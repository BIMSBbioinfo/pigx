# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('convert_loom_to_singleCellExperiment')



# ---------------------------------------------------------------------------- #
#' Import a loom file into SingleCellExperiment object
#'
#' This function was borrowed from the script 'make-data.R' in the repository
#' https://github.com/LTLA/HumanCellAtlasData/blob/master/inst/scripts/make-data.R
#'
#' @param path Path to the .loom file to import into SingleCellExperiment object
#' @return A SingleCellExperiment object
loom2sce = function(path) {
  message("Importing",path,"into SingleCellExperiment object")
  message("Extracting the count matrix")
  mat = HDF5Array(path, "matrix")
  mat = t(mat)

  message("Extracting the row and column metadata")
  col.attrs = h5read(path, "col_attrs")
  if (length(col.attrs)) {
    col.df = data.frame(col.attrs)
  } else {
    col.df = DataFrame(matrix(0, nrow(mat), 0))
  }

  row.attrs = h5read(path, "row_attrs")
  if (length(row.attrs)) {
    row.df = data.frame(row.attrs)
  } else {
    row.df = NULL
  }


  message("Extracting layers (if there are any)")
  optional = h5ls(path)
  is.layer = optional$group=="/layer"
  if (any(is.layer)) {
    layer.names = optional$name[is.layer,]
    other.layers = vector("list", length(layer.names))
    names(other.layers) = layer.names

    for (layer in layer.names) {
      current = HDF5Array(path, file.path("/layer", layer))
      other.layers[[layer]] = t(current)
    }
  } else {
    other.layers = list()
  }

  message("Returning SingleCellExperiment object.")
  sce = SingleCellExperiment(c(matrix=mat, other.layers), rowData=row.df, colData=col.df)
  return(sce)
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column),
#' calculate counts per million
#'
#' @param dm A DelayedMatrix object containing counts of genes per cell
#' @return A DelayedMatrix object of normalized counts per million of genes per
#'   cell in log2 scale
getCPM = function(dm) {
  cpm = t((t(dm)/DelayedMatrixStats::rowSums2(t(dm)) * 10^6))
  return(log2(cpm+1))
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column), scale the
#' matrix values
#'
#' @param dm A DelayedMatrix object
#' @return A DelayedMatrix object of scaled and centered columns
scaleDM = function(dm) {
  #center
  #subtract column means from corresponding columns
  cdm = t(t(dm) - DelayedMatrixStats::rowMeans2(t(dm))) #centred delayed matrix

  #scale
  #divide the centered columns by their standard deviations
  scdm = t(t(cdm) / DelayedMatrixStats::rowSds(t(cdm))) #scaled centered delayed matrix
  return(scdm)
}

# ---------------------------------------------------------------------------- #
#' runPrComp
#'
#' run stats::prcomp to obtain PCA rotations and loadings
#'
#' Given a SingleCellExperiment object with at least one assay matrix named `expr_values`,
#' compute the PCA and loadings, update the SingleCellExperiment object with 'PCA' in
#' reducedDimensions slot of the object.
#'
#'  @param x A SingleCellExperiment object
#'  @param features A character vector of features (should be subset of rowData(x)$Genes)
#'  @param center whether to center the data
#'  @param scale whether to scale the data
#'  @param ncomponents Max number of principle components to look for
#'  @param expr_values The name of the assay data matrix in `x`
#'
#'  @return A matrix of dimensions [ncol(x), ncomponents] with attributes matrix
#'    of [len(features), ncol(x)]
runPrComp = function(x, features = NULL, center = FALSE, scale = FALSE, ncomponents = 2, expr_values = 'scale') {
  dm = assays(x)[[expr_values]]
  if(!is.null(features)){
    select = match(features, rowData(x)$Genes)
    dm = assay(x, expr_values)[select,]
    rownames(dm) = features
  }
  results = stats::prcomp(x = dm,
                       center = center,
                       scale. = scale,
                       rank. = ncomponents)
  PCA = results$rotation
  attr(PCA, 'gene.loadings') = results$x
  reducedDim(x, 'PCA') = PCA
  return(x)
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column),
#' calculate basic statistics and return a DataFrame object
#'
#' @param dm A DelayedMatrix object
#' @return A DataFrame object (nrow = number of genes)
#' of basic stats of genes versus cells
getGeneStats = function(m) {
  #ndetected -> number of cells detected per gene
  nCells = DelayedMatrixStats::colSums2(t(m) > 0)
  #max_gene -> maximum gene expression per gene detected in any cell
  maxGene = DelayedArray::rowMaxs(m)
  #mean_gene -> average expression per gene in all cells
  meanGene = DelayedMatrixStats::colMeans2(t(m))
  #mean_expr -> average expression per gene only in detected cells
  meanExpr = ifelse(nCells > 0, DelayedMatrixStats::colSums2(t(m))/nCells, 0)
  return(DataFrame('ndetected' = nCells,
                   'max_gene' = maxGene,
                   'mean_gene' = meanGene,
                   'mean_expr' = meanExpr))
}

#' Given a DelayedMatrix object (genes on the row, cells on the column),
#' calculate basic statistics and return a DataFrame object
#'
#' @param dm A DelayedMatrix object
#' @return A DataFrame object (nrow = number of cells)
#' of basic stats of genes versus cells
getCellStats = function(m) {
  #nGene -> number of genes detected per cell
  nGene = DelayedMatrixStats::rowSums2(t(m) > 0)
  #max_gene -> maximum gene expression detected per cell
  maxGene = DelayedArray::colMaxs(m)
  #mean_gene -> average gene expression detected per cell
  meanGene = DelayedMatrixStats::rowMeans2(t(m))
  #mean_expr -> average non-zero gene expression per cell
  meanExpr = DelayedMatrixStats::rowSums2(t(m))/nGene
  return(DataFrame('nGene' = nGene,
                   'max_gene' = maxGene,
                   'mean_gene' = meanGene,
                   'mean_expr' = meanExpr))
}


# ---------------------------------------------------------------------------- #
loomToSingleCellExperiment = function(
    loom_file         = NULL,
    sample_sheet_file = NULL,
    gtf_file          = NULL,
    outfile           = NULL,
    genome_version    = NULL,
    script_path       = NULL

){

  # -------------------------------------------------------------------------- #
  suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(DelayedMatrixStats)
    library(DelayedArray)
    library(data.table)
    library(rtracklayer)
    library(stringr)
  })


  message(date()," Importing loom file into SingleCellExperiment object")
  sce = loom2sce(path = loom_file)

  outname = str_replace(outfile,'.RDS','')
  saveRDS(object = sce, file = paste0(outname, '.raw.RDS'))
  #4.1.1 Subset the sce object, remove cells that don't exist in the metaDataFile
  # warning(date()," Removing cells that don't exist in the meta data file",metaDataFile)
  sample_sheet = data.table::fread(sample_sheet_file)
  sample_sheet$barcode = NULL
  sample_sheet$reads   = NULL
  sample_sheet = unique(sample_sheet)

  #4.1.2 remove cells with zero expression for all genes
  warning(date()," Removing cells with zero expression for all genes")
  zeros = which(DelayedMatrixStats::rowSums2(t(assays(sce)[[1]])) == 0)

  if(length(zeros) > 0){
    #subset sce object to exclude those cells
    sce = sce[,-zeros]
  }

  message(date()," Updating colData")
  #4.1.1 update colData with additional meta data
  colData(sce) = merge(colData(sce),
                       DataFrame(sample_sheet),
                       by   = 'sample_name',
                       sort = FALSE)


  #count matrix
  counts = assays(sce)[[1]]

  #4.2 Normalize counts (get counts per million)
  message(date()," Normalizing counts")
  counts_per_million = round(getCPM(dm = counts), 2)

  ##4.2.1 update cell stats using cpm values
  colData(sce) = cbind(colData(sce),
                        getCellStats(m = counts_per_million))


  #4.2.2 update gene stats using cpm values
  rowData(sce) = cbind(rowData(sce),
                        getGeneStats(m = counts_per_million))


  message(date()," Scaling normalized counts")
  #4.3 scale normalized counts
  scaled_counts = round(scaleDM(dm = counts_per_million), 2)


  #4.4 save processed assay data
  assays(sce) = SimpleList(
    "cnts"  = counts,
    "cpm"   = counts_per_million,
    "scale" = scaled_counts
  )

  saveRDS(object = sce, file = paste0(outname, '.intermediate.RDS'))


  # 4.8 Map gene ids to gene names and update rowData
  #read mappings from the gtf file
  message(date()," Importing GTF data to map gene ids to gene names")
  gtf = rtracklayer::import.gff(gtf_file, format = 'gtf')
  #update rowData slot with gene name column
  rowData(sce)$gene_name = gtf[match(rowData(sce)$gene_id, gtf$gene_id),]$gene_name

  #5 save SingleCellExperiment object in .RDS format
  saveRDS(object = sce, file = outfile)

}

# -------------------------------------------------------------------------- #
loomToSingleCellExperiment(
    loom_file         = argv$input[['infile']],
    outfile           = argv$output[['outfile']],
    genome_version    = argv$params[['genome_version']],
    script_path       = argv$params[['script']],
    sample_sheet_file = argv$params[['sample_sheet_file']],
    gtf_file          = argv$params[['gtf_file']]
)
