# ---------------------------------------------------------------------------- #
#Author: BU
#Date: February, 2018
#This script takes a loom file, imports it into singleCellExperiment object.
#Then, it makes further processing on the matrix:
# - normalizes and scales the count matrix (DelayedMatrixStats package)
# - calculate reduced dimensions such as PCA, and t-SNE (scater package)
# - find variable genes  (scrna package)
# - assign cell cycle scores to each cell (scrna package)
#Finally, it saves everything as a single SingleCellExperiment object in .RDS format,
#   which should be the starting point for creating an HTML Report of the whole experiment.

# ---------------------------------------------------------------------------- #
#1. Collect arguments
args <- commandArgs(TRUE)

cat("arguments:", args,"\n")

## Default setting when no arguments passed
if(length(args) == 0) {
  args <- c("--help")
}

help_command = "
convert_loom_to_singleCellExperiment.R: Import a loom file into
SingleCellExperiment object and further processes the data to make it
ready for reporting.

Arguments:
--loomFile Path to .loom file from the scRNA-seq experiment
--sampleSheetFile Sample Sheet file used to initiate the pipeline
--gtfFile GTF genome annotation file used as input to the pipeline
--genomeBuild Genome build version (e.g. hg38, mm10)
--outFile Path to the writing location of the output .RDS file

Example:
Rscript convert_loom_to_singleCellExperiment.R --loomFile=foo.loom --metaDataFile bar.tsv --genomeBuild hg38 --outFile=foo.RDS"

## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

# ---------------------------------------------------------------------------- #
## Parse arguments (we expect the form --arg=value)
#parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
parseArgs <- function(x) {
  myArgs <- unlist(strsplit(x, "--"))
  myArgs <- myArgs[myArgs != '']
  #when no values are provided for the argument
  myArgs <- gsub(pattern = "=$", replacement = "= ", x = myArgs)
  myArgs <- as.data.frame(do.call(rbind, strsplit(myArgs, "=")))
  myArgs$V2 <- gsub(' ', '', myArgs$V2)
  return(myArgs)
}

argsDF <- parseArgs(args)
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(!("loomFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: loomFile. Provide the path to .loom file")
}

if(!("sampleSheetFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: sampleSheetFile. Provide the path to sampleSheetFile file")
}

if(!("gtfFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: gtfFile Provide the path to gtfFile file")
}

if(!("genomeBuild" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: genomeBuild. Provide the genome build version (e.g. hg38, mm10)")
}

if(!("outFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: outFile. Provide the path to the output RDS file")
}

loomFile     = argsL$loomFile
sampleSheetFile = argsL$sampleSheetFile
gtfFile = argsL$gtfFile
genomeBuild  = argsL$genomeBuild
outFile      = argsL$outFile

# ---------------------------------------------------------------------------- #
#2. load libraries
library(HDF5Array)
library(SingleCellExperiment)
library(scater)
library(scran)
library(DelayedMatrixStats)
library(DelayedArray)
library(data.table)
library(rtracklayer)

#3. Define functions
# ---------------------------------------------------------------------------- #
#' Import a loom file into SingleCellExperiment object
#'
#' This function was borrowed from the script 'make-data.R' in the repository
#' https://github.com/LTLA/HumanCellAtlasData/blob/master/inst/scripts/make-data.R
#'
#' @param path Path to the .loom file to import into SingleCellExperiment object
#' @return A SingleCellExperiment object
loom2sce <- function(path) {
  message("Importing",path,"into SingleCellExperiment object")
  message("Extracting the count matrix")
  mat <- HDF5Array(path, "matrix")
  mat <- t(mat)

  message("Extracting the row and column metadata")
  col.attrs <- h5read(path, "col_attrs")
  if (length(col.attrs)) {
    col.df <- data.frame(col.attrs)
  } else {
    col.df <- DataFrame(matrix(0, nrow(mat), 0))
  }

  row.attrs <- h5read(path, "row_attrs")
  if (length(row.attrs)) {
    row.df <- data.frame(row.attrs)
  } else {
    row.df <- NULL
  }


  message("Extracting layers (if there are any)")
  optional <- h5ls(path)
  is.layer <- optional$group=="/layer"
  if (any(is.layer)) {
    layer.names <- optional$name[is.layer,]
    other.layers <- vector("list", length(layer.names))
    names(other.layers) <- layer.names

    for (layer in layer.names) {
      current <- HDF5Array(path, file.path("/layer", layer))
      other.layers[[layer]] <- t(current)
    }
  } else {
    other.layers <- list()
  }

  message("Returning SingleCellExperiment object.")
  sce <- SingleCellExperiment(c(matrix=mat, other.layers), rowData=row.df, colData=col.df)
  return(sce)
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column),
#' calculate counts per million
#'
#' @param dm A DelayedMatrix object containing counts of genes per cell
#' @return A DelayedMatrix object of normalized counts per million of genes per
#'   cell in log2 scale
getCPM <- function(dm) {
  cpm <- t((t(dm)/DelayedMatrixStats::rowSums2(t(dm)) * 10^6))
  return(log2(cpm+1))
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column), scale the
#' matrix values
#'
#' @param dm A DelayedMatrix object
#' @return A DelayedMatrix object of scaled and centered columns
scaleDM <- function(dm) {
  #center
  #subtract column means from corresponding columns
  cdm <- t(t(dm) - DelayedMatrixStats::rowMeans2(t(dm))) #centred delayed matrix

  #scale
  #divide the centered columns by their standard deviations
  scdm <- t(t(cdm) / DelayedMatrixStats::rowSds(t(cdm))) #scaled centered delayed matrix
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
runPrComp <- function(x, features = NULL, center = FALSE, scale = FALSE, ncomponents = 2, expr_values = 'scale') {
  dm <- assays(x)[[expr_values]]
  if(!is.null(features)){
    select <- match(features, rowData(x)$Genes)
    dm <- assay(x, expr_values)[select,]
    rownames(dm) <- features
  }
  results <- stats::prcomp(x = dm,
                       center = center,
                       scale. = scale,
                       rank. = ncomponents)
  PCA <- results$rotation
  attr(PCA, 'gene.loadings') <- results$x
  reducedDim(x, 'PCA') <- PCA
  return(x)
}

# ---------------------------------------------------------------------------- #
#' Given a DelayedMatrix object (genes on the row, cells on the column),
#' calculate basic statistics and return a DataFrame object
#'
#' @param dm A DelayedMatrix object
#' @return A DataFrame object (nrow = number of genes)
#' of basic stats of genes versus cells
getGeneStats <- function(m) {
  #ndetected -> number of cells detected per gene
  nCells <- DelayedMatrixStats::colSums2(t(m) > 0)
  #max_gene -> maximum gene expression per gene detected in any cell
  maxGene <- DelayedArray::rowMaxs(m)
  #mean_gene -> average expression per gene in all cells
  meanGene <- DelayedMatrixStats::colMeans2(t(m))
  #mean_expr -> average expression per gene only in detected cells
  meanExpr <- ifelse(nCells > 0, DelayedMatrixStats::colSums2(t(m))/nCells, 0)
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
getCellStats <- function(m) {
  #nGene -> number of genes detected per cell
  nGene <- DelayedMatrixStats::rowSums2(t(m) > 0)
  #max_gene -> maximum gene expression detected per cell
  maxGene <- DelayedArray::colMaxs(m)
  #mean_gene -> average gene expression detected per cell
  meanGene <- DelayedMatrixStats::rowMeans2(t(m))
  #mean_expr -> average non-zero gene expression per cell
  meanExpr <- DelayedMatrixStats::rowSums2(t(m))/nGene
  return(DataFrame('nGene' = nGene,
                   'max_gene' = maxGene,
                   'mean_gene' = meanGene,
                   'mean_expr' = meanExpr))
}

#4. start analysis

# ---------------------------------------------------------------------------- #
#4.1 import loom into SingleCellExperiment object
message(date()," Importing loom file into SingleCellExperiment object")
sce <- loom2sce(path = loomFile)

saveRDS(object = sce, file = paste0(outFile, '.raw.RDS'))
#4.1.1 Subset the sce object, remove cells that don't exist in the metaDataFile
# warning(date()," Removing cells that don't exist in the meta data file",metaDataFile)
sample_sheet <- data.table::fread(sampleSheetFile)
sample_sheet$barcode = NULL
sample_sheet$reads   = NULL
sample_sheet = unique(sample_sheet)

#4.1.2 remove cells with zero expression for all genes
warning(date()," Removing cells with zero expression for all genes")
zeros <- which(DelayedMatrixStats::rowSums2(t(assays(sce)[[1]])) == 0)

if(length(zeros) > 0){
  #subset sce object to exclude those cells
  sce <- sce[,-zeros]
}

message(date()," Updating colData")
#4.1.1 update colData with additional meta data
colData(sce)$sample_name = sub('_[ACGT]+$','',colData(sce)$cell_id)

colData(sce) <- merge(colData(sce),
                     DataFrame(sample_sheet),
                     by   = 'sample_name',
                     sort = FALSE)
#count matrix
counts <- assays(sce)[[1]]

#4.2 Normalize counts (get counts per million)
message(date()," Normalizing counts")
counts_per_million <- round(getCPM(dm = counts), 2)

##4.2.1 update cell stats using cpm values
colData(sce) <- cbind(colData(sce),
                      getCellStats(m = counts_per_million))


#4.2.2 update gene stats using cpm values
rowData(sce) <- cbind(rowData(sce),
                      getGeneStats(m = counts_per_million))


message(date()," Scaling normalized counts")
#4.3 scale normalized counts
scaled_counts <- round(scaleDM(dm = counts_per_million), 2)


#4.4 save processed assay data
assays(sce) <- SimpleList("cnts"  = counts,
                          "cpm"   = counts_per_million,
                          "scale" = scaled_counts)

saveRDS(object = sce, file = paste0(outFile, '.intermediate.RDS'))


#4.5 compute gene variability across conditions
message(date()," Computing gene variability across cells")
fit    <- scran::trendVar(x = sce, assay.type = 'cpm', use.spikes = FALSE)
decomp <- scran::decomposeVar(x = sce, fit = fit, assay.type = 'cpm')
rowData(sce) <- cbind(rowData(sce),
                      DataFrame('fitted_variability' = decomp$bio))

max <- 500
if(max > nrow(sce)) {
  max <- nrow(sce)
}

rowvars <- DelayedMatrixStats::colVars(t(assay(sce, 'cpm')))
rowData(sce) <- cbind(rowData(sce),
                      DataFrame('row_variability' = rowvars))
topgenes <- rowData(sce)[order(rowData(sce)$row_variability, decreasing = T),][1:max,]$Genes

#4.6 get PCA results using genes with top row variance
message(date()," Computing PCA")
sce <- runPrComp(x = sce, expr_values = 'scale', ncomponents = 10,
            features = topgenes, center = FALSE, scale = FALSE)
saveRDS(object = sce, file =  paste0(outFile, '.intermediate.RDS'))

message(date()," Computing t-SNE")
#4.7 get t-SNE results
sce <- scater::runTSNE(object           = sce,
                       components       = 2,
                       use_dimred       = 'PCA',
                       rand_seed        = 1) # check_duplicates is now set to FALSE by the function itself

saveRDS(object = sce, file =  paste0(outFile, '.intermediate.RDS'))

#4.8 Assign cell cycle phase scores
skipCycleScore <- FALSE
if(grepl('hg[0-9]+', genomeBuild)) {
  organism <- 'human'
} else if(grepl('mm[0-9]+', genomeBuild)) {
  organism <- 'mouse'
} else {
  skipCycleScore <- TRUE
  warning("Skipping cell cycle score assignment as genome builds except for human and mouse are not supported")
}
if(skipCycleScore == FALSE) {
  cc.pairs <- readRDS(system.file("exdata", paste0(organism,"_cycle_markers.rds"), package="scran"))
  #check if the gene id namespace in the SingleCellExperimet object follows Ensembl style
  overlap <- sum(grepl('^ENS', rowData(sce)$Genes)) / nrow(sce)

  if(overlap > 0.5) { #require at least 50% overlap in the gene id overlap
    message(date()," Computing cell cycle phase scores")
    m = as.matrix(assays(sce)[['cpm']])
    assigned <- scran::cyclone(x = m,
                               gene.names = rowData(sce)$Genes,
                               pairs      = cc.pairs,
                               iter       = 100)
    phases <- assigned$phases
    phases[is.na(phases)] <- 'unknown'
    colData(sce) <- cbind(colData(sce),
                          DataFrame("CellCyclePhase" = phases))
  } else {
    warning("Skipping cell cycle score assignment as
            gene id namespace used in scran package doesn't match
            the gene id namespace found in the SingleCellExperiment object")
  }
}

# 4.8 Map gene ids to gene names and update rowData
#read mappings from the gtf file
message(date()," Importing GTF data to map gene ids to gene names")
gtf <- rtracklayer::import.gff(gtfFile, format = 'gtf')
#update rowData slot with gene name column
rowData(sce)$geneName <- gtf[match(rowData(sce)$Genes, gtf$gene_id),]$gene_name

#5 save SingleCellExperiment object in .RDS format
saveRDS(object = sce, file = outFile)
