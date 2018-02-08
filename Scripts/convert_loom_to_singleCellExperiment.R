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
--outFile Path to the writing location of the output .RDS file 

Example:
Rscript convert_loom_to_singleCellExperiment.R --loomFile=foo.loom --outFile=bar.RDS"

## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

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
  stop("Missing argument: loomFile Provide the path to .loom file")
}

if(!("outFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: outFile Provide the path to the output RDS file")
}


loomFile = argsL$loomFile
outFile = argsL$outFile

#2. load libraries
library(HDF5Array)
library(SingleCellExperiment)
library(scater)
library(scran)
library(DelayedMatrixStats)

#3. Define functions

#' Import a loom file into SingleCellExperiment object
#' 
#' This function was borrowed from the script 'make-data.R' in the repository 
#' https://github.com/LTLA/HumanCellAtlasData/blob/master/inst/scripts/make-data.R
#' 
#' @param path Path to the .loom file to import into SingleCellExperiment object
#' @return A SingleCellExperiment object 
loom2sce <- function(path) {
  # Extracting the count matrix.
  mat <- HDF5Array(path, "matrix")
  mat <- t(mat)
  
  # Extracting the row and column metadata.
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
  
  # Extracting layers (if there are any)
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
  
  # Returning SingleCellExperiment object.
  sce <- SingleCellExperiment(c(matrix=mat, other.layers), rowData=row.df, colData=col.df)
  return(sce)
}

#' Given a DelayedMatrix object (genes on the row, cells on the column), calculate 
#' counts per million
#' 
#' @param dm A DelayedMatrix object containing counts of genes per cell
#' @return A DelayedMatrix object of normalized counts per million of genes per cell
getCPM <- function(dm) {
  cpm <- t((t(dm)/DelayedMatrixStats::rowSums2(t(dm)) * 10^6))
  return(cpm)
}

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

#4. start analysis

#4.1 import loom into SingleCellExperiment object
sce <- loom2sce(path = loomFile)
#find cells with zero expression for all genes 
zeros <- which(colSums(assays(sce)[[1]]) == 0)
#subset sce object to exclude those cells
sce <- sce[,-zeros]

#count matrix
counts <- assays(sce)[[1]]

#4.2 Normalize counts (get counts per million)
counts_per_million <- round(getCPM(dm = counts), 2)

#4.3 scale normalized counts
scaled_counts <- round(scaleDM(dm = counts_per_million), 2)

#4.4 save processed assay data
assays(sce) <- SimpleList("cnts" = counts, 
                          "cpm" = counts_per_million, 
                          "scale" = scaled_counts)

#4.5 get PCA results
sce <- scater::runPCA(object = sce, ncomponents = 10, exprs_values = 'scale')

#4.6 get t-SNE results
sce <- scater::runTSNE(object = sce, ncomponents = 2, use_dimred = 'PCA')#,exprs_values = 'PCA')

#4.7 get highly variable genes 
#fit <- scran::trendVar(x = sce, parametric=TRUE, assay.type = 'cpm', use.spikes = FALSE)
#decomp <- scran::decomposeVar(x = sce, fit = fit, assay.type = 'cpm')
#top.hvgs <- rownames(decomp[order(decomp$bio, decreasing=TRUE),])[1:50]
#rowData(sce) <- cbind(DataFrame(VariableGenes = top.hvgs))

#5 save SingleCellExperiment object in .RDS format
saveRDS(object = sce, file = outFile)


  
