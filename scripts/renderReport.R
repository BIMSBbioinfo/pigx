# Author: BU
# Date: February, 2018
# This script takes as input an RDS file containing a SingleCellExperiment object 
# and renders the scrnaReport.Rmd script 

library(rmarkdown)

#1. Collect arguments
args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) == 0) {
  args <- c("--help")
}
help_command = "
renderReport.R: Run scrnaReport.Rmd script to render an HTML report

Arguments:
--reportFile Path to .loom file from the scRNA-seq experiment
--sceRdsFile Path to the RDS format file containing the SingleCellExperiment object
--covariates Comma-sepated list of factors to use when plotting PCA/t-SNE 
--workdir Path to working directory where the output files will be written
--prefix Prefix to use for output file

Example:
Rscript renderReport.R --reportFile=<path to scranReport.Rmd> \\\
                       --sceRdsFile=<path to sce.RDS> \\\
                       --covariates='Age,Sex,Treatment' \\\
                       --workdir=<path to working directory> \\\
                       --prefix=MyProjectName \n"
                        

## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
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

if(!("reportFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: reportFile. Provide the path to .Rmd file")
}

if(!("sceRdsFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: sceRdsFile. Provide the path to .Rds file 
      containing the SingleCellExperiment object.")
}

if(!("covariates" %in% argsDF$V1)) {
  cat(help_command, "\n")
  covariates <- 'sample_id'
  warning("No covariates provided. Will use 'sample_id' as the only covariate when plotting PCA/t-SNE\n")
} else {
  covariates <- argsL$covariates
}

if(!("prefix" %in% argsDF$V1)) {
  cat(help_command, "\n")
  prefix <- 'PiGx'
  warning("No prefix provided. Will use '",prefix,"' as the prefix to output files\n")
} else {
  prefix <- argsL$prefix
}

if(!("workdir" %in% argsDF$V1)) {
  workdir <- getwd()
  warning("No output folder provided. Setting working directory to current directory:\n",workdir,"\n")
} else {
  workdir <- argsL$workdir
  cat("setting working directory to ",workdir,"\n")
}

reportFile <- argsL$reportFile
sceRdsFile <- argsL$sceRdsFile
  
#2. Define report rendering function
#' runReport
#' 
#' Generate a PiGx-scRNAseq Report in a self-contained HTML file
#'
#' @param reportFile Path to .Rmd script to generate a HTML report
#' @param sceRdsFile Path to the RDS format file containing the
#'   SingleCellExperiment object
#' @param covariates Comma-sepated list of factors to use when plotting
#'   PCA/t-SNE
#' @param workdir Path to working directory where the output files will be
#'   written
#' @param prefix Prefix to be attached to the beginning of output files
#' @return An html generated using rmarkdown/knitr/pandoc that contains
#'   interactive figures, tables, and text that provide an overview of the
#'   experiment
runReport <- function(reportFile, 
                      sceRdsFile, 
                      covariates, 
                      workdir = getwd(),
                      prefix = 'PiGx') {
  
  outFile <- paste0(prefix, '.scRNA-Seq.report.html')
  
  htmlwidgets::setWidgetIdSeed(1234)
  rmarkdown::render(
    input = reportFile, 
    output_dir = workdir,
    intermediates_dir = file.path(workdir, prefix),
    clean = TRUE,
    output_file = outFile,
    output_format = rmarkdown::html_document(
      code_folding = 'hide', 
      depth = 2,
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = TRUE
    ),
    output_options = list(self_contained = selfContained),
    params = list(sceRdsFile = sceRdsFile,
                  covariates = covariates, 
                  outFile = outFile,
                  workdir = workdir),
    quiet = FALSE
  )
  
  if(dir.exists(file.path(workdir, prefix))) {
    unlink(file.path(workdir, prefix), recursive = TRUE)
  }
}

runReport(reportFile = reportFile, 
          sceRdsFile = sceRdsFile, 
          covariates = covariates,
          workdir = workdir,
          prefix = prefix)


