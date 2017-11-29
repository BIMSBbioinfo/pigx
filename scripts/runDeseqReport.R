
#' runReport
#' 
#' Generate a DESeq2 Report in a self-contained HTML file
#' 
#' 
#' @param reportFile Path to .Rmd script to generate a HTML report
#' @param countDataFile Path to count data file which contains raw read counts 
#'   per gene/transcript for each sample replicate
#' @param colDataFile Path to a tab-separated file with the experimental set-up 
#'   description. The row-names are the sample names, and the columns consist of
#'   meta-data such as sample group, batch, sequencing run, condition, treatment
#'   etc.
#' @param caseSampleGroups Comma separated list of sample group names (not 
#'   sample replicate names) that should be treated as 'case' groups (e.g. 
#'   mutant or treated samples)
#' @param controlSampleGroups Comma separated list of sample group names (not 
#'   sample replicate names) that should be treated as 'control' groups (e.g. 
#'   wild-type or untreated samples)
#' @param geneSetsFolder A folder with one or more .txt files, where each txt
#'   file contains a list of gene ids (subset of the gene ids used as row names
#'   in the count data table), one id per line. GAGE package is utilized to find
#'   out if the given gene set is collectively up/down regulated in the case
#'   samples compared to the control samples.
#' @param workdir Path to working directory where the output files will be
#'   written
#' @param organism The organism for which the analysis is done (e.g. hsapiens, 
#'   mmusculus, celegans) via g:Profiler. This argument only affects GO term
#'   analysis results. If the organism is not supported or there is no internet
#'   connection, GO results will not be displayed.
#' @param prefix Prefix to be attached to the beginning of output files
#' @param quiet boolean value (default: FALSE). If set to TRUE, progress bars 
#'   and chunk labels will be suppressed while knitting the Rmd file.
#' @param selfContained boolean value (default: TRUE). By default, the generated
#'   html file will be self-contained, which means that all figures and tables 
#'   will be embedded in a single html file with no external dependencies (See 
#'   rmarkdown::html_document)
#' @return An html generated using rmarkdown/knitr/pandoc that contains 
#'   interactive figures, tables, and text that provide an overview of the 
#'   experiment
runReport <- function(reportFile, 
                      countDataFile,
                      colDataFile,
                      caseSampleGroups,
                      controlSampleGroups,
                      geneSetsFolder,
                      workdir = getwd(),
                      organism, 
                      prefix, 
                      selfContained = TRUE, 
                      quiet = FALSE) {
  
  outFile <- paste0(prefix, '.deseq.report.html')
  
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
      theme = 'united',
      highlight = 'tango',
      number_sections = TRUE
    ),
    output_options = list(self_contained = selfContained),
    params = list(countDataFile = countDataFile,
                  colDataFile = colDataFile,
                  caseSampleGroups = caseSampleGroups,
                  controlSampleGroups = controlSampleGroups,
                  geneSetsFolder = geneSetsFolder, 
                  prefix = prefix,
                  workdir = workdir, 
                  organism = organism),
    quiet = quiet
  )
  
  if(dir.exists(file.path(workdir, prefix))) {
    unlink(file.path(workdir, prefix), recursive = TRUE)
  }
}

#1. Collect arguments
args <- commandArgs(TRUE)

cat("arguments:", args,"\n")

## Default setting when no arguments passed
if(length(args) == 0) {
  args <- c("--help")
}

help_command = "
runDeseqReport.R: Generate an HTML report based on the analysis of raw count data for differential expression analysis using DESeq2

Arguments:
--reportFile Path to .Rmd script to generate a HTML report
--countDataFile Path to count data file which contains raw read counts 
per gene/transcript for each sample replicate
--colDataFile Path to a tab-separated file with the experimental set-up 
description. The row-names are the sample names, and the columns consist of
meta-data such as sample group, batch, sequencing run, condition, treatment
etc.
--caseSampleGroups Comma separated list of sample group names (not 
sample replicate names) that should be treated as 'case' groups (e.g. 
mutant or treated samples)
--controlSampleGroups Comma separated list of sample group names (not 
sample replicate names) that should be treated as 'control' groups (e.g. 
wild-type or untreated samples)
--geneSetsFolder (Optional) A folder with one or more .txt files, where each txt
file contains a list of gene ids (subset of the gene ids used as row names
in the count data table), one id per line. GAGE package is utilized to find
out if the given gene set is collectively up/down regulated in the case
samples compared to the control samples.
--workdir (Optional, default: current working directory) Path to working directory 
where the output files will be written
--organism (Optional) The organism for which the analysis is done. Supported
organisms are 'human', 'mouse', 'worm', and 'fly'. This argument only
affects GO term analysis results. If the organism is not supported, GO
results will not be displayed.
--prefix (Optional, default: 'comparison1') Prefix to be attached to the beginning 
of output files
--selfContained boolean value (default: TRUE). By default, the generated
html file will be self-contained, which means that all figures and tables 
will be embedded in a single html file with no external dependencies (See 
markdown::html_document)

Example:
Rscript runDeseqReport.R --reportFile=./deseqReport.Rmd \\\
--countDataFile=./sample_data/counts.tsv \\\
--colDataFile=./sample_data/colData.tsv \\\
--caseSampleGroups='spt.16_N2_L4, hmg.4_N2_L4' \\\
--controlSampleGroups='ctrl_N2_L4' \\\
--geneSetsFolder='./sample_data/genesets' \\\
--workdir=`pwd` \\\
--organism='hsapiens' \\\
--prefix='spt-16_hmg-4_vs_ctrl' \\\
--selfContained=TRUE"

## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(!("reportFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: reportFile. Provide the path to .Rmd script")
}

if(!("countDataFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: countDataFile. Provide the path to counts data file")
}

if(!("colDataFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: colDataFile. Provide the path to colData file which defines the experimental design\n")
}

if(!("caseSampleGroups") %in% argsDF$V1) {
  cat(help_command, "\n")
  stop("Missing argument: caseSampleGroups Provide a comma separated list of sample group names that will be treated as 'case' samples")
}

if(!("controlSampleGroups") %in% argsDF$V1) {
  cat(help_command, "\n")
  stop("Missing argument: controlSampleGroups Provide a comma separated list of sample group names that will be treated as 'case' samples")
}

if(!("geneSetsFolder") %in% argsDF$V1) {
  cat(help_command, "\n")
  warning("Missing argument: geneSetsFolder. Will not do gene set enrichment analysis\n")
  geneSetsFolder <- ''
} else {
  geneSetsFolder <- argsL$geneSetsFolder
}

if(!("organism") %in% argsDF$V1) {
  cat(help_command, "\n")
  warning("Missing argument: organism Will skip GO term analysis\n")
  organism <- ''
} else {
  organism <- argsL$organism
}

if(!("prefix" %in% argsDF$V1)) {
  cat(help_command, "\n")
  prefix <- 'comparison1'
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

if(!("selfContained" %in% argsDF$V1)) {
  selfContained <- TRUE
} else {
  selfContained <- argsL$selfContained
}

reportFile = argsL$reportFile
countDataFile = argsL$countDataFile
colDataFile = argsL$colDataFile
caseSampleGroups = argsL$caseSampleGroups
controlSampleGroups = argsL$controlSampleGroups

runReport(reportFile = reportFile, 
          countDataFile = countDataFile, 
          colDataFile = colDataFile, 
          caseSampleGroups = caseSampleGroups, 
          controlSampleGroups = controlSampleGroups, 
          geneSetsFolder = geneSetsFolder,
          workdir = workdir, 
          organism = organism,
          prefix = prefix, 
          selfContained = selfContained)
