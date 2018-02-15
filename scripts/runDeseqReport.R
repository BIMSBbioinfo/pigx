# PiGx RNAseq Pipeline.
#
# Copyright © 2017 Bora Uyar <bora.uyar@mdc-berlin.de>
# Copyright © 2018 Jonathan Ronen <yablee@gmail.com>
# Copyright © 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx RNAseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
#' @param gtfFile Path to the GTF file that was used as reference to calculate
#' countDataFile
#' @param caseSampleGroups Comma separated list of sample group names (not 
#'   sample replicate names) that should be treated as 'case' groups (e.g. 
#'   mutant or treated samples)
#' @param controlSampleGroups Comma separated list of sample group names (not 
#'   sample replicate names) that should be treated as 'control' groups (e.g. 
#'   wild-type or untreated samples)
#' @param covariates Comma separated list of co-variates to control for in
#'   differential expression analysis (e.g. batch, age, temperature,
#'   sequencing_technology etc.). Must correspond to a column field in the
#'   colDataFile file).
#' @param workdir Path to working directory where the output files will be
#'   written
#' @param organism The organism for which the analysis is done (e.g. hsapiens, 
#'   mmusculus, celegans) via g:Profiler. This argument only affects GO term
#'   analysis results. If the organism is not supported or there is no internet
#'   connection, GO results will not be displayed.
#' @param prefix Prefix to be attached to the beginning of output files
#' @param logo Location of PiGx logo
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
                      gtfFile,
                      caseSampleGroups,
                      controlSampleGroups,
                      covariates,
                      geneSetsFolder,
                      workdir = getwd(),
                      organism, 
                      prefix,
                      logo,
                      selfContained = TRUE, 
                      quiet = FALSE) {
  
  outFile <- paste0(prefix, '.deseq.report.html')

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
    params = list(countDataFile = countDataFile,
                  colDataFile = colDataFile,
                  gtfFile = gtfFile, 
                  caseSampleGroups = caseSampleGroups,
                  controlSampleGroups = controlSampleGroups,
                  covariates = covariates,
                  prefix = prefix,
                  workdir = workdir,
                  logo = logo,
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
--gtfFile Path to the GTF file that was used as reference to calculate
countDataFile
--caseSampleGroups Comma separated list of sample group names (not 
sample replicate names) that should be treated as 'case' groups (e.g. 
mutant or treated samples)
--controlSampleGroups Comma separated list of sample group names (not 
sample replicate names) that should be treated as 'control' groups (e.g. 
wild-type or untreated samples)
--covariates Comma separated list of co-variates to control for in
differential expression analysis (e.g. batch, age, temperature,
sequencing_technology etc.). Must correspond to a column field in the 
colDataFile file).
--workdir (Optional, default: current working directory) Path to working directory 
where the output files will be written
--organism (Optional) The organism for which the analysis is done. Supported
organisms are 'human', 'mouse', 'worm', and 'fly'. This argument only
affects GO term analysis results. If the organism is not supported, GO
results will not be displayed.
--prefix (Optional, default: 'comparison1') Prefix to be attached to the beginning 
of output files
--workdir (Optional, default: 'current working directory')
--selfContained boolean value (default: TRUE). By default, the generated
html file will be self-contained, which means that all figures and tables 
will be embedded in a single html file with no external dependencies (See 
markdown::html_document)

Example:
Rscript runDeseqReport.R --reportFile=./deseqReport.Rmd \\\
--countDataFile=./sample_data/counts.tsv \\\
--colDataFile=./sample_data/colData.tsv \\\
--gtfFile=./Ensembl.Celegans.90.gtf \\\
--caseSampleGroups='spt.16_N2_L4, hmg.4_N2_L4' \\\
--controlSampleGroups='ctrl_N2_L4' \\\
--covariates='batch, temp' \\\
--workdir=`pwd` \\\
--organism='hsapiens' \\\
--prefix='spt-16_hmg-4_vs_ctrl' \\\
--logo='/usr/local/share/pigx_rnaseq/logo.png' \\\
--selfContained=TRUE"

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

if(!("gtfFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: gtfFile Provide the path to GTF file\n")
}

if(!("caseSampleGroups") %in% argsDF$V1) {
  cat(help_command, "\n")
  stop("Missing argument: caseSampleGroups Provide a comma separated list of sample group names that will be treated as 'case' samples")
}

if(!("controlSampleGroups") %in% argsDF$V1) {
  cat(help_command, "\n")
  stop("Missing argument: controlSampleGroups Provide a comma separated list of sample group names that will be treated as 'case' samples")
}

if(!("covariates") %in% argsDF$V1) {
  covariates <- ''
} else {
  covariates <- argsL$covariates
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

if(!("logo" %in% argsDF$V1)) {
  warning("No logo specified.\n")
} else {
  logo <- argsL$logo
}

reportFile = argsL$reportFile
countDataFile = argsL$countDataFile
colDataFile = argsL$colDataFile
gtfFile = argsL$gtfFile
caseSampleGroups = argsL$caseSampleGroups
controlSampleGroups = argsL$controlSampleGroups



runReport(reportFile = reportFile, 
          countDataFile = countDataFile, 
          colDataFile = colDataFile, 
          gtfFile = gtfFile,
          caseSampleGroups = caseSampleGroups, 
          controlSampleGroups = controlSampleGroups, 
          covariates = covariates,
          workdir = workdir, 
          organism = organism,
          prefix = prefix,
          logo = logo,
          selfContained = selfContained)
