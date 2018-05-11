#' runReport
#' 
#' Generate a CRISPR mutation analysis report in a self-contained HTML file
#'
#'
#' @param reportFile Path to .Rmd script to generate a HTML report
#' @param sampleSheetFile Path to sample sheet used to run the pigx_crispr
#'   pipeline
#' @param ampliconFastaFile Path to the fasta file for the targeted amplicon
#' @param ampliconName Name of the amplicon (should be same as written
#' in sampleSheetFile amplicon column, and the header in the fasta file)
#' @param indelsFolder Path to the folder containing output of getIndelStats.R
#' script, which includes bedgraph, bed, and tsv files related to indels
#' @param cutSitesFile Path to the tab separated two column file containing
#' cut sites (target sites of the designed sgRNAs)
#' @param workdir Path to working directory where the output files will be
#'   written
#' @param prefix Prefix to be attached to the beginning of output files
#' @return An html generated using rmarkdown/knitr/pandoc that contains
#'   interactive figures, tables, and text that provide an overview of the
#'   experiment
runReport <- function(reportFile, 
                      sampleSheetFile,
                      ampliconFastaFile,
                      ampliconName, 
                      cutSitesFile,
                      indelsFolder, 
                      workdir,
                      prefix
                      ) {
  
  outFile <- paste0(prefix, '.report.html')
  tabsetDropdown <- file.path(dirname(reportFile), 'tabset-dropdown.html')
  
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
      number_sections = TRUE,
      includes = rmarkdown::includes(after_body = tabsetDropdown)
    ),
    output_options = list(self_contained = TRUE),
    params = list(ampliconFastaFile = ampliconFastaFile,
                  ampliconName = ampliconName, 
                  sampleSheetFile = sampleSheetFile,
                  indelsFolder = indelsFolder, 
                  cutSitesFile = cutSitesFile,
                  workdir = workdir,
                  prefix = prefix
                  ),
    quiet = TRUE
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
runReport: Generate a CRISPR mutation analysis report in a self-contained HTML file

Arguments:

--reportFile Path to .Rmd script to generate a HTML report
--sampleSheetFile Path to sample sheet (comma-separated format) used to run the pigx_crispr
pipeline
--ampliconFastaFile Path to the fasta file for the targeted amplicon
--ampliconName  Name of the amplicon (should be same as written
in sampleSheetFile amplicon column, and the header in the fasta file)
--cutSitesFile Path to the tab separated two column file containing
cut sites (target sites of the designed sgRNAs)
--indelsFolder Path to the folder containing output of getIndelStats.R
script, which includes bedgraph, bed, and tsv files related to indels
--prefix (Optional, default: 'comparison1') Prefix to be attached to the beginning 
of output files
--workdir (Optional, default: 'current working directory')

Example:
Rscript runReport.R \\\
--reportFile=./report.Rmd \\\
--sampleSheetFile=./sample_sheet.csv \\\
--ampliconFastaFile=./myGene.fasta \\\
--ampliconName=myGene \\\
--cutSitesFile=myGene.cutsites.gff \\\
--indelsFolder=./indels/ampliconName/ \\\
--workdir=`pwd` \\\
--prefix='myGene' \\\
"

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

if(!("sampleSheetFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: sampleSheetFile Provide the path to sample sheet file")
}

if(!("ampliconFastaFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: colDataFile. Provide the path to fasta file for the amplicon\n")
}

if(!("ampliconName" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: ampliconName Provide name of the analysed amplicon\n")
}

if(!("cutSitesFile" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: cutSitesFile Provide path to the cut sites file\n")
}

if(!("indelsFolder" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Missing argument: indelsFolder Provide the path to indels output folder\n")
}

if(!("prefix" %in% argsDF$V1)) {
  cat(help_command, "\n")
  prefix <- 'amplicon_analysis'
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


reportFile = argsL$reportFile
sampleSheetFile = argsL$sampleSheetFile
ampliconFastaFile = argsL$ampliconFastaFile
ampliconName = argsL$ampliconName
indelsFolder = argsL$indelsFolder
cutSitesFile = argsL$cutSitesFile

runReport(reportFile = reportFile, 
          sampleSheetFile = sampleSheetFile,
          ampliconFastaFile = ampliconFastaFile,
          ampliconName = ampliconName, 
          cutSitesFile = cutSitesFile, 
          indelsFolder = indelsFolder,
          prefix = prefix, 
          workdir = workdir)
