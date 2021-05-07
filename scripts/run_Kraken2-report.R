# Title  : run Kraken2 report
# author: mfaxel
# Based on the runDeseqReport.R from Bora Uyar
# Created on: 2.05.21

#library(rmarkdown)

#function to run report
runReport <- function(reportFile, 
                      sample_name, 
                      kraken_output,
                      output_file,
                      selfContained = TRUE,
                      quiet = FALSE) {
  
#   outFile <- paste0(output_dir,"/",sample_name,"kraken_report.html")
  
  # htmlwidgets::setWidgetIdSeed(1234)
  rmarkdown::render(
    input = reportFile, 
    # output_dir = sample_dir,
    #  intermediates_dir = file.path(workdir, prefix),
    clean = TRUE,
    output_file = output_file,
    output_format = rmarkdown::html_document(
      code_folding = 'hide', 
      depth = 2,
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = TRUE
    ),
    output_options = list(self_contained = selfContained),
    params = list(kraken_output = kraken_output,
                  sample_name = sample_name
                  ),
    quiet = quiet
  )
  
#   if(dir.exists(file.path(sample_dir, sample_name))) {
#     unlink(file.path(sample_dir, sample_name), recursive = TRUE)
#   }
  
}

#1. Collect arguments
args <- commandArgs(TRUE)

cat("arguments:", args,"\n")

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

reportFile <- argsL$reportFile
sample_name <- argsL$sample_name
kraken_output <- argsL$kraken_output
output_file <- argsL$output_file

runReport(reportFile = reportFile,
          sample_name = sample_name,
          kraken_output = kraken_output,
          output_file = output_file)