# Title     : runVariantReport_SarsCov2
# Objective : run variant reports
# Created by: vfs
# Based on the runDeseqReport.R from Bora Uyar
# Created on: 29.04.21

# files should be given with full pathes?
# ! signature mutation sources are hardcoded

runReport <- function(reportFile,
                      coverage_file,
                      sample_name,
                      outFile,
                      # logo,
                      selfContained = TRUE, 
                      quiet = FALSE) {
  
  outFile <- outFile
  
 # htmlwidgets::setWidgetIdSeed(1234)
  rmarkdown::render(
    input = reportFile,
  #  intermediates_dir = file.path(workdir, prefix),
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
    params = list(coverage_file = coverage_file,
                  outFile = outFile,
                  sample_name = sample_name
               #  logo = logo, 
                   ),
    quiet = quiet
  )
  
#   if(dir.exists(file.path(sample_dir, sample_name))) {
#   unlink(file.path(sample_dir, sample_name), recursive = TRUE)
#   }
# }

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

#print(class(myArgs <- unlist(strsplit(args, "--"))))

argsDF <- parseArgs(args)
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

reportFile <- argsL$reportFile
coverage_file <- argsL$coverage_file
outFile <- argsL$outFile
sample_name <- argsL$sample_name

runReport(reportFile = reportFile,
          coverage_file = coverage_file,
          sample_name = sample_name,
          outFile = outFile,
          # logo = logo,
          selfContained = selfContained)