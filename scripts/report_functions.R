#!/usr/bin/env Rscript



# Parse Command Line Arguments --------------------------------------------


## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Render to report
 
      Arguments:
        --reportFile report template
        --outFile output file
        --outDir output directory
        --finalReportDir final report location
        --report.params parameters to report, write in 
                        the form of 'p1:v1;p2:v2' 
        --logFile file to print the logs to
      --help              - print this text
 
      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
parseList <- function(x) unlist(strsplit(x,split = ";"))
parseListArgs <- function(x) strsplit(x,split = ":")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## get deeper list elements 
if(!is.null(argsL$report.params)) {
  
  ## construct the list of report params 
  argsDF2 <- as.data.frame(do.call("rbind",(parseListArgs(parseList(argsL$report.params)))))
  argsL2 <- as.list(as.character(argsDF2$V2))
  names(argsL2) <- argsDF2$V1
  
  ## convert numbers from string to integer
  argsL2 <- lapply(argsL2,
                   FUN = function(x){
                     ai <- suppressWarnings(as.integer(x))
                     ifelse(is.na(ai),x,ai)
                     }
                   )
  
  argsL$report.params <- argsL2
  
}

# print(argsL)

# Function Definitions ----------------------------------------------------



## modified from https://github.com/rstudio/bookdown/blob/master/R/utils.R#L205-L208
Rscript_render2 = function(file, ..., log.file=NULL) {
  args = shQuote(c(bookdown:::bookdown_file('scripts', 'render_one.R'), file, ...))
  ## we append the stderr and stdout to the log file
  if(!is.null(log.file)) args = paste(args,">>",log.file,"2>&1")
  if (bookdown:::Rscript(args)!= 0) stop('Failed to compile ', file)
}

## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included
render2Markdown <- function(reportFile,
                            outFile,
                            outDir,
                            finalReportDir,
                            report.params=NULL,
                            self.contained=FALSE,
                            logFile = NULL)
{
  
  output_format = rmarkdown::all_output_formats(reportFile, 'UTF-8')
  
  if(!dir.exists(finalReportDir)) dir.create(finalReportDir,recursive = FALSE)
  
  if(is.null(report.params)) report.params <- list()

  ## render single report
  message("render single report")

  ## make independent intermediate dirs
  interDir <- paste0(outDir,"/",outFile,"_tmp")

  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    intermediates_dir = interDir,
    output_file = outFile,
    knit_root_dir = outDir,
    output_format = rmarkdown::html_document(
      code_folding = "hide",
      code_download = TRUE,
      self_contained = self.contained
    ),
    params=c(report.params,
             list("sessioninfo"=TRUE,
                  "references"=TRUE)),
    quiet = FALSE,
    clean = TRUE
  )

  on.exit(unlink(interDir,recursive = TRUE),add = TRUE)

  ## render for multireport
  message("render for multireport")
  
  ## save a copy of render arguments in a temp file
  render_args = tempfile('render', tmpdir = finalReportDir, '.rds')
  on.exit(unlink(render_args), add = TRUE)
  saveRDS(
    list(output_format = "rmarkdown::html_notebook",
         params=c(report.params,
                  list("sessioninfo"=FALSE,
                       "references"=FALSE)),
         output_file = outFile,
         output_dir = finalReportDir,
         intermediates_dir = finalReportDir, 
         knit_root_dir = finalReportDir,
         clean = FALSE, 
         quiet = FALSE,
         envir = parent.frame()
    ),
    render_args
  )
  
  ## store metadata for Final Report
  render_meta = paste0(finalReportDir,"/knitr_meta.rds")
  
  Rscript_render2(file = reportFile,render_args,render_meta,log.file = logFile)
  ## move the sessioninfo to final folder
  session_file <- list.files(path = outDir,pattern = "session",full.names = TRUE)
  session_file <- session_file[endsWith(session_file,".rds")]
  if(length(session_file)!=0) {
    file.rename(from = session_file,
                to = paste0(finalReportDir,"/",basename(session_file)))
  }
}



# render2Report <- function(reportFile,
#                           outFile,
#                           outDir,
#                           report.params)
# {
#   
#   #print(getwd())
#   
#   
#   ## write stdout to log file
#   # sink(snakemake@log[[1]])
#   
#  
#   ## the logo is stored in the template directory
#   pathToLogo <- paste0(normalizePath(dirname(reportFile)),"/pigx_bsseq_logo.html")
#   
#   ## we set the knitr root dir to be the base directory,
#   ## such that all paths are relative from there
#   rootDir <- dirname(dirname(reportFile))
# 
#   interDir <- paste0(outDir,"/inter")
# 
#   
#   rmarkdown::render(
#     input = reportFile,
#     output_file = outFile,
#     output_dir = outDir,
#     # intermediates_dir = paste0(outDir,"/tmp"),
#     intermediates_dir = interDir,
#     knit_root_dir =  outDir,#rootDir
#     output_format = rmarkdown::html_notebook(
#       toc = TRUE,
#       toc_float = TRUE,
#       theme = 'lumen',
#       number_sections = FALSE,
#       code_folding = "hide",
#       self_contained = TRUE,
#       includes = list(in_header = pathToLogo)
#     ),
#     params = report.params,
#     quiet = FALSE,
#     clean = TRUE,
#     envir = new.env()
#   )
#   #unlink(paste0(outDir,"/tmp"),recursive = TRUE)
#   #unlink(list.files(path = outDir,pattern = "knit|utf8|nb_files"),recursive = TRUE)
#
#
#
#
# }


# Call Functions ----------------------------------------------------------



## catch output and messages into log file
# out <- file(snakemake@log[[1]], open = "wt")
out <- file(argsL$logFile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")

#
# Rscript = function(args, ...) {
#   system2(file.path(R.home('bin'), 'Rscript'), args, ...,stdout = out, stderr = out)
# }


## debugging
# save.image(file = "snakemakeObj.RData")

## check for filepaths to be normalized
# snakeParams <- snakemake@params[nchar(names(snakemake@params)) > 0]
# snakeParams <- lapply(snakeParams, function(x) {
#   if(class(x) == "character") return(normalizePath(x))
#   else return(x) })

# cat(paste(
#   Sys.time(),"\n\n",
#   "Rendering report:",basename(snakemake@output[["report"]]),"\n",
#   "from template:",basename(snakemake@input[["template"]]),"\n",
#   "into directory:",normalizePath(dirname(snakemake@output[["report"]])),"\n\n"
#   ))

cat(paste(
  Sys.time(),"\n\n",
  "Rendering report:",basename(argsL$reportFile),"\n",
  "from template:",basename(argsL$outFile),"\n",
  "into directory:",normalizePath(dirname(argsL$outFile)),"\n\n"
))


# render2Markdown(reportFile = normalizePath(snakemake@input[["template"]]),
#                 outFile = basename(snakemake@output[["report"]]),
#                 outDir = normalizePath(dirname(snakemake@output[["report"]])),
#                 finalReportDir = normalizePath(dirname(snakemake@output[["knitr_meta"]])),
#                 report.params = snakemake@params[nchar(names(snakemake@params)) > 0],
#                 logFile = snakemake@log[[1]])

render2Markdown(reportFile = normalizePath(argsL$reportFile),
                outFile = basename(argsL$outFile),
                outDir = normalizePath(argsL$outDir),
                finalReportDir = normalizePath(argsL$finalReportDir),
                report.params = argsL$report.params,
                logFile = argsL$logFile)

## remove empty intermediate file
# on.exit(unlink(snakemake@output[["knitr_meta"]]))


#load("snakemakeObj.RData")
