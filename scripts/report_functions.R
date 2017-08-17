
## modified from http://github.com/rstudio/bookdown/blob/master/R/utils.R#L201-L203
# Rscript = function(args, ...) {
#   message("using own version")
#   message(...)
#   system2(file.path(R.home('bin'), 'Rscript'), args, ...)
# }

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
      # toc = TRUE,
      # toc_float = TRUE,
      # theme = 'lumen',
      # number_sections = FALSE,
      code_folding = "hide",
      code_download = TRUE,
      self_contained = self.contained#,
      #includes = list(in_header = pathToLogo)
    ),
    params=c(report.params,
             list("sessioninfo"=TRUE,
                  "references"=TRUE)),
    quiet = FALSE,
    clean = TRUE#,
    #envir = new.env()
  )

  on.exit(unlink(interDir,recursive = TRUE),add = TRUE)

  ## render for multireport
  message("render for multireport")
  
  ## save a copy of render arguments in a temp file
  render_args = tempfile('render', tmpdir = finalReportDir, '.rds')
  on.exit(unlink(render_args), add = TRUE)
  saveRDS(
    list(output_format = "rmarkdown::html_notebook",# rmarkdown::html_notebook(
         #   toc = TRUE,
         #   toc_float = TRUE,
         #   theme = 'lumen',
         #   number_sections = FALSE,
         #   code_folding = "hide",
         #   self_contained = TRUE,
         #   includes = list(in_header = pathToLogo)
         # ), 
         #knit_root_dir = outDir,
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
  render_meta = paste0(finalReportDir,"/knitr_meta.rds")#,bookdown:::with_ext(outFile, '.rds'))
  # render_meta = paste0(finalReportDir,"/",bookdown:::with_ext(outFile, '.rds'))
  
  
  Rscript_render2(file = reportFile,render_args,render_meta,log.file = logFile)
  ## move the sessioninfo to final folder
  session_file <- list.files(path = outDir,pattern = "session",full.names = TRUE)
  session_file <- session_file[endsWith(session_file,".rds")]
  if(length(session_file)!=0) {
    #if(file.exists(session_file)) {
      file.rename(from = session_file,
                to = paste0(finalReportDir,"/",basename(session_file)))
    #}
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



## catch output and messages into log file
out <- file(snakemake@log[[1]], open = "wt")
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

cat(paste(  
  Sys.time(),"\n\n",
  "Rendering report:",basename(snakemake@output[["report"]]),"\n",
  "from template:",basename(snakemake@input[["template"]]),"\n",
  "into directory:",normalizePath(dirname(snakemake@output[["report"]])),"\n\n"
  ))


render2Markdown(reportFile = normalizePath(snakemake@input[["template"]]),
                outFile = basename(snakemake@output[["report"]]),
                outDir = normalizePath(dirname(snakemake@output[["report"]])),
                finalReportDir = normalizePath(dirname(snakemake@output[["knitr_meta"]])),
                report.params = snakemake@params[nchar(names(snakemake@params)) > 0],
                logFile = snakemake@log[[1]])

## remove empty intermediate file
on.exit(unlink(snakemake@output[["knitr_meta"]]))


#load("snakemakeObj.RData")
