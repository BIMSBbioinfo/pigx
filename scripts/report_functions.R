
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

## catch output and messages into log file
out <- file(snakemake@log[[1]], open = "wt")
sink(out,type = "output")
sink(out, type = "message")

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
