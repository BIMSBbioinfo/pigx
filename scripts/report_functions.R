# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
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
#
# This file incorporates code by Yihui Xie <xie@yihui.name> under the
# terms of the GPL version 3.  The "Rscript_render2" procedure was
# adapted from
# https://github.com/rstudio/bookdown/blob/master/R/utils.R

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

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## get deeper list elements 
if(!is.null(argsL$report.params)) {
  argsL$report.params <- jsonlite::fromJSON(argsL$report.params)
}


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
  
  if(!dir.exists(finalReportDir)) dir.create(finalReportDir, recursive = TRUE)
  
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




# Call Functions ----------------------------------------------------------


## catch output and messages into log file
out <- file(argsL$logFile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")

cat(paste(
  format(as.POSIXct(if ("" != Sys.getenv("SOURCE_DATE_EPOCH")) { Sys.getenv("SOURCE_DATE_EPOCH") } else { Sys.time() }, origin="1970-01-01"), "%Y-%m-%d %H:%M:%S"),"\n\n",
  "Rendering report:",basename(argsL$reportFile),"\n",
  "from template:",basename(argsL$outFile),"\n",
  "into directory:",normalizePath(dirname(argsL$outFile)),"\n\n"
))

render2Markdown(reportFile = normalizePath(argsL$reportFile),
                outFile = basename(argsL$outFile),
                outDir = normalizePath(dirname(argsL$outFile)),
                finalReportDir = normalizePath(argsL$finalReportDir),
                report.params = argsL$report.params,
                logFile = argsL$logFile)
