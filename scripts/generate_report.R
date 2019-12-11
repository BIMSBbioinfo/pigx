# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
# Copyright © 2017, 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017, 2018, 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
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


## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included
#' Title
#'
#' @param reportFile rmarkdown file to render
#' @param outFile name of the output file
#' @param workdir working directory
#' @param report.params
#' @param logo Location of PiGx logo
#' @param quiet boolean value (default: FALSE). If set to TRUE, progress bars 
#'   and chunk labels will be suppressed while knitting the Rmd file.
#' @param selfContained boolean value (default: TRUE). By default, the generated
#'   html file will be self-contained, which means that all figures and tables 
#'   will be embedded in a single html file with no external dependencies (See 
#'   rmarkdown::html_document)
#' @param bibTexFile path to bibTex biblographie file
#' @param prefix prefix for report
#'   
#' @return An html generated using rmarkdown/knitr/pandoc that contains 
#'   interactive figures, tables, and text that provide an overview of the 
#'   experiment
#' @export
#'
#' @examples
render2HTML <- function(reportFile,
                        outFile,
                        workdir = getwd(),
                        report.params,
                        logo,
                        bibTexFile,
                        prefix,
                        selfContained=TRUE,
                        quiet = FALSE)
{
  
  
  if (is.null(report.params)) report.params <- list()
  
  ## render single report
  message("rendering report from template: ", reportFile)
 

  interdir <- file.path(workdir,prefix,"_tmp")
  # dir.create(interdir)

  htmlwidgets::setWidgetIdSeed(1234)
  rmarkdown::render(
    input = reportFile,
    output_dir = workdir,
    intermediates_dir = interdir,
    knit_root_dir = interdir,
	output_file = outFile,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      depth = 2,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = TRUE,
      code_folding = "hide",
      self_contained = selfContained,
      bibliography = bibTexFile
    ),
    params = c(report.params,
               logo = logo,
               prefix = prefix,
               workdir = workdir),
    quiet  = quiet,
    clean  = FALSE
  )

  # on.exit(unlink(interdir, recursive = TRUE),add = TRUE)
}




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
        --outFile    output file
        --workdir    working directory
        --logo       path to PiGx Logo
        --bibTexFile path to bibTex biblographie file
        --prefix     prefix for report
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

## catch output and messages into log file
if(!is.null(argsL$logFile)) {
  out <- file(argsL$logFile, open = "wt")
  sink(out,type = "output")
  sink(out, type = "message")
} else {
  sink()
}

# export a copy of the argument list for this rendering attempt
saveRDS( argsL, file=paste0(argsL$outFile,".RenderArgs.rds") )
message("Parameters")
print(argsL)

## get deeper list elements of report params
if(!is.null(argsL$report.params)) {
  argsL$report.params <- jsonlite::fromJSON(argsL$report.params)
}


## check wether all required report params are given 
paramsList <- knitr::knit_params(readLines(argsL$reportFile),
                                 evaluate = TRUE)

## exclude standard report params
paramsList <- paramsList[!names(paramsList) %in% c("logo","prefix","workdir")]

givenParams <- names(paramsList) %in% names(argsL$report.params)
if( !all(givenParams ))  {
  warning("Missing values for parameters: ",
          paste(names(paramsList)[!givenParams],collapse = ", "))
}


cat(paste(
    format(as.POSIXct(if ("" != Sys.getenv("SOURCE_DATE_EPOCH")) {
                          as.numeric(Sys.getenv("SOURCE_DATE_EPOCH"))
                      } else {
                          Sys.time()
                      }, origin="1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    "\n\n",
    "Rendering report:",basename(argsL$reportFile),"\n",
    "from template:",basename(argsL$outFile),"\n",
    "into directory:",argsL$workdir ,"\n\n"
))


render2HTML(reportFile = normalizePath(argsL$reportFile),
            outFile = basename(argsL$outFile),
            workdir = argsL$workdir,
            report.params = argsL$report.params,
            logo = argsL$logo,
            bibTexFile = argsL$bibTexFile,
            prefix = argsL$prefix)
