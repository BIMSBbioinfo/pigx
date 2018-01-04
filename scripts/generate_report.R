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
        --scriptsDir location of R scripts
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

source(paste0(argsL$scriptsDir, "/report_functions.R"))

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
