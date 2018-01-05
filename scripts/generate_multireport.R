# PiGx BSseq Pipeline.
#
# Copyright © 2017, 2018 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017, 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017, 2018 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017, 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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
# terms of the GPL version 3.  The "merge_chapters2" procedure was
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
      Render multiple reports to one

      Arguments:
        --scriptsDir location of R scripts
        --index index template
        --finalOutput output file
        --finalReportDir final report location
        --references references template
        --sessioninfo sessioninfo template
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

source(paste0(argsL$scriptsDir, "/report_functions.R"))

## catch output and messages into log file
out <- file(argsL$logFile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")


cat(paste(
    format(as.POSIXct(if ("" != Sys.getenv("SOURCE_DATE_EPOCH")) {
                          as.numeric(Sys.getenv("SOURCE_DATE_EPOCH"))
                      } else {
                          Sys.time()
                      }, origin="1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    "\n\n",
    "Rendering report:",basename(argsL$finalOutput),"\n",
    "into directory:",normalizePath(dirname(argsL$finalReportDir)),"\n\n"
))



render2multireport(final_output = normalizePath(argsL$finalOutput),
                   finalreportdir = normalizePath(argsL$finalReportDir),
                   index = argsL$index,
                   references = argsL$references,
                   sessioninfo = argsL$sessioninfo)

