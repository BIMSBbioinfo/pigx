# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
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
      --rds name of the input RDS file containting the methylRaw object
      --grds name of output RDS file containing Segments as GRanges object
      --outBed name of output BED file containing Segments
      --png name of file to save diagnostic plots to   
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
out <- file(argsL$logFile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")



# Run Functions -----------------------------------------------------------

## Segmentation

## load methylKit
suppressPackageStartupMessages(library("methylKit"))

input     <- argsL$rds
output    <- argsL$outBed
grFile    <- argsL$grds
pngFile   <- argsL$png

## read input methylRaw
methRaw <- readRDS(input)

## convert to GRanges
methRaw.gr= as(methRaw,"GRanges")
## calculate methylation score 
mcols(methRaw.gr)$meth=100*methRaw.gr$numCs/methRaw.gr$coverage
##destrand
strand(methRaw.gr) <- "*"
##sort 
methRaw.gr <- sort(methRaw.gr[,"meth"]) 


### Segmentation of methylation profile

err <- tryCatch(expr = {
                ## try to run the code
                    png(filename = pngFile,
                        units = "in",width = 8,
                        height = 4.5,res=300)
                    res.gr = methSeg(methRaw.gr,
                                     diagnostic.plot=TRUE)
                    dev.off()},
                error = function(x) {
                ## if it fails still generate empty output
                    file.create(grFile)
                    file.create(output)
                    on.exit(file.remove(pngFile))
                    stop(paste("error occured!!",x))
                })

## Saving object
saveRDS(res.gr,file=grFile) 


### Export

## export segments to bed file
methSeg2bed(segments = res.gr,
            trackLine = paste0("track name='meth segments ' ",
                               "description='meth segments of ",
                               methRaw@sample.id,
                               " mapped to ",
                               methRaw@assembly,
                               "' itemRgb=On"),
            colramp=colorRamp(c("gray","green", "darkgreen")),
            filename = output)

