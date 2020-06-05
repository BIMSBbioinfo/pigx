# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
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
      Segment methylation profile using methylKit

      Arguments:
      --tabix     name of the input tabix file containting the methylRaw object
      --outBed    name of output BED file containing Segments
      --png       name of file to save diagnostic plots to
      --sample.id sample name
      --assembly  genome assembly
      --logFile   file to print the logs to
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
out <- file(argsL$logFile, open = "at")
sink(out,type = "output")
sink(out, type = "message")



# Run Functions -----------------------------------------------------------

## Segmentation

## load methylKit
suppressPackageStartupMessages(library("methylKit"))

input     <- argsL$tabix
output    <- argsL$outBed
pngFile   <- argsL$png
sample.id <- argsL$sample.id
assembly  <- argsL$assembly


message("Reading tabix file.")
## read input tabix to methylRawDB
methRawDB <- methRead(location=input,
                      sample.id =sample.id,
                      assembly = assembly,
                      dbtype ="tabix")

# ## catch a possible error and touch empty files
# ## to trigger successful run
# err <- tryCatch(
#   expr = {
    ## try to run the code
    png(filename = pngFile,
        units = "in",width = 8,
        height = 4.5,res=300)
    
    message("Performing segmentation...")
    ### Segmentation of methylation profile
    res.gr = methSeg(methRawDB,
                     diagnostic.plot=TRUE)
    
    dev.off()

    ### Export

    message("Exporting segmentation...")
    ## export segments to bed file
    methSeg2bed(segments = res.gr,
                trackLine = paste0("track name='meth segments ' ",
                                   "description='meth segments of ",
                                   methRawDB@sample.id,
                                   " mapped to ",
                                   methRawDB@assembly,
                                   "' itemRgb=On"),
                colramp=colorRamp(c("gray","green", "darkgreen")),
                filename = output)
#   },
#   error = function(x) {
#     ## if it fails still generate empty output
#     file.create(output)
#     message(paste("error occured!!",x))
#   }
# )



