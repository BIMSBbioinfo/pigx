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
      --inBam location of input bam file
      --assembly assembly used to map the reads
      --mincov minimum coverage (default: 10)
      --minqual minimum base quality (default: 20)
      --rds name of the RDS output file
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


### Methylation Calling

## load methylKit
library("methylKit")


input     <- argsL$inBam
assembly  <- argsL$assembly
mincov    <- as.numeric(argsL$mincov)
minqual   <- as.numeric(argsL$minqual)
rdsfile   <- argsL$rds

### Extract Methylation Calls

## extract the sample id from sample file 
sample_id <- gsub(".bam","",basename(input))

## define the location to save intermediate file
save_folder <- dirname(rdsfile)

## read bam file into methylKit object
methRaw = processBismarkAln(location = input,
                            sample.id = sample_id,
                            assembly = assembly,
                            mincov = mincov,
                            minqual = minqual,
                            save.context = "CpG",
                            save.folder = save_folder)

## Saving object
saveRDS(methRaw,file=normalizePath(rdsfile)) 
