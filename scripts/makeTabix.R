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
      Create tabix files from methylDackel .methylKit files
      
      Arguments:
      --location location of input file
      --sample.id unique name of input sample
      --assembly assembly used to map the reads
      --treatment treatment of input sample
      --context context of methylation
      --mincov minimum coverage (default: 10)
      --dbdir name of the output folder
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
  out <- file(argsL$logFile, open = "at")
  sink(out,type = "output")
  sink(out, type = "message")
}



# Run Functions -----------------------------------------------------------


### Create tabix file from Methylation Calls

## load methylKit
# require("methylKit")

location  <- argsL$location 
sample.id <- argsL$sample.id     
assembly  <- argsL$assembly    
treatment <- argsL$treatment   
context   <- argsL$context
mincov    <- as.numeric(argsL$mincov)
dbdir     <- argsL$dbdir   

message("read file <",location,"> into methylKit object")

## read file into methylKit object
methRaw = methylKit::methRead(location = location,
                    sample.id = sample.id,
                    assembly = assembly,
                    treatment = treatment,
                    mincov = mincov,
                    context = context,
                    dbtype = 'tabix',
                    dbdir = dbdir)

