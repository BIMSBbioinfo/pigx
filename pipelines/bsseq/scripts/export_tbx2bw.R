# PiGx BSseq Pipeline.
#
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

#export_tbx2bw.R - takes a methylraw tabix file and a tabulated list of chromosome lengths and outputs a bigwig file
# ---last updated Sep. 2019 by A. Blume

#-------------------------------------------------------------------------

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Export tabix files to methylation bigwig
      
      Arguments:
      --filepath        path to tabix file
      --seqlengths_path path to chrominfo file
      --assembly        genome assembly
      --destrand        wether to merge methCalls from both strands 
      --out_path        output path of bigwig file
      --logFile         file to print the logs to
      --help            help 
      
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
} else {
  sink()
}

# Run Functions -----------------------------------------------------------


suppressPackageStartupMessages(expr = {
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(methylKit)
  library(rtracklayer)
  library(data.table)
})

data.table::setDTthreads(8)

# filepath        <-  "/data/local/agosdsc/projects/Animesh_bsseq/pigx-bsseq-results/06_methyl_calls/methylDackel/tabix_CpG/RP6.deduped_CpG.txt.bgz"
# seqlengths_path <-  "/data/local/agosdsc/projects/Animesh_bsseq/pigx-bsseq-results/04_mapping/Refgen_mm9_chromlengths.csv"
# assembly        <-  "mm9"
# out_path        <-  "/clusterhome/agosdsc/projects/pigx/pigx_bsseq/my.bw"
# destrand        <-  TRUE

filepath        <- argsL$filepath 
seqlengths_path <- argsL$seqlengths_path 
assembly        <- argsL$assembly 
out_path        <- argsL$out_path 
destrand        <- argsL$destrand


destrand <- ifelse(tolower(destrand) %in% c("true","yes"),TRUE,FALSE)


# ---------------------------------------------------

message("Parsing chromosome lengths.")
# load reference sequence info containing  
seqdat_temp = read.csv(seqlengths_path, sep="\t", header=FALSE)
Sinfo <- GenomeInfoDb::Seqinfo(seqnames   = as.character(seqdat_temp[,1]),
                 seqlengths = seqdat_temp[,2],
                 genome     = assembly)

message("Extracting methylation value from Tabix file.")

SinfoList <- split(as(Sinfo,"GRanges"),seqnames(Sinfo))


methList <- lapply(SinfoList, 
       FUN = function(gr) {

         message("Processing chromosome ",seqnames(gr),"...")

            # read directly from tabix file and process in chunks
            dt <- methylKit:::applyTbxByOverlap(
              tbxFile = filepath,
              ranges = gr, 
              return.type = "data.table", 
              chunk.size = 1e9, 
              FUN = function(dt){
                options(scipen = 999)
                methylKit:::.setMethylDBNames(dt)
                # merge strands if destrand==TRUE
                if(destrand) dt <- setDT(methylKit:::.CpG.dinuc.unify(dt))
                dt[,score := numCs/coverage]
                dt[,c("coverage","numCs","numTs") := NULL] 
                return(dt)    
                })
            return(GenomicRanges::makeGRangesFromDataFrame(
              dt,seqinfo = Sinfo,  keep.extra.columns=TRUE)
              )
        })

methList <- unlist(GRangesList(methList))

rtracklayer::export.bw( object = methList, con = out_path )

# bigwig exported. Program complete.
