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

suppressPackageStartupMessages(expr = {
  require("AnnotationHub")
  require("rtracklayer")
})


## this function tries to fetch the reference genes for the given assembly
fetchRefGene <- function(refgenes.loc = NULL,
                         assembly) {
  
  if(is.null(refgenes.loc)) refgenes.loc <- paste0("refseq.genes.",assembly,".bed")
  
  ## import local bed file if available
  if( file.exists(refgenes.loc) ) {
    ## parse it 
    message("Found RefSeq track at:")
    return(refgenes.loc)
    
  } else {
    message("Trying to fetch from AnnotationHub.\n")
    ## else query it from AnnotationHub 
    ah = AnnotationHub()
    ## query refseq genes for assembly
    refseq.q <- query(ah,c("refseq","genes",assembly))
    ## either there is exactly one record, so fetch it
    if(length(refseq.q) == 1) {
      message("Found single RefSeq track, downloading...\n")
      refGenes <- ah[[names(refseq.q)]]
      ## and write it to BED file
      export.bed(object = refGenes,
                 con = refgenes.loc,
                 trackLine=FALSE)
      message("Written the RefSeq track to:")
      return(refgenes.loc)
      
    } else if ( length(refseq.q) == 0 ) { 
      message("Trying to fetch from UCSC table browser.\n")
      ## or there is none, 
      ## so we check with rtracklayer for the latest ucsc data
      mySession = browserSession("UCSC")
      genome(mySession) <- assembly
      track.names <- trackNames(ucscTableQuery(mySession))
      # I am interested in the refGene track 
      if("refGene" %in% track.names) {
        message("Found single RefSeq track, downloading...\n")
        # fetch it as a GRanges object 
        targetTrack <- track(mySession,"refGene")
        ## and write it to BED file
        export.bed(object = targetTrack,
                   con = refgenes.loc,
                   trackLine=FALSE)
        
        message("Written the RefSeq track to:")
        return(refgenes.loc)
      } 
    } else {
      stop(paste("Could not find reference gene set for the given assembly <'",assembly,"'>." ))
    }
    
  }
}

args <- commandArgs(trailingOnly = TRUE)
logFile     <- args[1]
refgenesLoc <- args[2]
assembly    <- args[3]

## catch output and messages into log file
out <- file(logFile, open = "wt")
sink(out, type = "output")
sink(out, type = "message")

fetchRefGene(refgenes.loc = refgenesLoc, assembly = assembly)
