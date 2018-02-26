# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
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

suppressPackageStartupMessages(expr = {
  require("AnnotationHub")
  require("rtracklayer")
})


#' Fetch a table from the UCSC browser and save it as a BED file
#'
#' The windows are color coded based on their score (methylation or differential
#' methylation value).
#'
#' @param table.name a character indicating name of a UCSC table
#' @param table.loc a character indicating path to the output BED file
#' @param assembly a character indicating a genome version, e.g. "ce10""
#'
#' @return location of a output BED file which can be visualized in the UCSC browser
#'
#'
#' @export
#' @docType methods
#' @rdname fetchTablefromUCSC
fetchTableFromUCSC <- function (table.name, table.loc=NULL, assembly) {
    mySession = browserSession("UCSC")
    genome(mySession) <- assembly
    track.names <- trackNames(ucscTableQuery(mySession))

    if (table.name %in% track.names) {
        message(paste("Downloading", table.name, "..."))

        targetTrack <- track(mySession, table.name)
        ## and write it to BED file
        export.bed(object = targetTrack,
                   con = table.loc,
                   trackLine = FALSE)
        message(paste("Wrote track to ", table.loc))
        return(table.loc)
    } else {
        stop(paste("Could not find", table.name, "for the given assembly <'",assembly,"'>."))
    }
}

lookupBedFile <- function (type, filename, assembly, webfetch) {
  ## import local bed file if available
  if (file.exists(filename)) {
    return(filename)
  }
  gzipped <- paste0(filename, ".gz")
  if (file.exists(gzipped)) {
    return(gzipped)
  }
  
  # can't fine file locally: now check if we should try to download it:
  if( webfetch )
  {
    message(paste0("Could not find ", filename, ".  Fetching from Internet."))
    if (type == "refGene") {
      message("Trying to fetch from AnnotationHub.\n")
      hub = AnnotationHub()
      
      ## query refseq genes for assembly
      refseq.q <- query(hub, c("refseq", "genes", assembly))
      
      ## If there is exactly one record: fetch it
      if(length(refseq.q) == 1) {
        message("Found single RefSeq track, downloading...\n")
        refGenes <- hub[[names(refseq.q)]]
        ## and write it to BED file
        export.bed(object = refGenes,
                   con = filename,
                   trackLine=FALSE)
        message(paste("Wrote RefSeq track to:", filename))
        return(filename)
      }
    }
    
    tryCatch({
        return(fetchTableFromUCSC(type, filename, assembly))
    }, error = function (msg) {
        message(paste0("Error while downloading from UCSC browser: ", msg))
    })
  }
  else
  {
  # print( paste("Failed to find reference annotation file",type," for <'",assembly,"'> (see settings:general in settings file.)." ))
  return('')
  }
}
