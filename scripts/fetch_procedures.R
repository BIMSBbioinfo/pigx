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

#===================================================================
# Define functions:
#===================================================================


# Prepare Variable -----------------------------------------------------------

# load libraries
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
fetchTableFromUCSC <- function(table.name, filename, assembly) {
  require("rtracklayer")
  mySession = rtracklayer::browserSession("UCSC")
  rtracklayer::genome(mySession) <- assembly
  track.names <- rtracklayer::trackNames(rtracklayer::ucscTableQuery(mySession))
  
  if (table.name %in% track.names) {
      message(paste("Downloading", table.name, "..."))
  
      targetTrack <- rtracklayer::track(mySession, table.name)
      return(targetTrack)
  } else {
      warning(paste("Could not find", table.name, " track for the given assembly <'",assembly,"'>."))
      return(FALSE)
  }
}


#' Fetch a table from the Annotation hub
#'
#' @param table.name a character indicating name of a UCSC table
#' @param table.loc a character indicating path to the output BED file
#' @param assembly a character indicating a genome version, e.g. "ce10""
#'
#' @return
#' @export
#'
#' @examples
fetchTableFromAnnotationHub <- function(table.name, filename, assembly) {
  require("AnnotationHub")
  hub = AnnotationHub()
  
  ## query track for assembly
  track.q <- query(hub, c(table.name, "genes", assembly))
  
  ## If there is exactly one record: fetch it
  if(length(track.q) == 1) {
    message("Found single ", table.name ," track, downloading...\n")
    track <- hub[[names(track.q)]]
    return(track)
  } else {
    warning(paste("Could not find single", table.name, "track for the given assembly <'",assembly,"'>.\n"))
    return(NULL)
  }
}

#' Lookup annotation files
#' 
#' Check if annotation file already exist, if not download from UCSC
#'
#' @param type 
#' @param filename 
#' @param assembly 
#' @param webfetch 
#'
#' @return
#' @export
#'
#' @examples
lookupBedFile <- function(type, filename, assembly, webfetch) {
  
  if (file.exists(filename)) {
    return(filename)
  }
  
  if(!(type %in% c("cpgIslandExt", "knownGene"))) {
    stop('ERROR: type can be either "cpgIslandExt" or "knownGene".')
  }

  message("Trying to fetch from AnnotationHub.\n")
  
  track <- fetchTableFromAnnotationHub(type,filename,assembly)
  if(is.null(track)) {
    
    message("Trying to fetch from UCSC directly.\n")
    track <-fetchTableFromUCSC(type,filename,assembly)
    if(is.null(track)) {
      
      warning( "Failed to find reference annotation file",type," for assembly <'",assembly,"'>.",
               "Make sure you used a valid UCSC assembly for settings:general:assembly in settings file." )
      return(NULL)
    }
  }
  
  
  ## write gzipped file ?
  if (grepl(".gz$",filename)) {
    fileCon <- gzfile(filename)
  } else {
    fileCon <- file(filename)
  }
  
  ## else write it to BED file
  export.bed(object = track,
             con = fileCon,
             trackLine = FALSE)
  message(paste("Wrote ", table.name, " track to:", filename))
  
  return(filename)
    
  }
