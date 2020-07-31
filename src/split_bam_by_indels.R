# crispr-DART pipeline
#
# Copyright © 2017-2020 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the crispr-DART pipeline
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

library(GenomicAlignments)
library(Rsamtools)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

bamFile <- args[1]
sampleName <- args[2]
outDir <- args[3]

get_reads_with_indels <- function(aln) {
  dt <- data.table::data.table('cigar' = cigar(aln), 
                               'readID' = mcols(aln)$qname)
  #find reads with at least one insertions or deletions  
  indel_reads <- dt[which(stringi::stri_count(regex = "D|I", dt$cigar) > 0),]$readID

  return(indel_reads)
}


#parse alignments from bam file 
aln <- GenomicAlignments::readGAlignments(bamFile, param = ScanBamParam(what=c("qname", "seq")))

#get all reads with (insertions/deletions/substitions)
indelReads <- get_reads_with_indels(aln) 

filter_indels <- FilterRules(list(subset_reads=function(x) x$qname %in% indelReads))
indexBam(filterBam(file = bamFile,
                   destination = BamFile(file.path(outDir, paste0(sampleName, ".with_indels.bam"))),
                   filter=filter_indels,
                   param=ScanBamParam(what="qname")))


filter_no_indels <- FilterRules(list(subset_reads=function(x) ! x$qname %in% indelReads))
indexBam(filterBam(file = bamFile,
                   destination = BamFile(file.path(outDir, paste0(sampleName, ".without_indels.bam"))),
                   filter=filter_no_indels,
                   param=ScanBamParam(what="qname")))

