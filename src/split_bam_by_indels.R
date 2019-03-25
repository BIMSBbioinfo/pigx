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

