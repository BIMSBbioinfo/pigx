library(data.table)
library(GenomicAlignments)
library(ggplot2)

#read in output of samtools mpileup and parse_mpileup.py scripts and
#create a bedgraph file of per-base deletion scores
#also, print out some diagnostic plots
args = commandArgs(trailingOnly=TRUE)

bamFile <- args[1]
sampleName <- args[2]
outDir <- args[3]
cutSitesFile <- args[4] #path to sgRNA cutting sites for the target genome.
sgRNA_list <- args[5] # column (:) separated list of sgRNA ids (must match ids in cutSitesFile)
                      # that were used (or desired to be profiles) for the given sample.

cat("Running getIndelStats with arguments:",args,"\n")
#print bedgraph file of deletion ratios per base


#get the actual sequences that are inserted in alignments 
#' @param aln GAlignments object to extract inserted sequences  
#' @param ins insertion coordinates in the genome 
#' @return DNAStringSetList object; one DNAStringSet object per 
#'  read with at least one insertion ' 
getInsertedSequences <- function(aln, ins) {
  #get coordinates of insertions in the reads
  insertions <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar = cigar(aln), 
                                                              ops = c('I'),
                                                              with.ops = T)
  names(insertions) <- mcols(aln)$qname
  insertions <- insertions[lengths(insertions) > 0]
  
  sequences <- mcols(aln)$seq[match(names(insertions), mcols(aln)$qname)]

  insertedSequences <- Biostrings::extractAt(x = sequences, at = insertions)

  df <- data.frame('seqname' = as.character(seqnames(aln[match(as.character(mcols(ins)$name), 
                                                               mcols(aln)$qname)])), 
                   #here we use the genomic coordinate of the insertion rather than the 
                   #position in the query (read)
                   'start' = start(ins), 
                   'name' = as.character(mcols(ins)$name),
                   'insertedSequence' = paste(unlist(insertedSequences)),
                   'insertionWidth' = nchar(paste(unlist(insertedSequences))))
                   
  return(df)
}

getIndelReads <- function(aln) {
  indelReads <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(aln), 
                                                              ops = c('I', 'D'),
                                                              with.ops = T, 
                                                              pos = start(aln))
  names(indelReads) <- mcols(aln)$qname
  indelReads <- stack(indelReads)
  end(indelReads[which(names(indelReads) == 'I')]) <- start(indelReads[which(names(indelReads) == 'I')]) 
  
  #get seqnames fields 
  indelReads <- cbind(as.data.table(indelReads), as.data.table(mcols(indelReads)))
  indelReads <- merge(indelReads, 
                  data.table("seqnames" = as.data.table(seqnames(aln))$value, 
                             "name" = mcols(aln)$qname),
                  by = 'name')
  colnames(indelReads)[5] <- 'indelType'
  indelReads <- GenomicRanges::makeGRangesFromDataFrame(indelReads, keep.extra.columns = TRUE)
  
  return(indelReads)
}
#aln: GenomicAlignments object
#indelReads: GenomicRanges object (ins or del extracted from getIndelReads function)
getIndelScores <- function(aln, indelReads) {
  alnCov <- GenomicAlignments::coverage(aln)
  indelCov <- GenomicAlignments::coverage(indelReads)
  
  #for each chromosome, calculate per base score: indel count divided by alignment coverage
  scores <- sapply(names(alnCov), function(chr) {
    ac <- alnCov[[chr]] #alignment coverage
    ic <- indelCov[[chr]] #indel coverage
    
    #in case indel coverage doesn't cover the whole alignment 
    #fill in the remaining bases with 0 values
    ic <- c(ic, rep(0, length(ac) - length(ic)))
    scores <- as.numeric(ic)/as.numeric(ac)
    scores[is.na(scores)] <- 0
    return(Rle(scores))
  })
  return(RleList(scores))
}

#' @param cutStart expected cutting site start pos for the sgRNA
#' @param cutEnd expected cutting site end pos for the sgRNA
#' @param bamFile path to bam file
#' @param extend (integer, default 3 bp) extend the searching area from cutting
#'   sites in either direction
countEventsAtCutSite <- function(seqName, cutStart, cutEnd, bamFile, readsWithInDels, extend = 3) {
  if(cutEnd < cutStart) {
    stop("End position of cutting site must be larger than start position\n")
  }
  if(cutStart < 0 | cutEnd < 0) {
    stop("Start/End positions of cutting sites must be positive values")
  }

  cutStartExt <- cutStart - extend
  cutEndExt <- cutEnd + extend

  # find the total number of reads whose alignments overlap the extended region of the cut-site
  aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname",
                                                       which = GRanges(seqnames = seqName,
                                                                       ranges = IRanges(start = cutStartExt,
                                                                                        end = cutEndExt),
                                                                       strand = '*')))
  
  # find the number of reads with deletions that start or end within the
  # extended region of the cut-site

  indelsAtCutSites <- readsWithInDels[(start(readsWithInDels) >= cutStartExt & 
                                         start(readsWithInDels) <= cutEndExt) | 
                                        (end(readsWithInDels) >= cutStartExt & 
                                           end(readsWithInDels) <= cutEndExt),]

  stats <- data.frame('indel' = length(indelsAtCutSites), 
                      'del' = sum(names(indelsAtCutSites) == 'D'),
                      'ins' = sum(names(indelsAtCutSites) == 'I'),
                      'coverage' = length(aln))
  return(stats)
}

#parse alignments from bam file 
aln <- GenomicAlignments::readGAlignments(bamFile, param = ScanBamParam(what=c("qname", "seq")))
#export to bigwig 
rtracklayer::export.bw(object = GenomicAlignments::coverage(aln), 
                       con = file.path(outDir, paste0(sampleName, ".alnCoverage.bigwig")))

#get all reads with (insertions/deletions/substitions)
indelReads <- getIndelReads(aln) 
#get indel scores 
indelScores <- getIndelScores(aln, indelReads)
#export to bigwig
rtracklayer::export.bw(object = indelScores, 
                       con = file.path(outDir, paste0(sampleName, ".indelScores.bigwig")))

#split insertions and deletions and 
#get indel scores (indel count divided by alignment coverage) per base
ins <- indelReads[indelReads$indelType == 'I']
insScores <- getIndelScores(aln, ins)
#export to bigwig
rtracklayer::export.bw(object = insScores, 
                       con = file.path(outDir, paste0(sampleName, ".insertionScores.bigwig")))

del <- indelReads[indelReads$indelType == 'D']
delScores <- getIndelScores(aln, del)
#export to bigwig
rtracklayer::export.bw(object = delScores, 
                       con = file.path(outDir, paste0(sampleName, ".deletionScores.bigwig")))

#extract insertion sequences and print to file 
insertedSequencesFile <- file.path(outDir, paste0(sampleName,'.insertedSequences.tsv'))
if(length(ins) > 0) {
  write.table(x = getInsertedSequences(aln, ins), file = insertedSequencesFile, 
              sep = '\t', quote = F, row.names = F, col.names = T)
} else {
  write(x = paste('seqname', 'start', 'name', 'insertedSequence', 'insertionWidth', collapse  = '\t'), 
        file =  insertedSequencesFile, sep = '\t')
}


# cutSites <- read.table(cutSitesFile, stringsAsFactors = F)
# 
# cutSiteStats <- do.call(rbind, lapply(1:nrow(cutSites), function(i) {
#   x <- cutSites[i,]
#   guide <- x[[1]]
#   cutStart <- as.numeric(x[[2]])
#   cutEnd <- cutStart + 1
# 
#   stats <- countEventsAtCutSite(seqName = seqName,
#                              cutStart = cutStart,
#                              cutEnd = cutEnd,
#                              bamFile = bamFile,
#                              readsWithInDels = readsWithInDels,
#                              extend = 3)
#   
#   stats <- cbind(data.frame("sample" = sampleName,
#                          "sgRNA" = guide,
#                          "cutStart" = cutStart,
#                          "cutEnd" = cutEnd), stats)
#   
#   #efficiency : number of indels that originate around the cut-site
#   #             divided by the read coverage around the cut-site
#   stats$indelEfficiency <- round(stats$indel / stats$cov * 100, 2)
#   
#   return(stats)
# }))
# 
# write.table(x = cutSiteStats,
#             file = file.path(outDir, paste0(sampleName, '.indel_stats_at_cutsites.tsv')),
#             quote = F, sep = '\t', row.names = FALSE
#             )
# 


# filter out reads with at least 1 deletion and 1 insertion
# filter out reads that have substitutions adjacent  to deletions (suggest complex events)
readsToFilter <- function(aln) {
  dt <- data.table::data.table('cigar' = cigar(aln), 
                               'readID' = mcols(aln)$qname)
  #find reads with at least one insertion and one deletion
  ins_del <- dt[grepl(pattern = '(I.*D)|(D.*I)', x = dt$cigar)]$readID
  
  #find reads with substitutions adjacent to a deletion
  adj_del_sub <- dt[grepl(pattern = '(X[0-9]+D)|(D[0-9]+X)', x = dt$cigar)]$readID
  
  #take the union
  return(union(ins_del, adj_del_sub))
}

#TODO use readsToFilter function depending on the sequencing technology (don't apply to pacbio)
# indelReads <- indelReads[!indelReads$name %in% readsToFilter(aln)]

# collapse indel reads by coordinates and print them into BED files
indels <- as.data.table(indelReads)[,length(name),by = c('seqnames', 'start', 'end', 'indelType')][order(V1, decreasing = T)]
colnames(indels)[5] <- 'ReadSupport'
indels$name <- paste(indels$seqnames, indels$start, indels$end, indels$indelType, sep = ':')
#get number of reads that overlap each deletion
overlapCounts <- countOverlaps(split(GenomicRanges::GRanges(indels), indels$name), aln)
indels$coverage <- as.numeric(overlapCounts[indels$name])

#add a score column for visualization on IGV
indels$score <- indels$ReadSupport/max(indels$ReadSupport)*1000

#print summarized deletions to BED file
rtracklayer::export.bed(object = GenomicRanges::GRanges(indels[indelType == 'D']), 
                        con = file.path(outDir, paste0(sampleName, '.deletions.bed')), 
                        trackLine = new("BasicTrackLine", useScore = TRUE, name = 'deletions'))

#print summarized insertions to BED file
rtracklayer::export.bed(object = GenomicRanges::GRanges(indels[indelType == 'I']), 
                        con = file.path(outDir, paste0(sampleName, '.insertions.bed')), 
                        trackLine = new("BasicTrackLine", useScore = TRUE, name = 'insertions'))

#write collapsed indels to file
write.table(x = indels,
            file = file.path(outDir, paste0(sampleName, '.indels.tsv')),
            quote = F, sep = '\t', col.names = T, row.names = F)

#write indel reads to file
write.table(x = indelReads,
            file = file.path(outDir, paste0(sampleName, '.reads_with_indels.tsv')),
            quote = F, sep = '\t', col.names = T, row.names = F)

