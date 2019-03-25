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
tech <- args[5] #the sequencing platform (illumina or pacbio)

outDir <- file.path(outDir, sampleName)
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
  
  #remove insertions that don't exist in "ins", which contains the genomic coordinate of all insertions
  #ins object may be missing some insertion events due to the readsToFilter function. 
  insertions <- insertions[unique(ins$name)]
  
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
getIndelScores <- function(alnCov, indelReads) {
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

# filter out reads with at least 1 deletion and 1 insertion
# filter out reads that have substitutions adjacent  to deletions (suggest complex events)
readsToFilter <- function(aln) {
  dt <- data.table::data.table('cigar' = cigar(aln), 
                               'readID' = mcols(aln)$qname)
  #find reads with at least one insertion and one deletion
  #ins_del <- dt[grepl(pattern = '(I.*D)|(D.*I)', x = dt$cigar)]$readID
  
  #find reads with more than one insertions or deletions  
  mul_indel <- dt[which(stringi::stri_count(regex = "D|I", dt$cigar) > 1),]$readID
    
  #find reads with substitutions adjacent to a deletion
  adj_del_sub <- dt[grepl(pattern = '(X[0-9]+D)|(D[0-9]+X)', x = dt$cigar)]$readID
  
  #take the union of all types of reads to remove
  toFilter <- union(mul_indel, adj_del_sub)
  
  return(toFilter)
}


#parse alignments from bam file 
aln <- GenomicAlignments::readGAlignments(bamFile, param = ScanBamParam(what=c("qname", "seq")))
alnCoverage <- GenomicAlignments::coverage(aln)
#export to bigwig 
rtracklayer::export.bw(object = alnCoverage, 
                       con = file.path(outDir, paste0(sampleName, ".alnCoverage.bigwig")))

#get all reads with (insertions/deletions/substitions)
indelReads <- getIndelReads(aln) 

#use readsToFilter function depending on the sequencing technology (don't apply to pacbio)
if(tech == 'illumina') {
  toFilter <-  readsToFilter(aln)
  message("Filtering out ",length(toFilter)," out of ",length(unique(indelReads$name)), 
          " indel reads (",round(length(toFilter) / length(unique(indelReads$name)) * 100, 2),"%)",
          " \n\tdue to multiple indel events per read or mismatches adjacent to deletions.")
  
  indelReads <- indelReads[!indelReads$name %in% toFilter]
}

#get indel scores 
indelScores <- getIndelScores(alnCoverage, indelReads)
#export to bigwig
rtracklayer::export.bw(object = indelScores, 
                       con = file.path(outDir, paste0(sampleName, ".indelScores.bigwig")))

#split insertions and deletions and 
#get indel scores (indel count divided by alignment coverage) per base
ins <- indelReads[indelReads$indelType == 'I']
insScores <- getIndelScores(alnCoverage, ins)
#export to bigwig
rtracklayer::export.bw(object = insScores, 
                       con = file.path(outDir, paste0(sampleName, ".insertionScores.bigwig")))

del <- indelReads[indelReads$indelType == 'D']
delScores <- getIndelScores(alnCoverage, del)
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

# calculate guide RNA efficiencies (mean indel score at cut site +/- extension)
cutSites <- rtracklayer::import(cutSitesFile, 'bed')
cutSites <- GenomicRanges::flank(cutSites, width = 5, both = TRUE)

cutSiteStats <- do.call(rbind, lapply(split(cutSites, cutSites$name), 
                       function(cs) {
                         chr <- as.character(seqnames(cs))
                         if(chr %in% names(indelScores)) {
                            scores <- as.numeric(indelScores[[chr]][start(cs):end(cs)])
                            df <- data.table('sgRNA' = cs$name, 
                                             'scores' = round(max(scores)*100,2))
                          } else {
                            return(NULL)
                          }
}))
cutSiteStats$sample <- sampleName

write.table(x = cutSiteStats,
            file = file.path(outDir, paste0(sampleName, '.sgRNA_efficiency.tsv')),
            quote = F, sep = '\t', row.names = FALSE
            )


# collapse indel reads by coordinates and print them into BED files
indels <- as.data.table(indelReads)[,length(name),by = c('seqnames', 'start', 'end', 'indelType')][order(V1, decreasing = T)]
colnames(indels)[5] <- 'ReadSupport'
indels$name <- paste(indels$seqnames, indels$start, indels$end, indels$indelType, sep = ':')
#get coverage values across the indel boundary (pick the maximum value) 
indels <- do.call(rbind, pbapply::pbsapply(simplify = F, USE.NAMES = T, 
                  X = names(split(indels, indels$seqnames)), 
                  FUN = function(chr) {
  alnCov <- as.vector(alnCoverage[[chr]])
  dt <- indels[seqnames == chr]
  #for each indel in chromosome:chr, find the max coverage value across boundaries
  if(nrow(dt) > 0){
    dt$coverage <- pbapply::pbapply(dt, 1, function(x) {
      return(max(alnCov[x[['start']]:x[['end']]], na.rm = TRUE))
    })
    return(dt)
  } else {
    return(NULL)
  }
}))

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

