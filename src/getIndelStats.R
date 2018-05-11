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


printBedGraphFile <- function(filepath, 
                              trackInfo,
                              scores) {
  #write bedgraph file for deletions
  trackDefinition <- paste0("track type=bedGraph name=",trackInfo)
  writeLines(text = trackDefinition, con = filepath)
  
  #copy base position to have start/end equal to base position
  scores <- scores[,c(1,2,2,3)] 

  #convert to zero-based index
  scores$bp <- as.numeric(scores$bp) - 1
  
  #bedgraph file
  write.table(x = scores, file = filepath, append = T,
              sep = '\t', quote = F, row.names = F, col.names = F)
}

printCoverageStats <- function(bamFile, sampleName, outDir = getwd()) {
  aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname"))
  
  indels <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(aln), 
                                                                       ops = c('I', 'D'),
                                                                       with.ops = T, pos = start(aln))
  names(indels) <- mcols(aln)$qname
  
  indels <- stack(indels)
  end(indels[which(names(indels) == 'I')]) <- start(indels[which(names(indels) == 'I')]) + 1
  seqinfo(indels) <- seqinfo(aln)
  
  del <- indels[which(names(indels) == 'D')]
  ins <- indels[which(names(indels) == 'I')]
  
  alnCoverage <- GenomicAlignments::coverage(aln)[[1]]
  
  delCoverage <- GenomicAlignments::coverage(del)
  #in case del coverage deosn't cover the whole alignment 
  #fill in the remaining bases with 0 values
  delCoverage <- c(delCoverage, rep(0, length(alnCoverage) - length(delCoverage)))

  insCoverage <- GenomicAlignments::coverage(ins)
  #fill in the remaining bases with 0 values
  insCoverage <- c(insCoverage, rep(0, length(alnCoverage) - length(insCoverage)))
  
  indelCoverage <- GenomicAlignments::coverage(indels)
  #fill in the remaining bases with 0 values
  indelCoverage <- c(indelCoverage, rep(0, length(alnCoverage) - length(indelCoverage)))
  

  df <- data.frame('seqname' = levels(seqnames(aln))[1], 
                   'sample' = sampleName, 
                   'bp' = 1:length(alnCoverage), 
                   'cov' = as.vector(alnCoverage),
                   'del' = as.vector(delCoverage),
                   'ins' = as.vector(insCoverage),
                   'indel' = as.vector(indelCoverage))
  
  #add some more stats: 
  df$delRatio <- ifelse(df$cov > 0, df$del/df$cov, 0)
  df$insRatio <- ifelse(df$cov > 0, df$ins/df$cov, 0)
  df$indelRatio <- ifelse(df$cov > 0, df$indel/df$cov, 0)

  #print bedgraph file 
  indelScoresFile <- file.path(outDir, paste0(sampleName,'.indelScores.bedgraph'))
  printBedGraphFile(file = indelScoresFile, 
                    trackInfo = paste(sampleName, 'indel score (insertions + deletions / coverage per base)'), 
                    scores = df[,c('seqname', 'bp', 'indelRatio')])
  
  #print coverage stats to file
  statsOutputFile <- file.path(outDir, paste0(sampleName,'.coverageStats.tsv'))
  write.table(x = df, file = statsOutputFile, append = T,
              sep = '\t', quote = F, row.names = F, col.names = T)
  return(indels)
} 

summarizeInDels <- function(readsWithInDels) {
  indelCoords <- data.table::data.table('qname' = as.vector(mcols(readsWithInDels)$name),
                                      'start' = start(readsWithInDels),
                                      'end' = end(readsWithInDels), 
                                      'type' = names(readsWithInDels))
  
  #for insertions, end should be the same as start
  indelCoords[type == 'I']$end <- indelCoords[type == 'I']$start
  
  indelCoords$ID <- paste(indelCoords$start, indelCoords$end, sep = ':')
  dt <- indelCoords[,length(qname), by = c('ID', 'start', 'end', 'type')]
  colnames(dt)[5] <- 'ReadSupport'
  dt$width <- dt$end - dt$start + 1
  dt <- dt[order(ReadSupport, decreasing = T)]
  return(dt)
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

printBedFile <- function(outDir, sampleName, df, tracktype) {
  outfile <- file.path(outDir, paste0(sampleName, '.', tracktype, ".bed"))
  writeLines(text = paste0("track name=\"",sampleName," top 100 ",tracktype," useScore=1"),
             con = outfile)
  write.table(x = df,
              file = outfile,
              quote = F, sep = '\t', col.names = F, row.names = F, append = T)
}

readsWithInDels <- printCoverageStats(bamFile, sampleName, outDir)

seqName <- seqnames(seqinfo(readsWithInDels))[1]

cutSites <- read.table(cutSitesFile, stringsAsFactors = F)

cutSiteStats <- do.call(rbind, lapply(1:nrow(cutSites), function(i) {
  x <- cutSites[i,]
  guide <- x[[1]]
  cutStart <- as.numeric(x[[2]])
  cutEnd <- cutStart + 1

  stats <- countEventsAtCutSite(seqName = seqName,
                             cutStart = cutStart,
                             cutEnd = cutEnd,
                             bamFile = bamFile,
                             readsWithInDels = readsWithInDels,
                             extend = 3)
  
  stats <- cbind(data.frame("sample" = sampleName,
                         "sgRNA" = guide,
                         "cutStart" = cutStart,
                         "cutEnd" = cutEnd), stats)
  
  #efficiency : number of indels that originate around the cut-site
  #             divided by the read coverage around the cut-site
  stats$indelEfficiency <- round(stats$indel / stats$cov * 100, 2)
  
  return(stats)
}))

write.table(x = cutSiteStats,
            file = file.path(outDir, paste0(sampleName, '.indel_stats_at_cutsites.tsv')),
            quote = F, sep = '\t', row.names = FALSE
            )

#TODO print summarized indels both as raw and as a summary to file.
#use summarizeIndels function 
indels <- summarizeInDels(readsWithInDels)
indels$seqname <- seqName

#print to BED file only top indels based on read support
#because it doesn't make sense to load every deletion  
topN <- 100 

#print summarized deletions to BED file
deletions <- indels[ type == 'D', c('seqname', 'start', 'end', 'ID', 'ReadSupport')]
if(nrow(deletions) > topN) {
  deletions <- deletions[1:100,]
}
printBedFile(outDir = outDir, sampleName = sampleName, 
             df = deletions, 
             tracktype = 'deletions')

#print summarized insertions to BED file
insertions <- indels[ type == 'I', c('seqname', 'start', 'end', 'ID', 'ReadSupport')]
if(nrow(insertions) > topN) {
  insertions <- insertions[1:100,]
}
printBedFile(outDir = outDir, sampleName = sampleName, 
             df = insertions, 
             tracktype = 'insertions')

#write all unfiltered indels to file 
dt <- cbind(as.data.table(readsWithInDels), as.data.table(mcols(readsWithInDels)))
dt$seqname <- seqName
dt$sample <- sampleName
dt <- dt[,c('seqname', 'sample', 'start', 'end', 'width', 'names', 'name')]
colnames(dt)[6:7] <- c('indelType', 'readID')

write.table(x = dt,
            file = file.path(outDir, paste0(sampleName, '.indels.unfiltered.tsv')),
            quote = F, sep = '\t', col.names = T, row.names = F, append = T)


