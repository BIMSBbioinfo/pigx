library(data.table)
library(GenomicAlignments)
library(ggplot2)

#read in output of samtools mpileup and parse_mpileup.py scripts and
#create a bedgraph file of per-base deletion scores
#also, print out some diagnostic plots
args = commandArgs(trailingOnly=TRUE)

bamFile <- args[1]
mpileupOutput <- args[2]
parseMpileupOutput <- args[3]
sampleName <- args[4]
outDir <- args[5]
cutSitesFile <- args[6] #path to sgRNA cutting sites for the target genome.
sgRNA_list <- args[7] # column (:) separated list of sgRNA ids (must match ids in cutSitesFile)
                      # that were used (or desired to be profiles) for the given sample.

cat("Running extractDeletionProfiles.R with arguments:",args,"\n")
#print bedgraph file of deletion ratios per base

printIndelStats <- function (mpileupOutput, parseMpileupOutput, outDir = getwd(), sampleName) {
  deletionScoresFile <- file.path(outDir, paste0(sampleName,'.deletionScores.bedgraph'))
  insertionScoresFile <- file.path(outDir, paste0(sampleName,'.insertionScores.bedgraph'))
  statsOutputFile <- file.path(outDir, paste0(sampleName,'.coverageStats.tsv'))
  
  #create a bedgraph file that contains deletion ratios for each base

  dt1 <- fread(input = paste0("cut -f 1-4 ",mpileupOutput), header = F)
  colnames(dt1) <- c('seqname', 'bp', 'base', 'cov')
  dt2 <- fread(parseMpileupOutput, select = c(1,6,7))
  dt3 <- merge(dt1[,c(2,4)], dt2, by = 'bp')
  dt3$sample <- sampleName
  dt3$indel <- dt3$del + dt3$ins

  #write bedgraph file for deletions
  trackDefinition <- paste0("track type=bedGraph name=",sampleName," deletion read support")
  writeLines(text = trackDefinition, con = deletionScoresFile)
  deletionScores <- ifelse(dt3$cov > 0, dt3$del/dt3$cov, 0)
  bg <- data.frame(cbind(dt1$seqname,
                         dt3$bp, dt3$bp, deletionScores), stringsAsFactors = F)
  #convert to 0-based index
  bg$V2 <- as.numeric(bg$V2) - 1
  #bedgraph file
  write.table(x = bg, file = deletionScoresFile, append = T,
              sep = '\t', quote = F, row.names = F, col.names = F)


  #write bedgraph file for insertions
  trackDefinition <- paste0("track type=bedGraph name=",sampleName," insertion read support")
  writeLines(text = trackDefinition, con = insertionScoresFile)
  insertionScores <- ifelse(dt3$cov > 0, dt3$ins/dt3$cov, 0)
  bg <- data.frame(cbind(dt1$seqname,
                         dt3$bp, dt3$bp, insertionScores), stringsAsFactors = F)
  #convert to 0-based index
  bg$V2 <- as.numeric(bg$V2) - 1
  #bedgraph file
  write.table(x = bg, file = insertionScoresFile, append = T,
              sep = '\t', quote = F, row.names = F, col.names = F)
  
  #coverage stats file
  dt3$delRatio <- ifelse(dt3$cov > 0, dt3$del/dt3$cov, 0)
  dt3$insRatio <- ifelse(dt3$cov > 0, dt3$ins/dt3$cov, 0)
  dt3$indelRatio <- ifelse(dt3$cov > 0, dt3$indel/dt3$cov, 0)
  
  write.table(x = dt3, file = statsOutputFile, append = T,
              sep = '\t', quote = F, row.names = F, col.names = T)
  return(dt3)
}

getReadsWithInDels <- function(bamFile) {
  aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname"))

  readsWithInDels <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(aln), 
                                                                              ops = c('I', 'D'),
                                                                     with.ops = T, pos = start(aln))
  names(readsWithInDels) <- mcols(aln)$qname

  readsWithInDels <- stack(readsWithInDels)
  seqinfo(readsWithInDels) <- seqinfo(aln)
  return(readsWithInDels)
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
#' @param extend (integer, default 5 bp) extend the searching area from cutting
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

readsWithInDels <- getReadsWithInDels(bamFile = bamFile)
seqName <- seqnames(seqinfo(readsWithInDels))[1]
inDels <- summarizeInDels(readsWithInDels)

indelStats <- printIndelStats(mpileupOutput, parseMpileupOutput, outDir, sampleName)

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
            quote = F, sep = '\t'
            )






