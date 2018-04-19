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

printDeletionProfiles <- function (mpileupOutput, parseMpileupOutput, outDir = getwd(), sampleName) {
  bedgraphOutputFile <- file.path(outDir, paste0(sampleName,'.deletionScores.bedgraph'))
  #create a bedgraph file that contains deletion ratios for each base

  dt1 <- fread(input = paste0("cut -f 1-4 ",mpileupOutput), header = F)
  colnames(dt1) <- c('seqname', 'bp', 'base', 'cov')
  dt2 <- fread(parseMpileupOutput, select = c(1,6))
  dt3 <- merge(dt1[,c(2,4)], dt2, by = 'bp')
  dt3$sample <- sampleName

  #write bedgraph file
  trackDefinition <- paste0("track type=bedGraph name=",sampleName," deletion read support")
  writeLines(text = trackDefinition, con = bedgraphOutputFile)
  deletionScores <- ifelse(dt3$cov > 0, dt3$del/dt3$cov, 0)
  bg <- data.frame(cbind(dt1$seqname,
                         dt3$bp, dt3$bp, deletionScores), stringsAsFactors = F)
  #convert to 0-based index
  bg$V2 <- as.numeric(bg$V2) - 1
  #bedgraph file
  write.table(x = bg, file = bedgraphOutputFile, append = T,
              sep = '\t', quote = F, row.names = F, col.names = F)

  dt3$delRatio <- ifelse(dt3$cov > 0, dt3$del/dt3$cov, 0)

  return(dt3)
}

getReadsWithDeletions <- function(bamFile) {
  aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname"))

  readsWithDeletionsList <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(aln), ops = 'D',
                                                                     with.ops = T, pos = start(aln))
  names(readsWithDeletionsList) <- mcols(aln)$qname

  readsWithDeletions <- stack(readsWithDeletionsList)
  seqinfo(readsWithDeletions) <- seqinfo(aln)
  return(readsWithDeletions)
}

summarizeDeletions <- function(readsWithDeletions) {
  delCoords <- data.table::data.table('qname' = as.vector(mcols(readsWithDeletions)$name),
                                      'start' = start(readsWithDeletions),
                                      'end' = end(readsWithDeletions))

  delCoords$ID <- paste(delCoords$start, delCoords$end, sep = ':')
  dt <- delCoords[,length(qname), by = c('ID', 'start', 'end')]
  colnames(dt)[4] <- 'ReadSupport'
  dt$width <- dt$end - dt$start + 1
  dt <- dt[order(ReadSupport, decreasing = T)]
  return(dt)
}

#' @param cutStart expected cutting site start pos for the sgRNA
#' @param cutEnd expected cutting site end pos for the sgRNA
#' @param bamFile path to bam file
#' @param extend (integer, default 5 bp) extend the searching area from cutting
#'   sites in either direction
calculateGuideCuttingEfficiency <- function(seqName, cutStart, cutEnd, bamFile, readsWithDeletions, extend = 3) {
  if(cutEnd < cutStart) {
    stop("End position of cutting site must be larger than start position\n")
  }
  if(cutStart < 0 | cutEnd < 0) {
    stop("Start/End positions of cutting sites must be positive values")
  }

  cutStartExt <- cutStart - extend
  cutEndExt <- cutEnd + extend

  # find the number of reads with deletions that start or end within the
  # extended region of the cut-site
  deletionCount <- length(readsWithDeletions[
    (start(readsWithDeletions) >= cutStartExt & start(readsWithDeletions) <= cutEndExt) |
      (end(readsWithDeletions) >= cutStartExt & end(readsWithDeletions) <= cutEndExt),])

  # find the total number of reads whose alignments overlap the extended region of the cut-site
  aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname",
                                                       which = GRanges(seqnames = seqName,
                                                                       ranges = IRanges(start = cutStartExt,
                                                                                        end = cutEndExt),
                                                                       strand = '*')))

  #efficiency : number of deletions that originate around the cut-site
  #             divided by the read coverage around the cut-site

  efficiency <- round(deletionCount/length(aln) * 100, 2)

  return(efficiency)
}

plotDeletionSizeDistributionAtCutSites <- function(sampleName, guide, cutStart, cutEnd, deletions, extend = 3, maxDelSize = 100) {
  if(cutEnd < cutStart) {
    stop("End position of cutting site must be larger than start position\n")
  }
  if(cutStart < 0 | cutEnd < 0) {
    stop("Start/End positions of cutting sites must be positive values")
  }

  cutStartExt <- cutStart - extend
  cutEndExt <- cutEnd + extend

  cutSiteDeletions <- readsWithDeletions[
    (start(readsWithDeletions) >= cutStartExt & start(readsWithDeletions) <= cutEndExt) |
      (end(readsWithDeletions) >= cutStartExt & end(readsWithDeletions) <= cutEndExt),]

  p <- ggplot(data.frame("width" = width(cutSiteDeletions)), aes(x = width)) +
    geom_histogram(aes(fill = width %% 3 == 0), binwidth = 1, show.legend = F) + xlim(c(0, maxDelSize+1)) +
    labs(title = paste0('Deletion size distribution (for deletions <= ',maxDelSize,'bp)'),
         subtitle = paste0("Sample:",sampleName, "; sgRNA:",guide,"; coords:",cutStart,"-",cutEnd),
         x = 'Deletion Size (bp)', y = "Count") +
    theme_bw(base_size = 14)

  return(p)
}

readsWithDeletions <- getReadsWithDeletions(bamFile = bamFile)
seqName <- seqnames(seqinfo(readsWithDeletions))[1]
deletions <- summarizeDeletions(readsWithDeletions)

deletionProfile <- printDeletionProfiles(mpileupOutput, parseMpileupOutput, outDir, sampleName)

cutSites <- read.table(cutSitesFile, stringsAsFactors = F)

cutEfficiencies <- do.call(rbind, lapply(1:nrow(cutSites), function(i) {
  x <- cutSites[i,]
  guide <- x[[1]]
  cutStart <- as.numeric(x[[2]])
  cutEnd <- cutStart + 1
  cutEfficiency <- calculateGuideCuttingEfficiency(seqName = seqName,
                                                   cutStart = cutStart,
                                                   cutEnd = cutEnd,
                                                   bamFile = bamFile,
                                                   readsWithDeletions = readsWithDeletions,
                                                   extend = 3)
  return(data.frame("sample" = sampleName,
                    "sgRNA" = guide,
                    "cutStart" = cutStart,
                    "cutEnd" = cutEnd,
                    "cutEfficiency" = cutEfficiency))
}))

write.table(x = cutEfficiencies,
            file = file.path(outDir, paste0(sampleName, '.sgRNA.efficiency.tsv')),
            quote = F, sep = '\t'
            )

#for the guides that were used for the given sample, draw deletion size histograms
# guides <- gsub(' ', '', unlist(strsplit(x = sgRNA_list, split = ':')))
# 
# pdf(file = file.path(outDir, paste0(sampleName, ".deletionSizeHistogram.pdf")))
# for(guide in guides) {
#   x <- cutSites[cutSites$V1 == guide,]
#   cutStart <- as.numeric(x[[2]])
#   cutEnd <- cutStart + 1
#   p1 <- plotDeletionSizeDistributionAtCutSites(sampleName = sampleName,
#                                          guide = guide,
#                                          cutStart = cutStart,
#                                          cutEnd = cutEnd,
#                                          deletions = deletions,
#                                          maxDelSize = 100,
#                                          extend = 3)
#   p2 <- plotDeletionSizeDistributionAtCutSites(sampleName = sampleName,
#                                                guide = guide,
#                                                cutStart = cutStart,
#                                                cutEnd = cutEnd,
#                                                deletions = deletions,
#                                                maxDelSize = 10,
#                                                extend = 3)
# 
#   print(p1)
#   print(p2)
# }
# dev.off()




