library(GenomicAlignments)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

bamFile <- args[1]
sampleName <- args[2]
outDir <- args[3]

cat(date(), "Running extractDeletionCoordinates.R with arguments:",args,"\n")

aln <- readGAlignments(bamFile, param = ScanBamParam(what="qname"))

delCoordsList <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(aln), ops = 'D',
                                                                   with.ops = T, pos = start(aln))
names(delCoordsList) <- mcols(aln)$qname

delCoords <- stack(delCoordsList)

delCoords <- data.table::data.table('qname' = as.vector(mcols(delCoords)$name),
                                    'start' = start(delCoords),
                                    'end' = end(delCoords))


delCoords$ID <- paste(delCoords$start, delCoords$end, sep = ':')
dt <- delCoords[,length(qname), by = c('ID')]
dt <- dt[order(V1, decreasing = T)]
dt$width <- sapply(strsplit(dt$ID, ':'), function(x) {
  x <- as.numeric(x)
  x[2] - x[1] + 1
})
colnames(dt)[2] <- 'ReadSupport'
dt <- rbind(dt, data.frame("ID" = 'NoDeletion', "ReadSupport" = sum(lengths(delCoordsList) == 0), 'width' = 0))
write.table(x = dt[,2],
            file = file.path(outDir, paste0(sampleName,'.deletion.counts.tsv')),
            append = F, quote = F, sep = '\t', row.names = dt$ID, col.names = sampleName)

scoreThreshold <- quantile(dt$ReadSupport, c(1:100)/100)[[95]]

delCoords$score <- dt[match(delCoords$ID, dt$ID)]$ReadSupport
delCoords$seqname <-  as.character(unique(seqnames(aln)))

delCoords <- delCoords[order(score, decreasing = T)]

#convert to 0-based index
delCoords$start <- delCoords$start - 1

# write the coordinates of deletions with at least scoreThreshold number of reads. 
scoreThreshold <- 10
outfile <- file.path(outDir, paste0(sampleName, ".deletions.bed"))
writeLines(text = paste0("track name=\"",sampleName," deletions (with read support > ",scoreThreshold,")\" useScore=1"),
           con = outfile)
write.table(x = unique(delCoords[,c(6,2:5)][score > scoreThreshold]),
            file = outfile,
            quote = F, sep = '\t', col.names = F, row.names = F, append = T)

