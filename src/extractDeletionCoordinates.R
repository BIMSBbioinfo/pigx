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

pdf(file = file.path(outDir, paste0(sampleName,'.deletionSize_readSupport_distribution.pdf')))
scoreThreshold <- quantile(dt$ReadSupport, c(1:100)/100)[[95]]
print(ggplot(dt, aes(x = width, y = ReadSupport)) +
        geom_point(aes(color = ReadSupport > scoreThreshold)) +
        geom_hline(yintercept = scoreThreshold, color = 'red') +
        annotate("label",
                 x = max(dt$width)-100, y = scoreThreshold,
                 label = paste0('95th\npercentile =',scoreThreshold),
                 color = 'blue',
                 alpha = 0.4) +
        xlim(0, max(dt$width) + 100) +
        scale_y_log10() + theme(legend.position = 'bottom',
                                text = element_text(size = 16)) +
        labs(y = 'Read Support for Each Deletion\n (log10 scale)',
             x = 'Deletion Size',
             title = 'Distribution of Read Support per deletion sizes') +
        scale_color_brewer(palette = 'Set1'))
dev.off()

delCoords$score <- dt[match(delCoords$ID, dt$ID)]$ReadSupport
delCoords$seqname <-  as.character(unique(seqnames(aln)))

delCoords <- delCoords[order(score, decreasing = T)]

#convert to 0-based index
delCoords$start <- delCoords$start - 1

# write the coordinates of frequent deletions to a BED file
outfile <- file.path(outDir, paste0(sampleName, ".deletions.95th_percentile.BED"))
writeLines(text = paste0("track name=\"",sampleName," deletions (with read support > ",scoreThreshold,")\" useScore=1"),
           con = outfile)
write.table(x = unique(delCoords[,c(6,2:5)][score > scoreThreshold]),
            file = outfile,
            quote = F, sep = '\t', col.names = F, row.names = F, append = T)


setwd('/data/akalin/buyar/collaborations/jonathan/data/180409/annotation/')
s <- read.table('./sample_info.txt', header = T, sep = '\t')
s$sgRNAs <- gsub(', *', ':', s$sgRNAs)

samplesheet <- fread('../SampleSheet.bcl2fastq.csv', skip = 17)
colnames(samplesheet)[1:2] <- c('sample_id', 'sample_name')

readFiles <- dir(path = '../demultiplexed', pattern = '.gz$', full.names = T)

reads <- sapply(samplesheet$sample_name, function(x) {
  grep(pattern = paste0(x, '_S'), x = readFiles, value = T)
})

table(lengths(reads))

reads[lengths(reads) > 1]


ss <- merge(s[,c('sample_id', 'amplicon', 'sgRNAs')], samplesheet[,c('sample_id', 'sample_name')])

ss$reads <- unlist(reads[match(ss$sample_name, names(reads))])

ss$reads <- basename(ss$reads)

ss <- ss[,c(1,4,5,2,3)]
colnames(ss)[5] <- 'sgRNA_list'

write.csv(ss, file = './sample_sheet.csv', quote = FALSE, row.names = F)
