# PiGx RNAseq Pipeline.
#
# Copyright © 2019 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the PiGx RNAseq Pipeline.
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

# The script is used for counting number of reads per feature using 
# GenomicAlignments::summarizeOverlaps
# see reference workflow: 
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#read-counting-step

args <- commandArgs(trailingOnly = TRUE)

sampleName <- args[1]
bamFile <- args[2] 
gtfFile <- args[3]
singleEnd <- as(args[4], "logical") #whether reads are single/paired end
counting_mode <- args[5] #see "mode" in summarizeOverlaps
count_nonunique <- as(args[6], "logical") #see "inter.feature" in summarizeOverlaps
strandedness <- args[7] #unspecific, forward, reverse
feature <- args[8] #which feature to count
group_feature_by <- args[9] #group features by "transcript_id" or 'gene_id'? 
yieldSize <- as.numeric(args[10]) # number of reads to process at a time

require(GenomicAlignments)

# define BAM file connection 
bamfile <- Rsamtools::BamFile(file = bamFile,
                              asMates = !singleEnd,  
                              yieldSize = 2000000)

# check counting mode
if(!counting_mode %in% c('Union', 'IntersectionStrict', 
                         'IntersectionNotEmpty')) {
  stop("Counting mode must be one of 'Union', 'IntersectionStrict', 
       or 'IntersectionNotEmpty'")
}

# check count_nonunique argument
if(class(count_nonunique) != 'logical') {
  stop("count_nonunique must be set to TRUE or FALSE")
}

#check strandedness
if(!strandedness %in% c('unspecific', 'forward', 'reverse')) {
  stop("Acceptable strandedness options are 'unspecific, forward, 
       or reverse.")
}
strand_ignore <- FALSE
if(strandedness == 'unspecific') {
  strand_ignore <- TRUE
}

# get annotations
gtfData <- rtracklayer::import.gff(gtfFile, format = 'gtf')

feat <- gtfData[gtfData$type == feature,]

if(length(feat) == 0) {
  stop("GTF file doesn't seem to contain such a feature:",feature,
       "in the column 'type'")
}

if(!group_feature_by %in% colnames(mcols(feat))) {
  stop("The chose feature grouping factor (",group_feature_by,")",
      "for grouping features doesn't seem to 
      exist in the GTF file")
}

# group features by 
feat <- split(feat, mcols(feat)[group_feature_by][1])

if(strandedness == 'reverse') {
  se <- GenomicAlignments::summarizeOverlaps(features = feat, 
                          reads = bamfile,
                          mode = counting_mode,
                          singleEnd = singleEnd,
                          ignore.strand = strand_ignore, 
                          inter.feature = count_nonunique, 
                          preprocess.reads = invertStrand)
} else {
  se <- GenomicAlignments::summarizeOverlaps(features = feat, 
                          reads = bamfile,
                          mode = counting_mode,
                          singleEnd = singleEnd,
                          ignore.strand = strand_ignore, 
                          inter.feature = count_nonunique)
}

counts <- assay(se, 'counts')

colnames(counts) <- sampleName

outFile <- file.path(dirname(bamFile), 
                     paste0(sampleName, ".read_counts.csv"))

write.csv(counts, outFile, quote = FALSE)
