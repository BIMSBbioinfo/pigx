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

# R script takes as input a BAM file and exports a coverage track
# in bigwig format

args <- commandArgs(trailingOnly = TRUE)

countsFile <- args[1]
colDataFile <- args[2]
outDir <- args[3] 

counts <- read.table(countsFile, check.names = FALSE)
colData <- read.table(colDataFile)

common <- intersect(colnames(counts), rownames(colData)) 
if(length(common) < ncol(counts)) {
  stop("Not all samples in the counts table can be found in the 
       colData table")
}

# update counts and colData to respect the order of common samples
colData <- colData[common,]
counts <- counts[,common]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = colData, 
                                      design = ~1)

dds <- DESeq2::estimateSizeFactors(dds)

size_factor <- DESeq2::sizeFactors(dds)

normCountsFile <- file.path(outDir, "deseq_normalized_counts.tsv")
sizeFactorsFile <- file.path(outDir, "deseq_size_factors.txt")

write.table(x = as.data.frame(size_factor), file = sizeFactorsFile, 
            quote = FALSE, sep = '\t')
write.table(x = DESeq2::counts(dds, normalized = TRUE), file = normCountsFile, 
            quote = FALSE, sep = '\t')




