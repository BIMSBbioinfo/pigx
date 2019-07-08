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
# in bigwig format. Uses DESeq2 estimated size factors to normalize
# the coverage tracks by size factors computed across all 
# samples available in the sample sheet. 
# (see DESeq2::estimateSizeFactors)

args <- commandArgs(trailingOnly = TRUE)

bamFile <- args[1] #STAR alignment file 
sampleName <- args[2] 
size_factors_file <- args[3] #deseq size factors for all samples
outDir <- args[4] # where to write the bigwig files

#' @param cov RLE-list object 
#' @param size_factor single numeric value corresponding 
#' to the size factor for the sample 
#' (see DESeq2::estimateSizeFactors)
scale_coverage <- function(cov, size_factor) {
  cov_scaled <- lapply(cov, function(x) {
    S4Vectors::runValue(x) <- round(S4Vectors::runValue(x) / size_factor, 1)
    return(x)
  })
  return(as(cov_scaled, "SimpleRleList"))
}

message(date()," ... Reading alignments from bam file: \n",bamFile,"\n")
aln <- GenomicAlignments::readGAlignments(file = bamFile, index = bamFile)

size_factors <- read.table(size_factors_file)

message(date()," ... Getting strand-specific coverage data")
cov_pos <- scale_coverage(cov = GenomicRanges::coverage(aln[GenomicRanges::strand(aln) == '+',]), 
                          size_factor = size_factors[sampleName,])
cov_neg <- scale_coverage(cov = GenomicRanges::coverage(aln[GenomicRanges::strand(aln) == '-',]), 
                          size_factor = size_factors[sampleName,])

out_pos <- file.path(outDir, paste0(sampleName, ".forward.bigwig"))
out_neg <- file.path(outDir, paste0(sampleName, ".reverse.bigwig"))

message(date(), " ... exporting bigwig files")

rtracklayer::export.bw(cov_pos, con = out_pos, format = 'bw')
rtracklayer::export.bw(cov_neg, con = out_neg, format = 'bw')




