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

bamFile <- args[1]
outFile <- args[2] 

aln <- GenomicAlignments::readGAlignments(bamFile)

cov <- GenomicRanges::coverage(aln)

rtracklayer::export.bw(object = cov, con = outFile)





