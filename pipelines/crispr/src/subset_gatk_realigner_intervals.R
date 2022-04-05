# crispr-DART pipeline
#
# Copyright © 2017-2020 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the crispr-DART pipeline
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

library(GenomicRanges)
args = commandArgs(trailingOnly=TRUE)

sample_sheet_file <- args[1]
sample <- args[2]
all_intervals_file <- args[3]

sample_sheet <- read.csv(sample_sheet_file)
sample_sheet <- sample_sheet[sample_sheet$sample_name == sample,]

#subset all_intervals to only keep the ones that overlap target_region
all_intervals <- as(readLines(all_intervals_file), "GRanges")
target <- as(sample_sheet$target_region, "GRanges")

subset_intervals <- IRanges::subsetByOverlaps(all_intervals, target)

outfile <- gsub(".intervals$", '.target.intervals', all_intervals_file)

writeLines(text = paste(subset_intervals), con = outfile)
