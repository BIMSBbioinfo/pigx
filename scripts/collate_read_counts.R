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

# Collate read counts into one big matrix 

args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
out_file <- args[2]

count_files <- dir(input_dir, pattern = ".read_counts.csv$", full.names = TRUE)

# get read counts for each sample into a list
counts <- lapply(count_files, function(f) {
  data.table::fread(f)
})

# merge list of data frames
counts_all <- as.data.frame(Reduce(function(dtf1, dtf2) 
  merge(dtf1, dtf2, by = "V1", all.x = TRUE),
       counts))
rownames(counts_all) <- counts_all$V1
counts_all$V1 <- NULL

# save results to out file
write.table(counts_all, out_file, quote = FALSE, 
            sep = '\t')









