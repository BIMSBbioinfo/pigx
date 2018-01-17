# PiGx RNAseq Pipeline.
#
# Copyright © 2017 Bora Uyar <bora.uyar@mdc-berlin.de>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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

# R script takes sample_sheet.csv as input and creates a colData.tsv file

args <- commandArgs(trailingOnly = TRUE)

sample_sheet = args[1]
s = read.csv(file = sample_sheet)
rownames(s) = s$name
s$group = s$sample_type
s = s[colnames(s)[-grep("name|reads", colnames(s))]]
write.table(s, "colData.tsv", sep="\t")

