# PiGx RNAseq Pipeline.
#
# Copyright © 2017 Jonathan Ronen <yablee@gmail.com>
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

# R script takes outputs from STAR when counting readsPerGene and creates a reads matrix
# that is genes x samples

args <- commandArgs(trailingOnly = TRUE)

mapped_files_dir = args[1]
reads_per_gene_files = list.files(mapped_files_dir, 'ReadsPerGene')
sample_names = unlist(lapply(strsplit(reads_per_gene_files, "_Reads"), function(x) x[1]))

readspergene = lapply(reads_per_gene_files, function(f) read.table(file.path(mapped_files_dir, f), skip=4))
genes = readspergene[[1]]$V1

readsmatrix = matrix(nrow=length(genes), ncol=length(sample_names))
for(i in 1:length(sample_names)){
  readsmatrix[,i] = readspergene[[i]]$V2
}
colnames(readsmatrix) = sample_names
rownames(readsmatrix) = genes

write.table(readsmatrix, args[2], sep='\t')
