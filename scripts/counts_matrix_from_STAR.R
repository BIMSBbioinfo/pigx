# R script takes outputs from STAR when counting readsPerGene and creates a reads matrix
# that is genes x samples

mapped_files_dir = snakemake@params[["mapped_files_dir"]]
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

write.table(readsmatrix, snakemake@output[[1]], sep='\t')
