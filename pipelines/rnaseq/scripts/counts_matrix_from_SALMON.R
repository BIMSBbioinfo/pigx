# R script takes outputs from salmon and uses tximport package to create a matrix
# that is genes x samples

args <- commandArgs(trailingOnly = TRUE)
salmon_output_folder <- args[1]
colDataFile <- args[2]

colData <- read.table(colDataFile)
files <- file.path(salmon_output_folder, rownames(colData), "quant.sf")
names(files) <- rownames(colData)
txi <- tximport::tximport(files, type = "salmon", txOut = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, colData, ~group)

write.table(x = DESeq2::counts(dds), 
            file = file.path(salmon_output_folder, "counts_from_SALMON.tsv"), 
            quote = F, sep = '\t')
