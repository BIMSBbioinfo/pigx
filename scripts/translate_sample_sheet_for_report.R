# R script takes sample_sheet.csv as input and creates a colData.tsv file

args <- commandArgs(trailingOnly = TRUE)

sample_sheet = args[1]
s = read.csv(file = sample_sheet)
rownames(s) = s$name
s$group = s$sample_type
s = s[colnames(s)[-grep("name|reads", colnames(s))]]
write.table(s, "colData.tsv", sep="\t")

