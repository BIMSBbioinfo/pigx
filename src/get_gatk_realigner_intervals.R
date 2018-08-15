library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

fastaFile <- args[1]
outFile <- args[2]

seq <- Biostrings::readDNAStringSet(fastaFile)

#print start and end location of the amplicon in the format: e.g. "dpy-10:1-1320".
write(paste0(names(seq), ":1-", width(seq)), file = outFile)