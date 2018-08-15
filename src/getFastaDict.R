library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

fastaFile <- args[1]
outFile <- args[2]

seq <- Biostrings::readDNAStringSet(fastaFile)

# print a header line in the style of sam file specification
# the line should contain the length of the amplicon and its name
# e.g. "@SQ	SN:dpy-10	LN:1320"
write(x = paste0("@SQ\t", 
                 "SN:", names(seq), "\t",
                 "LN:", width(seq)), 
      file = outFile)