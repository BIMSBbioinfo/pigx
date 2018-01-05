#export_bw.R - takes a methylraw RDS file and a tabulated list of chromosome lengths and outputs a bigwig file
# ---last updated Jan. 2018 by B. Osberg

#-------------------------------------------------------------------------

library(GenomicRanges)
library(stringr)
library(methylKit)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

RDS_filepath    <- args[1]
seqlengths_path <- args[2]
assembly        <- args[3]
out_path        <- args[4]

# ---------------------------------------------------

m1          = readRDS(file = RDS_filepath ) # import the methylRaw object from file.
seqdat_temp = read.csv(seqlengths_path, sep="\t", header=FALSE)
Sinfo <- Seqinfo(seqnames   = as.character(seqdat_temp[,1]),
                 seqlengths = seqdat_temp[,2],
                 genome     = assembly)


G1          = as(m1 , "GRanges")            # convert it to a GRanges object
seqinfo(G1) = Sinfo

G1$score = G1$numCs/G1$coverage

G1$coverage = NULL
G1$numCs    = NULL
G1$numTs    = NULL

export.bw( object = G1, con=out_path  )

# bigwig exported. Program complete.
