# PiGx BSseq Pipeline.
#
# Copyright © 2018 Bren Osberg <brendan.osberg@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
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

#export_bw.R - takes a methylraw RDS file and a tabulated list of chromosome lengths and outputs a bigwig file
# ---last updated Jan. 2018 by B. Osberg

#-------------------------------------------------------------------------

suppressPackageStartupMessages(expr = {
  library(GenomicRanges)
  library(stringr)
  library(methylKit)
  library(rtracklayer)
})

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


G1            <- as(m1 , "GRanges")            # convert it to a GRanges object

seqlevels(G1) <- seqlevels(Sinfo)              # ensure the full set of seqnames 
                                               # from the ref-genome are included
                                               # (even if this data set is low-
                                               # coverage and missing chrom's)
seqinfo(G1)   <- Sinfo
G1$score = G1$numCs/G1$coverage

G1$coverage = NULL
G1$numCs    = NULL
G1$numTs    = NULL

export.bw( object = G1, con=out_path  )

# bigwig exported. Program complete.
