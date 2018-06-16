#!/usr/bin/env Rscript

# Author: EU
# Date: June, 2018
# This script takes as input a GTF file and returns annotation for genes

# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('get_GTF_genes')
# -------------------------------------------------------------------------- #
library(rtracklayer)
library(GenomicRanges)
get_GTF_genes <- function(gtfFile = NULL,
                          output = NULL) {
    if (is.null(gtfFile))
        stop("gtfFile not specified")
    if (is.null(output))
        stop("output file not specified")
    gtf <- import.gff2(gtfFile)
    g <- subset(gtf, type == "exon")
    gl <- split(g, g$gene_id)
    gl <- unlist(range(gl))
    
    export(gl, output, "gtf")
}

# -------------------------------------------------------------------------- #
get_GTF_genes(
    gtfFile = argv$input[['infile']],
    output = argv$output[['outfile']]
)