#!/usr/bin/env Rscript

# Author: EU
# Date: June, 2018
# This script takes as input a GTF file and returns annotation for genes

# 1. Collect Arguments ----------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

## Default setting when no arguments passed
if(length(args) == 0) {
    args = c("--help")
}
help_command = "
get_GTF_genes.R: Takes as input a GTF file and returns annotation for genes
Arguments:
--gtfFile Path to the GTF format file containing the genome annotation
--output Path to the output GTF file containing annotation for genes
Example:
Rscript get_GTF_genes.R --gtfFile=<path to genome.gtf> --output=<path to output.gtf>\n"

## Help section
if("--help" %in% args) {
    cat(help_command, "\n")
    q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs = function(x) {
    myArgs = unlist(strsplit(x, "--"))
    myArgs = myArgs[myArgs != '']
    #when no values are provided for the argument
    myArgs = gsub(pattern = "=$", replacement = "= ", x = myArgs)
    myArgs = as.data.frame(do.call(rbind, strsplit(myArgs, "=")))
    myArgs$V2 = gsub(' ', '', myArgs$V2)
    return(myArgs)
}

argsDF = parseArgs(args)
argsL = as.list(as.character(argsDF$V2))
names(argsL) = argsDF$V1

if(!("gtfFile" %in% argsDF$V1)) {
    cat(help_command, "\n")
    stop("Missing argument: gtfFile. Provide the path to .gtf file 
         containing the genome annotation.")
}

if(!("output" %in% argsDF$V1)) {
    cat(help_command, "\n")
    stop("Missing argument: output Provide the path to output .gtf file 
         that will contain annotation for genes.")
}

gtfFile = argsL$gtfFile
output = argsL$output

# 2. Get genes annotation ----------------------------------------------
library(rtracklayer)
library(GenomicRanges)

gtf <- import.gff2(gtfFile)
g <- subset(gtf, type == "exon")
gl <- split(g, g$gene_id)
gl <- unlist(range(gl))

export(gl, output, "gtf")
