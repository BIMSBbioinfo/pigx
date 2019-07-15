# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')

# ---------------------------------------------------------------------------- #
GRangesToBigWig = function(g, outfile, scale.fac){
    cov = coverage(g)
    cov = round(cov*scale.fac,2)
    rtracklayer::export.bw(cov, outfile)
}

# ---------------------------------------------------------------------------- #
BamToBigWig = function(
    bamfile            = NULL,
    bwfile             = NULL

){
    if(is.null(bamfile))
        stop('bamfile not specified')

    if(is.null(bwfile))
        stop('bwfile not specified')


    suppressPackageStartupMessages({
        library(Rsamtools)
        library(GenomicRanges)
        library(GenomicFiles)
        library(GenomicAlignments)
        library(yaml)
        library(dplyr)
        library(stringr)
    })


    # ---------------------------------------------------------------------------- #
    YIELD = function(bfl) {
        param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
        aln = readGAlignments(bfl, param=param)
        g   = unlist(grglist(aln, use.mcols=TRUE))
    }

    bf = BamFile(bamfile, yieldSize=1000000)
    result = reduceByYield(bf, YIELD, MAP=identity, REDUCE=c)

    if(length(result) == 0)
        stop('All the cells were filtered')

    scale.fac = round(1e8/length(result),3)

    g       = result[strand(result) == '-']
    outfile = str_replace(bwfile,'bw','m.bw')
    GRangesToBigWig(g, outfile, scale.fac)

    g       = result[strand(result) == '+']
    outfile = str_replace(bwfile,'bw','p.bw')
    GRangesToBigWig(g, outfile, scale.fac)

    GRangesToBigWig(result, bwfile, scale.fac)
}

# -------------------------------------------------------------------------- #
BamToBigWig(
    bamfile            = argv$input[['bamfile']],
    bwfile             = argv$output[['bwfile']]
)
