# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Signal_Peaks')

# ---------------------------------------------------------------------------- #
#' Extract_Signal_Peaks - given a bigWig file and a set of regions, extracts
#' the signal around the peak regions
#'
#' @param bed         - GRanges of peak regions
#' @param wig         - bigWig signal profile
#' @param outfile     - rds output file
#' @param expand.peak - regions extensions for the peaks
#' @param bin.num     - number of beans for binning
#' @param scriptdir   -  location of the R scripts and functions
#'
#' @return saves RDS object with a list of ScoreMatrix and summarized profiles
Extract_Signal_Peaks = function(
    bed         = NULL,
    wig         = NULL,
    outfile     = './',
    expand.peak = 2000,
    bin.num     = 50,
    scriptdir   = NULL,
    peakname    = 'Peak'
){

    if(is.null(scriptdir))
        stop('Please specify the script directory')

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)
    # --------------------------------------------------------------- #
    if(is.null(annotation))
        stop('Annotation is not specified')

    if(is.null(bed))
        stop('bed file is not specified')

    if(is.null(wig))
        stop('wig file is not specified')

    # --------------------------------------------------------------- #
    message('Annotation ...')
        annot = readRDS(annotation)
        annot$genomic_annotation$gene = subset(annot$gtf$gtf, type='gene')

    message('Data ...')
        peaktype = ifelse(grepl('narrow', basename(bed)),'readNarrowPeak','readBroadPeak')
        bed      = match.fun(peaktype)(bed)

    message('Extract Profiles ...')
        bed_expand = resize(bed, width = expand.peak, fix='center')
        sml  = ScoreMatrixBin(wig, bed_expand, bin.num = bin.num)

    message('Sumarize Profiles ...')
        profiles = data.table(
            sample  = peakname,
            signal = rowMeans(sml))


    lout = list(
        sml      = ssml,
        profiles = profiles
    )
    saveRDS(lout, outfile)
}


# ---------------------------------------------------------------------------- #
Extract_Signal_Peaks(
    bed         = argv$input[['peaks']],
    wig         = argv$input[['wig']],
    outfile     = argv$output[['outfile']],
    expand.peak = argv$params[['expand_peak']],
    bin.num     = argv$params[['bin_num']],
    scriptdir   = argv$params[['scriptdir']],
    peakname    = argv$params[['peakname']]
)
