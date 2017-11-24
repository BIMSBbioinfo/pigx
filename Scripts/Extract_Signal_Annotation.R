# ---------------------------------------------------------------------------- #
#' Extract_Signal - given a bigWig file and a set of regions, extracts
#' the signal for downstream analysis
#'
#' @param annotation  - list - processed annotation as given by Prepare_Annotation.R
#' @param bed         - GRanges of peak regions
#' @param wig         - bigWig signal profile
#' @param outfile     - rds output file
#' @param expand.peak - regions extensions for the peaks
#' @param bin.num     - number of beans for binning
#' @param scriptdir   -  location of the R scripts and functions
#'
#' @return saves RDS object with a list of ScoreMatrix and summarized profiles
Extract_Signal = function(
    annotation  = NULL, 
    bed         = NULL,
    wig         = NULL,
    outfile     = './',
    expand.peak = 2000,
    bin.num     = 50,
    scriptdir   = NULL
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
        nams = names(annot$genomic_annotation)
        lsml = lapply(setNames(nams, nams), function(x)
            ScoreMatrixBin(wig, x[[x]], bin.num = bin.num, strand.aware=TRUE))
    
        bed_expand = resize(bed, width = expand.peak, fix='center')
        lsml$peak   = ScoreMatrixBin(wig, bed_expand, bin.num = bin.num)
   
    message('Sumarize Profiles ...')    
        profiles = lapply(lsml, function(x){
            d = data.table(sample = peakname, 
                           signal = rowMeans(x))
        })
    
    lout = list(
        sml      = lsml,
        profiles = profiles
    )
    saveRDS(lout, outfile)
}


# ---------------------------------------------------------------------------- #
extract_Signal(
    annotation  = snakemake@input[['annotation']],
    bed         = snakemake@input[['peaks']],
    wig         = snakemake@input[['wig']],
    outfile     = snakemake@output[['outfile']],
    expand.peak = snakemake@params[['expand_peak']],
    bin.num     = snakemake@params[['bin_num']],
    scriptdir   = snakemake@params[['scriptdir']],
    peakname    = snakemake@params[['peakname']]
)