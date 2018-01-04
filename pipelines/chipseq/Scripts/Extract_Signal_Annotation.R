# ---------------------------------------------------------------------------- #
#' Extract_Signal_Annotation - given a bigWig file and a set of regions, extracts
#' the signal for downstream analysis
#'
#' @param annotation  - list - processed annotation as given by Prepare_Annotation.R
#' @param wig         - bigWig signal profile
#' @param outfile     - rds output file
#' @param scriptdir   -  location of the R scripts and functions
#'
#' @return saves RDS object with a list of ScoreMatrix and summarized profiles
Extract_Signal_Annotation = function(
    annotation  = NULL, 
    wig         = NULL,
    outfile     = './PEAKS.rds',
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
    
   
    if(is.null(wig))
        stop('wig file is not specified')
    
    # --------------------------------------------------------------- #
    message('Annotation ...')
        annot = readRDS(annotation)
        annot$genomic_annotation$gene = subset(annot$gtf$gtf, type='gene')
        
    message('Extract Profiles ...')
        nams = names(annot$genomic_annotation)
        lsml = lapply(setNames(nams, nams), function(x)
            ScoreMatrixBin(wig, x[[x]], bin.num = bin.num, strand.aware=TRUE))
    
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
Extract_Signal_Annotation(
    annotation  = snakemake@input[['annotation']],
    wig         = snakemake@input[['wig']],
    outfile     = snakemake@output[['outfile']],
    scriptdir   = snakemake@params[['scriptdir']],
    peakname    = snakemake@params[['peakname']]
)