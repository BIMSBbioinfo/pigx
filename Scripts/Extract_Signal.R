# ---------------------------------------------------------------------------- #
extract_Signal = function(
    annotation  = NULL, 
    bed         = NULL,
    wig         = NULL,
    outfile     = './',
    expand.peak = 2000,
    bin.num     = 50,
    scriptdir
){
    
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
    
    message('Data ...')
        peaktype = ifelse(grepl('narrow', basename(bed)),'readNarrowPeak','readBroadPeak')
        bed      = match.fun(peaktype)(bed)
        
    message('Profiles ...')
        nams = names(annot$genomic_annotation)
        lsml = lapply(setNames(nams, nams), function(x)
            ScoreMatrixBin(wig, x[[x]], bin.num = bin.num))
    
        bed_expand = resize(bed, width = expand.peak, fix='center')
        lsml$peak   = ScoreMatrixBin(wig, bed_expand, bin.num = bin.num)
   
    saveRDS(lsml, outfile)
}


# ---------------------------------------------------------------------------- #
# extract_Signal(
#     annotation  = snakemake@input[['annotation']], 
#     bed         = snakemake@input[['peaks']],
#     wig         = snakemake@input[['wig']],
#     outfile     = snakemake@output[['outfile']],
#     expand.peak = snakemake@params[['peaks_width']],
#     bin.num     = snakemake@params[['bin_width']],
#     scriptdir   = snakemake@params[['scriptdir']]
# )