# ---------------------------------------------------------------------------- #
Annotate_Peaks = function(
    annotation = NULL,
    peaks      = NULL,
    outfile    = NULL
){
    
    if(is.null(annotation))
        stop('Annotation is not defined')
    
    if(is.null(peaks))
        stop('Peaks are not defined')
    
    if(is.null(outfile))
        stop('output file path is not defined')
    
    # --------------------------------------------------------------- #
    message('Annotation ...')
        annot = readRDS(annotation)
        
    message('Reading peaks ...')
        peaktype = ifelse(grepl('narrow', basename(bed)),'readNarrowPeak','readBroadPeak')
        bed      = match.fun(peaktype)(bed)
    
        
    ### FINISH THIS STUFF WITH PROPER FUNCTIONS
        # - add the function
        # - add conversion to data table
    message('Annotating peaks ...')
       peaks_annot = Annotate_Ranges(bed, annot$genomic_annotation)
    
    
    message('Saving annotated peaks ...')
        saveRDS(peaks_annot, outfile)
    
}



# ---------------------------------------------------------------------------- #
Annotate_Peaks(
    annotation  = snakemake@input[['annotation']],
    peaks       = snakemake@input[['peaks']],
    outfile     = snakemake@input[['outfile']]
)