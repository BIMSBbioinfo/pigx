# ---------------------------------------------------------------------------- #
#' Annotate_Peaks
#'
#' @param annotation - list - processed annotation as given by Prepare_Annotation.R
#' @param peaks      - narrowPeaks/broadPeaks in macs2 output
#' @param outfile    - location of the output file
#' @param peakname   - Name of the peak file
#'
#' @return saves RDS object with an annotated GRanges object
Annotate_Peaks = function(
    annotation = NULL,
    peaks      = NULL,
    outfile    = NULL,
    peakname   = 'Peak',
){
    
    # --------------------------------------------------------------- #
    if(is.null(annotation))
        stop('Annotation is not defined')
    
    if(is.null(peaks))
        stop('Peaks are not defined')
    
    if(is.null(outfile))
        stop('output file path is not defined')
    
    # --------------------------------------------------------------- #
    library(data.table)
    library(stringr)
    
    # --------------------------------------------------------------- #
    message('Annotation ...')
        annot = readRDS(annotation)
        
    message('Reading peaks ...')
        peaktype = ifelse(grepl('narrow', basename(bed)),'readNarrowPeak','readBroadPeak')
        bed      = match.fun(peaktype)(bed)
    

    message('Annotating peaks ...')
        peaks$annot = Annotate_Ranges(bed, annot$genomic_annotation)
    
    message('CpGi ...')
        if(!is.null(annot$cpg)){
           
            peaks$cpgi = countOverlaps(peaks, annot$cpgi) > 0
        }else{
            peaks$cpgi = NA
        }
    
    message('Saving annotated peaks ...')
        saveRDS(peaks, outfile)
    
}



# ---------------------------------------------------------------------------- #
Annotate_Peaks(
    annotation  = snakemake@input[['annotation']],
    peaks       = snakemake@input[['peaks']],
    outfile     = snakemake@input[['outfile']],
    peakname    = snakemake@params[['name']]
)