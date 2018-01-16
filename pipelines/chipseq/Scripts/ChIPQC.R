# ---------------------------------------------------------------------------- #
#' ChIPQC
#'
#' @param annotation - list - processed annotation as given by Prepare_Annotation.R
#' @param peaks      - narrowPeaks/broadPeaks in macs2 output
#' @param outfile    - location of the output file
#' @param peakname   - Name of the peak file
#'
#' @return saves RDS object with an annotated GRanges object
ChIPQC = function(
    bamfiles   = NULL,
    peaks_path = NULL,
    outfile    = NULL,
    peakname   = NULL,
    scriptdir  = NULL
){
    
    # --------------------------------------------------------------- #
    # checks for default arugments
    deflist = as.list(formals(Annotate_Peaks))
    arglist = as.list(match.call)
    arg.ind = names(deflist) %in% names(arglist)
    if(any(arg.ind))
        stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))
    
    # --------------------------------------------------------------- #
    library(ChIPQC)
    library(stringr)
    source(file.path(scriptdir, 'Functions_Helper.R'), local=TRUE)
    
    # --------------------------------------------------------------- #
    bamReads <- grep(pattern = "mrg1.AB",x = aln, value = TRUE)
    bamControl <- grep(pattern = "no.AB",x = aln, value = TRUE)
    Tissue <- gsub("_.*","",basename(bamReads))
    SampleID <- paste0(Tissue,seq(1,3))
    Factor <- rep("mrg1",6)
    Replicate <- rep(seq(1,3),2)
    ControlID <- paste0(Tissue,"c",seq(1,3))
    Peaks <- list.files(projectDir,pattern = "_peaks.narrowPeak$",recursive = TRUE,include.dirs = TRUE)
    PeakCaller <- rep("narrowPeak",6)
    
    sampleSheet <- data.frame(SampleID,Tissue,Factor,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)
    
    
    experiment = ChIPQC(sampleSheet,
                        consensus   = TRUE,
                        bCount      = TRUE,
                        summits     = 250,
                        chromosomes = NULL)
    
    saveRDS(experiment,file = outfile)
}



# ---------------------------------------------------------------------------- #
ChIPQC(
    annotation  = snakemake@input[['annotation']],
    peaks_path  = snakemake@input[['peaks']],
    outfile     = snakemake@output[['outfile']],
    peakname    = snakemake@params[['name']],
    scriptdir   = snakemake@params[['scriptdir']]
)


annotation = '/home/vfranke/Tmp/pigxtest/Annotation/Processed_Annotation.rds'
peaks      = '/home/vfranke/Tmp/pigxtest/Peaks/MACS2/Peaks1/Peaks1_qsort.bed'
outfile    = '~/Tmp/file.rds'
peakname   = 'peak1'
scriptdir  = '/home/vfranke/Projects/AAkalin_PIX/Scripts'








