# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Scripts/Argument_Parser.R'))
argv = Parse_Arguments('Make_UCSC_HUB')


# ---------------------------------------------------------------------------- #
#' ChIPQC
#'
#' @param config     - yaml formatted file used as snakemake config
#' @param outfile    - location of the output file
#' @param scriptdir  - location of R scripts
#'
#' @return saves RDS object with an annotated GRanges object
ChIPQC = function(
    config      = NULL,
    outfile     = NULL,
    scriptdir   = NULL,
    path_mapped = NULL,
    path_peaks  = NULL
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
    sample_sheet = data.frame(
        SampleID = names(config$peak_calling)
    )
    sample_sheet$bamReads = file.path(path_mapped)

    bamReads <- grep(pattern = "mrg1.AB",x = aln, value = TRUE)
    bamControl <- grep(pattern = "no.AB",x = aln, value = TRUE)
    Tissue <- gsub("_.*","",basename(bamReads))

    Factor <- rep("mrg1",6)
    Replicate <- rep(seq(1,3),2)
    ControlID <- paste0(Tissue,"c",seq(1,3))
    Peaks <- list.files(projectDir,pattern = "_peaks.narrowPeak$",recursive = TRUE,include.dirs = TRUE)
    PeakCaller <- rep("narrowPeak",6)

    sampleSheet <- data.frame(SampleID,Tissue,Factor,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)


    SampleID Tissue Factor Condition Treatment Replicate

    bamReads ControlID bamControl Peaks PeakCaller

    experiment = ChIPQC(sampleSheet,
                        consensus   = TRUE,
                        bCount      = TRUE,
                        summits     = 250,
                        chromosomes = NULL)

    saveRDS(experiment,file = outfile)
}



# ---------------------------------------------------------------------------- #
ChIPQC(
    config    = snakemake@input[['peaks']],
    outfile   = snakemake@output[['outfile']],
    scriptdir = snakemake@params[['scriptdir']]
)
