# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Signal_Annotation')

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
    outfile     = NULL,
    scriptdir   = NULL,
    bin.num     = 10
){

    # --------------------------------------------------------------- #
    # checks for default arugments
    deflist = as.list(formals(Extract_Signal_Annotation))
    arglist = as.list(match.call)
    arg.ind = names(deflist) %in% names(arglist)
    if(any(arg.ind))
        stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))

    # --------------------------------------------------------------- #
    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('stringr'))
    suppressPackageStartupMessages(library('tibble'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    # --------------------------------------------------------------- #
    message('Annotation ...')
        annot = readRDS(annotation)$full_annotation


    message('Extract Profiles ...')
        nams = names(annot)
        lsml = lapply(setNames(nams, nams), function(x)
            ScoreMatrixBin(wig, annot[[x]], bin.num = bin.num, strand.aware=TRUE))


    message('Summarize Profiles ...')
        sample_name = str_replace(basename(wig),'.bw','')
        profiles = lapply(names(lsml), function(x){
              dat = lsml[[x]]
              dm  = colMeans(dat)
              d = tibble(genomic_location = x,
                         sample_name      = sample_name,
                         position         = seq(dm),
                         value            = dm)
          })

    lout = list(
        lsml     = lsml,
        profiles = profiles
    )
    saveRDS(lout, outfile)
}


# ---------------------------------------------------------------------------- #
Extract_Signal_Annotation(
    annotation  = argv$input[['annotation']],
    wig         = argv$input[['wig']],
    outfile     = argv$output[['outfile']],
    scriptdir   = argv$params[['scriptdir']],
    bin.num     = argv$params[['number_of_bins']]
)
