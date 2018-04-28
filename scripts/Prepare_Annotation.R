# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Prepare_Annotation')


# ---------------------------------------------------------------------------- #
Prepare_Annotation = function(
    gtf_path     = NULL,
    outfile      = NULL,
    scriptdir    = './scripts',
    width_params
){

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    message('Annotation...')
    if(is.null(gtf_path))
        stop('Annotation is not specified')

    annotation = ReadGTFAnnotation(gtf_path)

    # selects the genomic regions used for the profiles
    annotation$full_annotation    = GTFGetAnnotation(annotation$gtf,
                                                     width_params=width_params)

    # selects the genomic regions used for annotation
    annotation$genomic_annotation = annotation$full_annotation[c('tss','tts','exon','intron')]


    # --------------------------------------------------------------- #
    message('RDS ...')
    saveRDS(annotation, outfile)
}


# ---------------------------------------------------------------------------- #
# function call
Prepare_Annotation(
    gtf_path     = argv$input[['gtf_path']],
    outfile      = argv$output[['outfile']],
    scriptdir    = argv$params[['scriptdir']],
    width_params = argv$params[['width_params']]
)
