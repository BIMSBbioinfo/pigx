# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/scripts/Argument_Parser.R'))
argv = Parse_Arguments('Prepare_Annotation')

# ---------------------------------------------------------------------------- #
Prepare_Annotation = function(
    gtf_path=NULL,
    outfile,
    scriptdir
){

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    message('Annotation...')
    if(is.null(gtf_path))
        stop('Annotation is not specified')

    annotation = ReadGTFAnnotation(gtf_path)

    annotation$genomic_annotation = GTFGetAnnotation(annotation$gtf)

    # --------------------------------------------------------------- #
    message('RDS ...')
    saveRDS(annotation, outfile)
}


# ---------------------------------------------------------------------------- #
# function call
Prepare_Annotation(
    gtf_path  = argv$input[['gtf_path']],
    outfile   = argv$output[['outfile']],
    scriptdir = argv$params[['scriptdir']]
)
