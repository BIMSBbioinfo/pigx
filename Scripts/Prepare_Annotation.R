# ---------------------------------------------------------------------------- #
Prepare_Annotation = function(
    annotation=NULL, 
    scriptdir, 
    outfile
){
    
    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)
    
    print('Annotation...')
    if(is.null(annotation))
        stop('Annotation is not specified')
    
    annotation = read_Annotation(annotation)
    
    annotation$genomic_annotation = GTFGetAnnotation(annotation$gtf$gtf)
   
    # --------------------------------------------------------------- #
    message('RDS ...')
    saveRDS(annotation, outfile)
}


# ---------------------------------------------------------------------------- #
# function call
prepare_Annotation(
    annotation  = snakemake@input[['annotation']],
    outfile     = snakemake@output[['outfile']],
    scriptdir   = snakemake@params[['scriptdir']]
)