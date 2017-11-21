# ---------------------------------------------------------------------------- #
prepare_Annotation = function(
    annotation=NULL, 
    scriptdir, 
    outpath
){
    
    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)
    
    print('Annotation...')
    if(is.null(annotation))
        stop('Annotation is not specified')
    
    annotation         = read_Annotation(annotation)
    
    genomic_annotation = GTFGetAnnotation(annotation$gtf$gtf)
   
    # --------------------------------------------------------------- #
    message('RDS ...')
    ldat = list(annotation         = annotation, 
                genomic_annotation = genomic_annotation)
    saveRDS(ldat, outpath)
}


# ---------------------------------------------------------------------------- #
# function call
# prepare_Annotation(
#     annotation  = snakemake@params[['annotation']],
#     outpath     = snakemake@params[['outpath']],
#     scriptdir   = snakemake@params[['scriptdir']]
# )