# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extend_Regions')

# ---------------------------------------------------------------------------- #
Extend_Regions = function(
    inpath,
    outpath,
    extend       = NULL,
    scale_index  = FALSE,
    scale_factor = 1e6,
    chunk_size   = 1e7,
    library      = 'single'
){

    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(GenomicAlignments)
        library(rtracklayer)
        library(Rsamtools)
    })
    if(!is.numeric(extend))
        stop('Extend_Regions: extend parameter needs to be a number')

    if(!file.exists(inpath))
        stop('Extend_Regions: input file does not exist')
    
    if(library == 'paired'){
        param   = ScanBamParam(flag=scanBamFlag(isPaired = TRUE, isProperPair=TRUE))
    }else{
        param   = ScanBamParam(flag=scanBamFlag(isPaired = FALSE))
    }
    bamfile = BamFile(inpath, yieldSize=chunk_size)
    open(bamfile)
    lcov = NULL
    total = 0
    while (length(chunk <- readGAlignments(bamfile, param=param))) {
        gchunk = resize(granges(chunk), width=extend)
        cov = coverage(gchunk)
        if(is.null(lcov)){
            lcov = cov   
        }else{
            lcov = lcov + cov
        }
        total = total+length(gchunk)
    }
    close(bamfile)
    
    if(scale_index == 'yes' || scale_index == TRUE)
        lcov = round(lcov * (scale_factor/total),2)

    export.bw(lcov, outpath)
}

# ---------------------------------------------------------------------------- #
Extend_Regions(
  inpath        = argv$input[['file']],
  outpath       = argv$output[['outfile']],
  extend        = argv$params[['extend']],
  scale_index   = argv$params[['scale']],
  library       = argv$params[['library']] 
)
