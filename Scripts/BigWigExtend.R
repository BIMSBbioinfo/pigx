# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Scripts/Argument_Parser.R'))
argv = Parse_Arguments('Extend_Regions')

# ---------------------------------------------------------------------------- #
Extend_Regions = function(
    inpath,
    outpath,
    extend       = NULL,
    scale_index  = FALSE,
    scale_factor = 1e6
){

    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(rtracklayer))
    if(!is.numeric(extend))
        stop('Extend_Regions: extend parameter needs to be a number')

    if(!file.exists(inpath))
        stop('Extend_Regions: input file does not exist')


    f = fread(inpath)
    setnames(f,c('chr','start','end','name','x','strand'))
    f[strand == '+', end := start[strand == '+'] + as.integer(extend)]
    f[strand == '-', start := end[strand == '-'] - as.integer(extend)]
    if(any(f$start < 0)){
        warning('Removing reads with start < 0')
        f = f[f$start >= 0]
    }

    g = makeGRangesFromDataFrame(as.data.frame(f))
    cov = coverage(g)

    if(scale_index == 'yes' || scale_index)
        cov = round(cov * (scale_factor/nrow(f)),2)

    export.bw(cov, outpath)
}

# ---------------------------------------------------------------------------- #
Extend_Regions(
  inpath        = argv$input[['file']],
  outpath       = argv$output[['outfile']],
  extend        = argv$params[['extend']],
  scale_index   = argv$params[['scale']])
