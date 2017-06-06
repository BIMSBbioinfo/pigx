Extend_Regions = function(inpath, outpath, extend=NULL){

  if(!is.numeric(extend))
    stop('Extend_Regions: extend parameter needs to be a number')

  if(!file.exists(inpath))
    stop('Extend_Regions: input file does not exist')

  library(data.table)
  f = fread(inpath)
  setnames(f,c('chr','start','end','name','x','strand'))
  f[strand == '+', end := start[strand == '+'] + as.integer(extend)]
  f[strand == '-', start := end[strand == '-'] - as.integer(extend)]
  if(any(f$start > 0)){
      warning('Removing reads with start < 0')
      f = f[f$start > 0]
  }
  write.table(f, outpath, row.names=F, col.names=F, quote=F, sep="\t")
}


Extend_Regions(
  inpath  = snakemake@input[['file']],
  outpath = snakemake@output[['file']],
  extend  = snakemake@params[['extend']])
