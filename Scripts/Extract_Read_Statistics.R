# -------------------------------------------------------------------------- #
Extract_Read_Statistics = function(
    bamfile = NULL,
    outfile = NULL,
    outname = 'Sample'

){
  if(is.null(bamfile))
    stop('bamfile not specified')

  if(is.null(outfile))
      stop('outfile not specified')

  suppressPackageStartupMessages(library(dropbead))
  suppressPackageStartupMessages(library(Rsamtools))

    message('Extracting read statistics ...')
      stat = extractReadStatistics(bamfile)
      rownames(stat) = outname

    message('Writing read statistics ...')
      write.table(stat, outfile,
        row.names=TRUE, col.names=TRUE,
        sep='\t',quote=FALSE)
}

# -------------------------------------------------------------------------- #
Extract_Read_Statistics(
      bamfile = snakemake@input[['bamfile']],
      outfile = snakemake@output[['outfile']],
      outname = snakemake@params[['outname']]
  )
