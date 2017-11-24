# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff = function(
    infile  = NULL,
    outfile = NULL
){
  if(is.null(infile))
    stop('infile not specified')

  if(is.null(outfile))
      stop('outfile not specified')

  suppressPackageStartupMessages(library(dropbead))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(yaml))

  reads_by_cell = read.table(infile)

  message('Plot inflection point ...')
    png(str_replace(outfile,'yaml','pdf'),width=400, height=300)
      plotCumulativeFractionOfReads(reads_by_cell,
                                    cutoff = 10000,
                                    draw.infl.point = TRUE)
    dev.off()


  message('Print output yaml ...')
    infl_point = findInflectionPoint(reads_by_cell[,1], max.cells=10000)
    lout       = list(reads_cutoff = reads_by_cell[infl_point,1])
    cat(as.yaml(lout), file=outfile)

}

# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff(
      infile  = snakemake@input[['infile']],
      outfile = snakemake@output[['outfile']]
  )
