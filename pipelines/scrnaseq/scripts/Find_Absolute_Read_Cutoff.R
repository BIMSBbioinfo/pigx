# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Find_Absolute_Read_Cutoff')

# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff = function(
    infile  = NULL,
    outfile = NULL,
    cutoff  = 50000
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
    png(str_replace(outfile,'yaml','png'), width=400, height=300)
      p = plotCumulativeFractionOfReads(reads_by_cell,
                                    cutoff = cutoff,
                                    draw.knee.point = TRUE)
      print(p)
    dev.off()


  cutoff = min(cutoff, nrow(reads_by_cell))
  message('Print output yaml ...')
    infl_point = estimateCellNumber(reads_by_cell[,1], max.cells=cutoff)
    if(length(infl_point) == 0){
      message('knee cell selection did not succeed: including all cells')
      infl_point = min(reads_by_cell[,1])
    }
    lout       = list(reads_cutoff = reads_by_cell[infl_point,1])
    cat(as.yaml(lout), file=outfile)

}


# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff(
      infile  = argv$input[['infile']],
      outfile = argv$output[['outfile']],
      cutoff  = argv$params[['cutoff']]
  )
