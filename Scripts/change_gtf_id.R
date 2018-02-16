# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')


# -------------------------------------------------------------------------- #
change_gtf_id = function(
  infile  = NULL, 
  outfile = NULL
){

  library(GenomicRanges)
  library(rtracklayer)

  g = import.gff2(infile)
  g$gene_name = g$gene_id
  g$transcript_name = g$transcript_id
  export.gff2(g, outfile)
}

# -------------------------------------------------------------------------- #
change_gtf_id(
  infile  = argv$input[['infile']],
  outfile = argv$output[['outfile']]
  )
