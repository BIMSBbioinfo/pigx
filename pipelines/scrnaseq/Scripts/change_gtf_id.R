change_gtf_id = function(infile, outfile){

  library(GenomicRanges)
  library(rtracklayer)

  g = import.gff2(infile)
  g$gene_name = g$gene_id
  g$transcript_name = g$transcript_id
  export.gff2(g, outfile)
}


change_gtf_id(
  infile  = snakemake@input[['infile']],
  outfile = snakemake@output[['outfile']]
  )
