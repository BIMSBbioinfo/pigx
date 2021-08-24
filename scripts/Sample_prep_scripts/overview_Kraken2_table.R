library(dplyr)

concat_kraken_output <- function ( kraken_dir ){
  kraken_files <- list.files ( path = kraken_dir, 
                        pattern = "_classified_unaligned_reads.txt", 
                        full.names = TRUE, 
                        recursive = FALSE)
  df <- do.call(bind_rows, lapply(kraken_files, read_kraken))
  df <- as_tibble(lapply(df, function (x) {gsub("[[:space:]]", "", x)}))
}

read_kraken <- function ( kraken_file ){
  rd_tbl <- read.table( kraken_file, sep = "\t" )
  colnames (rd_tbl) <- c( "proportion", "num_cov_fragments", "num_assigned_fragments", "rank", "ncbi_taxID", "tax_class" )
  
  data.frame( proportion = rd_tbl$proportion,
              taxon_class = rd_tbl$tax_class,
              sample = strsplit( basename( kraken_file ), "_classified_unaligned_reads.txt")[[1]]
              )
}

args <- commandArgs (trailingOnly=TRUE)
kraken_dir <- args[1]
output <- args[2]

df <- concat_kraken_output( kraken_dir )

write.csv(df, output, row.names = FALSE)