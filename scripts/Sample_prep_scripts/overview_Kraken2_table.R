library(dplyr)

concat_kraken_output <- function ( kraken_dir ){
    #' docstring missing
    #' 
    require(dplyr)
    
    kraken_files <- list.files ( path = kraken_dir,
                                 pattern = "_classified_unaligned_reads.txt",
                                 full.names = TRUE,
                                 recursive = FALSE)
    # concat kraken results of all samples together
    df <- do.call( bind_rows, lapply(kraken_files, read_kraken ))
    # sort for decreasing number of assigned fragments
    df <- df[ order( df$num_assigned_fragments, decreasing = TRUE),]
    # remove indented tax_class spacing
    df <- as_tibble( lapply( df, function (x) {gsub("[[:space:]]", "", x)} ))
    # sum number of fragments per tax_class across all samples 
    df_summed <- df %>% dplyr::select(-sample) %>% 
                 group_by(taxon_class, rank ) %>% 
                 summarize( num_cov_fragments = sum(as.numeric(num_cov_fragments)),
                            num_assigned_fragments = sum(as.numeric(num_assigned_fragments)))
    return( list( df, df_summed ))
}

read_kraken <- function ( kraken_file ){
    #' docstring missing
    #' 
    
    rd_tbl <- read.table( kraken_file, sep = "\t" )
    colnames (rd_tbl) <- c( "proportion", "num_cov_fragments", "num_assigned_fragments", "rank", "ncbi_taxID", "tax_class" )

    data.frame( sample = strsplit( basename( kraken_file ), "_classified_unaligned_reads.txt")[[1]], 
                num_cov_fragments = rd_tbl$num_cov_fragments,
                num_assigned_fragments = rd_tbl$num_assigned_fragments,
                rank = rd_tbl$rank,
                taxon_class = rd_tbl$tax_class
    )
}

args <- commandArgs (trailingOnly=TRUE)
kraken_dir <- args[1]
output <- args[2]

Kraken_summary <- concat_kraken_output( kraken_dir )
summary_all <- Kraken_summary[[1]]
summary_summed <- Kraken_summary[[2]]

write.csv( summary_all, output, row.names = FALSE )
write.csv( summary_summed, paste0( 'summed_', output ), row.names = FALSE )