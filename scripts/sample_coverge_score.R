get_genome_cov <- function ( coverage_dir ) {
  # TODO check assumptions e.g. if the coverage files exists in dir and if the dir exists
  require(dplyr)
  
  files <- list.files(path = coverage_dir,
                      pattern = "_coverage.csv",
                      full.names = TRUE, 
                      recursive = FALSE) # why false?
  genome_cov.df <- data.frame(samplename = c(),
                              coverage = c())
  genome_cov.df <- bind_rows(genome_cov.df, lapply( files, function (x){
    rd_tbl <- read.table(x, sep = '\t') # ignores the row starting with # directly TODO: how to get column names?
    # TODO check if file has the correct header format
    data <- c( samplename = strsplit( basename (x), "\\_coverage.csv" )[[1]],
               coverage = rd_tbl$V6 )
    #df <- bind_rows(df, data)
    return(data)
  }))
  
  return(genome_cov.df) 
}