get_genome_cov <- function ( coverage_dir ) {
  # TODO check assumptions e.g. if the coverage files exists in dir and if the dir exists
  require(dplyr)
  
  files <- list.files(path = coverage_dir,
                      pattern = "_coverage.csv",
                      full.names = TRUE, 
                      recursive = FALSE) # why false?
  if (length(files) > 0) { 
    genome_cov.df <- data.frame(samplename = c(),
                                ref_genome_coverage = c())
    genome_cov.df <- bind_rows(genome_cov.df, lapply( files, function (x){
      rd_tbl <- read.table(x, sep = '\t') # ignores the row starting with # directly TODO: how to get column names?
      # TODO check if file has the correct header format
      data <- c( samplename = strsplit( basename (x), "\\_coverage.csv" )[[1]],
                 ref_genome_coverage = rd_tbl$V6 )
      #df <- bind_rows(df, data)
      return(data)
    }))
  }else{ stop("coverage directory is empty or does not exist")}
  
  return(genome_cov.df) 
}

get_mutation_cov <- function ( coverage_dir ) {
  # TODO check assumptions e.g. if the coverage files exists in dir and if the dir exists
  require(dplyr)
  require(stringr)
  
  files <- list.files(path = coverage_dir,
                      pattern = "_merged_covs.csv",
                      full.names = TRUE, 
                      recursive = FALSE) # why false?
  mutation_cov.df <- data.frame(samplename = c(),
                                total_num_muts = c(),
                                total_muts_cvrd = c(),
                                drop_out_muts = c())
  mutation_cov.df <- bind_rows(mutation_cov.df, lapply( files, function (x){
    rd_tbl <- read.table(x, sep = '\t', header = TRUE)
    # TODO check if file has the correct header format
    data <- c( samplename = strsplit( basename (x), "\\_merged_covs.csv" )[[1]],
               total_num_muts = rd_tbl$Total.number.of.tracked.mutations,
               total_muts_cvrd = rd_tbl$Total.number.of.mutations.covered,
               drop_out_muts = str_remove_all(rd_tbl$Number.of.mutations.not.covered,"\\[|\\]|\\s"))
    return(data)
  }))
  
  return(mutation_cov.df)
}