
library(dplyr)

# read in mutation sheet
# TODO: look up proper handling, probably somewhat similar to sample_sheet reading
args <- commandArgs(trailingOnly = TRUE)

sampleName <- args[1]
mutation_sheet <- args[2]

createSigMatrix <- function ( mutations.vector, sig_mutations.df ) {
  #' for making the signature matrix based on the signature mutations found in the sample (given as input as a vector)
  #' TODO: The data.frame should be builded dynamically based on a vector with variants, the filter step could be a function
  #' for it self
  #' returns simple signature matrix as data.frame without frequency values
  #' 
  # read in provided mutation sheet
  mutations.df <- read.csv(mutation_sheet)
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  
  msig <- setNames( data.frame( matrix( ncol = ncol(mutations.df)+1, nrow = 0 )), c("WT", colnames(mutations.df)))
  msig <- bind_rows(tibble(muts=mutations.vector), msig)
  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant to see how many characterising mutations
  # where found by variant
  
  for ( variant in names(mutations.df) ){
    msig[[ variant ]] <- msig$muts %in% mutations.df[[ variant ]]
  }
  msig[is.na(msig)] <- 0
  
  return( msig[,-match('muts', names(msig))]*1 ) # use the *1 to turn true/false to 0/1
}

simulateWT <- function ( mutations.vector, bulk_freq.vector, simple_sigmat.dataframe) {
  # for the deconvolution to work we need the "wild type" frequencies too. The matrix from above got mirrored, 
  # wild type mutations are simulated the following: e.g. T210I (mutation) -> T210T ("wild type")
  
  # 1. make "WT mutations" 
  muts_wt <- lapply(mutations.vector,function(x) str_replace(x,regex(".$"), str_sub(x, 1,1)))
  muts_wt.df <- data.frame(muts = unlist(muts_wt))
  # 2. make frequency values, subtract the observed freqs for the real mutations from 1
  bulk_wt <- lapply(bulk_freq.vector, function (x) {1-x})
  
  # 3. make matrix with wt mutations and inverse the values and wild type freqs
  msig_inverse <- bind_cols(muts_wt.df, as.data.frame(+(!simple_sigmat.dataframe[,-1])))
  
  # fixme: not sure if this really is a nice way to concat those things...
  muts_all <- c(muts_wt,muts)
  muts_all.df <- data.frame(muts = unlist(muts_all))
  
  bulk_all <- c(bulk_wt, bulk)
  bulk_all.df <- data.frame(freq = unlist(bulk_all))
  
  msig_all <- rbind(msig_inverse[,-1],msig[,-1])
  
  # 4. concat the data frames
  
  # without bulk freq for building the signature matrix
  msig_stable <- bind_cols(muts_all.df,msig_all)
  
  # with bulk freq for export and overview
  msig_stable_complete <- bind_cols(muts_all.df,msig_all,bulk_all.df)
  
  return ( msig_stable )
}