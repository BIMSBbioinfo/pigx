# Title     : TODO
# Objective : TODO
# Created by: vfs
# Created on: 23.06.21

createSigMatrix <- function ( mutations.vector, sig_mutations.df ) {
  #' for making the signature matrix based on the signature mutations found in the sample (given as input as a vector)
  #' TODO: The data.frame should be builded dynamically based on a vector with variants, the filter step could be a function
  #' for it self
  #' returns simple signature matrix as data.frame without frequency values
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  msig <- data.frame(mutations.vector,WT=0,b117=0,b1351=0,b1427=0,b1429=0,b1526=0,p1=0, b16172=0)
  names(msig)[names(msig) == "mutations.vector"] <- "muts"
  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant to see how many characterising mutations
  # where found by variant 
  
  b117 <- filter(sig_mutations.df, name == "b117")$value
  msig[msig$muts %in% b117,"b117"]=1
  
  b1351 <- filter(sig_mutations.df, name == "b1351")$value
  msig[msig$muts %in% b1351,"b1351"]=1
  
  b1427 <- filter(sig_mutations.df, name == "b1427")$value
  msig[msig$muts %in% b1427,"b1427"]=1
  
  b1429 <- filter(sig_mutations.df, name == "b1429")$value
  msig[msig$muts %in% b1429,"b1429"]=1
  
  b1526 <- filter(sig_mutations.df, name == "b1526")$value
  msig[msig$muts %in% b1526,"b1526"]=1
  
  p1 <- filter(sig_mutations.df, name == "p1")$value
  msig[msig$muts %in% p1,"p1"]=1
  
  b16172 <- filter(sig_mutations.df, name == "b16172")$value
  msig[msig$muts %in% b16172, "b16172"]=1
  
  msig[1,2] <- 0 # put the WT signature, here it's 0 and not 1 bc it will be inversed in the next step
  
  return( msig )
}
# TODO write test about the expected structure of the matrix

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