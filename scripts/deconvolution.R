# Title     : TODO
# Objective : TODO
# Created by: vfs
# Created on: 23.06.21

createSigMatrix <- function ( mutations.vector ) {
  #' for making the signature matrix based on the signature mutations found in the sample (given as input as a vector)
  #' TODO: The data.frame should be builded dynamically based on a vector with variants, the filter step could be a function
  #' for it self
  #' returns simple signature matrix as data.frame without frequency values
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  msig <- data.frame(muts,WT=0,b117=0,b1351=0,b1427=0,b1429=0,b1526=0,p1=0, b16172=0)
  
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

