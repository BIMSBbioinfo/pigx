library(tidyverse)

# -----------TEST---------------
mutation_freq <- read_delim( "/mnt/beast/pathogenomics/vpipe/akalinlab_pathogenomics/misc_src/pigx_sarscov2_ww/tests/output/variants/data_mutation_plot.csv", "\t", 
                             escape_double = FALSE, trim_ws = TRUE )
# -----------TEST---------------

#open the table
# TODO file input
# mutation_freq <- read_delim( "C:/Users/rafae/Desktop/mutation_freq.csv", "\t", escape_double = FALSE, trim_ws = TRUE )

#  get names of mutations (everything after the meta info from the sample sheet)
mutations <- names ( 
            mutation_freq[ (which( names(mutation_freq) %in% "coordinates_long")+1) : length( names( mutation_freq ))])

#enforcing date type for column dates ...it will get rid of the time
mutation_freq$dates <- as.Date ( mutation_freq$dates )

#remove mutations with NA for all rows and create a new dataframe
# TODO what is x?
not_all_na <- function( x ) any( !is.na( x ) )
lm_df <- mutation_freq %>% select( c( dates, all_of(mutations) )) %>% select( where( not_all_na ))

#get a list of mutations that still present after removing those with only NA values
mutations_with_values <- names( lm_df[ 2 : length( names( lm_df ))])

#initialize 3 lists for results
results_lm <- list ()
summaries <- list ()
pvalues <- list ()

#loop the mutations, doing the linear model, and extracting the pvalues
for ( i in mutations_with_values ){
  tmp <- lm_df %>% select( dates, all_of(i) )
  test <- lm( formula = tmp[[i]] ~ tmp$dates )
  results_lm[[i]] <- test
  summaries[[i]] <- summary( test )
  pvalues[[i]] <- as.data.frame( summary( test )$coefficients[,4]) 
  colnames( pvalues[[i]] ) <- "pvalues"
}

#generate a proper dataframe with pvalues, filtering by significance
pvalues_df <- do.call( rbind, pvalues ) %>% 
              tibble::rownames_to_column( "VALUE" ) %>% 
              filter( str_detect( VALUE, "dates" ) ) %>% 
              filter( pvalues < 0.05 ) 

#names of mutations get strange pattern, doing some split and getting the names correct
pvalues_df$mutation <- str_split_fixed( pvalues_df$VALUE, "[.]",2 )[,1] # TODO this assumes values

#making a list with mutations that are significant
mutations_sig <- str_split_fixed( pvalues_df$VALUE, "[.]", 2 )[,1]

#fixing the dataframe with mutations and pvalues
pvalues_df <- pvalues_df %>% select( mutation, pvalues )

# TODO put this in a set of functions you call then...probably from the index report.