
parsing_mutation_plot_data <- function ( mutation_plot_data ){
  #' taking csv as input which will be made by the variant reports
  #' returns df containing only mutations with frequency values
  
  require(tidyverse)
  # TODO check file format assumptions
  
  #  get names of mutations (everything after the meta info from the sample sheet)
  mutations <- names ( 
            mutation_plot_data[ (which( names(mutation_plot_data) %in% "coordinates_long")+1) : length( names( mutation_plot_data ))])
  
  # enforcing date type for column dates ...it will get rid of the time
  mutation_plot_data$dates <- as.Date ( mutation_plot_data$dates )
  
  # remove mutations with NA for all rows and create a new dataframe
  # fixme Vic I think there is a more compact way for this
  not_all_na <- function( x ) any( !is.na( x ) )
  lm_df <- mutation_plot_data %>% select( c( dates, all_of(mutations) )) %>% select( where( not_all_na ))
  
  return ( lm_df )
}

refined_lm_model <- function( mutations.df ){
 #' takes data frames with mutations, frequency values over time
 #' returns dataframe with related pvalues
  
  require(tidyverse)
  require(stringr)
 
  # TODO check file format assumptions
  
  # for mutation frequency - a missing value can be assumed as "not found", so NA can be set to 0
  mutations.df <- mutations.df %>% replace(is.na(.), 0)
  
  #initialize 3 lists for results
  results_lm <- list ()
  summaries <- list ()
  pvalues <- list () 
  coeff <- list()
  
  #loop the mutations, doing the linear model, and extracting the pvalues
  for ( i in names( mutations.df[,-which(names(mutations.df) %in% "dates")]) ){
    # take only mutations that were found in at least 5 samples (equals 5 values above 0) to ensure robust model calculation
    if ( nrow( filter( tmp, values > 0 ) ) >= 5 ){
      # only select mutation columns
      tmp <- mutations.df %>% select( dates, all_of(i) )
      # perform linear regression
      test <- lm( formula = tmp[[i]] ~ tmp$dates )
      # extract results from model
      coeff[[i]] <- as.data.frame(test["coefficients"]$coefficients["tmp$dates"])
      results_lm[[i]] <- test
      summaries[[i]] <- summary( test )
      pvalues[[i]] <- as.data.frame( summary( test )$coefficients[,4])
      }
    }
  
  # generate results data frame with all values
  coeff_df <- as.data.frame(do.call(rbind, coeff)) %>% tidyr::drop_na() %>%
              # generate a proper dataframe from coeff values
              tibble::rownames_to_column( "VALUE" ) %>% 
              # proper column names
              rename(mutation = VALUE, coefficients = "test[\"coefficients\"]$coefficients[\"tmp$dates\"]" )
  
  lm_res.df <-  do.call( rbind, pvalues ) %>%
                # make proper df from p-values
                tibble::rownames_to_column( "VALUE" ) %>% 
                # only take the actual p-value
                filter( stringr::str_detect( VALUE, "dates" ) ) %>%
                # split the suffix from the model away
                mutate(VALUE = str_split(VALUE, ".tmp", simplify = TRUE)[,1]) %>%
                # proper column names
                rename(mutation = VALUE, pvalues = "summary(test)$coefficients[, 4]" ) %>%
                # join coeffs and pvalues togehter
                left_join(coeff_df, by = "mutation")
  
  return ( lm_res.df )
}

filter_lm_res_top20 <- function ( lm_res.df, pvalue_cutoff){
  #' lm_res.df: a dataframe with variables, pvalues and regression coefficients 
  #' pvalue_cutoff: a numeric value to use as p-value cutoff for filtering
  
  # generate dataframe with significant results only
  lm_res_sig.df <- lm_res.df %>% 
                    # filter for increasing trends only
                    filter( coefficients > 0) %>%
                    # filter for significance
                    filter( pvalues < pvalue_cutoff) %>% 
                    # sort for decreasing coeffs
                    arrange( desc(coefficients)) %>% 
                    # only take the 20 strongest trends
                    slice_head(n = 20)
  
  return ( lm_res_sig.df )
}