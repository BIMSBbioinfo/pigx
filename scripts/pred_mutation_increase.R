
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


lin_reg_mutation_change <- function( mutations.df ){
  #' takes data frames with mutations, frequency values over time
  #' returns A) list with mutations with significant change, B) dataframe with related pvalues
  
  require(tidyverse)
  require(stringr)
  # TODO check file format assumptions
  
  #initialize 3 lists for results
  results_lm <- list ()
  summaries <- list ()
  pvalues <- list ()
  
  #loop the mutations, doing the linear model, and extracting the pvalues
  for ( i in names(mutations.df[,-which(names(mutations.df) %in% "dates")]) ){
    if ( length(na.omit(mutations.df[,i])) >= 3 ){
      tmp <- mutations.df %>% select( dates, all_of(i) )
      test <- lm( formula = tmp[[i]] ~ tmp$dates )
      # only write the p-values for positive coefficients
      if (test["coefficients"]$coefficients["tmp$dates"] > 0){
          results_lm[[i]] <- test
          summaries[[i]] <- summary( test )
          pvalues[[i]] <- as.data.frame( summary( test )$coefficients[,4])
      }
    }
  }
  
  #generate a proper dataframe with pvalues, filtering by significance
  pvalues_df <- do.call( rbind, pvalues ) %>% 
                tibble::rownames_to_column( "VALUE" ) %>% 
                filter( stringr::str_detect( VALUE, "dates" ) ) %>%
                filter (`summary(test)$coefficients[, 4]` < 0.01) 
  
  #names of mutations get strange pattern, doing some split and getting the names correct
  pvalues_df$mutation <- str_split_fixed( pvalues_df$VALUE, "[.]",2 )[,1]
  
  #fixing the dataframe with mutations and pvalues
  pvalues_df <- pvalues_df %>% select( mutation, `summary(test)$coefficients[, 4]` )
  colnames(pvalues_df) <- c("mutation", "pvalues")
  return ( pvalues_df)
}