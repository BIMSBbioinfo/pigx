write_lm_results <- function ( df, output ){
  #' takes a dataframe with 3 col: samplename, p-values, sigmutflag
  #' writes them to a csv
  # this function is here because I don't know yet if we need to do smth else with it too
  write.csv(df, output)
}

write_mutations_count <- function ( mutation_plot_data, mutation_sheet, mutations_sig, output ){
  #' takes data_mut_plot.csv df, mutation_sheet.df, mutations_sig.df as input
  #' counts mutations and return them as a dataframe
  #'
  
  fields <- c("sample", "total_muts", 
              "total_sigmuts", "sigmuts_b117", "sigmuts_p1", "sigmuts_b16172", "sigmuts:b1351",
              "non-sigmuts", "tracked_muts_after_lm" )
  
  samples <- mutation_plot_data$samplenames
  
  #  get names of mutations (everything after the meta info from the sample sheet)
  mutations <- names ( mutation_plot_data[ (which( names(mutation_plot_data) %in% "coordinates_long")+1) : length( names( mutation_plot_data ))])
  total_muts <- length(mutations)
  
  # from mutation_plot_data select columns which are in mutation_sheet
  sigmuts.df <- mutation_plot_data %>% select(all_of(mutations) %in% mutation_sheet) # TODO maybe str_detect? or grepl?
  total_sigmuts <- length(sigmuts.df)
  
  non-sigmuts <- ( total_muts - total_sigmuts )
  
  # check if mutations_sig$mutations is in mutation_plot_data
    # if yes, count to "trakced_muts_after_lm" 
  
  # for each sample
    # check if value in "mutations col" contains value (not NA)
      # if yes, count to mutations_per_sample
        # check if it's also in "mutations_sig" 
          # if yes, count to "tracked_muts_after_lm"_per_sample
    # check if value in sigmuts.df contains value
      # if yes, count to sigmuts_per_sample
      # if yes and value in b117 col from mutation_sheet, count to sigmut_b117_per_sample (for all of them, probab. write a loop)
    #  "non-sigmuts" <- "mutations_per_sample" - "sigmuts_per_samples"
}