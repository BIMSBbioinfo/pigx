group_by_date <- function( df ){
  require(dplyr)
  df %>% 
        # discard location
        dplyr::select ( -location_name, -coordinates_long, -coordinates_lat) %>%
        # discard time only keep day
        mutate(dates = as.Date( dates )) %>%
        # pool samples per date and calc. the mean for every mutation column
        group_by(dates)
}

pool_by_mean <- function(df, na_handling) {
  #' docstring missing
  #' 
    df_pooled <- group_by_date(df) %>%
                summarize(across(where(is.numeric), mean, na.rm = na_handling)) %>%
                # rename samples to indicated that they were pooled
                mutate(samplename = paste0(dates, "_pooled")) %>%
                # put names first again
                relocate(samplename) %>%
                ungroup()
    return (df_pooled)
}

read_files <- function ( sample_sheet.df ){
  samples <- sample_sheet.df$name
  do.call( bind_rows, lapply(samples, apply_fun_lookup, sample_sheet = sample_sheet.df))
}

apply_fun_lookup <- function ( sample, sample_sheet.df ){
  sample_row <- which(sample_sheet.df$name == sample)
  # works only for paired end for now - I don't know how to make num of returned read files conditional based on sample_sheet columns
  # TODO make num of returend reads depended on paired or single end
  data.frame( samplename = sample,
              raw_reads1 = sample_sheet.df[sample_row,"reads"],
              raw_reads2 = sample_sheet.df[sample_row, "reads2"])
}

read_num_raw <- function ( raw_reads_vector, reads_dir){
  do.call( bind_rows, lapply(raw_reads_vector, apply_fun_get_read_num, reads_dir = reads_dir))
}

apply_fun_get_read_num <- function (read, reads_dir) {
  read_num <- as.numeric(
                      system2(command = "echo", 
                      args = c ("$(zcat ", file.path(reads_dir, read), "|wc -l)/4|bc"),
                      stdout = TRUE))
  data.frame( read_num )
}

get_num_raw_reads <- function (reads_dir, sample_sheet){
  
  sample_sheet.df <- read.csv(sample_sheet, header = TRUE)
  # get read files matching samples
  cat("get samples and reads from sample_sheet...\n")
  read_counts <- read_files ( sample_sheet.df )
  # fixme can do in one lline with dplyr and across I think
  read_counts$reads_r1 <- read_num_raw(read_counts$raw_reads1, reads_dir)$read_num
  read_counts$reads_r2 <- read_num_raw(read_counts$raw_reads2, reads_dir)$read_num
  read_counts <- read_counts %>% mutate( total_reads = reads_r1 + reads_r2 )
  
  return(read_counts)
}

pool_by_weighted_mean <- function(df, weights) {
  #' docstring missing
  #' weigths is a dataframe with minimum samplenames and total_reads as column 
  #' total reads is the number of reads used for alignment of one sample, should be the sum of read1 and read2 with 
  #' paired end data
  require(dplyr)

  weights <- weights  %>% 
              # only take weights from approved samples
              semi_join(df, by = "samplename")
  # get variants from data frame, being every column after the metadata
  variants <- names (
              df[ ( which( names(df) %in% "coordinates_long")+1) : length( names( df ))]
  )
  
  df_pooled <- group_by_date(df) %>%
                left_join(weights, by = "samplename")  %>%
                relocate(total_reads) %>%
                # summarize by calc weighted mean, with the num of raw reads as weight
                summarise_at(vars( all_of(variants) ), list(~ weighted.mean(., total_reads))) %>%
                # rename samples to indicated that they were pooled
                mutate(samplename = paste0(dates, "_pooled")) %>%
                # put names first again
                relocate(samplename) %>%
                ungroup()
}

pool_by_mean <- function(df, na_handling) {
  #' docstring missing
  #' 
    df_pooled <- group_by_date(df) %>%
                summarize(across(where(is.numeric), mean, na.rm = na_handling)) %>%
                # rename samples to indicated that they were pooled
                mutate(samplename = paste0(dates, "_pooled")) %>%
                # put names first again
                relocate(samplename) %>%
                ungroup()
    return (df_pooled)
}