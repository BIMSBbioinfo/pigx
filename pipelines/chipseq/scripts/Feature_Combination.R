# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Feature_Combination')

# ---------------------------------------------------------------------------- #
Feature_Combination = function(
  features,
  outfile,
  scriptdir
){

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    suppressPackageStartupMessages(library('stringr'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    # ------------------------------------------------------------------------ #
    message('Collapsing Peaks ...')
      feat.list = GRangesList(lapply(features, readGeneric))
      lnames = unique(str_replace(basename(unlist(features)),'.(narrow|broad)Peak',''))
      lnames = str_replace(lnames,'.bed$','')

      if(length(lnames) != length(feat.list))
          stop('feature_Combination: input features are duplicated')
      names(feat.list) = lnames
      print(lnames)

      if(length(feat.list) == 1){
        feat.comb = feat.list[[1]]
        feat.comb$class = names(feat.comb)
      }else{
          feat.comb = findFeatureComb(feat.list, use.names=TRUE)
      }
      feat.comb$peak_id = paste0('Peak',sprintf('%07d', seq_along(feat.comb)))

      # export RDS
      saveRDS(feat.comb, file = outfile)
      # export txt file with class
      write.table(x = as.data.frame(feat.comb),
                  file = gsub("\\.rds","\\.txt",outfile, ignore.case = TRUE),
                  sep = "\t", append = TRUE,
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      # export BED format file
      write.table(x = as.data.frame(feat.comb)[, c("seqnames", "start",
                                                   "end", "peak_id"
                                                   )],
                  file = gsub("\\.rds","\\.bed",outfile, ignore.case = TRUE),
                  sep = "\t", append = TRUE,
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ---------------------------------------------------------------------------- #
# function call
Feature_Combination(
    features         = argv$input,
    outfile          = argv$output[['outfile']],
    scriptdir        = argv$params[['scriptdir']]
)
