# ---------------------------------------------------------------------------- #
Feature_Combination = function(
  features,
  scriptdir,
  outfile
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
      feat.comb$peak_id = paste0('Peak',sprintf('%07d', 1:length(feat.comb)))
      

      saveRDS(feat.comb, file = outfile)

}




# ---------------------------------------------------------------------------- #
# function call
Feature_Combination(
    features         = snakemake@input,
    outfile          = snakemake@output[['outfile']],
    scriptdir        = snakemake@params[['scriptdir']]
)

# feature_Combination(
#     features    =args[['features']],
#     bw          =args[['bw']],
#     annotation  = args[['annotation']],
#     outpath     = args[['outpath']],
#     scriptdir   = args[['scriptdir']]
# )
