
prefix = snakemake@wildcards$prefix
assembly = snakemake@wildcards$assembly

diff.meth.reports = snakemake@input[['diffmeth']]
# there can be more than 1 diff meth report per sample
# because the sample can be compared with more than
# one for diff. meth. calling
treatment_pairs = sapply(diff.meth.reports, function(x){
  basex=basename(x)
  strsplit(basex, paste0(".sorted_",assembly,"_annotation.diff.meth.nb.html"))
  })

# 1. Copy session info from diff meth to the sample-specific directory
diffmeth_dir = lapply(treatment_pairs, function(x) paste0("final_Report/",x,"/"))
sesion.info.diffmeth.pairs = lapply(diffmeth_dir, function(x) list.files(path =x, pattern = "session", full.names = TRUE))

for(x in sesion.info.diffmeth.pairs){
  file.copy(x, paste0("final_Report/",prefix,"/"))
}

# 2. Merge knitr_meta.rds together
diffmeth_knitrmeta = lapply( treatment_pairs, function(x) readRDS( paste0("final_Report/",x,"/knitr_meta.rds") ) )
diffmeth_annot_knitrmeta = lapply( treatment_pairs, function(x) readRDS( paste0("final_Report/",x,".",assembly,"/knitr_meta.rds") ) )

diffmeth_knitrmeta = lapply(diffmeth_knitrmeta, function(x) x[[1]])
names(diffmeth_knitrmeta) = paste0(tools::file_path_sans_ext(as.character(diffmeth_knitrmeta)), ".Rmd")

diffmeth_annot_knitrmeta = lapply(diffmeth_annot_knitrmeta, function(x) x[[1]])
names(diffmeth_annot_knitrmeta) = paste0(tools::file_path_sans_ext(as.character(diffmeth_annot_knitrmeta)), ".Rmd")


final_knitrmeta = readRDS(paste0("final_Report/",prefix,"/knitr_meta.rds"))


# # If diff meth reports are already included in th final report, then dont include it again
if( !(any(   names(diffmeth_knitrmeta) %in% names(final_knitrmeta)  ) )){

  merged_final_report=c(final_knitrmeta, diffmeth_knitrmeta, diffmeth_annot_knitrmeta)

  saveRDS(merged_final_report, paste0("final_Report/",prefix,"/knitr_meta.rds"))

}



