# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

args <- commandArgs(trailingOnly = TRUE)
prefix            <- args[1]
assembly          <- args[2]
targetdir         <- args[3]
diff.meth.reports <- tail(args, n=-3)


# if a sample doesnt take part in diff. meth.
# then dont change knitr_meta.rds for final_report rule
if(!is.null(diff.meth.reports) && (length(diff.meth.reports) > 0)){
  
  # there can be more than 1 diff meth report per sample
  # because the sample can be compared with more than
  # one for diff. meth. calling
  treatment_pairs = sapply(diff.meth.reports, function(x){
    basex=basename(x)
    strsplit(basex, paste0(".sorted_",assembly,"_annotation.diff.meth.nb.html"))
  })
  
  # 1. Copy session info from diff meth to the sample-specific directory
  diffmeth_dir = lapply(treatment_pairs, function(x) paste0(targetdir,x,"/"))
  sesion.info.diffmeth.pairs = lapply(diffmeth_dir, function(x) list.files(path =x, pattern = "session", full.names = TRUE))
  
  for(x in sesion.info.diffmeth.pairs){
    file.copy(x, paste0(targetdir,prefix,"/"))
  }
  
  # 2. Merge knitr_meta.rds together
  diffmeth_knitrmeta = lapply( treatment_pairs, function(x) readRDS( paste0(targetdir,x,"/knitr_meta.rds") ) )
  diffmeth_annot_knitrmeta = lapply( treatment_pairs, function(x) readRDS( paste0(targetdir,x,".",assembly,"/knitr_meta.rds") ) )
  
  diffmeth_knitrmeta = lapply(diffmeth_knitrmeta, function(x) x[[1]])
  names(diffmeth_knitrmeta) = paste0(tools::file_path_sans_ext(as.character(diffmeth_knitrmeta)), ".Rmd")
  
  diffmeth_annot_knitrmeta = lapply(diffmeth_annot_knitrmeta, function(x) x[[1]])
  names(diffmeth_annot_knitrmeta) = paste0(tools::file_path_sans_ext(as.character(diffmeth_annot_knitrmeta)), ".Rmd")
  
  
  final_knitrmeta = readRDS(paste0(targetdir,prefix,"/knitr_meta.rds"))
  
  
  # # If diff meth reports are already included in th final report, then dont include it again
  if( !(any(   names(diffmeth_knitrmeta) %in% names(final_knitrmeta)  ) )){
    
    merged_final_report=c(final_knitrmeta, diffmeth_knitrmeta, diffmeth_annot_knitrmeta)
    
    saveRDS(merged_final_report, paste0(targetdir,prefix,"/knitr_meta.rds"))
    
  }
}
