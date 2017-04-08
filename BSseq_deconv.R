#== BSseq_deconv_funcs.R -deconvolves BSseq pipeline output  (.cov file) into cell fractions ======
# ---last updated on  Sat Apr 8 23:11:46 CEST 2017  by  blosberg  at location  , brenosbsmacbook.fritz.box
#
#  changes from  Sat Apr 8 23:11:46 CEST 2017 : Made first draft template: script executes bug free but output inconsistent with expectations. Further testing needed
#
# By Bren Osberg, MDC Berlin, started October 2016.
#==================================================================================================

# install.packages("stringr", repos='http://cran.us.r-project.org')
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
#biocLite("methylKit")
## install.packages("devtools")
# library("devtools")
# install_github("al2na/methylKit", build_vignettes=FALSE, repos=BiocInstaller::biocinstallRepos(),ref="1.1.7BioC3.4", dependencies=TRUE)
rm(list=ls()) # CLEAN UP EVERYTHING

# install.packages("methylKit")
# methylKit::methRead() 

library(GenomicRanges)
library("stringr")
library(methylKit)
source('~/postdoc_work/bs/pigx_bsseq/BSseq_deconv_funcs.R');

#================================================================================
#---------- IMPORT SIGNATURE MATRIX AND LIST OF BIOMARKER LOCATIONS -------------

SIGMAT_PATH="/Users/blosberg/postdoc_work/bs/pigx_out/HK_Sun_data/"
Sundat = get_refdat(SIGMAT_PATH)

#================================================================================
#----------- READ IN YOUR SAMPLE DATA:
Cond_dat          = list()
deconv_out        = list()
deconv_subset_out = list()

mincounts  = 10;  

#PATH_DATA  ="/home/bren/Desktop/Science/postdoc_MDC/bs/pigx_out/HK_Sun_data/Cond_read_data/"; #--- home machine:
PATH_DATA  = "/Users/blosberg/postdoc_work/bs/pigx_out/HK_Sun_data/Cond_read_data/"; #--- laptop:
prefix     = "Cond"
suffix     = ".read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov"

NCT=14;
epsilon=0.01
conditions <- get_conditions(NCT, 1, 0, epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
ui_COND = conditions$ui_COND; 
ci_COND = conditions$ci_COND; 
w       = conditions$w

i=1
# for (i in c(1:18))
# {
#  print(paste(" ---------- Working on dataset i = ",as.character(i), " ----------------" ) )
   fin=paste(PATH_DATA, prefix, as.character(i), suffix, sep="")

   Cond_dat[[i]] = get_Cond_dat(fin, refdat=Sundat , mincounts=1 )

   temp=dim(Cond_dat[[i]]$Sigmat_red);
   L=temp[1]

   #============     Start deconvolving:  ==================
   
   CTsubset=c(1,13,14)
   Sigmat_red_CTsubset = Cond_dat[[i]]$Sigmat_red[,CTsubset]
   
   NCT_subset =length(CTsubset)
   conditions_subset <- get_conditions(NCT_subset, 1, 0, epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
   ui_COND_subset = conditions_subset$ui_COND; 
   ci_COND_subset = conditions_subset$ci_COND; 
   w_subset       = conditions_subset$w
   
   deconv_out[[i]]        = constrOptim(w, function(w){func_SQ(w, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in= Cond_dat[[i]]$Sigmat_red)}, ui=ui_COND, ci=ci_COND, mu = 1e-04, control = list(), method = "Nelder-Mead") 
   deconv_subset_out[[i]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_red_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = "Nelder-Mead") 

   methods=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",   "SANN", "Brent")
   #--- for other methods you have to supply the gradient.
     for( j in c(1:length(methods)) ) 
     {
     m=methods[j]
     deconv_subset_out[[j]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_red_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = m) 
     deconv_subset_out[[j]]$par
     }

   
   # }

N=matrix(0,NCT,1)
for( k in c(1:NCT))
{
  N[k,1] =  norm(as.matrix(Cond_dat[[i]]$ROI_meth_profile -  Cond_dat[[i]]$Sigmat_red[,k] ) ,"F")
}

