#== BSseq_deconv.R -deconvolves BSseq pipeline output  (.cov file) into cell fractions ======
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
library(rtracklayer)

# whereami="MDC_server"
 whereami="laptop"

if(whereami=="MDC_server")
  {  
  WORKDIR          ="/home/bosberg/bs/pigx_bsseq"
  SIGMAT_PATH      ="/home/bosberg/bs/pigx_out/HK_Sun_data"
  ATLAS_WIG_FOLDER ="/data/akalin/bosberg/genome_atlas"
  }
elseif(whereami=="laptop")
  {
  ATLAS_WIG_FOLDER   ="/Users/blosberg/Desktop/Science/postdoc_MDC/bs/ref/genome_atlas"
  WORKDIR="/Users/blosberg/postdoc_work/bs/pigx_bsseq/"
  SIGMAT_PATH="/Users/blosberg/postdoc_work/bs/pigx_out/HK_Sun_data/"
  PATH_DATA  = "/Users/blosberg/postdoc_work/bs/pigx_out/HK_Sun_data/Cond_read_data/"; #--- laptop:
  }
elseif(whereami=="home_linuxbox")
  {
  PATH_DATA  ="/home/bren/Desktop/Science/postdoc_MDC/bs/pigx_out/HK_Sun_data/Cond_read_data/"; #--- home machine:
  }


setwd(WORKDIR)
source('./BSseq_deconv_funcs.R');
#================================================================================
#---------- IMPORT SIGNATURE MATRIX AND LIST OF BIOMARKER LOCATIONS -------------

Sundat = get_refdat(SIGMAT_PATH)

#================================================================================
#----------- READ IN YOUR TEST SAMPLE (bigWIGFILE IN THIS CASE):
# the genome atlas data are in .wig format, so
# First you need to convert your wig file to bigwig using: wigToBigWig()
# it should be included in one of the above packages
# to do this, you'll need to generate a Seq object, using the length and chromosome names provided in 
# hg19_std_chrom_sizes_ucsc_sorted.csv in the home folder on beast.
# seqdat=read.csv(file.path(ATLAS_WIG_FOLDER,"hg19_std_chrom_sizes_ucsc_sorted.csv"),sep="\t", header=FALSE ) 
# Seqobj <- Seqinfo(seqnames=as.character(seqdat[,1]), seqlengths=seqdat[,2])
# wigToBigWig(file.path(ATLAS_WIG_FOLDER, wigcov_file) , seqinfo=Seqobj )

# wigcov_file   = "UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.readcoverage.wig" # --> already converted to bigwig (binary)
# wig_file      = "UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.wig" # --> already converted to bigwig (binary)
bw_rc_dat_file  = "UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.readcoverage.bw"
bw_dat_file     = "UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.bw"

#--- now you can feed your bigwig files into the function

Atlas_dat = get_Atlas_dat(ATLAS_WIG_FOLDER, bw_dat_file, bw_rc_dat_file , Sundat$Biomarks_gr)



#================================================================================
#----------- READ IN YOUR SAMPLE DATA:
Cond_dat          = list()
deconv_out        = list()
deconv_subset_out = list()
mincounts  = 10;  
prefix     = "Cond"
suffix     = ".read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov"

#================================================================

NCT=14;
epsilon=0.01
conditions <- get_conditions(NCT, 1, 1, epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
ui_COND = conditions$ui_COND; 
ci_COND = conditions$ci_COND; 

CTsubset=c(1, 11, 12, 13, 14)
NCT_subset =length(CTsubset)
conditions_subset <- get_conditions(NCT_subset, 1, 1, epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
ui_COND_subset = conditions_subset$ui_COND; 
ci_COND_subset = conditions_subset$ci_COND; 


residue=matrix(0,14,1)
unit_CT_meth_profiles = as.matrix(Sundat$Sigmat) %*% diag(14) 
for (k in c(1:14))
  {
  residue[k,1] = norm( unit_CT_meth_profiles[,k] - Atlas_dat$ROI_meth_profile, "F")
  }

i=4
# for (i in c(1:18))
# {
#  print(paste(" ---------- Working on dataset i = ",as.character(i), " ----------------" ) )

   #============     Start deconvolving:  ==================
fin=paste(PATH_DATA, prefix, as.character(i), suffix, sep="")
Cond_dat[[i]]         = get_Cond_dat(fin, refdat=Sundat , mincounts=1 )
Sigmat_whits_CTsubset = Cond_dat[[i]]$Sigmat_whits[,CTsubset]

temp=dim(Cond_dat[[i]]$Sigmat_whits);
L=temp[1]

   w              = conditions$w
   w_subset       = conditions_subset$w
 
  # deconv_out        = constrOptim(w, function(w){func_SQ(w, y_target = Atlas_dat$ROI_meth_profile, SigMat_in= as.matrix(Sundat$Sigmat))  },
  #                                     ui=ui_COND, ci=ci_COND, mu = 1e-04, control = list(), method = "Nelder-Mead") 

   deconv_out[[i]]        = constrOptim(w, function(w){func_SQ(w, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in= Cond_dat[[i]]$Sigmat_whits)}, ui=ui_COND, ci=ci_COND, mu = 1e-04, control = list(), method = "Nelder-Mead") 
   deconv_subset_out[[i]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_whits_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = "Nelder-Mead") 

   #} --- end of for loop
   
   plot(c(1:NCT),deconv_out[[1]]$par, c(1:NCT),deconv_out[[2]]$par, c(1:NCT), deconv_out[[3]]$par, c(1:NCT), deconv_out[[4]]$par, color=c("black","red","blue","green") )
   
   plot(c(1:NCT),deconv_out[[1]]$par)
   
   plot(c(1:NCT_subset), deconv_subset_out[[1]]$par)
   plot(c(1:NCT_subset), deconv_subset_out[[2]]$par)
   plot(c(1:NCT_subset), deconv_subset_out[[3]]$par)
   plot(c(1:NCT_subset), deconv_subset_out[[4]]$par)
   
   methods=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",   "SANN", "Brent")
   #--- for other methods you have to supply the gradient.
     for( j in c(1:length(methods)) ) 
     {
     m=methods[j]
     deconv_subset_out[[j]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, y_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_whits_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = m) 
     deconv_subset_out[[j]]$par
     }

   
   # }

N=matrix(0,NCT,1)
for( k in c(1:NCT))
{
  N[k,1] =  norm(as.matrix(Cond_dat[[i]]$ROI_meth_profile -  Cond_dat[[i]]$Sigmat_whits[,k] ) ,"F")
}

