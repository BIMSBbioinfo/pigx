#== BSseq_deconv.R -deconvolves BSseq pipeline output  (.cov file) into cell fractions ======
# ---last updated on Wednesday April 19th 19:22
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

# whereami="home_linuxbox"
 whereami="MDC_server"
# whereami="laptop"

if(whereami=="MDC_server")
  {  
  WORKDIR          ="/home/bosberg/bs/pigx_bsseq"
  SIGMAT_PATH      ="/home/bosberg/bs/pigx_out/HK_Sun_data"
  ATLAS_WIG_FOLDER ="/data/akalin/bosberg/genome_atlas"
  PATH_DATA        ="/data/akalin/bosberg/Sun_PNAS_processed_output/MIX_FINAL_covfiles"
  }
elseif(whereami=="laptop")
  {
  #ATLAS_WIG_FOLDER   ="/Users/blosberg/Desktop/Science/postdoc_MDC/bs/ref/genome_atlas"
  ATLAS_WIG_FOLDER   ="/Volumes/Seagate Expansion Drive/Science_bak/postdoc_MDC/bs/ref/genome_atlas"
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
source('./BS_Selfconsist_check_deconv.R');

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

bw_file_list    = read.table(  file.path(ATLAS_WIG_FOLDER, "files_bw")   )
bwcov_file_list = read.table(  file.path(ATLAS_WIG_FOLDER, "files_bw_coverage")   )
bw_TT_names     = read.table(  file.path(ATLAS_WIG_FOLDER, "files_Tissue_types" )  )

temp = dim(bw_file_list)
Num_atlas_files=temp[1]

#--- now you can feed your bigwig files into the function
Atlas=list()

for (i in c(1:Num_atlas_files))
  {
  print(paste("-- Now importing Atlas data file ", as.character(i), " of ", as.character(Num_atlas_files) , " ---"))
  bw_file    = as.character(bw_file_list[i,1])
  bwcov_file = as.character(bwcov_file_list[i,1])
  
  Atlas[[i]] = get_Atlas_dat(ATLAS_WIG_FOLDER, bw_file, bwcov_file , Sundat$Biomarks_gr, as.character(bw_TT_names[i,1]) ) 
  }
# plot(Atlas[[1]]$ROI_meth_profile, Sundat$Sigmat[,8] )
# par(mfrow=c(2,2))
# plot(Atlas[[1]]$ROI_meth_profile, Sundat$Sigmat[,6],  xlab='Adipose -homemade', ylab = "Adrenal (Sun)", main="")
# plot(Atlas[[2]]$ROI_meth_profile, Sundat$Sigmat[,8],  xlab='Adrenal -homemade', ylab = "Adipose (Sun)", main="")
# plot(Atlas[[1]]$ROI_meth_profile, Sundat$Sigmat[,8],  xlab='Adipose -homemade', ylab = "Adipose (Sun)", main="" )
# plot(Atlas[[2]]$ROI_meth_profile, Sundat$Sigmat[,6],  xlab='Adrenal -homemade', ylab = "Adrenal (Sun)", main="" )

#================================================================================
#----------- READ IN YOUR SAMPLE DATA:
 Cond_dat   = list()
 mincounts  = 4;  
 prefix     = "Cond"
 suffix     = ".read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
 Num_samples= 18

 for (i in c(1:Num_samples))
   {
   print(paste("-- Now importing data from Condition mix ", as.character(i), " of ", as.character(Num_samples) , " ---")) 
   fin=file.path(PATH_DATA, paste(prefix, as.character(i), suffix, sep="") )
   print(paste(" using file: ", fin, "\n"))
   
   Cond_dat[[i]]         = get_Cond_dat(fin, refdat=Sundat , mincounts )   
   }
 
 

# Sigmat_whits_CTsubset = Cond_dat[[i]]$Sigmat_whits[,CTsubset]
# 
# temp=dim(Cond_dat[[i]]$Sigmat_whits);
# L=temp[1]
# deconv_out[[i]]        = constrOptim(w, function(w){func_SQ(w, yEXP_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in= Cond_dat[[i]]$Sigmat_whits)}, ui=ui_COND, ci=ci_COND, mu = 1e-04, control = list(), method = "Nelder-Mead") 
# deconv_subset_out[[i]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, yEXP_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_whits_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = "Nelder-Mead") 
# 
#================================================================

NCT=14;
 
epsilon=0.01


CT_subset=c(1, 11, 12, 13, 14)
NCT_subset =length(CTsubset)

residue=matrix(0,14,1)
unit_CT_meth_profiles = as.matrix(Sundat$Sigmat) %*% diag(14) 
unit_means=rowMeans(unit_CT_meth_profiles)

residue=matrix(0,Num_atlas_files,Num_samples)
for (i in c(1:Num_atlas_files))
  {
  for(j in c(1:Num_samples))
  {
  residue[i,j] = norm( as.matrix( Cond_dat[[j]]$ROI_meth_profile ) - as.matrix( Atlas[[i]]$ROI_meth_profile[ Cond_dat[[j]]$biomark_hits, 1 ])  , "F")
  }
}

#=========================================================================================


deconv_out          = list()
deconv_subset_out   = list()

p_init        =   matrix(1/NCT-(0.1*epsilon/NCT),NCT,1) #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------
p_init_subset =   matrix(1/NCT_subset-(0.1*epsilon/NCT_subset),NCT_subset,1) #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------

#============     Start deconvolving:  ==================

# for ( i in c(1: Num_atlas_files))
  for ( i in c(1: Num_samples))
  {
  print(paste(" ---------- Now deconvolving dataset i = ",as.character(i), " of ", as.character(Num_samples), " ----------------" ) )

  deconv_out_minbound[[i]]        = get_cell_fracs(initial_profile=p_init,        target_profile= Cond_dat[[i]]$ROI_meth_profile,   
                                            Sigmat=Cond_dat[[i]]$Sigmat_whits, epsilon=0.01, mu_in=1e-5, minbound=TRUE, maxbound=TRUE) 
  deconv_subset_out_minbound[[i]] = get_cell_fracs(initial_profile=p_init_subset, target_profile= Cond_dat[[i]]$ROI_meth_profile,   
                                            Sigmat=Cond_dat[[i]]$Sigmat_whits[,CT_subset], epsilon=0.01, mu_in=1e-5, minbound=TRUE, maxbound=TRUE) 
  } ## --- end of for loop running deconvolution on data sets.

par(mfrow=c(1,1))
plot(deconv_out[[i]]$par,deconv_out_2[[i]]$par)


#=========================================================================================
#==============----  SELF CONSISTENCY CHECK  --- =========================================
w_actual_original=c(0.004184101, 0.008368202, 0.016736401, 0.025104603, 0.033472803, 0.041841005, 0.050209204, 
                    0.066945605, 0.083682009, 0.100418411, 0.117154812, 0.133891213, 0.150627615, 0.167364016);

deconv_test = self_consist_check( w_actual_original, as.matrix(Sundat$Sigmat),  noise_level=0.1, metric="methylation", Numruns =10, epsilon=0.01 , lower_bound=TRUE)


labels=c("Lvr","Lng","Cln","sIn","Pnc","Adr","Eso","Adi","Hrt","Brn","Tcl","Bcl","Neu","Plc")
labels_subset = labels[CT_subset]

plot_deconv( deconv_out[[i]]$par, labels, paste("Condition : ", as.character(i)) ) 
plot_deconv( deconv_subset_out[[i]]$par, labels_subset, paste("Condition : ", as.character(i)) ) 

#=========================================================================================
#==============----  TRY TO MAKE SOME PIE CHARTS   ---====================================

Exp_ref=read.table(file.path(PATH_DATA,"info.mix.table"), stringsAsFactors = FALSE, header=TRUE)[,2:4]

  Liver_index=1
  Placenta_index=14
  Placenta_index_subset=5
  
  blood_indices=c(11:13)
  blood_indices_subset=c(2:4)
  
  Liver_index=1
  Placenta_index=14
  blood_indices=c(11:13)

  other_indices = setdiff(1:NCT,c(Liver_index, Placenta_index, blood_indices) )
  lbls <- c("blood cells (B-cell+T-cell+Neutrophils)", "Placenta", "Liver", "Other")
  
  # Simple Pie Chart
  par(mfrow=c(1,2))
  i=3
  cell_frac_calcd         <- c(sum(deconv_out[[i]]$par[blood_indices]), deconv_out[[i]]$par[Placenta_index], deconv_out[[i]]$par[Liver_index] , sum(deconv_out[[i]]$par[other_indices]) ) 
  cell_frac_expected <- cbind(Exp_ref[i,], 0)
  pie(labels=NA,  col=c("Green","red","blue","black"), border="black", cell_frac_calcd,     main="calculated")

  i=9
  cell_frac_calcd_subset  <- c(sum(deconv_subset_out[[i]]$par[blood_indices_subset]), deconv_subset_out[[i]]$par[Placenta_index_subset], deconv_subset_out[[i]]$par[Liver_index] ) 
  cell_frac_expected <- cbind(Exp_ref[i,], 0)
  pie(labels=NA,  col=c("Green","red","blue","black"), border="black", cell_frac_calcd_subset ,     main="calculated")
  pie(labels=NA,  col=c("Green","red","blue","black"), border="black", as.numeric(cell_frac_expected),  main="claimed")


  pie(cell_frac_calcd,
      labels=NA,
      clockwise=TRUE,
      col=c("Green","red","blue","black"),
      border="black",
      radius=0.7,
      cex=0.8,
      main="calculated")

#  legend("bottom", legend=lbls,
#       fill=c("Green","red","blue","brown"), 
#       ncol=1)
  
  
  
   # methods=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",   "SANN", "Brent")
   # #--- for other methods you have to supply the gradient.
   #   for( j in c(1:length(methods)) ) 
   #   {
   #   m=methods[j]
   #   deconv_subset_out[[j]] = constrOptim(w_subset , function(w_subset){func_SQ(w_subset, yEXP_target = Cond_dat[[i]]$ROI_meth_profile,  SigMat_in = Sigmat_whits_CTsubset)},     ui=ui_COND_subset, ci=ci_COND_subset, mu = 1e-04, control = list(), method = m) 
   #   deconv_subset_out[[j]]$par
   #   }
   # 
   # 
   # # }

   # N=matrix(0,NCT,1)
   # for( k in c(1:NCT))
   # {
   #   N[k,1] =  norm(as.matrix(Cond_dat[[i]]$ROI_meth_profile -  Cond_dat[[i]]$Sigmat_whits[,k] ) ,"F")
   # }

