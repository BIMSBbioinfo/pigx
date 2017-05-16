#!/usr/bin/env Rscript  ##== BSseq_deconv -deconvolves BSseq pipeline output from args (takes .bam file) into cell fractions ======
# ---last updated on May 5th.
#
# By Bren Osberg, MDC Berlin, started October 2016.
#==================================================================================================

# install.packages("devtools")
# install.packages("limSolve")
# install.packages("stringr", repos='http://cran.us.r-project.org')

# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# biocLite("methylKit")
# library("devtools")
# install_github("al2na/methylKit", build_vignettes=FALSE, repos=BiocInstaller::biocinstallRepos(),ref="1.1.7BioC3.4", dependencies=TRUE)
rm(list=ls()) # CLEAN UP EVERYTHING

# install.packages("methylKit")
# methylKit::methRead() 

library(GenomicRanges)
library(stringr)
library(methylKit)
library(rtracklayer)
library(limSolve)

CLargs = commandArgs(trailingOnly=TRUE)
if(length(CLargs) != 2)
  {
  print("command line arguments supplied are: ")
  for(i in length(CLargs))
    {
    print(CLargs[i])
    }
  
  stop("Incorrect number of command line arguments supplied. Exiting.")
  }


filename=CLargs[1]
sampleID=CLargs[2]


source('./Define_BSvars.R'); # ==== this script should contain the following definitions:
# WORKDIR          ="/home/bosberg/bs/pigx_bsseq"
# SIGMAT_PATH      ="/home/bosberg/bs/pigx_out/HK_Sun_data"
# (this one is omitted): ATLAS_WIG_FOLDER ="/data/akalin/bosberg/genome_atlas"
# PATH_DATA        ="/data/akalin/bosberg/Sun_PNAS_processed_output/MIX_FINAL_covfiles"
# setwd(WORKDIR)

source('./BSseq_deconv_funcs.R');

#===============================================================
# ----- GET THE SIGNATURE MATRIX AND MARKER LOCATIONS  ---------
refdat   = get_refdat(SIGMAT_PATH)
NCT_full = dim(refdat$Sigmat)[2];
  #------------------------------------
  # CT_subset=c(1, 11, 12, 13, 14)
  # NCT_subset =length(CT_subset)
  # refdat_s = refdat
  # refdat_s$Sigmat = refdat_s$Sigmat[,CT_subset]


#================================================================
#------------     IMPORT  YOUR SAMPLE DATA:      ----------------
mincounts  = 1;  

# fin=paste0(PATH_DATA,filename)
fin=filename #--- snakemake just feeds the whole filename into the command line.

type=".bam"

print(paste("Importing experimental data"))

Exp_dat = get_Exp_dat(fin, refdat=refdat , mincounts, filetype=type, ID=sampleID, genome="hg19")  

#================================================================
#---------------------    DECONVOLVE    -------------------------

epsilon=0.01  #---- tolerance parameter for the sum of cell-type fractions.
mu=1e-5       #---- tolerance parameter for Nelder mead convergence.

print(paste("Generating matrix conditions"))
p_init       =   matrix(1/NCT_full-(0.1*epsilon/NCT_full),NCT_full,1) #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------
conditions   =  get_conditions(NCT_full, 1, 0, epsilon )

print(paste("finished getting matrix conditions"))
#============     Start deconvolving:  ==================

E_in = matrix(1,1,NCT_full)
F_in = 1

# ===== starting my own convolution     =====
deconv_out = get_cell_fracs(initial_profile=p_init,        target_profile= Exp_dat$ROI_meth_profile,   
                                          Sigmat=Exp_dat$Sigmat_whits, epsilon=epsilon, mu_in=mu, minbound=FALSE, maxbound=TRUE) 

# ===== starting Peiyongs deconvolution =====
# this is the function recommended by Peiyong
Peis_fracs = lsei(A = Exp_dat$Sigmat_whits, B =  Exp_dat$ROI_meth_profile, E=E_in, F=F_in, 
                         G = conditions$ui_COND, H = conditions$ci_COND)
  
print(paste("saving data"))
save(deconv_out, Peis_fracs, refdat, Exp_dat, sampleID, file = paste0(PATHOUT, sampleID, "deconv_out.RData") )

print(paste("BSseq_deconv.R program complete."))
# par(mfrow=c(1,1))
