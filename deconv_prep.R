##------- genomic ranges management script; preps BSseq pipeline output for cell deconvolution -By Bren Osberg, MDC Berlin, started October 2016.
#-- last updated on never.
#
##--- Does some stuff to be elaborated on later...
## ==================================================================================================

# source("https://bioconductor.org/biocLite.R")
# biocLite("methylkit")
# library(devtools)
# install_github("al2na/methylKit", build_vignettes=FALSE, repos=BiocInstaller::biocinstallRepos(),ref="1.1.7BioC3.4", dependencies=TRUE)
# install.packages("GenomicRanges")

rm(list=ls()) # CLEAN UP EVERYTHING

library(methylKit)
library(GenomicRanges)
source('/Users/blosberg/Desktop/Science/Code_library/funcs_general.R')

#=====================================================================
#-------------- IMPORT SIGNATURE MATRIX ------------------------------
SIGMAT_PATH="/Users/blosberg/postdoc_work/bs/data/ref/"
Sun_biomark_file=paste(SIGMAT_PATH,"Sun_biomarkers_locs_sorted.csv", sep="");

# biomarkers_type1 = read.csv(paste(SIGMAT_PATH,"Sun_pnas_biomarker_1.csv", sep=""), sep="\t", header=TRUE, stringsAsFactors=TRUE)
# biomarkers_type2 = read.csv(paste(SIGMAT_PATH,"Sun_pnas_biomarker_2.csv", sep=""), sep="\t", header=TRUE, stringsAsFactors=TRUE)

Sun_biomark_locs=read.csv(Sun_biomark_file, sep="\t", header=FALSE, stringsAsFactors=TRUE, col.names=c("chrom","start","end"))

Sun_sigmat_gr=GRanges(seqnames=stringr::str_replace(as.character(Sun_biomark_locs[,1]),' ',''),
                      ranges=IRanges(start=Sun_biomark_locs[,2],end=Sun_biomark_locs[,3]))


#================================================================================
# ----  NOW GET THE AVERAGE METHYLATION OF THE SAMPLE AT THESE LOCATIONS --------


#----------- READ IN YOUR SAMPLE DATA:
PATH_DATA   ="/Users/blosberg/postdoc_work/bs/data/";

file7_cov     ="Cond7.read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov";
Cond7_methraw = methRead( paste(PATH_DATA, file7_cov, sep=""),"cond1","hg19", skip=1, pipeline="bismarkCoverage") 
Cond7_gr      = as(Cond7_methraw,"GRanges") # convert methylraw object ot GRanges
Cond7_gr      = sortSeqlevels(Cond7_gr) # sort internal chromosome ordering

seqlevels(Cond7_gr) # check chromosome ordering
Cond7_gr <- sort(Cond7_gr) # sort by chromosome level

#------ 

file8_cov     ="Cond8.read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov";
Cond8_methraw = methRead( paste(PATH_DATA, file7, sep=""),"cond1","hg19", skip=1, pipeline="bismarkCoverage") 
Cond8_gr      = as(Cond8_methraw,"GRanges") # convert methylraw object to GRanges
Cond8_gr      = sortSeqlevels(Cond8_gr) # sort internal chromosome ordering

seqlevels(Cond8_gr) # check chromosome ordering
Cond8_gr <- sort(Cond8_gr) # sort by chromosome level


#----@@ still need to figure out how to do it with zipped files; unfortunately OSX insists on appending a ".Z" to my filenames.

Cond7_gr[Cond7_gr %in% Sun_sigmat_gr] #--- find overlap between two GRanges objects
selectByOverlap(Cond7,Sun_sigmat_gr) #--- takes 1 methylraw object and 1 GRanges object

findOverlaps(Cond7_gr,Sun_sigmat_gr)
tmp = subsetByOverlaps(Cond7_gr,Sun_sigmat_gr)
sum(mcols(tmp)$coverage * mcols(tmp)$numCs) / ( sum(mcols(tmp)$coverage * (mcols(tmp)$numCs + mcols(tmp)$numTs) ) ) 

ovlps=findOverlaps(Cond7_gr,Sun_sigmat_gr)
subsetByOverlaps(Cond7_gr,Sun_sigmat_gr[2])
Sun_sigmat_gr[3]

#===============================================================================
reduced_bamfile="/Users/blosberg/postdoc_work/bs/data/Cond7.sorted.reduced.bam" 
temp=processBismarkAln(reduced_bamfile,"Cond7","hg19", mincov = 0)

#====================================

file9_cov     ="Cond9.read1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz";
Cond9_methraw = methRead( paste(PATH_DATA, file9_cov, sep=""),"cond1","hg19", skip=1, pipeline="bismarkCoverage") 
Cond9_gr      = as(Cond9_methraw,"GRanges") # convert methylraw object ot GRanges
Cond9_gr      = sortSeqlevels(Cond9_gr) # sort internal chromosome ordering

seqlevels(Cond9_gr) # check chromosome ordering
Cond9_gr <- sort(Cond9_gr) # sort by chromosome level

