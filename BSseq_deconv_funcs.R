#== BSseq_deconv_funcs.R -function definitions necessary for deconvolution script on output from the BSseq pipeline (.cov file) ======
# ---last updated on  Sat Apr 8 23:11:46 CEST 2017  by  blosberg  at location  , brenosbsmacbook.fritz.box

#  changes from  Sat Apr 8 23:11:46 CEST 2017 : Made first draft template: script executes bug free but output inconsistent with expectations. Further testing needed


#==================================================================
#--- plot deconvolution output: ----

plot_deconv <- function( yvals, name_list, Title)
  {
  L=length(yvals)
  
  plot(yvals, xaxt = "n", xlab='cell type', ylab = "cell fraction", main=Title)
  axis(side=1, at=seq(1,L,2), labels=name_list[seq(1,L,2)], cex.axis=0.85)
  axis(side=1, at=seq(2,L,2), labels=name_list[seq(2,L,2)], cex.axis=0.85)
  
  }

#==================================================================
#--- get experimental data ranges: ----

get_Atlas_dat  <- function( ATLAS_WIG_FOLDER, bw_dat_file, bw_rc_dat_file, Biomarks_gr, name="unnamed")
{ 
  #----- imports bigwig files:
bw_rc_gr    =   sort( sortSeqlevels( import(file.path(ATLAS_WIG_FOLDER, bw_rc_dat_file) ) ) )
bw_gr       =   sort( sortSeqlevels( import(file.path(ATLAS_WIG_FOLDER, bw_dat_file) ) ) )

Atlas_dat_gr = subsetByOverlaps(bw_gr, bw_rc_gr) 
names(mcols(Atlas_dat_gr))<-"meth_pc"

temp             = subsetByOverlaps(bw_rc_gr, bw_gr) 
Atlas_dat_gr$cov = temp$score  #----these are now all of the reads with >3 coverage. not necessarily found on a ROI

biomarkers_with_hits_gr = subsetByOverlaps(Biomarks_gr, Atlas_dat_gr) 
L                       = length(biomarkers_with_hits_gr)
if(L!=5820)
  {
  stop("not all the biomarkers in matrix are accounted for")
  }
ROI_meth_profile=matrix(0,L,1)

ranges_hits = findOverlaps(Atlas_dat_gr, Sundat$Biomarks_gr)
ranges_mat  = as.matrix(ranges_hits)


# in this loop:  ranges_mat[,2]==j is true for the block of reads corresponding to the j-th ROI 
# the left column then tells you the corresponding read indexes for that ROI -it should be ordered but NOT comprehensive. Thus: 
# ranges_mat[ ranges_mat[,2]==j ,1] then tells us the list of reads to grab from the total atlas data for this ROI 
# Atlas_dat_gr[ranges_mat[ ranges_mat[,2]==j ,1] ]  <== this is then the GRanges object of reads for the j-th ROI
# "used" just keeps track of which entries were used 

used=matrix(0, length(ranges_hits), 1)

window_min=1
for (j in c(1:L)) #---- now loop through each j biomarker being considered and grab the avg. meth.
  {
  # ROI_meth_profile[j,1] = sum( (Atlas_dat_gr[ranges_mat[ ranges_mat[,2]==j ,1] ]$meth_pc)*(Atlas_dat_gr[ranges_mat[ ranges_mat[,2]==j ,1] ]$cov) ) /sum(Atlas_dat_gr[ranges_mat[ ranges_mat[,2]==j ,1] ]$cov)

  #---- trying to find a faster method:
  temp=which(ranges_mat[,2]==j)  
  window_min=temp[1]
  window_max=tail(temp, n=1)

  ROI_meth_profile[j,1] = sum( (Atlas_dat_gr[ranges_mat[ window_min:window_max ,1] ]$meth_pc)*(Atlas_dat_gr[ranges_mat[ window_min:window_max ,1] ]$cov) ) /sum(Atlas_dat_gr[ranges_mat[ window_min:window_max ,1] ]$cov)

#  if(j %%100 ==0)
#    {  print(paste("--working on ROI ",as.character(j), " out of", as.character(L) ))  }

  used[ window_min:window_max ] =  used[ window_min:window_max ] +1

  #  window_min=window_max+1
  }

if(min(used)!=1 || max(used)!=1)
  {
  stop("at least one entry used a number of times different from 1")
  }
  
# ---the above is equivalent to this:
# range                 = Atlas_dat_gr[ranges_mat[ ranges_mat[,2]==j ,1] ]
# avg_meth_range[j]        = sum( (range$meth_pc)*(range$cov) ) / ( sum(  range$cov) ) 


return( list( "GR" = Atlas_dat_gr, "ROI_meth_profile"=ROI_meth_profile, "name"= name) )
}

#==================================================================
#--- get experimental data ranges: ----
get_Cond_dat <- function(fin, refdat, mincounts ) 
{    
  genome="hg19"
  label=paste("Cond",as.character(i),sep="")
  
  EXP_methraw           = methRead(fin, label, genome, skip=1, pipeline="bismarkCoverage", mincov = 1)  

  #---- now filter only the experimental reads that hit upon the selected biomarkers.
  EXP_methraw_filtered  = regionCounts( EXP_methraw, refdat$Biomarks_gr, strand.aware = FALSE, cov.bases=mincounts) 
  ROI_meth_profile      = EXP_methraw_filtered$numCs/EXP_methraw_filtered$coverage
  
  EXP_gr    = as(EXP_methraw_filtered,"GRanges")  # convert methylraw object to GRange
  EXP_gr    = sortSeqlevels(EXP_gr)                        # sort internal chromosome ordering
  EXP_gr    = sort(EXP_gr)                                 # sort by chromosome level
  #--- GRanges object containing all of the experimental data that hits upon RsOI
  
  biomarkers_with_hits_gr =  subsetByOverlaps(refdat$Biomarks_gr, EXP_gr) 
  #--- GRanges object containing RsOI that were hit upon in this experiment
  
  ovlps         = findOverlaps(biomarkers_with_hits_gr , EXP_gr)
  biomark_hits  = which(countLnodeHits(ovlps)>0)
  #--- list of RsOI that were hit upon
  
  Sigmat_whits      = as.matrix( refdat$Sigmat[biomark_hits,] )                #---- Matrix using just the biomarkers hit on in the experiment, 
  #--- target values for each tissue type are listed by col.
  

  #
  # My old way of getting avg_methylation --produces the same thing as ROI_meth_profile, above
  # the following is clumsy, but explicit; good to keep as a reference:
  #  for (j in c(1:L)) #---- now loop through each j biomarker being considered and grab the avg. meth.
  #  {
  #  range                 = subsetByOverlaps(EXP_gr, refdat$Biomarks_gr[j])
  #  avg_meth_range        = sum( mcols(range)$numCs) / ( sum(mcols(range)$coverage ) ) 
  #  ROI_meth_profile[j,1] = avg_meth_range
  #  }
  

  #----------------------------------
  result <- list("EXP_gr" = EXP_gr, "ROI_meth_profile" = ROI_meth_profile, "biomark_hits"= biomark_hits, "biomarkers_with_hits_gr" = biomarkers_with_hits_gr, "Sigmat_whits" = Sigmat_whits )
  return(result)
}

#==================================================================
#--- get experimental data ranges: ----

get_grange_from_wig <- function( WIG_PATH, wig_file, readcov_file)
{    
  wig    = import.wig(file.path(WIG_PATH, wig_file))
  wigcov = import.wig(file.path(WIG_PATH, readcov_file))
  
#  ?readWig
  
  w2g_out         = c(wig, wigcov)
  w2g_out$score = NULL
  w2g_out = unique(w2g_out)
  
  w2g_out$percmet = 0
  w2g_out$cov     = 0
  
  fow = findOverlaps(w2g_out, wig)
  w2g_out$percmet[queryHits(fow)] = wig$score[subjectHits(fow)]
  
  
  foc = findOverlaps(w2g_out, wigcov)
  w2g_out$cov[queryHits(foc)] = wigcov$score[subjectHits(foc)]
  
  return(w2g_out)
}

#==================================================================
#--- get ranges to compare against: ----

get_refdat <- function(SIGMAT_PATH) 
  {   

  Sun_biomark_file = file.path(SIGMAT_PATH,"Sun_biomarkers_locs_sorted.csv");
  
  Sun_biomark_locs = read.csv(Sun_biomark_file, sep="\t", header=FALSE, stringsAsFactors=TRUE, col.names=c("chrom","start","end"))
  
  Biomarks_gr      = GRanges(seqnames=stringr::str_replace(as.character(Sun_biomark_locs[,1]),' ',''),
                             ranges=IRanges(start=Sun_biomark_locs[,2],end=Sun_biomark_locs[,3]))
  
  Sun_sigmat_file  = file.path(SIGMAT_PATH,"Sun_sigmat_sorted.csv");
  Sigmat           = (1.0/100.0)*read.csv(Sun_sigmat_file, sep=" ", header=TRUE, stringsAsFactors=TRUE)
  
  result <- list("Sigmat" = Sigmat, "Biomarks_gr" = Biomarks_gr)
  return(result)
  }

#==================================================================
#--- DEFINE RESIDUAL FUNCTION TO SEEK A MINIMUM IN ----
func_SQ <- function(w, yEXP_target, SigMat_in) 
{   
  y_calc   = SigMat_in %*% w   
  Res      = (y_calc - yEXP_target)
  
   result = sqrt( sum( Res*Res ) )
  # result = sum(abs( Res ))
  
  return(result)
  
}
#==================================================================
#--- arrays to impose constraints on functional optimization ----
get_conditions <- function(NCT, impose_le_1, impose_ge_1, epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
{  

 if( !impose_le_1 && !impose_ge_1 ) #--- impose no constraints on the sum of w --only that each element is positive.
    {print("getting conditions without any imposition on sum")

    ui_COND         <-  matrix(0,NCT,NCT)
    ci_COND         <-  matrix(0,NCT,1)
    w <- matrix(1/NCT-(0.1*epsilon/NCT),NCT,1) #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------

    diag(ui_COND)   <-  1
    }
else if( impose_le_1 && !impose_ge_1 ) #--- impose constraint that sum(w)<=1 does not exceed one.
  { print("getting conditions with upper bound sum <=1+epsilon, but no lower bound")
  
  ui_COND         <-  matrix(0,NCT+1,NCT)
  ci_COND         <-  matrix(0,NCT+1,1)   #--- CONDITION MATRIX AND VECTOR.
  w <- matrix(1/NCT-(0.1*epsilon/NCT),NCT,1)     #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------

  diag(ui_COND)   <-  1
  ui_COND[NCT+1,] <- -1

  #--- NOW WE IMPOSE THAT ALL FRACTIONS ARE POSITIVE AND SUM TO <= UNITY.
  #--- PERHAPS WE SHOULD ADD A LOWER BOUND TO THE SUM AS WELL @@@ 
  ci_COND[NCT+1] <- -1-epsilon
  }
else if( impose_le_1 && impose_ge_1 ) #--- impose BOTH constraints that 1-eps <= sum(w) <= 1+eps --that sum is within eps of one
  { print("getting conditions with upper and lower bounds 1-epsilon <= sum(w) <= 1+epsilon")
    
    ui_COND         <-  matrix(0,NCT+2,NCT)
    ci_COND         <-  matrix(0,NCT+2,1)   #--- CONDITION MATRIX AND VECTOR.
    w <- matrix(1/NCT-(0.1*epsilon/NCT),NCT,1)     #--- INITIAL CONDITION OF CELL TYPE FRACTIONS ------
    
    diag(ui_COND)   <-  1
    ui_COND[NCT+1,] <- -1
    ui_COND[NCT+2,] <-  1
    
    #--- NOW WE IMPOSE THAT ALL FRACTIONS ARE POSITIVE AND SUM TO <= UNITY.
    #--- PERHAPS WE SHOULD ADD A LOWER BOUND TO THE SUM AS WELL @@@ 
    ci_COND[NCT+1] <- -1-epsilon
    ci_COND[NCT+2] <-  1-epsilon
  }
else if( !impose_le_1 && impose_ge_1 ) #--- impose BOTH constraints that 1-eps <= sum(w) <= 1+eps --that sum is within eps of one
  {
  stop("prompted for conditions with lower bound, but no upper bound. This doesn't make sense. Exiting.")
  }
  
#---------------------------------------------------------------

result <- list("ui_COND" = ui_COND, "ci_COND" = ci_COND, "w" = w)

return(result)
}


#==================================================================
get_dummy_CT_profile  <- function(SigMat_in, SigDEG_in) #=== calculates a "dummy" cell type where everything is "off"
{
  if(min(SigDEG_in<0) || max(SigDEG_in>1))
  {stop("signature digital matrix is not a binary 1/0 matrix")}
  if(!all(dim(SigDEG_in)==dim(SigMat_in) ))  
  {stop("dimensions of analogue and digital signature matrices don't agree.") }
  
  Avg_expression = matrix()
  
  num_CTs  = rowSums(SigDEG_in)   #---  each of these represent an array of genes 
  #---  cataloging how many cell types are "on" or "off" in each case.
  Avg_expression  = (rowSums(SigDEG_in * SigMat_in) ) / num_CTs
  
  return(Avg_expression );
}

