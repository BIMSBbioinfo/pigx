##-------------- QP_celltypes_deconv.R - By Bren Osberg, MDC Berlin, started October 2016.
#-- last updated on never.

##--- Analyzes genetic expression patterns and uses a matrix of expression trends from known cells to 
##--- infer what fraction of the collected sample originates from each cell type. Initially run using 
##--- artifically constructed data as a proof of concept.

## ==================================================================================================


# CLEAN UP EVERYTHING:

self_consist_check <- function( w_actual_original, SigMat, noise_level, metric="methylation", Numruns =10, epsilon=0.01, lower_bound=FALSE)
{
  
temp=dim(SigMat)
NG =temp[1] #-- should be 547 for Newman, 1013 for Sun
NCT=temp[2] #-- should be 22  for Newman, 14 for Sun

xmin=0; xmax=1.01*max( w_actual_original );
ymin=0; ymax=1.01*max( w_actual_original );

#--- SET UP ARTICIFIAL DATA AS A TEST---------------------------------------

alpha    = noise_level; #---- DEFINE NOISE LEVELS

w_actual = list()
ymeas    = list()

ynoisy  = list()

# y1=sample(w_actual_original), y1=sample(w_actual_original),y1=sample(w_actual_original)
  for (i in c(1:Numruns)){
        w_actual[i] = list(sample(w_actual_original))
        ymeas[i]    = list(SigMat %*% w_actual[[i]]) #--- "ACTUAL" cell expression data that would be measured exactly in this hypothetical example.

        ynoisy[i]   = list( abs(rnorm(NG, 0, alpha*sqrt(ymeas[[i]] ))+ymeas[[i]] ) )
        }

#========== PLOT GENE EXPRESSION NOISE ===============================
# par(mfrow=c(1,1))
# plot(ymeas[[1]], ynoisy1[[1]], main = "blurring level of expression data", xlab=paste("actual ",metric), ylab=paste("noisy ",metric, " supplied to algorithm"), col="red", lwd=1)
# for (i in c(1:Numruns)){ points(ymeas[[i]],ynoisy1[[i]],col="red",pch=1,lwd=1);  points(ymeas[[i]],ynoisy2[[i]],col="blue",pch=0,lwd=1); points(ymeas[[i]],ynoisy3[[i]],col="green",pch=2,lwd=1) }
# lines(c(0,100),c(0,100),col="black",lw=4)
# legend("topleft",inset=.05, title="Noise levels", c(paste("a=",alpha1), paste("a=",alpha2), paste("a=",alpha3)), horiz=TRUE, lty=c(1,1), lwd=c(1,1,1),  col=c("red","blue","green"),pch=c(1,0,2)) 
#---------------------------------------------------------------------------
#--- NOW SET UP WORKING PARAMETERS

conditions <- get_conditions(NCT, 1, as.integer(lower_bound), epsilon ) # DEFINITION: return ui_COND, ci_COND, w, s.t. ui_COND * w >= ci_COND
ui_COND = conditions$ui_COND; ci_COND = conditions$ci_COND; w = conditions$w


Results=list()

for (i in c(1:Numruns)){
    w <- matrix(1/NCT-(0.01/NCT),NCT,1) 
    Results[i] = list(constrOptim(w , function(w){func_SQ(w, yEXP_target = ynoisy[[i]],  SigMat_in= SigMat)}, ui=ui_COND, ci=ci_COND, mu = 1e-04, control = list(), method = "Nelder-Mead"))
    
    print(paste("completed run",i," of ", Numruns))
    }

#========== PLOT DECONVOLUTION OUTPUT ================================
plot(w_actual[[1]], Results[[1]]$par[1:NCT], main = "Cell fractions: input vs. deconvoluted ", xlab="actual fraction", ylab="deconvolution output", col="red", lwd=1, xlim=c(xmin,xmax), ylim=c(ymin,ymax))
lines(c(0,1),c(0,1),col="black",lwd=5)
for (i in c(2:Numruns)){
  points(w_actual[[i]], Results[[i]]$par[1:NCT],col="red",lwd=1)
}
legend("topleft",inset=.05, title="Noise level", paste("a=",alpha) , horiz=TRUE, lty=c(1,1), lwd=c(5,3,3),  col=c("red"), pch=c(1) ) 
lines(c(0,1),c(0,1),col="black",lwd=5)


return(list(w_actual, Results ))

}

# # ===== BELOW ARE SOME SNIPPETS THAT WERE USED TO SHOW THE NOISINESS LEVEL.
# #======================================================================
# 
# rot_temp=list()
# for (n in c(1:Numruns)){
#   rot_temp[[n]]=matrix(0,NCT,NCT)
#   }
# 
# for (n in c(1:Numruns))
#   {
#   rot_temp[[n]]= get_rot_mat(w_actual_original,w_actual[[n]])
#   }
# 
# # ----------------------------------------
# for (n in c(1:Numruns)){
#   print(paste("for value n=", n,", max val=",max(  (rot_temp[[n]]%*%w_actual[[n]] -w_actual_original) )) )
#   }
# # ----------------------------------------
# 
# #--- rot_temp is now definitely a transfer matrix that rotates the used data into the basis of the _original_ data.
# 
# #======================================================================
# 
# temp_noisy=matrix(0, NCT, Numruns)
# 
# for (n in c(1:Numruns)){
#     temp_noisy[,n] = (rot_temp[[n]]%*% Results_noisy1[[n]]$par )
# }
# 
# plot_noisy = list()
# 
# for (j in c(1:NCT)){
#   plot_noisy[j] =list(temp_noisy1[j,])
# }
# 
# library(beeswarm)
# beeswarm(plot_clean,  vertical = TRUE, log = TRUE, pch = 16, col = rainbow(8), main = 'components-clean',xlab="initial fraction index",ylab="output fraction")
# 
# beeswarm(plot_noisy1,  vertical = TRUE, pch = 16, col = rainbow(8), main = paste('components, noise\n a=',alpha1) ,xlab="initial fraction index",ylab="output fraction")
# beeswarm(plot_noisy2,  vertical = TRUE, pch = 16, col = rainbow(8), main = paste('components, noise\n a=',alpha2) ,xlab="initial fraction index",ylab="output fraction")
# beeswarm(plot_noisy3,  vertical = TRUE, pch = 16, col = rainbow(8), main = paste('components, noise\n a=',alpha3) ,xlab="initial fraction index",ylab="output fraction")
# 
# #beeswarm(plot_noisy3,  vertical = TRUE, log = TRUE, pch = 16, col = rainbow(8), main = paste('components, noise\n a=',alpha3) ,xlab="initial fraction index",ylab="output fraction")
# 
# Numruns_data=10
# Normfac_data=matrix(0,Numruns_data,4)
# for (n in c(1:Numruns_data)){
#   Normfac_data[n,1]  = sum(Results_clean[[n]]$par)
#   Normfac_data[n,2]  = sum(Results_noisy1[[n]]$par)
#   Normfac_data[n,3]  = sum(Results_noisy2[[n]]$par)
#   Normfac_data[n,4]  = sum(Results_noisy3[[n]]$par)
#   }
# 
# dres = data.frame(lapply(lres, function(x)x$par))
