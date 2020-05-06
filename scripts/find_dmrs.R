.calcFisherStat <- function(pval, jitter = 1e-16){
    if (any(pval==0)){
        pval[pval==0] = jitter
    }
    return(sum(-2*log(pval)))
}

#' @importFrom stats pchisq
.calcFisherPval <- function(pval, ...){
    stat <- .calcFisherStat(pval, ...)
    return(1-pchisq(stat, 2*length(pval)))
}

#' @importFrom stats pnorm qnorm
.calcStoufferPval <- function(pval, effect.size, jitter = 1e-16){
    pval[pval==1] <- 1-jitter
    pval[pval==0] <- jitter
    tmp <- qnorm(1-pval/2)*sign(effect.size)
    stat <- sum(tmp)/sqrt(length(tmp))
    return((1-pnorm(abs(stat)))*2)
}

#' @importFrom stats pnorm qnorm
.calcStoufferPvalOneSided <- function(pval, jitter = 1e-16){
    pval[pval==1] <- 1-jitter
    tmp <- qnorm(1-pval)
    stat <- sum(tmp)/sqrt(length(tmp))
    return(1-pnorm(stat))
}

# 
joinSegmentNeighbours <- function(res) {
    if (length(unique(seqnames(res))) > 1) {
        gr <- lapply(split(res, seqnames(res)), joinSegmentNeighbours)
        gr <- do.call(c, unlist(gr, use.names = FALSE))
        return(gr)
    }  else if (length(res) <= 1) {
        return(res)
    }   else {
        res_dt <- data.table::copy(data.table::as.data.table(res))
        ID = endRow = num.mark = seg.group = seg.mean = startRow = NULL
        group.neighbours <- rle(res$seg.group)
        N = length(group.neighbours$lengths)
        if (N == 1) {
            res_dt[, `:=`(
                seqnames = unique(seqnames),
                start = min(start),
                end = max(end),
                strand = "*",
                width = sum(as.numeric(width)),
                ID = unique(ID),
                num.mark = sum(as.numeric(num.mark)),
                seg.mean = mean(seg.mean),
                startRow = min(startRow),
                endRow = max(endRow),
                seg.group = unique(seg.group)
            )]
        } else {
            k <- numeric(N)
            k[1] <- 1
            l <- group.neighbours$lengths - 1
            for (i in 2:N) {
                k[i] = k[i - 1] + l[i - 1] + 1
            }
            for (i in which(l != 0)) {
                res_dt[k[i]:(k[i] + l[i]), `:=`(
                    seqnames = unique(seqnames),
                    start = min(start),
                    end = max(end),
                    strand = "*",
                    width = sum(as.numeric(width)),
                    ID = unique(ID),
                    num.mark = sum(as.numeric(num.mark)),
                    seg.mean = mean(seg.mean),
                    startRow = min(startRow),
                    endRow = max(endRow),
                    seg.group = unique(seg.group)
                )]
            }
        }
        res_dt <- unique(res_dt)
        return(GenomicRanges::makeGRangesFromDataFrame(res_dt, keep.extra.columns = TRUE))
    }
}

reg <- dmrs[61864]
reg.orig <- selectByOverlap(methylDiffDB,reg)
reg.pval <- reg.orig$pvalue

#' Find regions of differential methylation
#' 
#' This function performs unsupervised segmentation of 
#' the methylation differences calculated from calculateDiffMeth.
#' Then it classifies hyper and hypo methylated segments based 
#' on a methylation difference threshold and joins 
#' neighbouring segments with the same class. 
#'
#' @param methylDiff methylDiff or methylDiffDB object generated with calculateDiffMeth
#' @param meth.diff integer methylation difference threshold
#'
#' @return A GRanges object with segment classification and information. 
#' 'seg.mean' column shows the mean methylation per segment. 
#' 'seg.group' column shows the segment groups obtained by mixture modeling
#' @export
#'
#' @examples
find_DMR <- function(methylDiff, meth.diff = 25, regionStatMethod=c("Stouffer","Fisher")) {
    
    # segment first
    mDiffRegion <- methSeg(methylDiff)
    
    # assign clusters based on meth.diff cutoff
    mDiffRegion$seg.group <- "none"
    if( any(mDiffRegion$seg.mean < -meth.diff)) {
        mDiffRegion[mDiffRegion$seg.mean < -meth.diff]$seg.group <- "hypo"
    }
    if( any(mDiffRegion$seg.mean > meth.diff)) {
        mDiffRegion[mDiffRegion$seg.mean > meth.diff]$seg.group <- "hyper"
    }
    # join neighbouring with same groups
    mDiffRegion <- joinSegmentNeighbours(mDiffRegion)
    
    regionStatFun <- switch (regionStatMethod[1],
                             Stouffer = .calcStoufferPvalOneSided,
                             Fisher = .calcFisherPval
    )
    
    mDiffRegion$pvalue <- sapply(seq_along(mDiffRegion), 
               FUN = function(i,.fun) { 
                    mDiff = selectByOverlap(methylDiff,mDiffRegion[i])
                    return(.fun(mDiff$pvalue))
           },
           .fun = regionStatFun
           )
    
    mDiffRegion$qvalue <- p.adjust(p =mDiffRegion$pvalue, method = "fdr")
        
    
    return(mDiffRegion)
    
}
