# copied from https://gist.github.com/al2na/f002dafb1f3a5782f099#file-ideodmc-r

#' function for making ideogram for differential methylation values
#' requires methylKit, ggbio and GenomicRanges
#'
#' @example
#' library(BSgenome)
#' library("BSgenome.Hsapiens.UCSC.hg18")
#' chr.len = seqlengths(Hsapiens)  # get chromosome lengths
#' # remove X,Y,M and random chromosomes
#' chr.len = chr.len[grep("_|M|X|Y", names(chr.len), invert = T)] 
#' 
#' download.file("http://methylkit.googlecode.com/files/myDiff.rda", 
#'               destfile = "myDiff.rda")
#' load("myDiff.rda")
#' 
#' ideoDMC(myDiff, chrom.length = chr.len, difference = 25, qvalue = 0.01, 
#'        circos = TRUE, title = "test", hyper.col = "magenta", hypo.col = "green")
#'        
ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  hypo = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                        type = "hypo")
  hyper = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                         type = "hyper")
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.per)=(chrom.length)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.po)=(chrom.length)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
      scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
    
  } else {
    
    p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                         aes(x = midpoint, 
                             y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                          hypo.col)) + labs(title = title)
    # new alternative commented out
    #autoplot(c(g.po, g.per), layout = "karyogram", geom = "point", size = 0.65, 
    #aes(x = midpoint,y = meth.diff, color = id))  + scale_colour_manual(values = c(hyper.col, 
    #                                                                                        hypo.col)) + labs(title = title)
    
  }
}

ideoDMC_hyper_hypo <- function(methylDiff.hyper, methylDiff.hypo, chrom.length, 
                               circos = FALSE, title = "Differentially methylated cytosines", hyper.col = "magenta", 
                               hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  hypo = methylDiff.hypo
  hyper = methylDiff.hyper
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.per)=(chrom.length)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.po)=(chrom.length)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
      scale_colour_manual("Regions",values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
    
  } else {
    
    # p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    # p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
    #                      aes(x = midpoint, 
    #                          y = meth.diff, color = id)) + scale_colour_manual("Regions", values = c(hyper.col, 
    #                                                                                                  hypo.col)) + labs(title = title)
    # new alternative commented out
    d = c(g.per, g.po)
    p = autoplot(myIdeo, layout = "karyogram")
    p + layout_karyogram(d, geom = "point", size = 0.65,
                         aes(x = start,  y = meth.diff, color = id))+
      scale_colour_manual("", values = c(hyper.col, hypo.col))+
      labs(title = title) 
  }
}




