#!/usr/local/bin/Rscript


  suppressPackageStartupMessages(expr = {
    require("AnnotationHub")
    require("rtracklayer")
  })


## this function tries to fetch the reference genes for the given assembly
fetchRefGene <- function(refgenes.loc = NULL,
                         assembly) {
  
  if(is.null(refgenes.loc)) refgenes.loc <- paste0("refseq.genes.",assembly,".bed")
  
  ## import local bed file if available
  if( file.exists(refgenes.loc) ) {
    ## parse it 
    message("Found RefSeq track at:")
    return(refgenes.loc)
    
  } else {
    message("Trying to fetch from AnnotationHub.\n")
    ## else query it from AnnotationHub 
    ah = AnnotationHub()
    ## query refseq genes for assembly
    refseq.q <- query(ah,c("refseq","genes",assembly))
    ## either there is exactly one record, so fetch it
    if(length(refseq.q) == 1) {
      message("Found single RefSeq track, downloading...\n")
      refGenes <- ah[[names(refseq.q)]]
      ## and write it to BED file
      export.bed(object = refGenes,
                 con = refgenes.loc,
                 trackLine=FALSE)
      message("Written the RefSeq track to:")
      return(refgenes.loc)
      
    } else if ( length(refseq.q) == 0 ) { 
      message("Trying to fetch from UCSC table browser.\n")
      ## or there is none, 
      ## so we check with rtracklayer for the latest ucsc data
      mySession = browserSession("UCSC")
      genome(mySession) <- assembly
      track.names <- trackNames(ucscTableQuery(mySession))
      # I am interested in the refGene track 
      if("refGene" %in% track.names) {
        message("Found single RefSeq track, downloading...\n")
        # fetch it as a GRanges object 
        targetTrack <- track(mySession,"refGene")
        ## and write it to BED file
        export.bed(object = targetTrack,
                   con = refgenes.loc,
                   trackLine=FALSE)
        
        message("Written the RefSeq track to:")
        return(refgenes.loc)
      } 
    } else {
      stop(paste("Could not find reference gene set for the given assembly <'",assembly,"'>." ))
    }
    
  }
}


# ## debugging
# save.image(file = "snakemakeObj.RData")

## catch output and messages into log file
out <- file(snakemake@log[[1]], open = "wt")
sink(out,type = "output")
sink(out, type = "message")


## call with snakemake 
fetchRefGene(refgenes.loc      = snakemake@output[["refgenes"]],
             assembly = snakemake@params[["assembly"]])


