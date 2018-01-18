# ---------------------------------------------------------------------------- #
dtfindOverlaps = function(reg1, reg2, ignore.strand=FALSE){

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    as.data.table(findOverlaps(reg1,reg2, ignore.strand=ignore.strand))
}

# ---------------------------------------------------------------------------- #
GtfSelectTranscript = function(gtf){

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(stringr))
    if(!is.null(gtf$exon_id)){
        gtf = gtf[gtf$exon_id != 'CCDS']
        gtf = gtf[!str_detect(gtf$exon_id, 'mRNA')]
        gtf = gtf[!str_detect(gtf$exon_id, 'cdna')]
    }
   # gtf = keepStandardChromosomes(gtf, pruning.mode='coarse')

    d = data.table(as.data.frame(values(gtf)))
    d$strand = as.character(strand(gtf))
    d$width = width(gtf)
    d[d$strand == '+' , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == '+']]
    d[d$strand == '-' , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == '-']]
    d.cnt = d[, c('gene_id','transcript_id','COUNT','width'), with=F]
    d.cnt = d[, list(COUNT=unique(COUNT), width=sum(width)), by=c('gene_id', 'transcript_id')]
    d.cnt = d.cnt[order(-d.cnt$COUNT, -d.cnt$width),]
    d.cnt = d.cnt[!duplicated(d.cnt$gene_id)]

    gtf.t = gtf[gtf$transcript_id %in% d.cnt$transcript_id]
    return(gtf.t)
}


# ---------------------------------------------------------------------------- #
ReadGTFAnnotation = function(
  gtf.path
){

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(genomation))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(rtracklayer))


    if(!file.exists(gtf.path))
        stop('the gtf file does not exist')

    message('Importing gtf...')
    gtf          = import.gff(gtf.path)
    gtf.exon     = subset(gtf, type == 'exon')
    gtf.exon.ge  = split(gtf.exon, gtf.exon$gene_id)
    gtf.range.ge = unlist(range(gtf.exon.ge))
    gtf.exon.tr  = split(gtf.exon, gtf.exon$transcript_id)


    message('Selecting transcripts...')
    gtf.selected = GtfSelectTranscript(gtf.exon)
    gtf.selected = split(gtf.selected, gtf.selected$transcript_id)

    message('Making union...')
    gtf.union = reduce(gtf.exon.ge)

    message('Constructing annotation...')
    gtf.annot = unique(as.data.frame(DataFrame(gtf.exon))[,c('gene_id','transcript_id','gene_name','gene_biotype')])

    gtf.annot$gcoord  = as.character(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
    gtf.annot$gwidth  = width(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
    gtf.annot$tcoord  = as.character(unlist(range(gtf.exon.tr))[match(gtf.annot$transcript_id, names(gtf.exon.tr))])
    gtf.annot$twidth  = sum(width(gtf.exon.tr))[gtf.annot$transcript_id]


    return(list(gtf       = gtf.exon,
                gtf.sel   = gtf.selected,
                gtf.union = gtf.union,
                annot     = gtf.annot))
}


# ---------------------------------------------------------------------------- #
setGeneric("AnnotateRanges",
           function(region, annotation,
                    ignore.strand = FALSE,
                    type          = 'precedence',
                    null.fact     = 'None',
                    collapse.char = ':',
                    precedence    = NULL,
                    id.col        = NULL)
               standardGeneric("AnnotateRanges") )

setMethod("AnnotateRanges",signature("GRanges","GRangesList"),
          function(
              region, 
              annotation, 
              ignore.strand = FALSE, 
              type          = 'precedence', 
              null.fact     = 'None',
              collapse.char = ':'
            
            ){

              if(! class(region) == 'GRanges')
                  stop('Ranges to be annotated need to be GRanges')

              if(! all(sapply(annotation, class) == 'GRanges'))
                  stop('Annotating ranges need to be GRanges')

              if(!type %in% c('precedence','all'))
                  stop('type may only be precedence and all')

              suppressPackageStartupMessages(library(data.table))
              suppressPackageStartupMessages(library(GenomicRanges))
              cat('Overlapping...\n')
              if(any(names(is.null(annotation))))
                  stop('All annotations need to have names')

              if(class(annotation) != 'GRangesList')
                  annotation = GRangesList(lapply(annotation, function(x){values(x)=NULL;x}))

              a = suppressWarnings(data.table(as.matrix(findOverlaps(region, annotation, ignore.strand=ignore.strand))))
              a$id = names(annotation)[a$subjectHits]
              a$precedence = match(a$id,names(annotation))
              a = a[order(a$precedence)]

              if(type == 'precedence'){
                  cat('precedence...\n')
                  a = a[!duplicated(a$queryHits)]
                  annot = rep(null.fact, length(region))
                  annot[a$queryHits] = a$id
              }
              if(type == 'all'){
                  cat('all...\n')
                  a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
                  annot = rep(null.fact, length(region))
                  annot[a$queryHits] = a$id

              }
              return(annot)

})

# ---------------------------------------------------------------------------- #
GTFGetAnnotation = function(g, downstream=500, upstream=1000){

    g.exon = g[g$type=='exon']
    exon = unlist(g.exon)
    gene = unlist(range(split(g.exon, g.exon$gene_id)))
    tss  = promoters(gene, downstream=downstream, upstream=upstream)
    tts  = promoters(resize(gene, width=1, fix='end'), downstream=downstream,
                    upstream=upstream)
    intron = GenomicRanges::setdiff(gene, exon)

    values(exon) = NULL
    gl = GRangesList(tss    = tss,
                     tts    = tts,
                     exon   = exon,
                     intron = intron)

    return(gl)

}

# # ---------------------------------------------------------------------------- #
# annotates a bam file with a given annotation list
Annotate_Reads = function(
    infile        = NULL,
    annotation    = NULL,
    ignore.strand = FALSE
){

    suppressPackageStartupMessages(library(GenomicAlignments))
    reads = readGAlignments(infile, use.names=TRUE, param=ScanBamParam(which=w, tag='NH'))
    g = granges(reads, use.names=TRUE, use.mcols=TRUE)
    if(length(g) == 0)
        return(data.table(rname=NA, annot=NA, uniq=NA))

    g$annot = AnnotateRanges(g, annotation, ignore.strand=ignore.strand)
    g = g[order(match(g$annot, c(names(annotation),'None')))]
    g$uniq  = factor(ifelse(g$NH == 1,'Uniq','Mult'),levels=c('Uniq','Mult'))
    dg = as.data.table(values(g)[,c('annot','uniq')])
    dg$rname = names(g)
    dg = dg[!duplicated(dg$rname)]

    ldg = rbindlist(lchr)
    ldg = ldg[order(match(ldg$annot, c(names(annotation),'None')))]
    ldg = ldg[!duplicated(ldg$rname)]
    ldg = na.omit(ldg)

    sdg = data.table(experiment = BamName(infile),
                     ldg[,list(cnts=length(rname)), by=list(annot,uniq)])

    sdg[,freq:=round(cnts/sum(cnts),2)]
    return(sdg)
}
