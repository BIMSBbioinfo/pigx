# ---------------------------------------------------------------------------- #
dtfindOverlaps = function(reg1, reg2, ignore.strand=FALSE){

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    as.data.table(findOverlaps(reg1,reg2, ignore.strand=ignore.strand))
}


# ---------------------------------------------------------------------------- #
ReadGTFAnnotation = function(
  gtf.path,

  # defines annotation columns which will be subsetted from the gtf file
  gtf.colnames = c('gene_id','transcript_id','gene_name','gene_biotype'),
    
  # defines the required column names
  required_colnames = c('type','gene_id','transcript_id')
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

        # checks whether the gtf file contains required column names
        if(!all(required_colnames %in% colnames(values(gtf)))){
            missing_colnames = setdiff(required_colnames, colnames(values(gtf)))
            stop('required colnames are missing:', paste(missing_colnames, collapse=','))
        }    

        # checks whether type column contains exon element
        if(!any(gtf$type == 'exon'))
            stop('type column should contain exon element')

        gtf.exon     = subset(gtf, type == 'exon')
        gtf.exon.ge  = split(gtf.exon, gtf.exon$gene_id)
        gtf.range.ge = unlist(range(gtf.exon.ge))
        gtf.exon.tr  = split(gtf.exon, gtf.exon$transcript_id)


    message('Making union...')
        gtf.union = reduce(gtf.exon.ge)

    message('Constructing annotation...')
        gtf.colnames = intersect(gtf.colnames, colnames(values(gtf.exon)))
        if(length(gtf.colnames) == 0)
            stop('.gtf file does not contain required column names: gene_id')
    
        gtf.annot = unique(as.data.frame(DataFrame(gtf.exon))[,gtf.colnames])

        gtf.annot$gcoord  = as.character(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
        gtf.annot$gwidth  = width(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
        gtf.annot$tcoord  = as.character(unlist(range(gtf.exon.tr))[match(gtf.annot$transcript_id, names(gtf.exon.tr))])
        gtf.annot$twidth  = sum(width(gtf.exon.tr))[gtf.annot$transcript_id]


    return(list(gtf       = gtf.exon,
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
GTFGetAnnotation = function(
  g,
  width_params = list(
  tss_width            = 1000,
  tts_width            = 1000,
  tss_wide_width       = 10000,
  tts_wide_width       = 10000,
  tss_body_upstream    = 1000,
  tss_body_downstream  = 10000,
  tts_body_upstream    = 10000,
  tts_body_downstream  = 1000,
  splicing_donor_width = 200,
  splicing_accep_width = 200)
){

    g.exon = g[g$type=='exon']
    exon = g.exon
    values(exon) = NULL

    #------------------------------------------------------------------------- #
    gene = unlist(range(split(g.exon, g.exon$gene_id)))
    tss  = promoters(gene, downstream=width_params$tss_width/2, upstream=width_params$tss_width/2)
    tts  = promoters(resize(gene, width=1, fix='end'), downstream=width_params$tts_width/2,
                    upstream=width_params$tts_width/2)
    intron = GenomicRanges::setdiff(gene, exon)

    #------------------------------------------------------------------------- #
    tss_body = gene
    tss_body = resize(tss_body, width=1, fix='start')
    tss_body = resize(tss_body, width=width_params$tss_body_upstream, fix='end')
    tss_body = resize(tss_body, width=width(tss_body)+width_params$tss_body_downstream, fix='start')

    #------------------------------------------------------------------------- #
    tss_wide = gene
    tss_wide = resize(tss_wide, width = 1, fix='start')
    tss_wide = resize(tss_body, width = width_params$tss_wide_width, fix='center')

    #------------------------------------------------------------------------- #
    tts_body = gene
    tts_body = resize(tts_body, width=1, fix='end')
    tts_body = resize(tts_body, width=width_params$tts_body_upstream, fix='end')
    tts_body = resize(tts_body, width=width(tts_body)+width_params$tts_body_downstream, fix='start')

    #------------------------------------------------------------------------- #
    tts_wide = gene
    tts_wide = resize(tts_body, width = 1,     fix='end')
    tts_wide = resize(tts_wide, width = width_params$tts_wide_width, fix='center')

    #------------------------------------------------------------------------- #
    splicing_donor    = exon
    splicing_donor    = resize(splicing_donor, width = 1,   fix='end')
    splicing_donor    = resize(splicing_donor, width = width_params$splicing_donor_width, fix='center')

    #------------------------------------------------------------------------- #
    splicing_acceptor = exon
    splicing_acceptor = resize(splicing_acceptor, width = 1,   fix='start')
    splicing_acceptor = resize(splicing_acceptor, width = width_params$splicing_accep_width, fix='center')

    gl = list(tss      = tss,
              tts      = tts,
              exon     = exon,
              intron   = intron,
              gene     = gene,
              tss_wide = tss_wide,
              tss_body = tss_body,
              tts_wide = tts_wide,
              tts_body = tts_body,
              splicing_acceptor = splicing_acceptor,
              splicing_donor    = splicing_donor)
    gl = GRangesList(lapply(gl, function(x){values(x) = NULL;x}))
    gl = endoapply(gl, function(x)x[width(x) > 1])

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

# ---------------------------------------------------------------------------- #
#' .listBamFiles
#'
#' @param path - top directory with mapped samples
#' @param suffix - suffix to strip of from bam files
#'
#' @return a data.frame with location of bam files

.list_BamFiles = function(path,suffix='.sorted.bam'){

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(tibble))
    d = tibble(bam_file = list.files(path,
                                         full.names=TRUE,
                                         pattern='bam$',
                                         recursive=TRUE)) %>%
        mutate(bam_name = basename(bam_file))             %>%
        dplyr::filter(str_detect(bam_name,'sorted'))      %>%
        mutate(bam_name = str_replace(bam_name,suffix,''))
    return(d)
}

# ---------------------------------------------------------------------------- #
.list_qsortPeaks = function(path){

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(tibble))
    d = tibble(bed_file = list.files(path, full.names=TRUE, pattern='bed$', recursive=TRUE)) %>%
        mutate(chip_name = basename(bed_file)) %>%
        dplyr::filter(str_detect(chip_name,'sort')) %>%
        mutate(chip_name = str_replace(chip_name,'_qsort.bed',''))
    return(d)
}

# ---------------------------------------------------------------------------- #
# for bowtie2
MappingStats_Bowtie2 = function(path){

    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\t', quiet=TRUE)
    s = str_replace(s,'^ +','')
    s = str_replace(s,' .+','')
    s = str_replace(s,'%','')
    if(length(s) > 6)
        s = s[c(1:5, 15)]
    s = as.numeric(s)
    d = data.table(value = s)
    d$stat = c('reads.total','reads.unpaired','reads.unmapped','reads.uniq','reads.mult','alignment.rate')
    d = rbind(d, data.table(
        stat  ='mapped.total',
        value = subset(d,stat=='reads.uniq')$value + subset(d,stat=='reads.mult')$value))
    return(d)
}
