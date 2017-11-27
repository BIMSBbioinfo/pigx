# ---------------------------------------------------------------------------- #
feature_Combination = function(
  features,
  annotation = NULL,
  bw,
  scriptdir,
  outpath,
  tss.up = 1000
){

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    message('Reading Annotation ...')
    if(is.null(annotation))
        stop('Annotation is not specified')

    annotation = read_Annotation(annotation)

    if(is.null(annotation$gtf))
        stop('GTF file is not specified')

    message('Collapsing and Annotating Peaks ...')
    feat.list = GRangesList(lapply(features, readNarrowPeak))
    lnames = unique(str_replace(basename(unlist(features)),'.narrowPeak',''))
    if(length(lnames) != length(feat.list))
        stop('feature_Combination: input features are duplicated')
    names(feat.list) = lnames

    feat.comb = findFeatureComb(feat.list, use.names=TRUE)
    feat.comb$peak_id = paste0('Peak',sprintf('%07d', 1:length(feat.comb)))

    annot.l = GTFGetAnnotation(annotation$gtf$gtf, upstream = tss.up)

    # peak location
    feat.comb$location = suppressWarnings(AnnotateRanges(feat.comb, annot.l, type='precedence'))

    # overlapping gene
    gtf.gens = annotation$gtf$gtf
    gtf.gens = unlist(range(split(gtf.gens, gtf.gens$gene_id)))
    fog = dtfindOverlaps(feat.comb, resize(gtf.gens, width=width(gtf.gens) + tss.up, fix='end'))
    fog$gene_id = names(gtf.gens)[fog$subjectHits]

    fogm = merge(fog, unique(annotation$gtf$annot[,c('gene_id','gene_name','gene_biotype','gcoord')]), by='gene_id')
    fogm = fogm[,c('queryHits','gene_name','gene_id'),with=FALSE][,lapply(.SD, function(x)paste(unique(x), sep=':', collapse=':')),by='queryHits']
    fog$subjectHits=NULL

    feat.comb$gene_name = 'None'
    feat.comb$gene_name[fogm$queryHits] = fogm$gene_name

    feat.comb$gene_id = 'None'
    feat.comb$gene_id[fogm$queryHits]   = fogm$gene_id

    if(!is.null(annotation$cpg))
        feat.comb$cpgi = suppressWarnings(countOverlaps(feat.comb, annotation$cpg) > 0)


    # --------------------------------------------------------------- #
    message('Peak Scores...')
    # features = unlist(feat.list)
    # features$class = names(features)
    # fo.score = as.data.table(findOverlaps(feat.comb, features))
    # fo.score[, class:=features$class[fo.score$subjectHits]]
    # fo.score[, score:=features$score[fo.score$subjectHits]]
    # fo.score = fo.score[,list(score = mean(score)), by=list(queryHits, class)]
    # fo.score = dcast(fo.score, queryHits~class, fill=0, value.var='score')
    # fo.score = fo.score[order(fo.score$queryHits)]
    # setnames(fo.score, -1, paste(colnames(fo.score)[-1],'score',sep='_'))
    # values(feat.comb) = cbind(values(feat.comb), DataFrame(as.data.frame(fo.score[,-1,with=FALSE])))
    sml = ScoreMatrixList(bw, feat.comb, bin.num=1)
    names(sml) = str_replace(names(sml),'.bw','')
    values(feat.comb) = cbind(values(feat.comb), DataFrame(data.frame(sml)))

    # --------------------------------------------------------------- #
    message('Output table...')
    dat = data.frame(coord = as.character(feat.comb),
                     width = width(feat.comb),
                     as.data.frame(values(feat.comb)))

    outfile.csv  = file.path(outpath, 'Feature_Combination.csv')
    write.table(dat, outfile.csv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')


    outfile.tsv  = file.path(outpath, 'Feature_Combination.tsv')
    write.table(dat, outfile.tsv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}




# ---------------------------------------------------------------------------- #
# function call
feature_Combination(
    features    = snakemake@input[['features']],
    bw          = snakemake@input[['bw']],
    annotation  = snakemake@params[['annotation']],
    outpath     = snakemake@params[['outpath']],
    scriptdir   = snakemake@params[['scriptdir']]
)
