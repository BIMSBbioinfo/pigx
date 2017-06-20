# ---------------------------------------------------------------------------- #
# missing function
# GTFGetAnnotation
# AnnotateRanges
# dtfindOverlaps
#

feature_Combination = function(features, annot=NULL, outname, tss.up = 1000){
    
    # modules::import_package('modules',       attach=TRUE)
    # modules::import_package('genomation',    attach=TRUE)
    # modules::import_package('GenomicRanges', attach=TRUE)
    # modules::import_package('readxl',        attach=TRUE)
    # modules::import_package('data.table',    attach=TRUE)
    # helper = modules::import('Functions_Helper')
    
    library('modules')
    library('genomation')
    library('GenomicRanges')
    library('readxl')
    library('data.table')
    source('Functions_Helper.R', local=TRUE)
    
    
    
    if(is.null(annot))
        stop('Annotation is not specified')
    
    annot = read_Annotation(annot)
    
    if(is.null(annot$gtf))
        stop('GTF file is not specified')
    
    feat.list = GRangesList(lapply(features, readNarrowPeak))
    lnames = unique(str_replace(basename(unlist(features)),'.narrowPeak',''))
    if(length(lnames) != length(feat.list))
        stop('feature_Combination: input features are duplicated')
    names(feat.list) = lnames
    
    feat.comb = findFeatureComb(feat.list, use.names=TRUE)
    feat.comb$peak_id = paste0('Peak',sprintf('%07d', 1:length(feat.comb)))
    
    annot.l = GTFGetAnnotation(annot$gtf$gtf, upstream = tss.up)
    
    # peak location
    feat.comb$location = suppressWarnings(AnnotateRanges(feat.comb, annot.l, type='precedence'))
    
    # overlapping gene
    gtf.gens = annot$gtf$gtf
    gtf.gens = unlist(range(split(gtf.gens, gtf.gens$gene_id)))
    fog = dtfindOverlaps(feat.comb, resize(gtf.gens, width=width(gtf.gens) + tss.up, fix='end'))
    fog$gene_id = names(gtf.gens)[fog$subjectHits]
    
    fogm = merge(fog, unique(annot$gtf$annot[,c('gene_id','gene_name','gene_biotype','gcoord')]), by='gene_id')
    fogm = fogm[,c('queryHits','gene_name','gene_id'),with=FALSE][,lapply(.SD, function(x)paste(unique(x), sep=':', collapse=':')),by='queryHits']
    fog$subjectHits=NULL
    
    feat.comb$gene_name = 'None'
    feat.comb$gene_name[fogm$queryHits] = fogm$gene_name
    
    feat.comb$gene_id = 'None'
    feat.comb$gene_id[fogm$queryHits]   = fogm$gene_id
    
    if(!is.null(annot$cpg))
        feat.comb$cpgi = suppressWarnings(countOverlaps(feat.comb, annot$cpg) > 0)

    
    # --------------------------------------------------------------- #
    message('Peak Scores...')
    features = unlist(feat.list)
    features$class = names(features)
    fo.score = as.data.table(findOverlaps(feat.comb, features))
    fo.score[, class:=features$class[fo.score$subjectHits]]
    fo.score[, score:=features$score[fo.score$subjectHits]]
    fo.score = fo.score[,list(score = mean(score)), by=list(queryHits, class)]
    fo.score = dcast(fo.score, queryHits~class, fill=0, value.var='score')
    fo.score = fo.score[order(fo.score$queryHits)]
    setnames(fo.score, -1, paste(colnames(fo.score)[-1],'score',sep='_'))
    values(feat.comb) = cbind(values(feat.comb), DataFrame(as.data.frame(fo.score[,-1,with=FALSE])))
    
    # --------------------------------------------------------------- #
    message('Output table...')
    dat = data.frame(coord = as.character(feat.comb),
                     width = width(feat.comb),
                     as.data.frame(values(feat.comb)))
    
    outfile.xlsx = file.path(outpath, 'Feature_Combination.xlsx')
    write.xlsx(dat, outfile.xlsx)
    
    outfile.tsv  = file.path(outpath, 'Feature_Combination.tsv')
    write.table(dat, outfile.tsv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}




# ---------------------------------------------------------------------------- #
# function call
feature_Combination(
    features    = snakemake@input[['features']],
    annot  = snakemake@params[['annotation']],
    outname = snakemake@params[['outname']]
)

# annot = list(gtf = '/data/akalin/Base/Annotation/mm9/EnsemblGenes/150520_Ensembl_Mus_musculus.NCBIM37.67.form.gtf',
#              cpg = '/data/akalin/Base/Annotation/mm9/CpGi/mm9.cpg.islands.bed')
# 
# path = '/data/akalin/vfranke/AAkalin_PIX/ChIP'
# features = list()
# features$ChIP_IDR = file.path(path, 'IDR', 'ChIP_IDR/ChIP_IDR.narrowPeak')
# features$Peaks5 = file.path(path, 'Peaks/MACS2/Peaks5/Peaks5_qsort.narrowPeak')
# features$Peaks6 = file.path(path, 'Peaks/MACS2/Peaks6/Peaks6_qsort.narrowPeak')


