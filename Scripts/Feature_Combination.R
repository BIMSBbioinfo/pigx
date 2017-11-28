# ---------------------------------------------------------------------------- #
# suppressPackageStartupMessages(library(argparse))
#
# parser = ArgumentParser(description="Feature_Combination.R")
# parser$add_argument('--features',   type='character', help='feature', nargs='+')
# parser$add_argument('--annotation', type='character', help='feature', nargs='+')
# parser$add_argument('--bw',         type='character', help='feature', nargs='+')
# parser$add_argument('--scriptdir',  type='character', help='feature')
# parser$add_argument('--outpath',    type='character', help='feature')
# parser$add_argument('--tss.up',     type='character', help='numeric')
# args   = parser$parse_args()

# print(args)


# ---------------------------------------------------------------------------- #
Feature_Combination = function(
  features,
  annotation_path = NULL,
  bw,
  scriptdir,
  outfile,
  tss.up = 1000
){

    suppressPackageStartupMessages(library('genomation'))
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('data.table'))
    suppressPackageStartupMessages(library('stringr'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)


    # ------------------------------------------------------------------------ #
    if(is.null(annotation_path))
        stop('Annotation is not specified')

    message('Reading Annotation ...')
      annotation = readRDS(annotation_path)

    # ------------------------------------------------------------------------ #
    message('Collapsing Peaks ...')
      feat.list = GRangesList(lapply(features, readNarrowPeak))
      lnames = unique(str_replace(basename(unlist(features)),'.(narrow|broad)Peak',''))

      if(length(lnames) != length(feat.list))
          stop('feature_Combination: input features are duplicated')
      names(feat.list) = lnames

      feat.comb = findFeatureComb(feat.list, use.names=TRUE)
      feat.comb$peak_id = paste0('Peak',sprintf('%07d', 1:length(feat.comb)))

      # peak location
    message('Peak Location ...')
      feat.comb$genomic_annotation = suppressWarnings(AnnotateRanges(feat.comb, annotation$genomic_annotation, type='precedence'))

    message('Peak Annotation ...')
      gtf.gens = annotation$gtf
      gtf.gens = unlist(range(split(gtf.gens, gtf.gens$gene_id)))
      fog = dtfindOverlaps(feat.comb, GenomicRanges::resize(gtf.gens, width=width(gtf.gens) + as.numeric(tss.up), fix='end'))
      fog$gene_id = names(gtf.gens)[fog$subjectHits]

      fogm = merge(fog, unique(annotation$annot[,c('gene_id','gene_name','gene_biotype','gcoord')]), by='gene_id')
      fogm = fogm[,c('queryHits','gene_name','gene_id'),with=FALSE][,lapply(.SD, function(x)paste(unique(x), sep=':', collapse=':')),by='queryHits']
      fog$subjectHits=NULL

      feat.comb$gene_name = 'None'
      feat.comb$gene_name[fogm$queryHits] = fogm$gene_name

      feat.comb$gene_id = 'None'
      feat.comb$gene_id[fogm$queryHits]   = fogm$gene_id


    # ------------------------------------------------------------------------ #
    message('Peak Scores...')
    sml = ScoreMatrixList(bw, feat.comb, bin.num=1)
    names(sml) = str_replace(names(sml),'.bw','')
    values(feat.comb) = cbind(values(feat.comb), DataFrame(data.frame(sml)))

    # ------------------------------------------------------------------------ #
    message('Output table...')
    dat = data.frame(coord = as.character(feat.comb),
                     width = width(feat.comb),
                     as.data.frame(values(feat.comb)))

    outfile.tsv  = outfile
    write.table(dat, outfile.tsv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

    outfile.csv  = str_replace(outfile, 'tsv','csv')
    write.table(dat, outfile.csv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

}




# ---------------------------------------------------------------------------- #
# function call
Feature_Combination(
    features         = snakemake@input[['features']],
    annotation_path  = snakemake@input[['annotation']],
    bw               = snakemake@input[['bw']],
    outfile          = snakemake@output[['outfile']],
    scriptdir        = snakemake@params[['scriptdir']],
    tss.up           = snakemake@params[['tss_up']]
)

# feature_Combination(
#     features    =args[['features']],
#     bw          =args[['bw']],
#     annotation  = args[['annotation']],
#     outpath     = args[['outpath']],
#     scriptdir   = args[['scriptdir']]
# )
