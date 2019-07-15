# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')


bamfile = '/clusterhome/vfranke/Projects/AAkalin_pigx/pigx_scrnaseq/tests/out/Mapped/WT_HEK_0h_br1/hg19/WT_HEK_0h_br1.sorted.bam'
sample='WT_HEK_0h_br1'
# ------------------------------------------------------------------------ #
# for bowtie
MappingStats_STAR = function(path, name){
    
    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\n', quiet=TRUE)
    s = as.numeric(str_replace(s[c(5,8,23,25)],'^.+\\t',''))
    d = data.table(sample = name,
                   mapped = c('reads.total','map.uniq','map.mult','map.disc','map.total'),
                   cnts   = c(s,s[2]+s[3]))
    
    d[,freq := round(cnts/cnts[1],3)]
    d$type = 'mapped'
    return(d)
}

# ------------------------------------------------------------------------ #
Star_Solo_Stat = function(path, name, type=''){
    
   s = scan(path, what='character', sep='\n', quiet=TRUE)
   s = s[-c(1,6)]
   s = gsub('\\s+','\t',s)
   s = strsplit(s,'\t')
   d = do.call(rbind, lapply(s, function(x)data.frame(t(x), stringsAsFactors = FALSE)))
   d$X1 = NULL
   d$X3 = as.numeric(as.character(d$X3))
   d = d[c(2,3,4,6,7,13,14),]
   d = data.frame(sample = name, d)
   colnames(d) = c('sample','mapped','cnts')
   d$freq = NA
   d$type = type
   return(d)
   
}

# -------------------------------------------------------------------------- #
Extract_Read_Statistics = function(
    bamfile  = NULL,
    outfile  = NULL,
    sample   = 'Sample',
    mito_chr = 'chrM'

){
    if(is.null(bamfile))
        stop('bamfile not specified')

    if(is.null(outfile))
        stop('outfile not specified')

    suppressPackageStartupMessages({
      library(Rsamtools)
      library(data.table)
    })

    basedir = dirname(bamfile)
    
    message('STAR Mapping Statistics ...')
    stat_file     = file.path(basedir, paste(sample, 'Log.final.out', sep='_'))
    stats_mapping = MappingStats_STAR(stat_file, sample)
  
    message('% Mapping to mitochondrial genome ...')
    targets = scanBamHeader(bamfile)[[1]]$targets
    mito_count = NA
    if(mito_chr %in% names(targets)){
        mito_count = countBam(
            bamfile, 
            param = ScanBamParam(which = GRanges(mito_chr, IRanges(1, targets[mito_chr])))
        )
    }
    stats_mito = data.frame(
        sample = sample, 
        mapped = 'mito_count',
        cnts   = mito_count,
        freq   = NA,
        type   = 'mito'
    )
  
    message('Exon UMI statistics ...')
    solo_path = file.path(basedir, paste(sample,'Solo.out', sep='_'))
    
    exon_path  = file.path(solo_path, 'Gene.stats')
    stats_exon = Star_Solo_Stat(exon_path, sample, 'exon')
    
    message('Gene UMI statistics ...')
    gene_path  = file.path(solo_path, 'GeneFull.stats')
    stats_gene = Star_Solo_Stat(gene_path, sample, 'gene')
  
    message('Spliced statistics ...')
    sj_path  = file.path(solo_path, 'SJ.stats')
    stats_sj = Star_Solo_Stat(sj_path, sample, 'spliced')

    message('Writing read statistics ...')
    read_statistics = rbindlist(list(
        stats_mapping, 
        stats_mito,
        stats_exon,
        stats_gene,
        stats_sj
    ))
    read_statistics$freq = round(read_statistics$cnts / subset(read_statistics, mapped == 'reads.total')$cnts,2)

    write.table(read_statistics, outfile,
        row.names=TRUE, col.names=TRUE,
        sep='\t',quote=FALSE)
}

# -------------------------------------------------------------------------- #
Extract_Read_Statistics(
      bamfile  = argv$input[['bamfile']],
      outfile  = argv$output[['outfile']],
      sample   = argv$params[['sample']],
      mito_chr = argv$params[['mito_chr']]
  )
