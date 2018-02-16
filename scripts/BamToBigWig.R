# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')

# ---------------------------------------------------------------------------- #
GRangesToBigWig = function(g, outfile, scale.fac){
    cov = coverage(g)
    cov = round(cov*scale.fac,2)
    rtracklayer::export.bw(cov, outfile)
}

# ---------------------------------------------------------------------------- #
BamToBigWig = function(
    bamfile            = NULL,
    reads_by_cell_file = NULL,
    cell_cutoff_file   = NULL,
    bwfile             = NULL

){
    if(is.null(bamfile))
        stop('bamfile not specified')

    if(is.null(reads_by_cell_file))
        stop('reads_by_cell not specified')

    if(is.null(cell_cutoff_file))
        stop('cell_cutoff not specified')

    if(is.null(bwfile))
        stop('bwfile not specified')


    suppressPackageStartupMessages(library(Rsamtools))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(GenomicFiles))
    suppressPackageStartupMessages(library(GenomicAlignments))
    suppressPackageStartupMessages(library(yaml))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(stringr))

    cell_cutoff   = yaml.load_file(cell_cutoff_file)
    reads_by_cell = read.table(reads_by_cell_file) %>%
        dplyr::filter(V1 >= cell_cutoff[[1]])


    # ---------------------------------------------------------------------------- #
    YIELD = function(bfl) {
        param = ScanBamParam(tag="XC", flag=scanBamFlag(isUnmappedQuery=FALSE))
        aln = readGAlignments(bfl, param=param)
        g   = granges(aln, use.mcols=TRUE)
        g[g$XC %in% reads_by_cell$V2]
    }

    bf = BamFile(bamfile, yieldSize=1000000)
    result = reduceByYield(bf, YIELD, MAP=identity, REDUCE=c)

    if(length(result) == 0)
        stop('All the cells were filtered')

    scale.fac = round(1e8/length(result),3)

    g       = result[strand(result) == '-']
    outfile = str_replace(bwfile,'bw','m.bw')
    GRangesToBigWig(g, outfile, scale.fac)

    g       = result[strand(result) == '+']
    outfile = str_replace(bwfile,'bw','p.bw')
    GRangesToBigWig(g, outfile, scale.fac)

    GRangesToBigWig(result, bwfile, scale.fac)
}

# -------------------------------------------------------------------------- #
BamToBigWig(
    bamfile            = argv$input[['bamfile']],
    reads_by_cell_file = argv$input[['reads_by_cell_file']],
    cell_cutoff_file   = argv$input[['cell_cutoff_file']],
    bwfile             = argv$output[['bwfile']]
)
