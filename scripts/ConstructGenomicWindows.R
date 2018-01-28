# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('ConstructGenomicWindows.R')

# ---------------------------------------------------------------------------- #
ConstructGenomicWindows = function(
    chrlen_path = NULL,
    tilewidth   = 10000,
    outfile     = NULL
){
    
    suppressPackageStartupMessages(library('GenomicRanges'))
    if(is.null(chrlen_path))
        stop('Chromosome lenghts file is not defined')
    
    chrlen = read.table(chrlen_path, header=FALSE)
    tiles = unlist(tileGenome(unlist(split(chrlen[,2], chrlen[,1])), tilewidth=tilewidth))
    
    # ------------------------------------------------------------------------ #
    saveRDS(tiles, file = outfile)
}

# ---------------------------------------------------------------------------- #
# function call
ConstructGenomicWindows(
    chrlen_path    = argv$input[['infile']],
    tilewidth      = argv$params[['tilewidth']],
    outfile        = argv$output[['outfile']]
)
