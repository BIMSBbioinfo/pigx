# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Nucleotide_Frequency.R')

# ---------------------------------------------------------------------------- #
Extract_Nucleotide_Frequency = function(
    genome_fasta    = NULL,
    tilling_windows = NULL,
    outfile         = NULL
){
    
    suppressPackageStartupMessages(library('GenomicRanges'))
    suppressPackageStartupMessages(library('Rsamtools'))
    suppressPackageStartupMessages(library('Biostrings'))

    tiles   = readRDS(tilling_windows)
    fa      = FaFile(file=genome_fasta)
    refbase = getSeq(fa, tiles)
    ofreq   = round(oligonucleotideFrequency(refbase, width = 1) / width(tiles),2)
    dfreq   = round(oligonucleotideFrequency(refbase, width = 2) / width(tiles),2)
    values(tiles) = DataFrame(cbind(ofreq, dfreq))
    
    # ------------------------------------------------------------------------ #
    saveRDS(tiles, file = outfile)
}

# ---------------------------------------------------------------------------- #
# function call
Extract_Nucleotide_Frequency(
    genome_fasta    = argv$input[['genome_fasta']],
    tilling_windows = argv$input[['tilling_windows']],
    outfile         = argv$output[['outfile']]
)
