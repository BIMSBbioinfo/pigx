# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Downstream_Statistics')


# ---------------------------------------------------------------------------- #
#' Extract_Downstream_Statistics - extracts the statistics about the PCR duplication rates
#' the signal around the peak regions
#'
#' @return saves a tsv file with the statistics
Extract_Downstream_Statistics = function(
  umi_matrix    = NULL,
  reads_matrix  = NULL,
  reads_stats   = NULL,
  reads_cutoff  = NULL,
  outfile       = NULL,
  file_location = NULL
){

  
    # ------------------------------------------------------------------------ #
    suppressPackageStartupMessages(library('yaml'))
    suppressPackageStartupMessages(library('data.table'))
    
    # ------------------------------------------------------------------------ #
    min.umi    = yaml::yaml.load_file(reads_cutoff)$reads_cutoff
   
    dge        = data.frame(fread(umi_matrix), row.names = 1)
    
    reads.cell = data.frame(fread(reads_matrix), row.names = 1)
    
    read.stats = read.table(reads_stats, header=TRUE)
    
    # ------------------------------------------------------------------------ #
    dge_colsum = colSums(dge)
    dge.ind    = which(dge_colsum >= min.umi)
    dge        = dge[, names(dge)[dge.ind]]
    reads.cell = reads.cell[reads.cell[, 2] %in% names(dge)[dge.ind], ]
    
    df = data.frame('cells'   = dim(dge)[2],
                    'reads'   = round(median(reads.cell[, 1])),
                    'genes'   = round(median(colSums(dge[,names(dge)[dge.ind]] > 0))),
                    'umis'    = round(median(dge_colsum)),
                    'sum.umi' = round(sum(dge)/10^6, 1),
                    'PCR'     = round(median(reads.cell[, 1]/dge_colsum), 1))
    df$PCR = round(median((reads.cell[, 1]*(1-sum(read.stats[, 3:6])/read.stats[, 2]))/dge_colsum), 1)
    rownames(df) = sub('_UMI.Matrix.txt','',basename(umi_matrix))

    # ------------------------------------------------------------------------ #
    write.table(df, file = outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')
}


# ---------------------------------------------------------------------------- #
Extract_Downstream_Statistics(
  umi_matrix    = argv$input[['umi_matrix']],
  reads_matrix  = argv$input[['reads_matrix']],
  reads_stats   = argv$input[['reads_stats']],
  reads_cutoff  = argv$input[['reads_cutoff']],
  outfile       = argv$output[['outfile']],
  file_location = argv$params[['file_location']]
)

