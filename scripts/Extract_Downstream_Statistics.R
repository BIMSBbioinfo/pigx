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
  bam_tag_hist  = NULL,
  outfile       = NULL,
  file_location = NULL,
  scriptdir     = NULL
){

    if(is.null(scriptdir))
        stop('Please specify the script directory')

    # ---------------------------------------------------------------------------- #
    suppressPackageStartupMessages(library('dropbead'))
    suppressPackageStartupMessages(library('yaml'))
    source(file.path(scriptdir,'Functions_Helper.R'), local=TRUE)

    umi_name   = basename(umi_matrix)
    reads_name = basename(reads_matrix)
    stats_name = basename(stats_matrix)
    min.umi    = yaml::yaml.load_file(bam_tag_hist)$reads_cutoff
    stats = extractDownstreamStatistics(
      object        = file_location,
      dge.name      = umi_name,
      reads.by.cell = reads_name,
      min.umi       = min.umi,
      read.statd    = stats_name)

    write.table(stats, file = outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')
}


# ---------------------------------------------------------------------------- #
Extract_Downstream_Statistics(
  umi_matrix    = argv$input[['umi_matrix']],
  reads_matrix  = argv$input[['reads_matrix']],
  reads_stats   = argv$input[['reads_stats']],
  bam_tag_hist  = argv$input[['bam_tag_hist']],
  outfile       = argv$output[['outfile']],
  file_location = argv$params[['file_location']],
  scriptdir     = argv$params[['scriptdir']]

)
