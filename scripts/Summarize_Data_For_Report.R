# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Summarize_Data_For_Report')


#
# config = yaml::yaml.load_file('Tests/sample_sheet.yaml')
# scriptdir = 'scripts'
# path_mapped = 'Tests/out/Mapped'
# path_peaks = 'Tests/out/Peaks'
# ---------------------------------------------------------------------------- #
Summarize_Data_For_Report = function(
  analysis_path  = NULL,
  analysis_names = NULL,
  report_chunks  = NULL,
  script_path    = NULL,
  outfile        = NULL
){
  # --------------------------------------------------------------- #
  # checks for default arugments
  deflist = as.list(formals(Summarize_Data_For_Report))
  arglist = as.list(match.call)
  arg.ind = names(deflist) %in% names(arglist)
  if(any(arg.ind))
      stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))

  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(genomation))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  source(file.path(script_path, 'Functions_Sample_Report.R'))

  # --------------------------------------------------------------- #
  .allowed_analysis_names = report_chunks
  set_analysis = setdiff(analysis_names, .allowed_analysis_names)
  if(length(set_analysis))
    stop('The current analysis are not supported:', set_analysis)

    rds_files = get_RDS_Input_files(analysis_path)

    # some sections of the pipeline are conditional. Because of that, the
    # document looks at which sections were executed and parses only this input
    # analysis_ind contains an index which says which chunks should be included
    # in the markdown
    analysis_ind = report_chunks %in% analysis_names
      names(analysis_ind) = report_chunks

    if(length(set_analysis) == 0)
        set_analysis = analysis_names

    # this loads in the data and processes the samples
    lstats = list()
    lstats = lapply(setNames(set_analysis, set_analysis),
                     function(x){
                     message(x)
                     match.fun(paste('format', x, sep='_'))(subset(rds_files, type==x))
                     })

    saveRDS(lstats, outfile)
}


# ---------------------------------------------------------------------------- #
Summarize_Data_For_Report(
  analysis_path  = argv$params[['analysis_path']],
  analysis_names = argv$params[['analysis_names']],
  report_chunks  = argv$params[['report_chunks']],
  script_path    = argv$params[['script_path']],
  outfile        = argv$output[['outfile']]

)
