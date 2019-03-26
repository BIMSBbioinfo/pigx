# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Knit_Report')


# ---------------------------------------------------------------------------- #
#' Knit_Report
#'
#' @param report_template rmd file template
#' @param analysis_path 
#' @param analysis_names 
#' @param report_chunks a vector which defines which parts of the report will be
#' knit
#' @param script_path locations of function files
#' @param width_params a list of parameters which define the width of annotation regions
#' @param infile snakemake input files - defines the dependency tree for the report
#' these are basically all files that are created in the pipeline
#' @param outfile name of the output report
#' @param logo 
#' @param selfContained an indicator on whether to create a self contained report
#' @param quiet knitr parameter
#'
#' @return
#' @export
#'
#' @examples
Knit_Report = function(
    report_template = NULL,
    analysis_path   = NULL,
    analysis_names  = NULL,
    analysis_mode   = NULL,
    genome          = NULL,
    report_chunks   = NULL,
    script_path     = NULL,
    width_params    = NULL,
    infile,
    outfile,
    logo,
    selfContained = TRUE,
    quiet         = FALSE
){

    suppressPackageStartupMessages(library(htmlwidgets))
    suppressPackageStartupMessages(library(rmarkdown))
    workdir = dirname(outfile)
    tempdir = file.path(workdir,'Temp')
    htmlwidgets::setWidgetIdSeed(1234)
    rmarkdown::render(
        input             = report_template,
        output_dir        = workdir,
        intermediates_dir = tempdir,
        clean             = TRUE,
        output_file       = outfile,
        output_format     = rmarkdown::html_document(
            code_folding    = 'hide',
            depth           = 2,
            toc             = TRUE,
            toc_float       = TRUE,
            theme           = 'lumen',
            number_sections = TRUE
        ),
        output_options = list(self_contained = selfContained),
        params = list(analysis_path  = analysis_path,
                      analysis_names = analysis_names,
                      report_chunks  = report_chunks,
                      analysis_mode   = analysis_mode,
                      genome          = genome,
                      script_path    = script_path,
                      infile         = infile,
                      logo           = logo),
        quiet = quiet
    )

    if(dir.exists(tempdir)) {
        unlink(tempdir, recursive = TRUE)
    }
}


# ---------------------------------------------------------------------------- #
Knit_Report(
    report_template = argv$params[['report_template']],
    analysis_path   = argv$params[['analysis_path']],
    analysis_names  = argv$params[['analysis_names']],
    analysis_mode   = argv$params[['analysis_mode']],
    genome          = argv$params[['genome']],
    report_chunks   = argv$params[['report_chunks']],
    script_path     = argv$params[['script_path']],
    width_params    = argv$params[['width_params']],
    logo            = argv$params[['logo']],
    infile          = argv$input[['infile']],
    outfile         = argv$output[['outfile']]
)
