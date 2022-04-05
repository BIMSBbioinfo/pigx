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
    report_params = NULL,
    outfile,
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
        output_file       = basename(outfile),
        output_format     = rmarkdown::html_document(
            code_folding    = 'hide',
            depth           = 2,
            toc             = TRUE,
            toc_float       = TRUE,
            theme           = 'lumen',
            number_sections = TRUE
        ),
        output_options = list(self_contained = selfContained),
        params = report_params,
        quiet = quiet
    )

    if(dir.exists(tempdir)) {
        unlink(tempdir, recursive = TRUE)
    }
}


# ---------------------------------------------------------------------------- #
## check wether all required report params are given 
defined_params <- knitr::knit_params(readLines(argv$params[['report_template']]),
                                 evaluate = TRUE)

given_params <- names(defined_params) %in% names(c(argv[["params"]],argv[["input"]]))
if( !all(given_params ))  {
      warning("Missing values for parameters: ",
                        paste(names(defined_params)[!given_params],collapse = ", "))
}
report_params <- c(argv[["params"]],argv[["input"]])[names(c(argv[["params"]],argv[["input"]])) %in% names(defined_params)]

print(report_params)

# ---------------------------------------------------------------------------- #
Knit_Report(
    report_template = argv$params[['report_template']],
    report_params   = report_params,
    outfile         = argv$output[['outfile']]
)

warnings()
