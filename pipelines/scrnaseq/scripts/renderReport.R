# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('renderReport.R')



# ---------------------------------------------------------------------------- #
#2. Define report rendering function
#' runReport
#'
#' Generate a PiGx-scRNAseq Report in a self-contained HTML file
#'
#' @param report_file Path to .Rmd script to generate a HTML report
#' @param sceRds_file Path to the RDS format file containing the
#'   SingleCellExperiment object
#' @param outfile
#' @param workdir Path to working directory where the output files will be
#'   written
#' @param prefix Prefix to be attached to the beginning of output files
#' @return An html generated using rmarkdown/knitr/pandoc that contains
#'   interactive figures, tables, and text that provide an overview of the
#'   experiment
Render_Report = function(
    sceRds_file  = NULL,
    report_file  = NULL,
    output_file  = NULL,
    read_stats   = NULL,
    workdir      = NULL,
    prefix       = 'PiGx',
    path_mapped  = NULL

){
    # ------------------------------------------------------------------------- #
    suppressPackageStartupMessages({
        library(rmarkdown)
    })

    if(is.null(report_file))
        stop("Missing argument: reportFile. Provide the path to .Rmd file")

    if(is.null(sceRds_file))
        stop("Missing argument: sceRdsFile. Provide the path to .Rds file containing the SingleCellExperiment object.")

    if(is.null(output_file))
        stop('Outfile is not defined')

    if(is.null(read_stats))
        stop("Missing argument: read_stats Provide the pathss to read statistic files")

    if(is.null(workdir))
        workdir = getwd()

    # ------------------------------------------------------------------------- #


    htmlwidgets::setWidgetIdSeed(1234)
    rmarkdown::render(
        input             = report_file,
        output_dir        = workdir,
        intermediates_dir = file.path(workdir, prefix),
        clean             = TRUE,
        output_file       = output_file,
        output_format = rmarkdown::html_document(
            code_folding    = 'hide',
            depth           = 2,
            toc             = TRUE,
            toc_float       = TRUE,
            theme           = 'lumen',
            number_sections = TRUE
        ),
        output_options = list(self_contained = selfContained),
        params = list(sceRds_file = sceRds_file,
                      read_stats  = read_stats,
                      output_file = output_file,
                      workdir     = workdir,
                      path_mapped = path_mapped),
        quiet = FALSE
    )

    if(dir.exists(file.path(workdir, prefix))) {
        unlink(file.path(workdir, prefix), recursive = TRUE)
    }
}


# ---------------------------------------------------------------------------- #
Render_Report(
    sceRds_file = argv$input[['sceRds_file']],
    read_stats  = argv$input[['read_stats']],
    output_file = argv$output[['outfile']],
    report_file = argv$params[['report_file']],
    workdir     = argv$params[['workdir']],
    path_mapped = argv$params[['path_mapped']]
)
