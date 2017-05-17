


args = commandArgs(trailingOnly=TRUE)

reportFile = "./methCall.report.Rmd"
outFile=args[2]
outDir=args[5]

params = list(input=args[1],
     output=args[2],
     context = args[3],
     sampleid = args[4],
     pathout= args[5],
     assembly = args[5],
     mincov=args[6],
     minqual=args[7],
     reportFile = args[8])






render2Report <- function(reportFile,
                          reportLogo,
                          outFile,
                          outDir,
                          ...)
{
  
  new.params = list(...)
  
  rmarkdown::render(
    input = reportFile, 
    output_dir = outDir,
    intermediates_dir = paste0(outDir,"/tmp"),
    output_file = outFile,
    output_format = rmarkdown::html_notebook(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = FALSE,
      code_folding = "hide",
      self_contained = TRUE
      #includes = list(in_header = "/Users/agosdsc/Development/Snakemake/pigx_bsseq/pigx_bsseq_logo.html"
      #                )
    ),
    params = new.params,
    quiet = TRUE,
    clean = TRUE
  )
  unlink(paste0(outDir,"/Tmp"),recursive = TRUE)
}



render2Report(reportFile = reportFile,
              outFile = outFile,
              outDir = outDir,
              params
)

