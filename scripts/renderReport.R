library(rmarkdown)
library(jsonlite)

args <- commandArgs (trailingOnly=TRUE)

report <- args[1]
output <- args[2]
header <- args[3]
parameters <- fromJSON (args[4])

css <- tempfile (fileext="css")
cat ("body { padding-top: 60px }", file = css)

settings <- html_document(includes=includes(before_body=header),
                          css=css,
                          theme="lumen",
                          highlight="monochrome",
                          code_folding="hide",
                          depth="3",
                          toc=TRUE,
                          toc_float=TRUE,
                          number_sections=TRUE)
render(report,
       output_file=output,
       output_format=settings,
       params=parameters,
       intermediates_dir=tempfile(pattern="dir"))
