library(rmarkdown)
library(jsonlite)

args <- commandArgs (trailingOnly=TRUE)

report <- args[1]
output <- args[2]
header <- args[3]
parameters <- fromJSON (args[4])

css <- tempfile (fileext="css")
cat ("\
#logo { position: relative; } \
#logo img { \
  width: 125px; \
} \
div#TOC { margin-top: 100px; } \
div#TOC span.header-section-number { display: none } \
span.header-section-number { color: #93a1a1; } \
@media all and (min-width:768px){ \
  body { font-size: 20px; padding-top: 60px; position: absolute; } \
  div#TOC { font-size: smaller; margin-top: 175px; } \
  #logo img { \
    top: 75px; \
    left: 50px; \
    position: fixed; \
  } \
  span.header-section-number { \
    position: absolute; \
    text-align: right; \
    width: 2em; \
    margin-left: -2.5em } \
} \
p { line-height: 1.5 } \
h1 { margin-top: 2em } \
h1.title { margin-top: 20px } \
h2 { margin-top: 2em } \
h4.date { color: #93a1a1 } \
h4.date:before { content: \"Created: \" } \
sample {font-weight: bold; color: #2C1966} \
sample:before {content: \" Sample: \"} \
", file = css)

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
