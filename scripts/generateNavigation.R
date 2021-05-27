library(data.table)
library(yaml)
library(rmarkdown)

args <- commandArgs (trailingOnly=TRUE)
reportsScriptDir <- args[1]
sampleSheetFile <- args[2]
navbarFile <- args[3]

sampleSheet <- fread (sampleSheetFile)
samples <- unique (sampleSheet$name)

isSampleReport <- Negate (function (file) {
    basename (file) %in% c("overview.Rmd", "index.Rmd")
})
sampleReports <- Filter (isSampleReport,
                         dir (reportsScriptDir,
                              pattern = ".Rmd$",
                              full.names = T))


navbar <- list("title" = "Reports",
               "left" = list(list("text" = "Overview",
                                  "href" = "overview.html"),
                             list("text" = "Timecourse",
                                  "href" = "index.html")))

# Read the YAML front matter of all reports.  Use the "nav" or "title"
# fields for the link text, or fall back to the basename of the report
# files.
tabs <- lapply(sampleReports, function (file) {
    meta <- yaml_front_matter (file)
    text <- Find (Negate (is.null),
                  c(meta$nav, meta$title, basename (file)))
    list("href" = basename(file),
         "text" = text,
         "menu" = lapply(samples, function (name) {
             list("text" = name,
                  "href" = paste0(name, ".", gsub(".Rmd", ".html", basename (file))))
         }))
})

navbar$left <- c(navbar$left, tabs)

tempfile <- navbar_html (navbar)
file.copy (tempfile, navbarFile)
unlink (tempfile)
