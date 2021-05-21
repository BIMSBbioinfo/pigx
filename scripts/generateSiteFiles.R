#Title: generate SiteFiles
# Based on the generateSiteFiles.R from Bora Uyar : https://github.com/BIMSBbioinfo/crispr_DART/blob/master/src/generateSiteFiles.R
#not changed yet
# the target regions would be samples in our case



library(data.table)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

reportsScriptDir <- args[1] #folder that contains the Rmd reports that will be rendered into a site
sampleSheetFile <- args[2] #path to sample_sheet.csv file
krakenDir <- args[3]
coverage_dir <- args[4]
variants_dir <- args[5]
sigmut_db <- args[6]
siteDir <- args[7] #path to folder where the site will be generated aka report_dir ?
# necessary? - also see below l67
pipelineOutputDir <- "../"


#read sample sheet
sampleSheet <- data.table::fread(sampleSheetFile)

#1. for each sample in sample sheet, create a copy of all Rmd files
#replacing the sample_name in the Rmd file's yaml header.
rmd_files <- dir(reportsScriptDir, pattern = '.Rmd$', full.names = T)

for(f in rmd_files) {
  if(basename(f) == 'timecourse.Rmd') {
    file.copy(from = f, to = siteDir)
  } else if (basename(f) == 'index.Rmd'){
    file.copy(from = f, to = siteDir)
  } else if (grepl('^_.+.Rmd$', basename(f))){ #copy child Rmd files don't duplicate 
    file.copy(from = f, to = siteDir)
  } else {
    text <- readLines(f)
    for(sample_name in unique(sampleSheet$name)) {
      outFile <- file.path(siteDir, paste0(sample_name, ".", basename(f)))
      modified_text <- gsub(x = text,
                            pattern = ' @PARAMS_SAMPLE_NAME@',
                            replacement = paste0("sample_name: '", sample_name, "'"))
      writeLines(text = modified_text, con = outFile)
    }
  }
}

# 2. Create a config.yml file that is a global configuration file that is
# re-used by separate Rmd files in the target site
config_yml <- list('sample_sheet' = sampleSheetFile,
                   'kraken_dir' = krakenDir,
                   'coverage_dir'=coverage_dir,
                   'variants_dir'=variants_dir,
                   'sigmut_db'=sigmut_db,
                   'pipeline_output_dir' = pipelineOutputDir,
                   'site_dir' = siteDir)
yaml::write_yaml(config_yml, file = file.path(siteDir, "config.yml"))

# 3. Create a _site.yml file that determines the layout of the rendered html
# site the layout is determined by the existing samples and their target regions
# in the sample sheet
# Each tab in the main frame will be represented by a drop-down menu of one type
# of analysis, where the menu items correspond to different target regions

# a yaml file is basically is list of list of ... lists

#global options for the output
site_yml <- list("name" = "Results",
                 "output_dir" = ".",
                 "output" = list("html_document" = list("theme" = "lumen",
                                                        "highlight" = "monochrome",
                                                        "code_folding" = "hide",
                                                        "depth" = "3",
                                                        "toc" = TRUE,
                                                        "toc_float" = TRUE,
                                                        "number_sections" = TRUE)))

#navigation bar options for the output
#NOTE:navigation bar is the whole reason we use this script, we want this to be dynamically generated
#depending on the available samples in the pipeline's output
navbar_yml <- list('title' = 'Reports',
                   'left' = list(
                     list("text" = "timecourse", "href" = "timecourse.html")))

#for all rmd files except index.Rmd(timecourse) and child Rmd files (beginning with "_"), 
# create a navbar item (for each sample in sample sheet)
tabs <- grep("(^timecourse.Rmd$)|(^_.+.Rmd$)", basename(rmd_files), invert = T, value = T)
tabs <- grep("(^index.Rmd$)|(^_.+.Rmd$)", basename(tabs), invert = T, value = T)

navbar_yml$left <- c(navbar_yml$left, lapply(tabs,
                                             function(f) {
                                               L <- list("text" =  gsub(".Rmd", "", f),
                                                         "menu" = lapply(unique(sampleSheet$name), function(t) {
                                                           href <- paste0(t, ".", gsub(".Rmd", ".html", f))
                                                           list("text" = t, "href" = href)
                                                         }))
                                               return(L)
                                             }))

#update site_yml with navbar
site_yml$navbar <- navbar_yml

yaml::write_yaml(site_yml, file = file.path(siteDir, "_site.yml"))

