#prepare _site.yml and other Rmd files to be rendered into a html report (see renderSite rule)

#        "{RSCRIPT} {params.script} {params.report_scripts_dir} {SAMPLE_SHEET_FILE} {CUT_SITES_FILE} {OUTPUT_DIR} {REPORT_DIR} {RSCRIPT} > {log} 2>&1"
library(data.table)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

reportsScriptDir <- args[1] #folder that contains the Rmd reports that will be rendered into a site
sampleSheetFile <- args[2] #path to sample_sheet.csv file
cutSitesFile <- args[3] #path to cut_sites.bed file
comparisonsFile <- args[4] #path to tsv file containing list of pairs of samples to compare 
pipelineOutputDir <- args[5] #root folder where the pipeline is written to
siteDir <- args[6] #path to folder where the site will be generated

#read sample sheet
sampleSheet <- data.table::fread(sampleSheetFile)

#1. for each target region in sample sheet, create a copy of all Rmd files
#replacing the target_name parameter in the Rmd file's yaml header.
rmd_files <- dir(reportsScriptDir, pattern = '.Rmd$', full.names = T)

for(f in rmd_files) {
  if(basename(f) == 'index.Rmd') {
    file.copy(from = f, to = siteDir)
  } else if (grepl('^_.+.Rmd$', basename(f))){ #copy child Rmd files don't duplicate 
    file.copy(from = f, to = siteDir)
  } else {
    text <- readLines(f)
    for(target_name in unique(sampleSheet$target_name)) {
      outFile <- file.path(siteDir, paste0(target_name, ".", basename(f)))
      modified_text <- gsub(x = text,
                            pattern = ' @PARAMS_TARGET_NAME@:',
                            replacement = paste0("target_name: '", target_name, "'"))
      writeLines(text = modified_text, con = outFile)
    }
  }
}

# 2. Create a config.yml file that is a global configuration file that is
# re-used by separate Rmd files in the target site
config_yml <- list('sample_sheet' = sampleSheetFile,
                   'pigx_output_dir' = pipelineOutputDir,
                   'cut_sites_file' = cutSitesFile, 
                   'comparisons_file' = comparisonsFile)
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
                     list("text" = "About", "href" = "index.html")))

#for all rmd files except index.Rmd and child Rmd files (beginning with "_"), 
# create a navbar item (for each target_name in sample sheet)
tabs <- grep("(^index.Rmd$)|(^_.+.Rmd$)", basename(rmd_files), invert = T, value = T)

navbar_yml$left <- c(navbar_yml$left, lapply(tabs,
       function(f) {
         L <- list("text" =  gsub(".Rmd", "", f),
                   "menu" = lapply(unique(sampleSheet$target_name), function(t) {
                     href <- paste0(t, ".", gsub(".Rmd", ".html", f))
                     list("text" = t, "href" = href)
         }))
         return(L)
}))

#update site_yml with navbar
site_yml$navbar <- navbar_yml

yaml::write_yaml(site_yml, file = file.path(siteDir, "_site.yml"))
