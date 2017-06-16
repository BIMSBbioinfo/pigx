## Wrapper function to combine a set of precompiled Rmd scripts
## and render them as multiRmd report 




merge_chapters2 = function(files, to, before = NULL, after = NULL, orig = files) {
  ## in the preview mode, only use some placeholder text instead of the full Rmd
  # preview = opts$get('preview'); input = opts$get('input_rmd')
  content = unlist(mapply(files, orig, SIMPLIFY = FALSE, FUN = function(f, o) {
    x = bookdown:::readUTF8(f)
    x = # if (preview && !(o %in% input)) create_placeholder(x) else {
      bookdown:::insert_code_chunk(x, before, after)
    #}
    c(x, '', paste0('<!--chapter:end:', o, '-->'), '')
  }))
  #if (preview && !(files[1] %in% input))
  #content = c(fetch_yaml(readUTF8(files[1])), content)
  unlink(to)
  bookdown:::writeUTF8(content, to)
  Sys.chmod(to, '644')
}



render2multireport <- function(final_output,
                               finalreportdir,
                               index=NULL,
                               references=NULL,
                               workdir=NULL,
                               clean=FALSE) {
  
  
  meta.file <- list.files(path = finalreportdir,
                          pattern = "knitr_meta.rds",
                          full.names = TRUE)
  
  if(!is.null(index) || !is.null(references)){
    
    ## save a copy of render arguments in a temp file
    render_args = tempfile('render', tmpdir = finalreportdir, '.rds')
    on.exit(unlink(render_args), add = TRUE)
    saveRDS(
      list(output_format = "rmarkdown::html_document",
           intermediates_dir = finalreportdir, 
           quiet = TRUE,
           clean = FALSE),
      render_args
    )
  }
  
  if(!is.null(index)) {
    
    render_list = readRDS(render_args)
    render_list$params=list(wdir = workdir)
    saveRDS(
      render_list,
      render_args
    )
  
    bookdown:::Rscript_render(file = normalizePath(index),render_args,meta.file)
  
  }
  if(!is.null(references)){
    
    ## comment out all other reference sections
    meta <- readRDS(meta.file)
    lapply(unlist(meta),
           FUN = function(file){
             
             x <- bookdown:::readUTF8(file)
             ## check for reference section
             refInd <- grep("References",x)
             if(!length(refInd)==0 && grepl("^#", x[refInd]))
               x[refInd] = paste0('<!--keep:final:', x[refInd], '-->')
             
             bookdown:::writeUTF8(x,file)
           })
    
    bookdown:::Rscript_render(file = normalizePath(references),render_args,meta.file)
    
  }
  
  meta <- readRDS(meta.file)
  inputFiles <- unlist(meta)
  indexPos = match('index',bookdown:::with_ext(basename(unlist(meta)),""))
  if (!is.na(indexPos)) meta <- c(meta[indexPos], meta[-indexPos])
  
  ## this part combines all meta output from all files and
  ## sorts by input files
  # meta <- bookdown:::clean_meta(meta_file = list.files(path = finalreportdir,
  #                                                      pattern = "knitr_meta.rds",
  #                                                      full.names = TRUE),
  #                               files = inputFiles)
  
  
  
  # merge_chapters2(c(index,unlist(meta),references), 
  #                 to = paste0(finalreportdir,"/","finalreport",".Rmd"),  
  #                 orig = c(index,template.list,references))
  
  
  merge_chapters2(unlist(meta),
                  to = paste0(finalreportdir,"/",bookdown:::with_ext(basename(final_output),".md"))#,
                  #orig = template.list
  )
  
  
  
  
  knit_meta = unlist(lapply(meta, attr, 'knit_meta', exact = TRUE), recursive = FALSE)
  intermediates = unlist(lapply(meta, attr, 'intermediates', exact = TRUE))
  if (clean) on.exit(unlink(intermediates, recursive = TRUE), add = TRUE)
  
  rmarkdown::render(
    input = paste0(finalreportdir,"/",bookdown:::with_ext(basename(final_output),".md")),
    output_format = rmarkdown::html_notebook(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = FALSE,
      code_folding = "hide",
      self_contained = TRUE,
      includes = list(in_header = "pigx_bsseq_logo.html"),
      bibliography= "reports.bib"
    ),
    output_file = paste0(finalreportdir,"/",bookdown:::with_ext(basename(final_output),".html")),
    # output_file = bookdown:::with_ext(basename(final_output),".html"),
    knit_root_dir = finalreportdir,
    output_dir = finalreportdir,
    knit_meta = knit_meta,
    clean = clean,
    run_pandoc = TRUE)
  

  
  file.copy(from = paste0(finalreportdir,"/",bookdown:::with_ext(basename(final_output),".html")),
              to = final_output)
  
  if(clean) unlink(finalreportdir,recursive = TRUE)
  
}


## catch output and messages into log file
out <- file(snakemake@log[[1]], open = "wt")
sink(out,type = "output")
sink(out, type = "message")

cat(paste(  
  Sys.time(),"\n\n",
  "Rendering final report: ",basename(snakemake@output[["finalreport"]]),"\n",
  "into directory:",normalizePath(dirname(snakemake@output[["finalreport"]])),"\n",
  "\n"
))



render2multireport(final_output = normalizePath(snakemake@output[["finalreport"]]),
                   finalreportdir = normalizePath(snakemake@params[["finalreportdir"]]),
                   workdir = normalizePath(snakemake@params[["workdir"]]),
                   index = snakemake@input[["index"]],
                   references = snakemake@input[["references"]])

# finalReportDir = "Final_Report/"
# 
# render2Markdown(reportFile = normalizePath(snakemake@input[["template"]]),
#                 outFile = basename(snakemake@output[["report"]]),
#                 outDir = normalizePath(dirname(snakemake@output[["report"]])),
#                 finalReportDir = finalReportDir,
#                 report.params = snakemake@params[nchar(names(snakemake@params)) > 0] )
