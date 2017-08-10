## Wrapper function to combine a set of precompiled Rmd scripts
## and render them as multiRmd report 


#' Merge two sessionInfo objects
#' 
#' This function takes two sessionInfo objects and merges them into one.
#'
#' @param sessionX 
#' @param sessionY 
#'
#' @return a sessionInfo object
#'
#'
#' @examples
.mergeSessionInfo <- function(sessionX,sessionY) {
  
  ## check which of both has more entries (minimum 8, maximum 10)
  if( length(sessionX)>=length(sessionY ) ){
    z <- sessionX
  } else { z <- sessionY }
  
  ## iterate over entries
  for( i in names(z)) { 
    
    ## merge content if entry is shared
    if( (i %in% names(sessionX)) & (i %in% names(sessionY)) ){
      z[[i]] <- union(sessionX[[i]],sessionY[[i]])
      names(z[[i]]) <- union(names(sessionX[[i]]),names(sessionY[[i]]))
    }
  }
  
  ## remove attached packages from namespace
  if(all(c("loadedOnly","otherPkgs") %in% names(z)) ) { 
    z[["loadedOnly"]] <- z[["loadedOnly"]][! z$loadedOnly %in% z$otherPkgs] }
  
  class(z) <- "sessionInfo"
  
  return(z)
}


#' Merge multiple sessionInfo objects
#' 
#' This function takes a list of sessionInfo objects and merges them into one
#'
#' @param sessionX 
#' @param sessionY 
#'
#' @return
#' @export
#'
#' @examples
mergeSessionInfos <- function(sessions = list()) {
  
  if(! all(sapply(sessions,class)=="sessionInfo") ) stop("Not all objects of class sessionInfo")
  
  if(length(sessions)>2) {
    Reduce(.mergeSessionInfo,sessions)
  } else if( length(sessions)==2) {
    .mergeSessionInfo(sessionX = sessions[[1]], 
                      sessionY = sessions[[2]])
  } else {
    stop("Not enough sessions given")
  }
}


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
                               sessioninfo=NULL,
                               self_contained=TRUE,
                               clean=FALSE) {
  
  
  meta.file <- list.files(path = finalreportdir,
                          pattern = "knitr_meta.rds",
                          full.names = TRUE)
  
  if(!is.null(index) || !is.null(references) || !is.null(sessioninfo)){
    
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
  
  if(!is.null(sessioninfo)){
    
    finalSessionInfoFile <- list.files(path = finalreportdir,
               pattern = "finalsessioninfo.rds",
               full.names = TRUE)
    if(!file.exists(finalSessionInfoFile)) {
      
      sessionfiles <- list.files(path = finalreportdir, pattern = "session", full.names = TRUE)
      sessionfiles <- sessionfiles[endsWith(sessionfiles,".rds")]
      finalSessionInfo <- mergeSessionInfos(lapply(sessionfiles,readRDS))
      saveRDS(finalSessionInfo,file = paste0(finalreportdir,"/","finalsessioninfo.rds"))
      
      on.exit(unlink(sessionfiles),add = TRUE)
      
    }
    
    render_list = readRDS(render_args)
    # render_list$knitr_root_dir = finalreportdir
    render_list$params=list("sessioninfo" = finalSessionInfoFile)
    saveRDS(
      render_list,
      render_args
    )
    
    bookdown:::Rscript_render(file = normalizePath(sessioninfo),render_args,meta.file)
    
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
  ## change ordering of static files
  indexPos = match('index',bookdown:::with_ext(basename(unlist(meta)),""))
  if (!is.na(indexPos)) meta <- c(meta[indexPos], meta[-indexPos])
  sessioninfoPos = match('sessioninfo',bookdown:::with_ext(basename(unlist(meta)),""))
  if (!is.na(sessioninfoPos)) meta <- c(meta[-sessioninfoPos], meta[sessioninfoPos])
  referencesPos = match('references',bookdown:::with_ext(basename(unlist(meta)),""))
  if (!is.na(referencesPos)) meta <- c(meta[-referencesPos], meta[referencesPos])
  
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
      self_contained = self_contained,
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
                   references = snakemake@input[["references"]],
                   sessioninfo = snakemake@input[["sessioninfo"]])

# finalReportDir = "Final_Report/"
# 
# render2Markdown(reportFile = normalizePath(snakemake@input[["template"]]),
#                 outFile = basename(snakemake@output[["report"]]),
#                 outDir = normalizePath(dirname(snakemake@output[["report"]])),
#                 finalReportDir = finalReportDir,
#                 report.params = snakemake@params[nchar(names(snakemake@params)) > 0] )
