# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This file incorporates code by Yihui Xie <xie@yihui.name> under the
# terms of the GPL version 3.  The "Rscript_render2" procedure was
# adapted from
# https://github.com/rstudio/bookdown/blob/master/R/utils.R

## modified from https://github.com/rstudio/bookdown/blob/master/R/utils.R#L205-L208
Rscript_render2 = function(file, ..., log.file=NULL) {
  args = shQuote(c(bookdown:::bookdown_file('scripts', 'render_one.R'), file, ...))
  ## we append the stderr and stdout to the log file
  if(!is.null(log.file)) args = paste(args,">>",log.file,"2>&1")
  if (bookdown:::Rscript(args)!= 0) stop('Failed to compile ', file)
}

## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included
render2Markdown <- function(reportFile,
                            outFile,
                            outDir,
                            finalReportDir,
                            report.params=NULL,
                            self.contained=FALSE,
                            logFile = NULL)
{
  
  output_format = rmarkdown::all_output_formats(reportFile, 'UTF-8')
  
  if(!dir.exists(finalReportDir)) dir.create(finalReportDir, recursive = TRUE)
  
  if(is.null(report.params)) report.params <- list()

  ## render single report
  message("render single report")

  ## make independent intermediate dirs
  interDir <- paste0(outDir,"/",outFile,"_tmp")

  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    intermediates_dir = interDir,
    output_file = outFile,
    knit_root_dir = outDir,
    output_format = rmarkdown::html_document(
      code_folding = "hide",
      code_download = TRUE,
      self_contained = self.contained
    ),
    params=c(report.params,
             list("sessioninfo"=TRUE,
                  "references"=TRUE)),
    quiet = FALSE,
    clean = TRUE
  )

  on.exit(unlink(interDir,recursive = TRUE),add = TRUE)

  ## render for multireport
  message("render for multireport")
  
  ## save a copy of render arguments in a temp file
  render_args = tempfile('render', tmpdir = finalReportDir, '.rds')
  on.exit(unlink(render_args), add = TRUE)
  saveRDS(
    list(output_format = "rmarkdown::html_notebook",
         params=c(report.params,
                  list("sessioninfo"=FALSE,
                       "references"=FALSE)),
         output_file = outFile,
         output_dir = finalReportDir,
         intermediates_dir = finalReportDir, 
         knit_root_dir = finalReportDir,
         clean = FALSE, 
         quiet = FALSE,
         envir = parent.frame()
    ),
    render_args
  )
  
  ## store metadata for Final Report
  render_meta = paste0(finalReportDir,"/knitr_meta.rds")
  
  Rscript_render2(file = reportFile,render_args,render_meta,log.file = logFile)
  ## move the sessioninfo to final folder
  session_file <- list.files(path = outDir,pattern = "session",full.names = TRUE)
  session_file <- session_file[endsWith(session_file,".rds")]
  if(length(session_file)!=0) {
    file.rename(from = session_file,
                to = paste0(finalReportDir,"/",basename(session_file)))
  }
}

# Function Definitions ----------------------------------------------------

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
  } else if( length(sessions)==1) {
    return(sessions[[1]])
  } else {
    stop("Not enough sessions given")
  }
}

## http://github.com/rstudio/bookdown/blob/master/R/utils.R#L141-L156
merge_chapters2 = function(files, to, before = NULL, after = NULL, orig = files) {
  ## in the preview mode, only use some placeholder text instead of the full Rmd
  # preview = opts$get('preview'); input = opts$get('input_rmd')
  content = unlist(mapply(files, orig, SIMPLIFY = FALSE, FUN = function(f, o) {
    x = bookdown:::readUTF8(f)
    x = # if (preview && !(o %in% input)) create_placeholder(x) else {
      bookdown:::insert_code_chunk(x, before, after)
    #}
    c(x, '', paste0("********"),'',paste0('<!--chapter:end:', o, '-->'), '')
  }))
  #if (preview && !(files[1] %in% input))
  #content = c(fetch_yaml(readUTF8(files[1])), content)
  unlink(to)
  bookdown:::writeUTF8(content, to)
  Sys.chmod(to, '644')
}



#' Combine multiple markdowns to one
#'
#' @param final_output string name of the final single file (creates *.md and *.html)
#' @param finalreportdir string output directory 
#' @param index path to index file prepended to final markdown  
#' @param references path to reference file appended to final markdown
#' @param sessioninfo path to session info file appended to final markdown
#' @param workdir path to working 
#' @param clean 
#'
#' @return
#' @export
#'
#' @examples
render2multireport <- function (final_output,
                               finalreportdir,
                               index=NULL,
                               references=NULL,
                               sessioninfo=NULL,
                               self_contained=TRUE,
                               clean=TRUE) {
  

  file.copy(paste0(dirname(index), "/", "pigx_bsseq_logo.html"),
            paste0(finalreportdir, "/", "pigx_bsseq_logo.html"))

  meta.file <- list.files(path = finalreportdir,
                          pattern = "knitr_meta.rds",
                          full.names = TRUE)
  
  if(!is.null(index) || !is.null(references) || !is.null(sessioninfo)){
    
    ## save a copy of render arguments in a temp file
    render_args = tempfile('render', tmpdir = finalreportdir, '.rds')
    on.exit(unlink(render_args), add = TRUE)
    saveRDS(
      list(output_format = "rmarkdown::html_fragment",
           output_dir = finalreportdir,
           intermediates_dir = finalreportdir, 
           quiet = TRUE,
           clean = FALSE),
      render_args
    )
  }
  
  if(!is.null(index)) {
    
    render_list = readRDS(render_args)
    #render_list$params=list(wdir = workdir)
    saveRDS(
      render_list,
      render_args
    )
  
    index_ = paste0(finalreportdir, "/", basename(index))
    file.copy(index, index_)
    bookdown:::Rscript_render(file = normalizePath(index_),render_args,meta.file)
  
  }
  
  if(!is.null(sessioninfo)){
    
      
    sessionfiles <- list.files(path = finalreportdir, pattern = "session", full.names = TRUE)
    sessionfiles <- sessionfiles[endsWith(sessionfiles,".rds")]
    
    finalSessionInfoFile <- paste0(finalreportdir,"/","finalsessioninfo.rds")
    if( length(sessionfiles) >= 1 ) {
      finalSessionInfo <- mergeSessionInfos(lapply(sessionfiles,readRDS))
      saveRDS(finalSessionInfo,file = finalSessionInfoFile)
      on.exit(unlink(sessionfiles),add = TRUE)
    } else stop("no session info file found")


      
    render_list = readRDS(render_args)
    # render_list$knitr_root_dir = finalreportdir
    render_list$params=list("sessioninfo" = finalSessionInfoFile)
    saveRDS(
      render_list,
      render_args
    )
    
    sessioninfo_ = paste0(finalreportdir, "/", basename(sessioninfo))
    file.copy(sessioninfo, sessioninfo_)
    bookdown:::Rscript_render(file = normalizePath(sessioninfo_),render_args,meta.file)
    
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
    
    references_ = paste0(finalreportdir, "/", basename(references))
    file.copy(references, references_)
    bookdown:::Rscript_render(file = normalizePath(references_),render_args,meta.file)
    
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
  
  merge_chapters2(unlist(lapply(meta ,function(x) attr(x,"intermediates")[2])),
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
  

  
  file.link(from = paste0(finalreportdir,"/",bookdown:::with_ext(basename(final_output),".html")),
              to = final_output)
}

