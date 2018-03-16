# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Peak_Statistics')


#
# peak_dict = yaml::yaml.load_file('Tests/sample_sheet.yaml')
# scriptdir = 'scripts'
# path_mapped = 'Tests/out/Mapped'
# path_peaks = 'Tests/out/Peaks'
# ---------------------------------------------------------------------------- #
Peak_Statistics = function(
  readnumber   = NULL,
  peak_dict    = NULL,
  lib_type_dict= NULL,
  outfile      = NULL,
  scriptdir    = NULL,
  path_mapped  = NULL,
  path_peak    = NULL,
  peaks_resize = NULL,
  ncores       = 8
){
  # --------------------------------------------------------------- #
  # checks for default arugments
  deflist = as.list(formals(Peak_Statistics))
  arglist = as.list(match.call)
  arg.ind = names(deflist) %in% names(arglist)
  if(any(arg.ind))
      stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))


  # -------------------------------------------------------------------------- #
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(genomation))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(GenomicAlignments))
  suppressPackageStartupMessages(library(Rsamtools))
  suppressPackageStartupMessages(library(BiocParallel))
  register(MulticoreParam(workers=ncores))
  source(file.path(scriptdir, 'Functions_Helper.R'), local=TRUE)

  # -------------------------------------------------------------------------- #
   message('Formatting sample sheet ...')
   print(path_peak)
   print(path_mapped)
    peaks_sheet = do.call(rbind, lapply(c('ChIP', 'Cont'), function(x){
        tibble(
            sample_name = names(peak_dict),
            bam_name    = lapply(peak_dict, '[[', x)) %>%
            dplyr::filter(!sapply(bam_name, is.null))           %>%
            unnest(bam_name)
        })) %>%
    merge(.list_BamFiles(path_mapped), by = 'bam_name') %>%
    merge(.list_qsortPeaks(path_peak), by.x = 'sample_name', by.y='chip_name')  %>%
    mutate(sample_id = paste(sample_name, bam_name, sep='_'))                   %>%
    mutate(bw_files  = str_replace(bam_file,'sorted.bam','bw'))                 %>%
    as.data.frame()
    peaks_sheet$library = unlist(lib_type_dict[peaks_sheet$bam_name])


  peaks_uniq = peaks_sheet %>%
    dplyr::select(sample_name, bed_file, library) %>%
    distinct()

  # -------------------------------------------------------------------------- #
  mapped_reads = readRDS(readnumber)      %>%
    dplyr::filter(stat == 'mapped.total') %>%
    mutate(stat = NULL)                   %>%
    dplyr::rename(mapped_total = value)   %>%
    dplyr::rename(bam_name = sample_name)

  # -------------------------------------------------------------------------- #
  message('Reading peaks ...')
    peaks     = lapply(peaks_uniq$bed_file,
      function(x)try(readGeneric(x), silent=TRUE))
    peaks_ind = sapply(peaks, function(x)class(x)=='GRanges')
    peaks = GRangesList(peaks[peaks_ind])
    names(peaks) = peaks_uniq$sample_name[peaks_ind]

    # -------------------------------------------------------------------------- #
  
  message('Counting reads ...')
    cnts = lapply(c('single', 'paired'), function(x){
        
        sub = subset(peaks_sheet, library==x)
        if(nrow(sub) > 0){
            bam_files  = BamFileList(unique(sub$bam_file),  yieldSize = 100000)
            summarizeOverlaps(peaks, bam_files, singleEnd = x == 'single', ignore.strand=TRUE)
        }else{
            NULL
        }
      })
    cnts   = cnts[sapply(cnts, function(x)!is.null(x))] 
    cntmat = as.data.frame(do.call(cbind, lapply(cnts, function(x)assays(x)[[1]]))) %>%
      mutate(sample_name = names(peaks)) %>%
      merge(x=peaks_sheet) %>%
      mutate(peak_number = elementNROWS(peaks)[sample_name]) %>%
      merge(mapped_reads, by='bam_name')

  # ------------------------------------------------------------------------- #
  message('Views ...')
    peak_names = names(peaks)
    lms = lapply(setNames(peak_names, peak_names), function(x){
      sml = ScoreMatrixList(
        windows = resize(peaks[[x]], width=peaks_resize, fix='center'),
        target  = unique(subset(peaks_sheet, sample_name == x)$bw_files),
        type    = 'bigWig',
        bin.num = 20,
        cores   = ncores)
      })

  # ------------------------------------------------------------------------- #
  message('Profiles ...')
    mms = ListListToDataFrame(lms, colMeans, 'mean')
    qms = ListListToDataFrame(lms, matrixStats::colSds, 'IQR')
    mes = ListListToDataFrame(lms, matrixStats::colMedians, 'median')
    dms = MergeStats(mms, qms)
    dms = MergeStats(dms, mes)

    lout = list(
      peaks_sample      = cntmat,
      score_matrix_list = lms,
      peak_profiles     = dms
      )
    saveRDS(lout, outfile)
}

# ---------------------------------------------------------------------------- #
ListListToDataFrame = function(lms, fun, colname){
  ll = lapply(lms, lapply, fun)
  d = do.call(rbind, lapply(names(ll), function(x){
          z = ll[[x]]
          data.frame(
            peak     = x,
            bw_name  = rep(names(z), sapply(z,length)),
            position = unlist(lapply(z, function(y)seq(length(y)))),
            score    = unlist(z)
            )})
        )

    colnames(d)[4] = colname
    return(d)
}

# ---------------------------------------------------------------------------- #
MergeStats = function(d1, d2){

  d1 = d1 %>%
    mutate(key = paste(peak, bw_name, position))
  d2 = d2 %>%
      mutate(key = paste(peak, bw_name, position)) %>%
      dplyr::select(-(peak:position))
  d = merge(d1, d2, by='key') %>%
    arrange(bw_name, position) %>%
    dplyr::select(-key)
  return(d)
}

# ---------------------------------------------------------------------------- #
Peak_Statistics(
  readnumber   = argv$input[['readnumber']],
  outfile      = argv$output[['outfile']],
  peak_dict    = argv$params[['peak_dict']],
  lib_type_dict= argv$params[['lib_type_dict']],
  scriptdir    = argv$params[['scriptdir']],
  path_mapped  = argv$params[['path_mapped']],
  path_peak    = argv$params[['path_peak']],
  peaks_resize = argv$params[['peaks_resize']],
  ncores       = argv$params[['threads']]
  
)
