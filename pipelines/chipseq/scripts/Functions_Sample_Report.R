# ---------------------------------------------------------------------------- #
get_RDS_Input_files = function(path){

    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tibble))

    rds_files = tibble(file = list.files(path, recursive=TRUE, pattern='rds', full.names=TRUE)) %>%
        mutate(name = basename(file)) %>%
        dplyr::filter(!str_detect(name, 'FeatureCombination')) %>%
        mutate(
            type = case_when(
                str_detect(name, 'ChIPQC')              ~ "ChIPQC",
                str_detect(name, 'Annotate_Peaks')      ~ "Annotate_Peaks",
                str_detect(name, 'Peak_Statistics.rds') ~ "Peak_Statistics",
                str_detect(name, 'Extract_Signal_Annotation.rds') ~ "Extract_Signal_Annotation"
            )) %>%
        mutate(name = str_replace(basename(name),'.rds',''))            %>%
        mutate(name = str_replace(basename(name),'_ChIPQC',''))         %>%
        mutate(name = str_replace(basename(name),'.Annotate_Peaks','')) %>%
        mutate(name = str_replace(basename(name),'.Extract_Signal_Annotation','')) %>%
        dplyr::filter(!is.na(type))
    rds_files
}



# ---------------------------------------------------------------------------- #
format_Annotate_Peaks = function(tab){

    ld = list()
    for(i in 1:nrow(tab)){
        name = tab$name[i]
        message(name)
        granges = readRDS(tab$file[i])
        if(length(granges) > 0)
            ld[[name]] = data.frame(sample_name=name, annot=granges$annot)
    }
    dd = as.data.frame(do.call(rbind, ld)) %>%
        group_by(sample_name, annot)       %>%
        dplyr::summarize(counts = n())     %>%
        mutate(annot = tolower(annot))     %>%
        dplyr::group_by(sample_name) %>%
        mutate(peak_number = length(sample_name)) %>%
        mutate(freq = round(counts/peak_number,2))


    return(dd)
}

# ---------------------------------------------------------------------------- #
format_ChIPQC = function(tab){

    ld = list()
    lres = lapply(setNames(tab$file, tab$name), readRDS)

    lout = list()
    message('GC ...')
      lout$GC = do.call(rbind, lapply(names(lres), function(x){
          dat = lres[[x]]$tilling_windows
          tibble(sample_name = x, GC = dat$G + dat$C, counts = dat$samplecounts[[2]])
      }))

    message('Shift Correlation ...')
    lout$shift_correlation = do.call(rbind, lapply(names(lres), function(x){
        scor = lres[[x]]$ShiftsCorAv
        tibble(sample_name = x, shift = seq(scor), strand_correlation = scor)
    }))

    message('Shift Average ...')
    lout$shift_average = do.call(rbind, lapply(names(lres), function(x){
        scor = lres[[x]]$ShiftsAv
        tibble(sample_name = x, shift = seq(scor), strand_correlation = scor)
    }))

    message('Counts ...')
    lout$counts = Reduce(left_join, lapply(names(lres), function(x){
        cnts = lres[[x]]$tilling_windows$samplecounts[[1]]
        tab = tibble(window = seq(cnts), counts = cnts)
        colnames(tab)[2] = x
        tab
    })) %>%
        mutate(window=NULL)

    message('Normalized Counts ...')
    lout$norm.counts = Reduce(left_join, lapply(names(lres), function(x){
        cnts = lres[[x]]$tilling_windows$samplecounts[[2]]
        tab  = tibble(window = seq(cnts), counts = cnts)
        colnames(tab)[2] = x
        tab
    })) %>%
    mutate(window=NULL)


    message('Annotation ...')
    lout$Annotation = do.call(rbind, lapply(lres, '[[','annot'))
    #
    # lout$readlength = do.call(rbind,lapply(names(lres), function(x){
    #     tibble(sample_name = x, readlength = lres[[x]]$readlength)
    # }))

    return(lout)
}

# ---------------------------------------------------------------------------- #
format_Peak_Statistics = function(tab){
    readRDS(tab$file)
}


# ---------------------------------------------------------------------------- #
format_Extract_Signal_Annotation = function(tab){

    library(genomation)
    lrds = lapply(setNames(tab$file, tab$name), readRDS)

    lsml = lapply(lrds, '[[', 'lsml')
    snams = names(lsml[[1]])
    lsml = lapply(setNames(snams,snams), function(x)ScoreMatrixList(lapply(lsml, '[[', x)))
    lsml = lapply(lsml, intersectScoreMatrixList)
    profiles = do.call(rbind, do.call(rbind, lapply(lrds, '[[', 'profiles')))


    lout = list(
      profiles = profiles,
      lsml     = lsml
      )
}
