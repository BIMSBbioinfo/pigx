# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('runFeatureCounts', parseFlags = TRUE)

# ---------------------------------------------------------------------------- #
## Count reads using Rsubread
runFeatureCounts <- function(peaks,
                             bam_files,
                             outRDS,
                             tempDir = tempdir("featureCounts_temp"),
                             ..., 
                             blacklisted = NULL,
                             n.cores = 1
                             ) {

    suppressPackageStartupMessages({
        require(Rsubread)
        require(GenomicRanges)
        require(rtracklayer)
    })

    if(!all(file.exists(bam_files))){
        stop(ipaste("Given bam files do not exist.","\n\t",
                    head(bam_files[!file.exists(bam_files)])
                    ))
    }

    fc_args <- list(...)
    
    message(paste("Started at ",timestamp(quiet = TRUE)))

    ## load peaks if the argument is a file
    if(is.character(peaks)) {
        if(!file.exists(peaks)){
            stop(paste("Given peak file do not exist.","\n\t",peaks))
        }
        if(grepl("txt$",peaks)) {
            peaks <- read.delim(peaks)
            names(peaks) <- c("chr","start","end", "name","class")
            peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks,
            keep.extra.columns = TRUE)
        } else {
            peaks <- rtracklayer::import(peaks)
        }
    }

    ## load blackliste regions
    if (!is.null(blacklisted)) {
        if(is.character(blacklisted)) {
            blacklisted <- rtracklayer::import.bed(blacklisted)
            peaks <- peaks[!peaks %over% blacklisted]
        }
    }


    if(!"name" %in% names(mcols(peaks)))  {
        peaks$name <- paste0("cons_peak_",seq_len(nrow(peaks)))
    }
    peaks_saf <- as.data.frame(peaks)[c("name","seqnames","start","end","strand")]
    names(peaks_saf) <- c("GeneID","Chr","Start","End","Strand")

    message("Counting ...")
    message(timestamp(quiet = TRUE))
            
    if(!dir.exists(tempDir)) dir.create(tempDir)    

    ## start counting
    fc <- Rsubread::featureCounts(bam_files,
                                  annot.ext = peaks_saf,
                                  ...,
                                  tmpDir = tempDir,
                                  nthreads = n.cores)
    fc$targets <- basename(bam_files)
    saveRDS(fc,outRDS) 

    message("Counting done.")
    message(timestamp(quiet = TRUE))

    message(paste("Finished at ",timestamp(quiet = TRUE)))
}

# ---------------------------------------------------------------------------- #
# function call

## debug
# print(argv)
## END debug

runFeatureCounts(
                 peaks = argv$input[['bedfile']],
                 bam_files =  argv$input[['bamfile']],
                 outRDS = argv$output[["outfile"]],
                 tempDir = dirname(argv$output[["outfile"]]),
                 blacklisted = NULL,
                 n.cores = argv$params[["threads"]],
                 argv$params[["params_tool"]]
)
