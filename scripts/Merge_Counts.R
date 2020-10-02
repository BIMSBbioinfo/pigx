
# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('mergeFeatureCounts')

# ---------------------------------------------------------------------------- #
## merge count files produced by FeatureCounts
mergeFeatureCounts <- function(countFiles, outFile, statFile) {

    message(paste("Started at ",timestamp(quiet = TRUE)))

    ### merge counting results
    message("Merge counts ...")
    message(timestamp(quiet = TRUE))

    fc <- lapply(countFiles,readRDS)

    names(fc) <- sapply(fc, function(x) x$targets)

    counts <- data.frame(sapply(fc,function(x) x$counts),row.names = rownames(fc[[1]]$counts))

    write.table(x = counts,
                file = outFile,
                sep = "\t",
                quote = FALSE)

    countStats <- data.frame("Status"=fc[[1]]$stat$Status,sapply(fc,function(x) x$stat[,2]))
    write.table(x = countStats,
                file = statFile,
                sep = "\t",
                quote = FALSE)

    message("Merge counts done.")
    message(timestamp(quiet = TRUE))

    message(paste("Finished at ",timestamp(quiet = TRUE)))

}

# ---------------------------------------------------------------------------- #
# function call

## debugging
# print(argv)
## END debugging

mergeFeatureCounts(
                   countFiles = argv$input[["countfiles"]],
                   outFile    = argv$output[["outfile"]],
                   statFile   = argv$output[["statfile"]]
)
