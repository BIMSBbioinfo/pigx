library(GenomicRanges)
args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
target_region <- args[1]
all_intervals_file <- args[2]

#subset all_intervals to only keep the ones that overlap target_region
all_intervals <- as(readLines(all_intervals_file), "GRanges")
target <- as(target_region, "GRanges")

subset_intervals <- IRanges::subsetByOverlaps(all_intervals, target)

outfile <- gsub(".intervals$", '.target.intervals', all_intervals_file)

writeLines(text = paste(subset_intervals), con = outfile)
