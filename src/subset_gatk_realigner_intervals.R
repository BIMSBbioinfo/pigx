library(GenomicRanges)
args = commandArgs(trailingOnly=TRUE)

sample_sheet_file <- args[1]
sample <- args[2]
all_intervals_file <- args[3]

sample_sheet <- read.csv(sample_sheet_file)
sample_sheet <- sample_sheet[sample_sheet$sample_name == sample,]

#subset all_intervals to only keep the ones that overlap target_region
all_intervals <- as(readLines(all_intervals_file), "GRanges")
target <- as(sample_sheet$target_region, "GRanges")

subset_intervals <- IRanges::subsetByOverlaps(all_intervals, target)

outfile <- gsub(".intervals$", '.target.intervals', all_intervals_file)

writeLines(text = paste(subset_intervals), con = outfile)
