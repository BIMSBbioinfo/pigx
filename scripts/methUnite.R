# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>,
# Copyright © 2019 Alexander Blume <alexander.blume@mdc-berlin.de>,
# Copyright © 2019 Katarzyna Wreczycka <katarzyna.wreczycka@mdc-berlin.de>
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

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
  args <- c("--help")
}

## Help section
if ("--help" %in% args) {
  cat("
      Merge methylation samples using methylKit

      Arguments:
      --inputfiles comma separated input tabix files
      --sampleids comma separated sample name vector
      --treatments comma separated treatment vector
      --assembly genome assembly
      --context methylation context
      --destrand wether to merge strands
      --cores number of processing cores
      --outdir output directory
      --suffix suffix for merged file
      --logFile file to print the logs to
      --help              - print this text

      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save = "no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

# ## catch output and messages into log file
out <- file(argsL$logFile, open = "wt")
sink(out, type = "output")
sink(out, type = "message")

print(argsL)

# Run Functions -----------------------------------------------------------

suppressPackageStartupMessages(library("methylKit"))
data.table::setDTthreads(8)

## Load variables
inputs <- strsplit(argsL$inputfiles, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
sampleids <- strsplit(argsL$sampleids, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
cores <- as.numeric(argsL$cores)
assembly <- argsL$assembly
suffix <- argsL$suffix
outdir <- argsL$outdir
destrand <- ifelse(tolower(argsL$destrand) %in% c("true", "yes"), TRUE, FALSE)

message("Remapping Treatment Description to number.")
# split all treatment values, could be numeric or not
treatmentsStr <- strsplit(argsL$treatments, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
# convert into named numeric vector
treatments <- as.numeric(as.factor(treatmentsStr))
names(treatments) <- treatmentsStr

message(paste(
  sort(unique(treatmentsStr)),
  sort(unique(treatments)),
  sep = ": ", collapse = "\n"
))

# convert variables and perform checks
context <- switch(tolower(argsL$context),
  cpg = "CpG",
  chh = "CHH",
  chg = "CHG",
  NULL
)

if (is.null(context)) {
  stop(
    "The given context <", argsL$context,
    "> is not among the supported ones ('CpG','CHG','CHH')"
  )
}

if (length(inputs) <= 1) {
  stop("At least two samples required to perform merging")
}


# ---------------------------------------------------
st <- system.time({

  ## Hardlink Tabix files to avoid race conditions
  inputFiles <- basename(inputs)
  linkedFiles <- file.path(outdir, inputFiles)
  file.link(inputs, linkedFiles)
  file.link(paste0(inputs, ".tbi"), paste0(linkedFiles, ".tbi"))
  inputs <- linkedFiles


  if (length(inputs > 1)) {
    inputs <- as.list(inputs)
    sampleids <- as.list(sampleids)
  }

  ## Read data
  message("Importing Tabix files.")
  methylRawListDB <- methRead(
    location = inputs,
    sample.id = sampleids,
    assembly = assembly,
    dbtype = "tabix",
    context = context,
    treatment = treatments
  )

  ## Destrand if needed
  if (destrand) {
    message("Destranding first.")

    destrandFun <- function(obj, cores) {
      print(getSampleID(obj))

      if (obj@resolution == "base") {
        dir <- dirname(obj@dbpath)
        filename <- paste(gsub(".txt.bgz", "", obj@dbpath),
          "destrand.txt",
          sep = "_"
        )

        filename <- basename(filename)

        # need to use .CpG.dinuc.unifyOld because output needs to be ordered
        newdbpath <- methylKit:::applyTbxByChr(obj@dbpath,
          dir = dir, filename = filename,
          return.type = "tabix",
          FUN = function(x) {
            options(scipen = 999)
            methylKit:::.CpG.dinuc.unifyOld(methylKit:::.setMethylDBNames(
              x,
              "methylRawDB"
            ))
          },
          mc.cores = cores
        )

        obj <- methylKit:::readMethylRawDB(
          dbpath = newdbpath, dbtype = "tabix",
          sample.id = obj@sample.id,
          assembly = obj@assembly, context = obj@context,
          resolution = obj@resolution
        )
      }
      return(obj)
    }

    new.list <- suppressMessages(lapply(methylRawListDB, destrandFun, cores = cores))
    methylRawListDB <- new("methylRawListDB", new.list, treatment = methylRawListDB@treatment)
    destrandFiles <- getDBPath(methylRawListDB)
  }

  ## Unite
  message("Merging samples.")
  methylBaseDB <- unite(methylRawListDB,
    destrand = FALSE,
    suffix = suffix,
    dbdir = outdir,
    mc.cores = cores,
    chunk.size = 1e7,
  )

  ## FIXME: check wether result has more than 1 rows and fail if not

  if (destrand) {

    ## Remove temp destrand files
    unlink(c(destrandFiles, paste0(destrandFiles, ".tbi")))
  }

  ## Remove temp hardlinked files
  unlink(c(linkedFiles, paste0(linkedFiles, ".tbi")))
})

message("Done.")
message("Process finished in (seconds): \n")
print(st)
