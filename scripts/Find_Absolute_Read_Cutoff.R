# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Find_Absolute_Read_Cutoff')

# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff = function(
    infile  = NULL,
    outfile = NULL,
    cutoff  = 5000
){
  if(is.null(infile))
    stop('infile not specified')

  if(is.null(outfile))
      stop('outfile not specified')

  suppressPackageStartupMessages({
      library(stringr)
      library(yaml)
      library(ggplot2)
  })
  
  reads_by_cell = read.table(infile, stringsAsFactors = FALSE)
  reads_by_cell = reads_by_cell[order(-reads_by_cell$V2),]

    cutoff = min(c(cutoff,nrow(reads_by_cell)))
               
 
    message('Calculate cutoff ...')
        cumsum = cumsum(reads_by_cell[1:cutoff, 2])
        df = data.frame(
            "cum"   = cumsum/max(cumsum),
            "cells" = 1:cutoff
         )
        df$rank = 1:nrow(df)
        df$perc = df$rank/max(df$rank)
        df$dist = df$cum - df$perc
    
        knee.point = which.max(df$dist)
    
    message('Plot inflection point ...')
        g = (ggplot(df, aes(cells, cum)) + geom_line(col="steelblue", size=1.25) + theme_minimal()
                + scale_x_continuous(expand=c(0.015, 0))
                + scale_y_continuous(expand = c(0.01, 0)) 
                + ylab("Cumulative fraction of reads")
                + xlab("Cell barcodes (descending number of reads)")
                + theme(text=element_text(size=24),
                      plot.margin = unit(c(1, 1 , 0.5, 0.5), "cm"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_blank()) )
        g = (g 
            + geom_vline(xintercept = knee.point, col='red', size=1)
            + ggtitle(paste0('Number of STAMPS: ', knee.point))
            + theme(title = element_text(size=16)))
        
    png(str_replace(outfile,'yaml','png'), width=400, height=300)
        print(g)
    dev.off()


    cutoff = min(cutoff, nrow(reads_by_cell))
  message('Print output yaml ...')
    
    if(knee.point == cutoff){
      message('knee cell selection did not succeed: including all cells')
        knee.point = min(reads_by_cell[,2])
    }
    lout       = list(reads_cutoff = reads_by_cell[knee.point,2])
    cat(as.yaml(lout), file=outfile)
}


# -------------------------------------------------------------------------- #
Find_Absolute_Read_Cutoff(
      infile  = argv$input[['infile']],
      outfile = argv$output[['outfile']],
      cutoff  = argv$params[['cutoff']]
  )
