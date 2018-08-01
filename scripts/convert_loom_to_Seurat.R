# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')

# ---------------------------------------------------------------------------- #
read_CellCycle = function(
    path_cc  = NULL,
    path_gtf = NULL,
    organism = NULL
){
    
    if(is.null(organism))
        stop('Please specify an organism')
    
    if(organism %in% c('mm','hg')){
        suppressPackageStartupMessages(library(dplyr))
            gtf = rtracklayer::import.gff(path_gtf)
            annot = unique(as.data.frame(values(gtf)[,c('gene_id','gene_name')]))
            cc = read.table(path_cc) %>%
                mutate(cycle = rep(c('S','G2M'), times=c(43,55))) %>%
                `colnames<-`(c('gene_name','cycle'))
            
            if(organims == 'mm'){
                gname = strsplit(cc$gene_name)
                gname = lapply(gname, function(x){
                    s = tolower(x)
                    s[1] = toupper(s[1])
                    paste(s, collapse='')
                })
                cc$gene_name = unlist(gname)
            }
            
            cc = cc %>%
                left_join(annot, by='gene_name') %>%
                arrange(cycle)
    }else{
        cc = NULL
    }
    
    return(cc)
}

# ---------------------------------------------------------------------------- #
loomToSeurat = function(
    loom_file       = NULL,
    outfile         = NULL,
    genome_version  = NULL,
    cycle_genes     = NULL,
    script_path     = NULL
    
){

    # -------------------------------------------------------------------------#
    suppressPackageStartupMessages({   
        library(Seurat)
        library(Matrix)
        library(stringr)
    })
    source(file.path(script_path, 'loom_Functions.R'))
    
    # cell_cycle = read_CellCycle(
    #     path_cc  = cycle_genes,
    #     path_gtf = path_gtf,
    #     organism = str_replace(genomeBuild,'\\d.+','')
    # )
    # 
    # -------------------------------------------------------------------------#
    loom = loom2sce(loom_file)

    meta.data = as.data.frame(SummarizedExperiment::colData(loom))
    colnames(meta.data) = 'CellID'
    raw.data  = Matrix(assays(loom)$counts, nrow=nrow(loom), ncol=ncol(loom))
    rownames(raw.data) = rowData(loom)$gene_id
    colnames(raw.data) = meta.data$CellID
    seu = CreateSeuratObject(raw.data  = raw.data, 
                             meta.data = meta.data)
    
    message('Normalize ...')
        seu = NormalizeData(seu)
        
    message('Scale ...')    
        seu = ScaleData(object = seu)
    
    message('Variable genes ...')
        seu = FindVariableGenes(object = seu, do.plot = FALSE)
  
    # message('PCA ...')
    #     pcs.compute = min()
    #     seu = RunPCA(seu, do.print=FALSE)
    # 
    # message('tSNE ...')
    #     npcs = min(c(20, seu@dr$pca@cell.embeddings))
    #     seu = RunTSNE(object = seu, dims.use = 1:npcs, do.fast = TRUE, check_duplicates = FALSE)
    # 
    # message('Cell cycle ...')
    # if(!is.null(cell_cycle)){
    #     seu = CellCycleScoring(object   = seu,
    #                           s.genes   = s.genes,
    #                           g2m.genes = g2m.genes,
    #                           set.ident = FALSE)
    # }
    seu@meta.data$CellID = rownames(seu@meta.data)
    saveRDS(seu, outfile)
}



# -------------------------------------------------------------------------- #
loomToSeurat(
    loom_file         = argv$input[['infile']],
    outfile           = argv$output[['outfile']],
    genome_version    = argv$params[['genome_version']],
    script_path       = argv$params[['script']]
)






  
