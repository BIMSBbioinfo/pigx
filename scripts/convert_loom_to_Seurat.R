# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')


# ---------------------------------------------------------------------------- #
loomToSeurat = function(
    loom_file       = NULL,
    outfile         = NULL,
    genome_version  = NULL,
    script_path     = NULL

){

    # -------------------------------------------------------------------------#
    suppressPackageStartupMessages({
        library(Seurat)
        library(Matrix)
        library(stringr)
    })
    source(file.path(script_path, 'loom_Functions.R'))

    
    

    # -------------------------------------------------------------------------#
    loom = loom2sce(loom_file)

    meta.data = as.data.frame(SummarizedExperiment::colData(loom))
    colnames(meta.data) = 'CellID'
    raw.data  = Matrix(assays(loom)$counts, nrow=nrow(loom), ncol=ncol(loom))

    rownames(raw.data) = as.character(rowData(loom)$gene_id)
    colnames(raw.data) = as.character(meta.data$CellID)
    seu = CreateSeuratObject(counts    = raw.data, 
                             meta.data = meta.data)

    message('Normalize ...')
        seu = NormalizeData(seu)

    message('Scale ...')
        # seu = ScaleData(object = seu)

    message('Variable genes ...')
        seu = FindVariableFeatures(object = seu, do.plot = FALSE)
  
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
    if(class(rownames(seu)) == 'array'){
        rownames(seu) = as.character(rownames(seu))
    }  
        
    if(class(colnames(seu)) == 'array'){
        colnames(seu) = as.character(colnames(seu))
    }  
        
    seu@meta.data$CellID = colnames(seu)
    saveRDS(seu, outfile)
}



# -------------------------------------------------------------------------- #
loomToSeurat(
    loom_file         = argv$input[['infile']],
    outfile           = argv$output[['outfile']],
    genome_version    = argv$params[['genome_version']],
    script_path       = argv$params[['script']]
)
