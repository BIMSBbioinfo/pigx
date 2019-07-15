# -------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Extract_Read_Statistics')


# ---------------------------------------------------------------------------- #
convert_loom_to_seurat = function(
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
    raw.data  = Matrix(assays(loom)$counts, nrow=nrow(loom), ncol=ncol(loom))

    rownames(raw.data) = as.character(rowData(loom)$gene_id)
    colnames(raw.data) = as.character(meta.data$cell_id)
    seu = CreateSeuratObject(counts    = raw.data,
                             meta.data = meta.data)

    message('Normalize ...')
        seu = NormalizeData(seu)

    message('Scale ...')
        seu = ScaleData(object = seu)

    seu@meta.data$CellID = colnames(seu)
    saveRDS(seu, outfile)
}



# -------------------------------------------------------------------------- #
convert_loom_to_seurat(
    loom_file         = argv$input[['infile']],
    outfile           = argv$output[['outfile']],
    genome_version    = argv$params[['genome_version']],
    script_path       = argv$params[['script']]
)
