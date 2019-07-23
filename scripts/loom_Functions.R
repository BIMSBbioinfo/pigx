array2char = function(df){
    df = DataFrame(lapply(df, function(x){
        if(class(x) == 'array')
            x = as.character(x)
        return(x)
    }))
    return(df)
}

loom2sce = function(
    path = NULL
) {

    if(is.null(path))
        stop('loom file does not exist')

    suppressPackageStartupMessages({
        library(HDF5Array)
        library(SingleCellExperiment)
    })


    message("Importing",path,"into SingleCellExperiment object")
    message("Extracting the count matrix")
    mat = HDF5Array(path, "matrix")
    mat = t(mat)

    message("Extracting the row and column metadata")
    col.attrs = h5read(path, "col_attrs")
    if (length(col.attrs)) {
        col.df = data.frame(col.attrs)
    } else {
        col.df = DataFrame(matrix(0, nrow(mat), 0))
    }
    col.df = array2char(col.df)

    row.attrs = h5read(path, "row_attrs")
    if (length(row.attrs)) {
        row.df = data.frame(row.attrs)
    } else {
        row.df = NULL
    }
    row.df = array2char(row.df)
    rownames(mat) = row.df$gene_id
    colnames(mat) = col.df$cell_id

    message("Extracting layers (if there are any)")
    optional = h5ls(path)
    is.layer = optional$group=="/layers"
    if (any(is.layer)) {
        layer.names = optional$name[is.layer]
        other.layers = vector("list", length(layer.names))
        names(other.layers) = layer.names

        for (layer in layer.names) {
            current = HDF5Array(path, file.path("/layers", layer))
            current = t(current)
            rownames(current) = rownames(mat)
            colnames(current) = colnames(mat)
            other.layers[[layer]] = current
        }
    } else {
        other.layers = list()
    }
    message("Returning SingleCellExperiment object.")
    sce = SingleCellExperiment(c(counts=mat, other.layers), rowData=row.df, colData=col.df)
    colnames(rowData(sce))[1] = 'gene_id'
    return(sce)
}
