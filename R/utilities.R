#' Calculates median fluorescence intensity per group based on "data" slot in DefaultAssay(seu)
#'
#' @param seu Seurat object.
#' @param group_by Indicate a column in meta.data of Seurat Object by which the data should be grouped. Default: Idents(Seu)
#' @return Matrix with median expression values across different groups.
#' @export
#' @examples
#' sample_mfi <- median_expression(seu, group_by = "sample_id")
median_expression <- function(seu,
                              group_by = "ident") {
    # add idents as metadata column
    if(group_by == "ident") seu$ident <- Idents(seu)
    # retrieve transformed counts from Seurat object (DefaultAssay)
    counts <- GetAssayData(seu, slot = "data")
    # rename colnames in counts by group.by
    colnames(counts) <- seu@meta.data[[group_by]]
    # get the median per group
    median <- t(apply(counts, 1, function(x) tapply(x, colnames(counts), median)))
    # return median matrix
    return(median)
}


#' Exports FCS files from Seurat object.
#'
#' This function uses a specified assay and slot from a Seurat object, converts the Object to a flowFrame, and exports the flowFrame as a FCS file. This can be used, for example, to visualize specific clusters or cells in a traditional Flow Cytometry Software.
#'
#' @param seu A seurat object.
#' @param assay Assay from Seurat object to be used for conversion. If empty, uses DefaultAssay.
#' @param slot Slot from Seurat object to be used for conversion. If empty, uses "data" slot.
#' @param filename Path and filename where FCS file should be written to.
#' @export
#' @examples
#' # subset Seurat object to specific cluster
#' seu_sub <- subset(seu, subset = seurat_clusters == 1)
#' # export FCS file containing cells from this cluster using raw data
#' export_fcs(seu_sub, filename = "cluster_1.fcs", assay = "fcs", slot = "counts")
export_fcs <- function(seu,
                       filename,
                       assay = NULL,
                       slot = "data") {
    # convert to flowframe
    ff <- convert_seurat(seu, "FF", assay = assay, slot = slot)
    # export fcs
    write.FCS(ff, filename = filename)
    message(paste("FCS file", filename, "written successfully."))
    
}


#' Converts a Seurat object to a FlowFrame, FlowSet or SingleCellExperiment Object.
#'
#' Takes Seurat object as input and returns FlowFrame, FlowSet or SingleCellExperiment Object. These can be used for integration of third-party tools. By default uses DefaultAssay and Slot from Seurat object as Matrix. Only channels used in matrix and present in panel are used for conversion. Unused Channels will not be included.
#'
#' @param seu A seurat object.
#' @param to Can be either FF for FlowFrame, FS for FlowSet, or SCE for SingleCellExperiment.
#' @param assay Assay from Seurat object to be used for conversion. If empty, uses DefaultAssay.
#' @param slot Slot from Seurat object to be used for conversion. If empty, uses "data" slot.
#' @param split_by Indicate which seu@meta.data column to use to split data into individual FlowFrames to be used in a FlowSet.
#' @return Either a FlowFrame, FlowSet, or SingleCellExperiment Object.
#' @export
#' @examples
#' sce = convert_seurat(seu, assay = "comp", to = "SCE")
convert_seurat <- function(seu,
                           to = c("FF", "FS", "SCE"),
                           assay = NULL,
                           slot = "data",
                           split_by) {
    # get expression data
    if(is.null(assay)) assay <- DefaultAssay(seu)
    mtx <- as.matrix(GetAssayData(seu, assay = assay, slot = slot))
    # get panel data
    panel <- as.data.frame(seu@misc)
    # FLOWFRAME
    if(to == "FF"){
        # fix panel data
        row.names(panel) <- paste0("$P", row.names(panel))
        names(panel)[1:2] <- c("desc", "name")
        panel$minRange <- 0
        panel$maxRange <- 0
        panel$range <- 0
        # create annotated data frame
        an_df <- AnnotatedDataFrame(data = panel)
        # create flowframe
        ff <- flowFrame(t(mtx), an_df)
        # set identifier
        identifier(ff) <- "FlowFrame from Seurat"
	# set cell IDs
	row.names(ff@exprs) <- colnames(mtx)
        # return flowframe
        return(ff)
    }
    # FLOWSET
    if(to == "FS"){
        if(missing(split_by)) stop("If converting to flowSet, split_by parameter must be provided.")
        # split Seurat Object
        seu_split <- SplitObject(seu, split.by = split_by)
        # make list of flowframes using recursion
        ff_list <- sapply(names(seu_split), function(x) convert_seurat(seu_split[[x]], to = "FF", assay = assay, slot = slot))
        # generate flowset
        fs <- flowSet(ff_list)
        # return FlowSet
        return(fs)
    }
    # SCE
    if(to == "SCE"){
        # get cell-level metadata
        metadata <- as.data.frame(seu@meta.data)
        # get reductions
        reductions <- list()
        for(reduction in names(seu@reductions)){
            reductions[[reduction]] <- seu@reductions[[reduction]]@cell.embeddings
        }
        # construct SCE
        sce <- SingleCellExperiment(assays = list(counts = mtx), colData = metadata, rowData = panel, reducedDims = reductions)
        # return SingleCellExperiment
        return(sce)
    }
}