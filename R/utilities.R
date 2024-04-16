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


#' Export specific cells to FCS file
#'
#' This function uses the initially create flowSet and exports original FCS files based on cell_IDs in Seurat object using [flowCore::write.FCS()]. Can be used to view data with external software such as FlowJo, e.g. to verify cluster annotation.
#'
#' @param fcs_fs The flowCore FlowSet that was prepared in the first step during data prep. See [create_flowset()].
#' @param cell_ids Name of cells that should be exported (colnames((seu)).
#' @param filename Path and filename where FCS file should be written to.
#' @return If no filename is indicated, this function will return a flowFrame instead of saving a FCS file.
#' @export
#' @examples
#' # load flowSet
#' fcs_fs <- readRDS("fcs_fs.rds")
#' # subset Seurat object to specific cluster
#' seu_sub <- subset(seu, subset = seurat_clusters == 1)
#' # export FCS file containing cells from this cluster
#' export_fcs(fcs_fs, cell_ids = colnames(seu_sub), filename = "cluster_1.fcs")
export_fcs <- function(fcs_fs,
                       cell_ids = NULL,
                       filename) {
    # merge all expression matrices from flowFrames from flowSet
    merged <- lapply(fcs_fs@frames, function(x) exprs(x))
    merged <- do.call(rbind, merged)
    # get a references flowFrame from flowSet
    fcs_ff <- fcs_fs[[1]]
    # if desired, subset the merged flowSet using specific cellnames
    if(!is.null(cell_ids)) {
        merged <- merged[which(row.names(merged) %in% cell_ids),]
    }
    # replace exprs in this flowFrame by the merged exprs matrix
    fcs_ff@exprs <- merged
    # export FCS
    if(!is.null(filename)) {
        # write FCS (flowCore)
        write.FCS(fcs_ff, filename = filename)
        print(paste("FCS file", filename, "containing", nrow(merged), "cells written successfully."))
    }else{
        # if exporting not desired, just return created flowFrame
        return(fcs_ff)
    }
}
