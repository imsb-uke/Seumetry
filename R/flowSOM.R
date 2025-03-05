#' Computes FlowSOM- and metaclustering
#'
#' This function computes FlowSOM- and metaclustering using the FlowSOM packages and writes SOM cluster and metacluster information into meta.data of Seurat object. Data should be compensated and transformed.
#'
#' @param seu A seurat object.
#' @param assay Assay from Seurat object to be used for conversion. If empty, uses DefaultAssay.
#' @param slot Slot from Seurat object to be used for conversion. If empty, uses "data" slot.
#' @param features Indicate which features should be used for FlowSOM clustering. Ideally, the features are identical to the features used for dimensionality reduction. Default: NULL, which means that all features are used.
#' @param metaclusters The number of metaclusters to be computed. Default: 10.
#' @param xdim The X dimension of the self-organizing map. Default: 10.
#' @param ydim The Y dimension of the self organizing map. Default: 10.
#' @param ... Additional parameters passed on to FlowSOM::BuildSOM() function.
#' @return Seurat Object containing FlowSOM- and metaclusters in seu@meta.data.
#' @export
#' @examples
#' seu <- run_FlowSOM(seu)
run_FlowSOM <- function(seu,
                        assay = NULL,
                        slot = "data",
                        features = NULL,
                        metaclusters = 10,
                        xdim = 10,
                        ydim = 10,
                        ...) {
    # set assay
    if(is.null(assay)) assay <- DefaultAssay(seu)
    # prepare list of features
    if(is.null(features)) features <- row.names(seu)
    # convert to flowframe
    ff <- convert_seurat(seu, to = "FF", assay = assay, slot = slot)
    # convert to flowSOM
    fSOM <- ReadInput(ff)
    # build SOM
    fSOM <- BuildSOM(fSOM, colsToUse = features, xdim = xdim, ydim = ydim, silent = TRUE, ...)
    # metaclustering
    metacl <- as.character(metaClustering_consensus(fSOM$map$codes, k = metaclusters))
    # add data to Seurat object
    seu$SOM_cl <- factor(GetClusters(fSOM))
    seu$SOM_metacl <- factor(GetMetaclusters(fSOM, metacl), levels = 1:length(metacl))
    # return Seurat object
    return(seu)
}