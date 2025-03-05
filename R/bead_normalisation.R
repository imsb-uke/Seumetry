#' Cytof bead normalization via CATALYST
#'
#' Applies a bead normalisation on raw intensities (counts slot) of active Seurat assay using CATALYST::normCytof. A Seurat object is returned with a new assay "beadnorm" and 2 new meta.data columns for beads and potential bead/cell doublets. For details how this normalization works, please refer to https://www.bioconductor.org/packages/release/bioc/html/CATALYST.html.
#'
#' @param seu A seurat object.
#' @param assay Assay from Seurat object to be used for conversion. If empty, uses DefaultAssay.
#' @param beads Indicate which beads were used (dvs or beta) or a custom numeric vector of bead masses (see CATALYST function).
#' @param plot_res Plots quality control graphs. Default: TRUE
#' @param ... Additional parameters passed on to CATALYST::normCytof() function (see CATALYST function).
#' @return Seurat Object with additional assay (beadnorm) and 2 additional meta.data columns indicating beads and bead_doublets with cells that can be removed.
#' @export
#' @examples
#' seu <- bead_norm(seu, beads = "dvs", k = 50) # Run bead normalisation
#' seu <- subset(seu, subset = bead_doublet == FALSE) # Remove beads and bead/cell doublets
bead_norm <- function(seu,
                      assay = NULL,
                      beads = "dvs",
                      plot_res = TRUE,
                      ...) {
    # get expression data
    if(is.null(assay)) assay <- DefaultAssay(seu)
    # convert and add data assay as exprs
    sce <- convert_seurat(seu, "SCE", slot = "counts")
    assay(sce, "data") <- counts(convert_seurat(seu, "SCE", slot = "data"))
    # add unused to int_colData
    int_colData(sce) <- cbind(int_colData(sce),
			      t(as.data.frame(GetAssayData(seu, assay = "unused", layer = "counts"))))
    # fix fcs_colname and antigen rowData columns to match with CATALYST workflow
    colnames(rowData(sce))[1:2] <- c("channel_name", "marker_name")
    # compute CATALYST normCytof
    res <- normCytof(sce,
                     beads = beads,
                     assays = c("counts", "data"),
                     overwrite = FALSE,
                     transform = FALSE,
                     remove_beads = FALSE,
                     cofactor = NA,
                     ...)
    # write everything into Seurat object
    seu[["beadnorm"]] <- CreateAssayObject(res$data@assays@data$normcounts)
    DefaultAssay(seu) <- "beadnorm"
    seu$bead <- colData(res$data)$is_bead
    seu$bead_doublet <- colData(res$data)$remove
    # plot results
    if(plot_res) {
        print(res$lines)
        print(res$scatter)
    }
    # return Seurat Object
    return(seu)
    
}