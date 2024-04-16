#' Compensate data based on spillover matrix save in FCS files
#'
#' Applies compensation matrices for each channel and saves compensated values into "comp" assay in Seurat Object. Uses [flowCore::compensate()] function.
#'
#' @param fcs_fs FlowSet containing all fcs files and compensation matrices. Must be the same FlowSet that was used for [create_seurat()].
#' @param seu Seurat object used for [create_seurat()].
#' @param comp_matrix Provide a compensation matrix which is the same for all samples. Names of columns and rows must match channel names indicated in panel data.frame.
#' @param fcs_matrix Use compensation matrices from FCS files. Can be different for every file. This parameter indicates which compensation matrix to use if there are multiple in the FCS files. Default: 1.
#' @return Seurat object with new assay "comp", where compensation values were saved in slot "counts".
#' @export
#' @examples
#' ## Option 1) Use external spillover matrix directly: same compensation for all files
#' # colnames must match panel data.frame; row.names are optional
#' comp_matrix <- read.csv("spillover_matrix.csv", row.names = 1)
#' seu_1 <- compensate_data(fcs_fs, seu, comp_matrix = comp_matrix)
#' ## Option 2) Use spillover matrix saved in FCS files: same compensation for all files
#' comp_matrix <- spillover(fcs_fs[[1]])[[3]] # Matrix 3 from first FCS file
#' seu_2 <- compensate_data(fcs_fs, seu, comp_matrix = comp_matrix)
#' ## Option 3) Use spillover matrix saved in FCS files: compensation matrix used from individual FCS files, thus can be different for each file
#' # Check different matrices using:
#' names(spillover(fcs_fs[[1]]))
#' # Pick the correct spillover matrix for compensation and provide with fcs_matrix parameter
#' seu <- compensate_data(fcs_fs, seu, fcs_matrix = 3)
compensate_data <- function(fcs_fs,
                            seu,
                            comp_matrix = NULL,
                            fcs_matrix = 1) {
  # create compensation object using a spillover matrix directly
  if(!is.null(comp_matrix)) {
      comp <- compensation(comp_matrix)
  # or use spillover matrices from FCS files (can be different for every file)
  }else{
      # check some parameters
      if(is.null(spillover(fcs_fs[[1]])[[fcs_matrix]]))
        stop("Spillover matrix empty!", call. = FALSE)
      # create object
      comp <- fsApply(fcs_fs, function(x) spillover(x)[[fcs_matrix]], simplify = FALSE)
  }
  # compensate and save in new flowSet object
  fs_comp <- compensate(fcs_fs, comp)
  # generate Seurat object
  seu_comp <- create_seurat(fs_comp, data.frame(seu@misc))
  # add compensated matrix to new assay
  seu[["comp"]] <- CreateAssayObject(counts = GetAssayData(seu_comp),
                                     min.cells = 0,
                                     min.features = 0)
  # make compensated data as default assay
  DefaultAssay(seu) <- "comp"
  return(seu)
}
