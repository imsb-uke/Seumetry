#' Manually remove noisy events using positive and/or negative thresholds
#'
#' Remove cells below or above certain threshold for each channel. These can be artefacts with very bright or dim fluorescence.
#'
#' @param seu Seurat object.
#' @param thresholds Named vector of thresholds per marker.
#' @param negative Negative or positive thresholds, i.e. remove events below or above the indicated value. Default: TRUE (= remove events below values).
#' @return Seurat object with cells removed that are below or above threshold.
#' @export
#' @examples
#' # Set positive and negative threshold for each channel
#' thresholds_neg = c("CD3" = -3, "CD4" = -2)
#' thresholds_pos = c("CD3" = 4, "CD4" = 2)
#' # Remove events below indicated thresholds
#' seu_manual = remove_outliers_manual(seu, thresholds)
#' # Remove events above indicated thresholds
#' seu_manual = remove_outliers_manual(seu_manual, thresholds, negative = FALSE)
remove_outliers_manual <- function(seu,
                                   thresholds,
                                   negative = TRUE) {
  remove_cells <- c()
  for(marker in names(thresholds)) {
    # get cell names of cells not passing threshold
    if(isTRUE(negative)) cells <- names(which(seu@assays$comp@data[marker,] < thresholds[[marker]]))
    if(!isTRUE(negative)) cells <- names(which(seu@assays$comp@data[marker,] > thresholds[[marker]]))
    # add cellnames to the array
    remove_cells <- c(remove_cells, cells)
  }
  # make unique
  remove_cells <- unique(remove_cells)
  # actually remove cells
  seu <- seu[,!colnames(seu) %in% remove_cells]
  return(seu)
}


#' Labels potentially noisy or outlier events in Seurat object (seu$outlier)
#' 
#' Uses isolation forest to score events based on their probability to be an outlier. Scores and predicted outliers are stored in Seurat meta.data "outlier" and "outlier_score".
#' 
#' @param seu Seurat object.
#' @param score_threshold The threshold by which an event is called outlier based on isolation forest score. Usually between 0.6 - 0.8. Default: 0.7.
#' @param ndim Isolation tree number of dim. Also see ?isolation.forest. Default: NULL.
#' @param ntrees Isolation tree number of trees. Also see ?isolation.forest. Default: 100.
#' @param seed Set seed for reproducible results. Default: 42.
#' @param ... Additional parameters passed on to isotree::isolation.forest function.
#' @return Seurat object with new cell-level meta.data "outlier" and "outlier_score".
#' @export
#' @examples
#' # score and predict outliers based on threshold
#' seu <- detect_outliers(seu, threshold = 0.6)
#' # remove outliers from Seurat object
#' seu <- subset(seu, subset = outlier == FALSE)
detect_outliers <- function(seu,
                            score_threshold = 0.7,
                            ndim = NULL,
                            ntrees = 100,
                            seed = 42,
                            ...) {
    # get expression matrix
    iso_forest = t(GetAssayData(seu))
    # isolation forest model
    model <- isotree::isolation.forest(iso_forest, ndim = ndim, ntrees = ntrees, ...)
    # label cells in Seurat object
    seu$outlier_score <- predict(model, iso_forest)
    seu$outlier <- seu$outlier_score > score_threshold
    # feedback on how many outliers found
    message(paste("Found",
                  sum(seu$outlier),
                  "outliers, which is",
                  round((sum(seu$outlier) / ncol(seu)) * 100, 2),
                  "percent of all cells."))
    return(seu)
}
