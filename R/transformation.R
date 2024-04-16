#' Biexponential transformation based on flowWorkspace::flowjo_biexp
#'
#' Transforms values for every channel based on panel data frame. Biexp can be used with 3 parameters (columns) per channel based on FlowJo 10: biexp_pos, biexp_neg, biexp_width. If a column is not provided, function will use default values. See [flowWorkspace::flowjo_biexp()] function for default values. It's best to use the wrapper function [transform_data()] instead of using this function directly.
#'
#' @param matrix Provide matrix of flow data. E.g. "Seurat::GetAssayData(seu)".
#' @param panel Provide data.frame with panel and transformation information. See vignette for details how to generate a panel_df.
#' @return Transformed matrix that can be added to Seurat object.
#' @export
#' @examples
#' seu <- SetAssayData(seu, slot = "data", new.data = transform_biexp(matrix, panel))
transform_biexp <- function(matrix,
                            panel) {
  # function to read settings from panel
  do_transform <- function(values,
                           marker,
                           panel) {
    # panel subset
    panel_sub <- panel[which(panel$antigen == marker),]
    # get args from panel
    args <- as.list(formals(flowjo_biexp))
    if("biexp_pos" %in% names(panel_sub)) args["pos"] <- panel_sub$biexp_pos
    if("biexp_neg" %in% names(panel_sub)) args["neg"] <- panel_sub$biexp_neg
    if("biexp_width" %in% names(panel_sub)) args["widthBasis"] <- panel_sub$biexp_width
    if("biexp_max" %in% names(panel_sub)) args["maxValue"] <- panel_sub$biexp_max
    if("biexp_range" %in% names(panel_sub)) args["channelRange"] <- panel_sub$biexp_range
    # build trans function
    trans <- do.call(flowjo_biexp, args)
    # apply trans function
    new_values <- trans(values)
    return(new_values)
  }
  # transform each channel individually
  new_matrix <- t(sapply(row.names(matrix), FUN = function(marker) do_transform(matrix[marker,], marker, panel)))
  # add cellnames again
  colnames(new_matrix) <- colnames(matrix)
  return(new_matrix)
}


#' Arcsinh transformation using cofactors for each individual channel
#'
#' Transforms values for every channel based on cofactor provided in panel data.frame (arcsinh_cofactor). If column not present, uses default cofactor: 5. Function used: y = arcsinh(x/cofactor). It's best to use the wrapper function [transform_data()] instead of using this function directly.
#'
#' @param matrix Provide matrix of flow data. E.g. "Seurat::GetAssayData(seu)".
#' @param panel Provide data.frame with panel and transformation information. See vignette for details how to generate a panel_df.
#' @return Transformed matrix that can be added to Seurat object.
#' @export
#' @examples
#' seu <- SetAssayData(seu, slot = "data", new.data = transform_arcsinh(matrix, panel))
transform_arcsinh <- function(matrix,
                              panel) {
  # function to read settings from panel
  do_transform <- function(values,
                           marker,
                           panel) {
    # panel subset
    panel_sub <- panel[which(panel$antigen == marker),]
    # get arguments for asinh
    if("arcsinh_cofactor" %in% names(panel_sub)) {
      cofactor <- panel_sub$arcsinh_cofactor
    }else{
      cofactor <- 5
    }
    # apply arcsinh function
    new_values <- asinh(values/cofactor)
    return(new_values)
  }
  # transform each channel individually
  new_matrix <- t(sapply(row.names(matrix), FUN = function(marker) do_transform(matrix[marker,], marker, panel)))
  # add cellnames again
  colnames(new_matrix) <- colnames(matrix)
  return(new_matrix)
}


#' Wrapper for all transformation functions using Seurat object as input & output
#'                         
#' See transformation functions for detail: [transform_arcsinh()], [transform_biexp()]. Input Seurat Object uses "counts" slot from current assay, performs transformation and writes transformed counts into "data" slot from current assay using panel information stored in seu@misc.
#' 
#' @param seu Seurat object.
#' @param transformation Indicate which transformation to use (arcsinh, biexp)
#' @return Seurat object with transformed data in "data" slot of current assay.
#' @export
#' @examples
#' seu <- transform_data(seu, "arcsinh")
transform_data <- function(seu,
                           transformation) {
  # get untransformed data
  matrix <- GetAssayData(seu, slot = "counts")
  # get panel (stored in misc slot by prep_seurat function)
  panel <- data.frame(seu@misc)
  # transform data
  if(transformation == "arcsinh") transformed <- transform_arcsinh(matrix, panel)
  if(transformation == "biexp") transformed <- transform_biexp(matrix, panel)
  # write into Seurat object
  seu <- SetAssayData(seu, slot = "data", new.data = transformed)
  return(seu)
}
