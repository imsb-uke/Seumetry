#' Computes differential abundances between conditions using edgeR
#'
#' This function will use [edgeR::glmFit()] to model number of cells per group given attribute/feature x condition (e.g. seurat_clusters x treatment) and [edgeR::glmLRT()] to make contrast between all conditions. Returns the topTags results corrected for multiple comparison using BH and sorted by pvalue.
#'
#' @param cellcounts Matrix with cellcounts per attribute.
#' @param attribute Attribute used to calculate differential abundance. Must be present in meta.data of Seurat Object. Default: seurat_clusters.
#' @param group_by Attribute used to group cells. Default: sample_id.
#' @param formula Formula for GLM, e.g. "~condition1+condition2+condition3". Conditions must be present as attritubes in seu@meta.data.
#' @param contrast Contrast for coefficients, e.g. c(0,1,-1). To determine how to set the desired contrast, use check_coeff = TRUE.
#' @param check_coeff Return coefficients to pick contrasts.
#' @return Data.frame containing topTable results of contrast, corrected for multiple comparison using BH.
#' @export
#' @examples
#'' # first check the coefficiencs
#' differential_abundance(seu, attribute = "seurat_clusters", group_by = "sample_id", formula = as.formula("~0+condition"), check_coef = TRUE)
#' '# calculate differential abundance
#' res <- differential_abundance(seu, attribute = "seurat_clusters", group_by = "sample_id", formula = as.formula("~0+condition"), contrast = c(1, -1))
#' res
differential_abundance <- function(seu,
                                   attribute = "seurat_clusters",
                                   group_by = "sample_id",
                                   formula,
                                   contrast,
                                   check_coeff = FALSE) {
  # create matrix with cellcounts
  cellcounts <- as.matrix(table(seu@meta.data[, attribute], seu@meta.data[, group_by]))
  # create sample-level metadata
  metadata <- seu@meta.data[match(colnames(cellcounts), seu@meta.data[, group_by]),]
  # make model matrix
  model <- model.matrix(formula, data = metadata)
  # dge analysis
  dge <- edgeR::DGEList(counts = cellcounts, lib.size = colSums(cellcounts))
  dge <- edgeR::estimateDisp(dge, model, trend.method = "none")
  fit <- edgeR::glmFit(dge, model)
  # return coefficients to pick contrast
  if(isTRUE(check_coeff)) return(colnames(fit$coeff))
  else{
      # make comparison and return results
      res <- edgeR::glmLRT(fit, contrast = contrast)
      res <- as.data.frame(edgeR::topTags(res, adjust.method = "BH", sort.by = "PValue", n = Inf))
      return(res)
  }
}
