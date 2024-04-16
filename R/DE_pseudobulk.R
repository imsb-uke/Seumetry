#' Computes differentially expressed proteins based on median fluorescence intensity per group using limma
#'
#' This function calculates the median fluorescence intensity per group using [median_expression()] to create a "pseudobulk" matrix. To detect differentially expressed proteins, limma is used: [limma::lmFit()] to fit the linear model and [limma::contrasts.fit()], [limma::eBayes()] to make contrasts. Returns data.frame of [limma::topTable()] results.
#'
#' @param seu Seurat object.
#' @param group_by Indicate group which is used for pseudobulking using [median_expression()]. Default: sample_id.
#' @param fixed_vars Fixed variables for the linear model.
#' @param contrast Indicate which contrast to make, also see limma vignette. For example, attribute1condition1-attribute1condition2. Second mention is always the reference.
#' @param min.cells Minimum number of cells that should be in a group. If less cells, the group is dropped.
#' @return Data.frame containing [limma::topTable()] results.
#' @export
#' @examples
#' # compute differential expression
#' DE_res <- DE_pseudobulk(seu, fixed_vars = "condition", contrast = "conditionA-conditionB")
DE_pseudobulk <- function(seu,
                          group_by = "sample_id",
                          fixed_vars,
                          contrast,
                          min.cells = 3) {
    # get median expression per sample
    counts <- median_expression(seu, group_by = group_by)
    # create metadata data.frame based on seu@meta.data
    metadata <- seu@meta.data[match(colnames(counts), seu@meta.data[, group_by]),]
    # filter out samples with low number of cells
    keep <- names(which(table(seu@meta.data[, group_by]) >= min.cells))
    counts <- counts[,keep]
    metadata <- metadata[which(metadata[[group_by]] %in% keep),]
    # create formula
    formula <- formula(paste0("~0+", paste(fixed_vars, collapse="+")))
    # create design
    design <- model.matrix(formula, metadata)
    # fit model
    fit <- limma::lmFit(counts, design)
    # make contrast
    contrast.matrix <- do.call(limma::makeContrasts, args = list(contrast, levels = design))
    # fit contrast
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    # ebayes test
    fit2 <- limma::eBayes(fit2)
    # toptable
    res <- limma::topTable(fit2, n = Inf)
    # return
    return(res)
}
