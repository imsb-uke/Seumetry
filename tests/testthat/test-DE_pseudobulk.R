test_that("DE_pseudobulk returns a data.frame with expected limma columns", {
  seu <- make_da_seurat()
  res <- suppressWarnings(DE_pseudobulk(seu, fixed_vars = "condition",
                                        contrast = "conditionctrl-conditiontreat"))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(res)))
})

test_that("DE_pseudobulk drops samples below min.cells without erroring", {
  seu <- make_da_seurat(n_cells = 200)
  # shrink S1 to only 1 cell (below min.cells) by reassigning the rest to S2
  idx <- which(seu@meta.data$sample_id == "S1")
  seu@meta.data$sample_id[idx[-1]] <- "S2"
  res <- suppressWarnings(DE_pseudobulk(seu, fixed_vars = "condition",
                                        contrast = "conditionctrl-conditiontreat",
                                        min.cells = 3))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 6)  # 6 features, all samples present bar the dropped S1
})

test_that("DE_pseudobulk errors when fixed_vars is missing", {
  seu <- make_da_seurat()
  expect_error(DE_pseudobulk(seu, contrast = "conditionctrl-conditiontreat"),
               "fixed_vars must be provided")
})

test_that("DE_pseudobulk errors when contrast is missing", {
  seu <- make_da_seurat()
  expect_error(DE_pseudobulk(seu, fixed_vars = "condition"),
               "contrast must be provided")
})

test_that("DE_pseudobulk errors when group_by column does not exist", {
  seu <- make_da_seurat()
  expect_error(DE_pseudobulk(seu, group_by = "nonexistent", fixed_vars = "condition",
                             contrast = "conditionctrl-conditiontreat"),
               "group_by column 'nonexistent' not found")
})

test_that("DE_pseudobulk errors when fixed_vars column does not exist", {
  seu <- make_da_seurat()
  expect_error(DE_pseudobulk(seu, fixed_vars = "nonexistent",
                             contrast = "conditionctrl-conditiontreat"),
               "not found in seu@meta.data: nonexistent")
})

test_that("DE_pseudobulk errors when min.cells filtering leaves fewer than 2 groups", {
  seu <- make_da_seurat()
  expect_error(DE_pseudobulk(seu, fixed_vars = "condition",
                             contrast = "conditionctrl-conditiontreat",
                             min.cells = 1000),
               "Fewer than 2 groups")
})
