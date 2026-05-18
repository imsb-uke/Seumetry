test_that("differential_abundance returns coefficient names with check_coeff = TRUE", {
  skip_if_not_installed("edgeR")
  # Needs ≥3 samples so estimateDisp has residual df; use the 4-sample fixture
  seu <- make_da_seurat()
  coeff <- suppressWarnings(differential_abundance(
    seu,
    attribute   = "seurat_clusters",
    group_by    = "sample_id",
    formula     = as.formula("~0 + condition"),
    check_coeff = TRUE
  ))
  expect_type(coeff, "character")
  expect_true(length(coeff) >= 1)
})

test_that("differential_abundance returns a data.frame with expected columns", {
  skip_if_not_installed("edgeR")
  seu <- make_da_seurat()
  result <- suppressWarnings(differential_abundance(
    seu,
    attribute = "seurat_clusters",
    group_by  = "sample_id",
    formula   = as.formula("~0 + condition"),
    contrast  = c(1, -1)
  ))
  expect_s3_class(result, "data.frame")
  expect_true("logFC"  %in% colnames(result))
  expect_true("FDR"    %in% colnames(result))
  expect_true("PValue" %in% colnames(result))
  expect_equal(nrow(result), 4)   # one row per cluster
})
