test_that("detect_outliers adds outlier_score and outlier metadata columns", {
  seu <- make_test_seurat(n_cells = 100)
  result <- suppressMessages(detect_outliers(seu))
  expect_s4_class(result, "Seurat")
  expect_true("outlier_score" %in% colnames(result@meta.data))
  expect_true("outlier"       %in% colnames(result@meta.data))
})

test_that("detect_outliers outlier column is logical", {
  seu <- make_test_seurat(n_cells = 100)
  result <- suppressMessages(detect_outliers(seu))
  expect_type(result$outlier, "logical")
})

test_that("detect_outliers outlier_score values are between 0 and 1", {
  seu <- make_test_seurat(n_cells = 100)
  result <- suppressMessages(detect_outliers(seu))
  expect_true(all(result$outlier_score >= 0 & result$outlier_score <= 1))
})

test_that("detect_outliers flags fewer cells with a higher score_threshold", {
  seu <- make_test_seurat(n_cells = 200, seed = 7)
  strict <- suppressMessages(detect_outliers(seu, score_threshold = 0.95))
  lenient <- suppressMessages(detect_outliers(seu, score_threshold = 0.50))
  expect_lte(sum(strict$outlier), sum(lenient$outlier))
})

test_that("detect_outliers result has same number of cells as input", {
  seu <- make_test_seurat(n_cells = 80)
  result <- suppressMessages(detect_outliers(seu))
  expect_equal(ncol(result), 80)
})

test_that("remove_outliers_manual removes cells below negative threshold", {
  seu <- make_test_seurat(n_cells = 100)
  # Force CD1 values for first 10 cells to be very low so they are removed
  mat <- as.matrix(SeuratObject::GetAssayData(seu, layer = "data"))
  mat["CD1", 1:10] <- -10
  seu[["fcs"]] <- SeuratObject::SetAssayData(seu[["fcs"]], layer = "data",
                                             new.data = mat)
  result <- remove_outliers_manual(seu, thresholds = c(CD1 = 0), negative = TRUE)
  expect_lt(ncol(result), ncol(seu))
  expect_equal(ncol(result), 90)
})
