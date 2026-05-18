test_that("median_expression returns matrix with correct dimensions", {
  seu <- make_test_seurat(n_cells = 100, n_features = 8)
  result <- median_expression(seu, group_by = "sample_id")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 8)   # features
  expect_equal(ncol(result), 2)   # 2 groups (S1, S2)
  expect_setequal(colnames(result), c("S1", "S2"))
})

test_that("median_expression computes correct medians", {
  set.seed(42)
  n_cells <- 40
  mat <- matrix(1, nrow = 2, ncol = n_cells,
                dimnames = list(c("CD1", "CD2"), paste0("cell_", seq_len(n_cells))))
  mat["CD1", 1:20] <- 3   # S1 group
  mat["CD1", 21:40] <- 7  # S2 group
  panel <- data.frame(fcs_colname = c("CD1", "CD2"),
                      antigen     = c("CD1", "CD2"),
                      stringsAsFactors = FALSE)
  seu <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = mat, assay = "fcs",
                                                           min.cells = 0, min.features = 0))
  seu[["fcs"]] <- SeuratObject::SetAssayData(seu[["fcs"]], layer = "data",
                                             new.data = mat)
  seu@misc <- panel
  seu@meta.data$sample_id <- rep(c("S1", "S2"), each = 20)
  result <- median_expression(seu, group_by = "sample_id")
  expect_equal(result["CD1", "S1"], 3)
  expect_equal(result["CD1", "S2"], 7)
})

test_that("median_expression column names match the group_by values", {
  seu <- make_test_seurat(n_cells = 60, n_features = 5)
  seu@meta.data$condition <- rep(c("ctrl", "treat", "ctrl"), 20)
  result <- median_expression(seu, group_by = "condition")
  expect_setequal(colnames(result), c("ctrl", "treat"))
})
