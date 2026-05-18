test_that("convert_seurat to FF returns a flowFrame", {
  seu <- make_test_seurat()
  ff  <- convert_seurat(seu, to = "FF")
  expect_s4_class(ff, "flowFrame")
})

test_that("convert_seurat to FF has correct dimensions", {
  seu <- make_test_seurat(n_cells = 80, n_features = 6)
  ff  <- convert_seurat(seu, to = "FF")
  expect_equal(nrow(ff@exprs), 80)  # cells
  expect_equal(ncol(ff@exprs), 6)   # features
})

test_that("convert_seurat to FF preserves expression values", {
  seu <- make_test_seurat(n_cells = 40, n_features = 4)
  ff  <- convert_seurat(seu, to = "FF")
  original <- as.matrix(SeuratObject::GetAssayData(seu, layer = "data"))
  # flowFrame is cells × features; Seurat is features × cells
  expect_equal(t(ff@exprs), original, ignore_attr = TRUE, tolerance = 1e-5)
})

test_that("convert_seurat to SCE returns a SingleCellExperiment", {
  seu <- make_test_seurat()
  sce <- convert_seurat(seu, to = "SCE")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("convert_seurat to SCE carries over metadata as colData", {
  seu <- make_test_seurat(n_cells = 60)
  sce <- convert_seurat(seu, to = "SCE")
  expect_true("sample_id" %in% colnames(SummarizedExperiment::colData(sce)))
})

test_that("convert_seurat to FS returns a flowSet split by the given column", {
  seu <- make_test_seurat(n_cells = 60)
  fs  <- convert_seurat(seu, to = "FS", split_by = "sample_id")
  expect_s4_class(fs, "flowSet")
  expect_equal(length(fs), 2)  # S1 and S2
})

test_that("convert_seurat to FS errors when split_by is missing", {
  seu <- make_test_seurat()
  expect_error(convert_seurat(seu, to = "FS"),
               "split_by parameter must be provided")
})
