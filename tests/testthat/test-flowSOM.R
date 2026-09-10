test_that("run_FlowSOM adds SOM_cl and SOM_metacl to meta.data", {
  seu <- make_test_seurat(n_cells = 200, n_features = 6)
  seu2 <- suppressWarnings(suppressMessages(
    run_FlowSOM(seu, xdim = 5, ydim = 5, metaclusters = 3)
  ))
  expect_true(all(c("SOM_cl", "SOM_metacl") %in% colnames(seu2@meta.data)))
  expect_equal(ncol(seu2), ncol(seu))
  expect_true(is.factor(seu2$SOM_cl))
  expect_true(is.factor(seu2$SOM_metacl))
  # every cell must be assigned to one of the requested metaclusters
  expect_true(all(as.integer(as.character(seu2$SOM_metacl)) %in% seq_len(3)))
})

test_that("run_FlowSOM errors when seu is not a Seurat object", {
  expect_error(run_FlowSOM("not_a_seurat"), "seu must be a Seurat object")
})

test_that("run_FlowSOM errors when features are not present in the Seurat object", {
  seu <- make_test_seurat(n_cells = 50, n_features = 6)
  expect_error(run_FlowSOM(seu, features = c("CDX")),
               "not present in the Seurat object")
})

test_that("run_FlowSOM errors when assay does not exist", {
  seu <- make_test_seurat(n_cells = 50, n_features = 6)
  expect_error(run_FlowSOM(seu, assay = "nope"),
               "Assay 'nope' not found")
})

test_that("run_FlowSOM errors on non-positive xdim, ydim, or metaclusters", {
  seu <- make_test_seurat(n_cells = 50, n_features = 6)
  expect_error(run_FlowSOM(seu, xdim = 0), "xdim must be a single positive whole number")
  expect_error(run_FlowSOM(seu, ydim = -2), "ydim must be a single positive whole number")
  expect_error(run_FlowSOM(seu, metaclusters = -1), "metaclusters must be a single positive whole number")
})
