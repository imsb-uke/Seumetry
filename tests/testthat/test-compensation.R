# Fixture pairing a flowSet with a matching Seurat object (same cell names
# and channels), built the same way create_seurat() expects.
make_comp_fixture <- function() {
  fs <- make_test_flowset()
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- suppressWarnings(create_seurat(fs, panel))
  comp_matrix <- diag(6)
  colnames(comp_matrix) <- rownames(comp_matrix) <- paste0("CD", 1:6)
  list(fs = fs, seu = seu, comp_matrix = comp_matrix)
}

test_that("compensate_data returns Seurat object with comp assay as default", {
  fx <- make_comp_fixture()
  seu2 <- suppressWarnings(compensate_data(fx$fs, fx$seu, comp_matrix = fx$comp_matrix))
  expect_s4_class(seu2, "Seurat")
  expect_true("comp" %in% SeuratObject::Assays(seu2))
  expect_equal(SeuratObject::DefaultAssay(seu2), "comp")
})

test_that("compensate_data preserves cell and feature dimensions", {
  fx <- make_comp_fixture()
  seu2 <- suppressWarnings(compensate_data(fx$fs, fx$seu, comp_matrix = fx$comp_matrix))
  expect_equal(ncol(seu2), ncol(fx$seu))
  expect_equal(nrow(seu2[["comp"]]), nrow(fx$seu))
})

test_that("compensate_data errors when fcs_fs is not a flowSet", {
  fx <- make_comp_fixture()
  expect_error(compensate_data("not_a_flowset", fx$seu, comp_matrix = fx$comp_matrix),
               "fcs_fs must be a flowSet")
})

test_that("compensate_data errors when seu is not a Seurat object", {
  fx <- make_comp_fixture()
  expect_error(compensate_data(fx$fs, "not_a_seurat", comp_matrix = fx$comp_matrix),
               "seu must be a Seurat object")
})

test_that("compensate_data errors when comp_matrix channels do not match panel", {
  fx <- make_comp_fixture()
  bad_matrix <- fx$comp_matrix
  colnames(bad_matrix) <- rownames(bad_matrix) <- paste0("XX", 1:6)
  expect_error(compensate_data(fx$fs, fx$seu, comp_matrix = bad_matrix),
               "None of the colnames of comp_matrix match")
})

test_that("compensate_data errors when comp_matrix is not a matrix or data.frame", {
  fx <- make_comp_fixture()
  expect_error(compensate_data(fx$fs, fx$seu, comp_matrix = list(a = 1)),
               "comp_matrix must be a matrix or data.frame")
})
