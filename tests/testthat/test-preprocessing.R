test_that("create_seurat returns a Seurat object", {
  fs    <- make_test_flowset()
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  expect_s4_class(seu, "Seurat")
})

test_that("create_seurat has correct total cell count", {
  fs    <- make_test_flowset(n_cells_per_sample = 50, n_samples = 2)
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  expect_equal(ncol(seu), 100)  # 50 cells × 2 samples
})

test_that("create_seurat feature names match panel$antigen", {
  fs    <- make_test_flowset()
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  expect_equal(sort(rownames(seu)), sort(panel$antigen))
})

test_that("create_seurat stores the panel data.frame directly in @misc", {
  fs    <- make_test_flowset()
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  misc_df <- as.data.frame(seu@misc)
  expect_equal(misc_df[, c("fcs_colname", "antigen")],
               panel[, c("fcs_colname", "antigen")],
               ignore_attr = TRUE)
})

test_that("create_seurat replaces underscores with dashes in antigen names", {
  fs    <- make_test_flowset(n_features = 2)
  panel <- data.frame(fcs_colname = c("CD1", "CD2"),
                      antigen     = c("CD_1", "CD_2"),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  expect_true("CD-1" %in% rownames(seu))
  expect_true("CD-2" %in% rownames(seu))
  expect_false("CD_1" %in% rownames(seu))
})

test_that("create_seurat joins sample metadata into meta.data", {
  fs    <- make_test_flowset(n_samples = 2)
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  meta  <- data.frame(sample_id = c("S1", "S2"),
                      condition  = c("ctrl", "treat"),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel, metadata = meta)
  expect_true("condition" %in% colnames(seu@meta.data))
  expect_equal(
    unique(seu@meta.data[seu@meta.data$sample_id == "S1", "condition"]),
    "ctrl"
  )
})

test_that("create_seurat puts channels absent from panel into 'unused' assay", {
  fs    <- make_test_flowset(n_features = 6)
  # panel only references 4 of the 6 channels; CD5 and CD6 should go to unused
  panel <- data.frame(fcs_colname = paste0("CD", 1:4),
                      antigen     = paste0("CD", 1:4),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel)
  expect_true("unused" %in% names(seu@assays))
  expect_equal(nrow(seu[["unused"]]), 2)
})

test_that("create_seurat with derandomize = TRUE produces ceiling-rounded counts", {
  fs    <- make_test_flowset()
  panel <- data.frame(fcs_colname = paste0("CD", 1:6),
                      antigen     = paste0("CD", 1:6),
                      stringsAsFactors = FALSE)
  seu <- create_seurat(fs, panel, derandomize = TRUE)
  counts <- as.matrix(SeuratObject::GetAssayData(seu, layer = "counts"))
  expect_true(all(counts == ceiling(counts)))
})
