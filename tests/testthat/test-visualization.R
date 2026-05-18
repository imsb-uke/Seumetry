test_that("plot_pca returns a ggplot object", {
  seu <- make_test_seurat(n_cells = 60)
  p <- plot_pca(seu, group_by = "sample_id", color = "sample_id")
  expect_s3_class(p, "ggplot")
})

test_that("plot_frequency returns a ggplot when return_table = FALSE", {
  seu <- make_test_seurat(n_cells = 80)
  seu@meta.data$cluster <- rep(c("C1", "C2", "C3", "C4"), 20)
  p <- plot_frequency(seu, attribute_1 = "cluster", attribute_2 = "sample_id")
  expect_s3_class(p, "ggplot")
})

test_that("plot_frequency returns a data.frame when return_table = TRUE", {
  seu <- make_test_seurat(n_cells = 80)
  seu@meta.data$cluster <- rep(c("C1", "C2", "C3", "C4"), 20)
  tbl <- plot_frequency(seu, attribute_1 = "cluster", attribute_2 = "sample_id",
                        return_table = TRUE)
  expect_s3_class(tbl, "data.frame")
  expect_true("freq" %in% colnames(tbl))
  expect_true("perc" %in% colnames(tbl))
})

test_that("plot_frequency percentages sum to 100 per attribute_2 group", {
  seu <- make_test_seurat(n_cells = 80)
  seu@meta.data$cluster <- rep(c("C1", "C2", "C3", "C4"), 20)
  tbl <- plot_frequency(seu, attribute_1 = "cluster", attribute_2 = "sample_id",
                        return_table = TRUE)
  totals <- tapply(tbl$perc, tbl$attribute_2, sum)
  expect_true(all(abs(totals - 100) < 1e-6))
})

test_that("plot_cellnumber returns a ggplot when return_table = FALSE", {
  seu <- make_test_seurat()
  p <- plot_cellnumber(seu, group_by = "sample_id")
  expect_s3_class(p, "ggplot")
})

test_that("plot_cellnumber returns a data.frame when return_table = TRUE", {
  seu <- make_test_seurat(n_cells = 100)
  tbl <- plot_cellnumber(seu, group_by = "sample_id", return_table = TRUE)
  expect_s3_class(tbl, "data.frame")
  expect_true("Freq" %in% colnames(tbl))
  expect_equal(sum(tbl$Freq), 100)
})

test_that("plot_cyto returns a ggplot for 2d_density style", {
  seu <- make_test_seurat(n_cells = 50)
  p <- plot_cyto(seu, x = "CD1", y = "CD2", style = "2d_density")
  expect_s3_class(p, "ggplot")
})

test_that("plot_cyto returns a ggplot for density style", {
  seu <- make_test_seurat(n_cells = 50)
  p <- plot_cyto(seu, x = "CD1", style = "density")
  expect_s3_class(p, "ggplot")
})

test_that("plot_cyto errors on invalid style", {
  seu <- make_test_seurat(n_cells = 50)
  expect_error(plot_cyto(seu, x = "CD1", style = "invalid"), "valid style")
})
