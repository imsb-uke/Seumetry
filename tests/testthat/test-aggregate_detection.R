test_that("detect_aggregate_channels errors when no correlated channels found", {
  seu <- make_test_seurat(n_cells = 200, n_features = 6, seed = 1)
  # Use a very high threshold so nothing passes
  expect_error(detect_aggregate_channels(seu, threshold = 0.999),
               "No channels containing aggregates found")
})

test_that("detect_aggregate_channels returns list of 3 elements for correlated data", {
  seu <- make_correlated_seurat()
  result <- detect_aggregate_channels(seu, threshold = 0.7)
  expect_type(result, "list")
  expect_length(result, 3)
})

test_that("detect_aggregate_channels [[1]] is a square correlation matrix", {
  seu <- make_correlated_seurat()
  result <- detect_aggregate_channels(seu, threshold = 0.7)
  mat <- result[[1]]
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), ncol(mat))
  expect_equal(rownames(mat), colnames(mat))
})

test_that("detect_aggregate_channels [[2]] contains only 0 and 1", {
  seu <- make_correlated_seurat()
  result <- detect_aggregate_channels(seu, threshold = 0.7)
  binary <- result[[2]]
  expect_true(all(binary %in% c(0, 1)))
})

test_that("detect_aggregate_channels [[3]] data.frame has expected columns", {
  seu <- make_correlated_seurat()
  result <- detect_aggregate_channels(seu, threshold = 0.7)
  df <- result[[3]]
  expect_s3_class(df, "data.frame")
  expect_true(all(c("Channel_1", "Channel_2", "Pearson_R") %in% colnames(df)))
})

test_that("detect_aggregate_channels detects CD1-CD2 correlation in correlated data", {
  seu <- make_correlated_seurat()
  result <- detect_aggregate_channels(seu, threshold = 0.7)
  df <- result[[3]]
  pairs <- paste(df$Channel_1, df$Channel_2, sep = "-")
  expect_true(any(pairs %in% c("CD1-CD2", "CD2-CD1")))
})
