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

# --- run_modified_ransac ---------------------------------------------------

make_ransac_data <- function(n = 200, slope_ok = TRUE, seed = 1) {
  set.seed(seed)
  base <- abs(rnorm(n, mean = 4, sd = 1))
  if (slope_ok) {
    dat <- data.frame(x = base + rnorm(n, sd = 0.05), y = base + rnorm(n, sd = 0.05))
  } else {
    # unrelated x/y: no linear model should satisfy the slope constraint
    dat <- data.frame(x = abs(rnorm(n, mean = 4)), y = abs(rnorm(n, mean = 4)))
  }
  row.names(dat) <- paste0("cell_", seq_len(n))
  dat
}

test_that("run_modified_ransac fits a model with slope near 1 for x=y-like data", {
  dat <- make_ransac_data()
  model <- run_modified_ransac(formula = y ~ x, data = dat)
  expect_s3_class(model, "lm")
  expect_equal(unname(coef(model)["x"]), 1, tolerance = 0.1)
})

test_that("run_modified_ransac returns NULL with a warning when too few events", {
  dat <- make_ransac_data(n = 5)
  expect_warning(result <- run_modified_ransac(formula = y ~ x, data = dat, min_to_fit = 10),
                 "Not enough events")
  expect_null(result)
})

test_that("run_modified_ransac errors when formula is not a formula object", {
  dat <- make_ransac_data()
  expect_error(run_modified_ransac(formula = "y ~ x", data = dat),
               "formula must be a formula object")
})

test_that("run_modified_ransac errors when formula variables are missing from data", {
  dat <- make_ransac_data()
  expect_error(run_modified_ransac(formula = z ~ x, data = dat),
               "not found in data: z")
})

test_that("run_modified_ransac errors on invalid min_to_fit", {
  dat <- make_ransac_data()
  expect_error(run_modified_ransac(formula = y ~ x, data = dat, min_to_fit = 1),
               "min_to_fit must be")
})

# --- detect_aggregates ------------------------------------------------------

test_that("detect_aggregates flags correlated channels as aggregates", {
  seu <- make_correlated_seurat()
  channels <- detect_aggregate_channels(seu, threshold = 0.7)
  df <- channels[[3]]
  seu2 <- suppressMessages(detect_aggregates(seu, df, score_threshold = 1))
  expect_true("aggregate" %in% colnames(seu2@meta.data))
  expect_true("aggregate_score" %in% colnames(seu2@meta.data))
  expect_true(any(seu2$aggregate))
})

test_that("detect_aggregates respects score_threshold", {
  seu <- make_correlated_seurat()
  channels <- detect_aggregate_channels(seu, threshold = 0.7)
  df <- channels[[3]]
  # only one channel pair -> max possible score is 1, so score_threshold = 2 flags nothing
  seu2 <- suppressMessages(detect_aggregates(seu, df, score_threshold = 2))
  expect_false(any(seu2$aggregate))
})

test_that("detect_aggregates errors when seu is not a Seurat object", {
  seu <- make_correlated_seurat()
  channels <- detect_aggregate_channels(seu, threshold = 0.7)
  expect_error(detect_aggregates("not_a_seurat", channels[[3]]),
               "seu must be a Seurat object")
})

test_that("detect_aggregates errors when aggregate_channels is not a valid data.frame", {
  seu <- make_correlated_seurat()
  expect_error(detect_aggregates(seu, list(a = 1)),
               "aggregate_channels must be a data.frame")
})

test_that("detect_aggregates errors when aggregate_channels has no rows", {
  seu <- make_correlated_seurat()
  channels <- detect_aggregate_channels(seu, threshold = 0.7)
  empty_df <- channels[[3]][0, ]
  expect_error(detect_aggregates(seu, empty_df),
               "aggregate_channels has no rows")
})

test_that("detect_aggregates errors when channels are not found in the Seurat object", {
  seu <- make_correlated_seurat()
  channels <- detect_aggregate_channels(seu, threshold = 0.7)
  bad_df <- channels[[3]]
  bad_df$Channel_1 <- "NOPE"
  expect_error(detect_aggregates(seu, bad_df),
               "not found in the Seurat object")
})
