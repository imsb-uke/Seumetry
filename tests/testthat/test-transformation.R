test_that("transform_arcsinh applies asinh(x/cofactor) correctly", {
  mat <- matrix(c(0, 5, 50, 500), nrow = 1,
                dimnames = list("X", paste0("c", 1:4)))
  panel <- data.frame(fcs_colname = "X", antigen = "X",
                      arcsinh_cofactor = 5, stringsAsFactors = FALSE)
  result <- transform_arcsinh(mat, panel)
  expect_equal(result[1, 1], asinh(0 / 5))
  expect_equal(result[1, 2], asinh(5 / 5))
  expect_equal(result[1, 3], asinh(50 / 5))
  expect_equal(result[1, 4], asinh(500 / 5))
})

test_that("transform_arcsinh uses default cofactor of 5 when column absent", {
  mat <- matrix(5, nrow = 1, dimnames = list("X", "c1"))
  panel <- data.frame(fcs_colname = "X", antigen = "X",
                      stringsAsFactors = FALSE)
  result <- transform_arcsinh(mat, panel)
  # result[1,1] is a named numeric; ignore_attr drops the name for comparison
  expect_equal(result[1, 1], asinh(1), ignore_attr = TRUE)   # 5/5 = 1
})

test_that("transform_arcsinh handles per-channel cofactors independently", {
  # Must use ≥2 columns: sapply over 2 rows returns a 2×n matrix; t() then
  # gives n×2. With only 1 column, t(sapply(...)) produces a 1×2 result but
  # colnames<- tries to assign 1 name to 2 columns → dimension error.
  mat <- matrix(10, nrow = 2, ncol = 2,
                dimnames = list(c("A", "B"), c("c1", "c2")))
  panel <- data.frame(fcs_colname = c("A", "B"), antigen = c("A", "B"),
                      arcsinh_cofactor = c(5, 10), stringsAsFactors = FALSE)
  result <- transform_arcsinh(mat, panel)
  expect_equal(result["A", "c1"], asinh(10 / 5),  ignore_attr = TRUE)
  expect_equal(result["B", "c1"], asinh(10 / 10), ignore_attr = TRUE)
})

test_that("transform_data returns Seurat with transformed data layer", {
  seu <- make_test_seurat()
  raw <- as.matrix(SeuratObject::GetAssayData(seu, layer = "counts"))
  result <- transform_data(seu, transformation = "arcsinh")
  expect_s4_class(result, "Seurat")
  tfm <- as.matrix(SeuratObject::GetAssayData(result, layer = "data"))
  expect_false(identical(raw, tfm))
  # arcsinh output must be non-negative for non-negative input
  expect_true(all(tfm >= 0))
})

test_that("transform_data does not modify the counts layer", {
  seu <- make_test_seurat()
  before <- as.matrix(SeuratObject::GetAssayData(seu, layer = "counts"))
  result <- transform_data(seu, transformation = "arcsinh")
  after  <- as.matrix(SeuratObject::GetAssayData(result, layer = "counts"))
  expect_equal(before, after)
})
