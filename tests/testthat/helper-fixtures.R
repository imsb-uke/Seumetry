# Shared synthetic data fixtures — auto-loaded by testthat before all test files.

# Standard Seurat fixture: n_features rows (markers) × n_cells cols
make_test_seurat <- function(n_cells = 100, n_features = 10, seed = 1) {
  set.seed(seed)
  mat <- matrix(
    abs(rnorm(n_features * n_cells, mean = 3, sd = 1)),
    nrow     = n_features,
    ncol     = n_cells,
    dimnames = list(paste0("CD", seq_len(n_features)),
                    paste0("cell_", seq_len(n_cells)))
  )
  panel <- data.frame(
    fcs_colname = paste0("CD", seq_len(n_features)),
    antigen     = paste0("CD", seq_len(n_features)),
    stringsAsFactors = FALSE
  )
  # suppressWarnings: Seurat v5 emits a "Coercing to dgCMatrix" note for dense input
  seu <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = mat, assay = "fcs",
                                     min.cells = 0, min.features = 0)
  )
  seu[["fcs"]] <- SeuratObject::SetAssayData(seu[["fcs"]], layer = "data",
                                             new.data = mat)
  # create_seurat stores the panel data.frame directly in @misc
  seu@misc <- panel
  seu@meta.data$sample_id <- rep(c("S1", "S2"), each = n_cells / 2)
  seu
}

# In-memory FlowSet fixture (no FCS files needed) — for create_seurat tests.
# Mirrors what create_flowset does: sets cell IDs directly on @exprs so that
# create_seurat's sub("_\\d+$", "", cell_id) parse returns the right sample_id.
make_test_flowset <- function(n_cells_per_sample = 50, n_features = 6,
                              n_samples = 2, seed = 1) {
  set.seed(seed)
  sample_names <- paste0("S", seq_len(n_samples))
  feature_names <- paste0("CD", seq_len(n_features))
  make_ff <- function(sid) {
    mat <- matrix(
      abs(rnorm(n_cells_per_sample * n_features, mean = 5)),
      nrow = n_cells_per_sample, ncol = n_features,
      dimnames = list(NULL, feature_names)
    )
    ff <- flowCore::flowFrame(exprs = mat)
    # Set row names via the proper accessor — flowFrame() does not preserve
    # matrix row names, so we must set them explicitly (same as create_flowset)
    e <- flowCore::exprs(ff)
    rownames(e) <- paste0(sid, "_", seq(0, n_cells_per_sample - 1))
    flowCore::exprs(ff) <- e
    ff
  }
  ff_list <- stats::setNames(lapply(sample_names, make_ff), sample_names)
  flowCore::flowSet(ff_list)
}

# Seurat fixture with a pair of strongly correlated channels (for aggregate tests)
make_correlated_seurat <- function(n_cells = 300, seed = 1) {
  set.seed(seed)
  base <- abs(rnorm(n_cells, mean = 4, sd = 1))
  mat <- rbind(
    CD1 = base + rnorm(n_cells, sd = 0.05),   # very correlated with CD2
    CD2 = base + rnorm(n_cells, sd = 0.05),
    CD3 = abs(rnorm(n_cells, mean = 3)),
    CD4 = abs(rnorm(n_cells, mean = 3))
  )
  colnames(mat) <- paste0("cell_", seq_len(n_cells))
  panel <- data.frame(
    fcs_colname = rownames(mat),
    antigen     = rownames(mat),
    stringsAsFactors = FALSE
  )
  seu <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = mat, assay = "fcs",
                                     min.cells = 0, min.features = 0)
  )
  seu[["fcs"]] <- SeuratObject::SetAssayData(seu[["fcs"]], layer = "data",
                                             new.data = mat)
  seu@misc <- panel
  seu@meta.data$sample_id <- rep(c("S1", "S2"), each = n_cells / 2)
  seu
}

# 4-sample Seurat (2 ctrl, 2 treat) for differential abundance tests.
# Needs ≥3 samples so edgeR's estimateDisp has residual degrees of freedom.
make_da_seurat <- function(n_cells = 200, n_features = 6, seed = 1) {
  set.seed(seed)
  mat <- matrix(
    abs(rnorm(n_features * n_cells, mean = 3, sd = 1)),
    nrow = n_features, ncol = n_cells,
    dimnames = list(paste0("CD", seq_len(n_features)),
                    paste0("cell_", seq_len(n_cells)))
  )
  panel <- data.frame(fcs_colname = paste0("CD", seq_len(n_features)),
                      antigen     = paste0("CD", seq_len(n_features)),
                      stringsAsFactors = FALSE)
  seu <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = mat, assay = "fcs",
                                     min.cells = 0, min.features = 0)
  )
  seu[["fcs"]] <- SeuratObject::SetAssayData(seu[["fcs"]], layer = "data",
                                             new.data = mat)
  seu@misc <- panel
  seu@meta.data$sample_id        <- rep(paste0("S", 1:4), each = n_cells / 4)
  seu@meta.data$condition         <- rep(c("ctrl", "treat"), each = n_cells / 2)
  seu@meta.data$seurat_clusters   <- factor(rep(1:4, n_cells / 4))
  seu
}
