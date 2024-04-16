#' Load multiple FCS files as flowSet Object
#'
#' Uses [flowCore::read.flowSet()] to generate a flowSet based on files within a certain directory.
#' In addition, this function will give each cell a unique name based on the sample_id so it can be identified later in downstream analyses.
#' This flowSet Object should be stored in order to export FCS files based on unique cell IDs later.
#'
#' @param fcs_dir Directory containing all FCS files.
#' @param metadata Data.frame with metadata for each FCS file (rows). Must contain columns "sample_id" and "file_name".
#' @return A flowCore FlowSet object.
#' @export
#' @examples
#' create_flowset("fcs_files_directory", metadata)
create_flowset <- function(fcs_dir,
                           metadata) {
  # check if all files in folder are present in metadata
  if(length(setdiff(list.files(fcs_dir), metadata$file_name)) != 0) {
    stop("File names provided in metadata do not match files in fcs folder.")
  }
  # use flowCore and create flowSet
  fcs_fs <- read.flowSet(
    path = fcs_dir,
    transformation = FALSE, # Do not transform data. Transformation happens later.
    truncate_max_range = FALSE, # Inhibits removal of very bright events!
    alter.names = TRUE # Change names to R compatible format (remove dashes etc)
  )
  # name cells (sample_id_n)
  for(x in 1:length(fcs_fs@frames)) {
    # get name of FCS file (flowFrame); corresponds to metadata$file_name
    filename <- fcs_fs[[x]]@description[["GUID"]]
    # get sample_id of FCS file according to metadata$sample_id column
    sample_id <- as.character(metadata[which(metadata$file_name == filename), "sample_id"])
    # get number of cells in flowFrame
    n_cells <- nrow(fcs_fs[[x]]@exprs)
    # make vector of unique cell names
    unique_names <- make.unique(rep(sample_id, n_cells), sep = "_")
    unique_names[1] <- paste0(sample_id, "_0")
    # add cellnames to flowFrame
    row.names(fcs_fs[[x]]@exprs) <- unique_names
  }
  # check if all files loaded successfully
  if(length(setdiff(metadata$file_name, names(fcs_fs@frames))) != 0) {
    stop("Some FCS files were not loaded properly.")
  }else{
    message("All files loaded successfully.")
    return(fcs_fs)
  }
}


#' Converts a flowCore flowSet containing multiple flowFrames into a Seurat object
#'
#' Uses sample metadata to generate cell-level metadata. Each row = one FCS file. Must contain "file_name" and "sample_id" columns. Uses panel data.frame to subset a flowSet to channels that should be analyzed, and channels that will not be used for downstream analysis. It's important that panel$fcs_colnames contains all channels that should be used. For additional explanation of panel_df, see vignette.
#'
#' Channel names can be checked in the flowFrame using: colnames(fcs_fs@frames$fcs_id@exprs). Raw FCS data is stored in seu@assays$fcs. Raw FCS data of unused channels is stored in seu@assays$unused. Panel data is stored in seu@misc slot.
#'
#' @param fcs_fs A flowCore flowSet generated with create_flowset.
#' @param metadata Data.frame with metadata for each FCS file (rows). Must contain columns "sample_id" and "file_name" and each row must correspond to the flowFrames in the flowSet.
#' @param panel A panel containing all information for the markers. See vignette for more information on how to generate a panel_df. Function will only use channels indicated in panel$fcs_colnames and will rename channels to panel$antigen.
#' @param ... Additional parameters that will be passed on to Seurat::CreateSeuratObject.
#' @return A Seurat Object of all merged FCS files.
#' @export
#' @examples
#' seu = create_seurat(fcs_fs, panel, metadata)
create_seurat <- function(fcs_fs,
                          panel,
                          metadata = NULL,
                          ...) {
  # create matrix of all fcs files in flowSet
  matrix <- fsApply(fcs_fs, exprs)
  # check if all channels in panel can be found in FCS files
  if(!all(panel$fcs_colname %in% colnames(matrix))) {
      stop("Not all channels provided in panel$fcs_colname are present in FCS files! Make sure names of channels overlap (names(fcs_fs[[1]])).")
  }
  # create matrix of unused channels
  matrix_unused <- matrix[, !(colnames(matrix) %in% panel$fcs_colname)]
  # subset to only contain channels present in panel
  matrix <- matrix[, colnames(matrix) %in% panel$fcs_colname]
  # rename channels to antigen in panel
  colnames(matrix) <- panel[match(colnames(matrix), panel$fcs_colname), ]$antigen
  # transpose to have Seurat format
  matrix <- t(matrix) 
  # create cell level metadata
  if(!is.null(metadata)) {
    # get all cell_ids
    cell_ids <- colnames(matrix)
    # make dataframe with cell_ids as rownames
    cell_meta <- data.frame(row.names = cell_ids)
    # extract sample_id from cell_id
    sample_id <- sub("_\\d+$", "", cell_ids)
    # add meta to each cell
    meta_df <- cbind(cell_meta, metadata[match(sample_id,metadata$sample_id),])
  }else{
    meta_df <- NULL
  }
  # create seurat object
  seu <- CreateSeuratObject(counts = matrix,
                            assay = "fcs",
                            meta.data = meta_df,
                            min.cells = 0,
                            min.features = 0,
                            ...)
  # add panel data to seu@misc slot
  seu@misc <- panel
  # add unused channels to new assay: "unused"
  seu[["unused"]] <- CreateAssayObject(counts = t(matrix_unused),
                                       min.cells = 0,
                                       min.features = 0)
  # set idents to sample_id if metadata was provided
  if(!is.null(metadata)) Idents(seu) = seu$sample_id
  # return Seurat object
  return(seu)
}
