#' Plot PCA plot based on attribute level, not single cells, using prcomp
#'
#' First calculates median fluorescence intensity for each sample and then runs a PCA using prcomp. Returns ggplot2 object.
#'
#' @param seu Seurat object.
#' @param group_by Indicate group which is used for pseudobulking using [median_expression()]. Default: sample_id.
#' @param color Indicate attribute used to color points. Default: sample_id.
#' @param shape Add a second level of visualization using shape.
#' @param label_points Indicate if points should be labelled.
#' @return ggplot2 plot.
#' @export
#' @examples
#' plot = plot_sample_pca(seu, color = "attribute1", shape = "attribute2", label_points = TRUE)
#' plot + ggtitle("This is a sample-level PCA plot")
plot_pca <- function(seu,
                     group_by = "sample_id",
                     color = "sample_id",
                     shape = NULL,
                     label_points = FALSE) {
  # get average expression
  ave_df <- median_expression(seu, group_by = group_by)
  # perform a PCA on the data
  pca <- prcomp(t(ave_df))
  # the contribution to the total variance for each component
  percentVar <- round(pca$sdev^2/sum(pca$sdev^2)*100)
  # create metadata data.frame based on seu@meta.data
  metadata <- seu@meta.data[match(rownames(pca$x), seu@meta.data[, group_by]),]
  # add meta data
  ggdata <- cbind(pca$x, metadata)
  # make list of additional aesthetics
  aes_list <- list()
  if (!is.null(color)) aes_list$color <- ggdata[[color]]
  if (!is.null(shape)) aes_list$shape <- ggdata[[shape]]
  # plot PCA using ggplot2
  plot <- ggplot(ggdata) +
          geom_point(size = 3, aes(x = PC1, y = PC2, !!!aes_list)) + 
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(color = color, shape = shape) +
          theme_linedraw() + theme(aspect.ratio=1,
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank())
  # add labels
  if(isTRUE(label_points)) plot <- plot + geom_text(aes(label = rownames(ggdata)))
  return(plot)
}
                     
                     
#' Stacked barplot with percentage or number of cells per attribute by attribute
#'
#' @param seu Seurat object.
#' @param attribute_1 Attribute used for y axis.
#' @param attribute_2 Attribute used for x axis.
#' @param n_cells Plot number of cells instead of frequencies. Default: FALSE.
#' @param return_table Return table with frequencies instead of ggplot2 plot. Default: FALSE.
#' @return ggplot2 plot, or table if return_table = TRUE
#' @export
#' @examples
#' plot_frequency(seu, "seurat_clusters", "sample_id")
plot_frequency <- function(seu,
                           attribute_1,
                           attribute_2,
                           n_cells = FALSE,
                           return_table = FALSE) {
  # create data.frame with cellcounts per attributes
  datatable <- data.frame(table(seu@meta.data[, attribute_1], seu@meta.data[, attribute_2]))
  # create ggplot dataframe
  names(datatable) <- c("attribute_1", "attribute_2", "freq")
  # calculate percentage of cells per clusters for the provided attributes
  for(i in 1:nrow(datatable)) {
    sum <- sum(datatable[datatable[,"attribute_2"] == datatable[i,"attribute_2"], "freq"])
    datatable[i,"perc"] <- (datatable[i,"freq"] / sum) * 100
  }
  # plot graph with frequencies or ncells
  if(isTRUE(n_cells))
      plot <- ggplot(datatable, aes(x = attribute_2, y = freq, fill = attribute_1)) + ylab("Number of cells")
  else
      plot <- ggplot(datatable, aes(x = attribute_2, y = perc, fill = attribute_1)) + ylab("Frequency (%)")
  # add other plotting visualisations to plot
  plot <- plot + geom_bar(stat = "identity", position = "stack") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = c(0, 0)) + xlab("") + labs(fill = attribute_1)
  # return table or ggplot
  if(isTRUE(return_table)) {
    return(datatable)
  }else{
    return(plot)
  }
}
                     
                     
#' Barplot cell numbers
#'
#' @param seu Seurat object
#' @param group Group to plot number of cells. Default: sample_id.
#' @param return_table Return table with frequencies instead of ggplot2 plot. Default: FALSE.
#' @return ggplot2 plot, or table if return_table = TRUE
#' @export
#' @examples
#' plot_cellnumber(seu)
plot_cellnumber <- function(seu,
                            group_by = "sample_id",
                            return_table = FALSE) {
  # retrieve cell number data.frame
  ncells_df <- data.frame(table(seu[[group_by]]))
  # compite stats
  max <- max(ncells_df$Freq)
  min <- min(ncells_df$Freq)
  # make ggplot
  plot <- ggplot(ncells_df, aes(x = .data[[group_by]], y = Freq)) +
    geom_bar(stat = "summary", width=0.75) +
    xlab(group_by) + ylab("Number of cells") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(caption = paste("Max:", max(ncells_df$Freq),
                         "// Min:", min(ncells_df$Freq),
                         "// Median:", median(ncells_df$Freq)))
  # return table or ggplot
  if(isTRUE(return_table)) {
    return(ncells_df)
  }else{
    return(plot)
  }
}


#' Plot Cytometry-Style plots on Seurat objects
#'
#' @param seu Seurat object.
#' @param x Which marker to plot on x axis.
#' @param y Which marker to plot on y axis (only required for style = 2d_density or points)
#' @param style Choose from: "2d_density", "point", or "density". Default: 2d_density.
#' @param slot Choose which Seurat slot to use. Default: data (transformed).
#' @param assay Choose which Seurat assay to use. Default: DefaultAssay(seu).
#' @param scale Choose which scale to use for x and y axis (linear, log, biexp). Default: linear.
#' @param limits List of x and y limits. Example: list(x = c(0,5), y = c(0,5)).
#' @param color Color cells based on meta.data column in Seurat object, e.g. "sample_id".
#' @param pt_size When using style = "point", this parameter adjusts the point size (Default: 0.2).
#' @param rasterize Rasterize plot. Useful for plotting large datasets (Default: FALSE).
#' @param alpha When using style = "point", this parameter adjusts the alpha (Default: 1).
#' @param bins When using style = "2d_density", this parameter adjusts the number of bins.
#' @param biexp_pos If scale = "biexp", indicate biexp_pos here (see [transform_biexp()]).
#' @param biexp_neg If scale = "biexp", indicate biexp_neg here (see [transform_biexp()]).
#' @param biexp_widthBasis If scale = "biexp", indicate biexp_widthBasis here (see [transform_biexp()]).
#' @return ggplot2 plot.
#' @export
#' @examples
#' plot_cyto(seu, x = "CD3", y = "CD4")
#' plot_cyto(seu, x = "CD3", y = "CD4", style = "point", color = "sample_id")
#' plot_cyto(seu, x = "CD3", style = "density", color = "sample_id")
plot_cyto <- function(seu,
                      x,
                      y,
                      style = "2d_density",
                      slot = "data",
                      assay = NULL,
                      scale = "linear",
                      limits = NULL,
                      color = NULL,
                      pt_size = 0.1,
                      rasterize = FALSE,
                      alpha = 1,
                      bins = 200,
                      biexp_pos = 4.5,
                      biexp_neg = 0,
                      biexp_widthBasis = -10) {
  # check input
  if(!style %in% c("2d_density", "density", "point")) stop("Indicate a valid style.", call. = FALSE)
  # get dataframe with intensities
  gg_df = as.data.frame(t(as.matrix(GetAssayData(seu, assay = assay, slot = slot))))
  # add color to gg data.frame and aesthetics
  aes_list = list()
  if(!is.null(color)) {
    if(color %in% colnames(seu@meta.data)) {
        # add color to gg data.frame
        gg_df[[color]] = seu@meta.data[[color]]
        # shuffel so that overlapping points are plotted at random
        gg_df = gg_df[sample(1:nrow(gg_df)), ]
        # add color to the aes list
        aes_list = list(color = gg_df[[color]])
    }else{
        stop("Color input is not part of Seurat metadata.", call. = FALSE)
    }
  }
  # plot data
  plot = ggplot(gg_df, aes(x = .data[[x]], !!!aes_list)) +
    theme_linedraw() + theme(aspect.ratio=1,
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())
  # add style of ggplot
  if(style == "2d_density") plot = plot + geom_bin2d(aes(y = .data[[y]]), bins = bins) +
    viridis::scale_fill_viridis()
  if(style == "density") plot = plot + geom_density()
  if(style == "point") plot = plot + geom_point(aes(y = .data[[y]]), size = pt_size, alpha = alpha)
  # change scale if provided
  if(scale == "log") plot = plot + scale_y_log10() + scale_x_log10()
  if(scale == "biexp") {
    plot = plot +
      ggcyto::scale_x_flowjo_biexp(widthBasis=biexp_widthBasis, pos=biexp_pos, neg=biexp_neg) +
      ggcyto::scale_y_flowjo_biexp(widthBasis=biexp_widthBasis, pos=biexp_pos, neg=biexp_neg)
  }
  # add limits if provided
  if(!is.null(limits$x)) plot = plot + xlim(limits$x)
  if(!is.null(limits$y)) plot = plot + ylim(limits$y)
  # add legend title if color was changed
  if("color" %in% names(aes_list)) plot = plot + labs(color = color) 
  # rasterize
  if(rasterize) plot <- ggrastr::rasterize(plot, dpi = 600)
  #return 
  return(plot)
}
