% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_pca}
\alias{plot_pca}
\title{Plot PCA plot based on attribute level, not single cells, using prcomp}
\usage{
plot_pca(
  seu,
  group_by = "sample_id",
  color = "sample_id",
  shape = NULL,
  label_points = FALSE
)
}
\arguments{
\item{seu}{Seurat object.}

\item{group_by}{Indicate group which is used for pseudobulking using \code{\link[=median_expression]{median_expression()}}. Default: sample_id.}

\item{color}{Indicate attribute used to color points. Default: sample_id.}

\item{shape}{Add a second level of visualization using shape.}

\item{label_points}{Indicate if points should be labelled.}
}
\value{
ggplot2 plot.
}
\description{
First calculates median fluorescence intensity for each sample and then runs a PCA using prcomp. Returns ggplot2 object.
}
\examples{
plot = plot_sample_pca(seu, color = "attribute1", shape = "attribute2", label_points = TRUE)
plot + ggtitle("This is a sample-level PCA plot")
}
