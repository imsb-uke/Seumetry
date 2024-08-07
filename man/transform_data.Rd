% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation.R
\name{transform_data}
\alias{transform_data}
\title{Wrapper for all transformation functions using Seurat object as input & output

See transformation functions for detail: \code{\link[=transform_arcsinh]{transform_arcsinh()}}, \code{\link[=transform_biexp]{transform_biexp()}}. Input Seurat Object uses "counts" slot from current assay, performs transformation and writes transformed counts into "data" slot from current assay using panel information stored in seu@misc.}
\usage{
transform_data(seu, transformation)
}
\arguments{
\item{seu}{Seurat object.}

\item{transformation}{Indicate which transformation to use (arcsinh, biexp)}
}
\value{
Seurat object with transformed data in "data" slot of current assay.
}
\description{
Wrapper for all transformation functions using Seurat object as input & output

See transformation functions for detail: \code{\link[=transform_arcsinh]{transform_arcsinh()}}, \code{\link[=transform_biexp]{transform_biexp()}}. Input Seurat Object uses "counts" slot from current assay, performs transformation and writes transformed counts into "data" slot from current assay using panel information stored in seu@misc.
}
\examples{
seu <- transform_data(seu, "arcsinh")
}
