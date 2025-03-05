#' @import ggplot2
#' @import Seurat
#' @importFrom flowCore read.flowSet exprs fsApply spillover fsApply compensate write.FCS compensation flowFrame flowSet identifier<-
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom FlowSOM ReadInput BuildSOM metaClustering_consensus GetClusters GetMetaclusters
#' @importFrom CATALYST normCytof
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @importFrom isotree isolation.forest
NULL