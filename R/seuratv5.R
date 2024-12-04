#' SeuratV5 convertion
#'
#' @param obj SingleCellExperiment objcet to convert into a Seurat V3
#'
#' @importFrom Seurat ScaleData CreateDimReducObject IntegrateLayers
#' @importFrom SeuratObject CreateAssay5Object JoinLayers CreateSeuratObject
#' @importFrom SingleCellExperiment counts logcounts colData reducedDim
#'
#' @export
#'
sce_to_seuratv5 <- function(obj){
  options(Seurat.object.assay.version = "v5")

  out <- CreateAssay5Object(counts = counts(obj), data = logcounts(obj))
  out <- CreateSeuratObject(out)
  out@meta.data <- cbind(out@meta.data, as.data.frame(colData(obj)[, c(batch, cell_type)]))

  out[["RNA"]] <- split(out[["RNA"]], f = out@meta.data[, batch])
  out@reductions[[tolower(reduction)]] <- CreateDimReducObject(embeddings = reducedDim(obj, reduction),
                                                               loadings = attr(reducedDim(obj, reduction), "rotation"),
                                                               key = paste0(reduction , "_"), assay = 'RNA')
  out <- ScaleData(out, verbose = FALSE)
  return(out)
}
