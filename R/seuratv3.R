#' SeuratV3 convertion
#'
#' @param obj SingleCellExperiment objcet to convert into a Seurat V3
#'
#' @importFrom Seurat SplitObject VariableFeatures<- CreateDimReducObject
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SingleCellExperiment counts logcounts colData reducedDim
#'
#' @export
#'
sce_to_seuratv3 <- function(obj){
  options(Seurat.object.assay.version = "v3")

  so <- CreateSeuratObject(counts = counts(obj))
  so@assays$RNA$data <- logcounts(obj)
  so@meta.data <- cbind(so@meta.data, as.data.frame(colData(obj)[, c(batch, cell_type)]))
  VariableFeatures(so) <- rownames(so)

  so@reductions[[tolower(reduction)]] <- CreateDimReducObject(embeddings = reducedDim(obj, reduction),
                                                              loadings = attr(reducedDim(obj, reduction), "rotation"),
                                                              key = paste0(tolower(reduction), '_'), assay = 'RNA')

  so_ll <- SplitObject(so, split.by = batch)
  return(so_ll)
}
