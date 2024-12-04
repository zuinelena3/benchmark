#' Batch effect correction methods
#'
#' A common interface for single-cell batch correction methods.
#'
#' @param obj a \linkS4class{SingleCellExperiment} object or list containing single-cell gene expression matrices
#' @param batch a string specifying the batch
#' @param assay a string specifying the assay to use for correction
#' @param hvgs a vector specifying which features to use for correction
#' @param dims number of dimensions to use for dimension reduction
#' @param reduction a string specifying the dimension reduction to use for correction
#' @param anchor a string specifying the anchor integration type for Seurat methods
#' @param k_anchor number of anchors
#' @param genelist list containing the genes in each batch
#' @param cell_type string specifying the cell-type labels
#' @param METHOD a \code{MethodParam} object specifying the batch correction method
#' @param alt_out alternative output: a \code{\link{AltOutput}} class
#'
#' @import methods
#' @rdname scIntegration
#' @importFrom limma removeBatchEffect
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData assays
#'
setMethod("scIntegration", "limmaMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                    dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                    genelist = NULL, cell_type = NULL, METHOD) {
  out <- removeBatchEffect(x = logcounts(obj), batch = colData(obj)[, batch])
  return(out)
})

#' @rdname scIntegration
#' @importFrom SummarizedExperiment assays colData
#' @importFrom sva ComBat_seq
#'
setMethod("scIntegration", "combatMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                    dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                    genelist = NULL, cell_type = NULL, METHOD) {
  out <- ComBat_seq(counts = as.matrix(counts(obj)), batch =  colData(obj)[, batch])
  return(out)
})

#' @rdname scIntegration
#' @param obj A SingleCellExperiment object
#' @param anchor a string specifying the anchors finding type (cca, rpca, jpca, rlsi)
#' @param k_anchor number of anchors (default: k_anchor = 5)
#'
#' @importFrom Seurat FindIntegrationAnchors IntegrateData GetAssayData SplitObject VariableFeatures<- CreateDimReducObject
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SingleCellExperiment counts logcounts colData reducedDim
#' @importFrom scater runPCA
#'
setMethod("scIntegration", "seuratv3Method", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                      dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                      genelist = NULL, cell_type = NULL, METHOD) {
  anchorset <- FindIntegrationAnchors(object.list = so_ll, reduction = anchor, anchor.features = hvgs, dims = 1:dims,  k.anchor = k_anchor, verbose = FALSE)
  out <- IntegrateData(anchorset = anchorset, verbose = FALSE)

  return(out)
})

#' @rdname scIntegration
#' @param anchor a string specifying the anchors finding type (CCAIntegration, RPCAIntegration, HarmonyIntegration, JointPCAIntegration)
#' @param k_anchor number of anchors (default: k_anchor = 5)
#'
#' @importFrom Seurat ScaleData CreateDimReducObject IntegrateLayers
#' @importFrom SeuratObject CreateAssay5Object JoinLayers CreateSeuratObject
#' @importFrom SingleCellExperiment reducedDim counts logcounts
#' @importFrom SummarizedExperiment colData
#'
setMethod("scIntegration", "seuratv5Method", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                      dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                      genelist = NULL, cell_type = NULL, METHOD) {
  out <- IntegrateLayers(object = out, method = anchor, orig.reduction = tolower(reduction),
                         new.reduction = "integrated", features = hvgs, k.anchor = k_anchor, verbose = FALSE)
  return(out)
})

#' @rdname scIntegration
#' @param obj a SingleCellExperiment object  or a list containing SingleCellExperiment objects
#' @param batch a string specifying the batch. batch = NULL when obj is a list
#'
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assays colData
#' @importFrom Seurat as.sparse
#'
setMethod("scIntegration", "fastMNNMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD) {
  out <- fastMNN(obj = obj, batch = colData(obj)[batch][, 1], subset.row = hvgs, d = dims)
  return(out)
})

#' @rdname scIntegration
#' @importFrom harmony RunHarmony
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assays colData
#' @importFrom Seurat as.sparse
#'
setMethod("scIntegration", "harmonyMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD) {
  out <- RunHarmony(obj, batch)
  return(out)
})

#' @rdname scIntegration
#' @param obj a list of "matrix" "array" objects
#' @param genelist a list of genes
#' @param hvgs number of highly variable genes
#'
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#'
setMethod("scIntegration", "scanoramaMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                       dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                       genelist = NULL, cell_type = NULL, METHOD) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))
  out <- basiliskRun(proc = proc, fun = function(obj, genelist, hvgs) {
    scanorama <- import("scanorama")
    method <- scanorama$correct(obj, genelist, return_dimred = TRUE, return_dense = TRUE, verbose = FALSE, hvg = hvgs, dimred = as.integer(dims))
  }, obj = obj, genelist = genelist, hvgs = hvgs)

  return(out)
})

#' @rdname scIntegration
#' @param obj A SingleCellExperiment
#'
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom Rtsne Rtsne_neighbors
#'
setMethod("scIntegration", "bbknnMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                   dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                   genelist = NULL, cell_type = NULL, METHOD) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(adata, batch, reduction) {
    bbknn <- import("bbknn")
    sc <- import("scanpy")
    anndata <- import("anndata")

    bbknn$bbknn(adata, batch_key = batch, neighbors_within_batch = as.integer(neighbors_within_batch))

    return(obj)
  }, obj = adata, batch = batch, reduction = reduction)
  return(out)
})

#' @rdname scIntegration
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom reticulate import
#'
setMethod("scIntegration", "scVIMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                  dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                  genelist = NULL, cell_type = NULL, METHOD) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(andata, batch, assay, dims) {
    scvi <- import("scvi")

    model_scvi <- scvi$model$SCVI(andata, n_latent = as.integer(dims))
    model_scvi$train()

    andata$obsm["X_scVI"] <- model_scvi$get_latent_representation()
    return(andata)
  }, obj = andata, batch = batch, assay = assay, dims = dims)

  return(out)
})

#' @rdname scIntegration
#' @param genelist negative controls
#'
#' @importFrom scMerge scMerge2
#' @importFrom SingleCellExperiment logcounts colData
#'
setMethod("scIntegration", "scMergeMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD) {
  if(is.null(cell_type)){
    out <- scMerge2(exprsMat = logcounts(obj), batch = colData(obj)[, batch], ctl = genelist,
                    verbose = FALSE)
  }
  else {
    out <- scMerge2(exprsMat = logcounts(obj), batch = colData(obj)[, batch], ctl = genelist,
                    cellTypes = colData(obj)[, cell_type], verbose = FALSE)
  }
  return(out)
})
