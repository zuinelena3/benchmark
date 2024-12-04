pre_bbknn <- function(obj, batch, reduction){
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, batch, reduction) {
    bbknn <- import("bbknn")
    sc <- import("scanpy")
    anndata <- import("anndata")

    adata <- SCE2AnnData(obj)
    adata$X <- adata$layers["logcounts"]
    adata$obs$batch <- colData(obj)[, batch]
    adata$obsm["X_pca"] <- reducedDim(obj, reduction)
    return(adata)
  }, obj = obj, batch = batch, reduction = reduction)
}
