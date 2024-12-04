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

    if (adata$n_obs > 100000) {
      neighbors_within_batch = 25
    } else neighbors_within_batch = 3
    return(obj)
  }, obj = obj, batch = batch, reduction = reduction)
}
