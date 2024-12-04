pre_bbknn <- function(obj, adata, batch, reduction){
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(adata, obj, batch, reduction) {
    bbknn <- import("bbknn")
    sc <- import("scanpy")
    anndata <- import("anndata")

    sc$tl$umap(adata)
    bbknn_umap <- adata$obsm[["X_umap"]]
    reducedDim(obj, "UMAP_bbknn") <- bbknn_umap

    bbknn <- knn_index_dist(dist = adata$obsp[["distances"]])
    perplexity <- ncol(x = bbknn$idx) - 1
    tsne <- Rtsne_neighbors(
      index = bbknn$idx,
      distance = bbknn$dist,
      perplexity = perplexity
    )$Y
    colnames(x = tsne) <- paste0("tSNE_", c(1, 2))
    rownames(x = tsne) <- colnames(obj)
    reducedDim(obj, "TSNE_bbknn") <- tsne
  }, obj = obj, adata =adata, batch = batch, reduction = reduction)
}
