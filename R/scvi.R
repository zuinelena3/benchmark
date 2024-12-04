#' scVi
#'
#' @param obj SingleCellExperiment object
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom reticulate import
#' @export
#'
scvi <- function(obj) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, batch, assay, dims) {
    scvi <- import("scvi")

    andata <- SCE2AnnData(obj)
    andata$layers["counts"] <- andata$X
    scvi$model$SCVI$setup_anndata(andata, layer = assay, batch_key = batch)
    return(andata)
  }, obj = obj, batch = batch, assay = assay, dims = dims)
  return(out)
}
