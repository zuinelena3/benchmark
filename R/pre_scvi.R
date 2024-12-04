#' scVi
#'
#' @param obj SingleCellExperiment object
#' @param batch batch
#' @param assay assay
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom reticulate import
#' @export
#'
pre_scvi <- function(obj, batch, assay) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, batch, assay) {
    scvi <- import("scvi")

    andata <- SCE2AnnData(obj)
    andata$layers["counts"] <- andata$X
    scvi$model$SCVI$setup_anndata(andata, layer = assay, batch_key = batch)
    return(andata)
  }, obj = obj, batch = batch, assay = assay)
  return(out)
}
