#' Object for scanorama
#'
#' @param obj SingleCellExperiment object
#' @param batch batch
#'
#' @importFrom SingleCellExperiment colData logcounts
#' @export
#'
scanorama <- function(obj, batch){

  ll <- lapply(unique(colData(obj)[, batch]), function(i) obj[, colData(obj)[, batch] == i])

  assaylist <- list()
  genelist <- list()
  for(i in seq_along(ll)) {
    assaylist[[i]] <- t(as.matrix(logcounts(ll[[i]])))
    genelist[[i]] <- rownames(ll[[i]])
  }

  return(c(assaylist, genelist))
}

