#' @export
nobs.MxModel <- function(object, ...) {
  as.integer(NROW(openmx_raw_data_ipcr(object)))
}
