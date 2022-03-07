iterated_ipcr <- function(x, ...) {
  UseMethod("iterated_ipcr")
}

iterated_ipcr.lm <- function(x, ...) {
  stop("Iterated IPC regression has not been implemented yet for lm-type models.")
}

iterated_ipcr.glm <- function(x, ...) {
  stop("Iterated IPC regression has not been implemented yet for glm-type models.")
}
