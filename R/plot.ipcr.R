#' @title Plot Correlation between Covariates and IPCs
#' @description This functions plots a heat map that visualizes the correlation between covariates and IPCs.
#' @param x an ipcr object.
#' @param method character string indicating whether the relationship between covariates and initial IPCs (method = “initial”, default option) or the iterated IPCs (method = “iterated”) is drawn.
#' @examples
#' # Plot the correlation between covariates and initial IPCs.
#' plot(x = IPC_reg)
#' @export plot.ipcr

plot.ipcr <- function(x, parameter = NULL, method = "iterated", ...) {

  # Check arguments
  if (method != "standard" & method != "iterated") {
    stop("'method' must be either 'initial' or 'iterated'")
  }

  if (method == "iterated" & is.null(x$iterated_regression)) {method <- "standard"}

  # Get covariates

  # Compute correlation matrix
  COR <- cor(x = x$IPCs, y = x$regression[[1]]$model[, 2:NCOL(x$regression[[1]])])

  # Transform data into longformatio
  longData <- melt(COR)

  # Get names
  if (method == "standard") {
    ipc_name <- "Standard IPCs"
  }
  if (method == "iterated") {
    ipc_name <- "Iterated IPCs"
  }

  ggplot(data = longData, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",  midpoint = 0,
                         limit = c(-1, 1), space = "Lab", name = "Corr.") +
    labs(x = "Covariates", y = "Parameters", title = ipc_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    coord_fixed()
}
