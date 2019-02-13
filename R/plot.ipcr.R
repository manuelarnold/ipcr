#' @title Plot Correlation between Covariates and IPCs
#' @description This functions plots a heat map that visualizes the correlation between covariates and IPCs.
#' @param x an ipcr object.
#' @param ipc_method character string indicating whether the relationship between covariates and initial IPCs (method = “initial”, default option) or the iterated IPCs (method = “iterated”) is drawn.
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman".
#' @examples
#' # Plot the correlation between covariates and initial IPCs.
#' plot(x = IPC_reg)
#' @export plot.ipcr

plot.ipcr <- function(x, ipc_method = "initial", cor_method = "pearson") {

  # Check arguments
  if (ipc_method != "initial" & ipc_method != "iterated") {
    stop("'ipc_method' must be either 'initial' or 'iterated'")
  }

  if (cor_method != "pearson" & cor_method != "kendall" &
      cor_method != "spearman") {
    stop("'cor_method' must be either 'pearson', 'kendall', or 'spearman'")
  }

  # Compute correlation matrix
  COR <- cor(x = x$InitialIPCs, y = x$Covariates, method = cor_method)

  # Transform data into longformatio
  longData <- melt(COR)

  # Get names
  if (ipc_method == "initial") {
    ipc_name <- "Initial IPCs"
  }
  if (ipc_method == "iterated") {
    ipc_name <- "Iterated IPCs"
  }

  if (cor_method == "pearson") {
    cor_name <- "Pearson\ncorrelation"
  }
  if (cor_method == "kendall") {
    cor_name <- "Kendall\ncorrelation"
  }
  if (cor_method == "spearman") {
    cor_name <- "Spearman\ncorrelation"
  }

  ggplot(data = longData, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = cor_name) +
    labs(x = "Covariates", y = "Parameters", title = ipc_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1)) +
    coord_fixed()
}
