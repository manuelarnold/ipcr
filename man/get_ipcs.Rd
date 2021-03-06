% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ipcs.R
\name{get_ipcs}
\alias{get_ipcs}
\title{Calcuate Individual Parameter Contributions}
\usage{
get_ipcs(x)
}
\arguments{
\item{x}{a fitted model object.}
}
\value{
A \code{data.frame} containing the of a fitted model. The dimensions of this
\code{data.frame} are \eqn{n * k}, where \eqn{n} denotes the number of observations and
\eqn{k} the number of parameters. The columns should be named as in
\code{\link[stats]{coef}}.
}
\description{
This functions calculates the individual parameter contributions (IPCs) of
 a fitted model.
}
\details{
\code{get_ipcs} is a convenience function for obtaining IPCs. The more
powerful \code{\link[ipcr]{ipcr}} function also provides IPCs and can be used to predict
differences in model parameters by regressing the IPCs on covariates.
}
\references{
Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2019). Identifying
heterogeneity in dynamic panel models with individual parameter contribution
regression. \emph{Structural Equation Modeling, 27}, 613-628. doi:
\href{https://doi.org/10.1080/10705511.2019.1667240}{10.1080/10705511.2019.1667240}
}
\seealso{
\code{\link[ipcr]{ipcr}}
}
