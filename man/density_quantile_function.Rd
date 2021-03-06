% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nuisance_quantities.R
\name{density_quantile_function}
\alias{density_quantile_function}
\title{Density Quantile Function}
\usage{
density_quantile_function(y, x, u, alpha, sparsity, bandwidth_estimator)
}
\arguments{
\item{y}{Vector of dependent data}

\item{x}{Matrix of covariates}

\item{u}{Quantile residuals}

\item{alpha}{Probability level}

\item{sparsity}{iid or ind}

\item{bandwidth_estimator}{Bofinger, Chamberlain or Hall-Sheather}
}
\description{
Estimate the density quantile function
}
\references{
For the iid and nid method, see Koenker (1994), and Hendricks and Koenker (1992).
For the bandwidth types, see Bofinger (1975), Chamberlain (1994), and Hall and Sheather(1988).
}
\keyword{internal}
