% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{G_vec}
\alias{G_vec}
\title{Vectorized call to the G1 / G2 functions}
\usage{
G_vec(z, g, type)
}
\arguments{
\item{z}{Vector}

\item{g}{String, either G1, G1_prime, G2_curly, G2 or G2_curly}

\item{type}{Numeric, for G1: 1-2; G2: 1-5
(see \link{G1_fun}, \link{G1_prime_fun}, \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})}
}
\description{
Vectorized call to the G1 / G2 functions
}
\keyword{internal}
