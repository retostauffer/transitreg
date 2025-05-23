% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{transitreg_get_number_of_cores}
\alias{transitreg_get_number_of_cores}
\title{Get number of cores for OpenMP}
\usage{
transitreg_get_number_of_cores(ncores = NULL, verbose = FALSE)
}
\arguments{
\item{ncores}{\code{NULL} or a positive integer.}

\item{verbose}{logical, if \code{TRUE} a message is shown.}
}
\value{
Number of cores to be used in OpenMP parallelization (integer).
}
\description{
Some parts of the package use C routines which allow for parallelization
using OpenMP. This function is used to specify how many cores to be used.
}
\details{
If \code{ncores} is \code{NULL} the number of available
cores is auto-detected and set to 'total number of cores - 2'.
If integer, it is checked if this number of cores is available,
else set tot he 'total number of cores available'.
}
\author{
Reto
}
