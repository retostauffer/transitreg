\name{tm_detect_cores}
\alias{tm_detect_cores}

\title{
    Detect availability cores available for OMP
}

\description{
    To increase performance, several functions are implemented in \code{C},
    using OpenMP for parallelization. This is a convenience function to
    (i) check if OpenMP is available (the package has been compiled with
    OpenMP), and (ii) is used internally to auto-detect the number of available
    cores.
}

\usage{
tm_detect_cores(verbose = TRUE)
}

\arguments{
    \item{verbose}{Logical, if \code{TRUE} a short message is printed.}
}

\value{
    If OpenMP is not available, \code{0L} is returned. Else, the maximum
    number of cores (processors) is returned.
}

\seealso{
  \code{\link{tm}}
}

\examples{
    tm_detect_cores()
    tm_detect_cores(verbose = FALSE)
}
