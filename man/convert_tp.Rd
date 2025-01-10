\name{convert_tp}
\alias{convert_tp}
\alias{cdf_to_tp}
\alias{tp_to_cdf}
\alias{tp_to_pdf}

\title{
    Convert between Density, Probabilities, and Transition Probabilities
}

\description{
    A transiton model consists of a series of transition probabilities (TP)
    which define the probability that an observation falls into the next hither
    bin. Based on these transition probabilities, the cummulative distribution
    function (CDF) as well as the density (PDF) can be calculated.
    This utility function allows to convert from CDF to TPs, and from TPs
    back to CDF/PDF.
}

\usage{
    convert_tp(x, from, to, width = NULL, drop = TRUE)
}

\arguments{
    \item{x}{numeric vector with density, probabilities, or transition probabilities.}
    \item{from}{character, current type of data in \code{x}.}
    \item{to}{character, type into which \code{x} should be converted.}
    \item{width}{\code{NULL} or numeric (either length 1 or same length as \code{x}).
        Width of the individual bins represented in \code{x}.}
    \item{drop}{if \code{TRUE} (default) the result is simplified if possible.}
}

\value{
    If \code{to} is a single character and \code{drop = TRUE}, the return is a numeric vector
    of the same length as the input argument \code{x}. Else the return is
    a \code{data.frame} with the same number of rows as \code{length(x)}.
}

\examples{
## For testing:
## Draw PDF and CDF from Poisson distribution
p <- ppois(0:10, lambda = 4)
d <- dpois(0:10, lambda = 4)

## Convert Poisson CDF to transition probabilities
tp <- convert_tp(p, from = "cdf", to = "tp")
round(tp, 2)

pd <- convert_tp(tp, from = "tp", to = c("pdf", "cdf"))

## Convert transition probabilities back to CDF
p2 <- convert_tp(tp, from = "tp", to = "cdf")
all.equal(p, p2) # Must be equal

## Convert transition probabilities to PDF
d2 <- convert_tp(tp, from = "tp", to = "pdf")
all.equal(d, d2) # Must be equal

## Or directly from tp to both, CDF and PDF
d2 <- convert_tp(tp, from = "tp", to = c("cdf", "pdf"))
all.equal(p, d2$cdf)
all.equal(d, d2$pdf)

## Quick visual representation
barplot(tp, col = "steelblue", main = "Transition Probabilities")
col <- c("gray80", "tomato")
barplot(rbind(p, p2), beside = TRUE, main = "CDF Comparison", col = col)
barplot(rbind(d, d2), beside = TRUE, main = "PDF Comparison", col = col)
}

\concept{transition modeling}

