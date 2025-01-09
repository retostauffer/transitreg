\name{tmutils}
\alias{tmutils}
\alias{cdf_to_tp}
\alias{tp_to_cdf}
\alias{tp_to_pdf}

\title{
    Transiton Model Utilities
}

\description{
    A transiton model consists of a series of transition probabilities (TP)
    which define the probability that an observation falls into the
    next hither bin. Based on these transition probabilities, the
    cummulative distribution function (CDF) as well as the density (PDF)
    can be calculated.
    These utilities allow to convert from CDF to TPs, and from TPs
    back to CDF/PDF. Internally used for testing.
}

\usage{
    cdf_to_tp(x)
    tp_to_cdf(tp)
    tp_to_pdf(tp)
}

\arguments{
    \item{x}{Numeric vector with probabilities (CDF) in \code{[0., 1.]},
        monotonically increasing.}
    \item{tp}{Numeric vector of transition probabilities. Must be
        in \code{[0., 1.]} and monotonically decreasing.}
}

\value{
    A numeric vector of the same length as \code{x}/\code{tp}.
}

\examples{
## For testing:
## Draw PDF and CDF from Poisson distribution
p <- ppois(0:15, lambda = 4)
d <- dpois(0:15, lambda = 4)

## Convert Poisson CDF to transition probabilities
tp <- TransitionModels:::cdf_to_tp(p) 
round(tp, 2)

## Convert transition probabilities back to CDF
p2 <- TransitionModels:::tp_to_cdf(tp)
all.equal(p, p2) # Must be equal

## Convert transition probabilities to PDF
devtools::load_all(); d2 <- tp_to_pdf(tp)
all.equal(d, d2) # Must be equal

## Quick visual representation
barplot(tp, col = "steelblue", main = "Transition Probabilities")
col <- c("gray80", "tomato")
barplot(rbind(p, p2), beside = TRUE, main = "CDF Comparison", col = col)
barplot(rbind(d, d2), beside = TRUE, main = "PDF Comparison", col = col)
}

\concept{transition modeling}

