\name{cdf_to_tp}
\alias{cdf_to_tp}

\title{
    Convert CDF to Transition Probabilities
}

\description{
    A transiton model consists of a series of transition probabilities
    which define the probability that an observation falls into the
    next hither bin. Based on these transition probabilities, the
    cummulative distribution function (CDF) can be calculated.
    This utility function does the reverse; it takes a CDF vector
    and calculates transition probabilities. Used internally for testing.
}

\usage{
    cdf_to_tp(x)
}

\arguments{
    \item{x}{Numeric vector with probabilities (CDF) in \code{[0., 1.]},
        monotonically increasing.}
}

\value{
    A numeric vector of the same length as \code{x} with
    transition probabilities based on the CDF vector \code{x}.
}

\examples{
## CDF from a binomial distribution
p <- pbinom(0:10, size = 10, prob = 0.1)

## Calculate transition probabilities
tp <- TransitionModels:::cdf_to_tp(p)

## Quick visual representation
barplot(p, main = "CDF")
barplot(tp, col = "steelblue", main = "Transition Probabilities")
}

\concept{transition modeling}

