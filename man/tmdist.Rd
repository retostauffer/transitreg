\name{tmdist}
\alias{tmdist}
\alias{random.tmdist}
\alias{pdf.tmdist}
\alias{log_pdf.tmdist}
\alias{cdf.tmdist}
\alias{quantile.tmdist}
\alias{median.tmdist}
\alias{mean.tmdist}
\alias{as.matrix.tmdist}
\alias{support.tmdist}
\alias{is_discrete.tmdist}
\alias{is_continuous.tmdist}
\title{Create an tmdist Distribution}
\usage{
    tmdist(z, newdata = NULL, newresponse = NULL)

    \method{pdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
    \method{log_pdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
    \method{cdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
    \method{quantile}{tmdist}(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
    \method{median}{tmdist}(x, na.rm = NULL, ncores = NULL, ...)
    \method{mean}{tmdist}(x, ncores = NULL, ...)
    \method{random}{tmdist}(x, n = 1L, drop = TRUE, ...)

    \method{is_discrete}{tmdist}(d, ...)
    \method{is_continuous}{tmdist}(d, ...)
    \method{support}{tmdist}(d, drop = NULL, ...)

    \method{as.matrix}{tmdist}(x, expand = FALSE, ...)
}
\arguments{
    \item{z}{an object that can be converted into a \code{c("tmdist", "distributions")}
        object. See section 'Details' for more information.}
    \item{x}{object of class \code{tmdist} (if main argument), else a numeric
        vector at which the CDF, PDF, or quantile should be evaluated.}
    \item{probs}{Probabilities at which to evaluate the function.}
    \item{...}{Currently not used.}
    \item{n}{The number of samples to draw. Defaults to `1L`.}
    \item{drop}{logical. Should the result be simplified to a vector if possible?}
    \item{d}{An `tmdist` object created by a call to [tmdist()].}
    \item{elementwise}{logical. Should each distribution in \code{x} be evaluated
        at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
        Or, if \code{x} and \code{probs} have the same length, should the evaluation be
        done element by element (\code{elementwise = TRUE}, yielding a vector)? The
        default of \code{NULL} means that \code{elementwise = TRUE} is used if the
        lengths match and otherwise \code{elementwise = FALSE} is used.}
    \item{ncores}{number of cores to use with OpenMP. If \code{NULL} (default)
        all but two available cores will be used if OpenMP is supported.}
    \item{expand}{Logical, if \code{FALSE} the result is a 'wide' matrix,
        if \code{TRUE} a 'long' format is returned.}
    \item{na.rm}{\code{NULL}, currently has no effect.}
    \item{newdata}{An optional data frame in which to look for variables with
          which to predict. If omitted, the fitted values are used.}
    \item{newresponse}{Numeric vector with new response, must be provided
        if \code{newdata} is specified, same length as \code{nrow(newdata)}.}
}
\value{
A `tmdist` object.

TODO(R): Missing.
}
\description{
TODO(R): Missing
}
\details{
TODO(R): Missing

Explain what 'z' can take.
}
\examples{
## Creating 'fake' Transition dist using Poisson distribution, converted to
## transition models using ## cdf_to_tp function.
p  <- ppois(0:30, lambda = 5) # CDF
tp <- TransitionModels:::cdf_to_tp(p) # Transition probabilities

## Create named list; one form to create a 'tmdist' object.
## Requires 'tp' (transition probabilities), as well as lower and
## upper boundaries of the bins.
tmp <- list(tp = tp,
            lo = seq(-1.25, by = 0.5, length.out = length(tp)),
            up = seq(-1.25, by = 0.5, length.out = length(tp)) + 0.5)
str(tmp)

## Convert to 'tmdist' distributions object
d <- tmdist(tmp)
print(d)

## Calculate PDF and CDF (on the center of the bins used to create the distribution)
x <- seq(-1, by = 0.5, length.out = length(tp)) # between 'lo' and 'hi'
d_cdf <- cdf(d, x)
d_pdf <- pdf(d, x)

barplot(rbind(d_cdf, d_pdf), beside = TRUE, col = c("tomato", "gray30"),
        las = 2, main = "CDF and PDF of tmdist")

## Mean and median
c("mean" = mean(d), "median" = median(d))

## Draw 1000 random values; compare empirical mean/median
## TODO(R): mean(d) and mean(r) differ quite a bit; is mean(d) correctly implemented?
r <- random(d, n = 1000L)
c("empirical mean" = mean(r), "empirical median" = median(r))

## Support and type of distribution
is_discrete(d)
is_continuous(d)
support(d)
}
\concept{tmdist distribution}
\concept{Normal distribution}
\concept{distributions}
