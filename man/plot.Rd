\name{plot.transitreg}
\alias{plot.transitreg}

\title{
  Plot Method for Transition Model Fits
}

\description{
  Provides diagnostic and visualization plots for transition models fitted
  using the \code{\link{transitreg}} function. The method supports plotting effects and
  quantile residual diagnostic plots.
}

\usage{
\method{plot}{transitreg}(x, which = "effects", spar = TRUE, k = 5, ...)
}

\arguments{
  \item{x}{An object of class \code{"transitreg"} resulting from a call to \code{\link{transitreg}}.}
  \item{which}{A character string or integer specifying the type of plot(s) to
    generate. Options include:
    \describe{
      \item{\code{"effects"}}{Plots the effects of the predictors on the response.
        Requires that the model is estimated by \code{\link[mgcv]{gam}}
        or \code{\link[mgcv]{bam}}.}
      \item{\code{"hist-resid"}}{Plots a histogram of the qauntile residuals.}
      \item{\code{"qq-resid"}}{Generates a Q-Q plot of the quantile residuals.}
      \item{\code{"wp-resid"}}{Creates a worm plot of the quantile residuals.}
    }
    Multiple options can be specified as a character vector or numeric indices.}
  \item{spar}{Logical. If \code{TRUE}, multiple plots are arranged in a
    single window. Default is \code{TRUE}.}
  \item{k}{Integer, TODO(N): Describe argument. Defaults to \code{5}.}
  \item{\dots}{Additional arguments passed to the underlying plotting functions.}
}

\details{
  The \code{plot.transitreg} method provides flexible visualization options for
  evaluating transition model fits. Users can choose to:
  \itemize{
    \item Visualize the effects of predictors on the response variable
      (if the model is a GAM, see \code{\link[mgcv]{gam}}).
    \item Evaluate quantile residuals through histograms, Q-Q plots, or worm plots.
  }

  The \code{which} argument controls the type of plots generated. By default,
  the \code{"effects"} plot is shown if the model supports it. Residual-based
  plots (\code{"hist-resid"}, \code{"qq-resid"}, \code{"wp-resid"}) provide
  insights into model calibration.
}

\value{
  Returns \code{NULL} invisibly. Generates plots as a side effect.
}

\seealso{
  \code{\link{transitreg}}, \code{link{transitreg_dist}}, \code{link{transitreg_data}}, \code{\link{predict.transitreg}}
}

\examples{
## Example: Fit a transition model and generate plots.
set.seed(123)
n <- 500
x <- runif(n, -3, 3)
y <- rpois(n, exp(2 + sin(x)))
b <- transitreg(y ~ s(theta) + s(x))

## Plot effects.
plot(b, which = "effects")

## Plot residuals.
plot(b, which = c("hist-resid", "qq-resid"))

## Custom plot layout.
par(mfrow = c(2, 1))
plot(b, which = 3:4, spar = FALSE)
}

\keyword{methods}
\keyword{models}
\keyword{visualization}


