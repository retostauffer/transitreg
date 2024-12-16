\name{predict.tm}
\alias{predict.tm}

\title{
  Predict Method for Transition Model Fits
}

\description{
  Provides predictions for transition models fitted using the \code{\link{tm}} function.
  Predictions can be generated for the probability density function (PDF), cumulative 
  distribution function (CDF), maximum probability, or specific quantiles of the response 
  distribution.
}

\usage{
\method{predict}{tm}(object, newdata = NULL, y = NULL, prob = NULL,
  type = c("pdf", "cdf", "pmax", "quantile"), ...)
}

\arguments{
  \item{object}{An object of class \code{"tm"} resulting from a call to \code{\link{tm}}.}
  \item{newdata}{Optional. A data frame containing new predictor values.
    If not provided, predictions are made for the data used in fitting the model.}
  \item{y}{Optional. A vector of response values for which the PDF or CDF should
    be computed. Required if \code{type} is \code{"pdf"} or \code{"cdf"}.}
  \item{prob}{Optional. A numeric value specifying the quantile to compute when
    \code{type = "quantile"}. Default is \code{0.5} (median). If provided,
    argument \code{type} is set to \code{type = "quantile"}.}
  \item{type}{Character. Specifies the type of prediction to return:
    \describe{
      \item{\code{"pdf"}}{The predicted probability density function (PDF).}
      \item{\code{"cdf"}}{The cumulative distribution function (CDF).}
      \item{\code{"pmax"}}{The expected value of the response (maximum probability).}
      \item{\code{"quantile"}}{The quantile of the response specified by \code{prob}.}
    }
  }
  \item{\dots}{Additional arguments passed to the prediction function.}
}

\details{
  The \code{predict.tm} method computes predictions based on the transition model fit. 
  Predictions can be made for the original training data or for new data provided via 
  \code{newdata}. The method also supports scaling of covariates if scaling was applied 
  during model fitting (argument \code{scale.x} in function \code{\link{tm}}).

  For \code{"pdf"} and \code{"cdf"}, the response values (\code{y}) must be provided
  unless the model was fit with those values already included. For \code{"quantile"},
  a specific quantile is computed based on \code{prob}.
}

\value{
  Returns predictions of the specified type:
  \itemize{
    \item For \code{"pdf"} and \code{"cdf"}, a vector or matrix of probabilities.
    \item For \code{"pmax"}, the expected value of the response.
    \item For \code{"quantile"}, the quantile of the response distribution at
      the specified \code{prob}.
  }
}

\seealso{
  \code{\link{tm}}, \code{\link{tm_data}}, \code{\link{tm_dist}}
}

\examples{
## Example: Predicting PDF and CDF.
set.seed(123)
n <- 500
x <- runif(n, -3, 3)
y <- rpois(n, exp(2 + sin(x)))
b <- tm(y ~ s(theta) + s(x))

## Predict PDF and CDF.
p <- list()
p$pdf <- predict(b, type = "pdf", y = 3)
p$cdf <- predict(b, type = "cdf", y = 3)

## Predict maximum probability (expected value).
p$pmax <- predict(b, type = "pmax")

## Predict quantiles.
p$qu95 <- predict(b, prob = 0.95)

print(head(as.data.frame(p)))

## Visualize predictions.
nd <- data.frame(x = seq(-3, 3, length = 100))

## Predict quantiles.
qu <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 0.95, 0.99)
p <- lapply(qu, function(prob) {
  predict(b, newdata = nd, prob = prob)
})

## Plot data and fitted quantiles.
plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
matplot(nd$x, do.call("cbind", p),
  type = "l", col = 4, lwd = 2, lty = 1,
  add = TRUE)
}
\keyword{methods}
\keyword{models}

