# -------------------------------------------------------------------
# Testing 'tm' for discrete response (count data).
# TODO(R): Remove useC and comparison against the R version
#          once we removed that.
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("TransitionModels"))


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Example from ?tm with 20 bins only for testing
set.seed(123)
n <- 1000
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = exp(-1 + cos(x)))

## Fit model with continuous response.
mod <- tm(y ~ s(theta) + s(x) + te(x, theta), breaks = 20, verbose = FALSE)
expect_inherits(mod, "tm", info = "Check model class")

## Setting up new data.frame and new response,
## predict CDF and PDF to compare them. The cumulative sum of the
## PDF should result in the CDF.
nd <- data.frame(x = rep(mean(x), length(mod$bins)))
ny <- seq_along(mod$bins) - 1

## Calulating CDF and PDF
expect_silent(mcdf <- predict(mod, newdata = nd, y = ny, type = "cdf"),
              info = "Calculating cdf (silent)")
expect_true(is.numeric(mcdf), info = "Return class")
expect_identical(length(mcdf), nrow(nd), info = "Length of cdf return")

expect_silent(mpdf <- predict(mod, newdata = nd, y = ny, type = "pdf"),
              info = "Calculating pdf (silent)")
expect_true(is.numeric(mpdf), info = "Return class")
expect_identical(length(mpdf), nrow(nd), info = "Length of pdf return")

## Comparing cumsum(pdf) to cdf
expect_equal(cumsum(mpdf), mcdf, tolerance = 1e-15,
             info = "Comparing cumulative PDF to its CDF")

