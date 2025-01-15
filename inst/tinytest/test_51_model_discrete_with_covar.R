# -------------------------------------------------------------------
# Testing 'transitreg' for discrete response (count data).
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Example from ?transitreg with 20 bins only for testing
set.seed(123)
n <- 1000
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = exp(-1 + cos(x)))

## Fit model with continuous response.
mod <- transitreg(y ~ s(theta) + s(x) + te(x, theta), breaks = 20, verbose = FALSE)
expect_inherits(mod, "transitreg", info = "Check model class")

## Setting up new data.frame and new response,
## predict CDF and PDF to compare them. The cumulative sum of the
## PDF should result in the CDF.
nd <- data.frame(x = rep(mean(x), length(mod$bins)))
ny <- seq_along(mod$bins) - 1

## TODO(R): Write more tests once the S3 method is back
## ## Calulating CDF and PDF
## expect_silent(mcdf <- predict(mod, newdata = nd, y = ny, type = "cdf"),
##               info = "Calculating cdf (silent)")
## expect_true(is.numeric(mcdf), info = "Return class")
## expect_identical(length(mcdf), nrow(nd), info = "Length of cdf return")
## 
## expect_silent(mpdf <- predict(mod, newdata = nd, y = ny, type = "pdf"),
##               info = "Calculating pdf (silent)")
## expect_true(is.numeric(mpdf), info = "Return class")
## expect_identical(length(mpdf), nrow(nd), info = "Length of pdf return")

