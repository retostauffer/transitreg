# -------------------------------------------------------------------
# Testing implementation of CDF, PDF, pmax
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# -------------------------------------------------------------------
# Drawing random values from a Normal distribution, estimate a
# Transition Model and compare the predicted CDF, PDF, and quantile
# against the true CDF/PDF/Quantiles.
# -------------------------------------------------------------------

set.seed(6020)
truth <- list(mu = 77, sd = 13, N = 1000)

# Creating data.frame with training data
data <- data.frame(y = rnorm(truth$N, truth$mu, truth$sd),
                   x = runif(truth$N)) # Completely uninformative

# Specify a searies of breaks for the binning of the continuous
# response. This covers a range wider than the range of data$y
# to be able to calculate the probabilities for the 'entire'
# distribution.
breaks <- seq(0, 200, length.out = 301)

expect_silent(mod <- transitreg(y ~ s(theta), data = data, breaks = breaks),
              info = "Estimating transitional regression model")
expect_inherits(mod, "transitreg", info = "Checking transitreg() return class")
expect_identical(mod$bins, breaks, info = "Testing if the breaks have been stored correctly")
# TODO(R): Currently it stores the bins on breaks _and_ bins?

# Calculating center of the bins
binmid <- (head(mod$bins, -1) + tail(mod$bins, -1)) / 2
expect_identical(length(binmid), length(mod$bins) - 1L) # <- safety measure

# Newdata for the prediction
nd <- data.frame(dummy_value = 42)

# -------------------------------------------------------------------
# Calculate true and modelled PDF
# -------------------------------------------------------------------

# Our ground truth
expect_silent(true_pdf <- dnorm(binmid, truth$mu, truth$sd))

# This uses 'elementwise = FALSE' by default.
expect_true(is.numeric(true_pdf) && length(true_pdf) == length(binmid),
               info = "Testing if 'true_pdf' has been calculated properly")
# Testing against the same result using elementwise = FALSE
expect_warning(mod_pdf  <- transitreg:::transitreg_predict(mod, newdata = nd, y = bm, type = "pdf"),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting PDF based on the transitreg model (elementwise = NULL)")
expect_warning(mod_pdf2 <- transitreg:::transitreg_predict(mod, newdata = nd, y = bm, type = "pdf", elementwise = FALSE),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting PDF based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_pdf, mod_pdf2, info = "Comparing result (PDF) with elementwise NULL/FALSE")
rm(mod_pdf2)

# Checking returned object
expect_inherits(mod_pdf, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_pdf), c(1L, length(binmid)),
              info = "Dimension of returned matrix")
expect_true(all(mod_pdf >= 0),
              info = "Quick check on range of PDF")
expect_identical(dimnames(mod_pdf), list(NULL, sprintf("d_%d", seq_along(binmid) - 1)),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true CDF and the predicted CDF should be small (elementwise)
expect_true(all(abs(mod_pdf - true_pdf) < 1e-3),
              info = "Checking elementwise difference between predicted PDF and ground truth.")

### matplot(x = bm, y = cbind(true = true_pdf, transitreg = as.vector(mod_pdf)),
###         col = c("gray", "tomato"), lwd = c(3, 2), type = "l",
###         xlab = "response 'y'", ylab = "Density", main = "Comparison")
### legend("topleft", c("true", "transitreg"), col = c("gray", "tomato"), lwd = c(3, 2), lty = 1:2)
rm(true_pdf, mod_pdf)




# -------------------------------------------------------------------
# Calculate true and modelled CDF
# -------------------------------------------------------------------

# Our ground truth
expect_silent(true_cdf <- pnorm(binmid, truth$mu, truth$sd))

# This uses 'elementwise = FALSE' by default.
expect_true(is.numeric(true_cdf) && length(true_cdf) == length(binmid),
               info = "Testing if 'true_cdf' has been calculated properly")
# Testing against the same result using elementwise = FALSE
expect_warning(mod_cdf  <- transitreg:::transitreg_predict(mod, newdata = nd, y = bm, type = "cdf"),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting CDF based on the transitreg model (elementwise = NULL)")
expect_warning(mod_cdf2 <- transitreg:::transitreg_predict(mod, newdata = nd, y = bm, type = "cdf", elementwise = FALSE),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting CDF based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_cdf, mod_cdf2, info = "Comparing result (CDF) with elementwise NULL/FALSE")
rm(mod_cdf2)

# Checking returned object
expect_inherits(mod_cdf, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_cdf), c(1L, length(binmid)),
              info = "Dimension of returned matrix")
expect_true(all(mod_cdf >= 0 & mod_cdf <= 1),
              info = "Quick check on range of CDF")
expect_identical(dimnames(mod_cdf), list(NULL, sprintf("p_%d", seq_along(binmid) - 1)),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true CDF and the predicted CDF should be small (elementwise)
expect_true(all(abs(mod_cdf - true_cdf) < 1e-1),
              info = "Checking elementwise difference between predicted CDF and ground truth.")

### matplot(x = bm, y = cbind(true = true_cdf, transitreg = as.vector(mod_cdf)),
###         col = c("gray", "tomato"), lwd = c(3, 2), type = "l",
###         xlab = "response 'y'", ylab = "Density", main = "Comparison")
### legend("topleft", c("true", "transitreg"), col = c("gray", "tomato"), lwd = c(3, 2), lty = 1:2)
rm(true_cdf, mod_cdf)




# -------------------------------------------------------------------
# Calculate true and modelled Quantile
# -------------------------------------------------------------------

# Probability vector
p <- seq(0.01, 0.99, by = 0.01)

# Our ground truth
expect_silent(true_q <- qnorm(p, truth$mu, truth$sd))

# This uses 'elementwise = FALSE' by default.
expect_true(is.numeric(true_q) && length(true_q) == length(p),
               info = "Testing if 'true_q' has been calculated properly")
# Testing against the same result using elementwise = FALSE
expect_warning(mod_q  <- transitreg:::transitreg_predict(mod, newdata = nd, prob = p, type = "quantile"),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting Quantile based on the transitreg model (elementwise = NULL)")
expect_warning(mod_q2 <- transitreg:::transitreg_predict(mod, newdata = nd, prob = p, type = "quantile", elementwise = FALSE),
               pattern = "Currently assuming 'discrete = FALSE'",
               info = "Predicting Quantile based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_q, mod_q2, info = "Comparing result (Quantile) with elementwise NULL/FALSE")
rm(mod_q2)

# Checking returned object
expect_inherits(mod_q, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_q), c(1L, length(binmid)),
              info = "Dimension of returned matrix")
expect_identical(dimnames(mod_q), list(NULL, trimws(paste0(format(p * 100, 0), "%"))),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true Quantile and the predicted Quantile should be small (elementwise)
expect_true(all(abs(mod_q - true_q) < 1.),
              info = "Checking elementwise difference between predicted Quantile and ground truth.")

### matplot(x = cbind(true = true_q, transitreg = as.vector(mod_q)), y = p,
###         col = c("gray", "tomato"), lwd = c(3, 2), type = "l",
###         xlab = "response 'y'", ylab = "Density", main = "Comparison")
### legend("topleft", c("true", "transitreg"), col = c("gray", "tomato"), lwd = c(3, 2), lty = 1:2)
rm(true_q, mod_q)






