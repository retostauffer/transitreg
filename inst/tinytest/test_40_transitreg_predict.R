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
expect_identical(mod$breaks, breaks, info = "Testing if the breaks have been stored correctly")
# TODO(R): Currently it stores the bins on breaks _and_ bins?

# Calculating center of the bins
binmid <- (head(mod$breaks, -1) + tail(mod$breaks, -1)) / 2
expect_identical(length(binmid), length(mod$breaks) - 1L) # <- safety measure
expect_identical(length(binmid), mod$bins) # <- safety measure

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
expect_silent(mod_pdf  <- transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "pdf"),
               info = "Predicting PDF based on the transitreg model (elementwise = NULL)")
expect_silent(mod_pdf2 <- transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "pdf", elementwise = FALSE),
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

### matplot(x = binmid, y = cbind(true = true_pdf, transitreg = as.vector(mod_pdf)),
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
expect_silent(mod_cdf  <- transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "cdf"),
               info = "Predicting CDF based on the transitreg model (elementwise = NULL)")
expect_silent(mod_cdf2 <- transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "cdf", elementwise = FALSE),
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

### matplot(x = binmid, y = cbind(true = true_cdf, transitreg = as.vector(mod_cdf)),
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
expect_silent(mod_q  <- transitreg:::transitreg_predict(mod, newdata = nd, prob = p, type = "quantile"),
               info = "Predicting Quantile based on the transitreg model (elementwise = NULL)")
expect_silent(mod_q2 <- transitreg:::transitreg_predict(mod, newdata = nd, prob = p, type = "quantile", elementwise = FALSE),
               info = "Predicting Quantile based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_q, mod_q2, info = "Comparing result (Quantile) with elementwise NULL/FALSE")
rm(mod_q2)

# Checking returned object
expect_inherits(mod_q, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_q), c(1L, length(p)),
              info = "Dimension of returned matrix")
expect_identical(dimnames(mod_q), list(NULL, trimws(paste0(format(p * 100, 0), "%"))),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true Quantile and the predicted Quantile should be small (elementwise)
expect_true(all(abs(mod_q - true_q) < 2.),
              info = "Checking elementwise difference between predicted Quantile and ground truth.")

### matplot(x = cbind(true = true_q, transitreg = as.vector(mod_q)), y = p,
###         col = c("gray", "tomato"), lwd = c(3, 2), type = "l",
###         xlab = "response 'y'", ylab = "Density", main = "Comparison")
### legend("topleft", c("true", "transitreg"), col = c("gray", "tomato"), lwd = c(3, 2), lty = 1:2)
rm(true_q, mod_q)


# Testing the prediction of transition probabilities.
expect_silent(tp <- transitreg:::transitreg_predict(mod, newdata = nd, type = "tp"),
               info = "Predicting transition probabilities.")
expect_inherits(tp, "matrix",
               info = "Return is a matrix array.")
expect_true(is.numeric(tp),
               info = "Matrix must be numeric")
expect_identical(dim(tp), c(nrow(nd), mod$bins),
               info = "Testing matrix dimension")
expect_true(all(grepl("^tp_[0-9]+$", colnames(tp))),
               info = "Checking if matrix column name are as expected.")
expect_true(is.null(rownames(tp)),
               info = "Checking if matrix row name are as expected.")

# Testing the values inside; we only predicted the tps for one data point
# (nrow(nd) == 1), thus we simply convert it into a vector.
tp <- as.vector(tp)
expect_true(all(diff(tp) < 0),
               info = "Transition probabilities must be monotonically decreasing")

# Calculate PDF, CDF
expect_silent(tp_pdf <- as.vector(transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "pdf")),
               info = "Predicting PDF for a test")
expect_silent(tp_cdf <- as.vector(transitreg:::transitreg_predict(mod, newdata = nd, y = binmid, type = "cdf")),
               info = "Predicting CDF for a test")
expect_equal(convert_tp(tp, from = "tp", to = "cdf"), tp_cdf,
               info = "Testing conversion from TP to CDF")
expect_equal(convert_tp(tp, from = "tp", to = "pdf", width = diff(mod$breaks)), tp_pdf,
               info = "Testing conversion from TP to PDF")



# Predicting TP for nrow(nd) > 1 to test if the result is correct.
nd2 <- data.frame(dummy_value = 1:2)
expect_silent(tp2 <- transitreg:::transitreg_predict(mod, newdata = nd2, type = "tp"),
               info = "Predicting transition probabilities for two data points")
expect_inherits(tp2, "matrix",
               info = "Return is a matrix array")
expect_true(is.numeric(tp2),
               info = "Matrix must be numeric")
expect_identical(dim(tp2), c(nrow(nd2), mod$bins),
               info = "Testing matrix dimension")
expect_true(all(grepl("^tp_[0-9]+$", colnames(tp2))),
               info = "Checking if matrix column name are as expected.")
expect_true(is.null(rownames(tp2)),
               info = "Checking if matrix row name are as expected.")
expect_identical(tp2[1, ], tp2[2, ],
               info = "As we have no covaraites, both should be identical")
rm(tp2, nd2)

