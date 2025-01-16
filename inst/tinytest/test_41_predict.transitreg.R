# -------------------------------------------------------------------
# Testing prediction method
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
w <- diff(range(data$y))
breaks <- seq(min(data$y) - 0.2 * w, max(data$y) + 0.2 * w, length.out = 201)

expect_silent(mod <- transitreg(y ~ s(theta), data = data, breaks = breaks),
              info = "Estimating transitional regression model")
expect_inherits(mod, "transitreg", info = "Checking transitreg() return class")
expect_identical(mod$breaks, breaks, info = "Testing if the breaks have been stored correctly")
expect_identical(mod$bins, length(breaks) - 1L, info = "Testing number of bins")


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
expect_silent(mod_pdf  <- predict(mod, newdata = nd, y = binmid, type = "pdf"),
               info = "Predicting PDF based on the transitreg model (elementwise = NULL)")
expect_silent(mod_pdf2 <- predict(mod, newdata = nd, y = binmid, type = "pdf", elementwise = FALSE),
               info = "Predicting PDF based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_pdf, mod_pdf2, info = "Comparing result (PDF) with elementwise NULL/FALSE")
rm(mod_pdf2)

# Checking returned object
expect_inherits(mod_pdf, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_pdf), c(nrow(nd), length(binmid)),
              info = "Dimension of returned matrix")
expect_true(all(mod_pdf >= 0),
              info = "Quick check on range of PDF")
expect_identical(dimnames(mod_pdf), list(NULL, sprintf("d_%d", seq_along(binmid) - 1)),
              info = "Checking dimnames of returned matrix")


# If the model and the prediction method work as expected, the difference
# between the true CDF and the predicted CDF should be small (elementwise)
expect_true(all(abs(mod_pdf - true_pdf) < 1e-3),
              info = "Checking elementwise difference between predicted PDF and ground truth.")
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
expect_silent(mod_cdf  <- predict(mod, newdata = nd, y = binmid, type = "cdf"),
               info = "Predicting CDF based on the transitreg model (elementwise = NULL)")
expect_silent(mod_cdf2 <- predict(mod, newdata = nd, y = binmid, type = "cdf", elementwise = FALSE),
               info = "Predicting CDF based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_cdf, mod_cdf2, info = "Comparing result (CDF) with elementwise NULL/FALSE")
rm(mod_cdf2)

# Checking returned object
expect_inherits(mod_cdf, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_cdf), c(nrow(nd), length(binmid)),
              info = "Dimension of returned matrix")
expect_true(all(mod_cdf >= 0 & mod_cdf <= 1),
              info = "Quick check on range of CDF")
expect_identical(dimnames(mod_cdf), list(NULL, sprintf("p_%d", seq_along(binmid) - 1)),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true CDF and the predicted CDF should be small (elementwise)
expect_true(all(abs(mod_cdf - true_cdf) < 1e-2),
              info = "Checking elementwise difference between predicted CDF and ground truth.")
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
expect_silent(mod_q  <- predict(mod, newdata = nd, prob = p, type = "quantile"),
               info = "Predicting Quantile based on the transitreg model (elementwise = NULL)")
expect_silent(mod_q2 <- predict(mod, newdata = nd, prob = p, type = "quantile", elementwise = FALSE),
               info = "Predicting Quantile based on the transitreg model (elementwise = FALSE)")
expect_identical(mod_q, mod_q2, info = "Comparing result (Quantile) with elementwise NULL/FALSE")
rm(mod_q2)

# Checking returned object
expect_inherits(mod_q, "matrix",
              info = "Testing return class")
expect_identical(dim(mod_q), c(nrow(nd), length(p)),
              info = "Dimension of returned matrix")
expect_identical(dimnames(mod_q), list(NULL, trimws(paste0(format(p * 100, 0), "%"))),
              info = "Checking dimnames of returned matrix")

# If the model and the prediction method work as expected, the difference
# between the true Quantile and the predicted Quantile should be small (elementwise)
expect_true(all(abs(mod_q - true_q) < 2.),
              info = "Checking elementwise difference between predicted Quantile and ground truth.")
rm(true_q, mod_q)


# Testing the prediction of transition probabilities.
expect_silent(tp <- predict(mod, newdata = nd, type = "tp"),
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
expect_silent(tp_pdf <- predict(mod, newdata = nd, y = binmid, type = "pdf"),
               info = "Predicting PDF for a test")
expect_silent(tp_cdf <- predict(mod, newdata = nd, y = binmid, type = "cdf"),
               info = "Predicting CDF for a test")
expect_equal(convert_tp(tp, from = "tp", to = "cdf"), unname(tp_cdf[1, ]),
               info = "Testing conversion from TP to CDF")
expect_equal(convert_tp(tp, from = "tp", to = "pdf", width = diff(mod$breaks)), unname(tp_pdf[1, ]),
               info = "Testing conversion from TP to PDF")



# -------------------------------------------------------------------
# Testing with 'nrow(newdata) == 2L' to check if the return is as expected
# -------------------------------------------------------------------
nd2 <- data.frame(dummy_value = 1:2)

# Transition probability
expect_silent(tp2 <- predict(mod, newdata = nd2, type = "tp"),
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
rm(tp2)

# CDF
expect_silent(cdf2 <- predict(mod, newdata = nd2, binmid, type = "cdf"),
               info = "Predicting transition probabilities for two data points")
expect_inherits(cdf2, "matrix",
               info = "Return is a matrix array")
expect_true(is.numeric(cdf2),
               info = "Matrix must be numeric")
expect_identical(dim(cdf2), c(nrow(nd2), mod$bins),
               info = "Testing matrix dimension")
expect_true(all(grepl("^p_[0-9]+$", colnames(cdf2))),
               info = "Checking if matrix column name are as expected.")
expect_true(is.null(rownames(cdf2)),
               info = "Checking if matrix row name are as expected.")
expect_identical(cdf2[1, ], cdf2[2, ],
               info = "As we have no covaraites, both should be identical")
rm(cdf2)

# PDF
expect_silent(pdf2 <- predict(mod, newdata = nd2, binmid, type = "pdf"),
               info = "Predicting transition probabilities for two data points")
expect_inherits(pdf2, "matrix",
               info = "Return is a matrix array")
expect_true(is.numeric(pdf2),
               info = "Matrix must be numeric")
expect_identical(dim(pdf2), c(nrow(nd2), mod$bins),
               info = "Testing matrix dimension")
expect_true(all(grepl("^d_[0-9]+$", colnames(pdf2))),
               info = "Checking if matrix column name are as expected.")
expect_true(is.null(rownames(pdf2)),
               info = "Checking if matrix row name are as expected.")
expect_identical(pdf2[1, ], pdf2[2, ],
               info = "As we have no covaraites, both should be identical")
rm(pdf2)

# Quantile
p <- seq(0.1, 0.9, by = 0.1)
expect_silent(q2 <- predict(mod, newdata = nd2, prob = p, type = "q"),
               info = "Predicting transition probabilities for two data points")
expect_inherits(q2, "matrix",
               info = "Return is a matrix array")
expect_true(is.numeric(q2),
               info = "Matrix must be numeric")
expect_identical(dim(q2), c(nrow(nd2), length(p)),
               info = "Testing matrix dimension")
expect_true(all(grepl("^[0-9]+\\%$", colnames(q2))),
               info = "Checking if matrix column name are as expected.")
expect_true(is.null(rownames(q2)),
               info = "Checking if matrix row name are as expected.")
expect_identical(q2[1, ], q2[2, ],
               info = "As we have no covaraites, both should be identical")
rm(q2)



# -------------------------------------------------------------------
# Testing functionality of 'elementwise'.
# By default it is NULL and auto-guesses what to do.
# If number of observations is equal to y/probs, elementwise is assumed.
# Else !elementwise -> resulting in a matrix. Testing this behavior.
# -------------------------------------------------------------------
nd2 <- data.frame(dummy = c(42, 42))
set.seed(123)
yy <- sample(data$y, 2)

# CDF
expect_silent(cdf1 <- predict(mod, newdata = nd2, y = yy, type = "cdf", elementwise = NULL),
              info = "CDF for 2 observations, elementwise = NULL")
expect_silent(cdf2 <- predict(mod, newdata = nd2, y = yy, type = "cdf", elementwise = TRUE),
              info = "CDF for 2 observations, elementwise = TRUE")
expect_silent(cdf3 <- predict(mod, newdata = nd2, y = yy, type = "cdf", elementwise = FALSE),
              info = "CDF for 2 observations, elementwise = FALSE")

expect_true(is.numeric(cdf1) && is.vector(cdf1) && length(cdf1) == nrow(nd2),
              info = "Result for CDF must be numeric vector of length nrow(nd)")
expect_identical(cdf1, cdf2,
              info = "Testing behavior of elementwise NULL/TRUE on CDF")
expect_inherits(cdf3, "matrix",
              info = "CDF with elementwise = FALSE should return a matrix")
expect_identical(dim(cdf3), c(2L, 2L),
              info = "Expecting matrix of dimension 2 x 2")
rm(cdf1, cdf2, cdf3)

# PDF
expect_silent(pdf1 <- predict(mod, newdata = nd2, y = yy, type = "pdf", elementwise = NULL),
              info = "CDF for 2 observations, elementwise = NULL")
expect_silent(pdf2 <- predict(mod, newdata = nd2, y = yy, type = "pdf", elementwise = TRUE),
              info = "CDF for 2 observations, elementwise = TRUE")
expect_silent(pdf3 <- predict(mod, newdata = nd2, y = yy, type = "pdf", elementwise = FALSE),
              info = "CDF for 2 observations, elementwise = FALSE")

expect_true(is.numeric(pdf1) && is.vector(pdf1) && length(pdf1) == nrow(nd2),
              info = "Result for CDF must be numeric vector of length nrow(nd)")
expect_identical(pdf1, pdf2,
              info = "Testing behavior of elementwise NULL/TRUE on CDF")
expect_inherits(pdf3, "matrix",
              info = "CDF with elementwise = FALSE should return a matrix")
expect_identical(dim(pdf3), c(2L, 2L),
              info = "Expecting matrix of dimension 2 x 2")
rm(pdf1, pdf2, pdf3)

# Quantile
p <- c(0.8, 0.9)
expect_silent(q1 <- predict(mod, newdata = nd2, prob = p, type = "quantile", elementwise = NULL),
              info = "CDF for 2 observations, elementwise = NULL")
expect_silent(q2 <- predict(mod, newdata = nd2, prob = p, type = "quantile", elementwise = TRUE),
              info = "CDF for 2 observations, elementwise = TRUE")
expect_silent(q3 <- predict(mod, newdata = nd2, prob = p, type = "quantile", elementwise = FALSE),
              info = "CDF for 2 observations, elementwise = FALSE")

expect_true(is.numeric(q1) && is.vector(q1) && length(q1) == nrow(nd2),
              info = "Result for CDF must be numeric vector of length nrow(nd)")
expect_identical(q1, q2,
              info = "Testing behavior of elementwise NULL/TRUE on CDF")
expect_inherits(q3, "matrix",
              info = "CDF with elementwise = FALSE should return a matrix")
expect_identical(dim(q3), c(2L, 2L),
              info = "Expecting matrix of dimension 2 x 2")
rm(q1, q2, q3)

# TP
expect_silent(tp1 <- predict(mod, newdata = nd2, type = "tp", elementwise = NULL),
              info = "CDF for 2 observations, elementwise = NULL")
expect_silent(tp2 <- predict(mod, newdata = nd2, type = "tp", elementwise = TRUE),
              info = "CDF for 2 observations, elementwise = TRUE")
expect_silent(tp3 <- predict(mod, newdata = nd2, type = "tp", elementwise = FALSE),
              info = "CDF for 2 observations, elementwise = FALSE")
expect_identical(tp1, tp2, info = "elementwise should have no effect on TP")
expect_identical(tp1, tp3, info = "elementwise should have no effect on TP")
rm(tp1, tp2, tp3)
