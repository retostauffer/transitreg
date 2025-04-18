# -------------------------------------------------------------------
# Basic tests of transitreg models
# -------------------------------------------------------------------

if (interactive()) { library("tinytest"); library("transitreg") }


set.seed(6020)
data <- data.frame(y = rpois(100, 3.5))

# -------------------------------------------------------------------
# Estimating a simple model
# -------------------------------------------------------------------
expect_silent(m <- transitreg(y ~ s(theta, k = 4), data = data),
              info = "Estmating simple model for testing S3methods on transitreg object")
expect_inherits(m, "transitreg", info = "Ensure that the return is a transitreg object")


# -------------------------------------------------------------------
# [] subsetting
# -------------------------------------------------------------------
expect_silent(d <- m[1], info = "Subsetting method, returns a Transition object (length 1)")
expect_inherits(d, "Transition", info = "Checking return class")
expect_identical(length(d), 1L, info = "Checking return length")

expect_silent(d <- m[5:6], info = "Subsetting method, returns a Transition object (length 2)")
expect_inherits(d, "Transition", info = "Checking return class")
expect_identical(length(d), 2L, info = "Checking return length")

idx <- sample(nrow(data), 12)
expect_silent(d <- m[idx], info = "Subsetting method, returns a Transition object (length 12)")
expect_inherits(d, "Transition", info = "Checking return class")
expect_identical(length(d), 12L, info = "Checking return length")
rm(d, idx)


# -------------------------------------------------------------------
# LogLik()
# -------------------------------------------------------------------
expect_silent(ll <- logLik(m),      info = "Extracting logLik from model")
expect_inherits(ll, "logLik",       info = "Checking logLik return class")
expect_equal(sum(log(m$probs$pdf)), as.numeric(ll), info = "Checking logLik value")

# Calculating logLik on newdata (same data) should result in the same log likelihood
expect_identical(logLik(m, newdata = data), ll, info = "Calculated logLik on training data")

# Calculating logLik on a data different to training data (head of training data)
expect_silent(lln <- logLik(m, newdata = head(data)),    info = "Calculating logLik for 'new data'")
expect_identical(sum(log(m$probs$pdf[1:6])), as.numeric(lln), info = "Compare to manual calculation")

expect_identical(logLik(m, newdata = data), ll, info = "Calculated logLik on training data")
rm(ll, lln)

# -------------------------------------------------------------------
# model.frame()
# -------------------------------------------------------------------
expect_silent(mf <- model.frame(m), info = "Extracting model.frame")
expect_equivalent(mf, data,         info = "Checking content of model.frame (data[, important vars])")


# -------------------------------------------------------------------
# formula()
# -------------------------------------------------------------------
expect_equivalent(formula(m), Y ~ s(theta, k = 4), info = "Testing formula")


# -------------------------------------------------------------------
# summary(), print()
# -------------------------------------------------------------------
s <- summary(m)
expect_inherits(summary(m), "summary.gam", info = "summary returns summary of the internal model")
expect_stdout(print(m),                    info = "print produces output, not yet tested in more detail")


# -------------------------------------------------------------------
# prodist(), newresponse()
# -------------------------------------------------------------------
expect_silent(pd <- prodist(m),            info = "calling prodist")
expect_inherits(pd, "Transition",          info = "Testing prodist return class")
expect_identical(length(pd), nrow(data),   info = "Testing prodist return length")
expect_identical(pd[1], m[1],              info = "[] is a shortcut for prodict()[]")

expect_silent(nr <- newresponse(m),        info = "Calling newresponse()")
expect_equivalent(nr, data,                info = "newresponse equivalent to 'data' in this case")

expect_silent(nr2 <- newresponse(m, newdata = head(data)), info = "Testing newresponse with newdata")
expect_equivalent(nr2, head(data),         info = "Testing return of newresponse with newdata")

# Provide newdata, but does not contain the response variable -> error
tmp <- head(data); names(tmp) <- "foo"
expect_error(newresponse(m, newdata = tmp),
             pattern = "Response variable 'y' missing in newdata!",
             info = "Testing newdata with missing response variable")
rm(tmp)


# -------------------------------------------------------------------
# Prediction methods
# -------------------------------------------------------------------


# Predict mode
expect_silent(p <- predict(m, type = "mode"),      info = "Predicting type = 'mode'")
expect_inherits(p, "numeric",                      info = "mode prediction return class")
expect_true(all(p == 3.5),                         info = "mode prediction return values")
rm(p)


# Predict mean
expect_silent(p <- predict(m, type = "mean"),      info = "Predicting type = 'mean'")
expect_inherits(p, "numeric",                      info = "mean prediction return class")
expect_true(all(p == mean(m[1])),                  info = "mean prediction return values")
expect_true(all.equal(p, rep(3.706188, length(p)), tol = 1e-7),
            info = "mean prediction return values")
rm(p)


# Predict quantile
expect_silent(p0 <- predict(m, type = "mode"))
expect_silent(p1 <- predict(m, type = "quantile"),
            info = "Predicting type = 'quantile'")
expect_inherits(p1, "numeric",                info = "quantile prediction return class")
expect_silent(p2 <- predict(m, type = "quantile", prob = 0.5),
            info = "Predicting type = 'quantile' prob = 0.5")
expect_inherits(p2, "numeric",                info = "quantile prediction return class")

# Compare result with median
expect_equal(p0, p1, info = "Testing default behavior of predict type quantile.")
expect_equal(p1, p2, info = "Testing default behavior of predict type quantile.")
rm(p0, p1, p2)


# Predict quantile
expect_silent(p <- predict(m, type = "quantile", prob = -0.0001),
              info = "Quantile prediction for prob < 0")
expect_true(all(is.na(p)), info = "Return vor negative quantile")
expect_silent(p <- predict(m, type = "quantile", prob = +1.0001),
              info = "Quantile prediction for prob > 1")
expect_true(all(is.na(p)), info = "Return vor quantile > 1")

# Quantile 0.0 must be mid of first bin (i.e., 0.5)
expect_silent(p <- predict(m, type = "quantile", prob = 0),
              info = "Quantile prediction for prob == 0")
expect_true(all(p == 0.5), info = "Return value quantile 0")
mid <- transitreg:::get_mids(m)[1]
expect_true(all(p == mid), info = "Return value quantile 0")
rm(mid, p)


# Quantile 1.0 must be end of range (last break)
expect_silent(p <- predict(m, type = "quantile", prob = 1),
              info = "Quantile prediction for prob == 1")
expect_true(all(p == m$bins), info = "Return value quantile 1")
bmax <- max(transitreg:::get_breaks(m))
expect_true(all(p == bmax), info = "Return value quantile 1")
rm(bmax, p)


# Quantile 1.0 must be end of range (last break)
expect_silent(p <- predict(m, type = "quantile", prob = c(0, 0.5, 1)),
              info = "Quantile prediction for prob = c(0, 0.5, 1)")
expect_inherits(p, "matrix",  info = "Checking return class")
expect_identical(dim(p), c(nrow(data), 3L), info = "Dimension of matrix returned")
expect_true(is.numeric(p),                  info = "Matrix must be numeric")
expect_identical(colnames(p), paste0(c(0, 50, 100), "%"), info = "Matrix column names")

expect_true(all(p[, 1L] == 0.5),    info = "Checking return value for lowest quantile")
expect_true(all(p[, 2L] == 3.5),    info = "Checking return value for median (q50)")
expect_true(all(p[, 3L] == m$bins), info = "Checking return value for highest quantile")
rm(p)


# CDF
mids <- transitreg:::get_mids(m)
expect_silent(p <- predict(m, type = "cdf", y = mids),
              info = "CDF prediction on transitreg object")
expect_inherits(p, "matrix",        info = "Return class")
expect_true(is.numeric(p),          info = "Return must be numeric")
expect_identical(dim(p), c(nrow(data), length(mids)), info = "Return dimension check")
tmp <- cdf(m[1], x = mids)
tmp_mat <- matrix(tmp, ncol = length(mids), nrow = nrow(data), byrow = T)
expect_equivalent(p, tmp_mat,       info = "Checking matrix values")
rm(mids, p, tmp, tmp_mat)


# PDF
mids <- transitreg:::get_mids(m)
expect_silent(p <- predict(m, type = "pdf", y = mids),
              info = "PDF prediction on transitreg object")
expect_inherits(p, "matrix",        info = "Return class")
expect_true(is.numeric(p),          info = "Return must be numeric")
expect_identical(dim(p), c(nrow(data), length(mids)), info = "Return dimension check")
tmp <- pdf(m[1], x = mids)
tmp_mat <- matrix(tmp, ncol = length(mids), nrow = nrow(data), byrow = T)
expect_equivalent(p, tmp_mat,       info = "Checking matrix values")
rm(mids, p, tmp, tmp_mat)


# TP
expect_silent(p <- predict(m, type = "tp"),
              info = "PDF prediction on transitreg object")
expect_inherits(p, "matrix",        info = "Return class")
expect_true(is.numeric(p),          info = "Return must be numeric")
expect_identical(dim(p), c(nrow(data), m$bins), info = "Return dimension check")
tmp <- as.vector(as.matrix(m[1]))
tmp_mat <- matrix(tmp, ncol = m$bins, nrow = nrow(data), byrow = T)
expect_equivalent(p, tmp_mat,       info = "Checking matrix values")
rm(p, tmp, tmp_mat)

