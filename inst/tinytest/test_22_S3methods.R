# -------------------------------------------------------------------
# Basic tests of transitreg models
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))


set.seed(6020)
data <- data.frame(y = rpois(100, 3))

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












