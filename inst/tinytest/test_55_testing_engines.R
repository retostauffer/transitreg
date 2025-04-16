# -------------------------------------------------------------------
# Testing 'transitreg' for discrete response (count data).
# -------------------------------------------------------------------


if (interactive()) { library("tinytest"); library("transitreg") }

data(cars)
cars <- na.omit(cars)

# -------------------------------------------------------------------
# Engine = bam
# -------------------------------------------------------------------
# Default behavior
bam1 <- transitreg(dist ~ s(theta, speed, k = 4), breaks = 10, data = cars)
expect_inherits(bam1$model, "bam",  info = "Check that engine = bam is used")

# Testing with engine = gam
bam2 <- transitreg(dist ~ s(theta, speed, k = 4), breaks = 10, data = cars, engine = "bam")
expect_inherits(bam2$model, "bam",  info = "Check that engine = bam is used")
expect_identical(logLik(bam1), logLik(bam2),      info = "Checking that 'bam' is our default engine")
expect_identical(coef(bam1), coef(bam2),          info = "Checking coefficients")

# -------------------------------------------------------------------
# Engine = gam
# -------------------------------------------------------------------
gam <- transitreg(dist ~ s(theta, speed, k = 4), breaks = 10, data = cars, engine = "gam")
expect_inherits(gam$model, "gam",  info = "Check that engine = gam is used")
expect_equal(as.numeric(logLik(gam)), as.numeric(logLik(bam1)), tol = 1e-4,
             info = "Testing logLik bam vs. gam (must be equal, not identical)")

# -------------------------------------------------------------------
# Engine = nnet
# -------------------------------------------------------------------
nnet <- transitreg(dist ~ theta + speed, size = 100, trace = F, breaks = 10, data = cars, engine = "nnet")
expect_inherits(nnet$model, "nnet",  info = "Check that engine = nnet is used")





