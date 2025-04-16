# -------------------------------------------------------------------
# Basic tests of transitreg models
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))


set.seed(6020)
data <- data.frame(y = rpois(100, 3))


# ---------------------------------------------------------------------
# Testing missing required args
# ---------------------------------------------------------------------
expect_silent(transitreg(y ~ s(theta, k = 3), data = data), info = "Test it works in general")
expect_error(transitreg(), info = "Missing 'formula' and 'data'")
expect_error(transitreg(formula = y ~ theta), info = "Missing 'data'")
expect_error(transitreg(data = data), info = "Missing 'formula'")

expect_error(transitreg(y ~ s(theta, k = 3), data = data, ncores = "foo"),
              info = "'ncores' not numeric")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, verbose = c(TRUE, FALSE)),
              info = "'verbose' not single logical")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, verbose = "TRUE"),
              info = "'verbose' not logical")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, breaks = numeric()),
              info = "'breaks', if set, must be numeric of length >= 1")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, breaks = -1.2),
              info = "'breaks', if single numeric, must be > 0")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, scale.x = c(TRUE, FALSE)),
              info = "'scale.x' not single logical")
expect_error(transitreg(y ~ s(theta, k = 3), data = data, scale.x = "TRUE"),
              info = "'scale.x' not logical")

bk_wrong <- seq(min(data$y), max(data$y) - 1L, length.out = 5)
expect_error(transitreg(y ~ s(theta, k = 3), data = data, breaks = bk_wrong),
              info = "'breaks' specified by user do not cover the range of the response")

# bam with intercept only
expect_error(transitreg(y ~ 1, engin = "bam", data = data),
             pattern = "Intercept only model.*with engine = \"bam\" not allowed\\.",
             info = "Intercept only model with bam must fail.")

# y is numeric data (not 'count data'), so breaks must be provided!
data_num <- data.frame(y = runif(100, 0, 30))
expect_error(transitreg(y ~ s(theta, k = 3), data = data_num),
             pattern = "response does not look like count data, in this case \\'breaks\\' must be specified",
             info = "Response numeric, but no breaks specified.")



# ---------------------------------------------------------------------
# Testing misspecified formulas
# ---------------------------------------------------------------------
demo_data <- data.frame(y = runif(50, 0, 3))

# Testing model estimation
expect_error(transitreg( ~ theta, data = demo_data, breaks = 10),
             pattern = "object 'theta' not found",
             info = "Invalid formula (no response)")
expect_error(transitreg(5 ~ theta, data = demo_data, breaks = 10),
             pattern = "invalid term in model formula",
             info = "Invalid formula (constant response)")
rm(demo_data)


# Testing 'thetaX', 'thetaY', ...
demo_data <- data.frame(y = sample(0:2, 30, replace = TRUE))
demo_bk   <- seq(-0.5, 3.5, by = 1)
expect_error(transitreg(y ~ theta0 + theta1 + theta2 + theta3, data = demo_data, breaks = demo_bk),
             pattern = "Formula contains 'theta3', but no observation falls into this bin",
             info = "Testing misspecified formula, theta3 but no observations in this bin.")
expect_error(transitreg(y ~ theta0 + theta1 + theta2, data = demo_data, breaks = demo_bk),
             pattern = "Formula contains 'theta0', .*, covering all bins observations fall into",
             info = "Testing overspecified formula")
rm(demo_data, demo_bk)




