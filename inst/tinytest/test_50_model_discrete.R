# -------------------------------------------------------------------
# Testing 'transitreg' for discrete response (count data).
# TODO(R): Remove useC and comparison against the R version
#          once we removed that.
# -------------------------------------------------------------------


suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))


# ===================================================================
# Simple intercept=only model; default options
# ===================================================================

suppressPackageStartupMessages(library("gamlss.data"))
data("CD4", package = "gamlss.data")

# -------------------------------------------------------------------
# Estimating simple model, testing return object.
# -------------------------------------------------------------------
m1   <- transitreg(cd4 ~ theta, data = CD4)

expect_inherits(m1, "transitreg", info = "Returnclass")
expect_true(is.list(m1))
expected <- c("new_formula", "model", "response", "model.frame",
              "maxcounts", "theta_vars", "factor", "probs", "bins", "ym", "yc_tab", "breaks")
expect_true(all(names(m1) %in% expected),
            info = "Checking if all expected elements are there")
# transitreg() auto-guesses ym/bins if the response looks like count data.
# check if it did what it should do.
m1_ym   <- seq.int(0, max(CD4$cd4))
m1_bins <- seq.int(0, max(CD4$cd4) + 1) - 0.5
expect_equal(m1_ym,   m1$ym,   info = "Testing if auto-guessed ym is correct")
expect_equal(m1_bins, m1$bins, info = "Testing if auto-guessed bins is correct")
rm(m1_ym, m1_bins)


# -------------------------------------------
# checking element $new_formula
# For comparison and convenience I am also testing
# the S3 method formula.transitreg here and compare $new_formula
# to the model formula.
# -------------------------------------------
expect_silent(f <- formula(m1), info = "Testing S3 method formula")
expect_inherits(f, "formula", info = "formula return class")
expect_identical(format(f), "Y ~ theta", info = "formula return value (text)")
expect_identical(length(f), 3L, info = "formula return value (Y ~ theta = 3)")

# Checking $new_formula, testing against model formula
expect_inherits(m1$new_formula, "formula", info = "formula return class")
expect_identical(format(f), format(m1$new_formula),
                 info = "Comparing model formula against $new_formula")
rm(f)

# -------------------------------------------
# Checking elements: $model, $model.frame,
#                    $response, $maxcounts, $theta_vars,
#                    $factor, $probs
# Note: Not checking for correctness, just a
# rudimentary check if everything is where it must be.
# -------------------------------------------
expect_inherits(m1$model, "bam", info = "Checking $model class")
expect_identical(dim(model.frame(m1$model)),
                 c(as.integer(sum(CD4$cd4)) + nrow(CD4), 2L),
                 info = "Dimension of model.frame($model)")

expect_inherits(m1$model.frame, "data.frame")
expect_identical(dim(m1$model.frame), c(nrow(CD4), 1L),
                 info = "dimension of $model.frame element")

expect_identical(m1$response, "cd4", info = "Value of $response (name of response variable)")

expect_identical(m1$maxcounts, max(CD4$cd4), info = "Value of $maxcounts")

expect_true(is.character(m1$theta_vars) && length(m1$theta_vars) == 0,
            info = "Value and length of $theta_vars")

expect_true(isFALSE(m1$factor), info = "Boolean value $factor")

expect_inherits(m1$probs, "data.frame", info = "Class of $probs")
expect_identical(dim(m1$probs), c(nrow(CD4), 2L),       info = "Dimension of $probs")
expect_identical(names(m1$probs), c("pdf", "cdf"),      info = "Names of $probs")
expect_true(all(m1$probs$cdf >= 0 & m1$probs$cdf <= 1), info = "All $probs$cdf in [0, 1]")
expect_true(all(m1$probs$pdf >= 0),                     info = "All $probs$pdf in [0, Inf]")


# -------------------------------------------
# Testing for correctness (against hardcoded values)
# -------------------------------------------

# (1) Compare coefficients
expect_silent(tmp <- coef(m1), info = "Testing S3 method 'coef'")
expect_true(is.double(tmp), info = "coef return type")
expect_identical(names(tmp), c("(Intercept)", "theta"),
                  info = "coef return names")
expect_equivalent(tmp, c(6.5583511036, -0.0004586303), tolerance = 1e-9,
                  info = "coef return values")
rm(tmp)

# (2) Compare overall logLik
expect_silent(tmp <- logLik(m1), info = "Testing S3 method 'logLik'")
expect_inherits(tmp, "logLik", info = "logLik return class")
expect_true(is.double(tmp), info = "logLik return type")
expect_identical(length(tmp), 1L, info = "logLik return length")
expect_equivalent(as.numeric(tmp), -4447.884, tolerance = 1e-2, info = "logLik return value")
rm(tmp)

# (3) Sum over all CDF and PDF values; that is just a quick'n'
#     dirty check if the content of $prob resembles what
#     we expect for this very specific model.
expect_equivalent(sum(m1$probs), 306.1717, tolerance = 1e-3, info = "sum($probs) as quick check")



# -------------------------------------------------------------------
# Estimate the same model, specifying all default options to
# test that they remain unchanged.
# -------------------------------------------------------------------
m2   <- transitreg(cd4 ~ theta, data = CD4,
           engine = "bam", scale.x = FALSE, breaks = NULL,
           model = TRUE, ncores = NULL, verbose = FALSE)

# Comparing m1 against m2 (rough check)
expect_equal(logLik(m1), logLik(m2), info = "Checking results w/ default args")
expect_equal(coef(m1), coef(m2),     info = "Checking results w/ default args")
rm(m2)



