# -------------------------------------------------------------------
# Basic tests of transitreg models
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# -------------------------------------------------------------------
# Simulating data; 100 random values from a Poisson distribution.
# -------------------------------------------------------------------
set.seed(6020)
ndp <- data.frame(y = rpois(100, 3))

# Estimating a simple model without 'breaks' specification.
# 'transitreg()' assumes that the response is count data, and the
# number of breaks is simply the max response. Breaks are not stored
# on the model object, this indicates that there were no user-defined
# breaks.
expect_warning(m1 <- transitreg(y ~ theta, data = ndp),
               pattern = "no smooths",
               info = "Estimating transitreg model without breaks")
expect_inherits(m1, "transitreg",
               info = "Checking return class")
expect_true("bins" %in% names(m1),
               info = "Element 'bins' must be stored in the model")
expect_identical(m1$bins, as.integer(max(ndp$y) * 3),
               info = "Checking number of bins (max response value * 3 as n <= 10)")
expect_identical(m1$breaks, NULL,
               info = "No user-defined breaks, thus 'breaks' must be NULL")
rm(m1)


# Let's do the same, but this time we are providing an integer
# for 'breaks'. This suggests that our response values are not counts but
# continuous. The 'bins' will span the entire range of the response.
nbk <- 10L
expect_warning(m2 <- transitreg(y ~ theta, data = ndp, breaks = nbk),
               pattern = "no smooths",
               info = "Estimating transitreg model with breaks = 10")
expect_inherits(m2, "transitreg",
               info = "Checking return class")
expect_true("bins" %in% names(m2),
               info = "Element 'bins' must be stored in the model")
expect_identical(m2$bins, nbk - 1L,
               info = "Checking number of bins")
expect_equal(m2$breaks, seq(-0.8, 8.8, length.out = m2$bins + 1),
               info = "Breaks are alculated automatically via make_breaks()")
rm(m2, nbk)


# Alternatively, we provide a numeric vector with pre-specified breaks.
bk <- seq(-1.5, 20.5, by = 1)
expect_warning(m3 <- transitreg(y ~ theta, data = ndp, breaks = bk),
               pattern = "no smooths",
               info = "Estimating transitreg model with breaks = vector")
expect_inherits(m3, "transitreg",
               info = "Checking return class")
expect_identical(m3$bins, length(bk) - 1L,
               info = "Checking number of bins")
expect_equal(m3$breaks, bk,
               info = "Test that the breaks are our user-defined breaks.")
rm(m3)



