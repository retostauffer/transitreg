# -------------------------------------------------------------------
# Basic tests of transitreg models
# -------------------------------------------------------------------


if (interactive()) { library("tinytest"); library("transitreg") }


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
# The function get_braks will return 'pseudo breaks' centered arou nd 0, ..., m1$bins - 1L
tmp <- seq.int(0L, m1$bins)
expect_equal(transitreg:::get_breaks(m1), tmp,
               info = "Checking pseudo-breaks for count data response")
# Mid points used get_breaks() and calculates mid points. In this model this
# should be equal to 0, 1, 2, ..., m2$bins - 1L
expect_equal(transitreg:::get_mids(m1), seq(0, m1$bins - 1L) + 0.5,
               info = "Checking return of get_mids() helper function")
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
expect_true(is.null(m2$breaks),
               info = "Count data, no breaks set")
# Breaks are stored on the object, get_breaks should return the exact same.
expect_silent(tmp <- transitreg:::get_breaks(m2),
               info = "Extracting breaks")
expect_inherits(tmp, "integer", info = "Integer breaks (count data model)")
expect_identical(tmp, 0:9,      info = "Checking break values")
# Checking get_mids helper function
tmp_mid <- (tmp[-1L] + tmp[-length(tmp)]) / 2.0
expect_equal(transitreg:::get_mids(m2), tmp_mid,
               info = "Checking return of get_mids() helper function")
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
               info = "Test that the breaks are our user-defined breaks")
# Checking helper function 'get_breaks()'
expect_equal(transitreg:::get_breaks(m3), bk,
               info = "Testing return of get_breaks() helper function")
# Checking get_mids helper function
tmp_mid <- (bk[-1L] + bk[-length(bk)]) / 2.0
expect_equal(transitreg:::get_mids(m3), tmp_mid,
               info = "Checking return of get_mids() helper function")
rm(m3)



