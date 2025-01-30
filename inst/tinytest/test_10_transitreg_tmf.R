# -------------------------------------------------------------------
# Testing transitreg_tmf, which sets up the 'transitreg model frame'
# for the (internal) binary logit model.
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# -------------------------------------------------------------------
# Taking this small data.set for testing
# -------------------------------------------------------------------
data <- data.frame(myresponse = c(5.3, 0, 2.1),
                   x = c(-1.2, 3.2, -0.5),
                   y = c(765, 731, 353),
                   z = as.factor(c("foo", "bar", "foo")))
f <- myresponse ~ x + y + z
mf <- model.frame(f, data)
bk <- seq(0, 8, by = 1)
idx <- transitreg:::num2bin(data$myresponse, bk) # Convert response to 'pseudo-index'



# -------------------------------------------------------------------
# (1) Testing wrong usage
# -------------------------------------------------------------------
expect_error(transitreg:::transitreg_tmf(), info = "All arguments missing")
expect_error(transitreg:::transitreg_tmf(mf), info = "Most arguments missing")
expect_error(transitreg:::transitreg_tmf(mf, "myresponse"), info = "Not enough required args")

# Most basic use; should work, including the version with all default
# arguments specified. Used to check the function works before testing
# the remaining sanity checks.
expect_silent(x1 <- transitreg:::transitreg_tmf(mf, "myresponse", bk),
              info = "Basic working test")

args_OK <- list(data = mf, response = "myresponse", breaks = bk,
                theta_vars = NULL, newresponse = NULL, scaler = NULL, verbose = FALSE)
expect_silent(x2 <- do.call(transitreg:::transitreg_tmf, args_OK),
              info = "Basic working test with defaults specified")
expect_identical(x1, x2, info = "Testing return w/o and w/ default inputs")

rm(x1, x2)

# Testing sanity
args <- args_OK; args$data <- list(a = 1:3, b = 1:3)
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "data is not a data.frame")
args <- args_OK; args$data <- data.frame(a = NULL, b = NULL)
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "data with 0-dimension")

args <- args_OK; args$response <- 3
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "response is not NULL or character")
args <- args_OK; args$response <- c("foo", "bar")
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "response not NULL or char length 1")

args <- args_OK; args$newresponse <- TRUE
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "newresponse not NULL or numeric")
args <- args_OK; args$newresponse <- numeric()
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "newresponse not NULL or numeric length > 1")
args <- args_OK; args$response <- "foo_bar"
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "'response' not found in 'data', and no newrsponse provided")


args <- args_OK; args$theta_vars <- TRUE
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "theta_vars not NULL or character")

args <- args_OK; args$breaks <- TRUE
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "breaks not numeric")
args <- args_OK; args$breaks <- numeric()
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "breaks must be of length > 0")

args <- args_OK; args$newresponse <- -3.5
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "newresponse not NULL or integer")
args <- args_OK; args$newresponse <- c(-5L, 3L)
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "newresponse integer but not in {0, 1, 2, ...}")

args <- args_OK; args$verbose <- logical()
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")
args <- args_OK; args$verbose <- c(TRUE, TRUE)
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")
args <- args_OK; args$verbose <- 3
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")

args <- args_OK; args$scaler <- "yes"
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")
args <- args_OK; args$scaler <- 1:5
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")
args <- args_OK; args$scaler <- rep(TRUE, 3)
expect_error(do.call(transitreg:::transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")

rm(args, args_OK)


# -------------------------------------------------------------------
# (1) Simple case without theta vars and scaler
# -------------------------------------------------------------------
expect_silent(tmf <- transitreg:::transitreg_tmf(mf, "myresponse", bk),  info = "Testing basic functionality")
expect_inherits(tmf, "data.frame",                          info = "Testing tmf return class")
expect_identical(attr(tmf, "response"), "myresponse",       info = "Testing tmf attribute 'response'")
expect_identical(names(tmf), c("index", "Y", "theta", "myresponse", "x", "y", "z"),
              info = "Testing tmf return variable names")

expect_identical(tmf$index, rep(1:3, idx + 1),              info = "Testing 'index' variable")
expect_identical(tmf$Y, do.call(c, lapply(idx, function(k) rep(1:0, c(k, 1))[seq_len(k + 1)])),
                 info = "Testing 'index' variable")
expect_identical(tmf$theta, do.call(c, lapply(idx, function(k) seq_len(k + 1) - 1L)),
                 info = "Testing 'theta' variable")
expect_identical(tmf$myresponse, rep(idx, idx + 1),          info = "Testing 'myresponse' variable")
expect_identical(tmf$x, rep(data$x, idx + 1),                info = "Testing 'x' variable")
expect_identical(tmf$y, rep(data$y, idx + 1),                info = "Testing 'y' variable")
expect_identical(tmf$z, rep(data$z, idx + 1),                info = "Testing 'z' variable")

# If 'response = NULL' it assumes the first variable is our response, thus 
# the following should result in the identical object.
expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, NULL, bk),  info = "Testing 'response = NULL'")
expect_identical(tmf2, tmf, info = "Testing return object when using transitreg_tmf with 'response = NULL'")
rm(tmf2)

# -------------------------------------------------------------------
# (2) Testing newresponse behavior and response misuse
# -------------------------------------------------------------------

# Doing the same as above, but providing the response via 'newresponse'
# and not via the data argument (mf) as vector (identical to current mf$myresponse)
expect_silent(tmf2 <- transitreg:::transitreg_tmf(subset(mf, select = -myresponse), "myresponse",
                                     bk, newresponse = idx),  info = "Testing 'newresponse'")
expect_identical(tmf2, tmf, info = "Testing newresponse: providing response via newrespose only")

# Providing newresponse in addition to mf$myresponse, will throw a warning and
# mf$myresponse will be overwritten with 'newresponse', same return.
expect_warning(tmf2 <- transitreg:::transitreg_tmf(mf, "myresponse",
                                     bk, newresponse = idx),  info = "Testing 'newresponse'")
expect_identical(tmf2, tmf, info = "Testing newresponse: providing response via data and newrespose")

# New 'newresponse' vector, will result in a new return, thus testing it again.
expect_silent(tmf2 <- transitreg:::transitreg_tmf(subset(mf, select = -myresponse), "myresponse",
                                     bk, newresponse = c(4L, 4L, 4L)),  info = "Testing new 'newresponse' (vector)")
tmf2
expect_identical(tmf2$index, rep(1:3, each = 5L), info = "Testing tmf2$index")
expect_identical(tmf2$Y, rep(rep(1:0, c(4, 1)), 3), info = "Testing tmf2$Y")
expect_identical(tmf2$myresponse, rep(4L, 15), info = "Testing tmf2$myresponse")
expect_identical(tmf2$x, rep(mf$x, each = 5), info = "Testing tmf2$x")
expect_identical(tmf2$y, rep(mf$y, each = 5), info = "Testing tmf2$y")
expect_identical(tmf2$z, rep(mf$z, each = 5), info = "Testing tmf2$z")

# Providing newresponse as a single value, will be re-used for each observation in 'mf',
# resulting in the same object as above (as we use 4L).
expect_silent(tmf3 <- transitreg:::transitreg_tmf(subset(mf, select = -myresponse), "myresponse",
                                     bk, newresponse = 4L),  info = "Testing new 'newresponse' (scalar)")
expect_identical(tmf2, tmf3, info = "Testing tmf with 'newresponse' vector vs. scalar")

rm(tmf2, tmf)

# Testing misuse
# Must fail if response contains missing values
mfx <- mf; mfx$myresponse[2L] <- NA
expect_error(transitreg:::transitreg_tmf(mfx, "myresponse", bk),
             pattern = "Response contains missing values \\(not allowed\\)\\.",
             info = "Testing error when response contains missing values")

expect_error(transitreg:::transitreg_tmf(mfx, "myresponse", bk, newresponse = integer()),
             info = "Testing newresponse must be integer with length > 0")
expect_error(transitreg:::transitreg_tmf(mfx, "myresponse", bk, newresponse = 1.3),
             info = "Testing newresponse must be integer with length > 0")
expect_error(transitreg:::transitreg_tmf(mfx, "myresponse", bk, newresponse = c(1L, NA_integer_, 3L)),
             pattern = "Response contains missing values \\(not allowed\\)\\.",
             info = "Testing error when response contains missing values")

# Must fail when newresponse does not contain positive integers
expect_error(transitreg:::transitreg_tmf(mfx, "myresponse", bk, newresponse = c(0L, -1L, 3L)),
             info = "Testing error when newresponse contains negative indizes")


# -------------------------------------------------------------------
# (3) Testing 'theta_vars' behavior.
# -------------------------------------------------------------------
# Adding binary flags for theta2, theta5, and theta10
tvars <- c("theta2", "theta5", "theta10")
expect_silent(tmf <- transitreg:::transitreg_tmf(mf, "myresponse", bk))
expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, "myresponse", bk, theta_vars = tvars), 
              info = "Testing theta_vars behavior")

# The first part should be identical to 'tmf' above.
expect_equivalent(subset(tmf2, select = c(index, Y, theta, myresponse, x, y, z)), tmf,
              info = "Compare 'basic' part against tmf without theta_vars")

# Now testing the 'tvars'.
expect_true(all(c("theta2", "theta5") %in% names(tmf2)), info = "Test that 'theta2', 'theta5' were added")
expect_true(!"theta10" %in% names(tmf2), info = "Text that 'theta10' was not created (outside range)")

# Testing content of theta2 and theta5
expect_identical(tmf2$theta2, ifelse(tmf2$theta == 2L, 1L, 0L), info = "Testing convent of 'theta2'")
expect_identical(tmf2$theta5, ifelse(tmf2$theta == 5L, 1L, 0L), info = "Testing convent of 'theta5'")

rm(tmf2)


# -------------------------------------------------------------------
# (4) Testing 'scaler' behavior
# - `scaler = NULL`;   no scaling applied to 'x'
# - `scaler = FALSE`:  no scaling applied to 'x'
# - `scaler = TRUE`:   scaler applied forst time, returned as attribute on return
# - `scaler = list`:   existing scaler is applied
# -------------------------------------------------------------------

# This is the default behavior; scaler = NULL
expect_silent(tmf <- transitreg:::transitreg_tmf(mf, NULL, bk), info = "Testing scaler, reference object")

# Scaler NULL must result in the same object
expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, NULL, bk), info = "Testing scaler = NULL explicitly")
expect_identical(tmf2, tmf,                info = "Testing return with scaler = NULL")

# Scaler FALSE should be identical as well (no scaling)
expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, NULL, bk, scaler = FALSE), info = "Testing scaler = FALSE explicitly")
expect_identical(tmf2, tmf,                info = "Testing return with scaler = FALSE")

# -------------------------------------------
# Scaler TRUE: Scaling is applied
expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, NULL, bk, scaler = TRUE), info = "Testing scaler = TRUE explicitly")
scaler <- attr(tmf2, "scaler")
expect_true(is.list(scaler),               info = "Checking that 'scaler' is returned as attribute")

# Expecting 'scaler' length 3, scaling theta, x, and y (our numeric covariates)
expect_identical(names(scaler), c("theta", "x", "y"), info = "Checking names of scaler list")

expect_equal(scaler$theta$mean, mean(tmf$theta), info = "Checking scaler$theta$mu")
expect_equal(scaler$theta$sd,   sd(tmf$theta),   info = "Checking scaler$theta$mu")

expect_equal(scaler$x$mean, mean(tmf$x),         info = "Testing scaler$x$mean")
expect_equal(scaler$x$sd,   sd(tmf$x),           info = "Testing scaler$x$sd")
expect_equal(scaler$y$mean, mean(tmf$y),         info = "Testing scaler$y$mean")
expect_equal(scaler$y$sd,   sd(tmf$y),           info = "Testing scaler$y$sd")

rm(tmf2, scaler)

# -------------------------------------------
# If scaler is a list, apply existing scaler. For testing we use some 'random'
# (artificial) means and sds to ensure the calculation works as expected.
scaler <- list(theta = list(mean = -5, sd = 10),
               x     = list(mean = 50, sd = 11),
               y     = list(mean = -100, sd = 2.2))

expect_silent(tmf2 <- transitreg:::transitreg_tmf(mf, NULL, bk, scaler = scaler), info = "Testing scaler = TRUE explicitly")
expect_identical(attr(tmf2, "scaler"), scaler,   info = "Test that 'scaler' is returned unchanged")

# Test calculations
expect_equal(tmf2$theta, (tmf$theta - scaler$theta$mean) / scaler$theta$sd,
            info = "Testing result after applying an existing scaler on theta")
expect_equal(tmf2$x, (tmf$x - scaler$x$mean) / scaler$x$sd,
            info = "Testing result after applying an existing scaler on x")
expect_equal(tmf2$y, (tmf$y - scaler$y$mean) / scaler$y$sd,
            info = "Testing result after applying an existing scaler on y")

rm(tmf2, scaler)

