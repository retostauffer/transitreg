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

# Making a copy for convenience
transitreg_tmf <- transitreg:::transitreg_tmf

# -------------------------------------------------------------------
# (1) Testing wrong usage
# -------------------------------------------------------------------
expect_error(transitreg_tmf(), info = "All arguments missing")
expect_error(transitreg_tmf(mf), info = "Most arguments missing")
expect_error(transitreg_tmf(mf, "myresponse"), info = "Not enough required args")

# Most basic use; should work, including the version with all default
# arguments specified. Used to check the function works before testing
# the remaining sanity checks.
expect_silent(x1 <- transitreg_tmf(mf, "myresponse", bk),
              info = "Basic working test")

args_OK <- list(data = mf, response = "myresponse", breaks = bk,
                theta_vars = NULL, newresponse = NULL, scaler = NULL, verbose = FALSE)
expect_silent(x2 <- do.call(transitreg_tmf, args_OK),
              info = "Basic working test with defaults specified")
expect_identical(x1, x2, info = "Testing return w/o and w/ default inputs")

rm(x1, x2)

# Testing sanity
args <- args_OK; args$response <- 3
expect_error(do.call(transitreg_tmf, args), info = "response is not NULL or character")
args <- args_OK; args$response <- c("foo", "bar")
expect_error(do.call(transitreg_tmf, args), info = "response not NULL or char length 1")

args <- args_OK; args$newresponse <- TRUE
expect_error(do.call(transitreg_tmf, args), info = "newresponse not NULL or numeric")
args <- args_OK; args$newresponse <- numeric()
expect_error(do.call(transitreg_tmf, args), info = "newresponse not NULL or numeric length > 1")
args <- args_OK; args$response <- "foo_bar"
expect_error(do.call(transitreg_tmf, args), info = "'response' not found in 'data', and no newrsponse provided")


args <- args_OK; args$theta_vars <- TRUE
expect_error(do.call(transitreg_tmf, args), info = "theta_vars not NULL or character")

args <- args_OK; args$breaks <- TRUE
expect_error(do.call(transitreg_tmf, args), info = "breaks not numeric")
args <- args_OK; args$breaks <- numeric()
expect_error(do.call(transitreg_tmf, args), info = "breaks must be of length > 0")

args <- args_OK; args$newresponse <- -3.5
expect_error(do.call(transitreg_tmf, args), info = "newresponse not NULL or integer")
args <- args_OK; args$newresponse <- c(-5L, 3L)
expect_error(do.call(transitreg_tmf, args), info = "newresponse integer but not in {0, 1, 2, ...}")

args <- args_OK; args$verbose <- logical()
expect_error(do.call(transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")
args <- args_OK; args$verbose <- c(TRUE, TRUE)
expect_error(do.call(transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")
args <- args_OK; args$verbose <- NULL
expect_error(do.call(transitreg_tmf, args), info = "verbose must be logical TRUE or FALSE")

args <- args_OK; args$scaler <- "yes"
expect_error(do.call(transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")
args <- args_OK; args$scaler <- 1:5
expect_error(do.call(transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")
args <- args_OK; args$scaler <- rep(TRUE, 3)
expect_error(do.call(transitreg_tmf, args), info = "scaler mut be NULL, TRUE, FALSE, or list")

rm(transitreg_tmf) # Convenience copy
rm(args, args_OK)


# -------------------------------------------------------------------
# (1) Simple case without theta vars and scaler
# -------------------------------------------------------------------
expect_silent(tmf <- transitreg:::transitreg_tmf(mf, "myresponse", bk),  info = "Testing basic functionality")
expect_inherits(tmf, "data.frame",                          info = "Testing tmf return class")
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


# Doing the same again, but providing the new response via 'newresponse' instead
# of 'data', should result in the very same object as above.

expect_silent(tmf2 <- transitreg:::transitreg_tmf(subset(mf, select = -myresponse), "myresponse",
                                     bk, newresponse = idx),  info = "Testing 'newresponse'")
expect_identical(tmf2, tmf)

rm(tmf2)



# -------------------------------------------------------------------
# (2) Testing 'theta_vars' behavior.
# -------------------------------------------------------------------
# Adding binary flags for theta2, theta5, and theta10
tvars <- c("theta2", "theta5", "theta10")
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
# (3) Testing 'scaler' behavior
# -------------------------------------------------------------------





