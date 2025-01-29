# -------------------------------------------------------------------
# Testing the functions in R/misc.R
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# # -------------------------------------------------------------------
# # Estimating simple model (from ?transitreg examples)
# # -------------------------------------------------------------------
# set.seed(123)
# n <- 1000
# x <- runif(n, -3, 3)
# y <- rpois(n, exp(2 + sin(x)))
# 
# # Fit transition count response model.
# b <- transitreg(y ~ s(theta) + s(x))


# -------------------------------------------------------------------
# 'grep2': Convenience function which allows to grep()
# for multiple patterns.
# -------------------------------------------------------------------
haystack <- c("Tallin", "Berlin", "Innsbruck", "London", "Dublin")
expect_identical(transitreg:::grep2("bruck", haystack),
                 grep("bruck", haystack))
expect_identical(transitreg:::grep2("lin$", haystack),
                 grep("lin$", haystack))
expect_identical(sort(transitreg:::grep2(c("Inns", "lin$"), haystack)),
                 sort(grep("(Inns|lin$)", haystack)))


# -------------------------------------------------------------------
# Testing num2bin and bin2num, converting from a numeric variable
# into pseudo-counts (indices) and back.
# -------------------------------------------------------------------
num2bin <- transitreg:::num2bin
bin2num <- transitreg:::bin2num

x <- runif(50, 10, 50)

bk <- seq(10, 50, by = 5)
expect_error(num2bin(), info = "All inputs missing")
expect_error(num2bin(x = x), info = "Argument 'breaks' missing")
expect_error(num2bin(breaks = breaks), info = "Argument 'x' missing")
expect_identical(num2bin(x, bk), num2bin(breaks = bk, x = x), info = "Testing default arguments order")
expect_error(num2bin(c(1, 2, NA, 3), bk), pattern = "Response contains missing values \\(not allowed\\)",
             info = "Testing error if input contains missing values")

# Check if the conversion returns what it should
expect_identical(num2bin(x, bk), cut(x, breaks = bk, labels = FALSE) - 1L,
                 info = "Testing conversion 'num2bin' via breaks")
expect_identical(num2bin(NULL, bk), NULL,
                 info = "Return NULL if input is NULL")

# Values below lowest break threshold: should result in -1
expect_identical(num2bin(c(-100, 0, min(bk) - sqrt(.Machine$double.eps)), bk), rep(-1L, 3),
                 info = "Testing return for values outside support (below min)")
expect_identical(num2bin(c(max(bk) + sqrt(.Machine$double.eps), 100, 100000), bk), rep(length(bk) - 1L, 3L),
                 info = "Testing return for values outside support (above max)")


# Now converting back from bins to numerics
idx <- num2bin(x, bk)

# First incomplete/wrong calls
expect_error(bin2num(), info = "All inputs missing")
expect_error(bin2num(idx = idx), info = "Argument 'breaks' missing")
expect_error(bin2num(breaks = breaks), info = "Argument 'idx' missing")
expect_identical(bin2num(idx, bk), bin2num(breaks = bk, idx = idx), info = "Testing default arguments order")
expect_identical(num2bin(NULL, bk), NULL, info = "Return NULL if input is NULL")

mids <- (head(bk, -1) + tail(bk, -1)) / 2.
expect_identical(bin2num(idx, bk), mids[idx + 1L])
rm(bin2num, num2bin, idx, x, bk)




# -------------------------------------------------------------------
# 'response_name': Convenience function to extract the response name
# of a formula object(s).
# -------------------------------------------------------------------
expect_identical(transitreg:::response_name(Y ~ x), "Y")
expect_identical(transitreg:::response_name(formula("Y ~ x")), "Y")
expect_identical(transitreg:::response_name(foo ~ a + b + c^2), "foo")

f <- list(formula(foo ~ a + b + c^2), formula(bar ~ poly(x)))
expect_identical(transitreg:::response_name(f), c("foo", "bar"))

# If input is a list with a $formula element, that one is taken
expect_identical(transitreg:::response_name(list(formula = Y ~ x, foo = 3)), "Y")


# -------------------------------------------------------------------
# 'fake_formula': Very basic checks
# -------------------------------------------------------------------

f <- formula(Y ~ a + b + c^2)
expect_identical(transitreg:::fake_formula(f), f)
rm(f)

f <- list(formula(foo ~ a + b + c^2), formula(bar ~ poly(x)))
expect_identical(lengths(f), c(3L, 3L))
f <- transitreg:::fake_formula(f)
expect_identical(length(f), c(2L, 2L))
expect_identical(length(f[[1]]), 1L)
expect_identical(length(f[[2]]), 3L)
expect_identical(format(f[[2]]), "a + b + c^2 | poly(x)")
rm(f)


f <- list(formula(foo ~ s(a) + te(b) + c^2), formula(bar ~ poly(x)))
expect_identical(lengths(f), c(3L, 3L))
f <- transitreg:::fake_formula(f)
expect_identical(length(f), c(2L, 2L))
expect_identical(length(f[[1]]), 1L)
expect_identical(length(f[[2]]), 3L)
expect_identical(format(f[[2]]), "c + a + b | poly(x)")


f <- list(formula(foo ~ s(a) + te(b) + c^2), formula(bar ~ poly(x)))

expect_identical(transitreg:::fake_formula(f),
                 transitreg:::fake_formula(f, specials = NULL,
                            nospecials = FALSE, onlyspecials = FALSE),
                 info = "Testing defaults of fake_formula")

fx <- transitreg:::fake_formula(f, nospecials = FALSE)
expect_identical(format(fx[[2]]), "c + a + b | poly(x)")
rm(fx)

fx <- transitreg:::fake_formula(f, nospecials = TRUE)
expect_identical(format(fx[[2]]), "c | poly(x)")
rm(fx)

fx <- transitreg:::fake_formula(f, onlyspecials = TRUE)
expect_true(is.list(fx) && length(fx) == 2L)
expect_identical(fx, list(c("s(a)", "te(b)"), character()))
rm(fx)


# -------------------------------------------------------------------
# 'ff_replace':
# -------------------------------------------------------------------

orig     <- formula("foo ~ a:b + bar*c + d + e + f:g:h")
##expected <- "foo ~ a + b + bar + c + d + e + f + g + h + bar:c"
expected <- "foo ~ a + b + bar + c + d + e + f + g + h"
expect_silent(tmp <- transitreg:::ff_replace(orig))
## FAILS ## expect_identical(format(tmp), expected)
# TODO(R): Niki, there is a remaining bar:c in the example above.
#          Is that intended?
rm(orig, expected, tmp)


orig     <- list(formula(foo ~ a:b), formula(bar ~ c*d*e))
##expected <- "foo ~ a + b + bar + c + d + e + f + g + h + bar:c"
expected <- "foo ~ a + b + bar + c + d + e + f + g + h"
expect_silent(tmp <- transitreg:::ff_replace(f))
# TODO(R): Niki ff_replace is intended for one single formula object I assume?
#          Given there are no sanity checks the call above
#          does not throw any error, but does not return the intended result right?


# -------------------------------------------------------------------
# 'get_elementwise_colnames': A helper function to return proper
# column names if cdf, pdf, quantile, ... are called with
# elementwise = FALSE.
# -------------------------------------------------------------------
fn <- transitreg:::get_elementwise_colnames
expect_true(is.function(fn), info = "Check if we can find 'get_elementwise_colnames'")

# --------------- quantiles --------------
# Testing return for quantiles; no prefix will convert the
# input 'p' into percentages. 
p <- seq(0.1, 0.9, by = 0.1)
expect_silent(cnames <- fn(p),
                info = "Create 'quantile' column names")
# Testing defaults
expect_silent(cnames2 <- fn(p, prefix = NULL),
                info = "Create column names with defaults")
expect_identical(cnames, cnames2,
                info = "Testing return using function defaults")
rm(cnames2)

# Checking expected names
expect_identical(cnames, paste0(seq.int(10, 90, by = 10), "%"),
                info = "Checking 'quantile' column names")
rm(cnames, p)

# ------------------ cdf -----------------
set.seed(6020)
digits <- pmax(3L, getOption("digits") - 3L)
y <- runif(100, -1000, 1000)

# Single digit
expect_silent(cnames <- fn(y, prefix = "x", digits = digits),
                info = "Create column names")
expect_identical(cnames, paste("x", trimws(format(y, digits = digits)), sep = "_"),
                info = "Checking column names")




# -------------------------------------------------------------------
# check_args_for_treg_predict and check_args_for_treg_predict_pdfcdf
# are two functions used to test if we have all the elements before
# calling the .C routines, as well as checking if the elements are
# of the correct class/type and/or correct length.
#
# Here we are testing a few special cases to cover the sanity
# checks, most is covered by other tests when running model estimation,
# prediction etc.
# -------------------------------------------------------------------

# check_args_for_treg_predict ---------------------------------------

# This list will succeed the tests (all fine)
args_OK <- list(uidx        = 1L,
                idx         = rep(1L, 5),
                tp          = runif(5),
                breaks      = c(0, 0.5, 1.0, 1.5),
                y           = c(5L, 3L),
                prob        = c(0.5, 0.6),
                type        = "pdf",
                ncores      = 10L,
                elementwise = TRUE,
                discrete    = FALSE,
                censored    = "left")

expect_silent(transitreg:::check_args_for_treg_predict(args_OK), info = "This should work like a charm")

# Testing 'class' errors
args <- args_OK; args$uidx <- 1.0
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$uidx not integer")
args <- args_OK; args$idx  <- 1.0
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$idx not integer")
args <- args_OK; args$tp   <- 1L
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$tp not double")
args <- args_OK; args$breaks <- 1L
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$breaks not double")
args <- args_OK; args$y    <- 5.0
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$y not integer")
args <- args_OK; args$prob <- 1L
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$prob not double")
args <- args_OK; args$type <- FALSE
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$type not character")
args <- args_OK; args$ncores <- 3.3
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$cores not integer")
args <- args_OK; args$elementwise <- 5
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$elementwise not logical")
args <- args_OK; args$discrete <- "foo"
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$discrete not logical")
args <- args_OK; args$discrete <- 1.234
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$censored no character")

# Testing content errors (mismatching lengths, ...)
args <- args_OK; args$elementwise <- logical()
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$elementwise not length 1")
args <- args_OK; args$idx <- rep(1L, 10L)
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$idx and args$tp length mismatch")
args <- args_OK; args$tp <- rnorm(3)
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$idx and args$tp length mismatch")
args <- args_OK; args$type <- c("pdf", "cdf")
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$type not of length 1")
args <- args_OK; args$uidx <- 1:3
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$uidx and args$discrete length mismatch")
args <- args_OK; args$discrete <- rep(TRUE, 10)
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$uidx and args$discrete length mismatch")
args <- args_OK; args$uidx <- NULL
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "Not all required elements in list")
args <- args_OK; args$censored <- c("left", "right")
expect_error(transitreg:::check_args_for_treg_predict(args, silent = TRUE), info = "args$censored must be of length 1")

rm(args, args_OK)


# This list will succeed the tests (all fine)
args_OK <- list(uidx        = 1:3,
                idx         = rep(1:3, 5),
                tp          = runif(15),
                y           = c(5L, 10L, 3L),
                breaks      = c(0, 0.5, 1.0, 1.5),
                ncores      = 10L,
                censored    = "left")

expect_silent(transitreg:::check_args_for_treg_predict_pdfcdf(args_OK), info = "This should work like a charm")

# Testing 'class' errors
args <- args_OK; args$uidx <- 1.0
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$uidx not integer")
args <- args_OK; args$idx  <- 1.0
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$idx not integer")
args <- args_OK; args$tp   <- 1L
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$tp not double")
args <- args_OK; args$breaks <- 1L
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$breaks not double")
args <- args_OK; args$y    <- 5.0
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$y not integer")
args <- args_OK; args$ncores <- 3.3
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$cores not integer")
args <- args_OK; args$censored <- TRUE
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$censored not character")

# Testing content errors (mismatching lengths, ...)
args <- args_OK; args$tp <- logical()
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$elementwise not length 1")
args <- args_OK; args$idx <- rep(1L, 10L)
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$idx and args$tp length mismatch")
args <- args_OK; args$tp <- rnorm(3)
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$idx and args$tp length mismatch")
args <- args_OK; args$censored <- c("left", "right")
expect_error(transitreg:::check_args_for_treg_predict_pdfcdf(args, silent = TRUE), info = "args$censored must be of length 1")

rm(args, args_OK)

