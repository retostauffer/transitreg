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










