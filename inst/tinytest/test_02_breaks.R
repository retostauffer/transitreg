# -------------------------------------------------------------------
# Testing the make_breaks function used in transitreg() to catch
# all possible responses, breaks, and censoring.
# -------------------------------------------------------------------


if (interactive()) { library("tinytest"); library("transitreg") }


# Drawing a number of count data {0, ..., 9}
set.seed(111)
xint <- pmin(rpois(500, lambda = 3), 9)
eps  <- sqrt(.Machine$double.eps) * 0.95
xint <- xint + runif(length(xint), -eps, eps)
expect_true(min(round(xint, eps * 2)) == 0L && max(round(xint, eps * 2)) == 9L) # Better safe than sorry

# -------------------------------------------------------------------
# Count data
# -------------------------------------------------------------------

###
msg <- "Integer response (small y), breaks NULL, uncensored"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = NULL, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b, "integer",
        info = paste(msg, "return class."))
expect_identical(b, 0:27,
        info = paste(msg, "values."))
expect_identical(b, transitreg:::make_breaks(xint),
        info = paste(msg, "comparing to default."))

###
msg <- "Integer response, breaks 9 (too small), uncensored"
expect_error(b <- transitreg:::make_breaks(xint, breaks = 9L, censored = "uncensored"),
        pattern = "response looks like count data, \\'breaks\\' \\(integer\\) must be >= max\\(response\\) \\+ 1",
        info = paste(msg, "throws error."))

###
msg <- "Integer response, breaks 10 (just enough), uncensored"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = 10, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b, "integer",
        info = paste(msg, "return class."))
expect_identical(b, 0:10,
        info = paste(msg, "values."))

###
msg <- "Integer response, breaks 555 (plenty), uncensored"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = 555, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b, "integer",
        info = paste(msg, "return class."))
expect_identical(b, 0:555,
        info = paste(msg, "values."))

# ---------------------
# Testing expansion

# Lower or equal 10L

##
msg <- "Integer response, max(y) == 5, uncensored, expecting max(y) * 3 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 5L)),      info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:15,                                    info = paste(msg, "values."))

##
msg <- "Integer response, max(y) == 10, uncensored, expecting max(y) * 3 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L)),     info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:30,                                    info = paste(msg, "values."))

# Lower or equal 100

##
msg <- "Integer response, max(y) == 99, uncensored, expecting max(y) * 1.5 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 99L)),     info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:148,                                   info = paste(msg, "values."))

##
msg <- "Integer response, max(y) == 100, uncensored, expecting max(y) * 1.5 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 100L)),    info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:150,                                   info = paste(msg, "values."))

# Lower or equal 1000

##
msg <- "Integer response, max(y) == 999, uncensored, expecting max(y) * 1.5 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 999L)),    info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:1248,                                  info = paste(msg, "values."))

##
msg <- "Integer response, max(y) == 1000, uncensored, expecting max(y) * 1.25 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 1000L)),   info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:1250,                                  info = paste(msg, "values."))

##
msg <- "Integer response, max(y) == 5000, uncensored, expecting max(y) * 1.25 as upper limit"
expect_silent(b <- transitreg:::make_breaks(c(0L, 5000L)),   info = paste(msg, "execution silent."))
expect_inherits(b, "integer",                                info = paste(msg, "return class."))
expect_identical(b, 0:6250,                                  info = paste(msg, "values."))
