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
# Testing wrong use (sanity and co)
# -------------------------------------------------------------------

# Misspecified argument 'y'
expect_error(transitreg:::make_breaks(vector("numeric"), NULL, "uncensored"),
        pattern = "'y' must be numeric length > 0L",
        info = "y of length 0L")
expect_error(transitreg:::make_breaks("foo", NULL, "uncensored"),
        pattern = "'y' must be numeric length > 0L",
        info = "y not numeric")
expect_error(transitreg:::make_breaks(TRUE, NULL, "uncensored"),
        pattern = "'y' must be numeric length > 0L",
        info = "y not numeric")

# Misspecified argument 'breaks'
expect_error(transitreg:::make_breaks(0:10, vector("numeric"), "uncensored"),
        pattern = "'breaks' must be 'NULL' or numeric vector of length 1L or length > 2L",
        info = "breaks of length 0L")
expect_error(transitreg:::make_breaks(0:10, c(0, 10), "uncensored"),
        pattern = "'breaks' must be 'NULL' or numeric vector of length 1L or length > 2L",
        info = "breaks of length 0L")
expect_error(transitreg:::make_breaks(0:10, "foo", "uncensored"),
        pattern = "'breaks' must be 'NULL' or numeric vector of length 1L or length > 2L",
        info = "breaks not numeric")
expect_error(transitreg:::make_breaks(0:10, TRUE, "uncensored"),
        pattern = "'breaks' must be 'NULL' or numeric vector of length 1L or length > 2L",
        info = "breaks not numeric")

expect_error(transitreg:::make_breaks(xint, breaks = 0:23, censored = 0),
        pattern = "if 'censored' is numeric, breaks must be NULL or a single numeric value",
        info = "Testing inproper use")
expect_error(transitreg:::make_breaks(xint, breaks = 0:23, censored = c(NA, 10)),
        pattern = "if 'censored' is numeric, breaks must be NULL or a single numeric value",
        info = "Testing inproper use")
expect_error(transitreg:::make_breaks(xint, breaks = 0:23, censored = c(0, NA)),
        pattern = "if 'censored' is numeric, breaks must be NULL or a single numeric value",
        info = "Testing inproper use")
expect_error(transitreg:::make_breaks(xint, breaks = 0:23, censored = c(0, 10)),
        pattern = "if 'censored' is numeric, breaks must be NULL or a single numeric value",
        info = "Testing inproper use")


# -------------------------------------------------------------------
# Count data (response looks like positive integers or are integers)
# -------------------------------------------------------------------


###
msg <- "Integer response (small y), breaks NULL, uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = NULL, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",         info = paste(msg, "return class."))
expect_identical(b$breaks, 0:27,             info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",   info = paste(msg, "censored return."))
expect_identical(b, transitreg:::make_breaks(xint, NULL, "uncensored"),
        info = "Testing default argument order.")

###
msg <- "Integer response, breaks 9 (too small), uncensored; checking"
expect_error(b <- transitreg:::make_breaks(xint, breaks = 9L, censored = "uncensored"),
        pattern = "response looks like count data, \\'breaks\\' \\(integer\\) must be >= max\\(response\\) \\+ 1",
        info = paste(msg, "throws error."))

###
msg <- "Integer response, breaks 10 (just enough), uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = 10, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",         info = paste(msg, "return class."))
expect_identical(b$breaks, 0:10,             info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",   info = paste(msg, "censored return."))

###
msg <- "Integer response, breaks 555 (plenty), uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(xint, breaks = 555, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:555,                info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))


# ---------------------
# Testing expansion

# Lower or equal 10L

##
msg <- "Integer response, max(y) == 5, uncensored, expecting max(y) * 3 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 5L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:15,                 info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Integer response, max(y) == 10, uncensored, expecting max(y) * 3 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:30,                 info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

# Lower or equal 100

##
msg <- "Integer response, max(y) == 99, uncensored, expecting max(y) * 1.5 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 99L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:148,                info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Integer response, max(y) == 100, uncensored, expecting max(y) * 1.5 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 100L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:150,                info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

# Lower or equal 1000

##
msg <- "Integer response, max(y) == 999, uncensored, expecting max(y) * 1.5 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 999L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:1248,               info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Integer response, max(y) == 1000, uncensored, expecting max(y) * 1.25 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 1000L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:1250,               info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Integer response, max(y) == 5000, uncensored, expecting max(y) * 1.25 as upper limit; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 5000L), NULL, "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, 0:6250,               info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))


# ---------------------
# Testing censoring: character input

##
msg <- "Integer response, censored = \"left\", breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, "left"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(0L, 0:20),          info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

##
msg <- "Integer response, censored = \"left\", breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, "right"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(0:20, 20L),         info = paste(msg, "values."))
expect_identical(b$censored, "right",            info = paste(msg, "censored return."))

##
msg <- "Integer response, censored = \"both\", breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, "both"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(0L, 0:20, 20L),     info = paste(msg, "values."))
expect_identical(b$censored, "both",             info = paste(msg, "censored return."))


# ---------------------
# Testing censoring: numeric input

##
msg <- "Integer response, censored = 5.1, breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, 5.1),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(5L, 5:25),          info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

##
msg <- "Integer response, censored = c(5.1, NA), breaks 30; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 30, c(5.1, NA)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(5L, 5:35),          info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))


# TODO(R): 'breaks' argument is completely obsolete when having right break?
##
msg <- "Integer response, censored = c(NA, 10.1), breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, c(NA, 10.1)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(0:10, 10L),         info = paste(msg, "values."))
expect_identical(b$censored, "right",            info = paste(msg, "censored return."))

##
msg <- "Integer response, censored = c(10, 20), breaks 20; checking"
expect_silent(b <- transitreg:::make_breaks(c(0L, 10L), 20, c(10, 20)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(10L, 10:20, 20L),   info = paste(msg, "values."))
expect_identical(b$censored, "both",             info = paste(msg, "censored return."))


# ---------------------
# Providing breaks vector
expect_error(transitreg:::make_breaks(0:30, breaks = 10:20, censored = "uncensored"),
             pattern = "breaks not spanning the range of the response.*",
             info = "Breaks not spanning range of data.")
expect_error(transitreg:::make_breaks(0:20, breaks = 10:20, censored = "uncensored"),
             pattern = "breaks not spanning the range of the response.*",
             info = "Breaks not spanning range of data.")
expect_error(transitreg:::make_breaks(0:20, breaks = 10:20, censored = "right"),
             pattern = "breaks not spanning the range of the response.*",
             info = "Breaks not spanning range of data.")
expect_error(transitreg:::make_breaks(10:30, breaks = 0:20, censored = "left"),
             pattern = "breaks not spanning the range of the response.*",
             info = "Breaks not spanning range of data.")

##
msg <- "Integer response, breaks provided, uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(10:20, breaks = 5:35, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(5:35),              info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Integer response, breaks provided, censored left; checking"
expect_silent(b <- transitreg:::make_breaks(10:20, breaks = 5:35, censored = "left"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(5L, 5:35),          info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

##
msg <- "Integer response, breaks provided, censored right; checking"
expect_silent(b <- transitreg:::make_breaks(10:20, breaks = 5:35, censored = "right"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "integer",             info = paste(msg, "return class."))
expect_identical(b$breaks, c(5:35, 35L),         info = paste(msg, "values."))
expect_identical(b$censored, "right",            info = paste(msg, "censored return."))



# -------------------------------------------------------------------
# Numeric response
# -------------------------------------------------------------------

set.seed(100)
xnum <- c(-5, 5, round(runif(100, -5, 5), 2)) # -5 to +5

##
expect_error(b <- transitreg:::make_breaks(xnum, breaks = NULL, censored = "uncensored"),
        pattern = "response does not look like count data, in this case 'breaks' must be specified",
        info = "Breaks must be specified for numeric response.")
expect_error(b <- transitreg:::make_breaks(xnum, breaks = 2, censored = "uncensored"),
        pattern = "if 'breaks' is a single numeric, it must be >= 3",
        info = "Breaks must be specified for numeric response.")

##
msg <- "Numeric response, breaks provided, uncensored; checking"
bk <- c(-5, 0, 5)
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = bk, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(b$breaks, bk,                   info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Numeric response, breaks single numeric, uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 10L,          info = paste(msg, "length of breaks vector"))
tmp <- seq(-6, 6, length.out = 10)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))

##
msg <- "Numeric response, breaks single numeric, uncensored; checking"
expect_silent(b <- transitreg:::make_breaks(xnum * 100, breaks = 100, censored = "uncensored"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 100L,         info = paste(msg, "length of breaks vector"))
tmp <- seq(-600, 600, length.out = 100)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "uncensored",       info = paste(msg, "censored return."))


# ---------------------
# Censoring via character

##
msg <- "Numeric response, breaks single numeric, censored left; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = "left"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 11L,          info = paste(msg, "length of breaks vector"))
tmp <- c(-6, seq(-6, 6, length.out = 10))
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

##
msg <- "Numeric response, breaks single numeric, censored right; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = "right"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 11L,          info = paste(msg, "length of breaks vector"))
tmp <- c(seq(-6, 6, length.out = 10), +6)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "right",            info = paste(msg, "censored return."))

##
msg <- "Numeric response, breaks single numeric, censored both; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = "both"),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 12L,          info = paste(msg, "length of breaks vector"))
tmp <- c(-6, seq(-6, 6, length.out = 10), +6)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "both",             info = paste(msg, "censored return."))


# ---------------------
# Censoring via numeric

## Left censoring ...

expect_error(transitreg:::make_breaks(xnum, breaks = 10, censored = +10),
        pattern = "censoring point > max\\(y\\)",
        info = "Invalid left censoring point")

##
msg <- "Numeric response, left censoring via numeric vector; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = +2),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 11L,          info = paste(msg, "length of breaks vector"))
tmp <- c(2, seq(2, 6, length.out = 10))
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

##
msg <- "Numeric response, left censoring via numeric vector; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = c(+2, NA)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 11L,          info = paste(msg, "length of breaks vector"))
tmp <- c(2, seq(2, 6, length.out = 10))
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "left",             info = paste(msg, "censored return."))

## Right censoring ...

expect_error(transitreg:::make_breaks(xnum, breaks = 10, censored = c(NA, -10)),
        pattern = "censoring point < min\\(y\\)",
        info = "Invalid left censoring point")

##
msg <- "Numeric response, right censoring via numeric vector; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = c(NA, +2)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 11L,          info = paste(msg, "length of breaks vector"))
tmp <- c(seq(-6, 2, length.out = 10), 2)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "right",            info = paste(msg, "censored return."))

## Both censoring ...

expect_error(transitreg:::make_breaks(xnum, breaks = 10, censored = c(10, -10)),
        pattern = "invalid censoring points",
        info = "Left censoring point > right censoring point")
expect_error(transitreg:::make_breaks(xnum, breaks = 10, censored = c(10, 10)),
        pattern = "invalid censoring points",
        info = "Left censoring point == right censoring point")
expect_error(transitreg:::make_breaks(xnum, breaks = 10, censored = c(-100, -50)),
        pattern = "censoring outside range of data",
        info = "No data within censoring range")

##
msg <- "Numeric response, left-right censoring via numeric vector; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = c(-10, 10)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 12L,          info = paste(msg, "length of breaks vector"))
tmp <- c(-10, seq(-10, 10, length.out = 10), 10)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "both",             info = paste(msg, "censored return."))

##
msg <- "Numeric response, left-right censoring via numeric vector; checking"
expect_silent(b <- transitreg:::make_breaks(xnum, breaks = 10, censored = c(-3, +2.5)),
        info = paste(msg, "execution silent."))
expect_inherits(b$breaks, "numeric",             info = paste(msg, "return class."))
expect_identical(length(b$breaks), 12L,          info = paste(msg, "length of breaks vector"))
tmp <- c(-3, seq(-3, 2.5, length.out = 10), 2.5)
expect_equal(b$breaks, tmp,                      info = paste(msg, "values."))
expect_identical(b$censored, "both",             info = paste(msg, "censored return."))



