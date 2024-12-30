# -------------------------------------------------------------------
# Testing implementation of CDF, PDF, pmax
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("TransitionModels"))
if (interactive()) library("tinytest")

# -------------------------------------------------------------------
# Drawing probability for x = 0:10 from a binomial distribution.
# This represents the probabilities for different outcomes
# for one specific observation y in Y.
# -------------------------------------------------------------------
Y  <- 0:10
pp <- dbinom(0:11, 7, 0.4)

# -------------------------------------------------------------------
# R implementation of calculating the CDF and PDF
# -------------------------------------------------------------------

#                        y                   y-1
# CDF: (1 - P(y = 0)) + sum (1 - P(y = i)) * prod P(y = j)
#                       i=1                  j=1
#
# Note that in the R implementation r is the index (r = 1 equals y = 0)
man_cdf <- function(r, p) {
    res <- (1 - p[1])
    if (r > 1) { for (i in 2:r) res <- res + (1 - p[i]) * prod(p[1:(i-1)]) } 
    return(res)
}

#          y-1
# PDF: 1 * prod P(y = i) * (1 - P(y = y))
#          i=1
#
# Note that in the R implementation r is the index (r = 1 equals y = 0)
man_pdf <- function(r, p) {
    res <- 1
    if (r > 1) { for (i in seq_len(r - 1)) res <- res * p[i] }
    return(res * (1 - p[r]))
}


# Calculate CDF and PDF for each y in Y
res_manual <- data.frame(pdf = sapply(Y + 1, man_pdf, p = pp),
                         cdf = sapply(Y + 1, man_cdf, p = pp),
                         y = Y)

# -------------------------------------------------------------------
# Doing the same but using the C implementation
# -------------------------------------------------------------------
# Again, r is the index not y itself, so r = 1 equals y = 0
tm_fun <- function(r, p) {
    # Index 42 is just a random index
    x <- as.data.frame(.Call("tm_predict_pdfcdf", 42L, rep(42L, r), p = p[1:r], 1L,
                   PACKAGE = "TransitionModels"))
    return(transform(x, y = r - 1))
}

# Testing if tm_predict_pdfcdf is callable and returns what it should return
expect_silent(tmp <- tm_fun(3, pp),
              info = "Calling tm_predict_pcfcdf, checking return")
expect_inherits(tmp, "data.frame")
expect_identical(dim(tmp), c(1L, 3L))
expect_identical(names(tmp), c("pdf", "cdf", "y"))
rm(tmp)

# Calculate PDF and CDF for all y in Y
expect_silent(res_tm <- do.call(rbind, lapply(Y + 1, tm_fun, p = pp)),
              info = "Calculating CDF/PDF for y in Y")
expect_inherits(res_tm, "data.frame")
expect_identical(dim(res_tm), c(11L, 3L))
expect_identical(names(res_tm), c("pdf", "cdf", "y"))
rm(tm_fun)

# -------------------------------------------------------------------
# When everything works as expected, res_manual and res_tm
# should be all equal.
# -------------------------------------------------------------------
expect_equal(res_tm, res_manual,
             info = "Comparing results of manual and C implementation")

# -------------------------------------------------------------------
# Above we've used tm_predict_pdfcdf which calcultes both CDF and PDF,
# now testing against tm_predict which returns either CDF or PDF,
# or pmax (tested later).
# -------------------------------------------------------------------
# prob = 42 is a dummy value (must be double; not used)
tm_fun2 <- function(r, p, type, prob = 42.0) {
    # Index 666 is just a random index
    .Call("tm_predict", 666L, rep(666L, r), p = p[1:r], type, prob, 1L,
          PACKAGE = "TransitionModels")
}

# First testing if the function tm_fun2 works as expected, and that
# the C function is callable as expected.
expect_silent(tmp <- tm_fun2(3L, pp, "cdf"), info = "Testing tm_fun2 for type CDF")
expect_true(is.numeric(tmp) && length(tmp) == 1L, info = "Returns one single numeric")
rm(tmp)
expect_silent(tmp <- tm_fun2(3L, pp, "pdf"), info = "Testing tm_fun2 for type CDF")
expect_true(is.numeric(tmp) && length(tmp) == 1L, info = "Returns one single numeric")
rm(tmp)

# Calculating PDF/CDF using tm_predict for all y in Y and compare
# against the results from above; should be all equal.
expect_silent(res_tm2 <- data.frame(pdf = sapply(Y + 1, tm_fun2, p = pp, type = "pdf"),
                                    cdf = sapply(Y + 1, tm_fun2, p = pp, type = "cdf"),
                                    y   = Y))
expect_equal(res_tm, res_tm2,
             info = "Comparing results from tm_predict_pdfcdf and tm_predict")


# -------------------------------------------------------------------
# Checking 'pmax'
# -------------------------------------------------------------------

# pmax; internally uses man_pdf and returns the position of maximum
# probability. Used to check correctness of tm_predict(..., type = "pmax")
# implemented in C.
man_pmax <- function(r, p) {
    x <- sapply(seq_len(r), man_pdf, p = p)
    return(which.max(x) - 1)
}
pmax_man <- sapply(Y + 1, man_pmax, p = pp)

# Re-using 'tm_fun2' defined above with type = 'pmax'
expect_silent(pmax_tm <- sapply(Y + 1, tm_fun2, p = pp, type = "pmax"),
              info = "Getting pmax via C")
expect_true(is.numeric(pmax_tm) && length(pmax_tm) == length(Y),
            info = "Checking return")
expect_identical(pmax_man, pmax_tm,
            info = "Comparing results for 'pmax'; manual implementation vs. C")

## TODO(R): Is pmax doing what it is intended to do and in any kind useful?
