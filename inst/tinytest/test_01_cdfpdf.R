# -------------------------------------------------------------------
# Testing implementation of CDF, PDF, pmax
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("TransitionModels"))

# -------------------------------------------------------------------
# Drawing probability for x = 0:15 from a Poisson distribution.
# Convert CDF to transition probabilities for testing.
# -------------------------------------------------------------------
pp <- ppois(0:15, lambda = 5)
dd <- dpois(0:15, lambda = 5)


# Testing misuse
expect_error(tm_convert(), info = "No input arguments")
expect_error(tm_convert(0.5, "cdf"), info = "Not enough input arguments")
expect_error(tm_convert(pp, from = 1, to = "cdf"), pattern = "'arg' should be one of", info = "Argument 'from' wrong")
expect_error(tm_convert(pp, from = "cdf", to = "foo"), pattern = "'arg' should be one of", info = "Argument 'to' wrong")
expect_error(tm_convert(LETTERS[1:3], "cdf", "tp"), info = "Argument 'x' not numeric")
expect_error(tm_convert(numeric(), "cdf", "tp"), info = "Argument 'x' of length 0")

# If input is tp: All values must be in [0, 1] and monotonically decreasing
expect_error(tm_convert(c(NA, 0.3, 0.5), "tp", "cdf"), info = "Missing value(s)")

expect_error(tm_convert(seq(1.2, 0.5, length.out = 3), "tp", "cdf"), info = "Values outside expected range")
expect_error(tm_convert(seq(0.5, -0.5, length.out = 3), "tp", "cdf"), info = "Values outside expected range")
expect_error(tm_convert(seq(0.2, 0.5, length.out = 3), "tp", "cdf"), info = "Values not monotonically decreasing")

expect_error(tm_convert(seq(1.2, 0.5, length.out = 3), "cdf", "tp"), info = "Values outside expected range")
expect_error(tm_convert(seq(0.5, -0.5, length.out = 3), "cdf", "tp"), info = "Values outside expected range")
expect_error(tm_convert(seq(0.5, 0.2, length.out = 3), "cdf", "tp"), info = "Values not monotonically increasing")


# Convert cdf -> tp
expect_silent(tp <- tm_convert(pp, from = "cdf", to = "tp"),
              info = "Converstion from cdf to tp (width = NULL)")
expect_identical(tm_convert(pp, "cdf", "tp", NULL, TRUE),
              tp, info = "Testing default arguments and argument order")

# Converting tp back to CDF, PDF (width = NULL)
expect_silent(p1 <- tm_convert(tp, from = "tp", to = "cdf"),  info = "Covnerting from tp to cdf")
expect_equal(p1, pp,                                          info = "Result of tp -> cdf")
expect_silent(d1 <- tm_convert(tp, from = "tp", to = "pdf"),  info = "Covnerting from tp to pdf")
expect_equal(d1, dd,                                          info = "Result of tp -> pdf")
rm(d1, p1)

# Checking evaluation of from/to
expect_equal(tp, tm_convert(pp, "CDF", "TP"))
expect_equal(tp, tm_convert(pp, "C", "T"))

# Convert to cdf/pdf all at once
expect_silent(pd1 <- tm_convert(tp, from = "tp", to = c("cdf", "pdf")),
              info = "Covnerting from tp to cdf and pdf")
expect_equal(pd1, data.frame(cdf = pp, pdf = dd),             info = "Testing result")
expect_silent(pd2 <- tm_convert(tp, from = "tp", to = c("cdf", "pdf"), drop = TRUE),
              info = "Testing that drop = TRUE has no effect here")
expect_equal(pd2, data.frame(cdf = pp, pdf = dd),             info = "Testing result")
rm(pd1, pd2)

# Keeping names
names(tp) <- LETTERS[seq_along(tp)]
expect_silent(p1 <- tm_convert(tp, "tp", "cdf"),              info = "Covnerting from tp to cdf (named)")
expect_equal(setNames(pp, names(tp)), p1,                     info = "Result of tp -> cdf (named)")
expect_silent(pd <- tm_convert(tp, "tp", c("cdf", "pdf")),    info = "Covnerting from tp to cdf and pdf (named)")
expect_equal(structure(data.frame(cdf = pp, pdf = dd), row.names = names(tp)), pd, info = "Result of tp to cdf and pdf (named)")
rm(p1, pd)

# "WIDTH" currently not implemented
expect_error(tm_convert(tp, "tp", "cdf", width = 1.5),       info = "Currently throws an error (TODO(R))")


# -------------------------------------------------------------------
# R implementation of calculating the CDF and PDF
# -------------------------------------------------------------------
fn <- function(i, tp) {
    tp <- tp[seq_len(i)] # Take first 'i' values
    res <- .Call("tm_predict_pdfcdf", 42L, rep(42L, length(tp)), tp = tp, 1L,
                 PACKAGE = "TransitionModels")
    return(data.frame(res))
}
tp <- tm_convert(pp, from = "cdf", to = "tp")
res <- do.call(rbind, lapply(seq_along(tp), fn, tp = tp))
expect_equal(res$pdf, dd)
expect_equal(res$cdf, pp)

# Calling tm_predict
lo <- 0:16 - 0.5
up <- 0:15 + 0.5
p3 <- .Call("tm_predict", uidx = 3L, idx = rep(3L, length(tp)), tp = tp,
            lower = lo, upper = up, y = as.numeric(0:15), type = "cdf",
            cores = 1L, elementwise = FALSE, discrete = TRUE)
expect_equal(pp, p3)
d3 <- .Call("tm_predict", uidx = 3L, idx = rep(3L, length(tp)), tp = tp,
            lower = lo, upper = up, y = as.numeric(0:15), type = "pdf",
            cores = 1L, elementwise = FALSE, discrete = TRUE)
expect_equal(dd, d3)



# -------------------------------------------------------------------
# Checking 'pmax'
# -------------------------------------------------------------------

# TODO(R): Update and check 'pmax', write tests.

