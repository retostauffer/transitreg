# -------------------------------------------------------------------
# Testing implementation of CDF, PDF, pmax
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# -------------------------------------------------------------------
# Drawing probability for x = 0:15 from a Poisson distribution.
# Convert CDF to transition probabilities for testing.
# -------------------------------------------------------------------
pp <- ppois(0:15, lambda = 5)
dd <- dpois(0:15, lambda = 5)


# Testing misuse
expect_error(convert_tp(), info = "No input arguments")
expect_error(convert_tp(0.5, "cdf"), info = "Not enough input arguments")
expect_error(convert_tp(pp, from = 1, to = "cdf"), pattern = "'arg' should be one of", info = "Argument 'from' wrong")
expect_error(convert_tp(pp, from = "cdf", to = "foo"), pattern = "'arg' should be one of", info = "Argument 'to' wrong")
expect_error(convert_tp(LETTERS[1:3], "cdf", "tp"), info = "Argument 'x' not numeric")
expect_error(convert_tp(numeric(), "cdf", "tp"), info = "Argument 'x' of length 0")

# If input is tp: All values must be in [0, 1] and monotonically decreasing
expect_error(convert_tp(c(NA, 0.3, 0.5), "tp", "cdf"), info = "Missing value(s)")

expect_error(convert_tp(seq(1.2, 0.5, length.out = 3), "tp", "cdf"), info = "Values outside expected range")
expect_error(convert_tp(seq(0.5, -0.5, length.out = 3), "tp", "cdf"), info = "Values outside expected range")
expect_error(convert_tp(seq(0.2, 0.5, length.out = 3), "tp", "cdf"), info = "Values not monotonically decreasing")

expect_error(convert_tp(seq(1.2, 0.5, length.out = 3), "cdf", "tp"), info = "Values outside expected range")
expect_error(convert_tp(seq(0.5, -0.5, length.out = 3), "cdf", "tp"), info = "Values outside expected range")
expect_error(convert_tp(seq(0.5, 0.2, length.out = 3), "cdf", "tp"), info = "Values not monotonically increasing")


# Convert cdf -> tp
expect_silent(tp <- convert_tp(pp, from = "cdf", to = "tp"),
              info = "Converstion from cdf to tp (width = NULL)")
expect_identical(convert_tp(pp, "cdf", "tp", NULL, TRUE),
              tp, info = "Testing default arguments and argument order")

# Converting tp back to CDF, PDF (width = NULL)
expect_silent(p1 <- convert_tp(tp, from = "tp", to = "cdf"),  info = "Covnerting from tp to cdf")
expect_equal(p1, pp,                                          info = "Result of tp -> cdf")
expect_silent(d1 <- convert_tp(tp, from = "tp", to = "pdf"),  info = "Covnerting from tp to pdf")
expect_equal(d1, dd,                                          info = "Result of tp -> pdf")
rm(d1, p1)

# Checking evaluation of from/to
expect_equal(tp, convert_tp(pp, "CDF", "TP"))
expect_equal(tp, convert_tp(pp, "C", "T"))

# Convert to cdf/pdf all at once
expect_silent(pd1 <- convert_tp(tp, from = "tp", to = c("cdf", "pdf")),
              info = "Covnerting from tp to cdf and pdf")
expect_equal(pd1, data.frame(cdf = pp, pdf = dd),             info = "Testing result")
expect_silent(pd2 <- convert_tp(tp, from = "tp", to = c("cdf", "pdf"), drop = TRUE),
              info = "Testing that drop = TRUE has no effect here")
expect_equal(pd2, data.frame(cdf = pp, pdf = dd),             info = "Testing result")
rm(pd1, pd2)

# Keeping names
names(tp) <- LETTERS[seq_along(tp)]
expect_silent(p1 <- convert_tp(tp, "tp", "cdf"),              info = "Covnerting from tp to cdf (named)")
expect_equal(setNames(pp, names(tp)), p1,                     info = "Result of tp -> cdf (named)")
expect_silent(pd <- convert_tp(tp, "tp", c("cdf", "pdf")),    info = "Covnerting from tp to cdf and pdf (named)")
expect_equal(structure(data.frame(cdf = pp, pdf = dd), row.names = names(tp)), pd, info = "Result of tp to cdf and pdf (named)")
rm(p1, pd)

# Testing width argument
expect_identical(convert_tp(tp, "tp", "cdf", width = NULL),
                 convert_tp(tp, "tp", "cdf", width = 1.234),
                 info = "Width must have no effect on CDF")

expect_identical(convert_tp(tp, "tp", "pdf", width = NULL),
                 convert_tp(tp, "tp", "pdf", width = 1),
                 info = "Result with width = NULL and width 1 must be identical")
expect_identical(convert_tp(tp, "tp", "pdf", width = NULL),
                 convert_tp(tp, "tp", "pdf", width = rep(1.0, length(tp))),
                 info = "Result with width = NULL and width 1 must be identical")

expect_identical(convert_tp(tp, "tp", "pdf", width = NULL) / 0.5,
                 convert_tp(tp, "tp", "pdf", width = 0.5),
                 info = "Width set to 0.5, must scale the CDF.")


# -------------------------------------------------------------------
# R implementation of calculating the CDF and PDF
# -------------------------------------------------------------------
fn <- function(i, tp, binwidth = 1) {
    tp <- tp[seq_len(i)] # Take first 'i' values
    # Fake bins: width 1
    bins <- seq(binwidth / 2, by = binwidth, length.out = i + 1)
    ynum <- (tail(bins, -1) + head(bins, -1)) / 2
    # Calculate numeric value falling into the last bin
    ynum <- mean(tail(bins, 2))
    # Convert numeric response to 'bin' representation (num2bin)
    y    <- transitreg:::num2bin(ynum, bins)
    res <- .Call("treg_predict_pdfcdf",
                 uidx = 42L, idx = rep(42L, length(tp)), tp = tp,
                 y = y, bins = bins, ncores = 1L, censored = "not-censored",
                 PACKAGE = "transitreg")
    return(data.frame(res))
}
tp <- convert_tp(pp, from = "cdf", to = "tp")
res <- do.call(rbind, lapply(seq_along(tp), fn, tp = tp))
expect_equal(res$pdf, dd)
expect_equal(res$cdf, pp)

# Calling transitreg_predict
bins <- seq(-0.5, by = 1, length.out = length(tp) + 1)
p3 <- .Call("treg_predict", uidx = 3L, idx = rep(3L, length(tp)), tp = tp,
            bins = bins, y = 0:15, prob = NA_real_, type = "cdf",
            cors = 1L, elementwise = FALSE, discrete = TRUE,
            censored = "not-censored")
expect_equal(pp, p3)
d3 <- .Call("treg_predict", uidx = 3L, idx = rep(3L, length(tp)), tp = tp,
            bins = bins, y = 0:15, prob = NA_real_, type = "pdf",
            cores = 1L, elementwise = FALSE, discrete = TRUE,
            censored = "not-censored")
expect_equal(dd, d3)



# -------------------------------------------------------------------
# Checking 'pmax'
# -------------------------------------------------------------------

# TODO(R): Update and check 'pmax', write tests.

