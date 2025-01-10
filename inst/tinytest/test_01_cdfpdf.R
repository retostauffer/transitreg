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
expect_silent(tp <- TransitionModels:::cdf_to_tp(pp))

# Testing reverse
expect_silent(p2 <- TransitionModels:::tp_to_cdf(tp))
expect_equal(pp, p2)
expect_silent(d2 <- TransitionModels:::tp_to_pdf(tp))
expect_equal(dd, d2)

# -------------------------------------------------------------------
# R implementation of calculating the CDF and PDF
# -------------------------------------------------------------------
fn <- function(i, tp) {
    tp <- tp[seq_len(i)] # Take first 'i' values
    res <- .Call("tm_predict_pdfcdf", 42L, rep(42L, length(tp)), tp = tp, 1L,
                 PACKAGE = "TransitionModels")
    return(data.frame(res))
}
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

