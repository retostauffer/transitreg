# -------------------------------------------------------------------
# Testing implementation distributions3 objects and methods
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("transitreg"))

# -------------------------------------------------------------------
# For testing, faking transition probabilities by drawing from
# a Poisson distribution and convert the CDF to TPs
# -------------------------------------------------------------------
# Converting CDF to transition probabilities

m <- rbind(A = convert_tp(ppois(0:15, lambda = 3.0), "cdf", "tp"),
           B = convert_tp(ppois(0:15, lambda = 8.0), "cdf", "tp"),
           C = convert_tp(ppois(0:15, lambda = 0.3), "cdf", "tp"))

# Fake breaks; width equal to 1
breaks <- seq(-0.5, by = 1, length.out = ncol(m) + 1)



# --------------------------------------------------------------------
# Testing constructor function
# --------------------------------------------------------------------

# ------------- distribution m[1, ] w/ fake data ---------------------

# One single distribution based on m[1, ]. Firstly, we test that
# the constructor function works and that the result is identical
# if we hand over a vector or a single-row matrix.
expect_silent(d1 <- Transition(as.vector(m[1, ]), breaks),
              info = "Calling constructor function (vector input)")
expect_silent(d1_m <- Transition(m[1, , drop = FALSE], breaks),
              info = "Constructor function with matrix input (nrow = 1)")
expect_identical(d1, setNames(d1_m, NULL),
               info = "Check that Transition is identical with vector/matrix input")
rm(d1_m)

# Testing the Transition object of length 1
expect_identical(class(d1), c("Transition", "distribution"), info = "Testing class")
expect_identical(length(d1), 1L,                             info = "Testing length")
pat <- "Transition_16\\([0-9\\.\\,\\ ]+)"
expect_stdout(print(d1), pattern = pat,                      info = "Testing standard representation")
expect_true(grepl(pat, format(d1)),                          info = "Testing format function")

# Checking attributes
expect_true("breaks" %in% names(attributes(d1)),             info = "Checking if attribute 'breaks' is available")
expect_identical(attr(d1, "breaks"), breaks,                 info = "Testing attribute 'breaks'")


expect_true("discrete" %in% names(attributes(d1)),           info = "Checking if attribute 'discrete' is available")
expect_identical(attr(d1, "discrete"), TRUE,                 info = "Testing attribute 'discrete'")
expect_identical(is_discrete(d1),   rep(TRUE, 1),            info = "Testing return of is_discrete")
expect_identical(is_continuous(d1), rep(FALSE, 1),           info = "Testing return of is_continuous")


# Convert to matrix
expect_silent(m1 <- as.matrix(d1),                          info = "Converting to matrix")
expect_identical(class(m1), c("Transitionmatrix", "matrix", "array"), info = "Testing matrix class")
expect_true(is.matrix(m1) && is.numeric(m1),                info = "Testing matrix class")
expect_identical(dim(m1), c(1L, ncol(m)),                   info = "Testing matrix dimension")
expect_true(all(grep("^tp_[0-9]+$", colnames(m1))),         info = "Matrix column names")
expect_true(is.null(rownames(m1)),                          info = "Matrix row names (unnamed)")
expect_identical(attr(m1, "breaks"), breaks,                info = "Testing 'breaks' attribute on matrix")
expect_identical(attr(m1, "discrete"), TRUE,                info = "Testing 'discrete' attribute on matrix")
rm(m1)

# Convert to extended (long) matrix
expect_silent(m1e <- as.matrix(d1, expand = TRUE),          info = "Converting to extended matrix")
expect_identical(class(m1e), c("Transitionmatrix", "matrix", "array"), info = "Testing matrix class")
expect_true(is.matrix(m1e) && is.numeric(m1e),              info = "Testing extended matrix class")
expect_identical(dim(m1e), c(ncol(m), 3L),                  info = "Testing extended matrix dimension")
expect_identical(colnames(m1e), c("index", "theta", "tp"),  info = "Matrix extended column names")
expect_true(is.null(rownames(m1e)),                         info = "Extended matrix row names (unnamed)")
expect_identical(attr(m1e, "breaks"), breaks,               info = "Testing 'breaks' attribute on matrix")
expect_identical(attr(m1e, "discrete"), TRUE,               info = "Testing 'discrete' attribute on matrix")

# Testing content ...
expect_true(all(m1e[, "index"] == 1L),                      info = "Checking matrix content (index)")
expect_identical(m1e[, "theta"], seq_len(ncol(m)) - 1.0,    info = "Checking matrix content (index)")
expect_equal(m1e[, "tp"], as.vector(m[1, ]),                info = "Checking matrix content (tp)")

# Testing that attributes survive when subsetting
# TODO(R): Not properly implemented, attemted it, but that broke
#          the format function as I need to support to access
#          specific elements using y[i] if dim(y) = c(1, N).
expect_true(is.null(attr(m1e[1, ], "breaks")))
expect_true(is.null(attr(m1e[1, ], "discrete")))
#expect_identical(attr(m1e[1, ], "breaks"), breaks,          info = "Testing 'breaks' attribute on matrix subset")
#expect_identical(attr(m1e[1, ], "discrete"), TRUE,          info = "Testing 'discrete' attribute on matrix subset")

rm(m1e)
rm(d1)


# ------------- distribution m, length 3 -----------------------------
# --------- testing non-discrete mode (breaks) -----------------------
# ---------------- continuous distribution ---------------------------

cbreaks <- seq(11.5, 34.1, length.out = length(breaks))
expect_silent(d1c  <- Transition(m, breaks = cbreaks),
              info = "Create continuous distribution, discrete = NULL")
expect_silent(d1c2 <- Transition(m, breaks = cbreaks, discrete = FALSE),
              info = "Create continuous distribution, discrete = FALSE")
expect_identical(d1c, d1c2, info = "Test that both objects are identical")
expect_true("breaks" %in% names(attributes(d1c)),         info = "Testing that attribute 'breaks' is available")
expect_identical(attr(d1c, "breaks"), cbreaks,            info = "Testing values of attribute 'breaks'")
expect_true("discrete" %in% names(attributes(d1c)),       info = "Testing that attribute 'discrete' is available")
expect_identical(attr(d1c, "discrete"), FALSE,            info = "Testing attribute 'discrete'")
expect_identical(is_continuous(d1c), rep(TRUE, nrow(m)),  info = "Testing return of is_continuous")
expect_identical(is_discrete(d1c),   rep(FALSE, nrow(m)), info = "Testing return of is_discrete")
rm(d1c, d1c2)

# ------------- distribution m, length 3 -----------------------------

# In contrast to the test above we name the distribution "A";
# should result in proper rownames on the matrices.
expect_silent(d3 <- Transition(m, breaks),                   info = "Calling constructor with three distributions")
expect_identical(class(d3), c("Transition", "distribution"), info = "Testing return class")
expect_identical(length(d3), nrow(m),                        info = "Testing length")
expect_identical(names(d3), rownames(m),                     info = "Testing names")
expect_identical(attr(d3, "breaks"), breaks,                 info = "Testing 'breaks' attribute on Transition")
expect_true(all(grepl(sprintf("^Transition_%d\\([0-9\\.\\,\\ ]+\\)$", ncol(m)), format(d3))), info = "Format")

# Convert to matrix
expect_silent(m3 <- as.matrix(d3),                           info = "Converting to matrix")
expect_true(is.matrix(m3) && is.numeric(m3),                 info = "Testing matrix class")
expect_identical(dim(m3), dim(m),                            info = "Testing matrix dimension")
expect_true(all(grep("^tp_[0-9]+$", colnames(m3))),          info = "Matrix column names")
expect_identical(rownames(m3), rownames(m),                  info = "Matrix row names")
expect_identical(attr(m3, "breaks"), breaks,                     info = "Testing 'breaks' attribute on matrix")

# Convert to extended (long) matrix
expect_silent(m3e <- as.matrix(d3, expand = TRUE),           info = "Converting to extended matrix")
expect_identical(class(m3e), c("Transitionmatrix", "matrix", "array"), info = "Testing matrix class")
expect_true(is.matrix(m3e) && is.numeric(m3e),               info = "Testing extended matrix class")
expect_identical(dim(m3e), c(3L * ncol(m), 3L),              info = "Testing extended matrix dimension")
expect_identical(colnames(m3e), c("index", "theta", "tp"),   info = "Matrix extended column names")
expect_true(is.null(rownames(m3e)),                          info = "Extended matrix row names (unnamed)")
expect_identical(attr(m3e, "breaks"), breaks,                    info = "Testing 'breaks' attribute on matrix")

# Testing content ...
expect_equal(m3e[, "index"], rep(seq_along(d3), each = ncol(m)),           info = "Checking matrix content (index)")
expect_identical(m3e[, "theta"], rep(seq_len(ncol(m)) - 1.0, length(d3)),  info = "Checking matrix content (index)")
expect_equal(m3e[, "tp"], as.vector(t(m)),                                 info = "Checking matrix content (tp)")
rm(m3e)
rm(d3)



# ------------- testing S3 methods for Transition --------------------

d3 <- Transition(m, breaks)

# Testing c(); We have already tested 'd3', so if the result
# of the c(...) call is identical we know the combined object is
# as it must be.
expect_identical(c(d3[1], d3[2], d3[3]), d3,                  info = "Combine Transition objects [1:3]")
expect_identical(c(d3[1], d3[2]), d3[1:2],                    info = "Combine Transition objects [1:2]")


expect_silent(d3df <- as.data.frame(d3),                      info = "Coerce to data.frame")
expect_true(is.data.frame(d3df),                              info = "Resulting class")
expect_identical(dim(d3df), c(length(d3), 1L),                info = "Data.frame dimension")
expect_identical(names(d3df), "d3",                           info = "Data.frame column names")
expect_identical(rownames(d3df), names(d3),                   info = "Data.frame row names")


# S3 method support
expect_silent(s3 <- support(d3),                              info = "Calling S3 method support")
expect_identical(s3, matrix(range(breaks), byrow = TRUE, ncol = 2L, nrow = length(d3), dimnames = list(names(d3), c("min", "max"))),
                 info = "Testing return of 'support()' method")


# S3 method is_discrete
expect_silent(is_discrete(d3),                                info = "Calling S3 method is_discrete")
expect_identical(is_discrete(d3), rep(TRUE, length(d3)),      info = "Testing return of is_discrete")
expect_silent(is_continuous(d3),                              info = "Calling S3 method is_continuous")
expect_identical(is_continuous(d3), rep(FALSE, length(d3)),   info = "Testing return of is_continuous")

# S3 method cdf
# Evaluate the distributions at all breaks (binmid)
binmid <- (head(breaks, -1) + tail(breaks, -1)) / 2
expect_silent(xcdf <- cdf(d3[1], binmid),                     info = "Calling S3 method cdf")
expect_identical(xcdf, convert_tp(m[1, ], "tp", "cdf"),       info = "Compare result from C/R")

# S3 method pdf
# Evaluate the distributions at all breaks (binmid)
i <- 1
expect_silent(dcdf <- cdf(d3[i], binmid),                     info = "Calling S3 method cdf")
expect_identical(dcdf, convert_tp(m[i, ], "tp", "cdf"),       info = "Compare result from C/R")
expect_silent(dpdf <- pdf(d3[i], binmid),                     info = "Calling S3 method pdf")
expect_identical(dpdf, convert_tp(m[i, ], "tp", "pdf"),       info = "Compare result from C/R")
i <- 2
expect_silent(dcdf <- cdf(d3[i], binmid),                     info = "Calling S3 method cdf")
expect_identical(dcdf, convert_tp(m[i, ], "tp", "cdf"),       info = "Compare result from C/R")
expect_silent(dpdf <- pdf(d3[i], binmid),                     info = "Calling S3 method pdf")
expect_identical(dpdf, convert_tp(m[i, ], "tp", "pdf"),       info = "Compare result from C/R")
i <- 3
expect_silent(dcdf <- cdf(d3[i], binmid),                     info = "Calling S3 method cdf")
expect_identical(dcdf, convert_tp(m[i, ], "tp", "cdf"),       info = "Compare result from C/R")
expect_silent(dpdf <- pdf(d3[i], binmid),                     info = "Calling S3 method pdf")
expect_identical(dpdf, convert_tp(m[i, ], "tp", "pdf"),       info = "Compare result from C/R")
rm(dpdf, dcdf)


# S3 method mean and median
dmean <- mean(d3)
expect_warning(dmean <- mean(d3),                             info = "Calling S3 method mean")
expect_true(is.double(dmean) && length(dmean) == length(d3),  info = "Checking return object")
expect_identical(names(dmean), names(d3),                     info = "Testing if names were carried along")
rm(dmean)

expect_silent(dmedian <- median(d3),                              info = "Calling S3 method median")
expect_true(is.double(dmedian) && length(dmedian) == length(d3),  info = "Checking return object")
expect_identical(names(dmedian), names(d3),                       info = "Testing if names were carried along")
expect_identical(median(d3), quantile(d3, 0.5),                   info = "Check that median == quantile_0.5")
rm(dmedian)


# S3 method quantile
# Just some quick tests for now
expect_silent(q <- quantile(d3, 0.5))
expect_true(is.double(q) && length(q) == length(d3) && all(names(q) == names(d3)))
rm(q)

expect_silent(q <- quantile(d3, c(0.1, 0.2, 0.3)))
expect_true(is.double(q) && length(q) == length(d3) && all(names(q) == names(d3)))
rm(q)

expect_silent(q <- quantile(d3, 0.5, elementwise = FALSE))
expect_true(is.double(q) && is.matrix(q) && identical(dim(q), c(3L, 1L)))
expect_true(identical(rownames(q), names(d3)) && identical(colnames(q), "50%"))
rm(q)

expect_silent(q <- quantile(d3, c(0.1, 0.2, 0.35), elementwise = FALSE))
expect_true(is.double(q) && is.matrix(q) && identical(dim(q), c(3L, 3L)))
expect_true(identical(rownames(q), names(d3)) && identical(colnames(q), c("10%", "20%", "35%")))
rm(q)

