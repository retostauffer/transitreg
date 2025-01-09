# -------------------------------------------------------------------
# Testing implementation distributions3 objects and methods
# -------------------------------------------------------------------

rm(list = objects())

suppressPackageStartupMessages(library("TransitionModels"))
if (interactive()) library("tinytest")

library("devtools")


# -------------------------------------------------------------------
# For testing, faking transition probabilities by drawing from
# a binomial distribution
# -------------------------------------------------------------------
# Converting CDF to transition probabilities
p2tp <- function(x) {
    tp <- numeric(length(x))
    tp[1] = 1 - x[1]
    prod <- 1
    for (i in seq.int(2, length(x))) {
        prod <- prod * tp[i - 1]
        tp[i] = (x[i - 1] - x[i] + prod) / prod
    }
    return(tp)
}

# First fake distribution. A data.frame, though
# the order of the columns is 'off'.
fake1 <- data.frame(tp = p2tp(pbinom(0:9, size = 10, prob = 0.10)),
                    lo = 0:9 - 0.5,
                    up = 0:9 + 0.5)
fake1 <- fake1[, c(3, 1, 2)] # Mixing

# Second fake distribution, different binning, different length.
# Just an unnamed list.
fake2 <- list(p2tp(pbinom(1:15, size = 15, prob = 0.20)), # transition prob
              1:15 * 2 - 3 - 1, # lower bound of bin
              1:15 * 2 - 3 + 1) # upper bound of bin




# --------------------------------------------------------------------
# Convert distribution to 'tmdist' object, testing return.
# --------------------------------------------------------------------

exp_colnames <- c("index", "tp", "lo", "up")

# ------------- distribution d1 (based on fake1) ---------------------
expect_silent(d1 <- tmdist(fake1), info = "Testing conversion")
expect_inherits(d1, "tmdist", info = "Testing return class")
expect_identical(length(d1), 1L, info = "Testing length")


# Convert to matrix
expect_silent(m1 <- as.matrix(d1),                          info = "Converting to matrix")
expect_true(is.matrix(m1) && is.numeric(m1),                info = "Testing matrix class")
expect_identical(dim(m1), c(1L, 3 * nrow(fake1)),           info = "Testing matrix dimension")
expect_true(all(grep("^(tp|lo|up)_[0-9]+$", colnames(m1))), info = "Matrix column names")
expect_true(is.null(rownames(m1)),                          info = "Matrix row names (unnamed)")
rm(m1)

# Convert to extended (long) matrix
expect_silent(m1e <- as.matrix(d1, expand = TRUE),           info = "Converting to extended matrix")
expect_true(is.matrix(m1e) && is.numeric(m1e),               info = "Testing extended matrix class")
expect_identical(dim(m1e), c(nrow(fake1), 4L),               info = "Testing extended matrix dimension")
expect_true(identical(exp_colnames, colnames(m1e)),          info = "Matrix extended column names")
expect_true(is.null(rownames(m1e)),                          info = "Extended matrix row names (unnamed)")

rm(m1e)
rm(d1)


# ------------- distribution d2 (based on fake2) ---------------------

# In contrast to the test above we name the distribution "A";
# should result in proper rownames on the matrices.
expect_silent(d2 <- tmdist(fake2), info = "Testing conversion")
expect_silent(names(d2) <- "A",    info = "Adding name")
expect_identical(names(d2), "A",   info = "Testing added name")
expect_inherits(d2, "tmdist",      info = "Testing return class")
expect_identical(length(d2), 1L,   info = "Testing length")


# Convert to matrix
expect_silent(m2 <- as.matrix(d2),                          info = "Converting to matrix")
expect_true(is.matrix(m2) && is.numeric(m2),                info = "Testing matrix class")
expect_identical(dim(m2), c(1L, 3 * nrow(fake1)),           info = "Testing matrix dimension")
expect_true(all(grep("^(tp|lo|up)_[0-9]+$", colnames(m2))), info = "Matrix column names")
expect_identical(rownames(m2), "A",                         info = "Matrix row name")
rm(m2)

# Convert to extended (long) matrix
expect_silent(m2e <- as.matrix(d2, expand = TRUE),           info = "Converting to extended matrix")
expect_true(is.matrix(m2e) && is.numeric(m2e),               info = "Testing extended matrix class")
expect_identical(dim(m2e), c(nrow(fake1), 4L),               info = "Testing extended matrix dimension")
expect_true(identical(exp_colnames, colnames(m2e)),          info = "Matrix extended column names")
expect_true(all(grepl("^A_[0-9]+$", rownames(m2e))),         info = "Extended matrix row names")

rm(m2e)
rm(d2)


# ------------- distribution d1 and d2 -------------------------------
# This results in a tmdist object of length 2; we're also adding
# row names for testing. As one distribution is longer than the other,
# the non-extended matrix will contain a series of missing values.
# In the extended version, they should be gone, and the number of
# rows corresponds to the sum of the length of both distributions
# (defined by fake1 and fake2).

expect_silent(d3 <- tmdist(list(fake1, fake2)), info = "Testing conversion (lenght 2)")
expect_inherits(d3, "tmdist", info = "Testing return class")
expect_identical(length(d3), 2L, info = "Testing length")
expect_silent(names(d3) <- c("FOO", "BAR"),     info = "Adding names (length 2)")
expect_identical(c("FOO", "BAR"), names(d3),    info = "Testing added names")


# Convert to matrix
expect_silent(m3 <- as.matrix(d3),                          info = "Converting to matrix")
expect_true(is.matrix(m3) && is.numeric(m3),                info = "Testing matrix class")
expect_identical(dim(m3), c(length(d3), 3L * max(nrow(fake1), length(fake2[[1]]))), info = "Testing matrix dimension")
expect_true(all(grep("^(tp|lo|up)_[0-9]+$", colnames(m3))), info = "Matrix column names")
expect_identical(c("FOO", "BAR"), rownames(m3),             info = "Matrix row names")

# As the two fake distributions are not identical in length,
# we are expecting a few missing values in the matrix.
exp_na <- abs(nrow(fake1) - length(fake2[[1]])) * 3L
expect_identical(sum(is.na(m3)), exp_na,                    info = "Checking number of expected missing values")

rm(m3)

# Convert to extended (long) matrix
expect_silent(m3e <- as.matrix(d3, expand = TRUE),           info = "Converting to extended matrix")
expect_true(is.matrix(m3e) && is.numeric(m3e),               info = "Testing extended matrix class")
expect_identical(dim(m3e), c(nrow(fake1) + length(fake2[[1]]), 4L), info = "Testing extended matrix dimension")
expect_true(identical(exp_colnames, colnames(m3e)),          info = "Matrix extended column names")
expect_true(all(grepl("^(FOO|BAR)_[0-9]+$", rownames(m3e))), info = "Extended matrix row names (unnamed)")

rm(m3e)
rm(d3)





# --------------------------------------------------------------------
# Testing other S3 methods, using the same fake distributions from above.
# --------------------------------------------------------------------
d1 <- tmdist(fake1)
d3 <- setNames(tmdist(list(fake1, fake2)), LETTERS[1:2])

# Just to ensure the two objects are available:
expect_inherits(d1, "tmdist"); expect_identical(length(d1), 1L)
expect_inherits(d3, "tmdist"); expect_identical(length(d3), 2L)
expect_identical(names(d3), c("A", "B"))

# Converting to data.frame
expect_silent(df1 <- as.data.frame(d1))
expect_inherits(df1, "data.frame")
expect_identical(names(df1), "d1")
expect_identical(dim(df1), c(1L, 1L))
expect_identical(rownames(df1), "1")


expect_silent(df3 <- as.data.frame(d3))
expect_inherits(df3, "data.frame")
expect_identical(names(df3), "d3")
expect_identical(dim(df3), c(2L, 1L))
expect_identical(rownames(df3), names(d3))


