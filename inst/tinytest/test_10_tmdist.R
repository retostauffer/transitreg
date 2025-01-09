# -------------------------------------------------------------------
# Testing implementation distributions3 objects and methods
# -------------------------------------------------------------------

rm(list = objects())

suppressPackageStartupMessages(library("TransitionModels"))
if (interactive()) library("tinytest")


# -------------------------------------------------------------------
# For testing, faking transition probabilities by drawing from
# a binomial distribution
# -------------------------------------------------------------------
# Converting CDF to transition probabilities
p2tp <- TransitionModels:::cdf_to_tp

# First fake distribution. A data.frame, though
# the order of the columns is 'off'.
fake1 <- data.frame(tp = p2tp(pbinom(0:9, size = 10, prob = 0.10)),
                    lo = 0:9 - 0.5,
                    up = 0:9 + 0.5)
fake1 <- fake1[, c(3, 1, 2)] # Mixing

# Second fake distribution, different binning, different length.
# Just an unnamed list.
fake2 <- list(p2tp(pbinom(1:15, size = 15, prob = 0.20)), # transition prob
              1:15 * 2 - 2.9 - 1, # lower bound of bin; -2.9 -> not discrete (not integer)
              1:15 * 2 - 2.9 + 1) # upper bound of bin; -2.9 (see above)




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

expect_true(all(m1e[, "index"] == 1L),                       info = "Checking matrix content (index)")
expect_equal(m1e[, "tp"], fake1$tp,                          info = "Checking matrix content (tp)")
expect_equal(m1e[, "lo"], fake1$lo,                          info = "Checking matrix content (lo)")
expect_equal(m1e[, "up"], fake1$up,                          info = "Checking matrix content (up)")

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

expect_true(all(m2e[, "index"] == 1L),                       info = "Checking matrix content (index)")
expect_equivalent(m2e[, "tp"], fake2[[1]],                   info = "Checking matrix content (tp)")
expect_equivalent(m2e[, "lo"], fake2[[2]],                   info = "Checking matrix content (lo)")
expect_equivalent(m2e[, "up"], fake2[[3]],                   info = "Checking matrix content (up)")

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

tmp_index <- rep(1:2, c(nrow(fake1), length(fake2[[1]])))
tmp_tp    <- c(fake1$tp, fake2[[1]])
tmp_lo    <- c(fake1$lo, fake2[[2]])
tmp_up    <- c(fake1$up, fake2[[3]])
expect_equivalent(m3e[, "index"], tmp_index,                 info = "Checking matrix content (tp)")
expect_equivalent(m3e[, "tp"],    tmp_tp,                    info = "Checking matrix content (tp)")
expect_equivalent(m3e[, "lo"],    tmp_lo,                    info = "Checking matrix content (lo)")
expect_equivalent(m3e[, "up"],    tmp_up,                    info = "Checking matrix content (up)")

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

# Print
p = "Transition Dist \\(.*\\)"
expect_stdout(print(d1),
              pattern = p,
              info = "Testing standard representation (format)")
expect_stdout(print(d3),
              pattern = paste(p, p, sep = ".*"), # <- twice
              info = "Testing standard representation (format)")

# Converting to data.frame
expect_silent(df1 <- as.data.frame(d1),                  info = "Convert distributions to data.frame")
expect_inherits(df1, "data.frame",                       info = "distributions data.frame")
expect_identical(names(df1), "d1",                       info = "Testing name on distributions data.frame")
expect_identical(dim(df1), c(1L, 1L),                    info = "Testing dimension of distributions data.frame")
expect_identical(rownames(df1), "1",                     info = "Testing rownames of distributions data.frame")

expect_silent(df3 <- as.data.frame(d3),                  info = "Convert distributions to data.frame")
expect_inherits(df3, "data.frame",                       info = "distributions data.frame")
expect_identical(names(df3), "d3",                       info = "Testing name on distributions data.frame")
expect_identical(dim(df3), c(2L, 1L),                    info = "Testing dimension of distributions data.frame")
expect_identical(rownames(df3), names(d3),               info = "Testing rownames of distributions data.frame")



# Checking support as well as is_discrete, is_continuous
exp_support <- rbind(c(min = min(fake1$lo),   max = max(fake1$up)),
                     c(min = min(fake2[[2]]), max = max(fake2[[3]])))
expect_silent(support(d3),                               info = "support() should be silent")
expect_equal(support(d3), exp_support,                   info = "Checking support (return class, names, and values")

expect_silent(is_discrete(d3),                           info = "is_discrete() should be silent")
expect_silent(is_continuous(d3),                         info = "is_continuous() should be silent")
expect_identical(is_discrete(d3), c(TRUE, FALSE),        info = "Testing return of is_discrete")
expect_identical(is_discrete(d3), !c(TRUE, FALSE),       info = "Testing return of is_continuous")
expect_identical(is_discrete(d3), !is_continuous(d3),    info = "Must be complementary")





