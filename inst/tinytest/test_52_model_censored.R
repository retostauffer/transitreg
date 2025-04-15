# -------------------------------------------------------------------
# Testing 'transitreg' for discrete response (count data).
# -------------------------------------------------------------------

if (interactive()) { library("tinytest"); library("transitreg") }

# Drawing some random values in [0.5, 1.5, ..., 9.5] for testing
# simple censoring (checking caluclations are correct).
set.seed(111)
x <- pmin(rpois(500, lambda = 3), 9) + 0.5
x <- sample(c(rep(0.5, 200), x)) # adding a lot of 0.5's

# Estinmating uncensored model w/ theta0
f1  <- x ~ theta0 + s(theta, k = 8)
bk1 <- seq(0.0, 10.0, by = 1.0)
expect_silent(m1 <- transitreg(f1, breaks = bk1),
              info = "Estmating uncensored model")
expect_identical(m1$censored, "uncensored",
              info = "Checking element $censored.")

# Estimating left censored model (w/ theta0, censoring point)
f2  <- x ~ theta0 + s(theta, k = 8)
bk2 <- bk1[-1] # All but the first element of 'bk1'
expect_error(transitreg(f2, breaks = bk1, censored = "left"),
             pattern = "Formula contains 'theta0', but no observation falls into this bin.*",
             info = "Using wrong breaks (bin outside censoring range).")
expect_silent(m2 <- transitreg(f2, breaks = bk2, censored = "left"),
              info = "Estmating left censored model")
expect_identical(m2$censored, "left",
              info = "Checking element $censored.")

# Estimating left-right censored model, but without theta9 (right censoring
# point) to check we get the same results; though not the 'correct' formula.
f3  <- x ~ theta0 + s(theta, k = 8)
bk3 <- bk2[-length(bk2)] # All but the last element of 'bk2'
expect_warning(m3 <- transitreg(f3, breaks = bk3, censored = "both"),
               pattern = "consider adding 'theta9' to your formula when using 'censored = \"both\"",
               info = "Estimating left-right censored model, warning expected.")
expect_identical(m3$censored, "both",
              info = "Checking element $censored.")

# -------------------------------------------------------------------
# Testing logLik; should be the very same for all models
# -------------------------------------------------------------------
expect_identical(logLik(m1), logLik(m2),
              info = "Checking logLik, identical (although censoring 'uncensored'/'left').")
expect_identical(logLik(m1), logLik(m3),
              info = "Checking logLik, identical (although censoring 'uncensored'/'both').")

expect_identical(m1$probs, m2$probs,
              info = "Checking probs, identical (although censoring 'uncensored'/'left').")
expect_identical(m1$probs, m3$probs,
              info = "Checking probs, identical (although censoring 'uncensored'/'both').")

# -------------------------------------------------------------------
# Prediction (for random dummy) should also be identical, testing
# transition probabilities here.
# -------------------------------------------------------------------
expect_silent(p1 <- predict(m1, newdata = data.frame(x = "dummy"), type = "tp"),
              info = "Predicting `tp` for uncenosred model.")
expect_silent(p2 <- predict(m2, newdata = data.frame(x = "dummy"), type = "tp"),
              info = "Predicting `tp` for left censored model.")
expect_silent(p3 <- predict(m3, newdata = data.frame(x = "dummy"), type = "tp"),
              info = "Predicting `tp` for left-right censored model.")

expect_identical(p1, p2,
              info = "Testing prediction censored 'uncensored' vs 'left'.")
expect_identical(p1, p3,
              info = "Testing prediction censored 'uncensored' vs 'both'.")

# The transition probability for the first bin (i.e., the transition
# probability to be in bin 1 or larger/above left censoring point(
# should be equal to the empirical fraction of observations in that
# bin/the censoring point, which is nothing else than mean(x > 1).
emp <- mean(x > 1)
expect_equal(p1[1L], emp,
    info = "Checking first transition probability vs. empirical value (uncensored model).")
expect_equal(p2[1L], emp,
    info = "Checking first transition probability vs. empirical value (left censored model).")
expect_equal(p3[1L], emp,
    info = "Checking first transition probability vs. empirical value (left-right censored model).")
