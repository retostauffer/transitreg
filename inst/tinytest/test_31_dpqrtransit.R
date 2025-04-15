# -------------------------------------------------------------------
# Testing implementation distributions3 objects and methods
# -------------------------------------------------------------------

if (interactive()) { library("tinytest"); library("transitreg") }

# -------------------------------------------------------------------
# For testing, faking transition probabilities by drawing from
# a Poisson distribution and convert the CDF to TPs
# -------------------------------------------------------------------
# Converting CDF to transition probabilities

N <- 10
d1 <- Transition(convert_tp(ppois(0:N, 3), "cdf", "tp"), seq.int(0, N+1))
d2 <- Transition(convert_tp(ppois(0:N, 3), "cdf", "tp"), as.double(seq.int(0, N + 1)))

# RETO(TODO): Just for testing; remove in the future
#k <- rtransit(500, d1); hist(k)
#k <- plot(d1, cdf = F, plot = TRUE, tp = TRUE)
#x <- seq(0, 5, by = 0.5)
#dtransit(x, d1)
#ptransit(x, d1)
#qtransit(seq(0.1, 0.9, by = 0.1), d1)


# Drawing from a Poisson distribution to check if the
# dpqr functions do what they promis.
N <- 10
true_cdf <- ppois(0:N, lambda = 3.2)
true_pdf <- dpois(0:N, lambda = 3.2)
breaks_int <- seq.int(0, N+1)
breaks_dbl <- as.double(seq.int(0, N+1))

# Setting up two distributions, one with integer breaks ('count data')
# and one with double breaks (pseudo-binned continuous distribution).
dint <- Transition(convert_tp(true_cdf, "cdf", "tp"), breaks_int)
ddbl <- Transition(convert_tp(true_cdf, "cdf", "tp"), breaks_dbl)


# --------------- function misuse ---------------------------------------
expect_error(dtransit(),                   info = "All inputs missing")
expect_error(dtransit(d = dint),           info = "Argument 'x' missing")
expect_error(dtransit(x = breaks_int),     info = "Argument 'd' missing")
expect_error(dtransit(vector("double"), dint),
                pattern = "'x' must be numeric of length > 0",
                info = "Incorrect argument on 'x'.")
expect_error(dtransit(LETTERS[1:5], dint),
                pattern = "'x' must be numeric of length > 0",
                info = "Incorrect argument on 'x'.")
expect_error(dtransit(c(1, 2, NA),, dint),
                pattern = "missing values in 'x' not allowed",
                info = "Missing data in 'x' not allowed.")
expect_error(dtransit(breaks_int, dint, log = c("foo", "bar")),
                pattern = "'log' must evaluate to TRUE/FALSE",
                info = "Incorrect argument on 'log'.")
expect_error(dtransit(breaks_int, dint, log = "YES PLEASE"),
                pattern = "'log' must evaluate to TRUE/FALSE",
                info = "Incorrect argument on 'log'.")
expect_error(dtransit(breaks_int, dint, log = vector("logical")),
                pattern = "'log' must evaluate to TRUE/FALSE",
                info = "Incorrect argument on 'log'.")
expect_error(dtransit(breaks_int, dint, log = list()),
                pattern = "'list' object cannot be coerced to type 'logical'",
                info = "Incorrect argument on 'log'.")







# Testing dtransit
library("tinytest")
expect_equal(dtransit(breaks_int, dint), true_pdf)
expect_equal(dtransit(breaks_dbl, dint), true_pdf)
expect_equal(ptransit(breaks_int, dint), true_cdf)
expect_equal(ptransit(breaks_dbl, dint), true_cdf)

expect_equal(dtransit(breaks_int, ddbl), true_pdf)
expect_equal(dtransit(breaks_dbl, ddbl), true_pdf)
expect_equal(ptransit(breaks_int, ddbl), true_cdf)
expect_equal(ptransit(breaks_dbl, ddbl), true_cdf)



# Check if the distribution of N random values matches the true PDF
# Note: The seed is important to ensure we draw random values from all
# bins, even if the probability is close to zero.
set.seed(112)
expect_silent(r <- rtransit(3000, dint))
tmp <- unname(as.vector(proportions(table(r))))
expect_true(all(abs(tmp - true_pdf) < 0.02))


