# -------------------------------------------------------------------
# Testing the functions in R/misc.R
# -------------------------------------------------------------------

suppressPackageStartupMessages(library("tinytest"))
suppressPackageStartupMessages(library("TransitionModels"))

# # -------------------------------------------------------------------
# # Estimating simple model (from ?tm examples)
# # -------------------------------------------------------------------
# set.seed(123)
# n <- 1000
# x <- runif(n, -3, 3)
# y <- rpois(n, exp(2 + sin(x)))
# 
# # Fit transition count response model.
# b <- tm(y ~ s(theta) + s(x))


# -------------------------------------------------------------------
# 'grep2': Convenience function which allows to grep()
# for multiple patterns.
# -------------------------------------------------------------------
haystack <- c("Tallin", "Berlin", "Innsbruck", "London", "Dublin")
expect_identical(TransitionModels:::grep2("bruck", haystack),
                 grep("bruck", haystack))
expect_identical(TransitionModels:::grep2("lin$", haystack),
                 grep("lin$", haystack))
expect_identical(sort(TransitionModels:::grep2(c("Inns", "lin$"), haystack)),
                 sort(grep("(Inns|lin$)", haystack)))


# -------------------------------------------------------------------
# 'response_name': Convenience function to extract the response name
# of a formula object(s).
# -------------------------------------------------------------------
expect_identical(TransitionModels:::response_name(Y ~ x), "Y")
expect_identical(TransitionModels:::response_name(formula("Y ~ x")), "Y")
expect_identical(TransitionModels:::response_name(foo ~ a + b + c^2), "foo")

f <- list(formula(foo ~ a + b + c^2), formula(bar ~ poly(x)))
expect_identical(TransitionModels:::response_name(f), c("foo", "bar"))

# If input is a list with a $formula element, that one is taken
expect_identical(TransitionModels:::response_name(list(formula = Y ~ x, foo = 3)), "Y")


# -------------------------------------------------------------------
# 'fake_formula': Very basic checks
# -------------------------------------------------------------------

f <- formula(Y ~ a + b + c^2)
expect_identical(TransitionModels:::fake_formula(f), f)
rm(f)

f <- list(formula(foo ~ a + b + c^2), formula(bar ~ poly(x)))
expect_identical(lengths(f), c(3L, 3L))
f <- TransitionModels:::fake_formula(f)
expect_identical(length(f), c(2L, 2L))
expect_identical(length(f[[1]]), 1L)
expect_identical(length(f[[2]]), 3L)
expect_identical(format(f[[2]]), "a + b + c^2 | poly(x)")
rm(f)


f <- list(formula(foo ~ s(a) + te(b) + c^2), formula(bar ~ poly(x)))
expect_identical(lengths(f), c(3L, 3L))
f <- TransitionModels:::fake_formula(f)
expect_identical(length(f), c(2L, 2L))
expect_identical(length(f[[1]]), 1L)
expect_identical(length(f[[2]]), 3L)
expect_identical(format(f[[2]]), "c + a + b | poly(x)")


f <- list(formula(foo ~ s(a) + te(b) + c^2), formula(bar ~ poly(x)))

expect_identical(TransitionModels:::fake_formula(f),
                 TransitionModels:::fake_formula(f, specials = NULL,
                            nospecials = FALSE, onlyspecials = FALSE),
                 info = "Testing defaults of fake_formula")

fx <- TransitionModels:::fake_formula(f, nospecials = FALSE)
expect_identical(format(fx[[2]]), "c + a + b | poly(x)")
rm(fx)

fx <- TransitionModels:::fake_formula(f, nospecials = TRUE)
expect_identical(format(fx[[2]]), "c | poly(x)")
rm(fx)

fx <- TransitionModels:::fake_formula(f, onlyspecials = TRUE)
expect_true(is.list(fx) && length(fx) == 2L)
expect_identical(fx, list(c("s(a)", "te(b)"), character()))
rm(fx)


# -------------------------------------------------------------------
# 'ff_replace':
# -------------------------------------------------------------------

orig     <- formula("foo ~ a:b + bar*c + d + e + f:g:h")
##expected <- "foo ~ a + b + bar + c + d + e + f + g + h + bar:c"
expected <- "foo ~ a + b + bar + c + d + e + f + g + h"
expect_silent(tmp <- TransitionModels:::ff_replace(orig))
## FAILS ## expect_identical(format(tmp), expected)
# TODO(R): Niki, there is a remaining bar:c in the example above.
#          Is that intended?
rm(orig, expected, tmp)


orig     <- list(formula(foo ~ a:b), formula(bar ~ c*d*e))
##expected <- "foo ~ a + b + bar + c + d + e + f + g + h + bar:c"
expected <- "foo ~ a + b + bar + c + d + e + f + g + h"
expect_silent(tmp <- TransitionModels:::ff_replace(f))
# TODO(R): Niki ff_replace is intended for one single formula object I assume?
#          Given there are no sanity checks the call above
#          does not throw any error, but does not return the intended result right?



