stump <- function(x, y, weights) {
  foo <- function(split) {
    mL <- mean(y[x <= split])
    mR <- mean(y[x > split])
    fit <- rep(mL, length(y))
    fit[x > split] <- mR
    err <- sum((y - fit)^2)
    return(err)
  }

  split <- optimize(foo, lower = min(x), upper = max(x))$minimum

  mL <- mean(y[x <= split])
  mR <- mean(y[x > split])
  fit <- rep(mL, length(y))
  fit[x > split] <- mR

  return(fit)
}

set.seed(1328)
n <- 1000
x <- sort(runif(n, -3, 3))
y <- sin(x) + rnorm(n, sd = 0.3)

plot(x, y)

e <- y
fit <- 0

for(i in 1:1000) {
  fit <- fit + 0.1 * stump(x, e)
  e <- y - fit
  lines(fit ~ x, lwd = 0.2, col = 2)
}

lines(fit ~ x, lwd = 3, col = 4)


