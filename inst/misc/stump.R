stump <- function(x, y, weights, split_vals) {
  foo <- function(split) {
    mL <- mean(y[x <= split])
    mR <- mean(y[x > split])
    fit <- rep(mL, length(y))
    fit[x > split] <- mR
    err <- sum((y - fit)^2)
    return(err)
  }

  # split <- optimize(foo, lower = min(x), upper = max(x))$minimum
  err <- numeric(length = length(split_vals))
  
  for(i in seq_along(split_vals)) {
    err[i] <- foo(split_vals[i])
  }

  split <- split_vals[which.min(err)]
  
  mL <- mean(y[x <= split])
  mR <- mean(y[x > split])
  fit <- rep(mL, length(y))
  fit[x > split] <- mR

  return(list(fit = fit, split = split))
}

set.seed(1328)
n <- 1000
x <- sort(runif(n, -3, 3))
y <- sin(x) + rnorm(n, sd = 0.3)

x11()
plot(x, y)

e <- y
fit <- 0

sv <- sort(unique(round(x, 4)))

for(i in 1:1000) {
  if(length(sv)) {
    s <- stump(x, e, split_vals = sv)
    sv <- sv[sv != s$split]
    fit <- fit + 1 * s$fit
    e <- y - fit
    lines(fit ~ x, lwd = 0.2, col = 2)
  } else {
    break
  }
}

points(x, y)
lines(fit ~ x, lwd = 3, col = 4)


