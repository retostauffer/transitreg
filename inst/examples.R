library("TransitionModels")
library("gamlss2")

## Load the 'storms' dataset from the 'WeatherGermany' package and fit a transition model.
data("storms", package = "WeatherGermany")
f <- counts ~ ti(alt) + ti(year) + ti(lon, lat) + ti(alt, year)
b <- tm(f, data = storms)

## Generate random count data and visualize its probability density function (PDF).
set.seed(123)
n <- 3000
y <- rpois(n, 10)
tm_dist(y)

## Fit a transition model to generated Poisson data and visualize predictions.
set.seed(123)
n <- 500
x <- runif(n, -3, 3)
y <- rpois(n, exp(2 + sin(x)))
b <- tm(y ~ s(x))

nd <- data.frame("x" = seq(-3, 3, length = 100))
fit <- cbind(
  predict(b, nd, type = "pmax"),
  predict(b, nd, prob = 0.05 / 2),
  predict(b, nd, prob = 0.5),
  predict(b, nd, prob = 1 - 0.05 / 2)
)

plot(y ~ x)
matplot(nd$x, fit, type = "l", lty = 1, lwd = 3, col = c(2, 4, 4, 4), add = TRUE)
topmodels::rootogram(b)

## Fit a transition model for continuous response data with a discretization step.
set.seed(123)
n <- 1000
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = exp(-1 + cos(x)))
b <- tm(y ~ s(x) + te(x, theta, k = 10), breaks = 201)

nd <- data.frame("x" = seq(min(x), max(x), length = 100))
fit <- cbind(
  predict(b, nd, type = "pmax"),
  predict(b, nd, prob = 0.05),
  predict(b, nd, prob = 0.5),
  predict(b, nd, prob = 0.95)
)

plot(y ~ x, xlim = range(nd$x), ylim = range(c(fit, y)))
matplot(nd$x, fit, type = "l", lty = 1, lwd = 3, col = c(2, 4, 4, 4), add = TRUE)

## Compare predictions from transition models and generalized additive models.
m <- gamlss2(y ~ s(x, k = 40) | s(x, k = 40), family = NO)

par <- predict(m, newdata = nd, drop = FALSE)
fit <- cbind(
  family(m)$q(0.05, par),
  family(m)$q(0.5, par),
  family(m)$q(0.95, par)
)

matplot(nd$x, fit, type = "l", lty = 1, lwd = 3, col = 3, add = TRUE)

## Load the 'rent' dataset, split it into training and testing data, and fit models.
data("rent", package = "gamlss.data")
set.seed(123)
i <- sample(1:2, size = nrow(rent), replace = TRUE, prob = c(0.6, 0.4))
dtrain <- rent[i == 1, ]
dtest <- rent[i == 2, ]

f <- R ~ loc + s(Fl) + s(A)
b <- tm(f, data = dtrain, breaks = 300)

m <- gamlss2(R ~ s(Fl) + s(A) + loc | s(Fl) + s(A) + loc, data = dtrain, family = GA)

p <- predict(b, newdata = dtest, prob = 0.5)
par <- predict(m, newdata = dtest)
pm <- family(m)$q(0.5, par)

err1 <- mean((dtest$R - p)^2)
err2 <- mean((dtest$R - pm)^2)

print(err1)
print(err2)

plot(p, dtest$R, xlim = range(c(p, pm), na.rm = TRUE))
points(pm, dtest$R, col = 2)
abline(0, 1, lwd = 2, col = 4)

## Fit a transition model to the 'film90' dataset and visualize results.
data("film90", package = "gamlss.data")
f <- lborev1 ~ s(lboopen) + te(theta, lboopen, k = 10)
b <- tm(f, data = film90, breaks = 200, scale.x = TRUE)

nd <- data.frame("lboopen" = seq(min(film90$lboopen), max(film90$lboopen), length = 100))
fit <- cbind(
  predict(b, nd, type = "pmax"),
  predict(b, nd, prob = 0.05),
  predict(b, nd, prob = 0.5),
  predict(b, nd, prob = 0.95)
)

plot(lborev1 ~ lboopen, data = film90, xlim = range(nd$lboopen, na.rm = TRUE), ylim = range(c(fit, film90$lborev1), na.rm = TRUE))
matplot(nd$lboopen, fit, type = "l", lty = 1, lwd = 3, col = c(2, 4, 4, 4), add = TRUE)

## Fit a model for BMI data and visualize quantile predictions.
data("dbbmi", package = "gamlss.data")

b <- tm(bmi ~ s(sqrt(age)) + te(theta,sqrt(age)), data = dbbmi, breaks = 300)

nd <- data.frame("age" = seq(min(dbbmi$age), max(dbbmi$age), length = 300))
fit <- NULL
for (p in c(0.02, 0.10, 0.25, 0.50, 0.75, 0.90, 0.98)) {
  fit <- cbind(fit, predict(b, nd, prob = p))
}

plot(bmi ~ age, data = dbbmi, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.1))
matplot(nd$age, fit, type = "l", lty = 1, lwd = 2, col = 4, add = TRUE)

## Compare predictions between transition and GAMLSS models for categorized data.
set.seed(123)
n <- 1000
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = exp(-1 + cos(x)))
yc <- cut(y, breaks = 30, labels = FALSE, include.lowest = TRUE) - 1

b <- tm(yc ~ x, engine = "nnet", size = 10, scale.x = TRUE)
m <- gamlss2(yc ~ s(x) | s(x) | s(x), family = NBI)

nd <- data.frame("x" = seq(min(x), max(x), length = 100))
fit <- cbind(
  predict(b, nd, type = "pmax"),
  predict(b, nd, prob = 0.05),
  predict(b, nd, prob = 0.5),
  predict(b, nd, prob = 0.95)
)

par <- predict(m, newdata = nd, drop = FALSE)
fit2 <- cbind(
  family(m)$q(0.05, par),
  family(m)$q(0.5, par),
  family(m)$q(0.95, par)
)

plot(yc ~ x, xlim = range(nd$x), ylim = range(c(fit, yc)))
matplot(nd$x, fit, type = "l", lty = 1, lwd = 3, col = c(2, 4, 4, 4), add = TRUE)
matplot(nd$x, fit2, type = "l", lty = 1, lwd = 3, col = 3, add = TRUE)

## Fit a transition model to the 'airquality' dataset and compare predictions.
d <- na.omit(airquality)
b <- tm(Ozone ~ s(Temp) + s(Wind) + ti(theta, Temp) + ti(theta, Wind), data = d)
p <- predict(b, prob = 0.5)

plot(p, d$Ozone)
abline(0, 1)

