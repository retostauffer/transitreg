library("TransitionModels")
library("gamlss2")
library("qgam")

## Pin loss function.
pinLoss <- function(u, tau = 0.5)
{
  u * (tau - 1 * (u < 0))
}

## CRPS function.
CRPS <- function(model, y, newdata)
{
  par <- NULL
  if(inherits(model, "gamlss2")) {
    par <- predict(model, newdata = newdata)
  }
  qu <- seq(0.01, 0.99, by = 0.01)
  err <- 0
  for(q in qu) {
    if(inherits(model, "gamlss2")) {
      pq <- model$family$q(q, par)
    }
    if(inherits(model, "tm")) {
      pq <- predict(model, newdata = newdata, prob = q)
    }
    err <- err + pinLoss(y - pq, q)
  }
  err <- err / length(qu)
  return(mean(err))
}

## Data generating processes.
dgp_NO <- function(n = 1000, probs = c(0.01, 0.1, 0.5, 0.9, 0.99), breaks = 20, ...)
{
  ## Covariate.
  x <- runif(n, -3, 3)

  ## Parameters.
  mu <- sin(x)
  sigma <- exp(-1 + cos(x))

  ## Response
  y_cont <- rnorm(n, mean = mu, sd = sigma)

  ## Discretize the normal data into count categories.
  yr <- range(y_cont)
  breaks <- seq(-4, 4, length.out = breaks)
  y_count <- cut(y_cont, breaks = breaks, labels = FALSE, include.lowest = TRUE) - 1

  ## Combine.
  d <- data.frame("x" = x, "counts" = y_count, "num" = y_cont)

  ## Add quantiles.
  qu <- quc <- NULL
  for(j in probs) {
    qj <- qNO(j, mu = mu, sigma = sigma)
    qu <- cbind(qu, qj)
    quc <- cbind(quc, cut(qj, breaks = breaks, labels = FALSE, include.lowest = TRUE) - 1)
  }
  colnames(qu) <- colnames(quc) <- paste0(probs * 100, "%")
  d$quantiles <- qu
  d$count.quantiles <- quc

  return(d)
}

dgp_ZINBI <- function(n = 1000, probs = c(0.01, 0.1, 0.5, 0.9, 0.99), ...)
{
  ## Covariate data.
  x <- runif(n, -2, 4)

  ## Parameters.
  mu <- exp(1 + 2 * sin(x))
  sigma <- exp(-0.5*x)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  nu <- sigmoid(-4 + cos(x))

  ## Sample response.
  y <- rZINBI(n, mu = mu, sigma = sigma, nu = nu)

  ## Combine.
  d <- data.frame("x" = x, "y" = y)

  ## Add quantiles.
  qu <- NULL
  for(j in probs) {
    qu <- cbind(qu, qZINBI(j, mu = mu, sigma = sigma, nu = nu))
  }
  colnames(qu) <- paste0(probs * 100, "%")
  d$quantiles <- qu

  return(d)
}

sim_NO <- function(n = 1000, breaks = NULL, counts = FALSE, family = NO, seed = 111,
  xlim = NULL, ylim = NULL, engine = "bam", pos = "topleft", ...)
{
  if(!is.null(seed))
    set.seed(seed)

  ## Simulate data.
  d <- dgp_NO(n, ...)

  ## Simulate new data.
  nd <- dgp_NO(1000, ...)

  if(counts) {
    d$num <- d$counts
    nd$num <- nd$counts
    d$quantiles <- d$count.quantiles
    nd$quantiles <- nd$count.quantiles
  }

  ## Set quantiles.
  qu <- colnames(nd$quantiles)
  qu <- as.numeric(gsub("%", "", qu, fixed = TRUE))/100

  ## Estimate transition model.
  if(engine != "nnet") {
    f <- num ~ s(x) + te(x,theta)
  } else {
    f <- num ~ x
  }
  b <- tm(f, data = d, breaks = breaks, engine = engine,
    scale.x = TRUE, size = 40, maxit = 1000, decay = 0.01)

  ## Corresponding GAMLSS.
  if(counts && (family()$type == "Continuous")) {
    m <- gamlss2(log(num + 1) ~ s(x) | s(x), family = family, data = d)
  } else {
    m <- gamlss2(num ~ s(x) | s(x), family = family, data = d)
  }

  ## Quantile regression.
  if(counts) {
    g <- mqgam(log(num + 1) ~ s(x), data = d, qu = qu)
  } else {
    g <- mqgam(num ~ s(x), data = d, qu = qu)
  }

  ## Predict quantiles.
  p <- do.call("cbind",
    lapply(qu, function(j) {
      predict(b, newdata = nd, prob = j)
  }))

  par <- predict(m, newdata = nd)
  pm <- do.call("cbind",
    lapply(qu, function(j) {
      m$family$q(j, par)
  }))
  if(counts && (family()$type == "Continuous"))
    pm <- exp(pm) - 1

  pg <- do.call("cbind",
    lapply(qu, function(j) {
      qdo(g, j, predict, newdata = nd)
  }))
  if(counts)
    pg <- exp(pg) - 1

  ## Compute loss.
  err_b <- sqrt(mean((p - nd$quantiles)^2))
  err_m <- sqrt(mean((pm - nd$quantiles)^2))
  err_g <- sqrt(mean((pg - nd$quantiles)^2))

  ## Plot data and fitted median.
  if(is.null(ylim))
    ylim <- range(d$num, nd$num, p, pm, pg)

  plot(num ~ x, data = d, xlim = xlim, ylim = ylim,
    col = rgb(0.1, 0.1, 0.1, alpha = 0.2), pch = 16,
    main = "", xlab = "x", ylab = "y", ...)

  i <- order(nd$x)

  matplot(nd$x[i], nd$quantiles[i,], type = "l",
    col = rgb(0.1, 0.1, 0.1, alpha = 0.3),
    lty = 1, lwd = 6, add = TRUE)

  matplot(nd$x[i], pm[i, ], type = "l",
    lty = 1, col = 2, lwd = 2, add = TRUE)
  matplot(nd$x[i], pg[i, ], type = "l",
    lty = 1, col = 3, lwd = 2, add = TRUE)
  matplot(nd$x[i], p[i, ], type = "l",
    lty = 1, col = 4, lwd = 2, add = TRUE)

  qe <- paste(c("TM", paste("GAMLSS", family()$family[1]), "QGAM"), "=", round(c(err_b, err_m, err_g), 2))

  legend(pos, qe, lwd = 2, col = c(4, 2, 3), bty = "n") ##title = "Quantile Error (RMSE)")
}

sim_ZINBI <- function(n = 1000, seed = 111, breaks = NULL,
  family = ZINBI, pos = "topleft", xlim = NULL, ylim = NULL, ...)
{
  if(!is.null(seed))
    set.seed(seed)

  ## Simulate data.
  d <- dgp_ZINBI(n, ...)

  ## Simulate new data.
  nd <- dgp_ZINBI(1000, ...)

  ## Set quantiles.
  qu <- colnames(nd$quantiles)
  qu <- as.numeric(gsub("%", "", qu, fixed = TRUE))/100

  ## Estimate transition model.
  f <- y ~ ti(theta) + ti(x) + ti(theta,x)
  b <- tm(f, data = d, breaks = breaks)

  ## Corresponding GAMLSS.
  m <- gamlss2(y ~ s(x) | s(x) | s(x), family = family, data = d)

  ## Quantile regression.
  g <- mqgam(log(y + 1) ~ s(x), data = d, qu = qu)

  ## Predict quantiles.
  p <- do.call("cbind",
    lapply(qu, function(j) {
      predict(b, newdata = nd, prob = j)
  }))

  par <- predict(m, newdata = nd)
  pm <- do.call("cbind",
    lapply(qu, function(j) {
      m$family$q(j, par)
  }))

  pg <- do.call("cbind",
    lapply(qu, function(j) {
      qdo(g, j, predict, newdata = nd)
  }))
  pg <- exp(pg) - 1

  ## Compute errors.
  err_b <- sqrt(mean((p - nd$quantiles)^2))
  err_m <- sqrt(mean((pm - nd$quantiles)^2))
  err_g <- sqrt(mean((pg - nd$quantiles)^2))

  ## Compute loss.
  err_b <- sqrt(mean((p - nd$quantiles)^2))
  err_m <- sqrt(mean((pm - nd$quantiles)^2))
  err_g <- sqrt(mean((pg - nd$quantiles)^2))

  ## Plot data and fitted median.
  if(is.null(ylim))
    ylim <- range(d$y, nd$y, p, pm, pg)

  plot(y ~ x, data = d, xlim = xlim, ylim = ylim,
    col = rgb(0.1, 0.1, 0.1, alpha = 0.2), pch = 16,
    main = "", xlab = "x", ylab = "y")

  i <- order(nd$x)

  matplot(nd$x[i], nd$quantiles[i,], type = "l",
    col = rgb(0.1, 0.1, 0.1, alpha = 0.3),
    lty = 1, lwd = 6, add = TRUE)

  matplot(nd$x[i], pm[i, ], type = "l",
    lty = 1, col = 2, lwd = 2, add = TRUE)
  matplot(nd$x[i], pg[i, ], type = "l",
    lty = 1, col = 3, lwd = 2, add = TRUE)
  matplot(nd$x[i], p[i, ], type = "l",
    lty = 1, col = 4, lwd = 2, add = TRUE)

  qe <- paste(c("TM", paste("GAMLSS", family()$family[1]), "QGAM"), "=", round(c(err_b, err_m, err_g), 2))

  legend(pos, qe, lwd = 2, col = c(4, 2, 3), bty = "n") ##title = "Quantile Error (RMSE)")
}


## (1)
if(FALSE) {
x11(width = 10, height = 6)

par(mfrow = c(2, 3), mar = rep(0, 4), oma = c(4, 4, 3, 3))


sim_NO(n = 1000, counts = TRUE, breaks = NULL, probs = 0.5,
  family = NBI, ylim = c(2, 20), axes = FALSE)
box()
axis(2)
mtext("Median", side = 3, line = 1, font = 2)

sim_NO(n = 1000, counts = TRUE, breaks = NULL, probs = c(0.1, 0.9),
  family = NBI, ylim = c(2, 20), axes = FALSE)
box()
mtext("10% and 90% Quantile", side = 3, line = 1, font = 2)

sim_NO(n = 1000, counts = TRUE, breaks = NULL, probs = c(0.01, 0.99),
  family = NBI, ylim = c(2, 20), axes = FALSE)
box()
mtext("1% and 99% Quantile", side = 3, line = 1, font = 2)


sim_NO(n = 1000, counts = TRUE, breaks = 200, probs = 0.5,
  family = NO, ylim = c(2, 20), axes = FALSE)
box()

sim_NO(n = 1000, counts = TRUE, breaks = 200, probs = c(0.1, 0.9),
  family = NO, ylim = c(2, 20), axes = FALSE) 
box()
axis(1)

sim_NO(n = 1000, counts = TRUE, breaks = 200, probs = c(0.01, 0.99),
  family = NO, ylim = c(2, 20), axes = FALSE)
box()
axis(4)

mtext("x", side = 1, line = 2.5, outer = TRUE)
mtext("y", side = 2 , line = 2.5, outer = TRUE)
}


## (2)
if(TRUE) {
x11(width = 10, height = 6)

par(mfrow = c(2, 3), mar = rep(0, 4), oma = c(4, 4, 3, 3))

ylim <- c(-2.7, 3.5)
seed <- 3

sim_NO(n = 1000, counts = FALSE, breaks = 30, probs = 0.5,
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()
axis(2)
mtext("Median", side = 3, line = 1, font = 2)

sim_NO(n = 1000, counts = FALSE, breaks = 30, probs = c(0.1, 0.9),
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()
mtext("10% and 90% Quantile", side = 3, line = 1, font = 2)

sim_NO(n = 1000, counts = FALSE, breaks = 30, probs = c(0.01, 0.99),
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()
mtext("1% and 99% Quantile", side = 3, line = 1, font = 2)


sim_NO(n = 1000, counts = FALSE, breaks = 200, probs = 0.5,
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()

sim_NO(n = 1000, counts = FALSE, breaks = 200, probs = c(0.1, 0.9),
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()
axis(1)

sim_NO(n = 1000, counts = FALSE, breaks = 200, probs = c(0.01, 0.99),
  family = NO, ylim = ylim, axes = FALSE, seed = seed)
box()
axis(4)

mtext("x", side = 1, line = 2.5, outer = TRUE)
mtext("y", side = 2 , line = 2.5, outer = TRUE)
}


if(FALSE) {
x11(width = 10, height = 4)

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

sim_ZINBI(1000, probs = 0.5)
sim_ZINBI(1000, probs = c(0.1, 0.9))
sim_ZINBI(1000, probs = c(0.01, 0.99))
}

