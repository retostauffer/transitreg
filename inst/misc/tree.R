## ------------------------------------------------------------------------------
## Title: Transition Probability Tree (TPT) Script
## Author: Nikolaus Umlauf
## Date: 2025-07-21
##
## Statement:
## This script contains the original implementation of transition probability trees
## as developed by Nikolaus Umlauf (2025) for conditional distribution estimation.
## Please cite appropriately if using or building on this work.
## ------------------------------------------------------------------------------

## Main ref: https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6729

## Hazard function from binned y.
hazard_estimates <- function(y_bin, m) {
  counts <- tabulate(y_bin, nbins = m)    ## Counts in each bin 1 to m.
  cum_counts <- rev(cumsum(rev(counts)))  ## At risk at r: sum(y_bin >= r).
  trans_counts <- c(cum_counts[-1], 0)    ## Transitions: sum(y_bin > r).
  h <- ifelse(cum_counts[-m] > 0, trans_counts[-m] / cum_counts[-m], 0)
  return(h)
}

## Transition probabilities.
transition_probs <- function(h) {
  m <- length(h) + 1
  ## Cumulative product of hazards.
  cumprod_h <- cumprod(h)
  ## Shifted cumulative product (prepend 1).
  cumprod_h_prev <- c(1, head(cumprod_h, -1))
  ## Probabilities.
  p <- cumprod_h_prev * (1 - h)
  ## Last bin: survival until the end.
  p[m] <- cumprod_h[length(h)]
  return(p)
}

## Find split using Hellinger distance between transition probabilities
find_split <- function(x, y_bin, m) {
  ## Candidate split points: equally spaced quantiles (excluding extremes).
  ux <- quantile(x, probs = c(0.001, 0.999))
  ux <- seq(ux[1], ux[2], length = 200)
  n <- length(x)
  out <- numeric(length(ux))

  for(i in seq_along(ux)) {
    split_val <- ux[i]
    j <- x <= split_val

    n_l <- sum(j)
    n_r <- n - n_l
    if(n_l < 30 || n_r < 30) next ## FIXME!

    ## Hazard estimation.
    h_left  <- hazard_estimates(y_bin[j], m)
    h_right <- hazard_estimates(y_bin[!j], m)

    ## Transition probabilities.
    p_left  <- transition_probs(h_left)
    p_right <- transition_probs(h_right)

    ## Hellinger distance.
    out[i] <- sqrt(sum((sqrt(p_left) - sqrt(p_right))^2))
  }

  data.frame(x = ux, err = out) |> na.omit()
}

set.seed(123)

## Simulated data: X uniformly in [-5, 5], split at X = 2.5.
n <- 1000
X <- runif(n, -5, 5)
y <- 10 + ifelse(X > 2, 1, -1) + rnorm(n, sd = 0.1)

## Discretize y into m bins.
m <- 100
breaks <- seq(min(y), max(y), length.out = m + 1)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
y_bin <- cut(y, breaks = breaks, labels = FALSE, include.lowest = TRUE)

## Find best split.
split_result <- find_split(X, y_bin, m)
split_est <- split_result$x[which.max(split_result$err)]
cat("Estimated split point:", split_est, "\n")

## Predictive transition probabilities.
j <- X <= split_est
h_left  <- hazard_estimates(y_bin[j], m)
h_right <- hazard_estimates(y_bin[!j], m)
p_left  <- transition_probs(h_left)
p_right <- transition_probs(h_right)

## Predict mean.
predict_means <- function(x, split, p_left, p_right, mids) {
  probs <- ifelse(x <= split, list(p_left), list(p_right))
  mu <- sapply(probs, function(p) sum(p * mids))
  return(mu)
}
mu_hat <- predict_means(X, split_est, p_left, p_right, midpoints)

## Predict quantiles.
predict_quantiles <- function(x, split, p_left, p_right, mids, probs = c(0.5)) {
  pred <- matrix(NA_real_, nrow = length(x), ncol = length(probs))
  colnames(pred) <- paste0("q", probs)

  cdf_left  <- cumsum(p_left)
  cdf_right <- cumsum(p_right)

  for(i in seq_along(x)) {
    cdf <- if(x[i] <= split) cdf_left else cdf_right
    for(j in seq_along(probs)) {
      idx <- which(cdf >= probs[j])[1]
      pred[i, j] <- mids[idx]
    }
  }

  return(pred)
}

q_hat <- predict_quantiles(X, split_est, p_left, p_right, midpoints, probs = c(0.15, 0.95))

print(round(mean(y <= q_hat[, 1]), 2))
print(round(mean(y >= q_hat[, 2]), 2))

## Plot: true vs. predicted mean + quantiles.
plot(y ~ X, col = adjustcolor(1, 0.3), pch = 16)
points(X, mu_hat, col = 2, pch = 16)
points(X, q_hat[, 1], col = 4, pch = ".")
points(X, q_hat[, 2], col = 4, pch = ".")
abline(v = split_est, col = "red", lty = 2)

## Full tree.
build_tree <- function(X, y_bin, m, depth = 1, max_depth = 100, minsize = 5) {
  n <- nrow(X)
  if(depth > max_depth || n < minsize) return(NULL)

  best_split <- NULL
  best_err <- -Inf
  best_var <- NULL

  for(var in names(X)) {
    xvar <- X[[var]]
    split_result <- find_split(xvar, y_bin, m)
    if(nrow(split_result) == 0) next
    i_max <- which.max(split_result$err)
    if(split_result$err[i_max] > best_err) {
      best_err <- split_result$err[i_max]
      best_split <- split_result$x[i_max]
      best_var <- var
    }
  }

  if(is.null(best_split)) return(NULL)

  ## Partition.
  left_idx <- which(X[[best_var]] <= best_split)
  right_idx <- which(X[[best_var]] > best_split)

  ## Stop if node is too small.
  if(length(left_idx) < minsize || length(right_idx) < minsize) return(NULL)

  ## Fit hazards and transition probabilities in each node.
  h_left  <- hazard_estimates(y_bin[left_idx], m)
  h_right <- hazard_estimates(y_bin[right_idx], m)
  p_left  <- transition_probs(h_left)
  p_right <- transition_probs(h_right)

  ## Recursively build child nodes.
  left  <- build_tree(X[left_idx, , drop = FALSE], y_bin[left_idx], m, depth + 1, max_depth, minsize)
  right <- build_tree(X[right_idx, , drop = FALSE], y_bin[right_idx], m, depth + 1, max_depth, minsize)

  return(list(
    var = best_var,
    split = best_split,
    p_left = p_left,
    p_right = p_right,
    left = left,
    right = right
  ))
}

predict_node <- function(node, xrow, type, mids, probs) {
  if(is.null(node$left) && is.null(node$right)) {
    ## Terminal node.
    if(type == "mean") return(sum(node$p_left * mids))
    if(type == "quantile") {
      cdf <- cumsum(node$p_left)
      return(sapply(probs, function(q) mids[which(cdf >= q)[1]]))
    }
  }

  ## Determine direction and recurse.
  direction <- if(xrow[[node$var]] <= node$split) "left" else "right"
  next_node <- node[[direction]]

  ## If child node is missing, fall back to current node's prediction.
  if(is.null(next_node)) {
    p <- if(direction == "left") node$p_left else node$p_right
    if(type == "mean") return(sum(p * mids))
    if(type == "quantile") {
      cdf <- cumsum(p)
      return(sapply(probs, function(q) mids[which(cdf >= q)[1]]))
    }
  }

  ## Recurse.
  predict_node(next_node, xrow, type, mids, probs)
}

predict_tree <- function(tree, Xnew, mids, type = c("mean", "quantile"), probs = 0.5) {
  type <- match.arg(type)
  
  if(type == "mean") {
    return(apply(Xnew, 1, function(x) predict_node(tree, as.list(x), type, mids, probs)))
  } else {
    qmat <- t(apply(Xnew, 1, function(x) predict_node(tree, as.list(x), type, mids, probs)))
    colnames(qmat) <- paste0("q", probs)
    return(qmat)
  }
}

n <- 1000

## Simulate 2 covariates.
d <- data.frame(
  x1 = runif(n, -6, 6),
  x2 = rnorm(n)
)
y <- 10 + sin(d$x1) + rnorm(n, sd = exp(-2 + 1.5 * cos(d$x1)))

m <- 100

## Binning.
breaks <- seq(min(y), max(y), length.out = m + 1)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
y_bin <- cut(y, breaks = breaks, labels = FALSE, include.lowest = TRUE)

## Build tree.
tree <- build_tree(d, y_bin, m, max_depth = 100)

## Predict mean.
mu_hat <- predict_tree(tree, d, mids = midpoints, type = "mean")

## Predict quantiles.
q_hat <- predict_tree(tree, d, mids = midpoints, type = "quantile", probs = c(0.05, 0.95))

## Plot predictions.
plot(y ~ d$x1, col = adjustcolor(1, 0.3), pch = 16)

i <- order(d$x1)

matplot(d$x1[i], cbind(mu_hat, q_hat)[i, ],
  type = "l", lty = 1, col = c(2, 4, 4), lwd = 3, add = TRUE)


transittree <- function(formula, data, subset, na.action = na.pass,
   breaks = NULL, ...)
{
  call <- match.call()
  if(missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) 
      names(Y) <- nm
  }

  y_con <- FALSE

  if(is.integer(Y)) {
    y_bin <- Y
  } else {
    if(is.factor(Y)) {
      y_bin <- as.integer(Y)
    } else {
      if(is.null(breaks))
        breaks <- seq(min(Y), max(Y), length = 50)
      y_bin <- cut(Y, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      y_con <- TRUE
      nbins <- length(breaks) - 1L
      y_mids <- (breaks[-(nbins + 1L)] + breaks[-1L]) / 2
    }
  }

  rval <- list("y" = Y, "model" = mf, "terms" = mt, "breaks" = breaks,
    "y_bin" = y_bin, "y_mids" = y_mids, "nbins" = nbins)

  return(rval)
}

b <- transittree(y ~ x1 + x2, data = d)

