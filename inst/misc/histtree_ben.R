rm(list = ls())
graphics.off()

set.seed(1328)
n <- 1000
x <- sort(runif(n, -3, 3))
y <- sin(x) + rnorm(n, sd = 0.3)

# plot(x, y, ylim = c(-2.5, 2.5))

breaks <- seq(min(y), max(y), length = 50)

y_bin <- cut(y, breaks = breaks, labels = FALSE, include.lowest = TRUE)
nbins <- length(breaks) - 1L
y_mids <- (breaks[-(nbins + 1L)] + breaks[-1L]) / 2

nbatches <- 100
batch_ids <- lapply(1:nbatches, function(i) sample(1:n, size = floor(0.63 * n)))
# batch_ids <- lapply(1:nbatches, function(i) sample(1:n, size = floor(0.2 * n)))
# batch_ids <- replicate(100, 1:n, simplify = FALSE)
# nbatches <- length(batch_ids)

prior_left <- rep(1/nbins, nbins)
prior_right <- rep(1/nbins, nbins)
# alpha <- 1
alpha <- .1


weighted_js_distance <- function(pL, pR, piL = 0.5, base = 2) {
  piR <- 1 - piL
  m <- piL * pL + piR * pR
  ## Avoid log(0).
  KL <- function(a, b) sum(ifelse(a > 0 & b > 0, a * (log(a / b)), 0))
  jsd <- piL * KL(pL, m) + piR * KL(pR, m)
  sqrt(jsd / log(base))
}

best_split_val <- NULL
best_score <- -Inf

fits <- list()
err <- rep(NA_real_, nbatches)

log_loss <- function(P, y, eps = 1e-12) {
  y <- pmin(pmax(y, 1L), ncol(P))
  p_obs <- P[cbind(seq_len(nrow(P)), y)]
  mean(-log(pmax(p_obs, eps)))
}

for(i in 1:nbatches) {
  
  ##############################################################################
  # i <- 1
  ##############################################################################

  y_bini <- y_bin[batch_ids[[i]]]
  xi <- x[batch_ids[[i]]]
  qxi <- unname(quantile(xi, probs = seq(0.01, 0.99, 0.01), na.rm = TRUE))
  
  for(split_val in qxi) {
    
    ############################################################################
    # split_val <- qxi[50]
    ############################################################################
    
    idx_left <- xi <= split_val
    idx_left[is.na(idx_left)] <- FALSE
    idx_right <- !idx_left
    
    n_left <- sum(idx_left)
    n_right <- sum(idx_right)
    
    pmf_left  <- (tabulate(y_bini[idx_left],  nbins = nbins) + 1e-9) / (n_left  + nbins * 1e-9)
    pmf_right <- (tabulate(y_bini[idx_right], nbins = nbins) + 1e-9) / (n_right + nbins * 1e-9)
    pmf_left  <- pmf_left  / sum(pmf_left)
    pmf_right <- pmf_right / sum(pmf_right)
    
    pmf_left <- alpha * pmf_left + (1 - alpha) * prior_left
    pmf_right <- alpha * pmf_right + (1 - alpha) * prior_right
    
    piL <- n_left / (n_left + n_right)
    wdis <- weighted_js_distance(pmf_left, pmf_right, piL = piL)
    bal <- 2 * piL * (1 - piL)
    current_score <- wdis * bal
    
    if(current_score > best_score) {
      best_score <- current_score
      best_split_val <- split_val
      best_pmf_left <- pmf_left
      best_pmf_right <- pmf_right
    }
  }
  
  prior_left <- best_pmf_left
  prior_right <- best_pmf_right
  
  m_left <- mean(y[x <= best_split_val])
  m_right <- mean(y[x > best_split_val])
  fit <- rep(m_left, length(y))
  fit[x > best_split_val] <- m_right
  fits[[i]] <- fit
  # err[i] <- sum((y - fit)^2)
  
}

par(mfrow = c(1, 2))
plot(x, y, ylim = c(-2.5, 2.5))
sapply(breaks, function(b) abline(h = b, col = 'lightgrey'))
lapply(fits, function(f) lines(f ~ x, lwd = 0.5, col = 2))
lines(fit ~ x, lwd = 4, col = 4)
matplot(pmf_left, 1:nbins, type = 's', col = 2)
matplot(pmf_right, 1:nbins, type = 's', col = 3, add = TRUE)
# plot(1:nbatches, err, type = 'l')


