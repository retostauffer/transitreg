histtree <- function(formula, data, subset, na.action = na.pass,
   breaks = NULL, model = TRUE, batch_ids = NULL, nu = 1, score = log_loss,
   combine = c("mean", "geom"), ...)
{
  combine <- match.arg(combine)

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
      if(length(breaks) < 2L) {
        breaks <- unique(seq(min(Y), max(Y), length = as.integer(breaks)))
      }
      y_bin <- cut(Y, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      y_con <- TRUE
      nbins <- length(breaks) - 1L
      y_mids <- (breaks[-(nbins + 1L)] + breaks[-1L]) / 2
    }
  }

  ## Ensure helper objects exist and labels start at 1.
  if(!exists("breaks"))
    breaks <- NULL
  if(!exists("y_mids"))
    y_mids <- NULL
  if(min(y_bin, na.rm = TRUE) < 1L)
    y_bin <- y_bin - min(y_bin, na.rm = TRUE) + 1L
  if(!exists("nbins"))
    nbins <- max(y_bin, na.rm = TRUE)

  ## Process batches.
  if(is.null(batch_ids)) {
    batch_ids <- list(1:nrow(mf))
  } else {
    if(!is.list(batch_ids))
      batch_ids <- lapply(1:batch_ids, function(i) sample(1:nrow(mf), size = floor(0.63 * nrow(mf))))
  }
  nbatches <- length(batch_ids)

  ## Storage for ensemble and diagnostics.
  trees <- list()
  oos_scores <- rep(NA_real_, nbatches)         ## per-tree OOS score on its j
  ens_scores <- rep(NA_real_, nbatches)         ## ensemble 1, ..., i on j
  ens_scores_prev <- rep(NA_real_, nbatches)    ## ensemble 1, ..., i-1 on j
  ens_improve <- rep(NA_real_, nbatches)        ## delta = ens_curr - ens_prev (negative is better)
  k <- 1L

  ## Cycle over batches.
  for(i in 1:length(batch_ids)) {
    yi <- y_bin[batch_ids[[i]]]
    Xi <- mf[batch_ids[[i]], -attr(mt, "response"), drop = FALSE]
    treei <- build_tree(Xi, yi, nbins)

    if(nbatches > 1L) {
      j <- if(i < 2) nbatches else i - 1
      yj <- y_bin[batch_ids[[j]]]
      Xj <- mf[batch_ids[[j]], -attr(mt, "response"), drop = FALSE]

      ## per-tree OOS (optional)
      pj <- predict_tree(treei, Xj)
      oos_scores[i] <- score(pj, yj)

      ## ensemble WITHOUT the last tree (if i > 1)
      if(i > 1L) {
        P_list_prev <- lapply(1:length(trees), function(t) predict_tree(trees[[t]], Xj))
        P_prev <- aggregate_pmfs(P_list_prev, combine = combine)
        ens_scores_prev[i] <- score(P_prev, yj)
      }

      ## ensemble WITH all trees 1, ..., i
      P_list_curr <- c(if(i > 1L) P_list_prev else list(), list(pj))
      P_curr <- aggregate_pmfs(P_list_curr, combine = combine)
      ens_scores[i] <- score(P_curr, yj)

      ## improvement due to adding tree i on the SAME OOS batch j
      if(i > 1L) {
        ens_improve[i] <- ens_scores[i] - ens_scores_prev[i]
      } else {
        trees[[k]] <- treei
        k <- k + 1L
      }

      if((i > 1L) && (ens_improve[i] < 0)) {
        trees[[k]] <- treei
        k <- k + 1L
      }

      ## quick trace plot of ensemble OOS (lower is better)
      plot(ens_scores[!is.na(ens_scores)], type = "l",
        xlab = "Iteration", ylab = "Ensemble OOS score", lwd = 2)
      abline(h = min(ens_scores, na.rm = TRUE), lty = 3)
    } else {
      trees[[k]] <- treei
      k <- k + 1L
    }
  }

  rval <- list("y" = Y, "terms" = mt, "breaks" = breaks,
    "y_con" = y_con, "y_bin" = y_bin, "y_mids" = y_mids, "nbins" = nbins,
    "trees" = trees, "combine" = combine,
    "oos_scores" = oos_scores, "ens_scores" = ens_scores,
    "ens_scores_prev" = ens_scores_prev, "ens_improve" = ens_improve)

  if(model)
    rval$model <- mf

  class(rval) <- "histtree"

  return(rval)
}

## Hellinger Distance.
hellinger_distance <- function(p, q) {
  sqrt(sum((sqrt(p) - sqrt(q))^2))
}

## Weighted Jensen-Shannon Divergence.
weighted_js_distance <- function(pL, pR, piL = 0.5, base = 2) {
  piR <- 1 - piL
  m <- piL * pL + piR * pR
  ## Avoid log(0).
  KL <- function(a, b) sum(ifelse(a > 0 & b > 0, a * (log(a / b)), 0))
  jsd <- piL * KL(pL, m) + piR * KL(pR, m)
  sqrt(jsd / log(base))
}

## Mean negative log-likelihood (log-loss).
log_loss <- function(P, y, eps = 1e-12) {
  y <- pmin(pmax(y, 1L), ncol(P))
  p_obs <- P[cbind(seq_len(nrow(P)), y)]
  mean(-log(pmax(p_obs, eps)))
}

## Mean multi-class Brier score. Lower is better.
brier_score <- function(P, y) {
  n <- nrow(P); K <- ncol(P)
  y <- pmin(pmax(y, 1L), K)
  Y1 <- matrix(0, n, K)
  Y1[cbind(seq_len(n), y)] <- 1
  mean(rowSums((P - Y1)^2))
}

## Aggregate helper.
aggregate_pmfs <- function(P_list, weights = NULL, combine = c("mean", "geom")) {
  combine <- match.arg(combine)
  n <- nrow(P_list[[1L]]); K <- ncol(P_list[[1L]])
  nb <- length(P_list)
  if(is.null(weights)) {
    w <- rep(1 / nb, nb)
  } else {
    w <- weights / sum(weights)
  }
  if(combine == "mean") {
    P <- matrix(0, n, K)
    for(t in seq_len(nb))
      P <- P + w[t] * P_list[[t]]
  } else {
    eps <- 1e-15
    logP <- matrix(0, n, K)
    for(t in seq_len(nb))
      logP <- logP + w[t] * log(pmax(P_list[[t]], eps))
    P <- exp(logP)
  }
  rs <- rowSums(P); rs[!is.finite(rs) | rs <= 0] <- 1
  P <- P / rs
  return(P)
}

## Function that builds one tree.
build_tree <- function(X, y, nbins,
  distance = weighted_js_distance,
  depth = 0, max_depth = 100, min_group = 30, min_prop = 0.01)
{
  n <- nrow(X)
  req <- max(min_group, ceiling(min_prop * n))
  
  ## Stopping criteria: if node is too small or max depth reached
  if(n < 2 * req || depth >= max_depth) {
    pmf <- tabulate(y, nbins = nbins) / n
    return(list(is_leaf = TRUE, pmf = pmf, count = n))
  }

  best_split_val <- NULL
  best_score <- -Inf
  best_var_name <- NULL
  best_is_factor <- FALSE
  best_is_ordered <- FALSE
  
  ## Loop through all predictor variables.
  for(j in names(X)) {
    x <- X[[j]]
    
    if(is.factor(x) || is.character(x)) {
      ## Factor variable splitting logic.
      u <- if(is.factor(x)) levels(x) else unique(x[!is.na(x)])
      if(is.ordered(x)) {
        ## Ordered factor.
        if(length(u) >= 2L) for(split_val in u[-length(u)]) {
          idx_left <- x <= split_val
          idx_left[is.na(idx_left)] <- FALSE
          idx_right <- !idx_left
          
          n_left <- sum(idx_left)
          n_right <- sum(idx_right)
          if(n_left < req || n_right < req) next
          
          pmf_left  <- (tabulate(y[idx_left],  nbins = nbins) + 1e-9) / (n_left  + nbins * 1e-9)
          pmf_right <- (tabulate(y[idx_right], nbins = nbins) + 1e-9) / (n_right + nbins * 1e-9)
          pmf_left  <- pmf_left  / sum(pmf_left)
          pmf_right <- pmf_right / sum(pmf_right)

          ## Use true pi_L and a simple balance penalty.
          piL <- n_left / (n_left + n_right)
          d <- if(identical(distance, weighted_js_distance))
                  distance(pmf_left, pmf_right, piL = piL)
               else distance(pmf_left, pmf_right)
          bal <- 2 * piL * (1 - piL)                      ## in [0, 0.5], favors balanced splits
          current_score <- d * bal
          
          if(current_score > best_score) {
            best_score <- current_score
            best_split_val <- split_val
            best_var_name <- j
            best_is_factor <- TRUE
            best_is_ordered <- TRUE
          }
        }
      } else {
        ## Unordered factor (or character).
        splits <- power_set(u)
        for(split_set in splits) {
          idx_left <- if(is.factor(x)) as.character(x) %in% as.character(split_set) else x %in% split_set
          idx_left[is.na(idx_left)] <- FALSE
          idx_right <- !idx_left
          
          n_left <- sum(idx_left)
          n_right <- sum(idx_right)
          if(n_left < req || n_right < req) next
          
          pmf_left  <- (tabulate(y[idx_left],  nbins = nbins) + 1e-9) / (n_left  + nbins * 1e-9)
          pmf_right <- (tabulate(y[idx_right], nbins = nbins) + 1e-9) / (n_right + nbins * 1e-9)
          pmf_left  <- pmf_left  / sum(pmf_left)
          pmf_right <- pmf_right / sum(pmf_right)

          piL <- n_left / (n_left + n_right)
          d <- if(identical(distance, weighted_js_distance))
                  distance(pmf_left, pmf_right, piL = piL)
               else distance(pmf_left, pmf_right)
          bal <- 2 * piL * (1 - piL)
          current_score <- d * bal
          
          if(current_score > best_score) {
            best_score <- current_score
            best_split_val <- split_set
            best_var_name <- j
            best_is_factor <- TRUE
            best_is_ordered <- FALSE
          }
        }
      }
    } else {
      ## Numeric variable splitting logic.
      ## ux <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
      ux <- sort(unique(x)) ##seq(ux[1], ux[2], length = 10)
      ux <- ux[is.finite(ux)]
      
      for(split_val in ux) {
        idx_left <- x <= split_val
        idx_left[is.na(idx_left)] <- FALSE
        idx_right <- !idx_left

        n_left <- sum(idx_left)
        n_right <- sum(idx_right)
        if(n_left < req || n_right < req) next

        pmf_left  <- (tabulate(y[idx_left],  nbins = nbins) + 1e-9) / (n_left  + nbins * 1e-9)
        pmf_right <- (tabulate(y[idx_right], nbins = nbins) + 1e-9) / (n_right + nbins * 1e-9)
        pmf_left  <- pmf_left  / sum(pmf_left)
        pmf_right <- pmf_right / sum(pmf_right)

        piL <- n_left / (n_left + n_right)
        d <- if(identical(distance, weighted_js_distance))
                distance(pmf_left, pmf_right, piL = piL)
             else distance(pmf_left, pmf_right)
        bal <- 2 * piL * (1 - piL)
        current_score <- d * bal
        
        if(current_score > best_score) {
          best_score <- current_score
          best_split_val <- split_val
          best_var_name <- j
          best_is_factor <- FALSE
          best_is_ordered <- FALSE
        }
      }
    }
  }

  ## If no optimal split found (e.g., all splits resulted in too small groups).
  if(is.null(best_var_name)) {
    pmf <- tabulate(y, nbins = nbins) / n
    return(list("is_leaf" = TRUE, "pmf" = pmf, "count" = n))
  }
  
  ## Determine the best split indices based on the type of variable.
  if(best_is_factor) {
    if(best_is_ordered) {
      idx_left <- X[[best_var_name]] <= best_split_val
    } else {
      xv <- X[[best_var_name]]
      if(is.factor(xv)) {
        idx_left <- as.character(xv) %in% as.character(best_split_val)
      } else {
        idx_left <- xv %in% best_split_val
      }
    }
  } else {
    idx_left <- X[[best_var_name]] <= best_split_val
  }
  idx_left[is.na(idx_left)] <- FALSE
  idx_right <- !idx_left
  
  ## Recursively build left and right children.
  left_child <- build_tree(X[idx_left, , drop = FALSE], y[idx_left], nbins,
    distance, depth + 1, max_depth, min_group, min_prop)
  right_child <- build_tree(X[idx_right, , drop = FALSE], y[idx_right], nbins,
    distance, depth + 1, max_depth, min_group, min_prop)

  ## Return the current node's information.
  rval <- list(
    "is_leaf" = FALSE,
    "var" = best_var_name,
    "is_factor" = best_is_factor,
    "is_ordered" = best_is_ordered,
    "split_val" = best_split_val,
    "score" = best_score,
    "left" = left_child,
    "right" = right_child,
    "count" = n
  )

  return(rval)
}

## Helper function to generate power set for unordered factor splits.
power_set <- function(x) {
  n <- length(x)
  if(n > 15) { ## Limit to prevent combinatorial explosion.
    stop("Too many levels for a factor variable. Please reduce the number of levels or treat as ordered.")
  }
  elements <- vector("list", 2^n - 2)
  counter <- 1
  for(i in 1:(2^n - 2)) {
    subset_indices <- which(intToBits(i)[1:n] == 1)
    elements[[counter]] <- x[subset_indices]
    counter <- counter + 1
  }
  return(elements)
}

predict_tree <- function(tree, newdata) {
  stopifnot(!is.null(tree))
  n <- nrow(newdata)
  if(n == 0) return(matrix(numeric(0), 0, 0))

  ## Find number of bins from the first leaf.
  get_k <- function(node) if(isTRUE(node$is_leaf)) length(node$pmf) else get_k(node$left)
  k <- get_k(tree)

  out <- matrix(NA_real_, n, k)

  fill <- function(node, idx) {
    if(length(idx) == 0) return(invisible(NULL))

    if(isTRUE(node$is_leaf)) {
      out[idx, ] <<- matrix(node$pmf, nrow = length(idx), ncol = k, byrow = TRUE)
      return(invisible(NULL))
    }

    x <- newdata[[node$var]][idx]

    if(isTRUE(node$is_factor)) {
      if(isTRUE(node$is_ordered)) {
        go_left <- x <= node$split_val
      } else {
        go_left <- as.character(x) %in% as.character(node$split_val)
      }
    } else {
      go_left <- x <= node$split_val
    }

    ## Route NAs/unseen levels to the right by default.
    go_left[is.na(go_left)] <- FALSE

    idxL <- idx[go_left]
    idxR <- idx[!go_left]

    fill(node$left,  idxL)
    fill(node$right, idxR)
  }

  fill(tree, seq_len(n))

  ## Sanity: make sure rows sum to 1 (and backfill if something slipped through).
  rs <- rowSums(out)
  if(any(!is.finite(rs))) {
    ## Backstop: send to leftmost leaf.
    leftmost <- (function(node) { while(!node$is_leaf) node <- node$left; node })(tree)
    miss <- which(!is.finite(rs))
    out[miss, ] <- matrix(leftmost$pmf, nrow = length(miss), ncol = k, byrow = TRUE)
  } else if(any(abs(rs - 1) > 1e-10)) {
    out <- out / rs
  }

  return(out)
}

## Prediction.
## Quantiles from histogram PMF.
pmf_quantiles <- function(P, breaks, probs = c(0.1, 0.5, 0.9), eps = 1e-15) {
  n <- nrow(P); K <- ncol(P)
  stopifnot(length(breaks) == K + 1L)
  Cum <- t(apply(P, 1L, cumsum))
  out <- matrix(NA_real_, n, length(probs))
  for(i in seq_len(n)) {
    for(p in seq_along(probs)) {
      q <- probs[p]
      j <- which(Cum[i, ] >= q)[1]
      if(is.na(j)) j <- K
      c0 <- if(j == 1L) 0 else Cum[i, j - 1L]
      w  <- (q - c0) / max(P[i, j], eps)
      out[i, p] <- breaks[j] + w * (breaks[j + 1L] - breaks[j])
    }
  }
  colnames(out) <- paste0("q", probs)
  out
}

## S3 method.
predict.histtree <- function(object, newdata = NULL,
  type = c("pmf", "class", "mode", "mean", "quantile", "cdf", "all"),
  probs = c(0.1, 0.5, 0.9), append = FALSE, ...)
{
  type <- match.arg(type)

  ## Prepare newdata.
  if(is.null(newdata)) {
    if(!is.null(object$model)) {
      newdata <- object$model[, -attr(object$terms, "response"), drop = FALSE]
    } else {
      stop("newdata must be supplied because the fit was created with model = FALSE.")
    }
  }

  ## Edge cases.
  if(length(object$trees) < 1L) {
    k <- object$nbins
    P <- matrix(1 / k, nrow(newdata), k)
  } else {
    ## Predict each tree and aggregate.
    P_list <- lapply(object$trees, function(tr) predict_tree(tr, newdata))
    P <- aggregate_pmfs(P_list, combine = object$combine)
  }

  ## Shared bits.
  K <- ncol(P)
  class_idx <- max.col(P, ties.method = "first")

  ## Helpers for continuous targets.
  y_con <- isTRUE(object$y_con)
  y_mids <- if(!is.null(object$y_mids)) object$y_mids else {
    if(y_con && !is.null(object$breaks)) {
      nb <- length(object$breaks) - 1L
      (object$breaks[-(nb + 1L)] + object$breaks[-1L]) / 2
    } else NULL
  }

  ## Dispatch by type.
  if(type == "pmf") {
    return(P)
  }

  if(type == "class") {
    if(!append) return(class_idx)
    preds <- data.frame(pred_class = class_idx)
    return(cbind(newdata, preds, row.names = rownames(newdata)))
  }

  if(type == "mode") {
    mode_val <- if(y_con && !is.null(y_mids)) y_mids[class_idx] else class_idx
    if(!append) return(mode_val)
    preds <- data.frame(pred_class = class_idx, pred_mode = mode_val)
    return(cbind(newdata, preds, row.names = rownames(newdata)))
  }

  if(type == "mean") {
    if(!y_con || is.null(y_mids))
      stop("type = 'mean' requires continuous targets with available y_mids/breaks.")
    mu <- as.numeric(P %*% y_mids)
    if(!append) return(mu)
    preds <- data.frame(pred_mean = mu)
    return(cbind(newdata, preds, row.names = rownames(newdata)))
  }

  if(type == "quantile") {
    if(!y_con || is.null(object$breaks))
      stop("type = 'quantile' requires continuous targets with available breaks.")
    Q <- pmf_quantiles(P, object$breaks, probs = probs)
    if(!append) return(Q)
    preds <- as.data.frame(Q)
    names(preds) <- paste0("pred_", colnames(preds))
    return(cbind(newdata, preds, row.names = rownames(newdata)))
  }

  if(type == "cdf") {
    C <- t(apply(P, 1L, cumsum))
    colnames(C) <- paste0("F", seq_len(K))
    return(C)
  }

  ## type == "all"
  res <- list(
    "pmf" = P,
    "cdf" = t(apply(P, 1L, cumsum)),
    "class" = class_idx,
    "mode" = if(y_con && !is.null(y_mids)) y_mids[class_idx] else class_idx,
    "mean" = if(y_con && !is.null(y_mids)) as.numeric(P %*% y_mids) else NULL,
    "quantiles" = if(y_con && !is.null(object$breaks)) pmf_quantiles(P, object$breaks, probs = probs) else NULL,
    "breaks" = object$breaks,
    "y_mids" = y_mids
  )

  if(!append) return(res)

  ## Append summaries to newdata (exclude heavy pmf/cdf).
  df <- data.frame(
    pred_class = res$class,
    pred_mode  = res$mode,
    pred_mean  = if(!is.null(res$mean)) res$mean else NA_real_
  )
  if(!is.null(res$quantiles)) {
    qdf <- as.data.frame(res$quantiles)
    names(qdf) <- paste0("pred_", names(qdf))
    df <- cbind(df, qdf)
  }
  cbind(newdata, df, row.names = rownames(newdata))
}


## Testing.
set.seed(1328)
n <- 3000

## Simulate 2 covariates.
d <- data.frame(
  x1 = round(runif(n, -6, 6), 2),
  x2 = round(rnorm(n), 2)
)
d$y <- 10 + sin(d$x1) + c(0, -6)[(d$x1 > 0) * 1 + 1] + rnorm(n, sd = exp(-2 + 2 * cos(d$x1)))

b <- histtree(y ~ x1 + x2, data = d, batch_ids = 50, breaks = 100)

Q  <- predict(b, d, type = "quantile", probs = c(0.05, 0.5, 0.95))

i <- order(d$x1)

plot(y ~ x1, data = d, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(d$x1[i], Q[i, ], type = "l", lty = 1, col = 4, lwd = 2, add = TRUE)

