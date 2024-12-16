## Main paper: https://link.springer.com/article/10.1007/s10260-021-00558-6

## Function to set up expanded data set.
tm_data <- function(data, response = NULL, verbose = TRUE) {
  ## Ensure data is a data frame.
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  ## Determine response column.
  if(is.null(response)) {
    response <- names(data)[1L]
  }
  if(!response %in% names(data)) {
    stop("The specified response column does not exist in the data!")
  }

  ## Initialize progress bar.
  if(verbose) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(data), style = 3)
  }

  ## Preallocate list.
  n <- nrow(data)
  if(n > 20)
    step <- floor(n / 20)
  else
    step <- 1
  response_values <- data[[response]]
  df_list <- vector("list", n)

  ## Process each row.
  for(i in seq_len(n)) {
    k <- response_values[i]
    if(is.na(k))
      stop("NA values in response data!")
    Y <- c(rep(1, k), 0)
    theta <- seq_len(k) - 1
    row_data <- data[i, , drop = FALSE]

    ## Create expanded data frame for the current row.
    di <- data.frame(
      "index" = i,
      Y = Y,
      theta = c(0:k),
      row_data[rep(1, length(Y)), ]
    )

    ## Add to list.
    df_list[[i]] <- di

    ## Update progress bar.
    if(verbose && (i %% step == 0)) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if(verbose) {
    close(pb)
  }

  ## Combine all rows into a single data frame.
  result <- do.call("rbind", df_list)

  ## Attach the response column as an attribute.
  attr(result, "response") <- response

  return(result)
}

## Predict function.
tm_predict <- function(object, newdata,
  type = c("pdf", "cdf", "quantile", "pmax"), 
  response = NULL, y = NULL, prob = 0.5, maxcounts = 1e+03,
  verbose = FALSE, theta_scaler = NULL, theta_vars = NULL,
  factor = FALSE)
{
  if(is.null(response))
    response <- names(newdata)[1L]

  if(is.null(newdata[[response]])) {
    newdata[[response]] <- maxcounts + floor(0.1 * maxcounts)
  }

  nd <- tm_data(newdata, response = response, verbose = verbose)
  if(factor)
    nd$theta <- as.factor(nd$theta)

  if(!is.null(theta_vars)) {
    for(j in theta_vars) {
      i <- as.integer(gsub("theta", "", j))
      nd[[j]] <- as.integer(nd$theta == i)
    }
  }

  if(!is.null(theta_scaler)) {
    nd$theta <- (nd$theta - theta_scaler$mean) / theta_scaler$sd
  }

  what <- switch(class(object)[1L],
    "bam" = "response",
    "gam" = "response",
    "nnet" = "raw"
  )

  p <- as.numeric(predict(object, newdata = nd, type = what))

  ## What to predict?
  type <- tolower(type)
  type <- match.arg(type)

  if(type == "quantile") {
    prob <- rep(prob, length.out = nrow(newdata))
  }

  ui <- unique(nd$index)

  probs <- numeric(length(ui))

  for(j in ui) {
    pj <- p[nd$index == j]
    k <- length(pj)

    if(type == "pdf") {
      probs[j] <- (1 - pj[k]) * prod(pj[-k])
    }

    if(type == "cdf") {
      cj <- 1 - pj[1]
      if(length(pj) > 1) {
        for(jj in 2:length(pj))
          cj <- cj + (1 - pj[jj]) * prod(pj[1:(jj - 1)])
      }
      probs[j] <- cj
    }

    if(type == "quantile") {
      cj <- 1 - pj[1]
      if(cj >= prob[j]) {
        probs[j] <- 0L
      } else {
        if(length(pj) > 1) {
          for(jj in 2:length(pj)) {
            cj <- cj + (1 - pj[jj]) * prod(pj[1:(jj - 1)])
            if(cj >= prob[j]) {
              probs[j] <- jj - 1L
              break
            }
          }
        } else {
          probs[j] <- 0L
        }
      }
    }

    if(type == "pmax") {
      cj <- numeric(k)
      cj[1] <- 1 - pj[1]
      if(length(pj) > 1) {
        for(jj in 2:length(pj))
          cj[jj] <- (1 - pj[jj]) * prod(pj[1:(jj - 1)])
      }
      probs[j] <- which.max(cj) - 1L
    }
  }

  return(probs)
}

tm_dist <- function(y, data = NULL, ...)
{
  if(is.null(y))
    stop("argument y is NULL!")

  is_f <- FALSE
  if(inherits(y, "formula")) {
    yn <- response_name(y)
    f <- y
    is_f <- TRUE
  } else {
    yn <- deparse(substitute(y), backtick = TRUE, width.cutoff = 500)
    f <- y ~ 1
    environment(f) <- environment(y)
    data <- list()
    data[["y"]] <- y
  }

  ## Estimate model.
  b <- tm(f, data = data, ...)

  if(inherits(y, "formula"))
    y <- model.response(b$model.frame)

  ## Predict probabilities.
  nd <- data.frame("y" = 0:b$maxcounts)
  if((yn != "y") & is_f)
    names(nd) <- yn
  pb <- predict(b, newdata = nd)

  nl <- NULL

  if(is.null(b$yc_tab)) {
    if(!is.null(data) & is_f)
      y <- data[[yn]]
    tab <- prop.table(table(y))
  } else {
    tab <- prop.table(b$yc_tab)
    nl <- format(b$ym, digits = 2)
  }

  tab2 <- numeric(b$maxcounts + 1L)
  names(tab2) <- as.character(0:b$maxcounts)
  tab2[names(tab)] <- tab
  tab <- tab2

  ## Set labels.
  ylim <- list(...)$ylim
  if(is.null(ylim)) 
    ylim <- range(c(0, tab, pb * 1.1), na.rm = TRUE)
  ylab <- list(...)$ylab
  if(is.null(ylab))
    ylab <- "Probability"
  xlab <- list(...)$xlab
  if(is.null(xlab)) {
    xlab <- if(is.null(b$yc_tab)) "#Counts" else yn
  }

  ## Plot.
  if(!is.null(nl)) {
    names(tab) <- nl[as.integer(names(tab)) + 1L]
  }
  x <- barplot(tab, xlab = xlab, ylab = ylab, ylim = ylim)
  lines(pb ~ x, col = 4, lwd = 2, type = "h")
  points(x, pb, col = 4, pch = 16)
  points(x, pb, col = rgb(0.1, 0.1, 0.1, alpha = 0.6))

  return(invisible(b))
}

## Wrapper function to estimate CTMs.
tm <- function(formula, data, subset, na.action,
  engine = "bam", scale.x = FALSE, breaks = NULL,
  model = TRUE, verbose = FALSE, ...)
{
  cl <- match.call()

  ## Evaluate the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ff <- as.formula(paste("~", paste(all.vars(fake_formula(formula)), collapse = "+")))
  environment(ff) <- environment(formula)
  ff2 <- formula
  ff2[3L] <- ff[2L]
  ff2 <- update(ff2, ~ . - theta)
  tv <- all.vars(ff2)
  tv <- grep("theta", tv, value = TRUE)
  tv <- tv[tv != "theta"]
  if(length(tv)) {
    for(j in tv)
      ff2 <- eval(parse(text = paste0("update(ff2, . ~ . -", j, ")")))
  }
  mf[["formula"]] <- ff2
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  ## Response name.
  rn <- response_name(formula)

  ## Discretize response?
  if(bin.y <- !is.null(breaks)) {
    if(length(breaks) < 2L) {
      dy <- diff(range(model.response(mf)))
      bins <- seq(min(model.response(mf)) - 0.1*dy,
        max(model.response(mf)) + 0.1*dy, length = breaks)
    } else {
      bins <- breaks
    }
    #bins[1] <- -Inf
    #bins[length(bins)] <- Inf

    ## Discretize numeric response into counts.
    yc <- cut(model.response(mf), breaks = bins, labels = FALSE,
      include.lowest = TRUE) - 1
    ym <- (bins[-1] + bins[-length(bins)]) / 2

    lower <- list(...)$lower
    upper <- list(...)$upper

    if(!is.null(lower))
      ym[ym < lower] <- lower
    if(!is.null(upper))
      ym[ym > upper] <- upper

    mf[[rn]] <- yc
  }

  ## Scaling data.
  scaler <- NULL
  if(scale.x & (ncol(mf) > 1L)) {
    scaler <- list()
    for(j in names(mf)[-1L]) {
      if(!is.factor(mf[[j]])) {
        scaler[[j]] <- list("mean" = mean(mf[[j]]), "sd" = sd(mf[[j]]))
        mf[[j]] <- (mf[[j]] - scaler[[j]]$mean) / scaler[[j]]$sd
      }
    }
  }

  ## Max. counts.
  ymax <- max(mf[[rn]], na.rm = TRUE)
  k <- min(c(ymax - 1L, 20L))

  ## Transform data.
  tmf <- tm_data(mf, response = rn, verbose = verbose)

  if(!is.null(scaler)) {
    scaler$theta <- list("mean" = mean(tmf$theta), "sd" = sd(tmf$theta))
    tmf$theta <- (tmf$theta - scaler$theta$mean) / scaler$theta$sd
  }

  if(length(tv)) {
    for(j in tv) {
      i <- as.integer(gsub("theta", "", j))
      tmf[[j]] <- as.integer(tmf$theta == i)
    }
  }

  ## Setup return value.
  rval <- list()

  ## New formula.
  if(engine %in% c("gam", "bam")) {
    if(isTRUE(list(...)$factor)) {
      tmf$theta <- as.factor(tmf$theta)
      rval$new_formula <- update(formula, as.factor(Y) ~ theta + .)
    } else {
      ##rval$new_formula <- eval(parse(text = paste0("update(formula, Y ~ s(theta,k=", k, ") + .)")))
      rval$new_formula <- update(formula, Y ~ .)
    }
  } else {
    rval$new_formula <- update(formula, as.factor(Y) ~ . + theta)
  }

  ## Estimate model.
  warn <- getOption("warn")
  options("warn" = -1)
  if(engine == "bam") {
    rval$model <- bam(rval$new_formula, data = tmf, family = binomial, discrete = TRUE)
  }
  if(engine == "gam") {
    rval$model <- gam(rval$new_formula, data = tmf, family = binomial, ...)
  }
  if(engine == "nnet") {
    rval$model <- nnet::nnet(rval$new_formula, data = tmf, ...)
  }
  options("warn" = warn)

  ## Additional info.
  rval$response <- rn
  rval$model.frame <- mf
  rval$scaler <- scaler
  rval$maxcounts <- max(mf[[rn]])
  rval$theta_vars <- tv
  rval$factor <- isTRUE(list(...)$factor)

  if(inherits(rval$model, "nnet")) {
    p <- predict(rval$model, type = "raw")
  } else {
    p <- predict(rval$model, type = "response")
  }

  ## Remove model frame.
  if(!model)
    rval$model$model <- NULL

  ## Compute probabilities.
  ui <- unique(tmf$index)
  probs <- cprobs <- numeric(length(ui))

  for(j in ui) {
    pj <- p[tmf$index == j]
    k <- length(pj)
    probs[j] <- (1 - pj[k]) * prod(pj[-k])

    cj <- numeric(k)
    cj[1] <- 1 - pj[1]
    if(length(pj) > 1) {
      for(jj in 2:length(pj))
        cj[jj] <- (1 - pj[jj]) * prod(pj[1:(jj - 1)])
    }
    cprobs[j] <- sum(cj)
  }

  probs[probs < 1e-15] <- 1e-15
  probs[probs > 0.999999] <- 0.999999
  cprobs[cprobs < 1e-15] <- 1e-15
  cprobs[cprobs > 0.999999] <- 0.999999

  rval$probs <- data.frame("pdf" = probs, "cdf" = cprobs)

  ## If binning.
  if(bin.y) {
    rval$bins <- bins
    rval$ym <- ym
    rval$yc_tab <- table(yc)
    rval$breaks <- breaks
  }

  ## Assign class.
  class(rval) <- "tm"

  return(rval)
}

## Plotting method.
plot.tm <- function(x, which = "effects", spar = TRUE, ...)
{
  ## What should be plotted?
  which.match <- c("effects", "hist-resid", "qq-resid", "wp-resid")
  if(!is.character(which)) {
    if(any(which > 4L))
      which <- which[which <= 4L]
    which <- which.match[which]
  } else which <- which.match[grep2(tolower(which), which.match, fixed = TRUE)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if(any("effects" %in% which) & inherits(x$model, "gam")) {
    plot(x$model, ...)
    return(invisible(NULL))
  } else {
    which <- which[which != "effects"]
    if(length(which) < 1L) {
      which <- c("hist-resid", "qq-resid", "wp-resid")
    }
  }

  resids <- residuals(x, newdata = list(...)$newdata)

  ## Number of plots.
  if(spar) {
    oma <- par(no.readonly = TRUE)
    par(mfrow = c(1, length(which)))
    on.exit(par(oma))
  }

  if("hist-resid" %in% which) {
    plot_hist(resids, ...)
  }

  if("qq-resid" %in% which) {
    plot_qq(resids, ...)
  }

  if("wp-resid" %in% which) {
    plot_wp(resids, ...)
  }
}

## Summary method.
summary.tm <- function(object, ...)
{
  summary(object$model)
}

## Model frame extractor.
model.frame.tm <- function(formula, ...)
{
  return(formula$model.frame)
}

## Printing method.
print.tm <- function(x, ...)
{
  cat("Count Transition Model\n---")
  print(x$model)
}

## Predict method.
predict.tm <- function(object, newdata = NULL,
  y = NULL, prob = NULL,
  type = c("pdf", "cdf", "pmax", "quantile"), ...)
{
  type <- tolower(type)
  type <- match.arg(type)

  if(!is.null(prob))
    type <- "quantile"

  if(is.null(prob) & (type == "quantile"))
    prob <- 0.5

  if(is.null(newdata)) {
    if((type %in% c("pdf", "cdf")) & is.null(y)) {
      return(object$probs[[type]])
    } else {
      newdata <- model.frame(object)
    }
  }

  if(!is.null(y))
    newdata[[object$response]] <- y

  if(type %in% c("pmax", "quantile"))
    newdata[[object$response]] <- NULL

  if(!is.null(object$scaler)) {
    for(j in names(object$scaler)) {
      if(j != "theta") {
        newdata[[j]] <- (newdata[[j]] - object$scaler[[j]]$mean) / object$scaler[[j]]$sd
      }
    }
  }

  pred <- tm_predict(object$model, newdata = newdata, 
    response = object$response, type = type,
    maxcounts = object$maxcounts, prob = prob,
    theta_scaler = object$scaler$theta,
    theta_vars = object$theta_vars,
    factor = object$factor, ...)

  if(!is.null(object$bins)) {
    if(type %in% c("quantile", "pmax")) {
      pred <- object$ym[pred + 1L]
    }
  }

  if(type %in% c("pdf", "cdf")) {
    pred[pred < 1e-15] <- 1e-15
    pred[pred > 0.999999] <- 0.999999
  }

  return(pred)
}

## logLik method.
logLik.tm <- function(object, newdata = NULL, ...)
{
  if(is.null(newdata)) {
    p <- object$probs$pdf
  } else {
    p <- predict(object, newdata = newdata, type = "pdf")
  }
  ll <- sum(log(p))
  attr(ll, "nobs") <- nrow(object$model.frame)
  attr(ll, "df") <- sum(object$model$edf)
  class(ll) <- "logLik"
  return(ll)
}

## Residuals method.
residuals.tm <- function(object, ...)
{
  p <- predict(object, ..., type = "cdf")
  p[p < 1e-15] <- 1e-15
  p[p > 0.999999] <- 0.999999
  return(qnorm(p))
}

## Rootogram method.
rootogram.tm <- function(object, newdata = NULL, plot = TRUE,
  width = 0.9, style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), expected = TRUE, confint = TRUE,
  ref = TRUE, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  if(is.null(newdata))
    newdata <- object$model.frame

  if(is.null(newdata[[object$response]]))
    stop("response missing in newdata!")

  y <- newdata[[object$response]]
  n <- length(y)

  counts <- 0:object$maxcounts
  p <- NULL; obs <- NULL
  for(j in counts) {
    if(isTRUE(list(...)$verbose))
      cat(j, "/", sep = "")
    p <- cbind(p, predict(object, newdata = newdata, type = "pdf", y = j))
    obs <- c(obs, sum(y == j))
  }

  if(isTRUE(list(...)$verbose))
    cat("\n")

  e <- colMeans(p) * n  

  rg <- data.frame("observed" = obs, "expected" = e,
    "mid" = counts, "width" = width)

  scale <- match.arg(scale)

  if(scale == "sqrt") {
    rg$observed <- sqrt(rg$observed)
    rg$expected <- sqrt(rg$expected)
  }

  p <- t(p)
  rownames(p) <- paste0("p_", counts + 0.5)
  colnames(p) <- as.character(1:n)
  rg$distribution <- p

  attr(rg, "style") <- match.arg(style)
  attr(rg, "scale") <- scale
  attr(rg, "expected") <- expected
  attr(rg, "confint") <- confint
  attr(rg, "ref") <- ref
  attr(rg, "xlab") <- if(is.null(xlab)) "#Counts" else xlab
  attr(rg, "ylab") <- if(is.null(ylab)) "sqrt(Frequency)" else ylab
  attr(rg, "main") <- if(is.null(main)) "Rootogram" else main

  class(rg) <- c("rootogram", "data.frame")

  if(plot)
    plot(rg, ...)

  return(invisible(rg))
}

