## Main paper: https://link.springer.com/article/10.1007/s10260-021-00558-6


#' Detect number of cores for OpenMP
#'
#' The calculation of CDFs and PDFs is implemented in C and allows
#' for parallelization using OpenMP. This function detects how may
#' cores are available in total (if OpenMP is available).
#'
#' @param verbose logical, if \code{TRUE} a message is shown.
#'
#' @return Number of available cores (integer). If OpenMP is not
#' available, \code{1L} is returned.
#'
#' @author Reto
tm_detect_cores <- function(verbose = TRUE) {
    stopifnot("'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose))
    ncores <- .Call("tm_detect_cores")
    if (verbose && ncores) {
        message("OMP available, number of cores detected: ", ncores)
    } else if (verbose) {
        message("OMP not available (not compiled with omp)")
    }
    return(ncores)
}

#' Get number of cores for OpenMP
#'
#' Some parts of the package use C routines which allow for parallelization
#' using OpenMP. This function is used to specify how many cores to be used.
#'
#' @param ncores \code{NULL} or a positive integer.
#' @param verbose logical, if \code{TRUE} a message is shown.
#'
#' @return Number of cores to be used in OpenMP parallelization (integer).
#'
#' @details If \code{ncores} is \code{NULL} the number of available
#' cores is auto-detected and set to 'total number of cores - 2'.
#' If integer, it is checked if this number of cores is available,
#' else set tot he 'total number of cores available'.
#'
#' @author Reto
tm_get_number_of_cores <- function(ncores = NULL, verbose = verbose) {
  ## Number of cores to be used for OpenMP. If
  ## - NULL: Guess cores (max cores - 2L)
  ## - Smaller or equal to 0: Set to 1L (single-core processing)
  ## - Else: Take user input; limited to maximum number of detected cores.
  ncores <- if (!is.null(ncores)) as.integer(ncores)[1L] else tm_detect_cores(verbose = FALSE) - 2L
  ncores <- if (ncores < 1L) 1L else pmin(ncores, tm_detect_cores(verbose = FALSE))
  if (verbose) message("Number of cores set to: ", ncores)
  return(ncores)
}


## Function to set up expanded data set.
tm_data <- function(data, response = NULL, verbose = TRUE) {
  ## Ensure data is a data frame.
  if (!is.data.frame(data))
    data <- as.data.frame(data)

  ## Determine response column.
  if (is.null(response))
    response <- names(data)[1L]

  if (!response %in% names(data))
    stop("The specified response column does not exist in the data!")

  ## Initialize progress bar.
  if (verbose)
    pb <- utils::txtProgressBar(min = 0, max = nrow(data), style = 3)

  ## Preallocate list.
  n    <- nrow(data)
  step <- if (n > 20) floor(n / 20) else 1

  response_values <- data[[response]]
  df_list <- vector("list", n)

  ## If any missing value in response: stop
  if (any(is.na(response_values)))
    stop("NA values in response data!")

  ## Setting up the new data.frame with (pseudo-)bins

  ## Length of vectors in list; names of list elements.
  nout <- sum(response_values) + length(response_values)
  names_out <- c("index", "Y", "theta", names(data))

  ## Building index vector; each observation 1:n gets its
  ## unique index (ID).
  result <- list()
  result$index <- rep(seq_len(n), response_values + 1L)

  ## Creating Y; always 1 except for the last entry per index.
  result$Y     <- rep(1L, nout)
  result$Y[cumsum(response_values + 1)] <- 0L

  ## Creating theta; a sequence from zero to the
  ## response_value for each index. The following
  ## Two lines create this sequence of sequences.
  reset_to_zero <- c(0, which(diff(result$index) > 0))
  result$theta <- seq_len(nout) - rep(reset_to_zero, response_values + 1) - 1

  ## Appending the remaining data from 'data'.
  for (n in names(data)) {
      ## If data[[n]] is a simple vector
      if (!is.matrix(data[[n]])) {
        result[[n]] <- rep(data[[n]], response_values + 1)
      ## Else create matrix
      } else {
        result[[n]] <- matrix(rep(data[[n]], rep(response_values + 1, ncol(data[[n]]))),
                              ncol = ncol(data[[n]]),
                              dimnames = list(NULL, colnames(data[[n]])))
      }
  }

  result <- as.data.frame(result)

  ## Attach the response column as an attribute.
  attr(result, "response") <- response

  return(result)
}

## Predict function.
tm_predict <- function(object, bins, newdata,
  type = c("pdf", "cdf", "quantile", "pmax"),
  response = NULL, y = NULL, prob = 0.5, maxcounts = 1e+03,
  verbose = FALSE, theta_scaler = NULL, theta_vars = NULL,
  factor = FALSE, ncores = NULL)
{

  ## Pre-processing inputs
  type <- tolower(type)
  type <- match.arg(type)

  ## TODO(R): I assume this check is not correct
  ##if (length(prob) > 1)
  ##  warning("Argument 'prob' has length > 1, only first element will be used.")

  ##if (type == "quantile") {
  ##  prob <- as.numeric(prob)[1L]
  ##  stopifnot("'prob' must be numeric in [0, 1]" = prob >= 0 & prob <= 1)
  ##}

  if (is.null(ncores))
    ncores <- tm_get_number_of_cores(ncores, verbose = verbose)

  ## Staying sane
  stopifnot(
    "'newdata' must be data.frame" = is.data.frame(newdata),
    "'response' must be character of length 1 (or NULL)" =
        (is.character(response) && length(response) == 1) || is.null(response),
    "'y' must be NULL or a vector" = is.null(y) || is.vector(y),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose),
    "'factor' must be logical TRUE or FALSE" = isTRUE(factor) || isFALSE(factor),
    "'ncores' must be NULL or numeric >= 1" = is.null(ncores) || (ncores >= 1)
  )
  ## TODO(R) Not all arguments are checked above

  ## If response is not specified explicitly, assume the first
  ## variable in the response; should be avoided if possible.
  if (is.null(response)) {
    response <- names(newdata)[1L]
  ## If response was specified but does not exist; set to maxcounts
  } else if (is.null(newdata[[response]])) {
    ## TODO(N): Reto is not sure when this is used and what
    ##          the reason for the + floor(0.1 * maxcounts) is?
    newdata[[response]] <- maxcounts + floor(0.1 * maxcounts)
  }

  ## Preparing data
  nd <- tm_data(newdata, response = response, verbose = verbose)
  if (factor)
    nd$theta <- as.factor(nd$theta)

  if (!is.null(theta_vars) && length(theta_vars) > 0L) {
    for (j in theta_vars) {
      i <- as.integer(gsub("theta", "", j))
      nd[[j]] <- as.integer(nd$theta == i)
    }
  }

  ## Scaling (standardizing) theta if requested
  if (!is.null(theta_scaler))
    nd$theta <- (nd$theta - theta_scaler$mean) / theta_scaler$sd

  ## Specify argument for generic prediction method called below
  what <- switch(class(object)[1L],
    "bam"  = "response",
    "gam"  = "response",
    "nnet" = "raw"
  )
  tp <- as.numeric(predict(object, newdata = nd, type = what))

  ## Extract unique indices
  ui <- unique(nd$index)
  prob <- rep(prob, length(ui))

  ## Ensure we hand over the correct thing to C
  if (type == "quantile") {
      stopifnot(is.numeric(prob), length(prob) == length(ui),
                all(!is.na(prob)), all(prob >= 0 & prob <= 1))
  } else if (is.null(prob)) {
      prob <- NA_real_ # dummy value for C (not used if type != 'quantile')
  }

  probs <- .Call("tm_predict",
                 uidx  = ui,                       # Unique distribution index (int)
                 idx   = nd$index,                 # Index vector (int)
                 tp    = tp,                       # Transition probabilities
                 lower = NA_real_,                 # Lower edge of the bin
                 upper = NA_real_,                 # Upper edge of the bin
                 y     = prob,                     # Where to evaluate the pdf
                 type  = type, ncores = ncores, elementwise = TRUE,
                 discrete = FALSE) # <- dummy value
  if (type == "quantile") print(summary(probs))

  return(probs)
}

tm_dist <- function(y, data = NULL, ...)
{
  if (is.null(y))
    stop("argument y is NULL!")

  is_f <- FALSE
  if (inherits(y, "formula")) {
    yn <- response_name(y)
    f <- y
    is_f <- TRUE
  } else {
    yn <- deparse(substitute(y), backtick = TRUE, width.cutoff = 500)
    f <- y ~ s(theta)
    environment(f) <- environment(y)
    data <- list()
    data[["y"]] <- y
  }

  ## Estimate model.
  b <- tm(f, data = data, ...)

  if (inherits(y, "formula"))
    y <- model.response(b$model.frame)

  ## Predict probabilities.
  nd <- data.frame("y" = 0:b$maxcounts)
  if ((yn != "y") & is_f)
    names(nd) <- yn
  pb <- predict(b, newdata = nd)

  nl <- NULL

  if (is.null(b$yc_tab)) {
    if (!is.null(data) & is_f)
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
  if (is.null(ylim)) 
    ylim <- range(c(0, tab, pb * 1.1), na.rm = TRUE)
  ylab <- list(...)$ylab
  if (is.null(ylab))
    ylab <- "Probability"
  xlab <- list(...)$xlab
  if (is.null(xlab)) {
    xlab <- if (is.null(b$yc_tab)) "#Counts" else yn
  }

  ## Plot.
  if (!is.null(nl)) {
    names(tab) <- nl[as.integer(names(tab)) + 1L]
  }
  x <- barplot(tab, xlab = xlab, ylab = ylab, ylim = ylim)
  lines(pb ~ x, col = 4, lwd = 2, type = "h")
  points(x, pb, col = 4, pch = 16)
  points(x, pb, col = rgb(0.1, 0.1, 0.1, alpha = 0.6))

  return(invisible(b))
}

## Function to create bins.
make_bins <- function(y, breaks = 30, scale = FALSE , ...) {
  if(scale) {
    my <- min(y)
    y <- sqrt(y - my + 0.01)
    dy <- diff(range(y))
    bins <- (seq(min(y) - 0.1*dy,
      max(y) + 0.1*dy, length = breaks))^2 - 0.01 + my
  } else {
    dy <- diff(range(y))
    bins <- seq(min(y) - 0.1*dy,
      max(y) + 0.1*dy, length = breaks)
#    bins <- c(min(y) - 0.5*dy,
#      quantile(y, probs = seq(0, 1, length = breaks)),
#      max(y) + 0.5*dy)
  }
  return(bins)
}

## Wrapper function to estimate CTMs.
tm <- function(formula, data, subset, na.action,
  engine = "bam", scale.x = FALSE, breaks = NULL,
  model = TRUE, ncores = NULL, verbose = FALSE, ...)
{
  if (!is.null(breaks)) breaks <- as.numeric(breaks)

  ## Staying sane
  stopifnot(
    "'ncores' must be NULL or numeric" = is.null(ncores) || is.numeric(ncores),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose)
  )
  ncores <- tm_get_number_of_cores(ncores, verbose = verbose)

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
  if (length(tv)) {
    for (j in tv)
      ff2 <- eval(parse(text = paste0("update(ff2, . ~ . -", j, ")")))
  }
  mf[["formula"]] <- ff2
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  ## Response name.
  rn <- response_name(formula)

  yscale <- NULL

  ## Discretize response?
  if (bin.y <- !is.null(breaks)) {
    if (length(breaks) < 2L) {
      bins <- make_bins(model.response(mf), breaks = breaks)
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

    if (!is.null(lower))
      ym[ym < lower] <- lower
    if (!is.null(upper))
      ym[ym > upper] <- upper

    mf[[rn]] <- yc
  }

  ## Scaling data.
  scaler <- NULL
  if (scale.x & (ncol(mf) > 1L)) {
    scaler <- list()
    for (j in names(mf)[-1L]) {
      if (!is.factor(mf[[j]])) {
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

  if (!is.null(scaler)) {
    scaler$theta <- list("mean" = mean(tmf$theta), "sd" = sd(tmf$theta))
    tmf$theta <- (tmf$theta - scaler$theta$mean) / scaler$theta$sd
  }

  if (length(tv)) {
    for (j in tv) {
      i <- as.integer(gsub("theta", "", j))
      tmf[[j]] <- as.integer(tmf$theta == i)
    }
  }

  ## Setup return value.
  rval <- list()

  ## New formula.
  if (engine %in% c("gam", "bam")) {
    if (isTRUE(list(...)$factor)) {
      tmf$theta <- as.factor(tmf$theta)
      rval$new_formula <- update(formula, as.factor(Y) ~ theta + .)
    } else {
      ##rval$new_formula <- eval(parse(text = paste0("update(formula, Y ~ s(theta,k=", k, ") + .)")))
      rval$new_formula <- update(formula, Y ~ .)
    }
  } else {
    rval$new_formula <- update(formula, as.factor(Y) ~ .)
  }

  ## Estimate model.
  warn <- getOption("warn")
  options("warn" = -1)
  if (engine == "bam") {
    rval$model <- bam(rval$new_formula, data = tmf, family = binomial, discrete = TRUE)
  } else if (engine == "gam") {
    rval$model <- gam(rval$new_formula, data = tmf, family = binomial, ...)
  } else if (engine == "nnet") {
    rval$model <- nnet::nnet(rval$new_formula, data = tmf, ...)
  } else if (engine == "glmnet") {
    rval$model <- tm_glmnet(rval$new_formula, data = tmf, ...)
  }
  options("warn" = warn)

  ## Additional info.
  rval$response <- rn
  rval$model.frame <- mf
  rval$scaler <- scaler
  rval$maxcounts <- max(mf[[rn]])
  rval$theta_vars <- tv
  rval$factor <- isTRUE(list(...)$factor)

  if (inherits(rval$model, "nnet")) {
    p <- predict(rval$model, type = "raw")
  } else {
    p <- predict(rval$model, type = "response")
  }

  ## Remove model frame.
  if (!model)
    rval$model$model <- NULL

  ## Compute probabilities.
  ui <- unique(tmf$index)
  probs <- cprobs <- numeric(length(ui))

  ## c_tm_predict_pdfcdf returns a list with PDF and CDF, calculating
  ## both simultanously in C to improve speed.
  tmp    <- .Call("tm_predict_pdfcdf", uidx = ui, idx = tmf$index, p = p, ncores = ncores)
  probs  <- tmp$pdf
  cprobs <- tmp$cdf
  rm(tmp)

  eps <- abs(.Machine$double.eps)
  probs[probs  < eps]      <- eps
  probs[probs  > 1 - eps]  <- 1 - eps
  cprobs[cprobs < eps]     <- eps
  cprobs[cprobs > 1 - eps] <- 1 - eps

  rval$probs <- data.frame("pdf" = probs, "cdf" = cprobs)

  ## If binning.
  if (bin.y) {
    rval$bins   <- bins
    rval$ym     <- ym
    rval$yc_tab <- table(yc)
    rval$breaks <- breaks
  }

  ## Assign class.
  class(rval) <- "tm"

  return(rval)
}

## Plotting method.
plot.tm <- function(x, which = "effects", spar = TRUE, k = 5, ...)
{
  ## What should be plotted?
  which.match <- c("effects", "hist-resid", "qq-resid", "wp-resid")
  if (!is.character(which)) {
    if (any(which > 4L))
      which <- which[which <= 4L]
    which <- which.match[which]
  } else which <- which.match[grep2(tolower(which), which.match, fixed = TRUE)]
  if (length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if (any("effects" %in% which) & inherits(x$model, "gam")) {
    plot(x$model, ...)
    return(invisible(NULL))
  } else {
    which <- which[which != "effects"]
    if (length(which) < 1L) {
      which <- c("hist-resid", "qq-resid", "wp-resid")
    }
  }

  resids <- NULL
  for (j in 1:k)
    resids <- cbind(resids, residuals(x, newdata = list(...)$newdata))
  resids <- apply(resids, 1, median)

  ## Number of plots.
  if (spar) {
    oma <- par(no.readonly = TRUE)
    par(mfrow = c(1, length(which)))
    on.exit(par(oma))
  }

  if ("hist-resid" %in% which) {
    plot_hist(resids, ...)
  }

  if ("qq-resid" %in% which) {
    plot_qq(resids, ...)
  }

  if ("wp-resid" %in% which) {
    plot_wp(resids, ...)
  }
}

## Summary method.
summary.tm <- function(object, ...)
{
  summary(object$model)
}

## formula method.
formula.tm <- function(x, ...)
{
  formula(x$model)
}


## Coef method.
coef.tm <- function(object, ...)
{
  coef(object$model)
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
  type = c("pdf", "cdf", "pmax", "quantile"), ncores = NULL, ...)
{
  type <- tolower(type)
  type <- match.arg(type)

  ## Get number of cores for OpenMP parallelization
  ncores <- tm_get_number_of_cores(ncores, FALSE)

  if (!is.null(prob))
    type <- "quantile"

  if (is.null(prob) && type == "quantile")
    prob <- 0.5

  if (is.null(newdata)) {
    if (type %in% c("pdf", "cdf") && is.null(y)) {
      ## Returning pdf/cdf stored on object
      return(object$probs[[type]])
    } else {
      ## Extracting model.frame used for model training as 'newdata'
      newdata <- model.frame(object)
    }
  }

  ## Overwriting response to evaluate pdf/cdf/pmax at specific y's,
  ## if type is either pmax or quantile: drop response altogether
  if (!is.null(y) && type %in% c("pdf", "cdf")) {
    newdata[[object$response]] <- y
  } else if (type %in% c("pmax", "quantile")) {
    newdata[[object$response]] <- NULL
  }

  ## Applying scaler if needed; standardize data (except theta)
  if (!is.null(object$scaler)) {
    for (j in names(object$scaler)) {
      if (j != "theta") {
        newdata[[j]] <- (newdata[[j]] - object$scaler[[j]]$mean) / object$scaler[[j]]$sd
      }
    }
  }

  ## Calling tm_predict to perform the actual prediction
  pred <- tm_predict(object$model,
                     bins         = object$bins,
                     newdata      = newdata,
                     ncores       = ncores,
                     response     = object$response,
                     type         = type,
                     maxcounts    = object$maxcounts,
                     prob         = prob,
                     theta_scaler = object$scaler$theta,
                     theta_vars   = object$theta_vars,
                     factor       = object$factor,
                     ...)

  if (type == "quantile") print(summary(pred))

  ## If binning is used (pseudo-counts), convert predicted
  ## bin to numeric (center of the bin)
  if (!is.null(object$bins)) {
    if (type %in% c("quantile", "pmax")) {
      pred <- object$ym[pred + 1L]
    }
  }

  ## Ensure PDF/CDF lie inside [0, 1]
  if (type %in% c("pdf", "cdf")) {
    eps <- abs(.Machine$double.eps)
    pred[pred < eps]     <- eps
    pred[pred > 1 - eps] <- 1 - eps
  }

  return(pred)
}

## logLik method.
logLik.tm <- function(object, newdata = NULL, ...)
{
  if (is.null(newdata)) {
    p <- object$probs$pdf
  } else {
    p <- predict(object, newdata = newdata, type = "pdf", ...)
  }
  ll <- sum(log(p))
  attr(ll, "nobs") <- nrow(object$model.frame)
  attr(ll, "df") <- sum(object$model$edf)
  class(ll) <- "logLik"
  return(ll)
}

## Residuals method.
residuals.tm <- function(object, newdata = NULL, y = NULL, ...)
{
  if (is.null(newdata)) {
    if (is.null(object$model.frame))
      stop("cannot compute residuals, no model.frame including the response!")
    newdata <- object$model.frame
  }

  if (is.null(y)) {
    y <- newdata[[object$response]]
  }

  if (is.null(y)) {
    stop("cannot compute residuals, response is missing!")
  }

  i <- y > 0
  p <- numeric(length(y))

  pL <- predict(object, newdata = newdata[i, , drop = FALSE], y = y[i] - 1L, type = "cdf", ...)
  pU <- predict(object, newdata = newdata[i, , drop = FALSE], y = y[i], type = "cdf", ...)

  p[i] <- runif(sum(i), pL, pU)

  if (any(!i)) {
    pU <- predict(object, newdata = newdata[!i, , drop = FALSE], y = y[!i], type = "cdf", ...)
    p[!i] <- runif (sum(!i), 0, pU)
  }

  eps <- abs(.Machine$double.eps)
  p[p < eps]     <- eps
  p[p > 1 - eps] <- 1 - eps

  return(qnorm(p))
}


## Rootogram method.
rootogram.tm <- function(object, newdata = NULL, plot = TRUE,
  width = 0.9, style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), expected = TRUE, confint = TRUE,
  ref = TRUE, K = NULL, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  if (is.null(newdata))
    newdata <- object$model.frame

  if (is.null(newdata[[object$response]]))
    stop("response missing in newdata!")

  y <- newdata[[object$response]]
  n <- length(y)

  if(!is.null(K))
    object$maxcounts <- K
  counts <- 0:object$maxcounts
  p <- NULL; obs <- NULL
  for (j in counts) {
    if (isTRUE(list(...)$verbose))
      cat(j, "/", sep = "")
    p <- cbind(p, predict(object, newdata = newdata, type = "pdf", y = j))
    obs <- c(obs, sum(y == j))
  }

  if (isTRUE(list(...)$verbose))
    cat("\n")

  e <- colMeans(p) * n  

  rg <- data.frame("observed" = obs, "expected" = e,
    "mid" = counts, "width" = width)

  scale <- match.arg(scale)

  if (scale == "sqrt") {
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
  attr(rg, "xlab") <- if (is.null(xlab)) "#Counts" else xlab
  attr(rg, "ylab") <- if (is.null(ylab)) "sqrt(Frequency)" else ylab
  attr(rg, "main") <- if (is.null(main)) "Rootogram" else main

  class(rg) <- c("rootogram", "data.frame")

  if (plot)
    plot(rg, ...)

  return(invisible(rg))
}

## Rootogram method.
rootogram.tmdist <- function(object, newdata = NULL, plot = TRUE,
  width = 0.9, style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), expected = TRUE, confint = TRUE,
  ref = TRUE, K = NULL, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  if (is.null(newdata))
    newdata <- object$model.frame

  if (is.null(newdata[[object$response]]))
    stop("response missing in newdata!")

  y <- newdata[[object$response]]
  n <- length(y)

  if(!is.null(K))
    object$maxcounts <- K
  counts <- 0:object$maxcounts
  p <- NULL; obs <- NULL
  for (j in counts) {
    if (isTRUE(list(...)$verbose))
      cat(j, "/", sep = "")
    p <- cbind(p, predict(object, newdata = newdata, type = "pdf", y = j))
    obs <- c(obs, sum(y == j))
  }

  if (isTRUE(list(...)$verbose))
    cat("\n")

  e <- colMeans(p) * n  

  rg <- data.frame("observed" = obs, "expected" = e,
    "mid" = counts, "width" = width)

  scale <- match.arg(scale)

  if (scale == "sqrt") {
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
  attr(rg, "xlab") <- if (is.null(xlab)) "#Counts" else xlab
  attr(rg, "ylab") <- if (is.null(ylab)) "sqrt(Frequency)" else ylab
  attr(rg, "main") <- if (is.null(main)) "Rootogram" else main

  class(rg) <- c("rootogram", "data.frame")

  if (plot)
    plot(rg, ...)

  return(invisible(rg))
}

