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
    ncores <- .Call(C_tm_detect_cores)
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
tm_data <- function(data, response = NULL, useC = FALSE, verbose = TRUE) {
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
  timer("[tm_data] extracted response, allocated list", verbose)

  ## If any missing value in response: stop
  if (any(is.na(response_values)))
    stop("NA values in response data!")

  ## TODO(R) Original implementation is used when useC == FALSE,
  ##         though the new implementation has nothing to do with C.
  if (!useC) {

    ## Process each row.
    for (i in seq_len(n)) {
      k <- response_values[i] # Count or pseudocount
      Y <- c(rep(1, k), 0)
      row_data <- data[i, , drop = FALSE]

      ## Create expanded data frame for the current row.
      di <- data.frame(
        "index" = i,
        Y = Y,
        theta = c(0:k),
        row_data[rep(1, length(Y)), , drop = FALSE]
      )

      ## Add to list.
      df_list[[i]] <- di

      ## Update progress bar.
      if (verbose && (i %% step == 0)) {
        utils::setTxtProgressBar(pb, i)
      }
    }
    timer("[tm_data] loop over n - R", verbose)

    if (verbose)
      close(pb)


    ## Combine all rows into a single data frame.
    result <- do.call("rbind", df_list)
    timer("[tm_data] results combined - loop", verbose)

  } else {
    ## TODO(R): This is the 'faster version' of the
    ##          original loop. Has nothing to do with C, but for
    ##          testing I let it run when useC = TRUE.

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

    ## Appending the remaining data from 'data'
    for (n in names(data)) result[[n]] <- rep(data[[n]], response_values + 1)

    result <- as.data.frame(result)
    timer("[tm_data] faster vectorized version", verbose)
  }

  ## Attach the response column as an attribute.
  attr(result, "response") <- response

  return(result)
}

## Predict function.
# TODO(R): Adding useC option for testing; must be removed in the future.
tm_predict <- function(object, newdata,
  type = c("pdf", "cdf", "quantile", "pmax"),
  response = NULL, y = NULL, prob = 0.5, maxcounts = 1e+03,
  verbose = FALSE, theta_scaler = NULL, theta_vars = NULL,
  factor = FALSE, ncores = NULL, useC = FALSE)
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
  nd <- tm_data(newdata, response = response, verbose = verbose, useC = useC)
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
  p <- as.numeric(predict(object, newdata = nd, type = what))

  ## Extract unique indices
  ui <- unique(nd$index)
  prob <- rep(prob, length(ui))

  if (useC) {
    ## Ensure we hand over the correct thing to C
    if (type == "quantile") {
        stopifnot(is.numeric(prob), length(prob) == length(ui),
                  all(!is.na(prob)), all(prob >= 0 & prob <= 1))
    } else if (is.null(prob)) {
        prob <- 42.0 # dummy value for C (not used if type != 'quantile')
    }
    probs <- .Call(C_tm_predict, ui, nd$index, p, type = type, prob = prob, ncores = ncores);
  }

  # -------------------
  # Original R version, TODO(R): May be removed in the future
  if (!useC) {
    probs <- numeric(length(ui))
    for (j in ui) {

      pj <- p[nd$index == j]
      k <- length(pj)

      if (type == "pdf") {
        probs[j] <- (1 - pj[k]) * prod(pj[-k])
      }

      if (type == "cdf") {
        cj <- 1 - pj[1]
        if (length(pj) > 1) {
          for (jj in 2:length(pj))
            cj <- cj + (1 - pj[jj]) * prod(pj[1:(jj - 1)])
        }
        probs[j] <- cj
      }

      if (type == "quantile") {
        cj <- 1 - pj[1]
        if (any(is.na(pj))) {
          probs[j] <- NA_integer_
        } else if (cj >= prob[j]) {
          probs[j] <- 0L
        } else {
          if (length(pj) > 1) {
            for (jj in 2:length(pj)) {
              cj <- cj + (1 - pj[jj]) * prod(pj[1:(jj - 1)])
              if (cj >= prob[j]) {
                probs[j] <- jj - 1L
                break
              }
            }
          } else {
            probs[j] <- 0L
          }
        }
      }

      if (type == "pmax") {
        # That is my count
        cj <- numeric(k)
        cj[1] <- 1 - pj[1]
        if (length(pj) > 1) {
          for (jj in 2:length(pj))
            cj[jj] <- (1 - pj[jj]) * prod(pj[1:(jj - 1)])
        }
        probs[j] <- if (all(is.na(cj))) NA_integer_  else which.max(cj) - 1L
      }
    }
  } ## end !useC (R version)

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

timer <- function(msg = NULL, verbose = TRUE) {
    if (is.null(msg) || !exists("ttotal")) ttotal <<- Sys.time()
    if (!is.null(msg) && "treto" %in% ls(envir = .GlobalEnv)) {
        t <- as.numeric(Sys.time() - treto, units = "secs")
        tt <- as.numeric(Sys.time() - ttotal, units = "secs")
        if (verbose)
            message(sprintf("[timing] Elapsed  %8.1f ms  - %8.1f total     (%s)",
                            t * 1000, tt * 1000, msg))
    }
    treto <<- Sys.time()
}
timer(NULL)

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
# TODO(R): Adding useC option for testing; must be removed in the future.
tm <- function(formula, data, subset, na.action,
  engine = "bam", scale.x = FALSE, breaks = NULL,
  model = TRUE, ncores = NULL, verbose = FALSE, useC = FALSE, ...)
{
  if (!is.null(breaks)) breaks <- as.numeric(breaks)

  ## Staying sane
  stopifnot(
    "'ncores' must be NULL or numeric" = is.null(ncores) || is.numeric(ncores),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose)
  )
  ncores <- tm_get_number_of_cores(ncores, verbose = verbose)

timer(NULL)
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

  timer("building model frame", verbose)

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
  timer("discretize response", verbose)

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
  timer("data scaling", verbose)

  ## Max. counts.
  ymax <- max(mf[[rn]], na.rm = TRUE)
  k <- min(c(ymax - 1L, 20L))
  timer("find min/max", verbose)

  ## Transform data.
  tmf <- tm_data(mf, response = rn, useC = useC, verbose = verbose)
  timer("transforming data (tm_df)", verbose)

  if (!is.null(scaler)) {
    scaler$theta <- list("mean" = mean(tmf$theta), "sd" = sd(tmf$theta))
    tmf$theta <- (tmf$theta - scaler$theta$mean) / scaler$theta$sd
    timer("scaling theta", verbose)
  }

  if (length(tv)) {
    for (j in tv) {
      i <- as.integer(gsub("theta", "", j))
      tmf[[j]] <- as.integer(tmf$theta == i)
    }
    timer(paste("loop over tv (length ", length(tv), ")"), verbose)
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
  timer("updated formula", verbose)

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
  timer("estimation", verbose)
  options("warn" = warn)

  ## Additional info.
  rval$response <- rn
  rval$model.frame <- mf
  rval$scaler <- scaler
  rval$maxcounts <- max(mf[[rn]])
  rval$theta_vars <- tv
  rval$factor <- isTRUE(list(...)$factor)

  if (inherits(rval$model, "nnet")) {
    p <- predict(rval$model, type = "raw", useC = useC)
  } else {
    p <- predict(rval$model, type = "response", useC = useC)
  }
  timer("prediction", verbose)

  ## Remove model frame.
  if (!model)
    rval$model$model <- NULL

  ## Compute probabilities.
  ui <- unique(tmf$index)
  probs <- cprobs <- numeric(length(ui))

  if (useC) {
    ## c_tm_predict_pdfcdf returns a list with PDF and CDF, calculating
    ## both simultanously in C to improve speed.
    tmp    <- .Call(C_tm_predict_pdfcdf, uidx = ui, idx = tmf$index, p = p, ncores = ncores)
    probs  <- tmp$pdf
    cprobs <- tmp$cdf
    rm(tmp)
    timer("calculating CDF - C", verbose)
  } else {
    for (j in ui) {
      pj <- p[tmf$index == j]
      k <- length(pj)
      probs[j] <- (1 - pj[k]) * prod(pj[-k])

      cj <- numeric(k)
      cj[1] <- 1 - pj[1]
      if (length(pj) > 1) {
        for (jj in 2:length(pj))
          cj[jj] <- (1 - pj[jj]) * prod(pj[1:(jj - 1)])
      }
      cprobs[j] <- sum(cj)
    }
    timer("calculating CDF - R", verbose)
  }

  eps <- abs(.Machine$double.eps)
  probs[probs  < eps]      <- eps
  probs[probs  > 1 - eps]  <- 1 - eps
  cprobs[cprobs < eps]     <- eps
  cprobs[cprobs > 1 - eps] <- 1 - eps

  rval$probs <- data.frame("pdf" = probs, "cdf" = cprobs)
  timer("building rval", verbose)

  ## If binning.
  if (bin.y) {
    rval$bins   <- bins
    rval$ym     <- ym
    rval$yc_tab <- table(yc)
    rval$breaks <- breaks
    timer("ended if bin.y", verbose)
  }

  ## Assign class.
  class(rval) <- "tm"

  return(rval)
}

## TODO(R): Added useC option, remove once we streamlined this. Both in the
##          method call as well as in the residuals() call further down.
## Plotting method.
plot.tm <- function(x, which = "effects", spar = TRUE, k = 5, useC = FALSE, ...)
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
    resids <- cbind(resids, residuals(x, newdata = list(...)$newdata, useC = useC))
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
formula.tm <- function(object, ...)
{
  formula(object$model)
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


## TODO(R): Remove useC option once we streamline that; as input to the
##          method as well as in the predict call further down.
## Rootogram method.
rootogram.tm <- function(object, newdata = NULL, plot = TRUE,
  width = 0.9, style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), expected = TRUE, confint = TRUE,
  ref = TRUE, K = NULL, xlab = NULL, ylab = NULL, main = NULL, useC = FALSE, ...)
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
    p <- cbind(p, predict(object, newdata = newdata, type = "pdf", y = j, useC = useC))
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

