## Main paper: https://link.springer.com/article/10.1007/s10260-021-00558-6





#' Transition Models Estimation
#'
#' Fits transition models to count or continuous response data. The method leverages
#' transition probabilities to construct flexible probabilistic models without assuming
#' a fixed distribution for the response variable. Estimation relies on \pkg{mgcv}
#' infrastructure for GLM-type binary modeling.
#'
#' @param formula A GAM formula, see [mgcv::gam()]. This formula should specify
#'        the covariates and their interactions, as well as any special terms involving `theta`
#'        for transition modeling.
#' @param data A data frame containing the data for modeling.
#' @param subset An optional vector specifying a subset of observations to be used in the
#'        fitting process.
#' @param na.action A function that determines how `NA` values in the data should be handled.
#' @param engine Character string specifying the estimation engine. Options include
#'        `"bam"`, `"gam"`, `"nnet"`, or `"glmnet"` (experimental).
#'        Default is `"bam"`.
#' @param scale.x Logical value indicating whether covariates should be scaled before estimation.
#' @param breaks Controls the splitting of continuous responses into intervals to mimic a
#'        count response. If a single number is provided, equidistant intervals are used. If a
#'        numeric vector is provided, it is used directly as the breakpoints. This argument is
#'        critical for continuous responses.
#' @param model Logical value indicating whether the model frame should be included in
#'        the return object.
#' @param ncores `NULL` (default) or single numeric. See 'OpenMP' for more information.
#' @param verbose Logical value indicating whether progress and diagnostic messages
#'        should be printed. Default is `FALSE`.
#' @param \dots Additional arguments to be passed to the estimation engine.
#' @param formula Object of class `transitreg`.
#'
#' @details
#' The function transforms the input data using [transitreg_data()] to a format
#' compatible with transition models. Estimation relies on binary GLM-type techniques
#' to model conditional transition probabilities. Note that the `theta` variable, representing
#' the current transition level, must be included in the model by the user, please
#' see the examples. Additional transition-specific covariates, such as those addressing excess
#' zeros, can be included by adding `theta0`, `theta1`, etc., to the formula.
#'
#' For continuous responses, the `breaks` argument specifies the intervals for
#' discretizing the response, enabling the application of count-based transition models.
#'
#' @section OpenMP:
#' For improved performance `"transitreg"` is partially implemented in `C`, making
#' use of the OpenMP library for parallelization. The argument `ncores` allows
#' the user to control how many cores to be used for the `C` functions shipped
#' with this package. By default, the number of cores will be detected
#' automatically (`ncores = NULL`). The following rules are applied:
#'
#' * If OpenMP is not available, the number of cores is set to `1L` in all cases.
#' * If `ncores = NULL` the total number of cores (`N`) is detected
#'         automatically, and `ncores` is set to `N - 2L`.
#' * If `ncores < 1L`, one core will be used (single-core processing).
#' * If `ncores > 1L` all available cores will be utilized.
#'
#' @return
#' An object of class `"transitreg"`, which includes the following components:
#' * Fitted model results compatible with \pkg{mgcv}-style outputs.
#' * Methods for `plot`, `summary`, `residuals`, and `predict`.
#' * Model diagnostics and transformation details.
#'
#' See [transitreg_detect_cores()] for some more details.
#'
#' @seealso [transitreg_data()], [transitreg_dist()], [mgcv::gam()]
#'
#' @examples
#' ## Example 1: Count data.
#' set.seed(123)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x)))
#'
#' ## Fit transition count response model.
#' b <- transitreg(y ~ s(theta) + s(x))
#'
#' ## GAM summary.
#' summary(b)
#'
#' ## GAM coefficients
#' coef(b)
#'
#' ## Effect plots.
#' plot(b)
#'
#' ## Quantile residuals.
#' plot(b, which = 2:4)
#'
#' ## Predictions and plotting.
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#' fit <- cbind(
#'   predict(b, nd, type = "pmax"),
#'   predict(b, nd, type = "quantile", p = 0.05/2),
#'   predict(b, nd, type = "quantile", p = 0.5),
#'   predict(b, nd, type = "quantile", p = 1 - 0.05/2)
#' )
#'
#' ## Plot data and fitted counts.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", lty = 1, lwd = 2,
#'   col = c(2, 4, 4, 4), add = TRUE)
#'
#' ## Example 2: Continuous response.
#' set.seed(123)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- sin(x) + rnorm(n, sd = exp(-1 + cos(x)))
#'
#' ## Fit model with continuous response.
#' b <- transitreg(y ~ s(theta) + s(x) + te(x, theta), breaks = 200)
#'
#' ## Predictions and plotting
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#' fit <- cbind(
#'   predict(b, nd, type = "pmax"),
#'   predict(b, nd, type = "quantile", p = 0.05/2),
#'   predict(b, nd, type = "quantile", p = 0.5),
#'   predict(b, nd, type = "quantile", p = 1 - 0.05/2)
#' )
#'
#' ## Plot data and fitted curves.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", lty = 1, lwd = 2,
#'   col = c(2, 4, 4, 4), add = TRUE)
#'
#' ## Example 3: Count response with neural network.
#' set.seed(123)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x)))
#'
#' ## Fit NN transition count response model.
#' b <- transitreg(y ~ theta + x, scale.x = TRUE, engine = "nnet",
#'                 size = 5, maxit = 1000, decay = 0.001)
#'
#' ## Predictions and plotting.
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#' fit <- cbind(
#'   predict(b, nd, type = "pmax"),
#'   predict(b, nd, type = "quantile", p = 0.05/2),
#'   predict(b, nd, type = "quantile", p = 0.5),
#'   predict(b, nd, type = "quantile", p = 1 - 0.05/2)
#' )
#'
#' ## Plot data and fitted counts.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", lty = 1, lwd = 2,
#'   col = c(2, 4, 4, 4), add = TRUE)
#'
#' @keywords models regression
#'
#' @importFrom stats as.formula binomial predict sd update
#' @importFrom mgcv gam bam
#'
#' @author Niki
#' @export
transitreg <- function(formula, data, subset, na.action,
                       engine = "bam", scale.x = FALSE, breaks = NULL,
                       model = TRUE, ncores = NULL, verbose = FALSE, ...) {

  if (!is.null(breaks)) breaks <- as.numeric(breaks)

  ## Staying sane
  stopifnot(
    "'ncores' must be NULL or numeric" = is.null(ncores) || is.numeric(ncores),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose)
  )
  ncores <- transitreg_get_number_of_cores(ncores, verbose = verbose)

  cl <- match.call()

  ## Evaluate the model frame
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
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
  response <- response_name(formula)

  yscale <- NULL

  ## Discretize response?
  if (bin.y <- !is.null(breaks)) {
    if (length(breaks) == 1L) {
      breaks <- make_breaks(model.response(mf), breaks = breaks)
    }
    bins <- length(breaks) - 1L

    ## Discretize numeric response into counts.
    yc <- num2bin(model.response(mf), breaks)
    ym <- (breaks[-1] + breaks[-length(breaks)]) / 2

    lower <- list(...)$lower
    upper <- list(...)$upper

    if (!is.null(lower))
      ym[ym < lower] <- lower
    if (!is.null(upper))
      ym[ym > upper] <- upper

    mf[[response]] <- yc
  } else {
    ## For the model we need an integer response; check if the response
    ## is integer. If not, break.
    if (!all(mf[[response]] %% 1 < sqrt(.Machine$double.eps)) &&
        all(mf[[response]] > sqrt(.Machine$double.eps)))
        stop("Response is not looking like count data (integers); binning via 'breaks' is required.")
    mf[[response]] <- as.integer(round(mf[[response]], 1))

    bins <- max(mf[[response]])
    ## Setting highest 'bin' to max(response) * mp (multiplier)
    mp     <- if (bins <= 10) { 2 } else if (bins <= 100) { 1.5 } else { 1.2 }
    bins   <- as.integer(ceiling(bins * mp))
    breaks <- seq_len(bins + 1) - 1.5 # 'Integer' bins
    if (verbose) message("Response considered to be count data, using max count ", bins - 1)
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
  ymax <- max(mf[[response]], na.rm = TRUE)
  k <- min(c(ymax - 1L, 20L))

  ## Transform data.
  tmf <- transitreg_data(mf, response = response, verbose = verbose)

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
    if (length(attr(terms(rval$new_formula), "term.labels")) == 0L)
      stop("Intercept only model (", format(formula), ") with engine = \"", engine, "\" not allowed.")
    rval$model <- bam(rval$new_formula, data = tmf, family = binomial, discrete = TRUE)
  } else if (engine == "gam") {
    rval$model <- gam(rval$new_formula, data = tmf, family = binomial, ...)
  } else if (engine == "nnet") {
    rval$model <- nnet::nnet(rval$new_formula, data = tmf, ...)
  } else if (engine == "glmnet") {
    rval$model <- transitreg_glmnet(rval$new_formula, data = tmf, ...)
  }
  options("warn" = warn)

  ## Additional info.
  rval$response    <- response # Name of the response variable
  rval$model.frame <- mf       # Model frame (with 'binned' response)
  rval$scaler      <- scaler
  rval$maxcounts   <- max(mf[[response]]) # Highest count in response
  rval$theta_vars  <- tv
  rval$factor      <- isTRUE(list(...)$factor)

  if (inherits(rval$model, "nnet")) {
    tp <- predict(rval$model, type = "raw")
  } else {
    tp <- predict(rval$model, type = "response")
  }

  ## Remove model frame.
  if (!model)
    rval$model$model <- NULL

  ## Compute probabilities.
  ui <- unique(tmf$index)
  probs <- cprobs <- numeric(length(ui))

  ## c_transitreg_predict_pdfcdf returns a list with PDF and CDF, calculating
  ## both simultanously in C to improve speed.
  args <- list(uidx = ui, idx = tmf$index,
               tp = tp, y = mf[[response]], breaks = breaks, ncores = ncores)

  ## Calling C
  args <- check_args_for_treg_predict_pdfcdf(args)
  tmp <- do.call(function(...) .Call("treg_predict_pdfcdf", ...), args)

  ## Fixing values close to 0/1
  tmp$pdf[tmp$pdf < 1e-15]    <- 1e-15
  tmp$cdf[tmp$cdf < 1e-15]    <- 1e-15
  tmp$cdf[tmp$cdf > 0.999999] <- 0.999999

  rval$probs <- as.data.frame(tmp)
  rm(tmp)

  ## We always store 'bins'. In case of bin.y this is he number
  ## of bins (i.e., length(breaks) - 1), if !bin.y this is the higest
  ## integer of the original response.
  rval$bins <- bins

  ## If binning.
  if (bin.y) {
    rval$breaks <- breaks
    rval$ym     <- ym
    rval$yc_tab <- table(yc)
  }

  ## Assign class.
  class(rval) <- "transitreg"

  return(rval)
}


#' Transition Model Data Preparer
#'
#' Transforms a data frame into the format required for estimating transition models.
#' The function generates binary response data for fitting GLM-type models,
#' specifically designed to represent transitions between counts or states in
#' a probabilistic framework.
#'
#' @param data A data frame containing the raw input data.
#' @param response Character string specifying the name of the response variable
#'        to be used for transition modeling. This variable must represent
#'        counts or categorical states.
#' @param newresponse `NULL` or integer. New response vector to overwrite the one
#'        in `data`.
#' @param verbose Logical value indicating whether information about the transformation
#'        process should be printed to the console. Default is `TRUE`.
#'
#' @details
#' Transition models focus on modeling the conditional probabilities of transitions
#' between states or counts. This function converts the input data into a long format
#' suitable for such models. Each row in the resulting data frame corresponds to a
#' binary transition indicator, representing whether a transition to a higher category
#' occurred. For details on the modeling framework, see Berger and Tutz (2021).
#'
#' @return
#'   A transformed data frame in the long format. Each row represents a binary transition
#'   indicator (`Y`) for the response variable. Additional columns in the output include:
#'
#' * `index`: The original row index from the input data.
#' * `Y`: The binary indicator for whether a transition to a higher
#'       category occurred.
#' * `theta`: The level corresponding to the current transition.
#'
#' This format is required for fitting transition models using GLM or GAM frameworks.
#' For instance, a response variable with a value of 3 will generate rows with
#' transitions up to its value (0, 1, 2, and 3).
#'
#' @seealso [transitreg()], [transitreg_dist()]
#'
#' @examples
#' ## Raw data frame.
#' d <- data.frame(
#'   "id" = 1:5,
#'   "counts" = c(1, 0, 2, 3, 1),
#'   "x" = 1:5 * 10
#' )
#'
#' ## Transformed data frame.
#' dtm <- transitreg_data(d, response = "counts", verbose = TRUE)
#' print(dtm)
#'
#' @keywords data transformation
#' @concept transition
#' @concept modeling
#'
#' @importFrom stats setNames
#'
#' @author Niki
transitreg_data <- function(data, response = NULL, newresponse = NULL, verbose = TRUE) {

  stopifnot(
    "'response' must be NULL or a character of length 1" =
        is.null(response) || (is.character(response) && length(response) == 1L),
    "'newresponse' must be NULL or a numeric vector with length > 0" =
        is.null(newresponse) || (is.numeric(newresponse) && length(newresponse > 0)),
    "'verbose' must be TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose),
    "'newresponse' must be NULL or bins (0, 1, 2, ...)" =
        is.null(newresponse) || (is.integer(newresponse) && all(newresponse >= 0))
  )

  ## Ensure data is a data frame.
  if (!is.data.frame(data))
    data <- as.data.frame(data)

  ## Determine response column if not specified by user.
  if (is.null(response))
    response <- names(data)[1L]

  ## If newresponse is give, it must be length == 1 or the same length
  ## as the number of observations in 'data'. Required as we 'blow up'
  ## the data for the binary model.
  if (!is.null(newresponse)) {
      if (length(newresponse) == 1L) {
          newresponse <- rep(newresponse, nrow(data))
      } else if (length(newresponse) != nrow(data)) {
          stop("'newresponse' must be of length one (recycled for all observations) ",
               "or a vector of the same length as number of observations (rows) ",
               "in 'data'.")
      }
  }

  ## If 'response' is provided in 'data' and an additional 'newresponse'
  ## is provided, we will replace 'data[[response]]' with the newresponse.
  if (response %in% names(data) & !is.null(newresponse)) {
    warning("'newresponse' will overwrite response provided in 'data'.")
    data[[response]] <- newresponse
  } else if (!response %in% names(data)) {
    ## 'newresponse' must be provided. If missing, stop, else append.
    if (is.null(newresponse))
      stop("'response' not contained in 'data', must be provided via 'newresponse'.")
    data <- cbind(setNames(data.frame(newresponse), response), data)
  }
  rm(newresponse) # Delete object, no longer needed

  ## If any missing value in response: stop
  if (any(is.na(data[[response]])))
      stop("NA values in response data!")

  ## response$values must all be bin indices, so integers >= 0
  check <- all(data[[response]] > -sqrt(.Machine$double.eps) |
               abs(data[[response]] %% 1) > sqrt(.Machine$double.eps))
  if (!check)
    stop("The response must be bin indices, so integers in the range of {0, Inf}.")
  data[[response]] <- as.integer(data[[response]])

  ## Setting up the new data.frame with (pseudo-)bins

  ## Define length of vectors in list
  ## Sum of response indices + nrow(data), the latter to account for the
  ## additional "0" bin.
  nout <- sum(data[[response]]) + nrow(data)

  ## ------ building transitreg data -------

  ## Names of list elements.
  names_out <- c("index", "Y", "theta", names(data))

  ## Building index vector; each observation 1:nrow(data) gets its unique index
  result <- list()
  result$index <- rep(seq_len(nrow(data)), times = data[[response]] + 1L)

  ## Creating Y; always 1 except for the last entry per index.
  fn_get_Y <- function(nout, resp) {
    res <- rep(1L, nout); res[cumsum(resp + 1)] <- 0L; return(res)
  }
  result$Y <- fn_get_Y(nout, data[[response]])

  ## Creating theta; a sequence from zero to the
  ## response_value for each index. The following
  ## Two lines create this sequence of sequences.
  fn_get_theta <- function(nout, resp, idx) {
    resettozero <- c(0, which(diff(idx) > 0))
    return(seq_len(nout) - rep(resettozero, resp + 1) - 1)
  }
  result$theta <- fn_get_theta(nout, data[[response]], result$index)

  ## Appending the remaining data from 'data'.
  for (n in names(data)) {
      ## If data[[n]] is a simple vector
      if (!is.matrix(data[[n]])) {
        result[[n]] <- rep(data[[n]], data[[response]] + 1)
      ## Else create matrix
      } else {
        result[[n]] <- matrix(rep(data[[n]], rep(data[[response]] + 1, ncol(data[[n]]))),
                              ncol = ncol(data[[n]]),
                              dimnames = list(NULL, colnames(data[[n]])))
      }
  }

  result <- as.data.frame(result)

  ## Attach the response column as an attribute.
  attr(result, "response") <- response

  return(result)
}





# Helper function for predictions on a transitreg model object.
transitreg_predict <- function(object, newdata = NULL,
        type = c("pdf", "cdf", "quantile", "pmax", "tp"), y = NULL, prob = NULL,
        elementwise = NULL, maxcounts = 1e+03,
        verbose = FALSE, theta_scaler = NULL, theta_vars = NULL,
        factor = FALSE, ncores = NULL) {

  ## This method is built on the 'transitreg' model as we need access
  ## to the original model.frame, breaks, ... although the prediction (TPs)
  ## will be done on object$model, the internal probability model.
  stopifnot(
    "'object' must be of class 'transitreg'" = inherits(object, "transitreg")
  )

  ## Pre-processing inputs
  type <- tolower(type)
  type <- match.arg(type)

  if (type == "pmax")
      stop("TODO(R): Partially implemented, but not yet tested.")

  if (is.null(ncores))
    ncores <- transitreg_get_number_of_cores(ncores, verbose = verbose)

  ## Staying sane
  stopifnot(
    "'newdata' must be NULL or data.frame" = 
        is.null(newdata) || is.data.frame(newdata),
    "'y' must be NULL or a vector" = is.null(y) || is.vector(y),
    "'elementwise' must be NULL, TRUE, or FALSE" =
        is.null(elementwise) || isTRUE(elementwise) || isFALSE(elementwise),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose),
    "'factor' must be logical TRUE or FALSE" = isTRUE(factor) || isFALSE(factor),
    "'ncores' must be NULL or numeric >= 1" = is.null(ncores) || (ncores >= 1)
  )
  ## TODO(R) Not all arguments are checked above

  ## Depending on 'type' we need to ensure y/prob is set correctly.
  ## Quantile:
  ## - prob must be not NULL, all values must be in [0, 1].
  ## PDF/CDF:
  ## - If 'newdata = NULL' we take the model frame, so y can be missing.
  ## - Else y must be not NULL, all values must be within support of 'breaks'.
  if (type == "quantile") {
    stopifnot(
      "for 'type = \"quantile\"' argument 'prob' must be set" = !is.null(prob),
      "'prob' must be numeric of length > 0" = is.numeric(prob) && length(prob) > 0,
      "missing values in 'prob' not allowed" = all(!is.na(prob)),
      "'prob' must be in [0, 1]" = all(prob >= 0.0) && all(prob <= 1.0)
    )
  } else if (type %in% c("cdf", "pdf") & !is.null(newdata)) {
    stopifnot(
      "for 'type = \"cdf\" or \"pdf\"' argument 'y' must be set if 'newdata' is given" =
          !is.null(y),
      "'y' must be numeric of length > 0" = is.numeric(y) && length(y) > 0,
      "missing values in 'y' not allowed" = all(!is.na(y))
    )
    # Discrete distribution? y must be in range [0, object$bins].
    if (is.null(object$breaks)) {
      y <- as.integer(y)
      if (!all(y >= 0 & y <= object$bins - 1L))
          stop("elements in 'y' must be in {0, ..., ", object$bins - 1L, "}")
    # Else continuous
    } else {
      if (!all(y >= min(object$breaks) & y <= max(object$breaks)))
          stop(sprintf("elements in 'y' must be in range [%s, %s]",
                       format(min(object$breaks), format(max(object$breaks)))))
      ## If 'y' was set by the user, we must convert all values in 'y' from
      ## it's original numeric value to it's pseudo-observation (bin; int).
      if (!is.null(y)) y <- num2bin(y, object$breaks)
    }
  }

  ## If newdata is empty, we are taking the model frame.
  ## In addition, we take the pseudo-response from the model.frame
  ## in case type %in% c("cdf", "pdf") and the user did not provide any 'y'.
  ## Thus, we will evaluate the cdf/pdf a the observation used for fitting the model.
  if (is.null(newdata)) {
      newdata <- model.frame(object)
      if (type %in% c("cdf", "pdf") && is.null(y))
          y <- newdata[[object$response]] ## Integer in 0, 1, 2, ... Nbreaks
  } else {
    ## Applying scaler if needed; standardize data (except theta)
    if (!is.null(object$scaler)) {
        for (j in names(object$scaler)) {
            if (j != "theta") {
                newdata[[j]] <- (newdata[[j]] - object$scaler[[j]]$mean) / object$scaler[[j]]$sd
            }
        }
    }
  }

  ## Scaling (standardizing) theta if requested
  if (!is.null(theta_scaler))
    newdata$theta <- (newdata$theta - theta_scaler$mean) / theta_scaler$sd

  ## Guessing 'elementwise' if is NULL
  ##
  ## TRUE:  If set TRUE, the prdiction is done elementwise. In other words:
  ##        element i in y/prob corresponds to observation newdata[i, ] and
  ##        for each observation in newdata we will get one result in return.
  ##
  ## FALSE: Each observation in newdata will be evaluated at every value in
  ##        y/prob. The return will therefore be a matrix of dimension
  ##        ncol(newdata) times (length(y) | length(probs)) depending on type.
  if (is.null(elementwise)) {
      if (type == "quantile") {
          elementwise <- length(prob) == 1L | length(prob) == nrow(newdata)
          ## Extending prob and sorting if !elementwise
          if (elementwise & length(prob) == 1L)
              prob <- rep(prob, length.out = nrow(newdata))
      } else if (type %in% c("cdf", "pdf")) {
          elementwise <- length(y) == 1L | length(y) == nrow(newdata)
          ## Extending y and sorting if !elementwise
          if (elementwise & length(y) == 1L)
              y <- rep(y, length.out = nrow(newdata))
      } else if (type == "pmax") {
          elementwise <- TRUE # for 'pmax' elementwise is always TRUE
      } else if (type == "tp") {
          NULL
      } else {
          stop("Mode for type = \"", type, "\" must be implemented!")
      }
  }

  ## Setting dummy values (required by C later on)
  if (type == "quantile") {
      ## Sorting 'prob'. This is important for the .C routine!
      if (!elementwise) prob <- sort(unique(prob))
      y    <- NA_integer_ ## Dummy value required for .C call
  } else if (type %in% c("cdf", "pdf")) {
      ## Sorting 'y'. This is important for the .C routine!
      if (!elementwise) y <- sort(unique(y))
      prob <- NA_real_    ## Dummy value required for .C call
  } else {
      y    <- NA_integer_ ## Dummy value required for .C call
      prob <- NA_real_    ## Dummy value required for .C call
  }


  ## Setting up 'newresponse'.
  ## Quantile, pmax:
  ##  - Setting the pseudo-response to the highest bin. This ensures
  ##    that the entire Transition distribution is evaluated.
  ## CDF/PDF:
  ##  - If elementwise = TRUE: We only need to evaluate each distribution
  ##    up to 'y[i]'.
  ##  - If elementwise = FALSE: We must evaluate each distribution up to
  ##    max(y).
  if (type %in% c("quantile", "pmax", "tp")) {
    ## object$bins - 1 as we start with bin '0' again.
    newresponse <- rep(object$bins - 1L, nrow(newdata))
  } else {
    newresponse <- if (elementwise) y else rep(max(y), nrow(newdata))
  }

  ## Preparing data
  newdata <- transitreg_data(newdata, response = object$response,
                             newresponse = newresponse, verbose = verbose)

  if (factor)
    newdata$theta <- as.factor(newdata$theta)

  if (!is.null(theta_vars) && length(theta_vars) > 0L) {
    for (j in theta_vars) {
      i <- as.integer(gsub("theta", "", j))
      newdata[[j]] <- as.integer(newdata$theta == i)
    }
  }

  ## Specify argument for generic prediction method called below
  what <- switch(class(object$model)[1L],
    "bam"  = "response",
    "gam"  = "response",
    "nnet" = "raw"
  )
  tp <- as.numeric(predict(object$model, newdata = newdata, type = what))

  ## If 'type = "tp"' (transition probabilities) we already have our
  ## result. Conver to matrix, and return.
  if (type == "tp") {
      tmp <- list(NULL, paste0("tp_", seq_len(object$bins) - 1))
      return(matrix(tp, byrow = TRUE, ncol = object$bins, dimnames = tmp))
  }

  ## Extract unique indices
  ui   <- unique(newdata$index)

  ## If object$breaks is NULL, we have discrete bins (e.g., count data).
  if (is.null(object$breaks)) {
    discrete <- rep(TRUE, length(ui))
    breaks   <- seq(-0.5, by = 1.0, length.out = object$bins + 1)
  } else {
    discrete <- rep(TRUE, length(ui))
    breaks   <- object$breaks
  }

  ## Setting up arguments for the .C call
  args <- list(uidx   = ui,                # int; Unique distribution index (int)
               idx    = newdata$index,     # int; Index vector (int)
               tp     = tp,                # num; Transition probabilities
               breaks = breaks,            # num; Point intersections of bins
               y      = y,                 # int; Response y, used for 'cdf/pdf'
               prob   = prob,              # num; Probabilities (used for 'quantile')
               type   = type,              # str; to predict/calculate
               ncores = ncores,            # int; Number of cores to be used (OpenMP)
               elementwise = elementwise,  # Elementwise (one prob or y per ui)
               discrete    = discrete)     # Discrete distribution?

  # Calling C
  args <- check_args_for_treg_predict(args)
  res  <- do.call(function(...) .Call("treg_predict", ...), args)

  # If 'ementwise = FALSE' we get length(y)/length(prob) results per
  # observation and have to glue them back together into a matrix.
  if (!elementwise && type == "quantile") {
    res <- matrix(res, byrow = TRUE, ncol = length(prob),
                  dimnames = list(NULL, get_elementwise_colnames(prob, NULL)))
  } else if (!elementwise) {
    prefix <- if (type == "pdf") "d" else "p"
    res <- matrix(res, byrow = TRUE, ncol = length(y),
                  dimnames = list(NULL, get_elementwise_colnames(y, prefix)))
  }

  return(res)
}



#' Transition Model Probability Density Visualization
#'
#' Visualizes the probability density function (PDF) and raw count or continuous data
#' based on transition models estimated using [transitreg()]. This function provides
#' an intuitive way to understand the distribution of the modeled response.
#'
#' @param y A response vector or a formula specifying the relationship between the
#'        response and covariates. For count data, this is typically a vector of counts.
#'        For continuous data, this can be paired with the `breaks` argument to
#'        discretize the response.
#' @param data Optional. If `y` is a formula, this specifies the data frame to
#'        be used for model fitting.
#' @param \dots Additional arguments to be passed to [transitreg()], including
#'        settings for the estimation engine, formula, and other relevant parameters.
#'
#' @details
#' This function estimates and visualizes the underlying probability density function
#' (PDF) for count or continuous response data using transition models. For continuous
#' data, the response is discretized based on the `breaks` argument passed through
#' \dots
#'
#' The function supports visualizations for raw counts, zero-inflated data, and transformed
#' distributions, providing insights into the modeled distribution of the response variable.
#'
#' @return
#' An object of class `"transitreg"`, as described in [transitreg()]. This includes:
#'
#' * Fitted transition model details.
#' * Model diagnostics and parameters.
#' * Visualization-ready data for plotting PDFs or transformed distributions.
#'
#' @seealso [transitreg()], [transitreg_data()].
#'
#' @examples
#' ## Example 1: Count data.
#' set.seed(123)
#' n <- 3000
#' y <- rpois(n, 10)
#'
#' # Visualize PDF for count data.
#' transitreg_dist(y)
#'
#' ## Example 2: Zero-inflated data.
#' y <- c(y, rep(0, 500))
#'
#' ## Include a zero-inflation term.
#' transitreg_dist(y ~ s(theta) + theta0)
#'
#' ## Example 3: Continuous data.
#' set.seed(123)
#' n <- 1000
#' y <- rgamma(n, shape = 10, rate = 0.1)
#'
#' ## Visualize PDF for continuous data with discretization.
#' transitreg_dist(y, breaks = 50)
#'
#' @keywords distribution visualization
#'
#' @importFrom stats model.response predict
#' @importFrom grDevices rgb
#' @importFrom graphics barplot lines points
transitreg_dist <- function(y, data = NULL, ...) {
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
  b <- transitreg(f, data = data, ...)

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




#' @exportS3Method "[" transitreg
#' @author Reto
`[.transitreg` <- function(x, i, ..., drop = TRUE) {
    tp <- transitreg_predict(x, newdata = model.frame(x)[i, , drop = FALSE],
                            type = "tp")
    breaks <- if (is.null(x$breaks)) seq_len(x$bins + 1) - 1.5 else x$breaks
    return(Transition(tp, breaks))
}


#' Plot Method for Transition Model Fits
#'
#' Provides diagnostic and visualization plots for transition models fitted using
#' the [transitreg()] function. The method supports plotting effects and quantile
#' residual diagnostic plots.
#'
#' @param x An object of class `"transitreg"` resulting from a call to [transitreg()]
#' @param which A character string or integer specifying the type of plot(s) to
#'        generate (See 'Details').
#' @param spar Logical. If `TRUE`, multiple plots are arranged in a
#'        single window. Default is `TRUE`.
#' @param k Integer, TODO(N): Describe argument. Defaults to `5`.
#' @param \dots Additional arguments passed to the underlying plotting functions.
#'
#' @details
#' The function allows to control what to plot via the argument `which`.
#' Options include:
#'
#' * `"effects"` Plots the effects of the predictors on the response.
#'   Requires that the model is estimated by [mgcv::gam()]
#'   or [mgcv::bam()].
#' * `"hist-resid"` Plots a histogram of the qauntile residuals.
#' * `"qq-resid"` Generates a Q-Q plot of the quantile residuals.
#' * `"wp-resid"` Creates a worm plot of the quantile residuals.
#'
#' Multiple options can be specified as a character vector or numeric indices.
#'
#' The [plot.transitreg()] method provides flexible visualization options for
#' evaluating transition model fits. Users can choose to:
#'
#' * Visualize the effects of predictors on the response variable
#'   (if the model is a GAM, see [mgcv::gam()]).
#' * Evaluate quantile residuals through histograms, Q-Q plots, or worm plots.
#'
#' The `which` argument controls the type of plots generated. By default, the
#' `"effects"` plot is shown if the model supports it. Residual-based plots
#' (`"hist-resid"`, `"qq-resid"`, `"wp-resid"`) provide insights into model
#' calibration.
#'
#' @return Returns `NULL` invisibly. Generates plots as a side effect.
#'
#' @seealso [transitreg()], [transitreg_dist()], [transitreg_data()], [predict.transitreg()].
#'
#' @examples
#' ## Example: Fit a transition model and generate plots.
#' set.seed(123)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x)))
#' b <- transitreg(y ~ s(theta) + s(x))
#'
#' ## Plot effects.
#' plot(b, which = "effects")
#'
#' ## Plot residuals.
#' plot(b, which = c("hist-resid", "qq-resid"))
#'
#' ## Custom plot layout.
#' par(mfrow = c(2, 1))
#' plot(b, which = 3:4, spar = FALSE)
#'
#' @keywords methods models visualization
#'
#' @exportS3Method plot transitreg
#' @importFrom stats residuals
#' @importFrom graphics par
plot.transitreg <- function(x, which = "effects", spar = TRUE, k = 5, ...)
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

#' @exportS3Method summary transitreg
summary.transitreg <- function(object, ...) {
  summary(object$model)
}

#' @importFrom stats formula
#' @exportS3Method formula transitreg
formula.transitreg <- function(x, ...) {
  formula(x$model)
}

#' @importFrom stats coef
#' @exportS3Method coef transitreg
coef.transitreg <- function(object, ...) {
  coef(object$model)
}


#' @author Reto
#' @exportS3Method model.frame transitreg
model.frame.transitreg <- function(formula, ...) {
  formula$model.frame
}

#' @author Niki
#' @exportS3Method print transitreg
print.transitreg <- function(x, ...) {
  cat("Count Transition Model\n---")
  print(x$model)
}


#' @importFrom distributions3 prodist
#' @importFrom stats setNames
#' @importFrom utils head tail
#'
#' @author Reto
#' @exportS3Method prodist transitreg
#' @rdname transitreg
prodist.transitreg <- function(object, newdata = NULL, ...) {
    n <- nrow(object$model.frame)
    object[seq_len(n)]
}


#' @importFrom topmodels newresponse
#'
#' @author Reto
#' @rdname transitreg
#' @exportS3Method newresponse transitreg
newresponse.transitreg <- function(object, newdata = NULL, ...) {
    ## Response name
    yn <- object$response

    if (is.null(newdata)) {
        newdata <- object$model.frame
        newdata[[yn]] <- object$breaks[newdata[[yn]]]
    }

    if (is.null(newdata[[object$response]]))
        stop("response missing in newdata!")

    y <- setNames(data.frame(newdata[[object$response]]), yn)
    return(y)
}



#' Predict Method for Transition Model Fits
#'
#' Provides predictions for transition models fitted using the [transitreg()] function.
#' Predictions can be generated for the probability density function (PDF), cumulative
#' distribution function (CDF), maximum probability, or specific quantiles of the response
#' distribution.
#'
#'
#' @param object An object of class `transitreg` resulting from a call to [transitreg()].
#' @param newdata Optional. A data frame containing new predictor values.
#'        If not provided, predictions are made for the data used in fitting the model.
#' @param y Optional. A vector of response values for which the PDF or CDF should
#'        be computed. Required if `type` is `"pdf"` or `"cdf"`.
#' @param prob Optional. A numeric value specifying the quantile to compute when
#'        `type = "quantile"`. Default is `0.5` (median). If provided,
#'        argument `type` is set to `type = "quantile"`.
#' @param type Character. Specifies the type of prediction to return (see Section 'Details').
#' @param ncores `NULL` (default) or single numeric. See section 'OpenMP'
#'        of the [transitreg()] man page for more details.
#' @param elementwise Logical. Should each distribution in `object` be evaluated at
#'        all elements of `prob`/`y`? (`elementwise = FALSE`, yielding a matrix)? Or, if
#'        `x` and `probs` have the same length, should the evaluation be done element
#'        by element (`elementwise = TRUE` yielding a vector)? The default of `NULL`
#'        means that `elementwise = TRUE` is used if the lengths match and otherwise
#'        `elementwise = FALSE` is used.
#' @param verbose Logical, if `TRUE` few messages will be printed.
#' @param \dots Additional arguments passed to the prediction function.
#'
#' @details
#' The `predict.transitreg` method computes predictions based on the transition
#' model fit. Predictions can be made for the original training data or for new
#' data provided via `newdata`. The method also supports scaling of covariates if
#' scaling was applied during model fitting (argument `scale.x` in function
#' [transitreg()]).
#'
#' The argument `type` controls the return, the following types are allowed:
#'
#' * `"pdf"`: The predicted probability density function (PDF).
#' * `"cdf"`: The cumulative distribution function (CDF).
#' * `"pmax"`: The expected value of the response (maximum probability).
#' * `"quantile"`: The quantile of the response specified by `prob`.
#'
#' For `"pdf"` and `"cdf"`, the response values (`y`) must be provided unless the
#' model was fit with those values already included. For `"quantile"`, a specific
#' quantile(s) are computed based on `prob`.
#'
#' @return
#' Returns predictions of the specified type:
#'
#' * For `"pdf"` and `"cdf"`, a vector or matrix of probabilities evaluated at `y`.
#' * For `"pmax"`, the expected value of the response.
#' * For `"quantile"`, the quantile of the response distribution at the specified `prob`.
#'
#' @seealso [transitreg()], [transitreg_data()], [transitreg_dist()].
#'
#' @examples
#' ## Example: Predicting PDF and CDF.
#' set.seed(123)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x)))
#' b <- transitreg(y ~ s(theta) + s(x))
#'
#' ## Predict PDF and CDF.
#' p <- list()
#' p$pdf <- predict(b, type = "pdf", y = 3)
#' p$cdf <- predict(b, type = "cdf", y = 3)
#'
#' ## Predict maximum probability (expected value).
#' p$pmax <- predict(b, type = "pmax")
#'
#' ## Predict quantiles.
#' p$qu95 <- predict(b, prob = 0.95)
#'
#' print(head(as.data.frame(p)))
#'
#' ## Visualize predictions.
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#'
#' ## Predict quantiles.
#' qu <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 0.95, 0.99)
#' p <- lapply(qu, function(prob) {
#'   predict(b, newdata = nd, prob = prob)
#' })
#'
#' ## Plot data and fitted quantiles.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, do.call("cbind", p),
#'   type = "l", col = 4, lwd = 2, lty = 1,
#'   add = TRUE)
#'
#' @keywords methods model
#'
#' @importFrom stats model.frame
#'
#' @exportS3Method predict transitreg
#' @author Niki
predict.transitreg <- function(object, newdata = NULL, y = NULL, prob = NULL,
        type = c("pdf", "cdf", "quantile", "pmax", "tp"), ncores = NULL,
        elementwise = NULL, verbose = FALSE, ...) {

  type <- tolower(type)
  type <- match.arg(type)

  if (!is.null(prob))
    type <- "quantile"

  if (is.null(prob) && type == "quantile")
    prob <- 0.5

  ## Calling transitreg_predict to perform the actual prediction
  args <- list(object       = object,
               newdata      = newdata,
               type         = type,
               y            = y,
               prob         = prob,
               elementwise  = elementwise,
               maxcounts    = object$maxcounts,
               verbose      = verbose,
               theta_scaler = object$scaler$theta,
               theta_vars   = object$theta_vars,
               factor       = object$factor,
               ncores       = ncores)

  return(do.call(transitreg_predict, args))
}

#' @importFrom stats runif
#'
#' @author Niki
#' @rdname transitreg
#' @exportS3Method logLik transitreg
logLik.transitreg <- function(object, newdata = NULL, ...) {
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


