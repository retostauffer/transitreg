
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
#' @param censored `NULL` or one of `"left"`, `"right"`, or `"both"`. Specifies
#'        if the distribution is censored on one or both ends.
#' @param model Logical value indicating whether the model frame should be included in
#'        the return object.
#' @param ncores `NULL` (default) or single numeric. See 'OpenMP' for more information.
#' @param verbose Logical value indicating whether progress and diagnostic messages
#'        should be printed. Default is `FALSE`.
#' @param \dots Additional arguments to be passed to the estimation engine.
#' @param formula Object of class `transitreg`.
#'
#' @details
#' The function transforms the input data using [transitreg_tmf()] to a format
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
#' @seealso [transitreg_tmf()], [transitreg_dist()], [mgcv::gam()]
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
                       censored = NULL,
                       model = TRUE, ncores = NULL, verbose = FALSE, ...) {

  if (!is.null(breaks)) breaks <- as.numeric(breaks)
  if (!is.null(censored)) {
      censored <- tolower(censored)
      censored <- match.arg(censored, c("left", "right", "both"))
  }

  ## Staying sane
  stopifnot(
    "'ncores' must be NULL or numeric" = is.null(ncores) || is.numeric(ncores),
    "'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose),
    "'breaks' must be NULL, single positive numeric or numeric vector" =
        is.null(breaks) || (is.numeric(breaks) && length(breaks) == 1L && breaks > 0) || is.numeric(breaks) && length(breaks) > 0,

    "'scale.x' must be FALSE or TRUE" = isFALSE(scale.x) || isTRUE(scale.x)
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

  theta_vars <- all.vars(ff2)
  theta_vars <- grep("theta", theta_vars, value = TRUE)
  theta_vars <- theta_vars[theta_vars != "theta"]
  if (length(theta_vars)) {
    for (j in theta_vars)
      ff2 <- eval(parse(text = paste0("update(ff2, . ~ . -", j, ")")))
  }

  mf[["formula"]] <- ff2
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  ## Response name.
  yscale <- NULL


  ## Setting up empty list for return value
  rval <- list()
  rval$censored <- censored
  rval$response <- response_name(formula)
  rval$ymax     <- max(mf[[rval$response]])

  ## Setting up 'breaks and bins'
  if (is.numeric(breaks) && length(breaks) == 1L) {
      # Create, and store breaks and bins
      breaks <- rval$breaks <- make_breaks(mf[[rval$response]], breaks = breaks)
      rval$bins <- length(breaks) - 1L ## 10 breaks = 9 bins
  ## User-specified breaks, check if they span the required range
  } else if (is.numeric(breaks)) {
      tmp_bk <- range(breaks)
      tmp_y  <- range(mf[[rval$response]])
      if (tmp_bk[[1L]] > tmp_y[[1L]] || tmp_bk[[2L]] < tmp_y[[2L]])
          stop("Breaks do not cover the full range of the response \"", rval$response, "\".")
      # Store breaks and bins
      rval$bins   <- length(breaks) - 1L
      rval$breaks <- breaks
  # No breaks specified? In this case the response must be count data
  } else {
      # Check that the response is positive integers only.
      if (!all(mf[[rval$response]] >= 0L) ||
          !all(abs(mf[[rval$response]] %% 1) < .Machine$doluble.eps)) {
          stop("Response not count data. Breaks must be specified for binning.")
      }
      # There are no breaks, but bins
      ymax <- as.integer(max(mf[[rval$response]], na.rm = TRUE))
      tmp  <- ceiling(ymax * if (ymax <= 10) { 3 } else if (ymax <= 100) { 1.5 } else { 1.25 })
      rval$bins <- as.integer(tmp)
      # Will not be stored on 'rval' but used to convert data
      breaks <- seq.int(0L, rval$bins) - 0.5
      rm(tmp, ymax)
  }

  ## Transform data.
  tmf <- transitreg_tmf(mf,
                        response   = rval$response,
                        breaks     = breaks, # <- note: not rval$breaks
                        theta_vars = theta_vars,
                        scaler     = scale.x, verbose = verbose, ...)

  ## Response (as bins)
  y    <- num2bin(mf[[rval$response]], breaks = breaks)
  ymax <- max(y, na.rm = TRUE)

  ## Store scaler, returned as attribute on 'tmf' if used.
  rval$scaler <- if (scale.x) attr(tmf, "scaler") else NULL

  ## New formula.
  if (engine %in% c("gam", "bam")) {
    if (isTRUE(list(...)$factor)) {
      tmf$theta <- as.factor(tmf$theta)
      rval$new_formula <- update(formula, as.factor(Y) ~ theta + .)
    } else {
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

  ## Storing additional info on return object.
  rval$model.frame <- mf       # Model frame (with 'binned' response)
  rval$theta_vars  <- theta_vars
  rval$factor      <- isTRUE(list(...)$factor)

  # Highest count in response
  rval$maxcounts   <- num2bin(max(mf[[rval$response]]), breaks = breaks)

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
  censored <- if (is.null(rval$censored)) "not-censored" else rval$censored
  args <- list(uidx = ui, idx = tmf$index,
               tp = tp, y = y, breaks = breaks, ncores = ncores, censored = censored)

  ## Calling C
  args <- check_args_for_treg_predict_pdfcdf(args)
  tmp  <- do.call(function(...) .Call("treg_predict_pdfcdf", ...), args)

  ## Fixing values close to 0/1
  tmp$pdf[tmp$pdf < 1e-15]    <- 1e-15
  tmp$cdf[tmp$cdf < 1e-15]    <- 1e-15
  tmp$cdf[tmp$cdf > 0.999999] <- 0.999999

  rval$probs <- as.data.frame(tmp)
  rm(tmp)

  ## Assign class.
  class(rval) <- "transitreg"

  return(rval)
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

  ## Getting breaks from object
  breaks <- get_breaks(object)

  ## If 'newdata' is not set we take the existing model.frame; here the response
  ## is already stored as 'bin indices', so we do not have to call numb2bin.
  mf <- if (is.null(newdata)) model.frame(object) else newdata

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
          elementwise <- length(prob) == 1L | length(prob) == nrow(mf)
          ## Extending prob and sorting if !elementwise
          if (elementwise & length(prob) == 1L)
              prob <- rep(prob, length.out = nrow(mf))
      } else if (type %in% c("cdf", "pdf")) {
          elementwise <- is.null(y) || length(y) == 1L || length(y) == nrow(mf)
      } else if (type == "pmax") {
          elementwise <- TRUE # for 'pmax' elementwise is always TRUE
      } else if (type == "tp") {
          NULL
      } else {
          stop("Mode for type = \"", type, "\" must be implemented!")
      }
  }

  ## Setting up 'newresponse'.
  ## Quantile, pmax, tp:
  ##  - Setting response to "max bin mid" to ensure we calculate the
  ##    transition probabilities for _all_ bins.
  ## CDF/PDF:
  ##  - If elementwise = TRUE: We only need to evaluate each distribution
  ##    up to 'y[i]'.
  ##  - If elementwise = FALSE: We must evaluate each distribution up to
  ##    max(y).
  if (type %in% c("quantile", "pmax", "tp")) {
    mf[[object$response]] <- max(get_mids(object))
  } else {
    ## Else highest bin specified on 'y' (if set) or highest
    ## highest ever seen observation (yc_tab).
    if (isTRUE(elementwise) && !is.null(y)) {
        mf[[object$response]] <- rep(y, length.out = nrow(mf))
    } else if (isFALSE(elementwise) && is.null(y)) {
        mf[[object$response]] <- max(as.integer(names(object$yc_tab)))
    } else if (!is.null(y)) {
        mf[[object$response]] <- max(y)
    }
  }

  ## Creating 'transition model frame' for the prediction of the
  ## transition probabilities using the object$model (binary response model).
  tmf <- transitreg_tmf(mf,
                        response   = object$response,
                        breaks     = breaks,
                        theta_vars = object$theta_vars,
                        scaler     = object$scaler, verbose = verbose)

  ## Specify argument for generic prediction method called below
  what <- switch(class(object$model)[1L],
    "bam"  = "response",
    "gam"  = "response",
    "nnet" = "raw"
  )
  tp <- as.numeric(predict(object$model, newdata = tmf, type = what))

  ## If 'type = "tp"' (transition probabilities) we already have our
  ## result. Convert to matrix and return.
  if (type == "tp") {
      tmp <- list(NULL, paste0("tp_", seq_len(object$bins) - 1))
      return(matrix(tp, byrow = TRUE, ncol = object$bins, dimnames = tmp))
  }

  ## Extract unique indices
  ui   <- unique(tmf$index)

  ## Setting dummy values (required by C later on)
  if (type == "quantile") {
      ## Sorting 'prob'. This is important for the .C routine!
      probC <- if (!elementwise) sort(unique(prob)) else prob
      yC   <- NA_integer_ ## Dummy value required for .C call
  } else if (type %in% c("cdf", "pdf")) {
      ## Sorting 'y'. This is important for the .C routine!
      if (elementwise) {
          yC <- num2bin(mf[[object$response]], breaks = breaks)
      } else {
          yC <- num2bin(sort(unique(y)), breaks = breaks)
      }
      probC <- NA_real_    ## Dummy value required for .C call
  } else {
      yC    <- NA_integer_ ## Dummy value required for .C call
      probC <- NA_real_    ## Dummy value required for .C call
  }

  ## If object$breaks is NULL, we have discrete bins (e.g., count data).
  discrete <- rep(is.null(object$breaks), length(ui))
  censored <- if (is.null(object$censored)) "not-censored" else object$censored

  ## Setting up arguments for the .C call
  args <- list(uidx   = ui,                # int; Unique distribution index (int)
               idx    = tmf$index,         # int; Index vector (int)
               tp     = tp,                # num; Transition probabilities
               breaks = breaks,            # num; Point intersections of bins
               y      = yC,                # int; Response y, used for 'cdf/pdf'
               prob   = probC,             # num; Probabilities (used for 'quantile')
               type   = type,              # str; to predict/calculate
               ncores = ncores,            # int; Number of cores to be used (OpenMP)
               elementwise = elementwise,  # Elementwise (one prob or y per ui)
               discrete    = discrete,     # Discrete distribution?
               censored    = censored)

  # Calling C
  args <- check_args_for_treg_predict(args)
  res  <- do.call(function(...) .Call("treg_predict", ...), args)

  # If 'ementwise = FALSE' we get length(y)/length(prob) results per
  # observation and have to glue them back together into a matrix.
  if (!elementwise && type == "quantile") {
    res <- matrix(res, byrow = TRUE, ncol = length(prob),
                  dimnames = list(NULL, get_elementwise_colnames(prob, NULL)))
  } else if (!elementwise) {
    ## Convert from predicted 'bin' back to the numeric value.
    prefix <- if (type == "pdf") "d" else "p"
    ym  <- bin2num(y, object$breaks)
    res <- matrix(res, byrow = TRUE, ncol = length(y),
                  dimnames = list(NULL, get_elementwise_colnames(ym, prefix)))
  }

  return(res)
}


# Helper function, returns breaks of a transitreg model.
# If the model is a 'count data model' the breaks are not stored
# on the object, but here calculcated on the fly.
get_breaks <- function(x) {
    stopifnot("'x' must be a transitreg model" = inherits(x, "transitreg"))
    if (!is.null(x[["breaks"]])) {
        res <- x$breaks
    } else {
        res <- seq.int(0, x$bins) - 0.5
    }
    return(res)
}

# Helper function to get bin mids of a transitreg model
get_mids <- function(x) {
    stopifnot("'x' must be a transitreg model" = inherits(x, "transitreg"))
    res <- get_breaks(x)
    res <- (res[-1L] + res[-length(res)]) / 2.0
    return(res)
}

## TODO(R) DELETE ME ## prediction_get_new_response <- function(newdata, y, response) {
## TODO(R) DELETE ME ##     ## Newdata provided (by user), y is null: Response must be in 'newdata'.
## TODO(R) DELETE ME ##     if (!is.null(newdata) && is.null(y)) {
## TODO(R) DELETE ME ##         if (!response %in% names(newdata))
## TODO(R) DELETE ME ##             stop("response \"", response, "\" not found in 'newdata'. ",
## TODO(R) DELETE ME ##                  "Must be in 'newdata' or provided via the extra 'y' argument.")
## TODO(R) DELETE ME ##         res <- newdata[[response]]
## TODO(R) DELETE ME ##     ## If both 'newdata' and 'y' are set, force to use 'y' but
## TODO(R) DELETE ME ##     ## throw a warning.
## TODO(R) DELETE ME ##     } else if (!is.null(newdata) && !is.null(y)) {
## TODO(R) DELETE ME ##         if (response %in% names(newdata)) {
## TODO(R) DELETE ME ##             warning("Response \"", response, "\" provided via 'newdata' as well
## TODO(R) DELETE ME ##                     as ", "via argument 'y'. 'y' will overwrite 'newdata$",
## TODO(R) DELETE ME ##                     object$response, "' (w/ recycling).")
## TODO(R) DELETE ME ##             res <- rep(y, length.out = nrow(newdata))
## TODO(R) DELETE ME ##         } else {
## TODO(R) DELETE ME ##             res <- rep(y, length.out = nrow(newdata))
## TODO(R) DELETE ME ##         }
## TODO(R) DELETE ME ##     }
## TODO(R) DELETE ME ##     return(res)
## TODO(R) DELETE ME ## }
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ## # Helper function setting up 'newdata' when using predict.transitreg
## TODO(R) DELETE ME ## get_newdata <- function(object, newdata, y, prob, type) {
## TODO(R) DELETE ME ##     # Keep a logical flag for whether or not the user provided 'newdata'
## TODO(R) DELETE ME ##     newdata_provided = !is.null(newdata)
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##     # Newdata not provided? Taking model.frame (outer model; not binary
## TODO(R) DELETE ME ##     # response model, see 'engine' in transitreg().
## TODO(R) DELETE ME ##     if (is.null(newdata)) newdata <- model.frame(object)
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##     ## If 'newdata' is provided by user, apply scaling if required.
## TODO(R) DELETE ME ##     ## Scales all covariates except the 'theta' variables.
## TODO(R) DELETE ME ##     if (newdata_provided && !is.null(object$scaler)) {
## TODO(R) DELETE ME ##         for (j in names(object$scaler)) {
## TODO(R) DELETE ME ##             if (j != "theta") {
## TODO(R) DELETE ME ##                 newdata[[j]] <- (newdata[[j]] - object$scaler[[j]]$mean) / object$scaler[[j]]$sd
## TODO(R) DELETE ME ##             }
## TODO(R) DELETE ME ##         }
## TODO(R) DELETE ME ##         ## Scaling (standardizing) theta if requested
## TODO(R) DELETE ME ##         if (!is.null(theta_scaler))
## TODO(R) DELETE ME ##           newdata$theta <- (newdata$theta - theta_scaler$mean) / theta_scaler$sd
## TODO(R) DELETE ME ##     }
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##     ## Depending on 'type' we need to ensure y/prob is set correctly.
## TODO(R) DELETE ME ##     ## Quantile:
## TODO(R) DELETE ME ##     ## - prob must be not NULL, all values must be in [0, 1].
## TODO(R) DELETE ME ##     ## PDF/CDF:
## TODO(R) DELETE ME ##     ## - If 'newdata = NULL' we take the model frame, so y can be missing.
## TODO(R) DELETE ME ##     ## - Else y must be not NULL, all values must be within support of 'breaks'.
## TODO(R) DELETE ME ##     if (type == "quantile") {
## TODO(R) DELETE ME ##       stopifnot(
## TODO(R) DELETE ME ##         "for 'type = \"quantile\"' argument 'prob' must be set" = !is.null(prob),
## TODO(R) DELETE ME ##         "'prob' must be numeric of length > 0" = is.numeric(prob) && length(prob) > 0,
## TODO(R) DELETE ME ##         "missing values in 'prob' not allowed" = all(!is.na(prob)),
## TODO(R) DELETE ME ##         "'prob' must be in [0, 1]" = all(prob >= 0.0) && all(prob <= 1.0)
## TODO(R) DELETE ME ##       )
## TODO(R) DELETE ME ##     } else if (type %in% c("cdf", "pdf")) {
## TODO(R) DELETE ME ##       ## 'newdata' is provided by user (newdata_provided) we need to check
## TODO(R) DELETE ME ##       ## 1. If the response is not included in 'newdata', the user must provide the new
## TODO(R) DELETE ME ##       ##    response via 'y' (length 1 or same length as 'newdata').
## TODO(R) DELETE ME ##       ## 2. If the response is in both, 'newdata' and 'y' the response in newdata
## TODO(R) DELETE ME ##       ##    will be replaced with 'y' with an additional warning.
## TODO(R) DELETE ME ##       ## 3. No missing values are allowed in 'y'.
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##       ## (1)
## TODO(R) DELETE ME ##       if (newdata_provided & !object$response %in% names(newdata)) {
## TODO(R) DELETE ME ##         if (is.null(y)) {
## TODO(R) DELETE ME ##           stop("If the response \"", object$response, "\" is not included in ",
## TODO(R) DELETE ME ##                "'newdata' it must be provided via the argument 'y'.")
## TODO(R) DELETE ME ##         } else if (length(y) != 1L || length(y) != nrow(newdata)) {
## TODO(R) DELETE ME ##           stop("'y' (if provided) must be of length `1` or the same length as the number of ",
## TODO(R) DELETE ME ##                "rows in 'newdata' (provided by user or taken from internal model.frame) ",
## TODO(R) DELETE ME ##                "which is nrow(newdata) = ", nrow(newdata), ".")
## TODO(R) DELETE ME ##         }
## TODO(R) DELETE ME ##         # Store new response
## TODO(R) DELETE ME ##         newdata[[object$response]] <- num2bin(y, object$breaks)
## TODO(R) DELETE ME ##       ## (2) If response is in 'newdata' but also provided via 'y', it must still
## TODO(R) DELETE ME ##       ##     be of the correct length, and we will overwrite newdata[[object$response]]
## TODO(R) DELETE ME ##       ##     with 'y' showing a warning.
## TODO(R) DELETE ME ##       } else if (!is.null(y)) {
## TODO(R) DELETE ME ##         if (length(y) != 1L || length(y) != nrow(newdata)) {
## TODO(R) DELETE ME ##           stop("'y' (if provided) must be of length `1` or the same length as the number of ",
## TODO(R) DELETE ME ##                "rows in 'newdata' (provided by user or taken from internal model.frame) ",
## TODO(R) DELETE ME ##                "which is nrow(newdata) = ", nrow(newdata), ".")
## TODO(R) DELETE ME ##         }
## TODO(R) DELETE ME ##         ## Else store
## TODO(R) DELETE ME ##         warning("Response provided via 'newdata' as well as 'y'. ",
## TODO(R) DELETE ME ##                 "Response in 'newdata' will be overwritten by the values provided on 'y'.")
## TODO(R) DELETE ME ##         newdata[[object$response]] <- num2bin(y, object$breaks)
## TODO(R) DELETE ME ##       }
## TODO(R) DELETE ME ##       ## (3) Any missing values?
## TODO(R) DELETE ME ##       if (any(is.na(newdata[[object$response]]))) {
## TODO(R) DELETE ME ##         stop("Missing values in response (in 'newdata' or provided via argument 'y') ",
## TODO(R) DELETE ME ##              "are not allowed.")
## TODO(R) DELETE ME ##       }
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##       # Discrete distribution? y must be in range [0, object$bins].
## TODO(R) DELETE ME ##       if (is.null(object$breaks)) {
## TODO(R) DELETE ME ##         y <- as.integer(y)
## TODO(R) DELETE ME ##         if (!all(y >= 0 & y <= object$bins - 1L))
## TODO(R) DELETE ME ##             stop("elements in 'y' must be in {0, ..., ", object$bins - 1L, "}")
## TODO(R) DELETE ME ##       # Else continuous
## TODO(R) DELETE ME ##       } else {
## TODO(R) DELETE ME ##         if (!all(y >= min(object$breaks) & y <= max(object$breaks)))
## TODO(R) DELETE ME ##             stop(sprintf("elements in 'y' must be in range [%s, %s]",
## TODO(R) DELETE ME ##                          format(min(object$breaks), format(max(object$breaks)))))
## TODO(R) DELETE ME ##         ## If 'y' was set by the user, we must convert all values in 'y' from
## TODO(R) DELETE ME ##         ## it's original numeric value to it's pseudo-observation (bin; int).
## TODO(R) DELETE ME ##         if (!is.null(y)) y <- num2bin(y, object$breaks)
## TODO(R) DELETE ME ##       }
## TODO(R) DELETE ME ##     }
## TODO(R) DELETE ME ## 
## TODO(R) DELETE ME ##     # Returning prepared object
## TODO(R) DELETE ME ##     return(newdata)
## TODO(R) DELETE ME ## }


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
#' @seealso [transitreg()], [transitreg_tmf()].
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
    return(Transition(tp, breaks, censored = x$censored))
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
#' @seealso [transitreg()], [transitreg_dist()], [transitreg_tmf()], [predict.transitreg()].
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
#' @importFrom topmodels qqrplot wormplot rootogram
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

  if ("hist-resid" %in% which) rootogram(x, ...)
  if ("qq-resid" %in% which) qqrplot(x, ...)
  if ("wp-resid" %in% which) wormplot(b, ...)
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

    if (is.null(newdata)) newdata <- model.frame(object)
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
#' @seealso [transitreg()], [transitreg_tmf()], [transitreg_dist()].
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

#' @param object Transition model, object of class `transitreg`.
#' @param newdata An optional data frame in which to look for variables with
#'        which to predict. If omitted, the fitted values are used.
#'
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


