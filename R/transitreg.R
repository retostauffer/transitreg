
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
#' @seealso [transitreg_tmf()], [mgcv::gam()]
#'
#' @examples
#' ## Example 1: Count data.
#' library("transitreg")
#' set.seed(123)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x))) + 12
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
#' ## Plotting hanging rootogram, quantile residuals, and
#' ## a probability integral transform (PIT) histogram.
#' plot(b, which = c("rootogram", "qqrplot", "pithist"))
#'
#' ## Predictions and plotting.
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#' fit <- cbind(
#'   "97.5%"  = predict(b, nd, type = "quantile", p = 1 - 0.05/2),
#'   "median" = predict(b, nd, type = "quantile", p = 0.5),
#'   "2.5%"   = predict(b, nd, type = "quantile", p = 0.05/2),
#'   "mode"   = predict(b, nd, type = "mode")
#' )
#'
#' ## Plot data and fitted counts.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", add = TRUE, lwd = 2,
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#' legend("topleft", legend = colnames(fit),
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#'
#' ## Visualizing three predicted distributions via Transitreg distributions
#' idx <- c(25, 50, 75)
#' d3 <- Transition(predict(b, nd[idx, , drop = FALSE], type = "tp"), breaks = seq.int(0, b$bins) - 0.5)
#' plot(d3, type = "pdf")
#' abline(v = predict(b, nd[idx, , drop = FALSE], type = "mode"), col = 1:3, lty = 3)
#'
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
#'   "97.5%"  = predict(b, nd, type = "quantile", p = 1 - 0.05/2),
#'   "median" = predict(b, nd, type = "quantile", p = 0.5),
#'   "2.5%"   = predict(b, nd, type = "quantile", p = 0.05/2),
#'   "mode"   = predict(b, nd, type = "mode")
#' )
#'
#' ## Plot data and fitted curves.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", add = TRUE, lwd = 2,
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#' legend("topleft", legend = colnames(fit),
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#'
#' ## Visualizing 5 randomly selected fitted distributions
#' idx <- sample(seq_len(nrow(model.frame(b))), 5L)
#' plot(b[idx], type = "pdf")
#' plot(b[idx], type = "cdf")
#' plot(b[idx], type = "tp")
#'
#'
#' ## Example 3: Count response with neural network.
#' set.seed(123)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- rpois(n, exp(2 + sin(x)))
#'
#' ## Fit NN transition count response model.
#' b <- transitreg(y ~ theta + x, engine = "nnet",
#'                 size = 5, maxit = 1000, decay = 0.001)
#'
#' ## Predictions and plotting.
#' nd <- data.frame(x = seq(-3, 3, length = 100))
#' fit <- cbind(
#'   "97.5%"  = predict(b, nd, type = "quantile", p = 1 - 0.05/2),
#'   "median" = predict(b, nd, type = "quantile", p = 0.5),
#'   "2.5%"   = predict(b, nd, type = "quantile", p = 0.05/2),
#'   "mode"   = predict(b, nd, type = "mode")
#' )
#'
#' ## Plot data and fitted counts.
#' plot(y ~ x, pch = 16, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
#' matplot(nd$x, fit, type = "l", add = TRUE, lwd = 2,
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#' legend("topleft", legend = colnames(fit),
#'        col = c(4, 4, 4, 2), lty = c(2, 1, 2, 1))
#'
#' @keywords models regression
#'
#' @importFrom stats as.formula binomial predict sd update
#' @importFrom mgcv gam bam
#'
#' @author Niki
#' @export
transitreg <- function(formula, data, subset, na.action,
                       engine = "bam", breaks = NULL,
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
    "'engine' must be character of length 1" = is.character(engine) && length(engine) == 1L
  )
  ncores <- transitreg_get_number_of_cores(ncores, verbose = verbose)

  ## Evaluate 'engine' argument
  engine <- tolower(engine)
  engine <- match.arg(engine, c("bam", "gam", "nnet", "glmnet"))

  ## Check if `scale.x` was specified via the `...` argument (hidden feature).
  ## If so, it must be TRUE or FALSE. Else it is set TRUE if `engine = "nnet"`
  ## and `FALSE` else.
  scale.x <- list(...)$scale.x
  if (!is.null(scale.x))
      stopifnot("'scale.x' (if specified via '...') must be TRUE or FALSE" =
                isTRUE(scale.x) || isFALSE(scale.x))
  if (is.null(scale.x)) scale.x <- engine == "nnet" ## defaults to TRUE if engine = "nnet"

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
  theta_vars <- theta_vars[grep("^theta[0-9]+$", theta_vars)]

  if (length(theta_vars) > 0L)
    ff2 <- eval(parse(text = paste0("update(ff2, . ~ . -",
                paste(theta_vars, collapse = " - "), ")")))

  mf[["formula"]] <- ff2
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

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
      if (any(mf[[rval$response]] < 0) ||
          any(abs(mf[[rval$response]] %% 1) > sqrt(.Machine$double.eps))) {
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

  ## Testing theta_vars. Will fail if:
  ## - We have no observations falling into bin X (thetaX)
  ## - There are dedicated thetaX, thetaY, thetaZ, ..., covering all bins
  ##   populated with observations (overspecified model).
  if (length(theta_vars) > 0L) {
    theta_int  <- as.integer(regmatches(theta_vars, regexpr("[0-9]+$", theta_vars)))
    # Check if formula contains thetaX but no observations fall into bin X
    tmp <- theta_int[!theta_int %in% unique(y)]
    if (length(tmp) > 0L)
        stop("Formula contains ", paste(sprintf("'theta%d'", tmp), collapse = ", "), ", ",
             "but no observation falls into this bin (misspecified model formula).")
    # Test if there is a thetaX for each bin we have observation (overspecified)
    # throw an error as well.
    if (all(!is.na(match(y, theta_int))))
        stop("Formula contains ", paste(sprintf("'theta%d'", theta_int), collapse = ", "), ", ",
             "covering all bins observations fall into (overspecified).")
  }

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

  ## Compute probabilities.
  ui <- unique(tmf$index)
  probs <- cprobs <- numeric(length(ui))

  ## c_transitreg_predict_pdfcdf returns a list with PDF and CDF, calculating
  ## both simultanously in C to improve speed.
  censored <- if (is.null(rval$censored)) "not-censored" else rval$censored
  args <- list(uidx = ui, idx = tmf$index,
               tp = tp, y = y, breaks = breaks,
               discrete = rep(is.null(rval$breaks), length(ui)),
               ncores = ncores, censored = censored)

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
        type = c("pdf", "cdf", "quantile", "mode", "tp"), y = NULL, prob = NULL,
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

  ## Setting number of cores for OMP
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
  ## Else we check if we have the variables needed, and extract them from
  ## the newdata object provided by the user.
  if (is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    tmp <- names(model.frame(object))
    tmp <- tmp[!tmp == object$response]
    if (!all(tmp[!tmp == object$response] %in% names(newdata)))
        stop("'newdata' does not provide all required variables, ",
             "missing: ", paste(tmp[!tmp %in% names(newdata)], collapse = ", "))
    mf <- newdata
    rm(tmp)
  }

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
      } else if (type == "mode") {
          elementwise <- TRUE # for 'mode' elementwise is always TRUE
      } else if (type == "tp") {
          NULL
      } else {
          stop("Mode for type = \"", type, "\" must be implemented!")
      }
  }

  ## Setting up 'newresponse'.
  ## Quantile, mode, tp:
  ##  - Setting response to "max bin mid" to ensure we calculate the
  ##    transition probabilities for _all_ bins.
  ## CDF/PDF:
  ##  - If elementwise = TRUE: We only need to evaluate each distribution
  ##    up to 'y[i]'.
  ##  - If elementwise = FALSE: We must evaluate each distribution up to
  ##    max(y).
  if (type %in% c("quantile", "mode", "tp")) {
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

  ## Get rows (row index) where we have missing data
  obs_na <- unname(apply(mf, MARGIN = 1, function(x) sum(is.na(x))) > 0)
  if (all(obs_na)) {
    stop("all observations (rows) contain missing data, prediction not possible")
    # TODO(R): Create tests for this
  }

  ## Specify argument for generic prediction method called below
  what <- switch(class(object$model)[1L],
    "bam"  = "response",
    "gam"  = "response",
    "nnet" = "raw"
  )


  ## ------------------------------------------------
  ## Calculating how long the tmf matrix will be.
  ## tmf_rc is the 'tmf data.frame row count' we expect.
  tmf_rc <- integer(nrow(mf))
  tmf_rc[!obs_na] <- num2bin(mf[!obs_na, object$response], get_breaks(object))
  tmf_rc <- cumsum(tmf_rc) # Cumulative sum

  tmf_maxrows <- 1e7
  blockindex  <- tmf_rc %/% tmf_maxrows + 1L

  tp <- list()
  for (block in seq_len(max(blockindex))) {
    idx <- which(blockindex == block)
    if (all(obs_na[idx])) next

    ## Creating 'transition model frame' for the prediction of the
    ## transition probabilities using the object$model (binary response model).
    tmf <- transitreg_tmf(mf[blockindex == block & !obs_na, , drop = FALSE],
                          response   = object$response,
                          breaks     = breaks,
                          theta_vars = object$theta_vars,
                          scaler     = object$scaler, verbose = verbose)

    tp[[block]] <- as.numeric(predict(object$model, newdata = tmf, type = what))
  }
  tp <- do.call(c, tp)


  ## ------------------------------------------------
  ## If 'type = "tp"' (transition probabilities) we already have our
  ## result. Convert to matrix and return.
  if (type == "tp") {
      tmp <- list(NULL, paste0("tp_", seq_len(object$bins) - 1))
      arr.ind <- cbind(row = rep(which(!obs_na), each = object$bins),
                       col = rep(seq_len(object$bins), times = sum(!obs_na)))
      # Initialize empty matrix, fill with tps
      x <- matrix(NA, ncol = object$bins, nrow = nrow(mf), dimnames = tmp)
      x[arr.ind] <- tp
      return(x)
  }

  ## Extract unique indices
  ui   <- unique(tmf$index)

  ## Setting dummy values (required by C later on)
  if (type == "quantile") {
      ## Sorting 'prob'. This is important for the .C routine!
      probC <- if (!elementwise) sort(unique(prob)) else prob[!obs_na]
      yC   <- NA_integer_ ## Dummy value required for .C call
  } else if (type %in% c("cdf", "pdf")) {
      ## Sorting 'y'. This is important for the .C routine!
      if (elementwise) {
          yC <- num2bin(mf[!obs_na, object$response], breaks = breaks)
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
  args <- list(uidx        = ui,           # int; Unique distribution index (int)
               idx         = tmf$index,    # int; Index vector (int)
               tp          = tp,           # num; Transition probabilities
               breaks      = breaks,       # num; Point intersections of bins
               y           = yC,           # int; Response y, used for 'cdf/pdf'
               prob        = probC,        # num; Probabilities (used for 'quantile')
               type        = type,         # str; to predict/calculate
               ncores      = ncores,       # int; Number of cores to be used (OpenMP)
               elementwise = elementwise,  # Elementwise (one prob or y per ui)
               discrete    = discrete,     # Discrete distribution?
               censored    = censored)

  # Calling C
  args <- check_args_for_treg_predict(args)
  res  <- do.call(function(...) .Call("treg_predict", ...), args)

  # If 'ementwise = FALSE' we get length(y)/length(prob) results per
  # observation and have to glue them back together into a matrix.
  if (!elementwise && type == "quantile") {
    arr.ind <- cbind(row = rep(which(!obs_na), each = length(probC)),
                     col = rep(seq_along(probC), times = sum(!obs_na)))
    x <- matrix(NA, nrow = nrow(mf), ncol = length(prob),
                  dimnames = list(NULL, get_elementwise_colnames(prob, NULL)))
    x[arr.ind] <- res
  } else if (!elementwise) {
    ## Convert from predicted 'bin' back to the numeric value.
    prefix <- if (type == "pdf") "d" else "p"
    ## Setting up empty matrix and populate with results where the
    ## observations did not contain missing data.
    arr.ind <- cbind(row = rep(which(!obs_na), each = length(y)),
                     col = rep(seq_along(y), times = sum(!obs_na)))
    x <- matrix(NA, nrow = nrow(mf), ncol = length(y),
                dimnames = list(NULL, get_elementwise_colnames(y, prefix)))
    print(mf)
    x[arr.ind] <- res
  ## Elementwise - building vector return
  } else {
    x <- rep(NA_real_, nrow(mf))
    x[!obs_na] <- res
  }

  return(x)
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
#' @param which A character or integer, specifying the type of plot(s) to
#'        generate (See 'Details').
#' @param ask Either `NULL`, `TRUE` or `FALSE`. If `NULL` it will evaluate to
#'        `TRUE` if more than one plot is requested via `which`.
#' @param \dots Additional arguments passed to the underlying plotting functions.
#'
#' @details
#' The function allows to control what to plot via the argument `which`.
#' Options include:
#'
#' * `"effects"` Effects of the predictors on the response.
#'   Requires that the model is estimated by [mgcv::gam()] or [mgcv::bam()].
#' * `"rootogram"`: Rootogram for assessing goodness of fit.
#' * `"qqrplot"`: Q-Q plot for quantile residuals.
#' * `"wormplot"`: Worm plot for quantile residuals.
#' * `"pithist"`: Probability integral transform (PIT) histogram.
#'
#' Multiple options can be specified as a character vector or numeric indices.
#'
#' The [plot.transitreg()] method provides flexible visualization options for
#' evaluating transition model fits. Users can choose to:
#'
#' * Visualize the effects of predictors on the response variable
#'   (if the model is a GAM, see [mgcv::gam()]).
#' * Evaluate quantile residuals through histograms, Q-Q plots, or worm plots
#'   employing the `topmodels` framework for model evaluation.
#'
#' The `which` argument controls the type of plots generated. By default, the
#' `"effects"` plot is shown if the model supports it. Residual-based plots
#' (`"rootogram"`, `"qqrplot"`, `"wormplot"`, `"pithist"`) provide insights into
#' the goodness of fit of the model and model calibration.
#'
#' @return Returns `NULL` invisibly. Generates plots as a side effect.
#'
#' @seealso [transitreg()], [predict.transitreg()].
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
#' @importFrom topmodels qqrplot wormplot rootogram pithist
#' @importFrom stats residuals
#' @importFrom graphics par
plot.transitreg <- function(x, which = "effects", ask = NULL, ...) {
  ## Staying sane
  stopifnot(
    "'effects' must be numeric or character" = is.character(which) || is.numeric(which),
    "'ask' must be NULL, TRUE, or FALSE" = is.null(ask) || isFALSE(ask) || isTRUE(ask)
  )

  # Helper function
  effects_gam <- function(x, ...) plot(x$model, ...)

  ## Defines the plotting functions, as well as what is available
  avail <- list(effects   = effects_gam,
                rootogram = rootogram,
                qqrplot   = qqrplot,
                pithist   = pithist,
                wormplot  = wormplot)

  ## Evaluate 'which'; can be numeric or string
  if (is.numeric(which)) {
      which <- as.integer(which)
      which <- names(avail)[which[which >= 1L & which <= length(avail)]]
  }
  which <- match.arg(unique(which), names(avail), several.ok = TRUE)

  ## 'effects' only for "gam" models
  if (!inherits(x$model, "gam")) {
    warning("'effects' plot only available for 'gam'-based transitreg models")
    which <- which[!which == "effects"]
    if (length(which) == 0L) {
      warning("nothing to plot, returning invisible NULL")
      invisible(NULL)
    }
  }

  if (length(which) > 1L && ((isTRUE(ask)) || is.null(ask))) {
    on.exit(par(ask = FALSE)); par(ask = TRUE)
  }

  # Plotting user requests in the order requested
  for (w in which) avail[[w]](x, ...)

  invisible(NULL)
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
        stop("Response variable '", object$response, "' missing in newdata!")

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
#' * `"mode"`: The expected value of the response (maximum probability).
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
#' * For `"mode"`, the expected value of the response.
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
#' ## Predict mode (expected value at highest probability)
#' p$mode <- predict(b, type = "mode")
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
        type = c("pdf", "cdf", "quantile", "mode", "tp"), ncores = NULL,
        elementwise = NULL, verbose = FALSE, ...) {

  type <- tolower(type)
  type <- match.arg(type)

  # TODO(R): Write test for this
  if (!is.null(prob))
    type <- "quantile"

  # TODO(R): Write test for this
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


