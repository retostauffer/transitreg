
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
#' @param theta_vars `NULL`  or character vector with the bin-specific theta variables
#'        to be set up (e.g., `c("theta0", "theta99")`). Can be an empty character
#'        vector (or `NULL`) if there are no `theta_vars`.
#' @param newresponse `NULL` or integer. New response vector to overwrite the one
#'        in `data`.
#' @param scaler Can be `FALSE`/`NULL`, `TRUE`, or a list. See section 'Scaler' for details.
#' @param verbose Logical value indicating whether information about the transformation
#'        process should be printed to the console. Default is `FALSE`.
#'
#' @details
#' Transition models focus on modeling the conditional probabilities of transitions
#' between states or counts. This function converts the input data into a long format
#' suitable for such models. Each row in the resulting data frame corresponds to a
#' binary transition indicator, representing whether a transition to a higher category
#' occurred. For details on the modeling framework, see Berger and Tutz (2021).
#'
#' @section 'Scaler':
#' When estimating new model ([transitreg()]), the user can specify `scale.x = TRUE`
#' which will standardize all covariates as well as the `theta` variable
#' using `(z - mean(z)) / sd(z)` with `z` being a covariate or `theta`. If the
#' input `scaler` to this function is set `TRUE` the return of this function
#' will provide an attribute `"scaler"` which contains a list, storing mean and
#' standard deviation used for scaling such that it can be applied to new data
#' (e.g., during prediction). If `FALSE` or `NULL` no scaling is applied, and no
#' attribute is set.
#'
#' If input `scaler` is a list, the scaling is applied given the content of that list.
#' No additional attribute will be set on the return value.
#'
#' @return
#' A transformed data frame in the long format. Each row represents a binary transition
#' indicator (`Y`) for the response variable. Additional columns in the output include:
#'
#' * `index`: The original row index from the input data.
#' * `Y`: The binary indicator for whether a transition to a higher
#'       category occurred.
#' * `theta`: The level corresponding to the current transition.
#' * `theta[0-9]+`: Special binary variables which are `1L` for all rows where
#'   (the unscaled) `theta` matches the integer `[0-9]+$`, else `0L`. Indicator
#'   that we are in bin `[0-9]+$`. Only set if `theta_vars` is specified accordingly.
#'
#' The return will contain the response name as attribute `"response"` and (potentially)
#' a `"scaler"` attribute (see section 'Scaler').
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
transitreg_tmf <- function(data, response, breaks, theta_vars = NULL,
                            newresponse = NULL, scaler = NULL, verbose = FALSE, ...) {

  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data) && all(dim(data) > 0L),
    "'response' must be NULL or a character of length 1" =
        is.null(response) || (is.character(response) && length(response) == 1L),
    "'newresponse' must be NULL or a numeric vector with length > 0" =
        is.null(newresponse) || (is.numeric(newresponse) && length(newresponse > 0)),
    "'theta_vars' must be NULL or character vector" =
        is.null(theta_vars) || is.character(theta_vars),
    "'breaks' must be numeric vector of length > 1L" =
        is.numeric(breaks) && length(breaks) > 1L,
    "'newresponse' must be NULL or bins (0, 1, 2, ...)" =
        is.null(newresponse) || (is.integer(newresponse) && length(newresponse) > 0L),
    "'verbose' must be TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose),
    "'scaler' must be either FALSE, TRUE, or a list" =
        isTRUE(scaler) || (isFALSE(scaler) || is.null(scaler)) || is.list(scaler)
  )


  ## Determine response column if not specified by user.
  if (is.null(response))
    response <- names(data)[1L]


  # -----------------------------------------------------------------
  # Converting response from original scale to (pseudo-)index
  # -----------------------------------------------------------------

  ## Discretize numeric response into counts.
  yc <- num2bin(data[[response]], breaks)
  ym <- (breaks[-1] + breaks[-length(breaks)]) / 2

  lower <- list(...)$lower
  upper <- list(...)$upper

  if (!is.null(lower)) {
      stopifnot("'lower' must be single numeric" = is.numeric(lower) && length(lower) == 1L)
      ym[ym < lower] <- lower
  }
  if (!is.null(upper)) {
      stopifnot("'upper' must be single numeric" = is.numeric(upper) && length(upper) == 1L)
      ym[ym > upper] <- upper
  }

  data[[response]] <- yc



  # -----------------------------------------------------------------
  # Preparing the rest
  # -----------------------------------------------------------------

  ## Using 'FALSE' instead of `NULL` for the rest of the function.
  if (is.null(scaler)) scaler <- FALSE

  ## If theta_vars is set: Expecting string of the form `theta0`, `theta1`, ...
  if (is.character(theta_vars) && length(theta_vars) == 0L)
      theta_vars <- NULL
  if (!is.null(theta_vars))
      stopifnot("unexpected entries in 'theta_vars'" = all(grepl("^theta[0-9]+$", theta_vars)))


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

  ## data[[response]] must all be bin indices, so integers >= -1
  check <- all(data[[response]] > -(1 + sqrt(.Machine$double.eps)) |
               abs(data[[response]] %% 1) > sqrt(.Machine$double.eps))
  if (!check)
    stop("The response must be bin indices, so integers in the range of {-1L, Inf}.")
  data[[response]] <- as.integer(data[[response]])

  ## WARNING: Why do we allow for -1L in data[[response]]?
  ## The pseudo index "-1" is for values falling below the lowest bound;
  ## required to allow predictions outside the range. To properly handle some of the
  ## methods we create a "pseudo_index" which is pmax(0, data[[response]]).
  pseudo_index <- pmax(0L, data[[response]])

  ## Setting up the new data.frame with (pseudo-)bins

  ## Define length of vectors in list Sum of response pseudo indices +
  ## nrow(data), the latter to account for the additional "0" bin.
  nout <- sum(pseudo_index) + nrow(data)

  ## ------ building transitreg data -------

  ## Names of list elements.
  names_out <- c("index", "Y", "theta", names(data))

  ## Named of all elements which MUST NOT BE standardized, even if a scaler is
  ## provided or scaler is set TRUE.
  names_dont_scale <- c("index", "Y", response, theta_vars)

  ## Building index vector; each observation 1:nrow(data) gets its unique index
  result <- list()
  result$index <- rep(seq_len(nrow(data)), times = pseudo_index + 1L)

  ## Creating Y; always 1 except for the last entry per index.
  fn_get_Y <- function(nout, resp) {
    res <- rep(1L, nout); res[cumsum(resp + 1)] <- 0L; return(res)
  }
  result$Y <- fn_get_Y(nout, pseudo_index)

  ## Creating theta; a sequence from zero to the response_value for each index.
  ## The following Two lines create this sequence of sequences.
  fn_get_theta <- function(nout, resp, idx) {
    resettozero <- c(0L, which(diff(idx) > 0))
    return(seq_len(nout) - rep(resettozero, resp + 1L) - 1L)
  }
  result$theta <- fn_get_theta(nout, pseudo_index, result$index)

  ## Adding theta_vars if needed.
  ## For 'theta99' in 'theta_vars' a new variable 'theta99' is generated which
  ## is 1L where 'theta == 99', else 0L. Same for all other required theta_vars.
  if (!is.null(theta_vars)) {
      theta_vars <- theta_vars[order(as.integer(regmatches(theta_vars, regexpr("[0-9]+$", theta_vars))))]
      for (tv in theta_vars) {
          tint <- as.integer(regmatches(tv, regexpr("[0-9]+$", tv)))
          ## If this theta does not exist, skip (do not add a new variable
          ## with constant zero values).
          if (sum(result$theta == tint) == 0) break
          ## Else setting up the new binary variable
          result[[tv]] <- integer(length(result$theta)) # Initialize 0s
          result[[tv]][result$theta == tint] <- 1L
      }

  }

  ## Appending the remaining data from 'data'.
  for (n in names(data)) {
      ## If data[[n]] is a simple vector
      if (!is.matrix(data[[n]])) {
        result[[n]] <- rep(data[[n]], pseudo_index + 1)
      ## Else create matrix
      } else {
        result[[n]] <- matrix(rep(data[[n]], rep(data[[response]] + 1, ncol(data[[n]]))),
                              ncol = ncol(data[[n]]),
                              dimnames = list(NULL, colnames(data[[n]])))
      }
  }

  result <- as.data.frame(result)

  ## scaler == TRUE (comes from `scale.x = TRUE` when calling `transitreg()`; initial
  ## standartizaion of all covariates as well as theta. Perform scaling and keep the
  ## used first and second order momentum. Will be added as an attribute to the return.
  if (isTRUE(scaler)) {
      scaler <- list()
      for (j in names(result)[!names(result) %in% names_dont_scale]) {
        if (!is.factor(result[[j]])) {
            scaler[[j]] <- list("mean" = mean(result[[j]]), "sd" = sd(result[[j]]))
            result[[j]] <- (result[[j]] - scaler[[j]]$mean) / scaler[[j]]$sd
        }
      }
  ## Existing 'scaler' is provided, typically used when performing predictions and
  ## new data ('newdata') are provided by the user on the original (unscaled) scale.
  ## Will apply the ame scaling as on the initial model estimation.
  } else if (is.list(scaler)) {
      for (j in names(scaler)) {
          result[[j]] <- (result[[j]] - scaler[[j]]$mean) / scaler[[j]]$sd
      }
  } else {
      scaler <- NULL
  }

  ## Attach the response column as an attribute.
  attr(result, "response") <- response
  attr(result, "scaler")   <- scaler

  return(result)
}

