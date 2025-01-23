

#' Building Model Response for Transition Model
#'
#' To estimate a model as well as to perform predictions the model frame needs
#' to be built. This function (for internal use only) sets up the model
#' response properly, scaling and discretizising the response (if requested).
#'
#' @param x model framme, a `data.frame` where the first variable is expected
#'        to be the model response.
#' @param response Character of length one, name of the response variable.
#' @param breaks numeric vector with breaks for discretizing the original response.
#'        Can be `NULL` if no discretization is required.
#'
#' @return Named list with the following three elements:
#'
#' * `mf`: Modified version of input `x`; altered response with (pseudo-)bins
#'   if required.
#' * `breaks`: Numeric vector with breaks, point intersection between (pseudo-)bins.
#' * `bins`: Number of bins, same as `length(breaks) - 1L`.
#' * `ym`: mid-point of the bins.
#' * `yc`: index vector, (pseudo-)bin index/number.
transitreg_response <- function(x, response, breaks, ...) {
    stopifnot(
        "'x' must be a data.frame with positive dimensions" =
            is.data.frame(x) && length(dim(x)) == 2L && all(dim(x) > 0),
        "all variables in 'x' must be numeric" = all(sapply(x, is.numeric)),
        "'response' must be character of length 1" = is.character(response) && length(response) == 1L,
        "'response' variable not found in 'x'!" = response %in% names(x),
        "'breaks' must be numeric or NULL" = is.numeric(breaks) || is.null(breaks)
    )

    # Store the name of the response for later
    response <- names(x)[1L]

    ## Discretize response?
    if (!is.null(breaks)) {
        ## Create breaks automatically (if needed)
        if (length(breaks) == 1L)
            breaks <- make_breaks(model.response(x), breaks = breaks)
        ## Number of bins (one less than 'breaks'; points of bin intersections)
        bins <- length(breaks) - 1L

        ## Discretize numeric response into counts.
        yc <- num2bin(x[[response]], breaks)
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

        x[[response]] <- yc

        result <- list(mf = x, breaks = breaks, bins = bins, ym = ym, yc = yc)
    } else {
        ## For the model we need an integer response; check if the response
        ## is integer. If not, break.
        if (!all(x[[response]] %% 1 < sqrt(.Machine$double.eps)) &&
             all(x[[response]]      > sqrt(.Machine$double.eps)))
            stop("Response is not looking like count data (integers); binning via 'breaks' is required.")
        x[[response]] <- as.integer(round(x[[response]], 1))

        bins <- max(x[[response]])
        ## Setting highest 'bin' to max(response) * mp (multiplier)
        mp     <- if (bins <= 10) { 3 } else if (bins <= 100) { 1.5 } else { 1.25 }
        bins   <- as.integer(ceiling(bins * mp))
        breaks <- seq_len(bins + 1) - 1.5 # 'Integer' bins
        if (verbose) message("Response considered to be count data, using max count ", bins - 1)

        result <- list(mf = x, breaks = NULL, bins = bins, ym = NULL, yc = NULL)
    }

    return(result)
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
#' @param theta_vars `NULL`  or character vector with the bin-specific theta variables
#'        to be set up (e.g., `c("theta0", "theta99")`). Can be an empty character
#'        vector (or `NULL`) if there are no `theta_vars`.
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
transitreg_data <- function(data, response, theta_vars = NULL,
                            newresponse = NULL, verbose = TRUE) {

  stopifnot(
    "'response' must be NULL or a character of length 1" =
        is.null(response) || (is.character(response) && length(response) == 1L),
    "'newresponse' must be NULL or a numeric vector with length > 0" =
        is.null(newresponse) || (is.numeric(newresponse) && length(newresponse > 0)),
    "'theta_vars' must be NULL or character vector" =
        is.null(theta_vars) || is.character(theta_vars),
    "'newresponse' must be NULL or bins (0, 1, 2, ...)" =
        is.null(newresponse) || (is.integer(newresponse) && all(newresponse >= 0)),
    "'verbose' must be TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose)
  )

  ## Ensure data is a data frame.
  if (!is.data.frame(data))
    data <- as.data.frame(data)

  ## Determine response column if not specified by user.
  if (is.null(response))
    response <- names(data)[1L]

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

  ## Creating theta; a sequence from zero to the response_value for each index.
  ## The following Two lines create this sequence of sequences.
  fn_get_theta <- function(nout, resp, idx) {
    resettozero <- c(0, which(diff(idx) > 0))
    return(seq_len(nout) - rep(resettozero, resp + 1) - 1)
  }
  result$theta <- fn_get_theta(nout, data[[response]], result$index)

  ## Adding theta_vars if needed.
  ## For 'theta99' in 'theta_vars' a new variable 'theta99' is generated which
  ## is 1L where 'theta == 99', else 0L. Same for all other required theta_vars.
  if (!is.null(theta_vars)) {
      theta_vars <- theta_vars[order(as.integer(regmatches(theta_vars, regexpr("[0-9]+$", theta_vars))))]
      for (tv in theta_vars) {
          tint <- as.integer(regmatches(tv, regexpr("[0-9]+$", tv)))
          result[[tv]] <- integer(length(result$theta)) # Initialize 0s
          result[[tv]][result$theta == tint] <- 1L
      }
  }

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


