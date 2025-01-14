
#' Creates a Transition Distribution
#'
#' A 'Transition' distrubiton consists of a series of \code{K} transition
#' probabilities (TP) for \code{K} 'bins' (counts or pseudo-counts).
#'
#' @param x numeric vector or a numeric matrix.
#' @param bins numeric vector of points of intersection of the bins.
#'        The length the vector must be of \code{length(x) + 1} (if \code{x} is
#'        a vector) or \code{ncol(x) + 1} if \code{x} is a matrix. Must be
#'        monotonically increasing.
#'
#' @return Returns an object of class \code{c("Transition", "distribution")}.
#'
#' @author Reto
Transition <- function(x, bins) {

    # Sanity checks
    stopifnot(
        "'x' must be numeric (vector or matrix)" = is.numeric(x) && is.atomic(x),
        "'x' must be a vector or a matrix" = is.vector(x) || is.matrix(x),
        "length of 'x' must be > 0" = length(x) > 0L,
        "'bins' must be a numeric vector" = is.atomic(bins) && is.numeric(bins),
        "missing values in 'bins' not allowed" = all(!is.na(bins)),
        "'bins' must be monotonically increasing" = all(diff(bins) > 0)
    )

    # If 'x' is a vector, convert to matrix
    if (is.vector(x)) x <- matrix(x, nrow = 1)

    # Checking 'bins' vector
    if (length(bins) != (ncol(x) + 1))
        stop("'bins' must be of length ", ncol(x) + 1)


    # Ensure to convert to double in case input is integer
    x[,] <- as.numeric(x)
    bins <- as.numeric(bins)

    res <- setNames(as.data.frame(x),
                    paste("tp", seq_len(ncol(x)) - 1, sep = "_"))

    structure(res, class = c("Transition", "distribution"), bins = bins)
}

# Combine Transition objects, only allowed if they have the very same bins
# (come from the same distribution; same number of transition probabilities).
c.Transition <- function(...) {
    x <- list(...)
    if (length(x) == 1) return(x[[1]])

    # Else check whether or not we can combine the objects
    for (i in seq.int(2, length(x))) {
        stopifnot("input not of class Transition" = inherits(x[[i]], "Transition"))
        if (!all.equal(attr(x[[1]], "bins"), attr(x[[2]], "bins")))
            stop("bins of the ", i, ifelse(i == 2, "nd", "th"),
                 "object not the same as for the first object. Can't be combined.")
    }

    # Combine and return
    res <- do.call(rbind, lapply(x, as.matrix))
    Transition(res, attr(x[[1]], "bins"))
}

prodist.transitreg <- function(object, newdata = NULL, ...) {
    # Extracting covariable names to create newdata
    covars <- attr(terms(fake_formula(formula(object))), "term.labels")
    covars <- covars[!covars == "theta"]

    # In-sample data
    if (is.null(newdata)) {
        res <- object$model.frame
    } else {
        res <- newdata; rm(newdata)
    }

    # Creating res
    nb <- length(object$ym) # Number of bins
    nd <- nrow(res)         # Number of observations

    expand_covar <- function(x, nb) rep(x, each = nb)
    res  <- lapply(res[, covars, drop = FALSE], expand_covar, nb = nb)
    res  <- data.frame(c(list(theta = rep(seq_len(nb) - 1, times = nd)), res))

    # TODO(R): Currently 'type = response' which differs
    # for different engines (see transitreg()).
    res  <- data.frame(tp = predict(object$model, newdata = res, type = "response"),
                       lo = rep(head(object$bins, -1), times = nd),
                       up = rep(tail(object$bins, -1), times = nd))
    # Split into individual data.frames
    res <- split(res, rep(seq_len(nd), each = nb))

    # Convert into distributions object, pack into a data.frame
    return(setNames(as.data.frame(Transition(res)), as.character(substitute(object))))
}

# TODO(R): Not needed! Remove once d3 is implemented properly.
####procast.transitreg <- function(object, newdata = NULL, na.action = na.pass,
####                       type = "distribution", at = 0.5, drop = FALSE, ...) {
####    object <- if (is.null(newdata)) {
####        prodist(object)
####    } else {
####        # TODO(R): na.action passed to 'procast.transitreg' but not yet implemented (no effect)
####        prodist(object, newdata = newdata, na.action = na.action)
####    }
####
####    return(if (drop) object[[1]] else object)
####}


#' Convert Transition Distributions to Matrix
#'
#' @param x object of class \code{c("Transition", ...)}.
#' @param expand logical, if FALSE (default) the wide format is
#'        returned, else the extended (long) form.
#' @param \dots unused.
#'
#' @return Numeric matrix. If \code{expand = FALSE} the return is of dimension
#' \code{c(length(x), <ncol>)}.
#'
#' If \code{expand = TRUE} the returned matrix is of dimension
#' \code{c(length(x) * <ncol>, 3L)} where the four columns contain 
#' \code{index} (1, ..., length(x)) where each index corresponds to the
#' row-index of the original input \code{x} (distribution identifier),
#' \code{theta} which is the 'bin' the transition probability belongs to,
#' followed by the transition probabilities \code{tp} itself.
#' This expanded version is used when calling the .C functions.
#'
#' @author Reto
as.matrix.Transition <- function(x, expand = FALSE, ...) {
    stopifnot("'expand' must be logical TRUE or FALSE" = isTRUE(expand) || isFALSE(expand))

    xnames <- names(x) # Keep for later
    bins   <- attr(x, "bins")

    # convert to data.frame -> matrix
    x <- as.matrix(structure(x, class = "data.frame"))
    rownames(x) <- xnames
    if (expand) {
        # Keep original dimension as we transpose(x) in a second
        nd <- nrow(x)          # Number of distributions
        nb <- length(bins) - 1 # Number of bins

        index <- rep(seq_len(nd), each = nb) # distribution index
        theta <- rep(seq_len(nb) - 1, times = nd) # 'bin' index or theta
        x <- cbind(index = index, theta = theta, tp = as.vector(t(x)))
    }

    structure(x, class = c("Transitionmatrix", class(x)), bins = bins)
}


# Format
format.Transition <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and bins
    fmtfun <- function(i) {
        y <- as.matrix(x[i], expand = FALSE)
        sprintf("n = %d", ncol(y))
    }
    f <- sapply(seq_along(x), fmtfun)
    f <- sprintf("Transition(%s)", f)
    setNames(f, xnames)
}


pdf.Transition <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(x) == 1L) x <- rep(x, length(d))
        elementwise <- length(d) == length(x)
    }
    if (elementwise & length(d) != length(x))
        stop("'elementwise = TRUE' but length of x does not match the number of distributions")

    # Store element names for return
    xnames <- names(d)
    bins <- attr(d, "bins")

    # Convert numeric values to corresponding 'bin indices' (int)
    x <- num2bin(x, bins)

    if (!elementwise) x <- sort(x) # Important
    ui  <- seq_along(d) # Unique index
    idx <- rep(ui, each = length(bins) - 1) # Index of distribution

    # Setting up arguments to call .C predict function
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(d)),          # Transition probabilities
                 bins  = bins,                     # Point intersection of bins
                 y     = x,                        # Where to evaluate the pdf
                 prob  = NA_real_,                 # Dummy, only used for 'quantile'
                 type  = "pdf", ncores = ncores, elementwise = elementwise,
                 discrete = FALSE) # <- dummy value

    # Calling C
    check_args_for_treg_predict(args)
    res  <- do.call(function(...) .Call("treg_predict", ...), args)

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Returning simplified result
    } else if (length(d) == 1L & drop) {
        return(res)
    # Else return a matrix with dimension length(ui) x length(x)
    } else {
        # Create and return matrix
        return(matrix(res, byrow = TRUE, ncol = length(x),
                      dimnames = list(xnames, get_elementwise_colnames(x, "d"))))
    }
}

## Just returns log(pdf(...))
log_pdf.Transition <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    return(log(pdf(d, x, drop = drop, elementwise = elementwise, ncores = ncores, ...)))
}

#' Create Column Names for Return Matrix
#'
#' Used in the S3 methods pdf, cdf, and quantile. If `elementwise = FALSE`
#' the return is a matrix; this helper function creates the names based
#' on the thresholds/probabilities used.
#'
#' @param x numeric, thresholds (pdf, cdf) or probabilities (quantile).
#' @param prefix If \code{NULL} (quantiles) the result is in percent,
#'        else this prefix is used for each 'x'.
get_elementwise_colnames <- function(x, prefix = NULL, digits = 3) {
    if (is.null(prefix)) {
        x <- paste0(format(1e2 * x), "%")
    } else {
        x <- paste(prefix, trimws(format(x, digits = digits)), sep = "_")
    }
    if (anyDuplicated(x)) {
        for (d in unique(x[duplicated(x)])) {
            idx <- which(x == d)
            x[idx] <- paste0(x[idx], c("", paste0("_", seq_len(length(idx) - 1))))
        }
    }
    return(x)
}

cdf.Transition <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(x) == 1L) x <- rep(x, length(d))
        elementwise <- length(d) == length(x)
    }
    if (elementwise & length(d) != length(x))
        stop("'elementwise = TRUE' but length of x does not match the number of distributions")

    # Store element names for return
    xnames <- names(d)
    bins   <- attr(d, "bins")

    # Convert numeric values to corresponding 'bin indices' (int)
    x <- num2bin(x, bins)

    if (!elementwise) x <- sort(x) # Important
    ui  <- seq_along(d) # Unique index
    idx <- rep(ui, each = length(bins) - 1) # Index of distribution

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(d)),          # Transition probabilities
                 bins  = bins,                     # Point intersection of bins
                 y     = x,                        # Where to evaluate the pdf
                 prob  = NA_real_,                 # Dummy, only used for 'quantile'
                 type  = "cdf", ncores = ncores, elementwise = elementwise,
                 discrete = FALSE) # <- dummy value

    # Calling C
    check_args_for_treg_predict(args)
    res  <- do.call(function(...) .Call("treg_predict", ...), args)

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Returning simplified result
    } else if (length(d) == 1L & drop) {
        return(res)
    # Else return a data.frame with dimension length(ui) x length(x)
    } else {
        get_colnames <- function(x) {
            x <- paste("p", trimws(format(x, digits = 3)), sep = "_")
            if (anyDuplicated(x)) {
                for (d in unique(x[duplicated(x)])) {
                    idx <- which(x == d)
                    x[idx] <- paste0(x[idx], c("", paste0("_", seq_len(length(idx) - 1))))
                }
            }
            return(x)
        }
        # Create and return matrix
        return(matrix(res, byrow = TRUE, ncol = length(x),
                      dimnames = list(xnames, get_elementwise_colnames(x, "p"))))
    }
}

quantile.Transition <- function(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(probs) == 1L) probs <- rep(probs, length(x))
        elementwise <- length(x) == length(probs)
    }
    if (elementwise & length(probs) == 1L) probs <- rep(probs, length(x))
    if (elementwise & length(probs) != length(x))
        stop("'elementwise = TRUE' but number of probs does not match the number of distributions")

    # Store element names for return
    xnames <- names(x)
    bins   <- as.numeric(attr(x, "bins"))

    # Discrete distributions?
    discrete <- is_discrete(x)

    # Number of probabilities
    nprobs <- length(probs)

    if (elementwise) probs <- sort(probs) # Important
    ui  <- seq_along(x) # Unique index
    idx <- rep(ui, each = length(bins) - 1) # Index of distribution

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(x)),          # Transition probabilities
                 bins  = bins,                     # Point intersection of bins
                 y     = NA_integer_,              # Dummy, only used for cdf/pdf
                 prob  = probs,                    # Probabilities where to evaluate the distribution
                 type  = "quantile", ncores = ncores, elementwise = elementwise,
                 discrete = as.logical(discrete))

    # Calling C
    check_args_for_treg_predict(args)
    res  <- do.call(function(...) .Call("treg_predict", ...), args)

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Returning simplified result
    } else if (length(x) == 1L & drop) {
        return(res)
    # Else return matrix
    } else {
        return(matrix(res, byrow = TRUE, ncol = length(probs),
                      dimnames = list(xnames, get_elementwise_colnames(probs, NULL))))
    }
}


median.Transition <- function(x, na.rm = NULL, ncores = NULL, ...) {
    quantile(x, probs = 0.5, ncores = ncores, ...)
}


mean.Transition <- function(x, ncores = NULL, ...) {

    ## TODO(R): Not correct if the distribution does not cover
    ## the entire range of the data. Throw a warning for now.
    warning("mean.Transition is only experimental. The result is incorrect ",
            "if the distribution provided does not span the full support/range ",
            "of the response distribution.")

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    bins <- attr(x, "bins")
    ui   <- seq_along(x) # Unique index
    idx  <- rep(ui, each = length(bins) - 1) # Index of distribution
    discrete <- rep(FALSE, length(ui))
    warning("TODO(R): Currently assuming discrete = TRUE in mean.Transition")

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(x)),          # Transition probabilities
                 bins  = bins,                     # Point intersection of bins
                 y     = NA_integer_,              # <- Dummy value
                 prob  = NA_real_,                 # <- Dummy value
                 type  = "mean", ncores = ncores,
                 elementwise = TRUE, discrete = discrete) # <- Dummy values

    # Calling C
    check_args_for_treg_predict(args)
    res  <- do.call(function(...) .Call("treg_predict", ...), args)

    setNames(res, names(x))
}


## Draw random values
random.Transition <- function(x, n = 1L, drop = TRUE, ...) {

    # Calculating 'bin mids'
    bins <- attr(x, "bins")
    i <- seq_len(length(bins) - 1)
    binmid <- (bins[i + 1] + bins[i]) / 2.
    binwidth <- diff(bins)

    # Logical vector, is the distribution discrete?
    discrete <- is_discrete(x)

    # Helper function, draw weighted sample of length 'n'.
    # Scoping 'x', 'n', 'binmid', 'binwidth'
    fn <- function(i) {
        y <- as.matrix(x[i], expand = TRUE)
        p <- as.numeric(pdf(x[i], binmid))
        r <- sample(binmid, size = n, prob = p, replace = TRUE)
        # Adding random uniform error
        if (!discrete[i])
            r <- r + runif(n, min = -binwidth[r] / 2, max = +binwidth[r] / 2)
        return(r)
    }
    res <- lapply(seq_along(x), fn)

    # Only one distribution: Return numeric vector
    if (length(res) == 1) {
        return(res[[1]])
    # Only one random value per distribution: Return vector as well
    } else if (n == 1) {
        return(unlist(res))
    # Else return matrix
    } else {
        res <- do.call(rbind, res)
        dimnames(res) <- list(names(x), paste("r", seq_len(ncol(res)), sep = "_"))
        return(res)
    }
}

## Check if distribution is discrete
is_discrete.Transition <- function(d, ...) {
    x <- attr(d, "bins")
    # Calculating mid of bins
    idx <- seq_len(length(x) - 1)
    x <- x[idx + 1] - x[idx]
    # If all 'bin mids' integer we assume it is a discrete distribution
    rep(all(abs(x %% 1) < sqrt(.Machine$double.eps)), length(d))
}

## Check if distribution is continuous
is_continuous.Transition <- function(d, ...) {
    return(!is_discrete(d))
}

## Support (bin range) of the distributions
support.Transition <- function(d, drop = NULL, ...) {
    x <- setNames(range(attr(d, "bins")), c("min", "max"))
    if (length(x) > 1) {
        x <- matrix(x, nrow = length(d), ncol = 2, byrow = TRUE,
                    dimnames = list(names(d), names(x)))
    }
    return(x)
}


newresponse.transitreg <- function(object, newdata = NULL, ...) {
    ## Response name
    yn <- object$response

    if (is.null(newdata)) {
        newdata <- object$model.frame
        newdata[[yn]] <- object$bins[newdata[[yn]]]
    }

    if (is.null(newdata[[object$response]]))
        stop("response missing in newdata!")

    y <- newdata[[object$response]]
    return(y)
}

