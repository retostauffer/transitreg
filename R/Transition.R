
#' Creates a Transition Distribution
#'
#' A 'Transition' distrubiton consists of a series of \code{K} transition
#' probabilities (TP) for \code{K} 'bins' (counts or pseudo-counts).
#'
#' @param x Numeric vector or a numeric matrix.
#' @param breaks numeric vector of points of intersection of the breaks.
#'        The length the vector must be of \code{length(x) + 1} (if \code{x} is
#'        a vector) or \code{ncol(x) + 1} if \code{x} is a matrix. Must be
#'        monotonically increasing.
#' @param discrete `NULL` or logical. If `NULL` and the breaks look like
#'        count data breaks (`c(-0.5, 0.5, 1.5, ...)`) discrete will be set
#'        `TRUE`, else `FALSE`.
#'
#' @return Returns an object of class \code{c("Transition", "distribution")}.
#' TODO(R): Missing.
#'
#' @details
#' TODO(R): Missing
#'
#' @examples
#' # TODO(R) Write example
#'
#' @concept Transition distribution
#' @concept distributions
#'
#' @importFrom stats setNames
#'
#' @author Reto
#' @rdname Transition
#' @export
Transition <- function(x, breaks, discrete = NULL) {

    # Sanity checks
    stopifnot(
        "'x' must be numeric (vector or matrix)" = is.numeric(x) && is.atomic(x),
        "'x' must be a vector or a matrix" = is.vector(x) || is.matrix(x),
        "length of 'x' must be > 0" = length(x) > 0L,
        "'breaks' must be a numeric vector" = is.atomic(breaks) && is.numeric(breaks),
        "missing values in 'breaks' not allowed" = all(!is.na(breaks)),
        "'breaks' must be monotonically increasing" = all(diff(breaks) > 0),
        "'discrete' must be NULL or logical TRUE/FALSE" = is.null(discrete) || (isTRUE(discrete) || isFALSE(discrete))
    )

    # If 'x' is a vector, convert to matrix
    if (is.vector(x)) x <- matrix(x, nrow = 1)

    # Checking 'breaks' vector
    if (length(breaks) != (ncol(x) + 1))
        stop("'breaks' must be of length ", ncol(x) + 1)

    # Ensure to convert to double in case input is integer
    x[,] <- as.numeric(x)
    breaks <- as.numeric(breaks)

    # Guessing 'discrete' option if needed.
    if (is.null(discrete)) {
        mids <- (head(breaks, -1) + tail(breaks, -1)) / 2
        if (all(abs(mids %% 1 < sqrt(.Machine$double.eps)) &
                isTRUE(all.equal(mids, seq_along(mids) - 1)))) {
            discrete <- TRUE
        } else {
            discrete <- FALSE
        }
    }

    res <- setNames(as.data.frame(x),
                    paste("tp", seq_len(ncol(x)) - 1, sep = "_"))

    structure(res, class = c("Transition", "distribution"), breaks = breaks, discrete = discrete)
}



#' @param log,log.p Logical, if `TRUE`, probabilities `p` are given as `log(p)`.
#'
#' @rdname Transition
#' @export
dtransit <- function(x, d, log = FALSE, ncores = NULL) {
    stopifnot(
        "'x' must be numeric of length > 0" = is.numeric(x) && length(x) > 0L,
        "missing values in 'x' not allowed" = all(!is.na(x)),
        "'d' must be of class 'Transition'" = inherits(d, "Transition")
    )
    log <- as.logical(log[1])
    stopifnot("'log' must evaluate to TRUE/FALSE" = isTRUE(log) || isFALSE(log))

    ## Evaluate cdf
    res <- dpq_get_results(x, d, ncores, type = "pdf")

    return(if (log) log(res) else res)
}


#' @param lower.tail Logical, if `TRUE` (default), probabilities are
#'        `P[X <= x]` otherwise, `P[X > x]`.
#'
#' @rdname Transition
#' @export
ptransit <- function(x, d, lower.tail = TRUE, log.p = FALSE, ncores = NULL) {
    stopifnot(
        "'x' must be numeric of length > 0" = is.numeric(x) && length(x) > 0L,
        "missing values in 'x' not allowed" = all(!is.na(x)),
        "'d' must be of class 'Transition'" = inherits(d, "Transition")
    )
    lower.tail <- as.logical(lower.tail[1])
    log.p <- as.logical(log.p[1])
    stopifnot(
        "'log.p' must evaluate to TRUE/FALSE" = isTRUE(log.p) || isFALSE(log.p),
        "'log.p' must evaluate to TRUE/FALSE" = isTRUE(log.p) || isFALSE(log.p)
    )

    ## Evaluate cdf
    res <- dpq_get_results(x, d, ncores, type = "cdf")

    if (!lower.tail) res <- 1.0 - res
    return(if (log.p) log(res) else res)
}

#' @param p Numeric, vector of probabilities.
#'
#' @rdname Transition
#' @export
qtransit <- function(p, d, lower.tail = TRUE, log.p = FALSE, ncores = NULL) {
    stopifnot(
        "'p' must be numeric of length > 0" = is.numeric(p) && length(p) > 0L,
        "missing values in 'p' not allowed" = all(!is.na(p)),
        "'d' must be of class 'Transition'" = inherits(d, "Transition")
    )

    # Evaluate quantiles
    res <- dpq_get_results(p, d, ncores, type = "quantile")
    return(res)
}


# Input 'z' is either 'x' (dtransit, ptransit), or 'p' (qtransit)
dpq_get_results <- function(z, d, ncores, type) {
    stopifnot(
        is.numeric(z), inherits(d, "Transition"),
        is.character(type) && length(type) == 1L
    )

    ## If type != 'quantile', z represents the numeric values where
    ## to evaluate the distribution(s). Convert to 'bins'
    if (type != "quantile") z <- num2bin(z, attr(d, "breaks"))

    ## If length(y) == 1: Recycle
    if (length(z) == 1L) z <- rep(z, length.out = length(d))

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    ## Store element names for return
    breaks <- attr(d, "breaks")

    ui  <- seq_along(d) # Unique index
    idx <- rep(ui, each = length(breaks) - 1) # Index of distribution

    # If length(y) == length(d): elementwise = TRUE
    elementwise <- length(z) == length(d)


    # Setting up arguments to call .C predict function
    args <- list(uidx    = ui,                       # Unique distribution index (int)
                 idx     = idx,                      # Index vector (int)
                 tp      = t(as.matrix(d)),          # Transition probabilities
                 breaks  = breaks,                    # Point intersection of breaks
                 type    = type, ncores = ncores, elementwise = elementwise,
                 discrete = rep(FALSE, length(ui))) # <- dummy value
    if (type == "quantile") {
        args$y    <- NA_integer_
        args$prob <- if (!elementwise) unique(sort(z)) else z # If !elementwise: take sorted unique
    } else {
        args$y    <- if (!elementwise) unique(sort(z)) else z # If !elementwise: take sorted unique
        args$prob <- NA_real_
    }

    # Calling C
    args <- check_args_for_treg_predict(args)
    res <- do.call(function(...) .Call("treg_predict", ...), args)

    # If !elementwise and length(d) != 1 we need to pick the correct elements
    if (!elementwise & length(d) != 1L) res <- dpq_get_elements(res, z, length(d))
    return(res)
}


# Helper function used in dtransit, ptransit, and qtransit.
# For performance, dpq are always calculated with elementwise = FALSE.
# Depending on length(x)|length(p) and length(d) we then need to
# extract the correct elements for the final return.
#
# As an exmaple: Imagine calling:
#   dtransit(c(1, 5, 2), d)
# .. where d is of length 2. We calculate the density at position 1, 5, 2 for
# all distributions in d, convert the result into a matrix, and then pick the
# correct elements afterwards. As the results are sorted in the matrix the
# order of 'x' is used. In this scenario the resulting matrix is:
#
#         xs_1    xs_2     xs_3
# d_1     r_11    r_12     r_13
# d_2     r_21    r_22     r_23
#
# ... with xs -> sort(x). We will create arr.ind as follows:
#
# arr.ind     row   col
#             1     1
#             2     3
#             1     2
#
# Thus, returning a numeric vector with elements c(r_11, r_23, r_12).
dpq_get_elements <- function(res, x, nd) {
    res  <- matrix(res, nrow = nd, byrow = TRUE)
    mtch <- match(x, sort(unique(x)))
    n    <- max(length(x), nd)
    arr.ind <- cbind(row = rep(seq_len(nd), length.out = n),
                     col = rep(mtch,        length.out = n))
    res <- res[arr.ind]
}


#' @rdname Transition
#' @export
rtransit <- function(n, d, ncores = NULL) {
    n <- as.integer(n)
    if (length(n) > 1) n <- length(n)
    stopifnot(
        "'n' must be of length > 0" = length(n) > 0L,
        "missing values in 'n' not allowed" = all(!is.na(n)),
        "'n' must be integer > 0" = all(n > 0L)
    )

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    # Calculating 'bin mids'
    breaks <- attr(d, "breaks")
    binmid <- (head(breaks, -1) + tail(breaks, -1)) / 2.
    binwidth <- diff(breaks)

    # Logical vector, is the distribution discrete?
    discrete <- is_discrete(d)

    # Calculating densities for all distributions
    ui  <- seq_along(d)
    idx <- rep(seq_along(binmid), length(d))
    y   <- seq_along(binmid) - 1L
    args <- list(uidx   = ui,                       # Unique distribution index (int)
                 idx    = idx,                      # Index vector (int)
                 tp     = t(as.matrix(d)),          # Transition probabilities
                 breaks = breaks,                     # Point intersection of breaks
                 y      = y,                        # Where to evaluate the pdf
                 prob   = NA_real_,                 # <- Dummy value
                 type   = "pdf", ncores = ncores,
                 elementwise = FALSE, discrete = discrete) # <- Dummy values

    args <- check_args_for_treg_predict(args)
    p <- matrix(do.call(function(...) .Call("treg_predict", ...), args),
                ncol = length(binmid), byrow = TRUE)

    # Helper function, draw weighted sample of length 'n'.
    # Scoping 'x', 'n', 'binmid', 'binwidth'
    fn <- function(i) {
        y <- as.matrix(d[i], expand = TRUE)
        p <- as.vector(p[i, ])
        r <- sample(binmid, size = n, prob = p, replace = TRUE)
        # Adding random uniform error
        if (!discrete[i])
            r <- r + runif(n, min = -binwidth[r] / 2, max = +binwidth[r] / 2)
        return(r)
    }
    res <- lapply(seq_along(d), fn)
    res <- do.call(rbind, res)
    return(as.vector(t(res)))
}


mean_transit <- function(x, ncores = NULL, ...) {
    ## TODO(R): Not correct if the distribution does not cover
    ## the entire range of the data. Throw a warning for now.
    warning("mean.Transition is only experimental. The result is incorrect ",
            "if the distribution provided does not span the full support/range ",
            "of the response distribution.")

    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    breaks <- attr(x, "breaks")
    ui   <- seq_along(x) # Unique index
    idx  <- rep(ui, each = length(breaks) - 1) # Index of distribution
    discrete <- rep(FALSE, length(ui))
    warning("TODO(R): Currently assuming discrete = TRUE in mean.Transition")

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(x)),          # Transition probabilities
                 breaks  = breaks,                     # Point intersection of breaks
                 y     = NA_integer_,              # <- Dummy value
                 prob  = NA_real_,                 # <- Dummy value
                 type  = "mean", ncores = ncores,
                 elementwise = TRUE, discrete = discrete) # <- Dummy values

    # Calling C
    args <- check_args_for_treg_predict(args)
    return(do.call(function(...) .Call("treg_predict", ...), args))
}


#' @param \dots objects to be concatenated. Must all come from the same
#'        transition distribution (i.e., be based on the same binning).
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method c Transition
c.Transition <- function(...) {
    x <- list(...)
    if (length(x) == 1) return(x[[1]])

    # Else check whether or not we can combine the objects
    for (i in seq.int(2, length(x))) {
        stopifnot("input not of class Transition" = inherits(x[[i]], "Transition"))
        if (!all.equal(attr(x[[1]], "breaks"), attr(x[[2]], "breaks")))
            stop("breaks of the ", i, ifelse(i == 2, "nd", "th"),
                 "object not the same as for the first object. Can't be combined.")
    }

    # Combine and return
    res <- do.call(rbind, lapply(x, as.matrix))
    Transition(res, attr(x[[1]], "breaks"))
}


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
#' @rdname Transition
#' @exportS3Method as.matrix Transition
as.matrix.Transition <- function(x, expand = FALSE, ...) {
    stopifnot("'expand' must be logical TRUE or FALSE" = isTRUE(expand) || isFALSE(expand))

    xnames <- names(x) # Keep for later
    breaks   <- attr(x, "breaks")
    discrete <- attr(x, "discrete")

    # convert to data.frame -> matrix
    x <- as.matrix(structure(x, class = "data.frame"))
    rownames(x) <- xnames
    if (expand) {
        # Keep original dimension as we transpose(x) in a second
        nd <- nrow(x)          # Number of distributions
        nb <- length(breaks) - 1 # Number of breaks

        index <- rep(seq_len(nd), each = nb) # distribution index
        theta <- rep(seq_len(nb) - 1, times = nd) # 'bin' index or theta
        x <- cbind(index = index, theta = theta, tp = as.vector(t(x)))
    }

    structure(x, class = c("Transitionmatrix", class(x)), breaks = breaks, discrete = discrete)
}

## TODO(R): Problem when subsetting a matrix issue
### #' @exportS3Method "[" Transitionmatrix
### `[.Transitionmatrix` <- function(x, i, j, drop = TRUE, ...) {
###     breaks   <- attr(x, "breaks")
###     discrete <- attr(x, "discrete")
###     xclass   <- class(x)
### 
###     # Convert to default matrix
###     class(x) <- class(x)[!class(x) == "Transitionmatrix"]
###     x <- x[i, j, drop = drop]
### 
###     if (!missing(j)) return(x)
###     # Re-create Transitionmatrix object
###     structure(x, class = xclass, breaks = breaks, discrete = discrete)
### }

#### @exportS3Method format Transitionmatrix
###format.Transitionmatrix <- function(x, ...) {
###    class(x) <- class(x)[!class(x) == "Transitionmatrix"]
###
###    # Extracting (and deleting) attributes
###    an <- c("breaks", "discrete")
###    att <- setNames(lapply(an, function(n) attr(x, n)), an)
###    for (n in an) attr(x, n) <- NULL
###
###    # Print
###    print(x)
###    cat("breaks: ", paste(att$breaks, sep = ", "), "\n")
###    cat("discrete: ", att$discrete, "\n")
###
###}

#' @importFrom stats setNames
#'
#' @exportS3Method format Transition
format.Transition <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and breaks
    fmtfun <- function(i) {
        y <- as.matrix(x[i], expand = FALSE)
        if (ncol(y) > 2L) {
            sprintf("Transition_%d(%s, ..., %s)", ncol(y),
                    format(y[1], digits = digits), format(y[ncol(y)], digits = digits))
        } else {
            # Typically unused, that is two bins only!
            sprintf("Transition_%d(%s, %s)", ncol(y),
                    format(y[1], digits = digits), format(y[2], digits = digits))
        }
    }
    f <- sapply(seq_along(x), fmtfun)
    setNames(f, xnames)
}


#' @param d An object of class `Transition`.
#' @param x Numeric, vector of quantiles.
#' @param drop Logical. Should the result be simplified to a vector if possible?
#' @param elementwise Logical. Should each distribution in `x` be evaluated at
#'        all elements of `probs` (`elementwise = FALSE`, yielding a matrix)? Or, if
#'        `x` and `probs` have the same length, should the evaluation be done element
#'        by element (`elementwise = TRUE` yielding a vector)? The default of `NULL`
#'        means that `elementwise = TRUE` is used if the lengths match and otherwise
#'        `elementwise = FALSE` is used.
#' @param ncores Number of cores to be used (see [transitreg()] for details).
#'
#' @importFrom distributions3 pdf
#' @importFrom stats setNames
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method pdf Transition
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
    breaks <- attr(d, "breaks")

    # Convert numeric values to corresponding 'bin indices' (int)
    x <- num2bin(x, breaks)

    if (!elementwise) x <- sort(x) # Important
    ui  <- seq_along(d) # Unique index
    idx <- rep(ui, each = length(breaks) - 1) # Index of distribution

    # Setting up arguments to call .C predict function
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(d)),          # Transition probabilities
                 breaks  = breaks,                     # Point intersection of breaks
                 y     = x,                        # Where to evaluate the pdf
                 prob  = NA_real_,                 # Dummy, only used for 'quantile'
                 type  = "pdf", ncores = ncores, elementwise = elementwise,
                 discrete = rep(FALSE, length(ui))) # <- dummy value

    # Calling C
    args <- check_args_for_treg_predict(args)
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

#' @importFrom distributions3 log_pdf
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method log_pdf Transition
log_pdf.Transition <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    return(log(pdf(d, x, drop = drop, elementwise = elementwise, ncores = ncores, ...)))
}

#' @importFrom distributions3 cdf
#' @importFrom stats setNames
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method cdf Transition
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
    breaks   <- attr(d, "breaks")

    # Convert numeric values to corresponding 'bin indices' (int)
    x <- num2bin(x, breaks)

    if (!elementwise) x <- sort(x) # Important
    ui  <- seq_along(d) # Unique index
    idx <- rep(ui, each = length(breaks) - 1) # Index of distribution

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(d)),          # Transition probabilities
                 breaks  = breaks,                     # Point intersection of breaks
                 y     = x,                        # Where to evaluate the pdf
                 prob  = NA_real_,                 # Dummy, only used for 'quantile'
                 type  = "cdf", ncores = ncores, elementwise = elementwise,
                 discrete = rep(FALSE, length(ui))) # <- dummy value

    # Calling C
    args <- check_args_for_treg_predict(args)
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

#' @param probs numeric vector of probabilities with values in `[0,1]`.
#'        (Values up to ‘2e-14’ outside that range are accepted and
#'        moved to the nearby endpoint.) TODO(R): SURE?
#'
#' @importFrom stats quantile setNames
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method quantile Transition
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
    breaks   <- as.numeric(attr(x, "breaks"))

    # Discrete distributions?
    discrete <- is_discrete(x)

    # Number of probabilities
    nprobs <- length(probs)

    if (elementwise) probs <- sort(probs) # Important
    ui  <- seq_along(x) # Unique index
    idx <- rep(ui, each = length(breaks) - 1) # Index of distribution

    ## Calling C to calculate the required values.
    args <- list(uidx  = ui,                       # Unique distribution index (int)
                 idx   = idx,                      # Index vector (int)
                 tp    = t(as.matrix(x)),          # Transition probabilities
                 breaks  = breaks,                     # Point intersection of breaks
                 y     = NA_integer_,              # Dummy, only used for cdf/pdf
                 prob  = probs,                    # Probabilities where to evaluate the distribution
                 type  = "quantile", ncores = ncores, elementwise = elementwise,
                 discrete = as.logical(discrete))

    # Calling C
    args <- check_args_for_treg_predict(args)
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


#' @importFrom stats median
#'
#' @param na.rm Unused.
#' @rdname Transition
#' @exportS3Method median Transition
median.Transition <- function(x, na.rm = NULL, ncores = NULL, ...) {
    quantile(x, probs = 0.5, ncores = ncores, ...)
}

#' @importFrom stats setNames
#'
#' @rdname Transition
#' @exportS3Method mean Transition
mean.Transition <- function(x, ncores = NULL, ...) {
    setNames(mean_transit(x, ncores = ncores), names(x))
}

#' @importFrom distributions3 random
#' @importFrom stats runif
#'
#' @param n Integer `>0`, number of random values to be drawn per distribution.
#' @rdname Transition
#' @exportS3Method random Transition
random.Transition <- function(x, n = 1L, drop = TRUE, ...) {

    # Calculating 'bin mids'
    breaks <- attr(x, "breaks")
    i <- seq_len(length(breaks) - 1)
    binmid <- (breaks[i + 1] + breaks[i]) / 2.
    binwidth <- diff(breaks)

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

#' @importFrom distributions3 is_discrete
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method is_discrete Transition
is_discrete.Transition <- function(d, ...) {
    rep(attr(d, "discrete"), length(d))
}

#' @importFrom distributions3 is_continuous
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method is_continuous Transition
is_continuous.Transition <- function(d, ...) {
    return(!is_discrete(d))
}

#' @importFrom distributions3 support
#' @importFrom stats setNames
#'
#' @author Reto
#' @rdname Transition
#' @exportS3Method support Transition
support.Transition <- function(d, drop = NULL, ...) {
    x <- setNames(range(attr(d, "breaks")), c("min", "max"))
    if (length(x) > 1) {
        x <- matrix(x, nrow = length(d), ncol = 2, byrow = TRUE,
                    dimnames = list(names(d), names(x)))
    }
    return(x)
}


#' @param type Character, type of plot to be created. Either
#'        `"tp"` (transition probabilities), `"cdf"` or `"pdf"`.
#' @param all Logical. If `TRUE` all distributions in `x` will be drawn.
#'        Else only the first N ones (default is 8).
#' @param n Integer. Maximum number of distributions to be plotted, defaults to 8.
#'
#' @importFrom utils head tail
#' @importFrom graphics matplot axis
#' @exportS3Method plot Transition
#' @rdname Transition
plot.Transition <- function(x, type = c("tp", "cdf", "pdf"), all = FALSE, n = 8L, ...) {
    stopifnot("'all' must be TRUE or FALSE" = isTRUE(all) || isFALSE(all))
    n <- as.integer(n)[1]
    stopifnot("'n' cannot be coerced to integer > 0L" =
              is.integer(n) && length(n) == 1L && n > 0L)

    type <- match.arg(type)
    titles <- c("tp"  = "Transition Probabilities",
               "cdf" = "Distribution",
               "pdf" = "Density")
    breaks <- attr(x, "breaks")

    # Take first 1:n distributions only
    if (length(x) > n & !all) x <- x[seq_len(n)]

    binmid <- (head(breaks, -1) + tail(breaks, -1)) / 2 # Mid of bin
    if (type == "tp") {
        m  <- as.matrix(x)
    } else if (type == "cdf") {
        m  <- cdf(x, binmid, elementwise = FALSE, drop = FALSE)
    } else {
        m  <- pdf(x, binmid, elementwise = FALSE, drop = FALSE)
    }

    matplot(x = binmid, y = t(m), type = "l",
            lty = 1, main = titles[type], ...)

    axis(side = 1, at = breaks, labels = FALSE, col = 1, tck = 0.015)
    invisible(NULL)
}







