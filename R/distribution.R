
#' Creates a Transition Distribution
#'
#' A 'Transition' distrubiton consists of a series of \code{K} transition
#' probabilities (TP) for \code{K} 'bins' (counts or pseudo-counts).
#'
#' @param x numeric vector or a numeric matrix.
#' @param breaks numeric vector of points of intersection of the breaks.
#'        The length the vector must be of \code{length(x) + 1} (if \code{x} is
#'        a vector) or \code{ncol(x) + 1} if \code{x} is a matrix. Must be
#'        monotonically increasing.
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
Transition <- function(x, breaks) {

    # Sanity checks
    stopifnot(
        "'x' must be numeric (vector or matrix)" = is.numeric(x) && is.atomic(x),
        "'x' must be a vector or a matrix" = is.vector(x) || is.matrix(x),
        "length of 'x' must be > 0" = length(x) > 0L,
        "'breaks' must be a numeric vector" = is.atomic(breaks) && is.numeric(breaks),
        "missing values in 'breaks' not allowed" = all(!is.na(breaks)),
        "'breaks' must be monotonically increasing" = all(diff(breaks) > 0)
    )

    # If 'x' is a vector, convert to matrix
    if (is.vector(x)) x <- matrix(x, nrow = 1)

    # Checking 'breaks' vector
    if (length(breaks) != (ncol(x) + 1))
        stop("'breaks' must be of length ", ncol(x) + 1)


    # Ensure to convert to double in case input is integer
    x[,] <- as.numeric(x)
    breaks <- as.numeric(breaks)

    res <- setNames(as.data.frame(x),
                    paste("tp", seq_len(ncol(x)) - 1, sep = "_"))

    structure(res, class = c("Transition", "distribution"), breaks = breaks)
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

#' @importFrom distributions3 prodist
#' @importFrom stats setNames
#' @importFrom utils head tail
#'
#' @author Reto
#' @exportS3Method prodist transitreg
#' @rdname transitreg
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
    nb <- length(object$ym) # Number of breaks
    nd <- nrow(res)         # Number of observations

    expand_covar <- function(x, nb) rep(x, each = nb)
    res  <- lapply(res[, covars, drop = FALSE], expand_covar, nb = nb)
    res  <- data.frame(c(list(theta = rep(seq_len(nb) - 1, times = nd)), res))

    # TODO(R): Currently 'type = response' which differs
    # for different engines (see transitreg()).
    res  <- data.frame(tp = predict(object$model, newdata = res, type = "response"),
                       lo = rep(head(object$breaks, -1), times = nd),
                       up = rep(tail(object$breaks, -1), times = nd))
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
#' @rdname Transition
#' @exportS3Method as.matrix Transition
as.matrix.Transition <- function(x, expand = FALSE, ...) {
    stopifnot("'expand' must be logical TRUE or FALSE" = isTRUE(expand) || isFALSE(expand))

    xnames <- names(x) # Keep for later
    breaks   <- attr(x, "breaks")

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

    structure(x, class = c("Transitionmatrix", class(x)), breaks = breaks)
}

#' @importFrom stats setNames
#'
#' @exportS3Method format Transition
format.Transition <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and breaks
    fmtfun <- function(i) {
        y <- as.matrix(x[i], expand = FALSE)
        sprintf("n = %d", ncol(y))
    }
    f <- sapply(seq_along(x), fmtfun)
    f <- sprintf("Transition(%s)", f)
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
                 discrete = rep(as.logical(discrete), length(ui)))

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
    check_args_for_treg_predict(args)
    res  <- do.call(function(...) .Call("treg_predict", ...), args)

    setNames(res, names(x))
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
    x <- attr(d, "breaks")
    # Calculating mid of breaks
    idx <- seq_len(length(x) - 1)
    x <- x[idx + 1] - x[idx]
    # If all 'bin mids' integer we assume it is a discrete distribution
    rep(all(abs(x %% 1) < sqrt(.Machine$double.eps)), length(d))
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


#' @importFrom utils head tail
#' @importFrom graphics matplot axis
#' @exportS3Method plot Transition
plot.Transition <- function(d, type = c("tp", "cdf", "pdf"), p = c(0.1, 99.9), len = 1000,
                            all = FALSE, ...) {


    type <- match.arg(type)
    titles <- c("tp"  = "Transition Probabilities",
               "cdf" = "Distribution",
               "pdf" = "Density")

    if (length(d) > 8 & !all) d <- d[1:8]

    breaks <- attr(d, "breaks")

    if (type == "tp") {
        x <- (head(breaks, -1) + tail(breaks, -1)) / 2 # Mid of bin
        m <- as.matrix(d)
    } else if (type == "cdf") {
        x <- seq(min(breaks), max(breaks), length.out = len)
        m <- cdf(d, x, elementwise = FALSE, drop = FALSE)
    } else {
        x <- seq(min(breaks), max(breaks), length.out = len)
        m <- pdf(d, x, elementwise = FALSE, drop = FALSE)
    }

    matplot(x = x, y = t(m), type = "l",
            lty = 1, main = titles[type], ...)

    axis(side = 1, at = attr(m, "breaks"), labels = FALSE, col = 1, tck = 0.015)
    invisible(NULL)
}








# TODO(R): Used?

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

    y <- newdata[[object$response]]
    return(y)
}

