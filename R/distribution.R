
#' Creates a tm Distribution
#'
#' A tm distrubiton consists of a series of \code{K} transition
#' probabilities for \code{K} 'bins' (counts or pseudo-counts) and
#' a series of \code{K} numeric values representing the (center of) the
#' corresponding bins.
#'
#' @param x one (or multiple) transitionmodel distributions. See section
#'        'Input' for more details.
#' @param probs numeric vector with transition probabilities (0, ..., K).
#' @param bins numeric vector with numeric value for bin 0, ..., K.
#'
#' @section Input:
#' In contrast to may parametric distributions, a 'tm' distribution
#' does not consist of a set of distribution parameters, but of a set
#' of transition probabilities (\code{tp}) and a set of 'bins' (same length
#' as \code{tp}). The \code{tp}s define the probability that that an observation
#' is larger than the current one (in terms of \code{bins}).
#'
#' This constructor function allows for different input formats for convenience,
#' and tries to accomodate everything that easy to detect and not ambiguous.
#' For all cases covered the following conditions must be met:
#'
#' * Both (\code{tp}, \code{bins}) must be numeric and of same length.
#' * \code{tp} must be in the range of \code{[0, 1]}
#'   (for now; TODO(R): Do we always need the full distribution or is it enough if
#'    we allow for 'partial' ones? If cut at the upper end it may be OK, but
#'    if cut at the lower that goes wrong at all times, right?).
#' * \code{bins} must be monotonically increasing (sorted).
#'
#' Input argument \code{x} can be:
#'
#' 1. Named list or data.frame with two elements: \code{tp} and \code{bins}.
#' 2. Unnamed list of length 2, expecting the first entry to be \code{tp},
#'    the second \code{bins}.
#'
#'
#' @return Returns an object of class \code{c("tmdist", "distribution")}.
#' @author Reto
tmdist <- function(x) {
    # Converting input into a list of data.frames. The lengt of the list
    # corresponds to the number of distributions, whereof each element is a
    # data.frame with the transition probabilities and bins.

    # If the input is a Transition Model object, create distributions
    # based on the data used for training.
    if (inherits(x, "tm")) {
        stop(" --- here call procast --- ")
        vars <- attr(terms(TransitionModels:::fake_formula(formula(b))), "term.labels")
        vars <- vars[!vars == "theta"]
        d    <- x$model.frame
        nb   <- length(x$bins)
        tmp  <- lapply(d[, vars, drop = FALSE], function(x, nb) rep(x, each = nb), nb = nb)
        print(str(tmp))
        print(nb)
        tmp  <- data.frame(c(list(theta = rep(seq_len(nb) - 1, times = nrow(d))), tmp))

        tmp  <- data.frame(tp   = predict(x$model, newdata = tmp, type = "response"),
                           bins = rep(x$bins, times = nrow(d)))
        # Split into individual data.frames
        x <- split(tmp, rep(seq_len(nrow(d)), each = nb))
        rm(d, tmp)
    # Else try to convert the input (different formats are possible)
    # into a distributions object.
    } else {
        x <- tmdist_convert_input(x)
        if (is.data.frame(x)) x <- list(x)
    }

    # Sanity check on the newly created list of data.frames
    # Sanity check my data.frames
    check_dfs <- function(y) {
        stopifnot(all(y$tp >= 0 & y$tp <= 1))
        stopifnot(all(diff(y$binmid) >= 0))
        stopifnot(all(sapply(y, is.numeric)))
    }
    lapply(x, check_dfs)

    # Find max length of the dfs
    nmax <- max(sapply(x, nrow))

    # Helper function to create the names for the matrix (and the matrix inserts)
    get_names <- function(n)
        c(sprintf("tp_%d",  seq_len(n) - 1L), sprintf("bin_%d", seq_len(n) - 1L))

    # Convert from long to wide matrix; scopes 'mat' and 'x'
    df_to_mat <- function(y, m) {
        y <- setNames(unlist(y), get_names(nrow(y))) # Create named vector
        m[, names(y)] <- unname(y) # Inserting values into matrix
        return(m)
    }
    m <- matrix(NA, 1, ncol = 2 * nmax, dimnames = list(NULL, get_names(nmax)))

    # Converting to data.frame, adding custom class, and return.
    res <- lapply(x, df_to_mat, m = m)
    res <- as.data.frame(do.call(rbind, res), expand = TRUE)

    class(res) <- c("tmdist", "distribution")
    return(res)
}


prodist.tm <- function(object, newdata = NULL, ...) {
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
    nb <- length(object$bins) # Number of bins
    nd <- nrow(res)           # Number of observations

    expand_covar <- function(x, nb) rep(x, each = nb)
    res  <- lapply(res[, covars, drop = FALSE], expand_covar, nb = nb)
    res  <- data.frame(c(list(theta = rep(seq_len(nb) - 1, times = nd)), res))

    # TODO(R): Currently 'type = response' which differs
    # for different engines (see tm()).
    res  <- data.frame(tp     = predict(object$model, newdata = res, type = "response"),
                       binmid = rep(object$bins, times = nd))
    # Split into individual data.frames
    res <- split(res, rep(seq_len(nd), each = nb))

    # Convert into distributions object, pack into a data.frame
    return(setNames(as.data.frame(tmdist(res)), as.character(substitute(object))))
}

procast.tm <- function(object, newdata = NULL, na.action = na.pass,
                       type = "distribution", at = 0.5, drop = FALSE, ...) {
    object <- if (is.null(newdata)) {
        prodist(object)
    } else {
        # TODO(R): na.action passed to 'procast.tm' but not yet implemented (no effect)
        prodist(object, newdata = newdata, na.action = na.action)
    }

    return(if (drop) object[[1]] else object)
}


tmdist_convert_input <- function(x) {
    # Unnamed list of length 2; where each entry in x is a vector
    if (is.list(x) && all(sapply(x, is.atomic)) && length(x) == 2L && is.null(names(x))) {
        x <- data.frame(tp = x[[1]], binmid = x[[2]])
    # Named list
    } else if (is.list(x) && all(c("tp", "binmid") %in% names(x))) {
        stopifnot(length(x$binmid) == length(x$tp))
        x <- as.data.frame(x[c("tp", "binmid")], expand = TRUE)
    # Else we expect to have gotten a series of distributions,
    # so we call this function again to convert them if possible
    } else if (is.list(x)) {
        # TODO(R): That is no safe idea!
        x <- lapply(x, tmdist_convert_input)
    # Else we don't really know what to do.
    } else {
        stop("Problems converting user input to 'tm distribution'")
    }
    #lapply(x, check_values)
    return(x)
}

# Convert one tmdist distributions into data.frame
as.data.frame.tmdist <- function(x, expand = FALSE, ...) {

    if (expand) {
        class(x) <- "data.frame"
        idx_tp <- grepl("^tp_", names(x))
        idx_bm <- grepl("^bin_", names(x))
        binid  <- names(x)[idx_tp]
        binid  <- as.integer(regmatches(binid, regexpr("[0-9]+$", binid)))

        x <- t(as.matrix(x))
        res <- data.frame(index  = rep(seq_len(ncol(x)), each = length(binid)),
                          bin    = rep(binid, ncol(x)),
                          tp     = as.numeric(x[idx_tp, ]),
                          binmid = as.numeric(x[idx_bm, ]))
        return(na.omit(res))
    } else {
        return(NextMethod())
    }
}

# TODO(R): Currently testing for one distribution only
format.tmdist <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and bins
    fmt <- function(i, x) {
        y <- as.data.frame(x[i], expand = TRUE)
        sprintf("tm distribution (%s:%d <--> %s:%d)",
                format(min(y$binmid), digits = digits), min(y$bin),
                format(max(y$binmid), digits = digits), max(y$bin))
    }

    f <- sapply(seq_along(x), fmt, x = x)
    setNames(f, xnames)
}


log_pdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    return(log(pdf(d, x, drop = drop, elementwise = elementwise, ncores = ncores, ...)))
}

pdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {

    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(x) == 1L) x <- rep(x, length(d))
        elementwise <- length(d) == length(x)
    }
    if (elementwise & length(d) != length(x))
        stop("'elementwise = TRUE' but length of x does not match the number of distributions")

    # Store element names for return
    xnames <- names(x)

    if (!elementwise) x <- sort(x) # Important
    ui <- seq_along(d) # Unique index
    d  <- as.data.frame(d, expand = TRUE) # convert distributions to data.frame

    ## Calling C to calculate the required quantiles
    ## binmid: Required as we need to know where to stop.
    ## y: threshold at which to evaluate the pdf.
    res <- .Call("tm_predict", uidx = ui, idx = d$index, tp = d$tp,
                 binmid = d$binmid, y = x, type = "pdf", ncores = ncores,
                 elementwise = elementwise, discrete = FALSE)

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Else return a matrix with dimension length(ui) x length(x)
    } else {
        # Create and return matrix
        return(matrix(res, byrow = TRUE, ncol = length(x),
                      dimnames = list(xnames, get_mat_colnames(x, "d"))))
    }
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
get_mat_colnames <- function(x, prefix = NULL, digits = 3) {
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

cdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {

    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(x) == 1L) x <- rep(x, length(d))
        elementwise <- length(d) == length(x)
    }
    if (elementwise & length(d) != length(x))
        stop("'elementwise = TRUE' but length of x does not match the number of distributions")

    # Store element names for return
    xnames <- names(x)

    if (!elementwise) x <- sort(x) # Important
    ui <- seq_along(d) # Unique index
    d  <- as.data.frame(d, expand = TRUE) # convert distributions to data.frame

    ## Calling C to calculate the required CDFs
    res <- .Call("tm_predict", uidx = ui, idx = d$index, tp = d$tp,
                 binmid = d$binmid, y = x, type = "cdf", ncores = ncores,
                 elementwise = elementwise, discrete = FALSE)


    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
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
                      dimnames = list(xnames, get_mat_colnames(x, "p"))))
    }
}

quantile.tmdist <- function(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) {
        if (length(probs) == 1L) probs <- rep(probs, length(x))
        elementwise <- length(x) == length(probs)
    }
    if (elementwise & length(probs) != length(x))
        stop("'elementwise = TRUE' but number of probs does not match the number of distributions")

    # Store element names for return
    xnames <- names(x)

    # Discrete distributions?
    discrete <- is_discrete(x)
    discrete <- TRUE
    discrete <- FALSE

    # Number of probabilities
    nprobs <- length(probs)

    if (elementwise) probs <- sort(probs) # Important
    ui <- seq_along(x) # Unique index
    x  <- as.data.frame(x, expand = TRUE) # convert distributions to data.frame

    ## Calling C to calculate the required quantiles
    res <- .Call("tm_predict", uidx = ui, index = x$index, tp = x$tp,
                 binmid = x$binmid, y = probs, type = "quantile", ncores = ncores,
                 elementwise = elementwise, discrete = discrete)

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Else return matrix
    } else {
        return(matrix(res, byrow = TRUE, ncol = length(probs),
                      dimnames = list(xnames, get_mat_colnames(probs, NULL))))
    }
}


median.tmdist <- function(x, ...) {
    quantile(x, probs = 0.5, ...)
}


mean.tmdist <- function(x, ncores = NULL, ...) {
    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    ui <- seq_along(x) # Number of distributions
    x  <- as.data.frame(x, expand = TRUE) # Convert to data.frame

    ## Calling C to calculate the required quantiles
    return(.Call("tm_predict", uidx = ui, index = x$index, tp = x$tp, binmid = x$binmid, y = NA_real_,
                 type = "mean", ncores = ncores, elementwise = TRUE))
}


## Draw random values
random.tmdist <- function(x, n = 1L, drop = TRUE, ...) {
    # Helper function, draw weighted sample of length 'n'.
    # Scoping 'x' and 'n'
    fn <- function(i) {
        tmp     <- as.data.frame(x[i], expand = TRUE)
        tmp$pdf <- as.numeric(pdf(x[i], tmp$binmid))
        delta   <- diff(tmp$binmid[1:2]) / 2 # bin width/2
        sample(tmp$binmid, size = n, prob = tmp$pdf, replace = TRUE) +
            runif(n, -delta, delta)
        # TODO(R): Currently adding +/- uniform random error
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
is_discrete.tmdist <- function(x, ...) {
    ## Guessing from 'bins'. If all bins are integers,
    ## we assume the original distribution has been discrete.
    ## Else continuous.
    class(x) <- "data.frame"
    x <- as.matrix(x[, grep("^bin_[0-9]+$", names(x))])
    x <- na.omit(unique(as.numeric(x)))
    ## Fraction of non-integer values
    tmp <- mean(x %% 1 > sqrt(.Machine$double.eps))
    ## If fraction of non-perfect integers is < 1 promille, we assume
    ## this is/was a discrete distribution (or all of them).
    return(tmp < 1e-3)
}
is_continuous.tmdist <- function(x, ...) {
    return(!is_discrete(x))
}

## TODO(R): Currently adding 'max difference' of any distribution
##          as delta to add on top of the support. This works
##          fine if all the distributions have the same bins (or
##          same diff(bins), any better solution? Question for Niki
support.tmdist <- function(x, ...) {
    class(x) <- "data.frame"
    x <- as.matrix(x[, grep("^bin_[0-9]+$", names(x))])
    delta <- max(apply(x, MARGIN = 1, function(y) max(diff(y), na.rm = TRUE)))
    return(range(na.omit(unique(as.numeric(x)))) + c(-0.5, 0.5) * delta)
}


newresponse.tm <- function(object, newdata = NULL, ...) {
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

