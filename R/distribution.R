
#' Creates a Transition Distribution
#'
#' A Transition distrubiton consists of a series of \code{K} transition
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
#' In contrast to may parametric distributions, a 'Transition' distribution
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
#' @return Returns an object of class \code{c("Transition", "distribution")}.
#' @author Reto
Transition <- function(z, newdata = NULL, newresponse = NULL) {
    # Converting input into a list of data.frames. The lengt of the list
    # corresponds to the number of distributions, whereof each element is a
    # data.frame with the transition probabilities and bins.

    # If the input is a Transition Model object, create distributions
    # based on the data used for training.
    if (inherits(z, "transitreg")) {
        warning("TODO(R): Here I sould call procast?")
        vars <- attr(terms(fake_formula(formula(x))), "term.labels")
        vars <- vars[!vars == "theta"]
        d    <- x$model.frame
        nb   <- length(x$bins)
        tmp  <- lapply(d[, vars, drop = FALSE], function(x, nb) rep(x, each = nb), nb = nb)
        tmp  <- data.frame(c(list(theta = rep(seq_len(nb) - 1, times = nrow(d))), tmp))

        tmp  <- data.frame(tp   = predict(x$model, newdata = tmp, type = "response"),
                           bins = rep(x$bins, times = nrow(d)))
        # Split into individual data.frames
        x <- split(tmp, rep(seq_len(nrow(d)), each = nb))
        rm(d, tmp)
    # Else try to convert the input (different formats are possible)
    # into a distributions object.
    } else {
        z <- Transition_convert_input(z)
        if (is.matrix(z)) z <- list(z)
    }

    # Sanity check on the newly created list of data.frames
    # Sanity check my data.frames
    check_matrices <- function(y) {
        stopifnot(is.matrix(y), is.numeric(y))
        stopifnot(all(y[, "tp"] >= 0 & y[, "tp"] <= 1))
        stopifnot(all(y[, "lo"] < y[, "up"]))
        stopifnot(all(y[, "up"] - y[, "lo"] > 0))
    }
    lapply(z, check_matrices)

    # Find max length of the dfs
    nmax <- max(sapply(z, nrow))

    # Helper function to create the names for the matrix (and the matrix inserts)
    get_names <- function(n)
        as.vector(outer(c("tp", "lo", "up"), seq_len(n) - 1, paste, sep = "_"))

    # Converting to data.frame, adding custom class, and return.
    to_matrix <- function(y, n) {
        m <- matrix(NA_real_, nrow = 1, ncol = 3L * n)
        m[seq_along(y)] <- as.vector(t(y))
        return(m)
    }
    res <- lapply(z, to_matrix, n = nmax)
    res <- as.data.frame(structure(do.call(rbind, res),
                                   dimnames = list(NULL, get_names(nmax))))

    class(res) <- c("Transition", "distribution")
    return(res)
}


Transition_convert_input <- function(x) {
    # Unnamed list of length 2; where each entry in x is a vector
    if (is.list(x) && all(sapply(x, is.atomic)) && length(x) == 3L && is.null(names(x))) {
        x <- cbind(tp = x[[1]], lo = x[[2]], up = x[[3]])
    # Named list
    } else if (is.list(x) && all(c("tp", "lo", "up") %in% names(x))) {
        stopifnot("length of input vectrs differ" = length(unique(sapply(x, length))) == 1L)
        x <- cbind(tp = x$tp, lo = x$lo, up = x$up)
    # Else we expect to have gotten a series of distributions,
    # so we call this function again to convert them if possible
    } else if (is.list(x)) {
        # TODO(R): Recursive call on level 1; safe option?
        x <- lapply(x, Transition_convert_input)
    # Else we don't really know what to do.
    } else {
        stop("Problems converting user input to 'Transition distribution'")
    }
    #lapply(x, check_values)
    return(x)
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
#'        returned, where the columns contain transition probabilities (\code{tp_})
#'        as well as the lower (\code{lo_}) and upper (\code{up}) bound of the
#'        corresponding bin.
#' @param \dots unused.
#'
#' @return Numeric matrix. If \code{expand = FALSE} the return is of dimension
#' \code{c(length(x), <ncol>)}.
#'
#' If \code{expand = TRUE} the returned matrix is of dimension
#' \code{c(length(x) * <ncol>, 4L)} where the four columns contain 
#' \code{index} (1, ..., length(x)) where each index corresponds to the
#' row-index of the original input \code{x} (distribution identifier),
#' the transition probabilities \code{tp}, as well as two columns containing
#' the lower and upper bound of the bin (\code{lo}, \code{up}).
#' This expanded version is used when calling the .C functions.
#'
#' @author Reto
as.matrix.Transition <- function(x, expand = FALSE, ...) {
    stopifnot("'expand' must be logical TRUE or FALSE" = isTRUE(expand) || isFALSE(expand))

    xnames <- names(x) # Keep for later

    # convert to data.frame -> matrix
    x <- as.matrix(structure(x, class = "data.frame"))
    rownames(x) <- xnames
    if (expand) {
        # Keep original dimension as we transpose(x) in a second
        nd <- nrow(x)      # Number of distributions
        nb <- ncol(x) / 3L # Number of bins

        x <- t(x) # Transpose 'x' to properly extract data

        res <- list(index = rep(seq_len(nd), each = nb))
        for (n in c("tp", "lo", "up")) {
            res[[n]] <- as.vector(x[grep(sprintf("^%s_[0-9]+$", n), rownames(x)), ])
        }
        # Create new rownames if there were any
        if (!is.null(xnames))
            xnames <- paste(rep(xnames, each = nb),
                            rep(seq_len(nb), times = nd), sep = "_")

        # Convert to matrix
        x <- as.matrix(do.call(cbind, res))
        if (!is.null(xnames)) rownames(x) <- xnames

        # Remove missing values
        x <- x[apply(is.na(x), MARGIN = 1, sum) == 0, ]
    }
    return(x)
}

# Format
format.Transition <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and bins
    fmtfun <- function(i) {
        y <- as.matrix(x[i], expand = TRUE)
        sprintf("%d; %s [%s,%s], ... , %s [%s,%s]", nrow(y),
                format(y[1,       "tp"], digits = digits),
                format(y[1,       "lo"], digits = digits),
                format(y[1,       "up"], digits = digits),
                format(y[nrow(y), "tp"], digits = digits),
                format(y[nrow(y), "lo"], digits = digits),
                format(y[nrow(y), "up"], digits = digits))
    }
    f <- sapply(seq_along(x), fmtfun)
    f <- sprintf("Transition Dist (%s)", f)
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

    if (!elementwise) x <- sort(x) # Important
    d  <- as.matrix(d, expand = TRUE) # convert distributions to matrix
    ui <- unique(d[, "index"])

    ## Calling C to calculate the required values.
    res <- .Call("treg_predict",
                 uidx  = as.integer(ui),           # Unique distribution index (int)
                 idx   = as.integer(d[, "index"]), # Index vector (int)
                 tp    = d[, "tp"],                # Transition probabilities
                 lower = d[, "lo"],                # Lower edge of the bin
                 upper = d[, "up"],                # Upper edge of the bin
                 y     = as.numeric(x),            # Where to evaluate the pdf
                 type  = "pdf", ncores = ncores, elementwise = elementwise,
                 discrete = FALSE) # <- dummy value

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
    xnames <- names(x)

    if (!elementwise) x <- sort(x) # Important
    d  <- as.matrix(d, expand = TRUE) # convert distributions to matrix
    ui <- unique(d[, "index"])

    ## Calling C to calculate the required values.
    res <- .Call("treg_predict",
                 uidx  = as.integer(ui),           # Unique distribution index (int)
                 idx   = as.integer(d[, "index"]), # Index vector (int)
                 tp    = d[, "tp"],                # Transition probabilities
                 lower = d[, "lo"],                # Lower edge of the bin
                 upper = d[, "up"],                # Upper edge of the bin
                 y     = as.numeric(x),            # Where to evaluate the pdf
                 type  = "cdf", ncores = ncores, elementwise = elementwise,
                 discrete = FALSE) # <- dummy value


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

quantile.Transition <- function(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

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

    # Number of probabilities
    nprobs <- length(probs)

    if (elementwise) probs <- sort(probs) # Important
    x  <- as.matrix(x, expand = TRUE) # convert distributions to matrix
    ui <- unique(x[, "index"])

    ## Calling C to calculate the required values.
    res <- .Call("treg_predict",
                 uidx  = as.integer(ui),           # Unique distribution index (int)
                 idx   = as.integer(x[, "index"]), # Index vector (int)
                 tp    = x[, "tp"],                # Transition probabilities
                 lower = x[, "lo"],                # Lower edge of the bin
                 upper = x[, "up"],                # Upper edge of the bin
                 y     = as.numeric(probs),        # Where to evaluate the pdf
                 type  = "quantile", ncores = ncores, elementwise = elementwise,
                 discrete = as.logical(discrete))

    # If elementwise: Return named vector
    if (elementwise) {
        return(setNames(res, xnames))
    # Else return matrix
    } else {
        return(matrix(res, byrow = TRUE, ncol = length(probs),
                      dimnames = list(xnames, get_mat_colnames(probs, NULL))))
    }
}


median.Transition <- function(x, na.rm = NULL, ncores = NULL, ...) {
    quantile(x, probs = 0.5, ncores = ncores, ...)
}


mean.Transition <- function(x, ncores = NULL, ...) {
    ## Get number of cores for OpenMP parallelization
    ncores <- transitreg_get_number_of_cores(ncores, FALSE)

    x  <- as.matrix(x, expand = TRUE) # convert distributions to matrix
    ui <- unique(x[, "index"])

    ## Calling C to calculate the required values.
    return(.Call("treg_predict",
                 uidx  = as.integer(ui),           # Unique distribution index (int)
                 idx   = as.integer(x[, "index"]), # Index vector (int)
                 tp    = x[, "tp"],                # Transition probabilities
                 lower = x[, "lo"],                # Lower edge of the bin
                 upper = x[, "up"],                # Upper edge of the bin
                 y     = NA_real_,                 # <- Dummy value
                 type  = "mean", ncores = ncores,
                 elementwise = TRUE, discrete = FALSE)) # <- Dummy values
}


## Draw random values
random.Transition <- function(x, n = 1L, drop = TRUE, ...) {
    # Helper function, draw weighted sample of length 'n'.
    # Scoping 'x' and 'n'
    fn <- function(i) {
        y        <- as.matrix(x[i], expand = TRUE)
        binmid   <- rowMeans(y[, c("lo", "up")])
        p        <- as.numeric(pdf(x[i], binmid))
        binwidth <- y[, "up"] - y[, "lo"]
        sample(binmid, size = n, prob = p, replace = TRUE) +
            runif(n, -binwidth / 2, binwidth / 2)
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
is_discrete.Transition <- function(d, ...) {
    # Calculating 'bin mid', if all bin mid points
    # are integer we assume it is a discrete dist.
    fn <- function(i) {
        y <- as.matrix(d[i], expand = TRUE)
        binmid <- (y[, "lo"] + y[, "up"]) / 2
        all(abs(binmid %% 1) < sqrt(.Machine$double.eps))
    }
    sapply(seq_along(d), fn)
}

## Check if distribution is continuous
is_continuous.Transition <- function(d, ...) {
    return(!is_discrete(d))
}

## Support (bin range) of the distributions
support.Transition <- function(d, drop = NULL, ...) {
    fn <- function(i) {
        y <- as.matrix(d[i], expand = TRUE)
        c(min = min(y[, "lo"]), max = max(y[, "up"]))
    }
    t(sapply(seq_along(d), fn))
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

