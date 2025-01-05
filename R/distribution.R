

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
    k <- tmdist_convert_input(x)
    print(k)
    stop()

    stopifnot(is.numeric(probs), is.numeric(bins), length(probs) == length(bins))
    n <- length(probs)
    x <- c(setNames(as.list(probs), sprintf("tp_%d", seq_len(n) - 1L)),
           setNames(as.list(bins),  sprintf("bin_%d", seq_len(n) - 1)))
    x <- as.data.frame(x)
    class(x) <- c("tmdist", "distribution")
    return(x)
}

tmdist_convert_input <- function(x) {
    # Unnamed list of length 2; where each entry in x is a vector
    if (is.list(x) && all(sapply(x, is.atomic)) && length(x) == 2L && is.null(names(x))) {
        x <- data.frame(tp = x[[1]], bins = x[[2]])
        stopifnot(all(sapply(x, is.numeric)))
        stopifnot(all(x$tp >= 0 & x$tp <= 1))
        stopifnot(all(diff(x$bins) >= 0))
    } else {
        stop("Problem converting user input to 'tm distribution'")
    }
    return(x)
}

# Convert one tmdist distributions into data.frame
as.data.frame.tmdist <- function(x, ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)
    x <- na.omit(unlist(x)) # Convert to named numeric vector
    p <- x[grepl("^tp_", names(x))] # probs
    names(p) <- gsub("^tp_", "", names(p))
    b <- x[grepl("^bin_", names(x))] # bins
    names(b) <- gsub("^bin_", "", names(b))
    stopifnot(identical(names(p), names(b)))
    x <- data.frame(bin    = as.numeric(names(p)),
                    tp     = as.numeric(p),
                    binmid = as.numeric(b))
    x[order(x$bin), ]
}

# TODO(R): Currently testing for one distribution only
format.tmdist <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)
    if (length(x) < 1L) return(character(0))
    n <- names(x) # Keep for later

    # Extracting probabilites and bins
    x <- as.data.frame(x)

    f <- sprintf("tm distribution (%s:%d <--> %s:%d)",
                 format(min(x$binmid), digits = digits), min(x$bin),
                 format(max(x$binmid), digits = digits), max(x$bin))
    setNames(f, n)
}


log_pdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    return(log(pdf(d, x, drop = drop, elementwise = elementwise, ncores = ncores, ...)))
}

pdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    stopifnot("currently only designed length 1" = length(d) == 1L)

    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) elementwise <- !length(x) == length(d)

    if (!elementwise & length(x) > length(d))
        stop("length of 'x' can't be larger than number of distributions in 'd'",
             "if elementwise = FALSE");

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (length(x) == 1 && length(d) > 1L) x <- rep(x, length(d))

    # Scopes 'd' and 'grd'
    fn <- function(i) {
        grd <- as.list(grd[i, ])
        tmp <- as.data.frame(d[grd$d])
        cutoff <- tmp$bin[which.min(abs(tmp$binmid - grd$p))]
        tmp$index <- grd$index
        tmp[tmp$bin <= cutoff, ]
    }

    # If there is one distribution and one point at which to evaluate this
    # distribution, or 'elementwise' is FALSE we can proceed setting up the
    # data.frame used to call the C function performing the calculations.
    # Note: If length(d) != length(x) the vectors will be recycled.
    if ((length(d) == 1L && length(x) == 1L) || isFALSE(elementwise)) {

        # Create 'grid' for evaluation
        if (length(d) > length(x)) {
            grd <- data.frame(d = seq_along(d), p = rep(x, length(d)))
        } else {
            grd <- data.frame(d = rep(seq_along(d), length(x)), p = x)
        }
        grd$index <- seq_len(nrow(grd)) # Pseudo-index

        ui <- seq_len(nrow(grd))
        d <- do.call(rbind, lapply(ui, fn))
    } else if (elementwise) {
        # We're 
        ui <- seq_along(d) # Number of distributions
        d  <- as.data.frame(d)
        d$index <- 1L ## TODO(R): Only works if length(d) == 1
    } else {
        stop(" --- TODO(R): Can we end up here, and if so, why? Work needed --- ")
    }

    ## Calling C to calculate the required PDFs
    res <- .Call("tm_predict", ui, d$index, p = d$tp, type = "pdf", prob = 42,
                 ncores = ncores, elementwise = elementwise)

    if (elementwise) {
        res <- matrix(res, ncol = length(x),
                      dimnames = list(NULL, paste("x", format(x, digits = 3), sep = "_")))
    }
    return(res)

  #FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), ...)
  #apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


cdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    stopifnot("currently only designed length 1" = length(d) == 1L)

    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (length(x) == 1 && length(d) > 1L) x <- rep(x, length(d))

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) elementwise <- !length(x) == length(d)

    if (!elementwise & length(x) > length(d))
        stop("length of 'x' can't be larger than number of distributions in 'd'",
             "if elementwise = FALSE");

    # Scopes 'd' and 'grd'
    fn <- function(i) {
        grd <- as.list(grd[i, ])
        tmp <- as.data.frame(d[grd$d])
        cutoff <- tmp$bin[which.min(abs(tmp$binmid - grd$p))]
        tmp$index <- grd$index
        subset(tmp, bin <= cutoff)
    }

    # If there is one distribution and one point at which to evaluate this
    # distribution, or 'elementwise' is FALSE we can proceed setting up the
    # data.frame used to call the C function performing the calculations.
    # Note: If length(d) != length(x) the vectors will be recycled.
    if ((length(d) == 1L && length(x) == 1L) || isFALSE(elementwise)) {

        # Create 'grid' for evaluation
        if (length(d) > length(x)) {
            grd <- data.frame(d = seq_along(d), p = rep(x, length(d)))
        } else {
            grd <- data.frame(d = rep(seq_along(d), length(x)), p = x)
        }
        grd$index <- seq_len(nrow(grd)) # Pseudo-index

        ui <- seq_len(nrow(grd))
        d <- do.call(rbind, lapply(ui, fn))
    } else if (elementwise) {
        # We're 
        ui <- seq_along(d) # Number of distributions
        d  <- as.data.frame(d)
        d$index <- 1L ## TODO(R): Only works if length(d) == 1
    } else {
        stop(" --- TODO(R): Can we end up here, and if so, why? Work needed --- ")
    }

    ## Calling C to calculate the required CDFs
    res <- .Call("tm_predict", ui, d$index, p = d$tp, type = "cdf", prob = 42,
                 ncores = ncores, elementwise = elementwise)

    if (elementwise) {
        res <- matrix(res, ncol = length(x),
                      dimnames = list(NULL, paste("x", format(x, digits = 3), sep = "_")))
    }
    return(res)

  #FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), ...)
  #apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

quantile.tmdist <- function(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)

    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) elementwise <- !length(x) == length(x)

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (length(probs) == 1 && length(x) > 1L) probs <- rep(probs, length(x))

    if (!elementwise & length(probs) > length(x))
        stop("length of 'probs' can't be larger than number of distributions in 'd'",
             "if elementwise = FALSE");

    if (elementwise) probs <- sort(probs) # Important

    ui <- seq_along(x)
    x <- as.data.frame(x)
    x$index <- 1L # TODO(R) Only works if length(d) == 1

    ## Calling C to calculate the required quantiles
    res <- .Call("tm_predict", ui, x$index, p = x$tp, type = "quantile", prob = probs,
                 ncores = ncores, elementwise = elementwise)

    # Translate quantile bins to numeric values (binmid)
    bin2num <- function(x, d, p, ewise) {
        idx <- if (ewise) rep(1L, length(p)) else 1L
        x   <- data.frame(bin = res, index = idx)
        res <- merge(d, x, by = c("bin", "index"), all.x = FALSE, all.y = TRUE)
        return(res[order(res$index, res$bin), "binmid", drop = TRUE])
    }
    res <- bin2num(res, x, probs, elementwise)

    # If elementwise: Return matrix of dimension length(d) x length(probs)
    if (elementwise) {
        res <- matrix(res, ncol = length(probs),
                      dimnames = list(NULL, paste("q", format(probs, digits = 3), sep = "_")))
    }

    return(res)
  #FUN <- function(at, d) qempirical(at, y = as.matrix(d), ...)
  #apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


median.tmdist <- function(x, ...) {
    quantile(x, probs = 0.5, ...)
}


## Draw random values
random.tmdist <- function(x, n = 1L, drop = TRUE, ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)

    ## Helper function, draw weighted sample of length 'n'.
    ## Scoping 'x' and 'n'
    fn <- function(i) {
        tmp     <- as.data.frame(x[i])
        tmp$pdf <- as.numeric(pdf(x[i], tmp$binmid))
        sample(tmp$binmid, size = n, prob = tmp$pdf, replace = TRUE)
    }
    res <- lapply(seq_along(x), fn)

    res <- if (length(res) > 1) rbind(res) else res[[1]]
    return(res)
}

mean.tmdist <- function(x, ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)
    d <- as.data.frame(x)
    d$pdf <- as.numeric(pdf(x, d$binmid))
    print(d)
    return(sum(d$pdf * d$binmid))
}
