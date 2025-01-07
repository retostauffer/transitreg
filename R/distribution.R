
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
    x <- tmdist_convert_input(x)
    if (is.data.frame(x)) x <- list(x)

    # Sanity check on the newly created list of data.frames
    # Sanity check my data.frames
    check_dfs <- function(y) {
        print(str(y))
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
    res <- as.data.frame(do.call(rbind, res))

    class(res) <- c("tmdist", "distribution")
    return(res)
}

tmdist_convert_input <- function(x) {
    # Unnamed list of length 2; where each entry in x is a vector
    if (is.list(x) && all(sapply(x, is.atomic)) && length(x) == 2L && is.null(names(x))) {
        x <- data.frame(tp = x[[1]], binmid = x[[2]])
    # Named list
    } else if (is.list(x) && all(c("tp", "binmid") %in% names(x))) {
        stopifnot(length(x$binmid) == length(x$tp))
        x <- as.data.frame(x[c("tp", "binmid")])
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
as.data.frame.tmdist <- function(x, ...) {

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
}

# TODO(R): Currently testing for one distribution only
format.tmdist <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
    if (length(x) < 1L) return(character(0))
    xnames <- names(x) # Keep for later

    # Extracting probabilites and bins
    fmt <- function(i, x) {
        y <- as.data.frame(x[i])
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

    xnames <- names(x) # Keep for later

    # Number of probabilities
    nx <- length(x)

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (nx == 1L && length(d) > 1L) x <- rep(x, length(d))

    ## Guessing elementwise if set NULL
    if (is.null(elementwise)) elementwise <- !length(x) == length(d)

    if (!elementwise & length(x) > length(d))
        stop("length of 'x' can't be larger than number of distributions in 'd'",
             "if elementwise = FALSE");

    if (elementwise) x <- sort(x) # Important

    ui <- seq_along(d) # Unique index
    d  <- as.data.frame(d) # convert distributions to data.frame

    ## Calling C to calculate the required quantiles
    ## binmid: Required as we need to know where to stop.
    ## y: threshold at which to evaluate the pdf.
    print(x)
    res <- .Call("tm_predict", uidx = ui, idx = d$index, tp = d$tp, binmid = d$binmid, y = x,
                 type = "pdf", ncores = ncores, elementwise = elementwise)
    print(res)
    print(head(d))
    stop('xxx')

    # Translate quantile bins to numeric values (binmid)
    bin2num <- function(x, d, np, ewise) {
        idx <- rep(ui, each = if (ewise) np else 1L)
        # 'order' is used to recreate the correct order before return
        x   <- data.frame(order = seq_along(idx), index = idx, bin = res)
        res <- merge(d, x, by = c("bin", "index"), all.x = FALSE, all.y = TRUE)
        return(res[order(res$order), "binmid", drop = TRUE])
    }
    res <- bin2num(res, x, nprobs, elementwise)

    # If elementwise: Return matrix of dimension length(d) x length(probs)
    if (elementwise) {
        res <- matrix(res, ncol = length(probs),
                      dimnames = list(xnames, paste("q", format(probs, digits = 3), sep = "_")))
        return(as.data.frame(res))
    # Else return named vector
    } else {
        return(setNames(res, xnames))
    }

#####    # Scopes 'd' and 'grd'
#####    fn <- function(i) {
#####        grd <- as.list(grd[i, ])
#####        tmp <- as.data.frame(d[grd$d])
#####        cutoff <- tmp$bin[which.min(abs(tmp$binmid - grd$p))]
#####        tmp$index <- grd$index
#####        tmp[tmp$bin <= cutoff, ]
#####    }
#####
#####    # If there is one distribution and one point at which to evaluate this
#####    # distribution, or 'elementwise' is FALSE we can proceed setting up the
#####    # data.frame used to call the C function performing the calculations.
#####    # Note: If length(d) != length(x) the vectors will be recycled.
#####    if ((length(d) == 1L && length(x) == 1L) || isFALSE(elementwise)) {
#####
#####        # Create 'grid' for evaluation
#####        if (length(d) > length(x)) {
#####            grd <- data.frame(d = seq_along(d), p = rep(x, length(d)))
#####        } else {
#####            grd <- data.frame(d = rep(seq_along(d), length(x)), p = x)
#####        }
#####        print(grd)
#####        grd$index <- seq_len(nrow(grd)) # Pseudo-index
#####
#####        ui <- seq_len(nrow(grd))
#####        d <- do.call(rbind, lapply(ui, fn))
#####    # If elementwise we build the data.frame differently, simply
#####    # Keeping the full distribution.
#####    } else if (elementwise) {
#####        ui <- seq_along(d) # Number of distributions
#####        d  <- as.data.frame(d)
#####    } else {
#####        stop(" --- TODO(R): Can we end up here, and if so, why? Work needed --- ")
#####    }
#####    print(d)
#####    stop(3)
#####
#####    ## Calling C to calculate the required PDFs
#####    res <- .Call("tm_predict", ui, d$index, p = d$tp, type = "pdf", prob = 42,
#####                 ncores = ncores, elementwise = elementwise)
#####
#####    if (elementwise) {
#####        res <- matrix(res, ncol = length(x),
#####                      dimnames = list(NULL, paste("x", format(x, digits = 3), sep = "_")))
#####    }
#####    return(res)

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
    ## Get number of cores for OpenMP parallelization
    ncores <- tm_get_number_of_cores(ncores, FALSE)

    # Guessing elementwise if set NULL
    if (is.null(elementwise)) elementwise <- !length(x) == length(x)

    xnames <- names(x) # Keep for later

    # Number of probabilities
    nprobs <- length(probs)

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (!elementwise & nprobs == 1L && length(x) > 1L) probs <- rep(probs, length(x))

    if (!elementwise & length(probs) > length(x))
        stop("length of 'probs' can't be larger than number of distributions in 'd'",
             "if elementwise = FALSE");

    if (elementwise) probs <- sort(probs) # Important

    ui <- seq_along(x) # Unique index
    x  <- as.data.frame(x) # convert distributions to data.frame

    ## Calling C to calculate the required quantiles
    ## binmid: Not required
    ## y: probabilities at which to evaluate the distributions
    res <- .Call("tm_predict", uidx = ui, index = x$index, tp = x$tp, binmid = x$binmid, y = probs,
                 type = "quantile", ncores = ncores, elementwise = elementwise)

    # If elementwise: Return matrix of dimension length(d) x length(probs)
    if (elementwise) {
        res <- matrix(res, ncol = length(probs), byrow = TRUE,
                      dimnames = list(xnames, paste0(format(probs * 1e2, digits = 3), "%")))
        return(as.data.frame(res))
    # Else return named vector
    } else {
        return(setNames(res, xnames))
    }
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
