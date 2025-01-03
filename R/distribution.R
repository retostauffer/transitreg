

#' Creates a tm Distribution
#'
#' A tm distrubiton consists of a series of \code{K} transition
#' probabilities for \code{K} 'bins' (counts or pseudo-counts) and
#' a series of \code{K} numeric values representing the (center of) the
#' corresponding bins.
#'
#' @param probs numeric vector with transition probabilities (0, ..., K).
#' @param bins numeric vector with numeric value for bin 0, ..., K.
#'
#' @return Returns an object of class \code{c("tmdist", "distribution")}.
#' @author Reto
tmdist <- function(probs, bins) {
    stopifnot(is.numeric(probs), is.numeric(bins), length(probs) == length(bins))
    n <- length(probs)
    x <- c(setNames(as.list(probs), sprintf("tp_%d", seq_len(n) - 1L)),
           setNames(as.list(bins),  sprintf("bin_%d", seq_len(n) - 1)))
    x <- as.data.frame(x)
    class(x) <- c("tmdist", "distribution")
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


#' tm Quantiles
#'
#' @param x object of class \code{"tmdist"}.
#' @param probs numeric, numeric vector of probabilities with values in [0,1].
#' @param drop logical TODO(R) explain.
#' @param elementwise TODO(R) not yet implemented as we force length(x) == 1
#' @param .. TODO(R) check and extend if needed.
#'
#' At the end there will be one man page for tmdist.
quantile.tmdist <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
    stopifnot("currently only designed length 1" = length(x) == 1L)
    # Convert input to data.frame; expand data.frame (make a copy
    # of the data.frame for each element in 'probs' such that we can
    # call tm_predict only once).
    ui <- seq_along(probs) # 'unique index'
    x  <- as.data.frame(x)
    x  <- lapply(ui, function(i, x) { x$index <- i; return(x) }, x = x)
    x  <- do.call(rbind, x)
    # TODO(R): Currently no OpenMP (ncores = 1L)
    warning("elementwise = NULL equals boolean elementwise = 1 on the C end; needs adaption");
    .Call("tm_predict", ui, x$index, p = x$tp, type = "quantile", prob = probs, 1L, elementwise = elementwise)
}

pdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
    stopifnot("currently only designed length 1" = length(d) == 1L)

    # TODO(R): Hack
    if (is.null(elementwise)) {
        elementwise <- FALSE
        warning("pdf.tmdist; setting elementwise FALSE if is NULL, requires some work")
    }

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (length(x) == 1 && length(d) > 1L) x <- rep(x, length(d))

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
    # TODO(R): Currently no OpenMP (ncores = 1L)
    res <- .Call("tm_predict", ui, d$index, p = d$tp, type = "pdf", prob = 42,
                 ncores = 1L, elementwise = elementwise)

    if (elementwise) {
        res <- matrix(res, ncol = length(x),
                      dimnames = list(NULL, paste("x", format(x, digits = 3), sep = "_")))
    }
    return(res)

  #FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), ...)
  #apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


cdf.tmdist <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
    stopifnot("currently only designed length 1" = length(d) == 1L)

    # TODO(R): Hack
    if (is.null(elementwise)) {
        elementwise <- FALSE
        warning("cdf.tmdist; setting elementwise FALSE if is NULL, requires some work")
    }

    # Setting up all unique combinations needed
    # (1) Use same 'x' for all distributions
    if (length(x) == 1 && length(d) > 1L) x <- rep(x, length(d))

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
    # TODO(R): Currently no OpenMP (ncores = 1L)
    res <- .Call("tm_predict", ui, d$index, p = d$tp, type = "cdf", prob = 42,
                 ncores = 1L, elementwise = elementwise)

    if (elementwise) {
        res <- matrix(res, ncol = length(x),
                      dimnames = list(NULL, paste("x", format(x, digits = 3), sep = "_")))
    }
    return(res)

  #FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), ...)
  #apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


