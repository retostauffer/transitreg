


# Allows to convert from transition probabilities to
# PDF/CDF and vice versa (not all combinations are allowed).
tm_convert <- function(x, from, to, width = NULL, drop = TRUE) {
    stopifnot(
        "'x' must be numeric of length > 0" = is.numeric(x) && length(x) > 0,
        "'x' does not allow for NAs" = all(!is.na(x)),
        "values in 'x' must be between 0.0 and 1.0" = all(x >= 0 & x <= 1)
    )
    from <- tolower(from); to <- tolower(to)
    from <- match.arg(from, c("cdf", "pdf", "tp"))
    to   <- match.arg(to,   c("cdf", "pdf", "tp"), several.ok = TRUE)

    # Evaluate input argument 'width'. If width = NULL no weighting
    # is preformed. This is correct if we are talking about count data,
    # or pdf/cdf/tp for continuous data where each bin has an exact
    # width of 1. Else we must perform weighting.
    # If 'width' is provided it must be a numeric vector of either
    # length one (same width for all 'bins') or a numeric vector
    # of the same length as input 'x'.
    if (!is.null(width)) {
        stopifnot(
            "'width' must be numeric" = is.numeric(width),
            "'width' must contain positive values" = all(width > 0),
            "'width' must have length 1 or the same length as 'x'" = 
                length(width) == 1L || length(width) == length(x)
        )
        width <- if (length(width) == 1L) rep(as.numeric(width), length(x)) else as.numeric(width)

        # TODO(R): Not yet implemented!
        stop("TODO(R): The function allows for input 'width', but logic/math not yet implemented")
    }

    # Allowed conversions (name corresponds to 'from', the value as 'to',
    # thus multiple elements in this vector can have the same name.
    allowed <- c("cdf" = "tp", "tp" = "cdf", "tp" = "pdf")

    # Check if the current conversion(s) are allowed
    do <- allowed[names(allowed) == from & unname(allowed) %in% to]
    if (!all(to %in% do)) {
        stop("Requested conversions (from -> to) were:\n",
             "       ", paste(paste(from, to, sep = " -> "), collapse = ", "), "\n",
             "  Not all allowed. Allowed conversions are:\n",
             "       ", paste(paste(names(allowed), allowed, sep = " -> "), collapse = ", "), "\n")
    }

    # Perform conversion
    res <- list()
    for (i in seq_along(do)) {
        fnname <- sprintf("%s_to_%s", names(do)[i], do[[i]])
        res[[do[[i]]]] <- get(fnname)(x)
    }

    # Return
    res <- as.data.frame(res)
    res <- res[, unname(do), drop = drop]

    # Append names (if any)
    if (!is.null(names(x)) && is.data.frame(res)) {
        rownames(res) <- names(x)
    } else if (!is.null(names(x))) {
        names(res) <- names(x)
    }
    return(res)
}

# Converts a cdf vector into transition probabilities.
# The reverse of what tm_predict(..., type = 'pdf') does.
cdf_to_tp <- function(x) {
    stopifnot(is.numeric(x), all(x >= 0 & x <= 1))
    stopifnot(all(diff(x) > 0))

    # Converting CDF to transition probabilities
    tp <- numeric(length(x))
    tp[1] = 1 - x[1]
    prod <- 1
    for (i in seq.int(2, length(x))) {
        prod <- prod * tp[i - 1]
        tp[i] = (x[i - 1] - x[i] + prod) / prod
    }
    return(tp)
}

# Converts transition probabilities (tp) to cdf
tp_to_cdf <- function(tp) {
    stopifnot(is.numeric(tp), all(tp >= 0 & tp <= 1), length(tp) >= 1L)
    stopifnot(all(diff(tp) <= 0))

    res <- (1 - tp[1])
    if (length(tp) > 1) {
        for (i in seq.int(2, length(tp))) {
            res <- c(res, res[i - 1] + (1 - tp[i]) * prod(tp[1:(i-1)]))
        }
    }
    return(res)
}

# Converts transition probabilities (tp) to pdf
tp_to_pdf <- function(tp) {
    stopifnot(is.numeric(tp), all(tp >= 0 & tp <= 1), length(tp) >= 1L)
    stopifnot(all(diff(tp) <= 0))

    res <- numeric()
    p <- 1.0
    for (i in seq_along(tp)) {
        res <- c(res, p * (1.0 - tp[i]))
        p   <- p * tp[i]
    }
    return(res)
}
