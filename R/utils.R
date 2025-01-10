

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
