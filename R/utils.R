

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

