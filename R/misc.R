grep2 <- function(pattern, x, ...) {
  i <- NULL
  for (p in pattern)
    i <- c(i, grep(p, x, ...))
  unique(i)
}

#' @importFrom stats formula
#' @importFrom Formula as.Formula
response_name <- function(formula) {
  if (is.list(formula)) {
    if (!is.null(formula$formula)) {
      formula <- formula$formula
    } else {
      formula <- do.call("as.Formula", formula)
    }
  }
  formula <- as.Formula(formula)
  formula <- formula(formula, rhs = 0, collapse = TRUE)
  rn <- all.vars(formula)
  rn <- rn[rn != "."]
  return(rn)
}

## Function takes a formula, Formula, or a list of formulas
## and extracts only the necessary parts to create a model.frame,
## as well as the parts that are needed to setup smooth term
## specification lists, or any type of special model terms
## that should be used for fitting.

#' @importFrom stats formula terms update
#' @importFrom Formula as.Formula
fake_formula <- function(formula, specials = NULL, nospecials = FALSE, onlyspecials = FALSE) {
  if (is.list(formula))
    formula <- do.call("as.Formula", formula)

  if (any(grepl("|", as.character(formula), fixed = TRUE)))
    formula <- as.Formula(formula)

  n <- length(formula)

  outer <- FALSE
  if (length(n) > 1) {
    if (n[2L] > 1)
      outer <- TRUE
  }

  if (outer) {
    fl <- list()
    lhs <- formula(formula, rhs = 0)
    for (i in 1:n[2L]) {
      fl[[i]] <- fake_formula(formula(formula, rhs = i, lhs = 0), specials = specials,
        nospecials = nospecials, onlyspecials = onlyspecials)
      if (!is.character(fl[[i]]))
        fl[[i]] <- formula(as.Formula(fl[[i]]), drop = TRUE, collapse = TRUE)
    }
    if (length(fl)) {
      if (!is.character(fl[[1L]])) {
        formula <- . ~ .
        fl <- do.call("as.Formula", fl)
        fl <- formula(fl, lhs = 0, drop = FALSE)
        formula <- as.Formula(fl)
        formula[[3L]] <- fl
        attr(formula, "lhs") <- attr(as.Formula(lhs), "lhs")
      } else {
        formula <- fl
      }
    }
  } else {
    stn <- c("s", "te", "t2", "sx", "s2", "rs", "ti",
      "tx", "tx2", "tx3", "tx4", "la", "lasso", "n", "lin",
      "pb", "pbc", "nn", "fk", "re", "ps", "pbz", "ga",
      "random", "ra", "lo", "tr", "tree", "cf", "NN", "pb2", "ct",
      "st", "ps2", "pdDiag", "user")
    stn <- unique(c(stn, specials))
    formula <- ff_replace(formula)

    mt <- terms(formula, specials = stn)

    os <- NULL
    st <- unlist(attr(mt, "specials"))
    if (!is.null(st)) {
      tls <- rownames(attr(mt, "factors"))[st]
      tls <- unique(tls)
      ff <- NULL
      for (j in tls) {
        p <- parse(text = j)
        if (as.character(p[[1]][[1]]) %in% c("la", "lasso")) {
          p <- p[[1]][1:2]
        }
        v <- all.vars(p)
        e <- eval(parse(text = paste0("quote(", j, ")")))
        for (i in seq_along(e)) {
          av <- all.vars(e[[i]])
          av <- av[av != "|"]
          if (length(av)) {
            if (all(av %in% v)) {
              ef <- try(eval(e[[i]]), silent = TRUE)
              if (!inherits(ef, "try-error")) {
                if (inherits(ef, "formula")) {
                  vf <- attr(terms(eval(ef)), "variables")
                  for (k in 2:length(vf)) {
                    if (as.character(e[1L]) == "lin") {
                      lv <- all.vars(vf[[k]])
                      for (l in lv)
                        ff <- c(ff, eval(parse(text = paste0("quote(", l, ")"))))
                    } else {
                      ff <- c(ff, vf[[k]])
                    }
                  }
                  next
                }
              }
              if (as.character(e[[i]])[1L] == "~") {
                vf <- attr(terms(eval(e[[i]])), "variables")
                for (k in 2:length(vf))
                  ff <- c(ff, vf[[k]])
              } else {
                if (as.character(e[1L]) %in% c("lin")) {
                  lv <- all.vars(e[[i]])
                  for (l in lv)
                    ff <- c(ff, eval(parse(text = paste0("quote(", l, ")"))))
                } else {
                  ff <- c(ff, e[[i]])
                }
              }
            }
          }
        }
        os <- c(os, j)
        eval(parse(text = paste0("formula <- update(formula, NULL ~ . -", j,")")))
        if (!nospecials) {
          for (i in ff) {
            eval(parse(text = paste0("formula <- update(formula, NULL ~ . +", deparse(i),")")))
          }
        }
      }
    }

    if (onlyspecials) {
      if (!is.character(os))
        os <- character(0)
      os <- gsub(" ", "", os)
      formula <- unique(os)
    } else {
      tf <- terms(formula, specials = stn)
      sj <- unlist(attr(tf, "specials"))
      if (!is.null(sj)) {
        formula <- fake_formula(formula)
      }
      tf <- terms(formula, specials = stn)
      if (length(j <- grep("list(", attr(tf, "term.labels"), fixed = TRUE, value = TRUE))) {
        fc <- paste("formula <- update(formula, . ~ . -", j, ")")
        eval(parse(text = fc))
      }
    }
  }

  return(formula)
}

## Replace * and : with +.
ff_replace <- function(formula) {
  n <- length(formula)
  if (length(n) > 1) {
    f <- formula[[max(n)]]
  } else {
    f <- formula[[n]]
  }
  f <- deparse(f)
  if (any(grepl(":", f, fixed = TRUE))) {
    f <- gsub(":", "+", f, fixed = TRUE)
    f <- as.call(parse(text = f))
    formula[[n]] <- f[[1L]]
    formula <- update(formula, . ~ .)
  }
  return(formula)
}


# -------------------------------------------------------------------
# Utility functions to convert between CDF, PDF, and Transition
# Probabilities (TP). Mainly used internally for testing as well
# as for demonstration purposes; not part of the rest of the software.
# -------------------------------------------------------------------

#' Converting between Densities, Probabilities, and Transition Probabilities
#'
#' A transition model consists of a series of transition probabilities (TP)
#' which define the probability that an observation falls into the next higher
#' bin. Based on these transition probabilities, the cumulative distribution
#' function (CDF) as well as the density (PDF) can be calculated.
#' This utility function allows to convert from CDF to TPs, and from TPs
#' back to CDF/PDF.
#'
#' @param x Numeric vector with density, probabilities, or transition probabilities.
#'        All non-missing `>=0`.
#' @param from Character, current type of data in `x`.
#'        Either `"tp"`, `"cdf"`, or `"pdf"`.
#' @param to  Character, type into which `x` should be converted (see `from`).
#'        Note that not all combinations of conversions are allowed. If not possible,
#'        an error will be thrown.
#' @param width `NULL` or numeric (either length `1` or same length as `x`.
#'        Width of the individual bins represented in `x`.
#' @param drop If `TRUE` (default) the result is simplified if possible.
#'
#' @details
#' If `from = "tp"` (transition probabilities) all values in `x` must be
#' in `[0, 1]` and monotonically decreasing. Transition probabilities
#' can be converted into both, `to = "cdf"` as well as `to = "pdf"` (or both
#' at the same time, setting `to = c("cdf", "pdf"`)`.
#'
#' If `from = "cdf"` (distribution) all elements in `x` must be in `[0,1]`.
#' Can be converted `to = "tp"`.
#'
#' Similarly, `from = "pdf"` can only be converted `to = "tp"`. This requires
#' all elements in `x` to be in `[0, Inf]`.
#'
#' @return If `to` is a single character and `drop = TRUE`, the return is a numeric vector
#' of the same length as the input argument `x`. Else the return is
#' a `data.frame` with the same number of rows as `length(x)`.
#'
#' @examples
#' ## For testing
#' ## Draw PDF and CDF from Poisson distribution
#' p <- ppois(0:10, lambda = 4)
#' d <- dpois(0:10, lambda = 4)
#'
#' ## Convert Poisson CDF to transition probabilities
#' tp <- convert_tp(p, from = "cdf", to = "tp")
#' round(tp, 2)
#'
#' pd <- convert_tp(tp, from = "tp", to = c("pdf", "cdf"))
#'
#' ## Convert transition probabilities back to CDF
#' p2 <- convert_tp(tp, from = "tp", to = "cdf")
#' all.equal(p, p2) # Must be equal
#'
#' ## Convert transition probabilities to PDF
#' d2 <- convert_tp(tp, from = "tp", to = "pdf")
#' all.equal(d, d2) # Must be equal
#'
#' ## Or directly from tp to both, CDF and PDF
#' d2 <- convert_tp(tp, from = "tp", to = c("cdf", "pdf"))
#' all.equal(p, d2$cdf)
#' all.equal(d, d2$pdf)
#'
#' ## Quick visual representation
#' barplot(tp, col = "steelblue", main = "Transition Probabilities")
#' col <- c("gray80", "tomato")
#' barplot(rbind(p, p2), beside = TRUE, main = "CDF Comparison", col = col)
#' barplot(rbind(d, d2$pdf), beside = TRUE, main = "PDF Comparison", col = col)
#'
#' @author Reto
#' @export
convert_tp <- function(x, from, to, width = NULL, drop = TRUE) {
    stopifnot(
        "'x' must be numeric of length > 0" = is.numeric(x) && length(x) > 0,
        "'x' does not allow for NAs" = all(!is.na(x)),
        "values in 'x' must be between 0.0 and 1.0" = all(x >= 0)
    )
    from <- tolower(from); to <- tolower(to)
    from <- match.arg(from, c("cdf", "pdf", "tp"))
    to   <- match.arg(to,   c("cdf", "pdf", "tp"), several.ok = TRUE)

    # Checking range
    if (from %in% c("tp", "cdf") & !all(x <= 1.0))
        stop("if from = \"", from, "\" all 'x' must be in [0, 1]")

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
        res[[do[[i]]]] <- get(fnname)(x, width = width)
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
# @param x numeric vector with CDFs
# @param \dots unused.
cdf_to_tp <- function(x, ...) {
    stopifnot(is.numeric(x), all(x >= 0 & x <= 1))
    stopifnot(all(diff(x) >= -.Machine$double.eps))

    # Converting CDF to transition probabilities
    tp <- numeric(length(x))
    tp[1] <- 1 - x[1]
    prod  <- 1
    for (i in seq.int(2, length(x))) {
        prod <- prod * tp[i - 1]
        if (prod < sqrt(.Machine$double.eps)) {
            tmp   <- 0.0
        } else {
            tp[i] <- (x[i - 1] - x[i] + prod) / prod
        }
    }
    return(tp)
}

# Converts transition probabilities (tp) to cdf
# @param x numeric vector with TPs.
# @param \dots unused.
tp_to_cdf <- function(tp, ...) {
    stopifnot(is.numeric(tp), all(tp >= 0 & tp <= 1), length(tp) >= 1L)

    res <- (1 - tp[1])
    if (length(tp) > 1) {
        for (i in seq.int(2, length(tp))) {
            res <- c(res, res[i - 1] + (1 - tp[i]) * prod(tp[1:(i-1)]))
        }
    }
    return(res)
}

# Converts transition probabilities (tp) to pdf
# @param x numeric vector with TPs.
# @param width `NULL` or a numeric vector of the same length as `x`.
#        Width of the bins, required to properly scale the CDF if
#        the width of the bins is not equal to `1` (as it is for count data).
# @param \dots unused.
tp_to_pdf <- function(tp, width = NULL, ...) {
    stopifnot(is.numeric(tp), all(tp >= 0 & tp <= 1), length(tp) >= 1L)

    if (!is.null(width)) {
        stopifnot(
            "'width' (if provided) must be numeric" = is.numeric(width),
            "'width' must be of the same length as 'tp'" = length(width) == length(tp),
            "'width' must be > 0" = all(width > 0)
        )
    }

    res <- numeric()
    p <- 1.0
    for (i in seq_along(tp)) {
        res <- c(res, p * (1.0 - tp[i]))
        p   <- p * tp[i]
    }
    return(if (is.null(width)) res else res / width)
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
#' @param digits Integer, number of significant digits for names.
#'
#' @details If `prefix = NULL` it is expected that `x` contains probabilities
#' in the range of `[0,1]`. The returned vector of names will therefore contain
#' `"10%"`, `"20%"`, `"99%`" etc.
#'
#' When a `prefix` is set, `x` is interpreted as numeric value and the returned
#' names will be a combination of the `prefix` and the numeric value of `x`.
#' This will result in e.g, `p_-12.4`, `p_0`, `p_11.302`, ... .
#'
#' @author Reto
#' @return A character vector with 'column names'.
get_elementwise_colnames <- function(x, prefix = NULL, digits = pmax(3L, getOption("digits") - 3L)) {
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
    return(trimws(x))
}



#' Convert Numerics to Pseudo-Bins
#'
#' Convert response on original scale to (pseudo-)bin index.
#'
#' @param x numeric vector.
#' @param breaks numeric vector with point intersections to create the bins.
#' @param censored character, defines if the (pseudo-)bins are censored or not
#'        (see 'Details' section).
#' @param verbose logical, if `TRUE` some messages are printed if needed
#'        (defaults to `FALSE`).
#'
#' @details
#' Converts the numeric values in `x` into (pseudo-)bins by cutting the data
#' into `length(breaks) - 1` segments. If `censored = "uncensored"` (no censoring),
#' data outside `range(break)` will be set to `-1L`
#' if `x < min(breaks)` or `length(breaks) - 1L` if `x > max(breaks)`.
#'
#' If `censored = "left"` all values `x <= min(breaks)` are assigned to as bin `0L`,
#' followed by bin `1L` for `(breaks[1], breaks[2])`, `2L` for `[breaks[2], breaks[3])` etc.
#' Thus, the highest bin containing data is `length(breaks)` as well.
#'
#' If `censored = "right"` all values `x >= max(breaks)` are set to the highest bin
#' `length(breaks)`, those `< max(breaks)` are still set `NA_integer_`.
#'
#' When `censored = "both"` the rules above are combined. Due to the necessity that
#' we need tone dedicated 'censored bin index' on the left (`0L`), the right most
#' upper bin is shifted by `+1` compared to when we only use `censored = "right"`.
#'
#' Less important to know for an end-user as this is mainly handelled internally.
#'
#' @examples
#' x      <- 0:10
#' breaks <- c(2, 6, 8)
#' transitreg:::num2bin(x, breaks = breaks)
#' transitreg:::num2bin(x, breaks = breaks, censored = "left")
#' transitreg:::num2bin(x, breaks = breaks, censored = "right")
#' transitreg:::num2bin(x, breaks = breaks, censored = "both")
#' transitreg:::num2bin(x, breaks = breaks, censored = "both")
#'
num2bin <- function(x, breaks = NULL, censored = c("uncensored", "left", "right", "both"),
                    verbose = FALSE) {
    if (is.null(x)) return(x)
    if (is.null(breaks))
        stop("Either breaks (continuous response) or bins (discrete response) must be set")
    if (any(is.na(x)))
        stop("Response contains missing values (not allowed).")
    censored <- match.arg(censored)

    # TODO(R): Add test cases where breaks are not unique at the lower
    #          and upper end; used for censoring.
    cens_left  <- censored == "left"  || censored == "both"
    cens_right <- censored == "right" || censored == "both"

    # Censored? ignore duplicated breaks on the left or right end if needed.
    if (cens_left)  breaks <- c(min(breaks), breaks[!breaks == min(breaks)])
    if (cens_right) breaks <- c(max(breaks), breaks[!breaks == max(breaks)])
    breaks <- sort(breaks) # Sorted breaks without duplicates on left/right (if possible)

    # 'Cut' data, limit to -1 to length(brekas) - 1.
    res <- cut(x, breaks = breaks, labels = FALSE, right = FALSE, include.lowest = TRUE) - 1L

    res[x < min(breaks)] <- -1L                 ## Outside
    res[x > max(breaks)] <- length(breaks) - 1L ## Outside

    # Adjusting the right side of the distribution in case censoring
    # has been set (MUST be done before 'left').
    if (cens_right) {
        res[x >= max(breaks)] <- length(breaks) - 1L
        # The rest is just for verbose output
        n <- if (verbose) sum(x > max(breaks)) else -999
        if (n > 0) message(sprintf("num2bin: %d observation%s been set ", n, ifelse(n == 1, " has", "s have")),
                           sprintf("to the right censoring point (bin %dL ", length(breaks) - 1L),
                           sprintf("at %s; censored = \"%s\").", format(min(breaks)), censored))
    }

    # Adjusting the left side of the distribution in case censoring
    # has been set (MUST be done after checking the right side).
    if (cens_left) {
        res <- res + 1L; res[x <= min(breaks)] <- 0L
        # The rest is just for verbose output
        n <- if (verbose) sum(x < min(breaks)) else -999
        if (n > 0) message(sprintf("num2bin: %d observation%s been set ", n, ifelse(n == 1, " has", "s have")),
                           "to the left censoring point (bin 0L ",
                           sprintf("at %s; censored = \"%s\").", format(min(breaks)), censored))
    }

    return(res)
}

# Converts bin index (idx) to numeric values, using the center
# of the bin (defined by breaks) the index falls into. For
# indices equal to -1, or indices > max bins the return will
# be set to the minimum/maximum break (i.e., the edge of the
# defined range).
# TODO(R): Write some tests for this
bin2num <- function(idx, breaks, censored = c("uncensored", "left", "right", "both")) {
    censored <- match.arg(censored)

    # TODO(R): Add test cases where breaks are not unique at the lower
    #          and upper end; used for censoring.
    cens_left  <- censored == "left"  || censored == "both"
    cens_right <- censored == "right" || censored == "both"

    # Censored? Create boule-breaks on left/right if needed
    if (cens_left)  breaks <- c(rep(min(breaks), 2L), breaks[!breaks == min(breaks)])
    if (cens_right) breaks <- c(rep(max(breaks), 2L), breaks[!breaks == max(breaks)])
    breaks <- sort(breaks) # Sorted breaks without duplicates on left/right (if possible)

    # Bin mid, draw bin-mid for each index to convert to numeric
    bm <- (head(breaks, -1) + tail(breaks, -1)) / 2
    return(bm[idx + 1L])
}


#' Create Breaks for (Pseudo-)bins
#'
#' Calculates the breaks (point intersection between breaks)
#' which span the range of `y` if the user did not specify 
#' a vector of breaks, but only a single numeric (integer).
#'
#' @param y a numeric vector with response data.
#' @param breaks a single numeric (integer) >= 1, number of breaks to be created.
#' @param censored character, one of `"uncensored"`, `"left"`, `"right"`.
#'
#' @return Returns a numeric vector with the breaks.
make_breaks <- function(y, breaks, censored) {
    stopifnot(
        "'y' must be numeric length > 0L" = is.numeric(y) && length(y) > 0L,
        "'breaks' must be 'NULL' or numeric vector of length 1L or length > 2L" =
            is.null(breaks) || (is.numeric(breaks) && length(breaks) == 1L) || (is.numeric(breaks) && length(breaks) > 2L)
    )
    cens_allowed <- c("uncensored", "left", "right", "both")
    if (is.character(censored) && length(censored) == 1L) {
        censored <- match.arg(censored, cens_allowed)
    } else if (is.numeric(censored) && length(censored) %in% 1:2) {
        # All fine
    } else {
        stop("'censored' must be or a numeric vector of length 1 or 2 or one of ",
             paste(cens_allowed, collapse = ", "))
    }

    # If 'censored' is numeric and length(breaks) > 1L, error.
    if (is.numeric(censored) && (is.numeric(breaks) && length(breaks) > 1L))
        stop("if 'censored' is numeric, breaks must be NULL or a single numeric value")


    # Check if all breaks are positive integers, as well as all response values.
    # This is used identify count data models.
    eps   <- sqrt(.Machine$double.eps)
    y_int <- all(abs(y - round(y)) <= eps) && all(round(y) >= 0L)

    # If response does not look like count data, 'breaks' must be set.
    if (!y_int && is.null(breaks))
        stop("response does not look like count data, in this case 'breaks' must be specified")

    # Number of distinct values in y, and unique y values
    if (y_int) ny <- max(y <- unique(round(y))) else ny <- max(y <- unique(y))

    # -------------------------------
    # (1) If all values in 'y' look like count data, and no breaks are specified,
    #     we expect this to be a count data model call.
    if (y_int && is.null(breaks)) {
        # Censored count data models are not allowed (for now)
        if (!isTRUE(censored == "uncensored"))
            stop("response looks like count data and no 'breaks' were specified. Censoring not allowed.")

        res <- NULL
        nbins <- round(max(y)) * if (ny <= 10) { 3 } else if (ny <= 100) { 1.5 } else { 1.25 }

    } else if (y_int && (is.numeric(breaks) && length(breaks) == 1L) && isTRUE(censored == "uncensored")) {
        if (breaks < max(y))
            stop("response looks like count data, breaks must be >= max(response)")

        res      <- NULL
        censored <- censored #"uncensored"
        nbins    <- breaks - 1L

    # -------------------------------
    # (2) If 'breaks' is numeric of length 1, the user just specified how many
    #     breaks they want. Automatically generate them in the range of y.
    } else if (is.numeric(breaks) && length(breaks) == 1L) {
        stopifnot("if 'breaks' is a single numeric, it must be >= 3" = breaks >= 3)
        breaks <- as.integer(breaks)
        dy     <- diff(range(y))

        if (is.character(censored)) {
            res     <- seq(min(y) - 0.1 * dy, max(y) + 0.1 * dy, length = breaks)
        } else if (length(censored) == 1L || (length(censored) == 2L && is.na(censored[[2L]]))) {
            if (censored[1L] > max(y)) stop("censoring point > max(y)")
            res     <- seq(censored[[1]], max(y) + 0.1 * dy, length = breaks)
            censored <- "left"
        } else if (length(censored) == 2L && is.na(censored[[1L]])) {
            if (censored[2L] < min(y)) stop("censoring point < min(y)")
            res     <- seq(min(y) - 0.1 * dy, censored[[2L]], length = breaks)
            censored <- "right"
        } else {
            if (censored[1L] >= censored[2]) stop("invalid censoring points provided")
            if (censored[2L] < min(y) || censored[1] > max(y)) stop("censoring outside range of data")
            res     <- seq(censored[[1L]], censored[[2L]], length = breaks)
            censored <- "both"
        }

    # -------------------------------
    # (3) If 'breaks' is a numeric vector, check if these breaks make sense
    } else if (is.numeric(breaks) && length(breaks) > 1L) {
        res <- as.numeric(breaks)

        # Breaks not spanning the range?
        if ((censored == "uncensored" || censored == "right") && min(y) < min(res))
            stop("breaks not spanning the range of the response (min(response) < min(breaks))")
        if ((censored == "uncensored" || censored == "left") && max(y) > max(res))
            stop("breaks not spanning the range of the response (max(response) > max(breaks))")

    # -------------------------------
    # (X) Uncatched case
    } else {
        tmp <- list(y = y, ny = ny, y_int = y_int, censored = censored, breaks = breaks)
        str(tmp)
        stop("whoops, case not catched yet")
    }

    # Censored left or right?
    if (!is.null(res) && (censored == "left"  || censored == "both")) res <- c(min(res), res)
    if (!is.null(res) && (censored == "right" || censored == "both")) res <- c(res, max(res))

    # Number of bins if breaks were specified
    if (!is.null(res)) nbins <- length(res) - 1L

    return(list(breaks = res, bins = as.integer(nbins), censored = censored))
}


# Helper functions to check the arguments we hand over to C. This
# helps identifying potential problems easier. Input 'x' is a named
# list with all the elements required when calling "treg_predict" in C.
#
# The function ensures that (i) all the objects are of the correct type,
# and that some elements do show the correct (expected) length. Else
# it is possible that C returns 'garbage' results as it tries to read
# elemnets outside of memory, resulting in either segfaults or werid
# results.
#
# Last but not least it ensures that the order of the arguments
# are as expected by .C!

#' @importFrom utils str
check_args_for_treg_predict <- function(x, silent = FALSE) {
    ## Expected elements (in the order expected by C)
    enames <- c("uidx", "idx", "tp", "breaks", "censored", "y", "prob",
                "type", "ncores", "elementwise", "discrete")

    ## Checking types first
    tmp <- list("integer"   = c("uidx", "idx", "y", "ncores"),
                "double"    = c("tp", "breaks", "prob"),
                "logical"   = c("elementwise", "discrete"),
                "character" = c("type", "censored"))

    debug_stop <- function(e) { if (!silent) { cat("\nDebugging output (str(args)):\n"); str(x) }; stop(e) }
    for (n in names(tmp)) {
        fn <- get(sprintf("is.%s", n))
        for (e in tmp[[n]]) {
            if (!fn(x[[e]])) {
                debug_stop(paste0("Element '", e, "' in args list is not \"", n, "\""))
            }
        }
    }

    ## Checking length of some of the elements
    tryCatch(
        {stopifnot(
            "'args$elementwise' must be of length 1" = length(x$elementwise) == 1L,
            "length of 'args$idx' and 'args$tp' must be identical" = 
                length(x$idx) == length(x$tp),
            "length of 'args$type' must be 1" = length(x$type) == 1L,
            "length of 'args$discrete' and 'args$uidx' must be identical" =
                length(x$uidx) == length(x$discrete),
            "not all required elements found in 'args'" =
                all(enames %in% names(x))
        )},
        error = function(e) debug_stop(e)
    )

    # Return (re-ordered) list
    return(x[enames])
}

#' @importFrom utils str
check_args_for_treg_predict_pdfcdf <- function(x, silent = FALSE) {
    ## Required elements in the order as expected by C
    enames <- c("uidx", "idx", "tp", "y", "breaks", "censored", "discrete", "ncores")

    ## Checking types first
    tmp <- list("integer"   = c("uidx", "idx", "y", "ncores"),
                "double"    = c("tp", "breaks"),
                "character" = "censored",
                "logical"   = "discrete") # TODO(R) Cover this in tests

    debug_stop <- function(e) { if (!silent) { cat("\nDebugging output (str(args)):\n"); str(x) }; stop(e) }
    for (n in names(tmp)) {
        fn <- get(sprintf("is.%s", n))
        for (e in tmp[[n]]) {
            if (!fn(x[[e]])) {
                debug_stop(paste0("Element '", e, "' in args list is not \"", n, "\""))
            }
        }
    }

    ## Checking length of some of the elements
    tryCatch(
        {stopifnot(
            "length of 'args$idx' and 'args$tp' must be identical" =
                length(x$idx) == length(x$tp),
            "length of 'args$y' and 'args$uidx' must be identical" =
                length(x$y) == length(x$uidx),
            "'args$discrete' must be logical with same length as args$uidx" = 
                length(x$discrete) == length(x$uidx),
            "'args$censored' must be one of \"uncensored\", \"left\", \"right\", \"both\"" =
                x$censored %in% c("uncensored", "left", "right", "both")
        )},
        error = function(e) debug_stop(e)
    )

    # Return (re-ordered) list
    return(x[enames])
}



#' Detect number of cores for OpenMP
#'
#' The calculation of CDFs and PDFs is implemented in C and allows
#' for parallelization using OpenMP. This function detects how may
#' cores are available in total (if OpenMP is available).
#'
#' @param verbose logical, if `TRUE` a message is shown.
#'
#' @return Number of available cores (integer). If OpenMP is not
#' available, `1L` is returned.
#'
#' @author Reto
transitreg_detect_cores <- function(verbose = TRUE) {
    stopifnot("'verbose' must be logical TRUE or FALSE" = isTRUE(verbose) || isFALSE(verbose))
    ncores <- .Call("treg_detect_cores")
    if (verbose && ncores) {
        message("OMP available, number of cores detected: ", ncores)
    } else if (verbose) {
        message("OMP not available (not compiled with omp)")
    }
    return(ncores)
}

#' Get number of cores for OpenMP
#'
#' Some parts of the package use C routines which allow for parallelization
#' using OpenMP. This function is used to specify how many cores to be used.
#'
#' @param ncores `NULL` or a positive integer.
#' @param verbose logical, if `TRUE` a message is shown.
#'
#' @return Number of cores to be used in OpenMP parallelization (integer).
#'
#' @details If `ncores` is `NULL` the number of available
#' cores is auto-detected and set to 'total number of cores - 2'.
#' If integer, it is checked if this number of cores is available,
#' else set tot he 'total number of cores available'.
#'
#' @author Reto
transitreg_get_number_of_cores <- function(ncores = NULL, verbose = FALSE) {
  ## Number of cores to be used for OpenMP. If
  ## - NULL: Guess cores (max cores - 2L)
  ## - Smaller or equal to 0: Set to 1L (single-core processing)
  ## - Else: Take user input; limited to maximum number of detected cores.
  ncores <- if (!is.null(ncores)) as.integer(ncores)[1L] else transitreg_detect_cores(verbose = FALSE) - 2L
  ncores <- if (ncores < 1L) 1L else pmin(ncores, transitreg_detect_cores(verbose = FALSE))
  if (verbose) message("Number of cores set to: ", ncores)
  return(ncores)
}


