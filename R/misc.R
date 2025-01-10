grep2 <- function(pattern, x, ...) 
{
  i <- NULL
  for (p in pattern)
    i <- c(i, grep(p, x, ...))
  unique(i)
}

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
fake_formula <- function(formula, specials = NULL, nospecials = FALSE, onlyspecials = FALSE)
{
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
ff_replace <- function(formula)
{
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

## Histogram and density plot.
plot_hist <- function(x, ...)
{
  x <- na.omit(x)
  h <- hist(x, breaks = "Scott", plot = FALSE)
  d <- density(x)
  ylim <- list(...)$ylim
  if (is.null(ylim))
    ylim <- range(c(h$density, d$y))
  main <- list(...)$main
  if (is.null(main))
    main <- "Histogram and Density"
  xlab <- list(...)$xlab
  if (is.null(xlab))
    xlab <- "Quantile Residuals"
  hist(x, breaks = "Scott", freq = FALSE, ylim = ylim,
    xlab = xlab, main = main, ...)
  lines(d, lwd = 2, col = 4)
  rug(x, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
}

## Q-Q plot.
plot_qq <- function(x, ...)
{
  z <- qnorm(ppoints(length(x)))
  pch <- list(...)$pch
  if (is.null(pch))
    pch <- 19
  col <- list(...)$col
  if (is.null(col))
    col <- rgb(0.1, 0.1, 0.1, alpha = 0.3)
  qqnorm(x, col = col, pch = pch)
  lines(z, z, lwd = 2, col = 4)
}

## Wormplot.
plot_wp <- function(x, ...)
{
  x <- na.omit(x)
  d <- qqnorm(x, plot = FALSE)
  probs <- c(0.25, 0.75)
  y3 <- quantile(x, probs, type = 7, na.rm = TRUE)
  x3 <- qnorm(probs)
  slope <- diff(y3)/diff(x3)
  int <- y3[1L] - slope * x3[1L]
  d$y <- d$y - (int + slope * d$x)

  xlim <- list(...)$xlim
  if (is.null(xlim)) {
    xlim <- range(d$x)
    xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
  }

  ylim <- list(...)$ylim
  if (is.null(ylim)) {
    ylim <- range(d$y)
    ylim <- ylim + c(-0.3, 0.3) * diff(ylim)
  }

  main <- list(...)$main
  if (is.null(main))
    main <- "Worm Plot"
  xlab <- list(...)$xlab
  if (is.null(xlab))
    xlab <- "Theoretical Quantiles"
  ylab <- list(...)$ylab
  if (is.null(ylab))
    ylab <- "Deviation"
  pch <- list(...)$pch
  if (is.null(pch))
    pch <- 19
  col <- list(...)$col
  if (is.null(col))
    col <- rgb(0.1, 0.1, 0.1, alpha = 0.3)

  plot(d$x, d$y, xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, main = main,
    col = col, pch = pch)
  grid(lty = "solid")

  dz <- 0.25
  z <- seq(xlim[1L], xlim[2L], dz)
  p <- pnorm(z)
  se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(d$y)))
  low <- qnorm((1 - 0.95)/2) * se
  high <- -low
  lines(z, low, lty = 2)
  lines(z, high, lty = 2)

  fit <- lm(d$y ~ d$x + I(d$x^2) + I(d$x^3))
  i <- order(d$x)
  lines(d$x[i], fitted(fit)[i], col = 4, lwd = 2)
}

# -------------------------------------------------------------------
# Utility functions to convert between CDF, PDF, and Transition
# Probabilities (TP). Mainly used internally for testing as well
# as for demonstration purposes; not part of the rest of the software.
# -------------------------------------------------------------------

# Allows to convert from transition probabilities to
# PDF/CDF and vice versa (not all combinations are allowed).
convert_tp <- function(x, from, to, width = NULL, drop = TRUE) {
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
