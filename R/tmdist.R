

#' @author Niki
#' @rdname tmdist
#' @exportS3Method rootogram tmdist
rootogram.tmdist <- function(object, newdata = NULL, plot = TRUE,
  width = 0.9, style = c("hanging", "standing", "suspended"),
  scale = c("sqrt", "raw"), expected = TRUE, confint = TRUE,
  ref = TRUE, K = NULL, xlab = NULL, ylab = NULL, main = NULL, ...) {
  if (is.null(newdata))
    newdata <- object$model.frame

  if (is.null(newdata[[object$response]]))
    stop("response missing in newdata!")

  y <- newdata[[object$response]]
  n <- length(y)

  if(!is.null(K))
    object$maxcounts <- K
  counts <- 0:object$maxcounts
  p <- NULL; obs <- NULL
  for (j in counts) {
    if (isTRUE(list(...)$verbose))
      cat(j, "/", sep = "")
    p <- cbind(p, predict(object, newdata = newdata, type = "pdf", y = j))
    obs <- c(obs, sum(y == j))
  }

  if (isTRUE(list(...)$verbose))
    cat("\n")

  e <- colMeans(p) * n  

  rg <- data.frame("observed" = obs, "expected" = e,
    "mid" = counts, "width" = width)

  scale <- match.arg(scale)

  if (scale == "sqrt") {
    rg$observed <- sqrt(rg$observed)
    rg$expected <- sqrt(rg$expected)
  }

  p <- t(p)
  rownames(p) <- paste0("p_", counts + 0.5)
  colnames(p) <- as.character(1:n)
  rg$distribution <- p

  attr(rg, "style") <- match.arg(style)
  attr(rg, "scale") <- scale
  attr(rg, "expected") <- expected
  attr(rg, "confint") <- confint
  attr(rg, "ref") <- ref
  attr(rg, "xlab") <- if (is.null(xlab)) "#Counts" else xlab
  attr(rg, "ylab") <- if (is.null(ylab)) "sqrt(Frequency)" else ylab
  attr(rg, "main") <- if (is.null(main)) "Rootogram" else main

  class(rg) <- c("rootogram", "data.frame")

  if (plot)
    plot(rg, ...)

  return(invisible(rg))
}

