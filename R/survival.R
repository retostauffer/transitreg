


#' Creating Survival Object
#'
#' Creates a survival object used as response variable in a transitreg
#' model formula.
#'
#' @param time the time to event, the varaible the binning is taking place on.
#' @param event must evaluate to logical, whether or not an event happened
#'        (`FALSE` equals survived, `TRUE` equals dead). If `event` is a single
#'        logical it is recycled for all observations, else the length of
#'        the two vectors `time` and `event` must match exactely.
#'
#' @return An object of class `survival` which can be used as response
#' in transitreg models.
#'
#' @author Reto
#' @export
survival <- function(time, event) {
    event <- as.logical(event)
    if (length(event) == 1L) event <- rep(event, length(time))
    stopifnot(
        "argument `event` must evaluate to logical" = is.logical(event),
        "length of `time` and `event` must match" = length(time) == length(event)
    )

    res <- cbind(time = time, event = event) |> structure(class = "survival")
    return(res)
}

#' @exportS3Method as.character survival
#' @rdname survival
as.character.survival <- function(x, ...) {
    sprintf("%s%s", format(x[, "time"]), c("+", "")[as.integer(x[, "event"]) + 1L])
}

#' @exportS3Method print survival
#' @rdname survival
print.survival <- function(x, quote = FALSE, ...) {
    invisible(print(as.character(x), quote = quote))
}

#' @exportS3Method format survival
#' @rdname survival
format.survival <- function(x, ...) {
    format(as.character(x))
}

#' @exportS3Method `[` survival
#' @rdname survival
`[.survival` <- function(x, i, ...) {
    cls <- class(x)
    if (!missing(i)) {
        x <- unclass(x)[i, , drop = FALSE]
        class(x) <- cls
    } else {
        x <- NextMethod("[")
    }
    return(x)
}






