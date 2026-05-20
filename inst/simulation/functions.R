
## Helper function to calculate CRPS of multiple models,
## will return the results (elementwise) in long format.
get_crps <- function(models, y) {
    res <- lapply(models, function(m) crps(prodist(m), y))
    data.frame(model = rep(names(models), each = length(y)), crps = do.call(c, res)) |>
        structure(row.names = seq_len(length(y) * 3))
}

## Calculating skill scores against 'transitreg'
calc_crps_skill_scores <- function(x, offset = -1) {
    stopifnot(
        is.data.frame(x),
        all(c("model", "crps") %in% names(x)),
        "transitreg" %in% x$model
    )
    # Long to wide first
    res <- list()
    for (m in unique(x$model)) res[[m]] <- x$crps[x$model == m]
    res <- data.frame(res)
    for (n in names(res)) {
        if (n == "transitreg") next
        res[[n]] <- res[[n]] / res$transitreg
    }
    res <- subset(res, select = -transitreg)
    reto <<- res
    # Convert wide to long again
    data.frame(model = rep(names(res), each = nrow(res)),
               crpss = unname(do.call(c, res)) + offset)
}
plot_crps <- function(x) {
    # Elementwise CRPS + mean(CRPS) labelled at the top
    boxplot(crps ~ model, data = x)
    tmp <- aggregate(crps ~ model, data = x, mean)
    axis(side = 3, at = seq_len(nrow(tmp)), format(tmp$crps))
    invisible(tmp)
}
plot_crpss <- function(x, offset = -1) {
    x <- calc_crps_skill_scores(x, offset)
    # Elementwise CRPSS + mean(CRPSS) labelled at the top
    boxplot(crpss ~ model, data = x); abline(h = 1 + offset, col = "tomato")
    tmp <- aggregate(crpss ~ model, data = x, mean)
    points(tmp$crpss, pch = 19, col = 2, cex = 2)
    axis(side = 3, at = seq_len(nrow(tmp)), format(tmp$crps))
    invisible(tmp)
}
