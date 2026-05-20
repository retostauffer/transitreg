# -------------------------------------------------------------------
# Simulation study
# Count data (marginal)
# -------------------------------------------------------------------

library("transitreg")
library("gamlss2")
# For evaluation
library("distributions3")
library("scoringRules")
library("topmodels")
library("ggplot2")

rm(list = objects())
source("functions.R")


# -------------------------------------------------------------------
# -------------------------------------------------------------------
simulate_data_negbin <- function(size, prob, n = 1000) {
    data.frame(y = rnbinom(n, size, prob))
}

data   <- simulate_data_negbin(5, 0.1, 1000)
gam    <- gamlss2(y ~ 1, data = data, family = NBI, trace = TRUE)
tra    <- transitreg(y ~ s(theta), data = data)
models <- list(gamlss_NBI = gam, transitreg = tra)

# Calculating CRPS
crps <- get_crps(models, data$y)

par(mfrow = c(1, 2))
(plot_crps(crps))
(plot_crpss(crps))

cbind(logLik = sapply(models, logLik),
      AIC    = sapply(models, AIC),
      BIC    = sapply(models, BIC))

## compute and plot wormplot as base graphic
w <- lapply(models, function(m) wormplot(m, plot = FALSE, simint = FALSE))
autoplot(do.call(c, w), single_graph = TRUE, col = seq_along(w), legend = TRUE)


## Comparing 'quantiles'?
sim <- function(models, n, size, prob, ncores = NULL) {
    stopifnot(length(size) == length(prob))
    result <- list()
    for (i in seq_along(size)) {
        d <- simulate_data_negbin(size[[i]], prob[[i]], n)

        tgam <- system.time(gam <- gamlss2(y ~ 1, data = d, family = NBI, trace = TRUE))
        ttra <- system.time(tra <- transitreg(y ~ s(theta), data = d, ncores = ncores))

        models <- list(gamlss_NBI = gam, transitreg = tra)
        # Time elapsed: Take care of order!
        time_elapsed <- c(gam = tgam["elapsed"], transitreg = ttra["elapsed"])

        crps <- get_crps(models, d$y)
        mean_crps  <- plot_crps(crps, plot = FALSE)

        tmp <- data.frame(n = n, size = size[i], prob = prob[i],
                          model  = names(models),
                          logLik = sapply(models, logLik),
                          AIC    = sapply(models, AIC),
                          BIC    = sapply(models, BIC),
                          time   = time_elapsed,
                          mean_crps = mean_crps[,2])
        result[[i]] <- tmp
        rm(tmp, d, tgam, gam, ttra, tra, models, crps, mean_crps, mean_crpss)
    }

    result <- do.call(rbind, result)
    return(result |> structure(row.names = seq_len(nrow(result))))
}

m      <- 20 # Number of repetitions
size   <- runif(m, 2, 5)
prob   <- runif(m, 0.05, 0.5)
res <- sim(models, 10000, size, prob, ncores = 10)
#head(res, 2)

par(mfrow = c(2, 2))
boxplot(logLik    ~ model, data = res)
boxplot(mean_crps ~ model, data = res)
boxplot(time      ~ model, data = res)
boxplot(BIC       ~ model, data = res)

# Some testing (m = 20, n = 10.000)
# on one single CPU (ncores = 1): avg. time for transitreg   1.4  seconds
# same on 10 cores (ncores = 10): avg. time for transitreg   0.40 seconds
# gamlss2 on average takes                                   0.46 seconds
aggregate(time ~ model, data = res, FUN = mean)
