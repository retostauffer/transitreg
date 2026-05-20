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
simulate_data_poisson <- function(lambda, n = 1000) {
    data.frame(y = rpois(n, lambda))
}

data   <- simulate_data_poisson(5.5, 1000)
glm    <- glm(y ~ 1, data = data, family = "poisson")
gam1   <- gamlss2(y ~ 1, data = data, family = PO, trace = TRUE)
gam2   <- gamlss2(y ~ 1, data = data, family = NBI, trace = TRUE)
tra    <- transitreg(y ~ s(theta), data = data)
models <- list(glm = glm, gamlss_PO = gam1, gamlss2_NBI = gam2, transitreg = tra)

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





