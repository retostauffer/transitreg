library("TransitionModels")
library("qgam")
library("gamlss2")
library("gamlss.cens")
library("bamlss")

library("transitreg")
devtools::load_all("../")

rm(list = objects())

d <- readRDS("ehydTirol_Tageschniederschlagssummen.rds")
d <- subset(d, date >= (max(date) - 365*30))
d$sqrt_pre <- sqrt(d$value)
d$day <- as.POSIXlt(d$date)$yday

stations <- unique(d$name)

df <- subset(d, name == "Kirchberg in Tirol")

set.seed(123)
i <- sample(1:2, size = nrow(df), prob = c(0.8, 0.4), replace = TRUE)
dtrain <- subset(df, i < 2)
dtest <- subset(df, i > 1)

breaks <- c(0, seq(0.2, floor(max(df$sqrt_pre)) + 1, by = 0.3))

m <- transitreg(sqrt_pre ~ theta0 + s(theta, k = 20), data = df, breaks = breaks, censored = "left")
library("devtools"); load_all("../"); m[1]

nd <- data.frame("sqrt_pre" = seq(0, max(df$sqrt_pre), length.out = 101))
py <- nd$sqrt_pre

pm <- predict(m, newdata = nd, y = nd$sqrt_pre, type = "pdf")

mids <- (head(breaks, -1) + tail(breaks, -1)) / 2

b <- bamlss(sqrt_pre ~ 1, data = df, family = cnorm_bamlss)
par <- predict(b, newdata = data.frame("sqrt_pre" = mids), type = "parameter")
db <- family(b)$d(mids, par)

f <- sqrt_pre ~ theta0 + s(theta, k = 20) + s(day, bs = "cc", k = 20) + te(theta, day, bs = c("cr", "cc"), k = 10)

# ----------------------------------------------
breaks2 <- c(0, seq(0.2, floor(max(df$sqrt_pre)) + 1, by = 0.05))
m2 <- transitreg(f, data = dtrain, breaks = breaks2, censored = "left")
#m2 <- transitreg(f, data = dtrain, breaks = breaks2)

f <- sqrt_pre ~ s(day, k = 20, bs = "cc") | s(day, k = 20, bs = "cc")
b2 <- bamlss(f, data = dtrain, family = cnorm_bamlss, binning = TRUE)

qu <- c(0.01, 0.1, 0.5, 0.9, 0.99)

g2 <- mqgam(sqrt_pre ~ s(day, k = 20, bs = "cc"), data = dtrain, qu = qu)

nd <- data.frame("day" = 0:365)
pm2 <- predict(m2, newdata = dtest, prob = qu, elementwise = FALSE)

par <- predict(b2, newdata = dtest, type = "parameter")
pb2 <- do.call("cbind",
  lapply(qu, function(j) {
    b2$family$q(j, par)
}))

pg2 <- do.call("cbind",
  lapply(qu, function(j) {
    qdo(g2, j, predict, newdata = dtest)
}))
pg2[pg2 < 0] <- 0 # <- lower limit

err_b <- err_m <- err_g <- NULL
for(j in 1:5) {
  err_b <- c(err_b, qgam::pinLoss(dtest$sqrt_pre, pb2[, j], qu[j]))
  err_m <- c(err_m, qgam::pinLoss(dtest$sqrt_pre, pm2[, j], qu[j]))
  err_g <- c(err_g, qgam::pinLoss(dtest$sqrt_pre, pg2[, j], qu[j]))
}
err_b <- sum(err_b)
err_m <- sum(err_m)
err_g <- sum(err_g)
c(bamlss_b2 = err_b, transitreg_m2 = err_m, mqgam_g2 = err_g)




grDevices::pdf("__test.pdf", width = 12, height = 7)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

hist(df$sqrt_pre, breaks = breaks, freq = FALSE,
  xlab = "sqrt(Precipitation)", main = "transitreg")

lines(pm ~ py, col = 4, lwd = 2)
lines(db ~ mids, col = 2, lwd = 2)
rug(df$sqrt_pre, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))

legend("center", c("TM", "CN"),
  lwd = 2, col = c(4, 2), bty = "n")

plot(sqrt_pre ~ day, data = dtest, type = "h", col = rgb(0.1, 0.1, 0.1, alpha = 0.4),
  xlab = "Day of the year", ylab = "sqrt(Precipitation)", ylim = c(0, 11))

j <- order(dtest$day)
matplot(dtest$day[j], pb2[j, ], type = "l", lty = 1, col = 2, add = TRUE)
matplot(dtest$day[j], pg2[j, ], type = "l", lty = 1, col = 3, add = TRUE)
matplot(dtest$day[j], pm2[j, ], type = "l", lty = 1, col = 4, add = TRUE)

err <- c("CN" = err_b, "QR" = err_g, "TM" = err_m)
col <- c(2, 3, 4)
i <- order(err)
err <- err[i]
col <- col[i]

legend("topleft", paste(paste(names(err), "PBL ="), round(err)),
  lwd = 2, col = col, bty = "n")

dev.off()
