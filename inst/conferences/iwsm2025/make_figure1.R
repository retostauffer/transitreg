
library("transitreg")
library("bamlss")
library("qgam")

PDF <- "stauffer-figure1.pdf"
if (!file.exists(PDF)) {
  data(Shannon)
  head(Shannon)
  Shannon <- transform(Shannon,
                       day = as.integer(format(date, "%j")),
                       year = as.integer(format(date, "%Y")))


  #set.seed(1)
  set.seed(666)
  idx    <- sample(1:2, size = nrow(Shannon), prob = c(2/3, 1/3), replace = TRUE)
  dtrain <- subset(Shannon, idx == 1L)
  dtest  <- subset(Shannon, idx == 2L)

  breaks <- seq(0, 12, by = 0.3)
  f1 <- sqrt(rain) ~ theta0 + ti(theta, k = 20)
  m1 <- transitreg(f1, data = Shannon, breaks = breaks, censored = "left")

  xx <- c(0, ((head(breaks, -1) + tail(breaks, -1L)) / 2)[-1])
  nd <- data.frame("rain" = xx^2)
  mids <- nd$rain

  py <- seq(0, 12, by = 0.01)
  pm <- as.vector(predict(m1, newdata = nd[1,,drop = F], y = py, type = "pdf"))

  b1file <- "__bamlss_b1.rds"
  if (!file.exists(b1file)) {
    b1 <- bamlss(sqrt(rain) ~ 1, data = Shannon, family = cnorm_bamlss)
    saveRDS(b1, b1file)
  } else {
    b1 <- readRDS(b1file)
  }

  xb <- seq(0, 12, by = 0.01)
  par <- predict(b1, newdata = data.frame("sqrt_rain" = xb), type = "parameter")
  db <- family(b1)$d(xb, par)


  # ----------------------------------------------
  f2 <- sqrt(rain) ~ theta0 + ti(theta, k = 20) + ti(theta, day, bs = c("cr", "cc"), k = c(20, 20))
  m2 <- transitreg(f2, data = dtrain, breaks = breaks, censored = "left")

  fb2 <- sqrt(rain) ~ s(day, k = 20, bs = "cc") | s(day, k = 20, bs = "cc")

  b2file <- "__bamlss_b2.rds"
  if (!file.exists(b2file)) {
    b2  <- bamlss(fb2, data = dtrain, family = cnorm_bamlss, binning = TRUE)
    saveRDS(b2, b2file)
  } else {
    b2 <- readRDS(b2file)
  }

  qu <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  g2file <- "__mqgam_g2.rds"
  if (!file.exists(g2file)) {
    g2 <- mqgam(sqrt(rain) ~ s(day, k = 20, bs = "cc"), data = dtrain, qu = qu)
    saveRDS(g2, g2file)
  } else {
    g2 <- readRDS(g2file)
  }

  nd <- data.frame("day" = 1:365)

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

  err_b <- err_m <- err_g <- NULL
  for(j in 1:5) {
    err_b <- c(err_b, qgam::pinLoss(sqrt(dtest$rain), pb2[, j], qu[j]))
    err_m <- c(err_m, qgam::pinLoss(sqrt(dtest$rain), pm2[, j], qu[j]))
    err_g <- c(err_g, qgam::pinLoss(sqrt(dtest$rain), pg2[, j], qu[j]))
  }
  print(rbind(bamlss_b2     = err_b,
              transitreg_m2 = err_m,
              mqgam_g2      = err_g))
  err_b <- sum(err_b)
  err_m <- sum(err_m)
  err_g <- sum(err_g)
  pinball_loss <- sort(c(bamlss_b2 = err_b, transitreg_m2 = err_m, mqgam_g2 = err_g))
  print(pinball_loss)
  print(pinball_loss[which.min(pinball_loss)])


  ## Create figure
  ##grDevices::pdf(file = "stauffer-figure1.pdf", width = 8, height = 2.5)
  grDevices::pdf(file = PDF, width = 8 * 0.9, height = 2.5 * 0.8)
    par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))

    # First subplot
    col    <- rep(c("white", "gray80"), times = c(1, length(breaks) - 2))
    border <- rep(c("white", "black"),  times = c(1, length(breaks) - 2))
    h <- hist(sqrt(Shannon$rain),
              breaks = breaks,
              col = col,
              border = border,
              freq = FALSE,
              xlab = "y; sqrt(Precipitation)", main = NA,
              xlim = c(0, 8), ylim = c(0, 0.4))

    # Adding observed frequency of zero
    obs_zero <- mean(Shannon$rain == 0)
    points(0, obs_zero, col = 1, pch = 19)

    lines(pm ~ py, col = 4, lwd = 2)
    points(py[1L], pm[1L], col = 4, pch = 1)
    lines(db ~ xb, col = 2, lwd = 2)
    points(xb[1L], db[1L], col = 2, pch = 1)
    rug(sqrt(Shannon$rain), col = rgb(0.1, 0.1, 0.1, alpha = 0.4))

    legend("right", c("TM", "CN"),
      lwd = 2, col = c(4, 2), bty = "n")
    legend("topright", legend = c(
                expression(paste("observed P(", y, "=", 0, ")")),
                expression(paste("modelled P(", y, "=", 0, ")"))
           ), pch = c(19, 1), col = "gray50", bty = "n")

    # Second subplot
    plot(sqrt(rain) ~ day, data = dtest, type = "h", col = rgb(0.1, 0.1, 0.1, alpha = 0.4),
      xlab = "Day of the year", ylab = "y; sqrt(Precipitation)", ylim = c(0, 8.5),
      xlim = c(1, 365), xaxs = "i", yaxs = "i")

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

    ####
    target_day <- as.integer(format(as.Date("2025-07-16"), "%j"))
    distr <- prodist(m2)[which(dtrain$day == target_day)[1]]
    plot(distr, cdf = TRUE, tp = TRUE, main = NULL, col = 4, ylab = NA,
         xlab = "y; sqrt(Precipitation)")
    legend("bottomleft", legend = c(
                expression(paste("CDF: P(", tilde(y) <= r, ")")),
                expression(paste("TP: P(", tilde(y) > r, "|", tilde(y) >= r, ")"))
           ),
           bty = "n",
           pch = c(19, NA),
           lwd = c(2, 1), col = c(4, 4), lty = c(1, 2))

  dev.off()
}

#if(!file.exists("premodel.png") & FALSE) {
#  d <- readRDS("../ehydTirol_Tageschniederschlagssummen.rds")
#  d <- subset(d, date >= (max(date) - 365*30))
#  d$sqrt_rain <- sqrt(d$value)
#  d$day <- as.POSIXlt(d$date)$yday
#
#  stations <- unique(d$name)
#
#  set.seed(123)
#  i <- sample(stations, size = floor(length(stations) * 0.8))
#
#  dtrain <- subset(d, name %in% i)
#  dtest <- subset(d, !(name %in% i))
#
#  breaks <- c(-0.05, seq(0.05, floor(max(d$sqrt_rain)) + 1, by = 0.1))
#
#  f1 <- sqrt_rain ~ theta0 +
#    ti(theta) +
#    ti(alt) +
#    ti(day,bs="cc",k=20) +
#    ti(lon,lat,bs="tp",d=2,k=30) +
#    ti(day,lon,lat,d=c(1,2),k=c(5,5,30),bs=c("cc","tp")) +
#    ti(theta,day,lon,lat,d=c(1,1,2),k=c(5,5,30),bs=c("cr","cc","tp")) +
#    ti(theta,day,bs=c("cr","cc")) +
#    ti(theta,alt) +
#    ti(theta,lon,lat,d=c(1,2),bs=c("cr","tp"),k=c(5,30))
#
#  m1 <- transitreg(f1, data = dtrain, breaks = breaks)
#
#  f2 <- sqrt_rain ~ s(day,bs="cc",k=20) + s(alt,k=5) + s(lon,lat) +
#    te(day,lon,lat,d=c(1,2),bs=c("cr","tp"),k=c(5,30)) |
#      s(day,bs="cc",k=8) + s(alt,k=5) + s(lon,lat) +
#      te(day,lon,lat,d=c(1,2),bs=c("cr","tp"),k=c(5,30))
#
#  m2 <- bamlss(f2, data = dtrain, family = cnorm_bamlss,
#    binning = TRUE, light = TRUE)
#
#  qu <- c(0.01, 0.1, 0.5, 0.9, 0.99)
#  err1 <- err2 <- NULL
#  par <- predict(m2, newdata = dtest, type = "parameter")
#  for(j in qu) {
#    pj1 <- predict(m1, newdata = dtest, prob = j)
#    pj2 <- family(m2)$q(j, par)
#    err1 <- c(err1, sum(qgam::pinLoss(dtest$sqrt_rain, pj1, j)))
#    err2 <- c(err2, sum(qgam::pinLoss(dtest$sqrt_rain, pj2, j)))
#  }
#
##R> err1
##[1]   2314.927  23149.270 112589.150  71447.505  11941.508
##R> err2
##[1]   2314.927  23149.270 112545.307  71472.810  12087.671
##R> round((err1 - err2) / err2 * 100, 2)
##[1]  0.00  0.00  0.04 -0.04 -1.21
#
#  download.file("http://bamlss.org/misc/precipitation_clim.tar.gz", "clim.tar.gz")
#  untar("clim.tar.gz", exdir = ".")
#
#  download.file("https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_AUT_1.json.zip", "AT-gadm.zip")
#  unzip("AT-gadm.zip")
#
#  library("raster")
#  library("stars")
#  library("sf")
#
#  dem <- st_as_stars(readRDS("dem.rds"))
#  dem <- st_transform(dem, crs = 4326)
#
#  Tyrol <- "gadm41_AUT_1.json" |>
#    read_sf() |>
#    subset(NAME_1 == "Tirol") |>
#    st_geometry()
#
#  dem <- mask(readRDS("dem.rds"), as(Tyrol, "Spatial"))
#  dem <- crop(dem, as(Tyrol, "Spatial"))
#
#  writeRaster(dem, "dem.tif", format = "GTiff", overwrite = TRUE)
#
#  nd <- as.data.frame(coordinates(dem))
#  names(nd) <- c("lon", "lat")
#  nd$alt <- values(dem)
#
#  days <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
#  names(days) <- month.name
#
#  for(j in names(days)) {
#    print(j)
#
#    nd$day <- days[j]
#    nd[[paste0("m1_", j)]] <- predict(m1, newdata = nd, prob = 0.99)
#
#    par <- predict(m2, newdata = nd, type = "parameter")
#    par <- as.data.frame(par)
#    p2 <- family(m2)$q(0.99, par)
#    nd[[paste0("m2_", j)]][!is.na(nd$alt)] <- p2
#  }
#
#  r_list <- lapply(names(days), function(j) {
#    a <- nd[, c("lon", "lat", paste0("m1_", j))]
#    a[[paste0("m1_", j)]] <- a[[paste0("m1_", j)]]^2
#    rasterFromXYZ(a)
#  })
#  pred <- stack(r_list)
#  names(pred) <- gsub("m1_", "", names(pred))
#  crs(pred) <- crs(dem)
#
#  library("ggplot2")
#  library("dplyr")
#  library("tidyr")
#
#  pred_df <- as.data.frame(pred, xy = TRUE)
#  pred_df <- pred_df[, c("x", "y", names(pred))]
#  
#  pred_long <- pivot_longer(
#    pred_df, 
#    cols = names(pred), 
#    names_to = "Season", 
#    values_to = "value"
#  )
#
#  pred_long$Season <- factor(pred_long$Season, levels = names(pred))
#
#  rain_colors <- function (n, h = c(-160, -38), c. = c(10, 80), l = c(86, 39), 
#    power = c(2.75806451612903, 1), fixup = FALSE, gamma = NULL, 
#    alpha = 1, ...) 
#  {
#    if (!is.null(gamma)) 
#        warning("'gamma' is deprecated and has no effect")
#    if (n < 1L) 
#        return(character(0L))
#    h <- rep(h, length.out = 2L)
#    c <- rep(c., length.out = 2L)
#    l <- rep(l, length.out = 2L)
#    power <- rep(power, length.out = 2L)
#    rval <- seq(1, 0, length = n)
#    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
#        C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) * 
#            rval), fixup = fixup, ...)
#    if (!missing(alpha)) {
#        alpha <- pmax(pmin(alpha, 1), 0)
#        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
#            width = 2L, upper.case = TRUE)
#        rval <- paste(rval, alpha, sep = "")
#    }
#    return(rval)
#  }
#
#  rc <- rain_colors(100)
#
#  png("stauffer-figure1.png", units = "in", res = 200, width = 13, height = 10)
#
#  xr <- range(pred_long$x)
#  xb <- round(seq(xr[1] + 0.1*abs(diff(xr)), xr[2] - 0.1*abs(diff(xr)), length = 5), 2)
#  yr <- range(pred_long$y)
#  yb <- round(seq(yr[1] + 0.1*abs(diff(yr)), yr[2] - 0.1*abs(diff(yr)), length = 5), 2)
#
#  ggplot() +
#    geom_tile(data = pred_long, aes(x = x, y = y, fill = value)) +
#    geom_sf(data = Tyrol, fill = NA, color = "black", size = 0.5) + # Overlay Tyrol border
#    facet_wrap(~ Season, nrow = 4, ncol = 3) + # Correct ordering
#    scale_fill_gradientn(colors = rc, na.value = "transparent") + # Custom color scale
#    scale_x_continuous(breaks = xb) + # Adjust x-axis ticks
#    scale_y_continuous(breaks = yb) + # Adjust y-axis ticks
#    coord_sf(crs = st_crs(Tyrol)) + # Respect CRS and aspect ratio
#    theme_minimal() +
#    labs(title = NULL,
#         x = "Longitude",
#         y = "Latitude",
#         fill = "Precipitation\n[mm]")
#
#   dev.off()
#}

