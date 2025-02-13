
library("transitreg")
library("bamlss")
library("qgam")

PNG <- "stauffer-figure1.png"
if (!file.exists(PNG)) {
  df <- readRDS("observations_stn4811.rds")

  set.seed(6020)
  idx    <- sample(1:2, size = nrow(df), prob = c(0.8, 0.4), replace = TRUE)
  dtrain <- subset(df, idx == 1L)
  dtest  <- subset(df, idx == 2L)

  breaks <- c(0, seq(0.3, 12, by = 0.2))
  f1 <- sqrt_rain ~ theta0 + ti(theta, k = 20)
  m1 <- transitreg(f1, data = df, breaks = breaks)


  xx <- c(0, 0.3, ((head(breaks, -1) + tail(breaks, -1L)) / 2)[-1])
  nd <- data.frame("sqrt_rain" = xx)
  mids <- nd$sqrt_rain

  py <- seq(0, 12, by = 0.01) #nd$sqrt_rain
  pm <- as.vector(predict(m1, newdata = nd[1,,drop = F], y = py, type = "pdf"))

  b1 <- bamlss(sqrt_rain ~ 1, data = df, family = cnorm_bamlss)
  par <- predict(b1, newdata = data.frame("sqrt_rain" = mids), type = "parameter")
  db <- family(b1)$d(mids, par)


  # ----------------------------------------------
  f2 <- sqrt_rain ~ theta0 + ti(theta, k = 20) + ti(theta, day, bs = c("cr", "cc"), k = c(20, 20))
  m2 <- transitreg(f2, data = dtrain, breaks = breaks)

  fb2 <- sqrt_rain ~ s(day, k = 20, bs = "cc") | s(day, k = 20, bs = "cc")
  b2  <- bamlss(fb2, data = dtrain, family = cnorm_bamlss, binning = TRUE)

  qu <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  g2 <- mqgam(sqrt_rain ~ s(day, k = 20, bs = "cc"), data = dtrain, qu = qu)

  nd <- data.frame("day" = 1:365)

  devtools::load_all("~/Software/transitreg")
  pm2 <- predict(m2, newdata = dtest, prob = qu, elementwise = FALSE)
  pm2[, 1:2] <- 0 # Reto: quick fix, must be properly investigated in the software!

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
    err_b <- c(err_b, qgam::pinLoss(dtest$sqrt_rain, pb2[, j], qu[j]))
    err_m <- c(err_m, qgam::pinLoss(dtest$sqrt_rain, pm2[, j], qu[j]))
    err_g <- c(err_g, qgam::pinLoss(dtest$sqrt_rain, pg2[, j], qu[j]))
  }
  print(rbind(bamlss_b2     = err_b,
              transitreg_m2 = err_m,
              mqgam_g2      = err_g))
  err_b <- sum(err_b)
  err_m <- sum(err_m)
  err_g <- sum(err_g)
  print(sort(c(bamlss_b2 = err_b, transitreg_m2 = err_m, mqgam_g2 = err_g)))

  ## Create figure
  png(PNG, units = "in", res = 200, width = 8, height = 4)
    par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

    # First subplot
    hist(df$sqrt_rain, breaks = breaks, freq = FALSE,
      xlab = "sqrt(Precipitation)", main = NA)

    lines(pm ~ py, col = 4, lwd = 2)
    lines(db ~ mids, col = 2, lwd = 2)
    rug(df$sqrt_rain, col = rgb(0.1, 0.1, 0.1, alpha = 0.4))

    legend("center", c("TM", "CN"),
      lwd = 2, col = c(4, 2), bty = "n")

    # Second subplot
    plot(sqrt_rain ~ day, data = dtest, type = "h", col = rgb(0.1, 0.1, 0.1, alpha = 0.4),
      xlab = "Day of the year", ylab = "sqrt(Precipitation)", ylim = c(0, 8))

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

