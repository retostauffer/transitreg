library("TransitionModels")
library("gamlss2")

library("ggplot2")
library("patchwork")
library("colorspace")
library("dplyr")
library("sf")
library("rnaturalearth")

# Loading WeatherGermay data set
data("WeatherGermany", package = "WeatherGermany")
data <- subset(WeatherGermany, !is.na(Wmax))
head(data)
range(data$date)

##PERIOD <- list(start = as.Date("2015-01-01"), end = as.Date("2022-12-01"))
PERIOD <- list(start = as.Date("2000-01-01"), end = as.Date("2022-12-01"))
# --- period 2000-01-01 to 2022-12-01 with 100 breaks took 1484 seconds (42 minutes)
x <- subset(data, date >= PERIOD$start & date <= PERIOD$end)
x$year <- as.integer(format(x$date, "%Y"))
counts <- x %>% group_by(id) %>% summarize(count = length(Wmax))

fn <- function(x) head(subset(x, select = c(alt, lon, lat)), 1)
counts <- counts %>% right_join(x %>% group_by(id) %>% group_modify(~fn(.x)))
counts <- na.omit(counts)
message(" Left with ", paste(dim(counts), collapse = " x "))

counts <- st_as_sf(counts, coords = c("lon", "lat"), crs = st_crs(4326))
germany <- st_geometry(subset(ne_countries("medium", returnclass = "sf"), iso_a3 == "DEU"))

# ---------------------------------------------------------------
# Geodata
# https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/digitale-gelandemodelle/digitales-gelandemodell-gitterweite-1000-m-dgm1000.html
# Using the 1000m gridded ASCII
# ---------------------------------------------------------------
library("stars")
dem <- read_stars("dgm1000_utm32s.asc")
st_crs(dem) <- st_crs(25832) # UTM32n
dem <- st_transform(dem, st_crs(4326)) # to WGS84
#install.packages("FNN")
dem <- st_warp(dem, st_as_stars(st_bbox(dem), dx = 0.05))
plot(dem)


# ---------------------------------------------------------------
# Quick plot
# ---------------------------------------------------------------
ggplot() +
    geom_stars(data = dem) +
    scale_fill_gradientn(colours = hcl.colors(20, "terrain")) +
    geom_sf(data = germany, fill = NA) +
    geom_sf(aes(color = count), data = counts, cex = 3) +
    labs(title = sprintf("Wmax obs counts (%s-%s)", min(x$date), max(x$date))) +
    theme_minimal()


# Estimating the model
f <- Wmax ~ ti(alt) + ti(lon, lat) + ti(year) + ti(alt, year)
system.time(b <- tm(f, data = x, breaks = 100, useC = TRUE))
saveRDS(b, "____b.rds")

dim(b$model.frame)
dim(model.frame(b$model))

nd <- as.data.frame(dem)
names(nd) <- c("lon", "lat", "alt")
nd$year <- 2022


# ---------------------------------------------------------------

p5  <- predict(b, newdata = nd, useC = TRUE, type = "pdf", y = 5)
p12 <- predict(b, newdata = nd, useC = TRUE, type = "pdf", y = 12)

dem12 <- dem5 <- dem
dem5[[1]] <- matrix(p5, ncol = ncol(dem5[[1]]))
names(dem5) <- "P(y>5)"
dem12[[1]] <- matrix(p12, ncol = ncol(dem12[[1]]))
names(dem12) <- "P(y>12)"
g5  <- ggplot() + geom_stars(data = dem5) +
    scale_fill_continuous_sequential("ag_sunset") +
    labs(title = "P(Wmax > 5)") + theme_minimal()
g12 <- ggplot() + geom_stars(data = dem12) +
    scale_fill_continuous_sequential("ag_sunset") +
    labs(title = "P(Wmax > 12)") + theme_minimal()
g5 + g12
ggsave(g5 + g12, filename = "__wind_p.png", width = 12, height = 6)

# ---------------------------------------------------------------

demq50 <- dem
names(demq50) <- "median"
q50 <- predict(b, newdata = nd, type = "quantile", prob = 0.5, useC = TRUE)
demq50[[1]] <- matrix(q50, ncol = ncol(demq50[[1]]))
gq50 <- ggplot() + geom_stars(data = demq50) +
           scale_fill_continuous_sequential("Hawaii") +
           labs(title = "Median") + theme_minimal()

demq75 <- dem
names(demq75) <- "q75"
q75 <- predict(b, newdata = nd, type = "quantile", prob = 0.9, useC = TRUE)
demq75[[1]] <- matrix(q75, ncol = ncol(demq75[[1]]))
gq75 <- ggplot() + geom_stars(data = demq75) +
           scale_fill_continuous_sequential("Hawaii") +
           labs(title = "75er Percentile") + theme_minimal()

gq50 + gq75
ggsave(gq50 + gq75, filename = "__wind_q.png", width = 12, height = 6)



