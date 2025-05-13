

library('transitreg')
data(weatherShannon)
d <- transform(weatherShannon, day = as.integer(format(date, "%j")))
get_sunshine_by_doy <- function(doy, latitude) {
    # Solar declination in radians
    decl <- 23.44 * pi / 180 * sin(2 * pi * (284 + doy) / 365)
    # Clamp for acos domain (latitudes near poles during solstices can exceed range)
    x <- pmax(-1, pmin(1, -tan(latitude * pi / 180) * tan(decl)))
    #x[x < -1] <- -1
    #x[x >  1] <-  1
    # Convert hour angle to daylight hours
    2 * acos(x) * 180 / pi / 15
}

# Example usage:
d$maxsun <- get_sunshine_by_doy(d$day, latitude = 52.69028)
head(d)
plot(maxsun ~ day, data = d, ylim = c(0, 17))
with(d, points(day, sun, col = 2))


