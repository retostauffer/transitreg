library("TransitionModels")
library("gamlss2")

if (!file.exists("__tmp.rda")) {
    library("dplyr")
    library("sf")
    library("rnaturalearth")
    
    rm(list = objects())
    
    # Loading WeatherGermay data set
    data("WeatherGermany", package = "WeatherGermany")
    data <- subset(WeatherGermany, !is.na(Wmax))
    head(data)
    
    STARTYEAR <- 2012
    #STARTYEAR <- 2019
    ENDYEAR   <- 2019
    x <- subset(data, id == 722 &
                      date >= as.Date(sprintf("%04d-01-01", STARTYEAR)) &
                      date <= as.Date(sprintf("%04d-12-31", ENDYEAR)))
    
    x <- subset(x, select = c(id, date, Wmax, name, alt, lat, lon))
    x$yday <- as.integer(format(x$date, "%j"))
    x$year <- as.integer(format(x$date, "%Y"))
    message("Dimension of training data: ", paste(dim(x), collapse = " x "))
    
    # Estimating the model
    f <- Wmax ~ s(yday, bs = "cc")
    #f <- Wmax ~ s(yday, bs = "cc")
    system.time(b <- tm(f, data = x, breaks = 100))
    
    p <- cbind(p10 = predict(b, newdata = x, type = "quantile", prob = .1),
               p50 = predict(b, newdata = x, type = "quantile", prob = .5),
               p90 = predict(b, newdata = x, type = "quantile", prob = .9))
    matplot(p, type = "l", lwd = 2, lty = 1,
            main = "Stn 722, year 2019, Percentiles 10/50/90")

    save(b, x, file = "__tmp.rda")
} else {
    load("__tmp.rda")
    print(objects())
}

#devtools::load_all("../"); dead_p <- predict(b, newdata = x, type = "quantile", prob = .5)
devtools::load_all("../")
#dead_p <- predict(b, newdata = x, type = "pdf", useC = FALSE)
#dead_p <- predict(b, newdata = x, type = "pdf", useC = TRUE)

dead_p <- predict(b, newdata = x, type = "cdf", useC = FALSE)
dead_p <- predict(b, newdata = x, type = "cdf", useC = TRUE)

#system.time(dead_p <- predict(b, newdata = x, type = "pdf"))
#system.time(dead_p <- predict(b, newdata = x, type = "quantile"))
#system.time(dead_p <- predict(b, newdata = x, type = "cdf"))

