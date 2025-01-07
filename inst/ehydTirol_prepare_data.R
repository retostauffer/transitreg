#!/usr/bin/env Rscript
# -------------------------------------------------------------------
# Downloading and prearing daily precipitation sums for Tyrol.
#
# https://ehyd.gov.at now has a 'download' button where you can
# download all daily observations at once. The file can, today,
# also be downloaded via the URL below (may change over time?).
# In this case make your download manually and call it 'ehydTirol.zip',
# the rest of the script should work if they do now change the format.
# -------------------------------------------------------------------

library("httr")

remove(list = objects()) # Just to be safe

# Some definitions
URL <- "https://ehyd.gv.at/eHYD/AreaSelection/download?cat=nlv&reg=7"
LOCALZIP <- "ehydTirol.zip"
RDSFILE <- "ehydTirol_Tageschniederschlagssummen.rds"

if (file.exists(RDSFILE)) {
    stop("File ", RDSFILE, " already exists, stopping here.")
}

# Trying to download the data if required
if (!file.exists(LOCALZIP)) {
    req <- POST(URL, write_disk(LOCALZIP))
    if (!status_code(req) == 200) stop("Whoops, problems downloading the ZIP file")
}

# -------------------------------------------------------------------
# Else we unzip the file and read the data
# -------------------------------------------------------------------
tmpdir <- tempfile()
files  <- unzip(LOCALZIP, exdir = tmpdir)
files  <- files[grepl("N-Tagessummen-[0-9]+\\.csv$", files)]
cat("Files to be read:", length(files), "\n")

get_data <- function(f, enc = "latin1", n = 100) {
    stopifnot(isTRUE(file.exists(f)))
    stopifnot(grepl("\\.csv$", f))

    # Reading the fist 100 lines of the CSV to check where 'Werte:' starts
    tmp <- readLines(f, n = n, encoding = enc)
    n   <- grep("Werte:$", tmp)

    # Reading header and data
    if (length(n) == 0)
        stop("Could not find 'Werte:' line or 'n' is not large enough!")
    header <- readLines(f, n = n, encoding = enc)

    # Extracting header information. We need
    # - Messstelle: Name
    #
    # Demo header:
    # HZB-Nummer:                ;111831
    # ...
    # Höhe:
    #  gültig seit:              ;Höhe [m ü.A.]:
    #   01.01.1902               ;710
    #   01.04.2012               ;753
    # Geographische Koordinaten (Referenzellipsoid: Bessel 1841):
    #  gültig seit:              ;Länge (Grad,Min,Sek):    ;Breite  (Grad,Min,Sek):
    #   01.04.2012               ;14 46 52                 ;47 24 47
    #   01.01.1902               ;14 49 25                 ;47 24 01
    name <- regmatches(header, regexpr("(?<=(Messstelle:)).*+", header, perl = TRUE))
    stopifnot(length(name) == 1L)
    name <- trimws(gsub("\\s+;", "", name))

    # This function extracts the 'altitude' definition which may 
    # have changed over the years, when the station was re-located.
    get_alt <- function(x) { 
        x <- regmatches(x, regexpr("\\s+[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}\\s+;\\s{0,}[0-9\\.]+$", x))
        stopifnot(length(x) > 0L)
        x <- setNames(read.table(text = x, sep = ";", strip.white = TRUE), c("date", "alt"))
        x$date <- as.Date(x$date, format = "%d.%m.%Y")
        return(x)
    }
    alt <- get_alt(header)

    # As 'get_alt' this function is extracting the geographical
    # location; may have been changed over time as well.
    get_loc <- function(x) {
        x <- regmatches(x, regexpr("\\s+[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}\\s+;[0-9\ ]+;[0-9\ ]+", x))
        stopifnot(length(x) > 0L)
        x <- setNames(read.table(text = x, sep = ";", strip.white = TRUE), c("date", "lon", "lat"))
        x$date <- as.Date(x$date, format = "%d.%m.%Y")
        # degrees, mins, secs to decimal
        str2num <- function(x) {
            fn <- function(y) {
                y <- as.numeric(y)
                y[1] + y[2] / 60 + y[3] / 3600
            }
            sapply(regmatches(x, gregexpr("[0-9]{2}", x)), fn)
        }
        for (n in c("lon", "lat")) x[[n]] <- str2num(x[[n]])
        return(x)
    }
    loc <- get_loc(header)

    # Reading data (setting 'kill' to NA)
    kill <- "Lücke"
    res   <- setNames(read.table(f, skip = n, header = FALSE, sep = ";",
                                 strip.white = TRUE,
                                 dec = ",", encoding = enc, na.strings = kill),
                       c("date", "value"))
    res$date <- as.Date(as.POSIXct(res$date, format = "%d.%m.%Y %H:%M:%S", tz = "CET"))
    res$value <- as.numeric(gsub(",", ".", res$value))
    res$name  <- name

    # Adding lon, lat, alt
    for (n in c("lon", "lat", "alt")) res[[n]] <- NA_real_
    for (i in seq_len(nrow(alt))) res$alt[res$date >= alt$date[i]] <- alt$alt[i]
    for (i in seq_len(nrow(alt))) {
        res$lon[res$date >= loc$date[i]] <- loc$lon[i]
        res$lat[res$date >= loc$date[i]] <- loc$lat[i]
    }

    # Returning data.frame
    invisible(res)
}
#x <- get_data(files[[50]])
#sapply(x, function(y) sum(is.na(y)))
#print(dim(x))

data <- lapply(files, get_data)
data <- do.call(rbind, data)
data <- na.omit(data)
head(data)
dim(data)
print(range(data$date))
barplot(table(as.integer(format(data$date, "%Y"))))


# Saving the data set
cat("Saving result into", RDSFILE, "\n")
saveRDS(data, RDSFILE)
