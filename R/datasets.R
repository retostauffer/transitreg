#' Rain Ireland
#'
#' adsf
#'
#' @details
#' Station information:
#' * Country/county: Ireland, Clare
#' * Station: Ardnacrusha (Gen.Stn.No.2; ID 4011)
#' * Altitude (m): `28`
#' * Location (deg): `52.70639`N/`-8.61444`W
#'
#' @format A data frame containing xyz
#'
#' * `date`: The observations are valid from
#'    0900 UTC of the previous day until 0900 UTC of the date provided.
#' * `rain`: 24 hour precipitation sum in mm.
#'
#' For details regarding the observations and recording period please
#' see 'Rainfall Stations' on <https://www.met.ie/climate/what-we-measure>.
#'
#' @source Met Éireann, the Irish Meteorological Service,
#' licensed under the Creative Commons Attribution 4.0 International
#' (CC BY 4.0) license via. Available at
#' <https://www.met.ie/climate/available-data/historical-data>.
#'
#' @examples
#' data(rainIreland)
#' hist(sqrt(rainIreland$rain), freq = FALSE,
#'      main = "Daily precipitation sums\nArdnacrusha, Ireland",
#'      xlab = expression(sqrt(rain)), col = "steelblue")
#'
#' @docType data
#' @author Reto
"rainIreland"
