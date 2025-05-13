
#' Weather Data Shannon Airport (Ireland)
#'
#' This data set contains daily meteorological observations from Shannon Airport,
#' located on the west coast of Ireland, covering the period from September 1,
#' 1945, to January 31, 2025 (29,008 days, approximately 79.5 years). The data,
#' provided by the Irish Meteorological Service (Met Éireann), includes daily
#' minimum and maximum temperatures, precipitation amounts, and sunshine duration.
#'
#' @details
#' **Station information**
#'
#' * Country/county: Ireland, Clare
#' * Station: Shannon Airport (ID 518)
#' * Altitude: 15 m a.m.s.l.
#' * Location: 52.69028°N, 8.91806°W
#'
#' For more details, see the [weather observations website](https://wow.met.ie/).
#'
#' @format A data frame containing historical meteorological observations from
#' Shannon Airport (Clare County, Ireland) spanning over 79 years, without gaps
#' or missing values. Observations are valid from 0000 UTC to 0000 UTC.
#'
#' The data frame contains the following variables:
#'
#' * `date`: Date of observation  
#' * `maxtp`: Daily maximum temperature \[°C\]
#' * `mintp`: Daily minimum temperature \[°C\]
#' * `rain`: Daily precipitation sum \[mm\]
#' * `sun`: Sunshine duration \[hours\]
#'
#' For more details about the observations and recording period,  
#' visit <https://www.met.ie/climate/what-we-measure>.
#'
#' @source Met Éireann, the Irish Meteorological Service,
#' licensed under the Creative Commons Attribution 4.0 International
#' (CC BY 4.0) license via. Available at
#' <https://www.met.ie/climate/available-data/historical-data>.
#'
#' @examples
#' data(Shannon)
#' summary(Shannon)
#'
#' @docType data
#' @author Reto
"Shannon"

