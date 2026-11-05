#' Solar elevation angle
#'
#' Low-cost solar position (NOAA / Astronomical Almanac approximation, accurate
#' to roughly 0.1 degrees, no refraction). Vectorised over all arguments. Used
#' to restrict radar profiles to darkness without needing extra packages.
#'
#' @param time POSIXct (UTC assumed if no time zone).
#' @param lon,lat Decimal degrees.
#' @return Solar elevation in degrees.
#' @export
sun_altitude <- function(time, lon, lat) {
  t <- as.numeric(as.POSIXct(time, tz = "UTC"))
  n <- t / 86400 + 2440587.5 - 2451545.0
  L <- (280.460 + 0.9856474 * n) %% 360
  g <- deg2rad((357.528 + 0.9856003 * n) %% 360)
  lambda <- deg2rad(L + 1.915 * sin(g) + 0.020 * sin(2 * g))
  eps <- deg2rad(23.439 - 0.0000004 * n)
  ra <- atan2(cos(eps) * sin(lambda), cos(lambda))
  dec <- asin(sin(eps) * sin(lambda))
  gmst <- (18.697374558 + 24.06570982441908 * n) %% 24
  ha <- deg2rad(gmst * 15 + lon) - ra
  phi <- deg2rad(lat)
  rad2deg(asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(ha)))
}

#' Approximate sunset time preceding a night
#'
#' Finds the time on `night_date` (UTC date of the evening) at which the sun
#' crosses `elev` degrees going down, by a 1-minute grid search between
#' 12:00 and 24:00 local solar time.
#'
#' @param night_date Date of the evening.
#' @param lon,lat Decimal degrees (scalars or vectors matching `night_date`).
#' @param elev Sun elevation defining sunset (default 0).
#' @return POSIXct (UTC).
#' @export
sunset_time <- function(night_date, lon, lat, elev = 0) {
  n <- max(length(night_date), length(lon), length(lat))
  night_date <- rep_len(as.Date(night_date), n)
  lon <- rep_len(lon, n); lat <- rep_len(lat, n)
  out <- rep(as.POSIXct(NA, tz = "UTC"), n)
  grid <- seq(0, 12 * 3600, by = 60)
  for (i in seq_len(n)) {
    noon <- as.POSIXct(night_date[i], tz = "UTC") + 12 * 3600 - lon[i] / 15 * 3600
    tt <- noon + grid
    a <- sun_altitude(tt, lon[i], lat[i])
    k <- which(a[-length(a)] >= elev & a[-1] < elev)[1]
    if (!is.na(k)) out[i] <- tt[k]
  }
  out
}
