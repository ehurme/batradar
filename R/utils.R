#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"

# Equal-area projection used for all distance work (ETRS89 / LAEA Europe).
crs_europe <- 3035

deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

#' Bearing (degrees clockwise from north) between two lon/lat points
#'
#' Initial great-circle bearing, returned in `[0, 360)`.
#' @param lon1,lat1,lon2,lat2 Coordinates in decimal degrees.
#' @export
bearing_deg <- function(lon1, lat1, lon2, lat2) {
  p1 <- deg2rad(lat1); p2 <- deg2rad(lat2); dl <- deg2rad(lon2 - lon1)
  y <- sin(dl) * cos(p2)
  x <- cos(p1) * sin(p2) - sin(p1) * cos(p2) * cos(dl)
  (rad2deg(atan2(y, x)) + 360) %% 360
}

#' Great-circle distance in km
#' @inheritParams bearing_deg
#' @export
haversine_km <- function(lon1, lat1, lon2, lat2) {
  p1 <- deg2rad(lat1); p2 <- deg2rad(lat2)
  a <- sin((p2 - p1) / 2)^2 + cos(p1) * cos(p2) * sin(deg2rad(lon2 - lon1) / 2)^2
  6371.0088 * 2 * asin(pmin(1, sqrt(a)))
}

#' Smallest absolute angular difference between two directions (degrees)
#' @param a,b Directions in degrees.
#' @export
angle_diff <- function(a, b) {
  d <- abs(((a - b) + 180) %% 360 - 180)
  d
}

# Weighted circular mean (degrees) and mean resultant length.
circ_mean <- function(theta_deg, w = rep(1, length(theta_deg))) {
  ok <- is.finite(theta_deg) & is.finite(w) & w > 0
  if (!any(ok)) return(c(mean = NA_real_, rho = NA_real_))
  th <- deg2rad(theta_deg[ok]); w <- w[ok]
  s <- sum(w * sin(th)); c <- sum(w * cos(th))
  c(mean = (rad2deg(atan2(s, c)) + 360) %% 360, rho = sqrt(s^2 + c^2) / sum(w))
}

#' Area under the ROC curve (Mann-Whitney form)
#'
#' @param score Numeric predictor.
#' @param label Logical or 0/1 outcome.
#' @return AUC in `[0, 1]`, `NA` if either class is empty.
#' @export
auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]; label <- as.logical(label[ok])
  n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Night date anchored at local solar noon: a timestamp between solar noon on
# day d and solar noon on day d + 1 belongs to night d.
solar_night_date <- function(time, lon) {
  as.Date(as.POSIXct(time, tz = "UTC") + lon / 15 * 3600 - 12 * 3600, tz = "UTC")
}

check_cols <- function(df, cols, arg = deparse(substitute(df))) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    cli::cli_abort("{.arg {arg}} is missing column{?s} {.field {miss}}.")
  }
  invisible(TRUE)
}
