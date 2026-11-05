# Synthetic fixtures for unit tests only (not real radar or tag data).
fake_vpts <- function(radar = "xxtst", lon = 9, lat = 47.7,
                      t0 = as.POSIXct("2024-09-10 12:00", tz = "UTC"),
                      t1 = as.POSIXct("2024-09-11 12:00", tz = "UTC"),
                      dens_fun = function(h, t) ifelse(h < 1000, 10, 0),
                      u = 0, v = 10) {
  times <- seq(t0, t1, by = 300)
  heights <- seq(0, 3000, by = 200)
  g <- expand.grid(height = heights, datetime = times)
  tibble::tibble(
    radar = radar, datetime = g$datetime, height = g$height,
    height_reference = "sea", u = u, v = v, w = 0,
    ff = sqrt(u^2 + v^2), dd = (atan2(u, v) * 180 / pi + 360) %% 360,
    dens = dens_fun(g$height, g$datetime), eta = 11 * dens_fun(g$height, g$datetime),
    radar_longitude = lon, radar_latitude = lat, radar_height = 0)
}

fake_track <- function() {
  df <- data.frame(
    individual_local_identifier = rep(c("a", "b"), each = 3),
    timestamp = as.POSIXct(c("2024-09-10 12:00", "2024-09-11 12:00", "2024-09-12 12:00",
                             "2024-09-10 12:00", "2024-09-11 12:00", "2024-09-12 12:00"), tz = "UTC"),
    lon = c(9, 9, 9.0, 9, 9.01, 9.02),
    lat = c(47.0, 47.7, 48.6, 47.7, 47.7, 47.7),
    sigfox_computed_location_radius = 5000,
    species = "Nyctalus noctula", tag_type = "Nanofox")
  x <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  attr(x, "track_id_column") <- "individual_local_identifier"
  attr(x, "time_column") <- "timestamp"
  x
}

fake_radars <- function() {
  standardise_radars(data.frame(radar = c("xxtst", "xxfar"), lon = c(9, 20), lat = c(47.7, 60)))
}
