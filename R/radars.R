#' European weather radar sites
#'
#' Thin wrapper around [getRad::get_weather_radars()] that keeps one row per
#' radar identifier (the most recent metadata entry) and standard columns.
#'
#' @param source Passed to `getRad::get_weather_radars()`.
#' @return sf POINT data frame with `radar`, `country`, `location`, `lon`,
#'   `lat`, `startyear`, `finishyear`, `status`.
#' @export
radar_sites <- function(source = "opera") {
  rlang::check_installed("getRad", reason = "to download radar metadata.")
  r <- getRad::get_weather_radars(source = source)
  standardise_radars(r)
}

#' Standardise a radar metadata table
#'
#' Accepts getRad output or any data frame with `radar` and `longitude` /
#' `latitude` (or `lon` / `lat`) columns.
#' @param r Data frame or sf.
#' @export
standardise_radars <- function(r) {
  df <- if (inherits(r, "sf")) sf::st_drop_geometry(r) else r
  df <- tibble::as_tibble(as.data.frame(df))
  if (!"lon" %in% names(df)) df$lon <- df$longitude
  if (!"lat" %in% names(df)) df$lat <- df$latitude
  check_cols(df, c("radar", "lon", "lat"), "r")
  for (k in c("country", "location", "startyear", "finishyear", "status")) {
    if (!k %in% names(df)) df[[k]] <- NA
  }
  df <- df |>
    dplyr::filter(!is.na(.data$radar), is.finite(.data$lon), is.finite(.data$lat)) |>
    dplyr::arrange(.data$radar, dplyr::desc(.data$startyear)) |>
    dplyr::distinct(.data$radar, .keep_all = TRUE) |>
    dplyr::select("radar", "country", "location", "lon", "lat",
                  "startyear", "finishyear", "status")
  sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

# Sparse within-distance join returning (row_x, row_y, dist_km).
near_pairs <- function(x, y, max_km) {
  xp <- sf::st_transform(sf::st_geometry(x), crs_europe)
  yp <- sf::st_transform(sf::st_geometry(y), crs_europe)
  idx <- sf::st_is_within_distance(xp, yp, dist = max_km * 1000)
  ix <- rep(seq_along(idx), lengths(idx))
  iy <- unlist(idx)
  if (!length(ix)) return(tibble::tibble(ix = integer(), iy = integer(), dist_km = numeric()))
  d <- as.numeric(sf::st_distance(xp[ix], yp[iy], by_element = TRUE)) / 1000
  tibble::tibble(ix = ix, iy = iy, dist_km = d)
}

#' Match nightly bat segments to nearby radars
#'
#' A segment matches a radar when the straight line between consecutive fixes
#' passes within `max_km` of the radar. vol2bird profiles describe the air
#' within roughly 5 to 35 km of the antenna, so the default of 50 km allows
#' for Sigfox location error. Distance from the radar to the segment start is
#' also returned, for exposure models that need the bat to begin the night
#' within range.
#'
#' @param segments Output of [bat_nights()].
#' @param radars Output of [radar_sites()] / [standardise_radars()].
#' @param max_km Maximum segment-radar distance.
#' @return Tibble with `segment_id`, `radar`, `radar_dist_km`, `dist_start_km`.
#' @export
match_segments_to_radars <- function(segments, radars, max_km = 50) {
  p <- near_pairs(segments, radars, max_km)
  out <- tibble::tibble(
    segment_id = segments$segment_id[p$ix],
    radar = radars$radar[p$iy],
    radar_dist_km = p$dist_km
  )
  out$dist_start_km <- haversine_km(
    segments$lon_start[match(out$segment_id, segments$segment_id)],
    segments$lat_start[match(out$segment_id, segments$segment_id)],
    radars$lon[match(out$radar, radars$radar)],
    radars$lat[match(out$radar, radars$radar)]
  )
  out
}

#' Match in-flight fixes to radars
#'
#' @param fixes Output of [bat_flight_fixes()].
#' @param radars Radar sites.
#' @param max_km Maximum fix-radar distance (default 35 km, the default
#'   outer range used by vol2bird).
#' @return Tibble with `fix_id`, `radar`, `radar_dist_km`.
#' @export
match_fixes_to_radars <- function(fixes, radars, max_km = 35) {
  if (!"fix_id" %in% names(fixes)) fixes$fix_id <- seq_len(nrow(fixes))
  p <- near_pairs(fixes, radars, max_km)
  tibble::tibble(fix_id = fixes$fix_id[p$ix], radar = radars$radar[p$iy],
                 radar_dist_km = p$dist_km)
}

#' Radar-days of VPTS data needed for the matches
#'
#' @param segments,seg_matches Segments and their radar matches.
#' @param fixes,fix_matches Flight fixes and their radar matches (optional).
#' @return Tibble of unique `radar`, `date` (UTC) pairs to download.
#' @export
radar_requests <- function(segments, seg_matches, fixes = NULL, fix_matches = NULL) {
  s <- dplyr::left_join(seg_matches,
                        sf::st_drop_geometry(segments)[, c("segment_id", "t_start", "t_end")],
                        by = "segment_id")
  days <- purrr::map2(as.Date(s$t_start), as.Date(s$t_end), seq, by = "day")
  req <- tibble::tibble(radar = rep(s$radar, lengths(days)),
                        date = as.Date(unlist(days), origin = "1970-01-01"))
  if (!is.null(fixes) && !is.null(fix_matches) && nrow(fix_matches)) {
    if (!"fix_id" %in% names(fixes)) fixes$fix_id <- seq_len(nrow(fixes))
    f <- dplyr::left_join(fix_matches, sf::st_drop_geometry(fixes)[, c("fix_id", "time")],
                          by = "fix_id")
    req <- dplyr::bind_rows(req, tibble::tibble(radar = f$radar, date = as.Date(f$time)))
  }
  dplyr::distinct(req) |> dplyr::arrange(.data$radar, .data$date)
}
