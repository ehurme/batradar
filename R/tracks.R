# Pull id / time column names from a move2 object, with fallbacks for plain sf.
track_cols <- function(x, id_col = NULL, time_col = NULL) {
  id_col <- id_col %||% attr(x, "track_id_column") %||% "individual_local_identifier"
  time_col <- time_col %||% attr(x, "time_column") %||% "timestamp"
  list(id = id_col, time = time_col)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Drop geometry, add lon/lat from geometry, standardise id/time names.
flatten_track <- function(x, id_col = NULL, time_col = NULL) {
  cols <- track_cols(x, id_col, time_col)
  if (inherits(x, "sf")) {
    xy <- sf::st_coordinates(sf::st_transform(x, 4326))
    df <- sf::st_drop_geometry(x)
    df <- tibble::as_tibble(as.data.frame(df))
    df$.lon <- xy[, 1]; df$.lat <- xy[, 2]
  } else {
    df <- tibble::as_tibble(x)
    check_cols(df, c("lon", "lat"), "x")
    df$.lon <- df$lon; df$.lat <- df$lat
  }
  check_cols(df, c(cols$id, cols$time), "x")
  df$.id <- as.character(df[[cols$id]])
  df$.time <- as.POSIXct(df[[cols$time]], tz = "UTC")
  df
}

opt_col <- function(df, col, default = NA) {
  if (col %in% names(df)) df[[col]] else rep(default, nrow(df))
}

#' Nightly movement segments from daily bat positions
#'
#' Converts one-position-per-day tracks (e.g. `bats_daily`) into
#' consecutive-position segments, each describing what one bat did between two
#' fixes. Segments longer than `max_gap_h` are dropped because they can span
#' several nights and cannot be matched to a single radar night.
#'
#' A segment is flagged `migrated` when the displacement exceeds both
#' `min_km` and `err_mult` times the combined Sigfox location error, so that
#' location noise is not mistaken for migration.
#'
#' @param x A move2 / sf object (or data frame with `lon`, `lat`).
#' @param min_km Minimum displacement (km) to call a migration night.
#' @param err_mult Multiplier on combined location error (km).
#' @param max_gap_h Maximum time between fixes (hours).
#' @param err_col Column with location radius in metres.
#' @param id_col,time_col Override track id / time columns.
#' @return An sf data frame of LINESTRING segments (EPSG:4326).
#' @export
bat_nights <- function(x, min_km = 20, err_mult = 2, max_gap_h = 36,
                       err_col = "sigfox_computed_location_radius",
                       id_col = NULL, time_col = NULL) {
  df <- flatten_track(x, id_col, time_col)
  df$.err_km <- suppressWarnings(as.numeric(opt_col(df, err_col))) / 1000
  df <- df[is.finite(df$.lon) & is.finite(df$.lat) & !is.na(df$.time), ]
  df <- df[order(df$.id, df$.time), ]

  seg <- df |>
    dplyr::group_by(.data$.id) |>
    dplyr::mutate(
      t_start = dplyr::lag(.data$.time),
      lon_start = dplyr::lag(.data$.lon),
      lat_start = dplyr::lag(.data$.lat),
      err_start = dplyr::lag(.data$.err_km)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$t_start))

  seg <- seg |>
    dplyr::transmute(
      individual = .data$.id,
      species = opt_col(seg, "species"),
      tag_type = opt_col(seg, "tag_type"),
      sex = opt_col(seg, "sex"),
      season = opt_col(seg, "season"),
      t_start = .data$t_start,
      t_end = .data$.time,
      gap_h = as.numeric(difftime(.data$.time, .data$t_start, units = "hours")),
      lon_start = .data$lon_start, lat_start = .data$lat_start,
      lon_end = .data$.lon, lat_end = .data$.lat,
      dist_km = haversine_km(.data$lon_start, .data$lat_start, .data$.lon, .data$.lat),
      bearing = bearing_deg(.data$lon_start, .data$lat_start, .data$.lon, .data$.lat),
      err_km = sqrt(dplyr::coalesce(.data$err_start, 0)^2 +
                      dplyr::coalesce(.data$.err_km, 0)^2)
    ) |>
    dplyr::filter(.data$gap_h > 0, .data$gap_h <= max_gap_h)

  mid <- seg$t_start + (seg$t_end - seg$t_start) / 2
  seg$night_date <- solar_night_date(mid, (seg$lon_start + seg$lon_end) / 2)
  seg$yday <- as.integer(format(seg$night_date, "%j"))
  seg$migrated <- seg$dist_km >= pmax(min_km, err_mult * seg$err_km)
  seg$bearing[!seg$migrated] <- NA_real_
  seg$segment_id <- seq_len(nrow(seg))

  geoms <- lapply(seq_len(nrow(seg)), function(i) {
    sf::st_linestring(matrix(c(seg$lon_start[i], seg$lon_end[i],
                               seg$lat_start[i], seg$lat_end[i]), ncol = 2))
  })
  sf::st_sf(seg, geometry = sf::st_sfc(geoms, crs = 4326))
}

#' Night-time in-flight fixes with altitude
#'
#' Extracts fixes flagged as flight, taken in darkness, with a pressure
#' altitude. These are the fixes that can be placed in a specific radar
#' height bin at a specific time. A heading is added from the next fix of the
#' same individual when it follows within `heading_max_h` hours.
#'
#' @param x A move2 / sf object with all locations (e.g. `bats_loc`).
#' @param flight_col Logical/0-1 column flagging flight.
#' @param alt_col Altitude column (metres above sea level).
#' @param sun_elev Maximum solar elevation for a fix to count as night.
#' @param heading_max_h Max hours to next fix for a heading.
#' @param min_agl Optional minimum height above ground (m). Uncalibrated
#'   pressure altitudes can sit near or below ground for low or roosting
#'   bats; set e.g. `min_agl = 50` to keep only fixes clearly airborne.
#' @param err_col Column with location radius in metres.
#' @param id_col,time_col Override track id / time columns.
#' @return sf POINT data frame (EPSG:4326).
#' @export
bat_flight_fixes <- function(x, flight_col = "is_flight", alt_col = "altitude_m",
                             sun_elev = -6, heading_max_h = 3, min_agl = NULL,
                             err_col = "sigfox_computed_location_radius",
                             id_col = NULL, time_col = NULL) {
  df <- flatten_track(x, id_col, time_col)
  df <- df[is.finite(df$.lon) & is.finite(df$.lat) & !is.na(df$.time), ]
  df <- df[order(df$.id, df$.time), ]
  df <- df |>
    dplyr::group_by(.data$.id) |>
    dplyr::mutate(
      .dt_next = as.numeric(difftime(dplyr::lead(.data$.time), .data$.time, units = "hours")),
      .heading = bearing_deg(.data$.lon, .data$.lat, dplyr::lead(.data$.lon), dplyr::lead(.data$.lat)),
      .d_next = haversine_km(.data$.lon, .data$.lat, dplyr::lead(.data$.lon), dplyr::lead(.data$.lat))
    ) |>
    dplyr::ungroup()
  fl <- opt_col(df, flight_col, FALSE)
  fl <- as.logical(as.numeric(fl)) %in% TRUE
  alt <- suppressWarnings(as.numeric(opt_col(df, alt_col)))
  sun <- sun_altitude(df$.time, df$.lon, df$.lat)
  keep <- fl & is.finite(alt) & sun < sun_elev
  df <- df[keep, ]
  out <- tibble::tibble(
    individual = df$.id,
    species = opt_col(df, "species"),
    tag_type = opt_col(df, "tag_type"),
    time = df$.time,
    lon = df$.lon, lat = df$.lat,
    altitude_asl = alt[keep],
    ground_elev = suppressWarnings(as.numeric(opt_col(df, "elevation"))),
    err_km = suppressWarnings(as.numeric(opt_col(df, err_col))) / 1000,
    heading = ifelse(!is.na(df$.dt_next) & df$.dt_next <= heading_max_h &
                       df$.d_next > 2 * dplyr::coalesce(
                         suppressWarnings(as.numeric(opt_col(df, err_col))) / 1000, 0),
                     df$.heading, NA_real_),
    night_date = solar_night_date(df$.time, df$.lon)
  )
  out$altitude_agl <- out$altitude_asl - out$ground_elev
  if (!is.null(min_agl)) out <- out[!is.na(out$altitude_agl) & out$altitude_agl >= min_agl, ]
  out$fix_id <- seq_len(nrow(out))
  sf::st_as_sf(out, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}
