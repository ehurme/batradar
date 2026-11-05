#' Download VPTS data for radar-days, with a per radar-month cache
#'
#' Calls [getRad::get_vpts()] once per radar and month and caches each result
#' as an `.rds` file so the pipeline can be re-run without downloading again.
#' Failures (no data for that radar-month, network errors) are recorded rather
#' than stopping the run.
#'
#' @param requests Output of [radar_requests()].
#' @param source VPTS source passed to `getRad::get_vpts()` (`"baltrad"`
#'   covers 2016 onwards for most OPERA radars, `"uva"` and `"ecog-04003"` are
#'   older archives).
#' @param cache_dir Directory for cached radar-month files.
#' @param coverage Optional output of [getRad::get_vpts_coverage()]; when
#'   given, requests for radar-days not in the archive are skipped.
#' @param quiet Suppress progress messages.
#' @return Tibble of VPTS rows (vpts-csv columns) with attribute `failures`.
#' @export
fetch_vpts <- function(requests, source = "baltrad", cache_dir = "cache/vpts",
                       coverage = NULL, quiet = FALSE) {
  rlang::check_installed("getRad", reason = "to download VPTS data.")
  check_cols(requests, c("radar", "date"))
  req <- requests
  if (!is.null(coverage)) {
    cov <- tibble::as_tibble(coverage)
    if (all(c("radar", "date") %in% names(cov))) {
      cov$date <- as.Date(cov$date)
      req <- dplyr::semi_join(req, cov, by = c("radar", "date"))
    }
  }
  req$month <- format(req$date, "%Y%m")
  jobs <- req |>
    dplyr::group_by(.data$radar, .data$month) |>
    dplyr::summarise(d0 = min(.data$date), d1 = max(.data$date), .groups = "drop")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  failures <- list()
  res <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    j <- jobs[i, ]
    f <- file.path(cache_dir, sprintf("%s_%s_%s.rds", source, j$radar, j$month))
    if (file.exists(f)) { res[[i]] <- readRDS(f); next }
    if (!quiet) cli::cli_inform("VPTS {j$radar} {j$month} ({i}/{nrow(jobs)})")
    dt <- if (j$d0 == j$d1) j$d0 else c(j$d0, j$d1 + 1)
    v <- tryCatch(
      getRad::get_vpts(radar = j$radar, datetime = dt, source = source,
                       return_type = "tibble"),
      error = function(e) e
    )
    if (inherits(v, "error")) {
      failures[[length(failures) + 1]] <- tibble::tibble(
        radar = j$radar, month = j$month, error = conditionMessage(v))
      next
    }
    v <- tibble::as_tibble(v)
    saveRDS(v, f)
    res[[i]] <- v
  }
  out <- dplyr::bind_rows(res)
  attr(out, "failures") <- dplyr::bind_rows(failures)
  out
}

#' Read cached VPTS files
#' @param cache_dir Directory used by [fetch_vpts()].
#' @export
read_vpts_cache <- function(cache_dir = "cache/vpts") {
  fs <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  dplyr::bind_rows(lapply(fs, readRDS))
}

#' Standardise VPTS rows
#'
#' Ensures POSIXct `datetime`, adds `height_asl` (bin lower bound above sea
#' level), `bin_km` (bin width), `dt_h` (profile interval, hours), and
#' `height_agl_radar` (height above the antenna). Rows from `height_reference
#' == "antenna"` are shifted by `radar_height`.
#'
#' @param v VPTS tibble (vpts-csv columns).
#' @return Tibble.
#' @export
vpts_standardise <- function(v) {
  check_cols(v, c("radar", "datetime", "height", "dens", "u", "v"), "v")
  v <- tibble::as_tibble(v)
  if (!inherits(v$datetime, "POSIXct")) {
    v$datetime <- as.POSIXct(v$datetime, tz = "UTC", tryFormats = c(
      "%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS"))
  }
  attr(v$datetime, "tzone") <- "UTC"
  if (!"radar_height" %in% names(v)) v$radar_height <- 0
  if (!"height_reference" %in% names(v)) v$height_reference <- "sea"
  v$height_asl <- ifelse(v$height_reference %in% "antenna",
                         v$height + v$radar_height, v$height)
  v$height_agl_radar <- v$height_asl - v$radar_height
  if (!"ff" %in% names(v)) v$ff <- sqrt(v$u^2 + v$v^2)
  if (!"dd" %in% names(v)) v$dd <- (rad2deg(atan2(v$u, v$v)) + 360) %% 360
  steps <- v |>
    dplyr::group_by(.data$radar) |>
    dplyr::summarise(
      bin_km = {
        h <- sort(unique(.data$height)); d <- diff(h)
        if (length(d)) min(d[d > 0]) / 1000 else 0.2
      },
      dt_h = {
        t <- sort(unique(.data$datetime))
        d <- as.numeric(diff(t), units = "hours")
        if (length(d)) min(stats::median(d), 1) else 5 / 60
      },
      .groups = "drop")
  v <- dplyr::left_join(dplyr::select(v, -dplyr::any_of(c("bin_km", "dt_h"))),
                        steps, by = "radar")
  dens0 <- dplyr::coalesce(v$dens, NA_real_)
  v$mtr_bin <- ifelse(!is.na(dens0) & dens0 == 0, 0, dens0 * v$ff * 3.6 * v$bin_km)
  v
}

# Rows in darkness, above the clutter floor, inside the height band.
vpts_filter <- function(v, sun_elev = -6, min_agl = 200, height_band = c(-Inf, Inf),
                        radar_lonlat = NULL) {
  lon <- if ("radar_longitude" %in% names(v)) v$radar_longitude else NA
  lat <- if ("radar_latitude" %in% names(v)) v$radar_latitude else NA
  if (!is.null(radar_lonlat)) {
    k <- match(v$radar, radar_lonlat$radar)
    lon <- dplyr::coalesce(lon, radar_lonlat$lon[k])
    lat <- dplyr::coalesce(lat, radar_lonlat$lat[k])
  }
  if (all(is.na(lon))) cli::cli_abort("Radar coordinates are needed (radar_longitude/latitude or {.arg radar_lonlat}).")
  ut <- unique(data.frame(radar = v$radar, datetime = v$datetime, lon = lon, lat = lat))
  ut$sun <- sun_altitude(ut$datetime, ut$lon, ut$lat)
  ut$night_date <- solar_night_date(ut$datetime, ut$lon)
  v <- dplyr::left_join(v, ut[, c("radar", "datetime", "sun", "night_date", "lon", "lat")],
                        by = c("radar", "datetime"))
  v[v$sun < sun_elev & v$height_agl_radar >= min_agl &
      v$height_asl >= height_band[1] & v$height_asl < height_band[2], ]
}

# Per-profile additive totals. Every period summary is built from these sums,
# so windows can be aggregated with cumulative sums instead of re-scanning rows.
profile_table <- function(vf) {
  w <- vf$dens * vf$bin_km
  okw <- is.finite(w) & w > 0 & is.finite(vf$dd)
  w0 <- ifelse(okw, w, 0)
  th <- deg2rad(ifelse(okw, vf$dd, 0))
  vf$.w <- w0; vf$.s <- w0 * sin(th); vf$.c <- w0 * cos(th)
  vf$.wff <- w0 * ifelse(okw, vf$ff, 0)
  vf$.wh <- w0 * (vf$height_asl + vf$bin_km * 500)
  vf$.vid <- ifelse(is.finite(w), w, 0)
  vf |>
    dplyr::group_by(.data$radar, .data$datetime, .data$night_date) |>
    dplyr::summarise(
      ok = any(!is.na(.data$dens)),
      mtr = sum(.data$mtr_bin, na.rm = TRUE),
      vid = sum(.data$.vid),
      dt_h = dplyr::first(.data$dt_h),
      sw = sum(.data$.w), ss = sum(.data$.s), sc = sum(.data$.c),
      swff = sum(.data$.wff), swh = sum(.data$.wh),
      .groups = "drop") |>
    dplyr::arrange(.data$radar, .data$datetime)
}

# Turn summed quantities into the reported metrics.
finish_period <- function(d) {
  d$coverage <- d$n_ok / d$n_profiles
  d$mt_scaled <- ifelse(d$coverage > 0, d$mt / d$coverage, NA_real_)
  d$mtr_mean <- ifelse(d$n_ok > 0, d$mtr_ok / d$n_ok, NA_real_)
  d$vid_mean <- ifelse(d$n_ok > 0, d$vid_ok / d$n_ok, NA_real_)
  d$dir_mean <- ifelse(d$sw > 0, (rad2deg(atan2(d$ss, d$sc)) + 360) %% 360, NA_real_)
  d$dir_rho <- ifelse(d$sw > 0, sqrt(d$ss^2 + d$sc^2) / d$sw, NA_real_)
  d$ff_mean <- ifelse(d$sw > 0, d$swff / d$sw, NA_real_)
  d$height_mean <- ifelse(d$sw > 0, d$swh / d$sw, NA_real_)
  d[, setdiff(names(d), c("n_ok", "mtr_ok", "vid_ok", "sw", "ss", "sc", "swff", "swh"))]
}

sum_cols <- function(prof) {
  data.frame(n_profiles = 1, n_ok = as.numeric(prof$ok), hours = prof$dt_h,
             mt = prof$mtr * prof$dt_h, mtr_ok = ifelse(prof$ok, prof$mtr, 0),
             vid_ok = ifelse(prof$ok, prof$vid, 0), sw = prof$sw, ss = prof$ss,
             sc = prof$sc, swff = prof$swff, swh = prof$swh)
}

# Sum profile totals over arbitrary time windows per radar.
# windows: data frame with radar, t_start, t_end (+ any key columns).
window_sums <- function(prof, windows) {
  out <- vector("list", 0)
  for (r in intersect(unique(windows$radar), unique(prof$radar))) {
    pr <- prof[prof$radar == r, ]
    wr <- windows[windows$radar == r, ]
    cs <- rbind(0, apply(as.matrix(sum_cols(pr)), 2, cumsum))
    tt <- as.numeric(pr$datetime)
    i0 <- findInterval(as.numeric(wr$t_start), tt, left.open = TRUE)  # profiles < t_start
    i1 <- findInterval(as.numeric(wr$t_end), tt)                     # profiles <= t_end
    tot <- cs[i1 + 1, , drop = FALSE] - cs[i0 + 1, , drop = FALSE]
    out[[r]] <- dplyr::bind_cols(tibble::as_tibble(wr), tibble::as_tibble(tot))
  }
  res <- dplyr::bind_rows(out)
  if (!nrow(res)) return(res)
  res <- res[res$n_profiles > 0, ]
  finish_period(res)
}

#' Radar migration per radar-night
#'
#' Integrates standardised VPTS over darkness for each radar and night
#' (night dates anchored at local solar noon, as for the tags).
#'
#' @param v Output of [vpts_standardise()].
#' @param sun_elev Solar elevation below which a profile counts as night.
#' @param min_agl Lowest height above the antenna to include (m), to avoid
#'   ground clutter.
#' @param height_band Height band (m above sea level) to integrate over.
#' @param radar_lonlat Radar coordinates if VPTS lacks `radar_longitude`.
#' @return Tibble with one row per `radar`, `night_date`: `mt` (migration
#'   traffic, animals/km), `mt_scaled` (corrected for missing profiles),
#'   `mtr_mean`, `vid_mean` (animals/km2), `dir_mean`, `dir_rho`, `ff_mean`,
#'   `height_mean`, `coverage`.
#' @export
vpts_nightly <- function(v, sun_elev = -6, min_agl = 200, height_band = c(-Inf, Inf),
                         radar_lonlat = NULL) {
  vf <- vpts_filter(v, sun_elev, min_agl, height_band, radar_lonlat)
  prof <- profile_table(vf)
  sc <- sum_cols(prof)
  sc$radar <- prof$radar; sc$night_date <- prof$night_date
  out <- sc |>
    dplyr::group_by(.data$radar, .data$night_date) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop")
  out$mtr_max <- as.numeric(tapply(prof$mtr, paste(prof$radar, prof$night_date), max)[
    paste(out$radar, out$night_date)])
  out <- finish_period(out)
  out$yday <- as.integer(format(out$night_date, "%j"))
  out
}

#' Radar migration during each bat segment
#'
#' Same metrics as [vpts_nightly()], but integrated only over the dark part
#' of each segment's `t_start`..`t_end` window at the matched radar.
#'
#' @param v Output of [vpts_standardise()].
#' @param segments Output of [bat_nights()].
#' @param seg_matches Output of [match_segments_to_radars()].
#' @inheritParams vpts_nightly
#' @return Tibble keyed by `segment_id`, `radar`.
#' @export
vpts_for_segments <- function(v, segments, seg_matches, sun_elev = -6, min_agl = 200,
                              height_band = c(-Inf, Inf), radar_lonlat = NULL) {
  vf <- vpts_filter(v, sun_elev, min_agl, height_band, radar_lonlat)
  win <- dplyr::left_join(seg_matches[, c("segment_id", "radar")],
                          sf::st_drop_geometry(segments)[, c("segment_id", "t_start", "t_end")],
                          by = "segment_id")
  res <- window_sums(profile_table(vf), win)
  if (!nrow(res)) return(tibble::tibble())
  res[, setdiff(names(res), c("t_start", "t_end"))]
}
