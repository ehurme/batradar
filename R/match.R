#' Place in-flight bat fixes inside radar profiles
#'
#' For each fix-radar pair, finds the nearest radar profile in time (within
#' `max_dt_min`) and the height bin containing the bat's pressure altitude,
#' then reports what the radar saw there and how that compares to the rest
#' of the column.
#'
#' @param fixes Output of [bat_flight_fixes()] (with `fix_id`).
#' @param fix_matches Output of [match_fixes_to_radars()].
#' @param v Output of [vpts_standardise()].
#' @param max_dt_min Maximum time offset between fix and profile (minutes).
#' @param min_agl Clutter floor above antenna (m), as in [vpts_nightly()].
#' @return Tibble, one row per matched fix-radar pair, with `dens_at_bat`,
#'   `col_vid`, `dens_rank_pct` (percentile of the bat's bin among bins in
#'   that profile), `cum_below` (share of column density below the bat),
#'   `height_mean_radar`, `dd_at_bat`, `ff_at_bat`, `heading_diff`.
#' @export
match_fixes_to_vpts <- function(fixes, fix_matches, v, max_dt_min = 15, min_agl = 200) {
  if (!"fix_id" %in% names(fixes)) fixes$fix_id <- seq_len(nrow(fixes))
  fx <- sf::st_drop_geometry(fixes)
  fm <- dplyr::left_join(fix_matches, fx[, c("fix_id", "individual", "species", "time",
                                             "altitude_asl", "heading")],
                         by = "fix_id")
  fm <- fm[fm$radar %in% unique(v$radar), ]
  if (!nrow(fm)) return(tibble::tibble())
  v <- v[v$height_agl_radar >= min_agl, ]
  # nearest profile time per pair
  pt <- dplyr::distinct(v, .data$radar, .data$datetime) |> dplyr::arrange(.data$radar, .data$datetime)
  fm$profile_time <- as.POSIXct(NA, tz = "UTC"); fm$dt_min <- NA_real_
  for (r in unique(fm$radar)) {
    tt <- as.numeric(pt$datetime[pt$radar == r]); ii <- which(fm$radar == r)
    k <- findInterval(as.numeric(fm$time[ii]), tt, all.inside = TRUE)
    k2 <- pmin(k + 1, length(tt))
    d1 <- abs(tt[k] - as.numeric(fm$time[ii])); d2 <- abs(tt[k2] - as.numeric(fm$time[ii]))
    best <- ifelse(d1 <= d2, k, k2)
    fm$profile_time[ii] <- as.POSIXct(tt[best], origin = "1970-01-01", tz = "UTC")
    fm$dt_min[ii] <- pmin(d1, d2) / 60
  }
  fm <- fm[fm$dt_min <= max_dt_min, ]
  if (!nrow(fm)) return(tibble::tibble())
  fm$pair <- seq_len(nrow(fm))
  vv <- dplyr::inner_join(
    fm[, c("pair", "radar", "profile_time", "altitude_asl")],
    v[, c("radar", "datetime", "height_asl", "bin_km", "dens", "dd", "ff")],
    by = c("radar", "profile_time" = "datetime"), relationship = "many-to-many")
  vv$d0 <- dplyr::coalesce(vv$dens, 0)
  vv$w <- vv$d0 * vv$bin_km
  vv$hit <- vv$altitude_asl >= vv$height_asl & vv$altitude_asl < vv$height_asl + vv$bin_km * 1000
  st <- vv |>
    dplyr::group_by(.data$pair) |>
    dplyr::summarise(
      col_vid = sum(.data$w),
      hit_n = sum(.data$hit),
      bin_height_asl = if (any(.data$hit)) .data$height_asl[.data$hit][1] else NA_real_,
      dens_at_bat = if (any(.data$hit)) .data$dens[.data$hit][1] else NA_real_,
      dens_rank_pct = if (any(.data$hit)) mean(.data$d0 <= .data$d0[.data$hit][1]) else NA_real_,
      cum_below = if (any(.data$hit) && sum(.data$w) > 0)
        sum(.data$w[.data$height_asl < .data$height_asl[.data$hit][1]]) / sum(.data$w) else NA_real_,
      height_mean_radar = if (sum(.data$w) > 0)
        sum(.data$w * (.data$height_asl + .data$bin_km * 500)) / sum(.data$w) else NA_real_,
      dd_at_bat = if (any(.data$hit)) .data$dd[.data$hit][1] else NA_real_,
      ff_at_bat = if (any(.data$hit)) .data$ff[.data$hit][1] else NA_real_,
      .groups = "drop")
  out <- dplyr::left_join(fm, st, by = "pair")
  out$in_profile <- !is.na(out$bin_height_asl)
  out$heading_diff <- angle_diff(out$heading, out$dd_at_bat)
  out$hit_n <- NULL; out$pair <- NULL
  tibble::as_tibble(out)
}
