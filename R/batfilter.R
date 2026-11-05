#' Radar cells (radar x time x height) in darkness, with bat-relevant features
#'
#' @param v Output of [vpts_standardise()].
#' @param wind Optional tibble with `radar`, `datetime`, `height_asl`,
#'   `u_wind`, `v_wind` (m/s), e.g. ERA5 interpolated to radar bins. When
#'   given, airspeed and heading are computed by subtracting wind from
#'   ground velocity.
#' @inheritParams vpts_nightly
#' @return Tibble of cells with `hrs_since_sunset`, `yday`, `rcs`, and
#'   (with wind) `airspeed`, `heading_air`.
#' @export
vpts_cells <- function(v, wind = NULL, sun_elev = -6, min_agl = 200, radar_lonlat = NULL) {
  vf <- vpts_filter(v, sun_elev, min_agl, radar_lonlat = radar_lonlat)
  rn <- dplyr::distinct(vf, .data$radar, .data$night_date, .data$lon, .data$lat)
  rn$sunset <- sunset_time(rn$night_date, rn$lon, rn$lat)
  vf <- dplyr::left_join(vf, rn[, c("radar", "night_date", "sunset")],
                         by = c("radar", "night_date"))
  vf$hrs_since_sunset <- as.numeric(difftime(vf$datetime, vf$sunset, units = "hours"))
  vf$yday <- as.integer(format(vf$night_date, "%j"))
  vf$rcs_est <- if ("eta" %in% names(vf)) ifelse(vf$dens > 0, vf$eta / vf$dens, NA_real_) else NA_real_
  if (!is.null(wind)) {
    check_cols(wind, c("radar", "datetime", "height_asl", "u_wind", "v_wind"), "wind")
    vf <- dplyr::left_join(vf, wind, by = c("radar", "datetime", "height_asl"))
    ua <- vf$u - vf$u_wind; va <- vf$v - vf$v_wind
    vf$airspeed <- sqrt(ua^2 + va^2)
    vf$heading_air <- (rad2deg(atan2(ua, va)) + 360) %% 360
  }
  vf
}

#' Flight envelope of tagged bats
#'
#' Summarises where and when tagged bats flew, per group, as quantile ranges
#' of altitude above the radar antenna-equivalent (above ground), hours since
#' sunset, and day of year. This is the tags' own description of bat flight,
#' used as a first, transparent "bat filter" on radar data.
#'
#' @param fixes Output of [bat_flight_fixes()].
#' @param by Grouping columns (default species).
#' @param probs Lower and upper quantiles.
#' @return Tibble of ranges per group.
#' @export
tag_envelope <- function(fixes, by = "species", probs = c(0.1, 0.9)) {
  fx <- sf::st_drop_geometry(fixes)
  fx$sunset <- sunset_time(fx$night_date, fx$lon, fx$lat)
  fx$hrs_since_sunset <- as.numeric(difftime(fx$time, fx$sunset, units = "hours"))
  fx$yday <- as.integer(format(fx$night_date, "%j"))
  fx |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_fixes = dplyr::n(),
      n_ind = dplyr::n_distinct(.data$individual),
      agl_lo = stats::quantile(.data$altitude_agl, probs[1], na.rm = TRUE),
      agl_hi = stats::quantile(.data$altitude_agl, probs[2], na.rm = TRUE),
      agl_median = stats::median(.data$altitude_agl, na.rm = TRUE),
      hss_lo = stats::quantile(.data$hrs_since_sunset, probs[1], na.rm = TRUE),
      hss_hi = stats::quantile(.data$hrs_since_sunset, probs[2], na.rm = TRUE),
      .groups = "drop")
}

#' Restrict radar cells to a tag flight envelope
#'
#' @param cells Output of [vpts_cells()].
#' @param env One row of [tag_envelope()] output.
#' @return Logical vector, `TRUE` for cells inside the envelope.
#' @export
in_envelope <- function(cells, env) {
  stopifnot(nrow(env) == 1)
  cells$height_agl_radar >= env$agl_lo & cells$height_agl_radar <= env$agl_hi &
    cells$hrs_since_sunset >= env$hss_lo & cells$hrs_since_sunset <= env$hss_hi
}

#' Label radar cells occupied by tagged bats
#'
#' Marks each cell (radar, profile time, height bin) that contains a matched
#' in-flight fix as `present = TRUE`, and samples `n_bg` background cells per
#' presence from the same radar-night, so the comparison asks where in the
#' radar column the tagged bats were, given that night.
#'
#' @param cells Output of [vpts_cells()].
#' @param fix_vpts Output of [match_fixes_to_vpts()].
#' @param n_bg Background cells per presence.
#' @param seed Random seed.
#' @return Tibble of presence and background cells.
#' @export
label_cells <- function(cells, fix_vpts, n_bg = 20, seed = 1) {
  fv <- fix_vpts[fix_vpts$in_profile %in% TRUE, ]
  if (!nrow(fv)) cli::cli_abort("No fixes fall inside radar profiles.")
  key_c <- paste(cells$radar, as.numeric(cells$datetime), cells$height_asl)
  pres_keys <- unique(paste(fv$radar, as.numeric(fv$profile_time), fv$bin_height_asl))
  pres <- cells[key_c %in% pres_keys, ]
  pres$present <- TRUE
  set.seed(seed)
  rn <- unique(paste(pres$radar, pres$night_date))
  pool <- cells[paste(cells$radar, cells$night_date) %in% rn & !key_c %in% pres_keys, ]
  np <- table(paste(pres$radar, pres$night_date))
  pool$.g <- paste(pool$radar, pool$night_date)
  pool$.r <- stats::runif(nrow(pool))
  pool <- pool[order(pool$.g, pool$.r), ]
  pool$.k <- stats::ave(seq_len(nrow(pool)), pool$.g, FUN = seq_along)
  bg <- pool[pool$.k <= as.numeric(np[pool$.g]) * n_bg, ]
  bg$.g <- NULL; bg$.r <- NULL; bg$.k <- NULL
  bg$present <- FALSE
  dplyr::bind_rows(pres, bg)
}

#' Fit and cross-validate a radar bat-cell model
#'
#' Binomial GAM of presence vs background radar cells. Because presence
#' cells are defined by where tagged bats flew, a model can score well just
#' by learning the bats' height and timing. The function therefore also fits
#' a null model (`formula_null`, height and timing only) on the same folds
#' and reports how much the radar-signal terms add (`auc_gain`). Only that
#' gain is evidence that the radar signal itself carries a bat signature.
#'
#' Cross-validation is grouped by `group` (default radar-night) so test
#' folds contain nights the model has not seen.
#'
#' @param labelled Output of [label_cells()].
#' @param formula Full model. Add `s(airspeed)` when wind is available.
#' @param formula_null Envelope-only model for comparison.
#' @param k_folds Number of folds.
#' @param group Column defining CV groups.
#' @return List with `model` (full model on all data), `cv` (fold AUCs for
#'   both models), `auc_full`, `auc_null`, `auc_gain`.
#' @export
fit_bat_filter <- function(labelled,
                           formula = present ~ s(height_agl_radar, k = 6) +
                             s(hrs_since_sunset, k = 6) + s(yday, k = 6) +
                             s(ff, k = 5) + s(log1p(dens), k = 5),
                           formula_null = present ~ s(height_agl_radar, k = 6) +
                             s(hrs_since_sunset, k = 6) + s(yday, k = 6),
                           k_folds = 5, group = NULL) {
  d <- labelled
  vars <- unique(c(all.vars(formula), all.vars(formula_null)))
  d <- d[stats::complete.cases(d[, vars]), ]
  g <- if (is.null(group)) paste(d$radar, d$night_date) else d[[group]]
  ug <- unique(g)
  set.seed(42)
  fold_of_group <- sample(rep_len(seq_len(k_folds), length(ug)))
  fold <- fold_of_group[match(g, ug)]
  fit1 <- function(f, dat) mgcv::bam(f, data = dat, family = stats::binomial(),
                                     method = "fREML", discrete = TRUE)
  cv <- purrr::map_dfr(seq_len(k_folds), function(k) {
    tr <- d[fold != k, ]; te <- d[fold == k, ]
    if (!any(te$present) || !any(tr$present)) return(NULL)
    pf <- stats::predict(fit1(formula, tr), te, type = "response")
    pn <- stats::predict(fit1(formula_null, tr), te, type = "response")
    tibble::tibble(fold = k, n_test = nrow(te), n_pres = sum(te$present),
                   auc_full = auc(pf, te$present), auc_null = auc(pn, te$present))
  })
  list(model = fit1(formula, d), cv = cv,
       auc_full = mean(cv$auc_full, na.rm = TRUE),
       auc_null = mean(cv$auc_null, na.rm = TRUE),
       auc_gain = mean(cv$auc_full - cv$auc_null, na.rm = TRUE))
}
