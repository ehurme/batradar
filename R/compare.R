#' Pair bat segments with radar migration at the nearest matched radar
#'
#' @param segments Output of [bat_nights()].
#' @param seg_matches Output of [match_segments_to_radars()].
#' @param seg_radar Output of [vpts_for_segments()].
#' @param min_coverage Minimum share of dark profiles with data.
#' @return Tibble, one row per segment (nearest radar with data), with
#'   `log_mtr` (log1p of the mean migration traffic rate over the dark part of
#'   the segment window, animals/km/h; a rate so that windows of different
#'   length are comparable), `dir_diff` (bat bearing vs radar mean direction,
#'   migrating segments only).
#' @export
pair_segments <- function(segments, seg_matches, seg_radar, min_coverage = 0.5) {
  s <- sf::st_drop_geometry(segments)
  p <- seg_matches |>
    dplyr::inner_join(seg_radar, by = c("segment_id", "radar")) |>
    dplyr::filter(.data$coverage >= min_coverage) |>
    dplyr::group_by(.data$segment_id) |>
    dplyr::slice_min(.data$radar_dist_km, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::inner_join(s, by = "segment_id")
  p$log_mtr <- log1p(p$mtr_mean)
  p$dir_diff <- ifelse(p$migrated, angle_diff(p$bearing, p$dir_mean), NA_real_)
  p
}

#' Agreement between tag-detected migration and radar migration
#'
#' Night-level summary, overall and by group: how well radar migration
#' traffic separates nights on which tagged bats migrated past the radar from
#' nights they stayed, and how closely directions agree on migration nights.
#'
#' @param paired Output of [pair_segments()].
#' @param by Grouping columns.
#' @return Tibble with `n_seg`, `n_migr`, `auc_mtr` (AUC of radar traffic rate for
#'   bat migration), `rho_spearman` (bat displacement vs radar traffic),
#'   `dir_diff_median`, `prop_dir_within_45`, `dir_n`.
#' @export
night_agreement <- function(paired, by = NULL) {
  f <- function(d) {
    tibble::tibble(
      n_seg = nrow(d),
      n_ind = dplyr::n_distinct(d$individual),
      n_radar = dplyr::n_distinct(d$radar),
      n_migr = sum(d$migrated),
      auc_mtr = auc(d$log_mtr, d$migrated),
      rho_spearman = suppressWarnings(stats::cor(d$dist_km, d$log_mtr,
                                                 method = "spearman", use = "complete.obs")),
      dir_n = sum(is.finite(d$dir_diff)),
      dir_diff_median = stats::median(d$dir_diff, na.rm = TRUE),
      prop_dir_within_45 = mean(d$dir_diff <= 45, na.rm = TRUE)
    )
  }
  if (is.null(by)) return(f(paired))
  paired |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::group_modify(~ f(.x)) |>
    dplyr::ungroup()
}

#' Does radar traffic predict whether a tagged bat migrates?
#'
#' Binomial GAM of tag migration on radar traffic at the matched radar, with
#' individual random effects and a seasonal smooth.
#'
#' @param paired Output of [pair_segments()].
#' @param formula Model formula.
#' @return A fitted `mgcv::bam` (fast REML; individual random effects make
#'   `gam()` slow with hundreds of bats).
#' @export
fit_detection_gam <- function(paired,
                              formula = migrated ~ s(log_mtr, k = 5) + s(yday, k = 6) +
                                s(individual, bs = "re")) {
  d <- paired
  d$individual <- factor(d$individual)
  d$species <- factor(d$species)
  mgcv::bam(formula, data = d, family = stats::binomial(), method = "fREML", discrete = TRUE)
}

#' Fit the same smooth structure to tag and radar migration
#'
#' Builds comparable GAMs: tag migration (binomial, per bat-segment) and
#' radar migration (Gaussian on `log1p(mt_scaled)`, per radar-night), with
#' identical smooths of each covariate. Covariates must be columns of both
#' tables (e.g. `yday`, or weather joined from ERA5 such as wind support,
#' temperature, pressure change).
#'
#' @param bat_df Tag table, e.g. [pair_segments()] output or all segments.
#' @param radar_df Radar table, e.g. [vpts_nightly()] output restricted to
#'   matched radars and the tagging season.
#' @param covariates Covariate names.
#' @param k Basis dimension for each smooth.
#' @param radar_binary If `TRUE`, the radar response is a binary
#'   "migration night" (traffic above each radar's `radar_q` quantile) fitted
#'   with a binomial GAM, which puts both models on the same scale.
#' @param radar_q Quantile used when `radar_binary = TRUE`.
#' @return List with `bat` and `radar` models and the formulas used.
#' @export
fit_paired_gams <- function(bat_df, radar_df, covariates = "yday", k = 6,
                            radar_binary = TRUE, radar_q = 0.75) {
  sm <- paste0("s(", covariates, ", k = ", k, ")", collapse = " + ")
  bat_df <- tibble::as_tibble(sf::st_drop_geometry(bat_df))
  bat_df$individual <- factor(bat_df$individual)
  radar_df$radar <- factor(radar_df$radar)
  fb <- stats::as.formula(paste("migrated ~", sm,
                                if (nlevels(bat_df$individual) > 1) "+ s(individual, bs = 're')"))
  rr <- if (nlevels(radar_df$radar) > 1) "+ s(radar, bs = 're')" else ""
  if (radar_binary) {
    radar_df <- radar_df |>
      dplyr::group_by(.data$radar) |>
      dplyr::mutate(radar_migr = .data$mt_scaled > stats::quantile(.data$mt_scaled, radar_q, na.rm = TRUE)) |>
      dplyr::ungroup()
    fr <- stats::as.formula(paste("radar_migr ~", sm, rr))
    mr <- mgcv::bam(fr, data = radar_df, family = stats::binomial(), method = "fREML", discrete = TRUE)
  } else {
    radar_df$log_mt <- log1p(radar_df$mt_scaled)
    fr <- stats::as.formula(paste("log_mt ~", sm, rr))
    mr <- mgcv::bam(fr, data = radar_df, method = "fREML", discrete = TRUE)
  }
  mb <- mgcv::bam(fb, data = bat_df, family = stats::binomial(), method = "fREML", discrete = TRUE)
  list(bat = mb, radar = mr, formulas = list(bat = fb, radar = fr), covariates = covariates)
}

#' Compare partial effects between the tag and radar models
#'
#' Predicts each covariate's smooth over the range where both data sets
#' overlap, standardises each curve (mean 0, sd 1 on the link scale) so
#' shapes can be compared regardless of scale, and reports their correlation.
#'
#' @param fits Output of [fit_paired_gams()].
#' @param n Grid size.
#' @return List: `curves` (long tibble for plotting) and `summary`
#'   (correlation and sign agreement per covariate).
#' @export
compare_partial_effects <- function(fits, n = 100) {
  get_curve <- function(m, cov, grid) {
    nd <- as.data.frame(lapply(m$model, function(x) {
      if (is.numeric(x)) rep(stats::median(x, na.rm = TRUE), length(grid))
      else rep(x[1], length(grid))
    }))
    nd[[cov]] <- grid
    tm <- stats::predict(m, nd, type = "terms", se.fit = TRUE,
                         exclude = grep("s\\((individual|radar)\\)", names(m$sp), value = TRUE))
    lab <- paste0("s(", cov, ")")
    list(fit = tm$fit[, lab], se = tm$se.fit[, lab])
  }
  curves <- list(); summ <- list()
  for (cov in fits$covariates) {
    rb <- range(fits$bat$model[[cov]], na.rm = TRUE)
    rr <- range(fits$radar$model[[cov]], na.rm = TRUE)
    lo <- max(rb[1], rr[1]); hi <- min(rb[2], rr[2])
    if (!(hi > lo)) next
    grid <- seq(lo, hi, length.out = n)
    b <- get_curve(fits$bat, cov, grid); r <- get_curve(fits$radar, cov, grid)
    z <- function(x) (x - mean(x)) / stats::sd(x)
    curves[[cov]] <- dplyr::bind_rows(
      tibble::tibble(covariate = cov, x = grid, model = "tags", fit = b$fit, se = b$se, z = z(b$fit)),
      tibble::tibble(covariate = cov, x = grid, model = "radar", fit = r$fit, se = r$se, z = z(r$fit)))
    summ[[cov]] <- tibble::tibble(
      covariate = cov, lo = lo, hi = hi,
      r_curves = stats::cor(b$fit, r$fit),
      sign_agreement = mean(sign(diff(b$fit)) == sign(diff(r$fit))),
      edf_tags = sum(fits$bat$edf[grep(paste0("s\\(", cov, "\\)"), names(fits$bat$edf))]),
      edf_radar = sum(fits$radar$edf[grep(paste0("s\\(", cov, "\\)"), names(fits$radar$edf))]))
  }
  list(curves = dplyr::bind_rows(curves), summary = dplyr::bind_rows(summ))
}

#' Does a tag-envelope bat filter improve agreement?
#'
#' Recomputes radar traffic per segment using only cells inside a tag flight
#' envelope, and compares the AUC for tag migration with that of unfiltered
#' traffic.
#'
#' @param cells Output of [vpts_cells()].
#' @param env One-row envelope from [tag_envelope()].
#' @param paired Output of [pair_segments()] (restricted to the same group).
#' @return Tibble with AUC for all-traffic and envelope-traffic.
#' @export
envelope_gain <- function(cells, env, paired) {
  ce <- cells[in_envelope(cells, env), ]
  win <- paired[, c("segment_id", "radar", "t_start", "t_end", "migrated", "log_mtr")]
  ws <- if (nrow(ce)) window_sums(profile_table(ce), win[, c("segment_id", "radar", "t_start", "t_end")])
        else tibble::tibble()
  mt_env <- if (nrow(ws)) ws$mtr_mean[match(paste(win$segment_id, win$radar),
                                              paste(ws$segment_id, ws$radar))] else rep(NA, nrow(win))
  mt_env[is.na(mt_env)] <- 0
  tibble::tibble(
    n = nrow(win), n_migr = sum(win$migrated),
    auc_all = auc(win$log_mtr, win$migrated),
    auc_envelope = auc(log1p(mt_env), win$migrated),
    share_traffic_in_envelope = stats::median(mt_env / expm1(win$log_mtr), na.rm = TRUE))
}
