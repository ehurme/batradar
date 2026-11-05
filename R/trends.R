#' Summarise tag-based migration by group
#'
#' One row per group (default species x season x year): sample sizes,
#' how often bats migrated, how far, when, and in which direction.
#'
#' @param segments Output of [bat_nights()].
#' @param by Grouping columns. `year` is derived from `night_date`.
#' @return Tibble.
#' @export
migration_summary <- function(segments, by = c("species", "season", "year")) {
  s <- tibble::as_tibble(sf::st_drop_geometry(segments))
  s$year <- as.integer(format(s$night_date, "%Y"))
  s |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_ind = dplyr::n_distinct(.data$individual),
      n_ind_migr = dplyr::n_distinct(.data$individual[.data$migrated]),
      nights = dplyr::n(),
      migr_nights = sum(.data$migrated),
      p_migr = mean(.data$migrated),
      km_median = stats::median(.data$dist_km[.data$migrated]),
      km_q90 = if (any(.data$migrated)) stats::quantile(.data$dist_km[.data$migrated], 0.9) else NA_real_,
      km_max = if (any(.data$migrated)) max(.data$dist_km[.data$migrated]) else NA_real_,
      yday_first = if (any(.data$migrated)) min(.data$yday[.data$migrated]) else NA_integer_,
      yday_median = stats::median(.data$yday[.data$migrated]),
      yday_last = if (any(.data$migrated)) max(.data$yday[.data$migrated]) else NA_integer_,
      dir_mean = circ_mean(.data$bearing[.data$migrated])[["mean"]],
      dir_rho = circ_mean(.data$bearing[.data$migrated])[["rho"]],
      .groups = "drop")
}

#' Weekly migration phenology
#'
#' Share of tracked bat-nights with migration, per week of year and group.
#'
#' @param segments Output of [bat_nights()].
#' @param by Grouping columns.
#' @param min_nights Weeks with fewer tracked nights are dropped.
#' @return Tibble with `week`, `nights`, `migr_nights`, `p_migr`, and a
#'   Wilson 95% interval (`lo`, `hi`).
#' @export
migration_phenology <- function(segments, by = c("species", "season"), min_nights = 5) {
  s <- tibble::as_tibble(sf::st_drop_geometry(segments))
  s$week <- as.integer(format(s$night_date, "%V"))
  out <- s |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "week")))) |>
    dplyr::summarise(nights = dplyr::n(), migr_nights = sum(.data$migrated),
                     n_ind = dplyr::n_distinct(.data$individual), .groups = "drop") |>
    dplyr::filter(.data$nights >= min_nights)
  p <- out$migr_nights / out$nights; n <- out$nights; z <- 1.96
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  out$p_migr <- p; out$lo <- centre - half; out$hi <- centre + half
  out
}

#' Per-individual migration metrics
#'
#' For each individual and season: number of migration nights, total and
#' net displacement, straightness, mean heading, and travel speed over the
#' migratory period (first to last migration night).
#'
#' @param segments Output of [bat_nights()].
#' @return Tibble.
#' @export
individual_migration <- function(segments) {
  s <- tibble::as_tibble(sf::st_drop_geometry(segments))
  s$year <- as.integer(format(s$night_date, "%Y"))
  s |>
    dplyr::arrange(.data$individual, .data$t_start) |>
    dplyr::group_by(.data$individual, .data$species, .data$season, .data$year) |>
    dplyr::summarise(
      tag_type = dplyr::first(.data$tag_type),
      sex = dplyr::first(.data$sex),
      nights = dplyr::n(),
      migr_nights = sum(.data$migrated),
      first_migr = if (any(.data$migrated)) min(.data$night_date[.data$migrated]) else as.Date(NA),
      last_migr = if (any(.data$migrated)) max(.data$night_date[.data$migrated]) else as.Date(NA),
      km_total = sum(.data$dist_km[.data$migrated]),
      km_net = if (any(.data$migrated)) {
        i0 <- which(.data$migrated)[1]; i1 <- utils::tail(which(.data$migrated), 1)
        haversine_km(.data$lon_start[i0], .data$lat_start[i0], .data$lon_end[i1], .data$lat_end[i1])
      } else 0,
      dir_mean = circ_mean(.data$bearing[.data$migrated])[["mean"]],
      .groups = "drop") |>
    dplyr::mutate(
      straightness = ifelse(.data$km_total > 0, .data$km_net / .data$km_total, NA_real_),
      migr_days = as.numeric(.data$last_migr - .data$first_migr) + 1,
      km_per_day = ifelse(.data$migr_days > 0, .data$km_net / .data$migr_days, NA_real_))
}

#' Phenology GAM of nightly migration probability
#'
#' Binomial `bam` of migration on day of year, with a separate smooth per
#' species x season, and random intercepts for individual and year.
#'
#' @param segments Output of [bat_nights()].
#' @param k Basis dimension of the day-of-year smooths.
#' @param min_migr Groups with fewer migration nights are dropped.
#' @return A fitted `mgcv::bam` with the modelled data in `$model`.
#' @export
fit_phenology_gam <- function(segments, k = 8, min_migr = 10) {
  s <- tibble::as_tibble(sf::st_drop_geometry(segments))
  s$year <- factor(format(s$night_date, "%Y"))
  s$group <- interaction(s$species, s$season, drop = TRUE, sep = " / ")
  keep <- names(which(tapply(s$migrated, s$group, sum) >= min_migr))
  s <- s[s$group %in% keep, ]
  s$group <- droplevels(s$group)
  s$individual <- factor(s$individual)
  f <- if (nlevels(s$group) > 1) {
    migrated ~ group + s(yday, by = group, k = k) + s(individual, bs = "re")
  } else {
    migrated ~ s(yday, k = k) + s(individual, bs = "re")
  }
  if (nlevels(s$year) > 1) f <- stats::update(f, . ~ . + s(year, bs = "re"))
  fit <- mgcv::bam(f, data = s, family = stats::binomial(), method = "fREML", discrete = TRUE)
  fit$pheno_data <- s
  fit
}

#' Predicted phenology curves from [fit_phenology_gam()]
#'
#' Population-level curves (random effects excluded) over the observed
#' day-of-year range of each group.
#'
#' @param fit Output of [fit_phenology_gam()].
#' @param n Points per curve.
#' @return Tibble with `group`, `yday`, `p`, `lo`, `hi` (95% interval).
#' @export
predict_phenology <- function(fit, n = 100) {
  d <- fit$pheno_data
  ex <- grep("s\\((individual|year)\\)", sapply(fit$smooth, `[[`, "label"), value = TRUE)
  purrr::map_dfr(levels(d$group), function(g) {
    r <- range(d$yday[d$group == g])
    nd <- data.frame(group = factor(g, levels = levels(d$group)),
                     yday = seq(r[1], r[2], length.out = n),
                     individual = d$individual[1])
    if ("year" %in% names(d)) nd$year <- d$year[1]
    pr <- stats::predict(fit, nd, type = "link", se.fit = TRUE, exclude = ex)
    tibble::tibble(group = g, yday = nd$yday, p = stats::plogis(pr$fit),
                   lo = stats::plogis(pr$fit - 1.96 * pr$se.fit),
                   hi = stats::plogis(pr$fit + 1.96 * pr$se.fit))
  })
}
