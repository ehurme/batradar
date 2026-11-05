test_that("paired GAMs recover a shared seasonal peak", {
  set.seed(3)
  yd <- sample(200:300, 600, TRUE)
  p <- plogis(-1 + 2 * dnorm(yd, 250, 15) / dnorm(0, 0, 15))
  bat <- tibble::tibble(migrated = rbinom(600, 1, p) == 1, yday = yd,
                        individual = sample(letters[1:8], 600, TRUE))
  rad <- tibble::tibble(radar = sample(c("r1", "r2", "r3"), 400, TRUE),
                        yday = sample(200:300, 400, TRUE))
  rad$mt_scaled <- exp(3 + 2 * dnorm(rad$yday, 250, 15) / dnorm(0, 0, 15) + rnorm(400, 0, .3))
  fits <- fit_paired_gams(bat, rad, covariates = "yday")
  cmp <- compare_partial_effects(fits)
  expect_equal(cmp$summary$covariate, "yday")
  expect_gt(cmp$summary$r_curves, 0.8)
  expect_true(all(c("tags", "radar") %in% cmp$curves$model))
})

test_that("bat filter model separates presence cells by height", {
  set.seed(5)
  cells <- tibble::tibble(
    radar = "r1", night_date = as.Date("2024-09-10") + rep(0:9, each = 300),
    datetime = as.POSIXct("2024-09-10 20:00", tz = "UTC") + rep(0:9, each = 300) * 86400 + rep(1:300, 10),
    height_asl = sample(seq(200, 3000, 200), 3000, TRUE), height_agl_radar = NA,
    hrs_since_sunset = runif(3000, 0, 10), yday = 254, ff = runif(3000, 5, 15), dens = rexp(3000, .1))
  cells$height_agl_radar <- cells$height_asl
  cells$present <- rbinom(3000, 1, plogis(-4 + 3 * (cells$height_asl < 800))) == 1
  fit <- fit_bat_filter(cells, formula = present ~ s(height_agl_radar, k = 6) + s(hrs_since_sunset, k = 4),
                        formula_null = present ~ s(hrs_since_sunset, k = 4), k_folds = 5)
  expect_gt(fit$auc_full, 0.7)
  expect_gt(fit$auc_gain, 0.1)  # height is the informative term here
})
