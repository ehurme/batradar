test_that("migration summaries count nights, distances and directions", {
  s <- bat_nights(fake_track())
  s$season <- "Fall"
  ms <- migration_summary(s, by = "season")
  expect_equal(ms$nights, 4)
  expect_equal(ms$migr_nights, 2)
  expect_equal(ms$n_ind_migr, 1)
  expect_equal(ms$dir_mean, 0, tolerance = 1)
  ph <- migration_phenology(s, by = "season", min_nights = 1)
  expect_true(all(ph$lo <= ph$p_migr & ph$p_migr <= ph$hi))
  im <- individual_migration(s)
  a <- im[im$individual == "a", ]
  expect_equal(a$migr_nights, 2)
  expect_equal(a$straightness, 1, tolerance = 0.01)
  expect_equal(a$km_net, haversine_km(9, 47, 9, 48.6), tolerance = 1e-6)
  expect_equal(im$km_total[im$individual == "b"], 0)
})

test_that("phenology GAM recovers a migration peak", {
  set.seed(9)
  n <- 3000
  d <- tibble::tibble(
    individual = sample(paste0("i", 1:40), n, TRUE),
    species = "N. noctula", season = "Fall",
    night_date = as.Date("2024-01-01") + sample(220:300, n, TRUE) - 1)
  d$yday <- as.integer(format(d$night_date, "%j"))
  d$migrated <- stats::rbinom(n, 1, stats::plogis(-3 + 3 * exp(-((d$yday - 260) / 12)^2))) == 1
  fit <- fit_phenology_gam(d)
  pr <- predict_phenology(fit)
  expect_equal(pr$yday[which.max(pr$p)], 260, tolerance = 6)
})
