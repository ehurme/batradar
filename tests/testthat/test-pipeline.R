test_that("bat_nights flags migration beyond location error", {
  s <- bat_nights(fake_track())
  expect_s3_class(s, "sf")
  expect_equal(nrow(s), 4)
  a <- s[s$individual == "a", ]
  expect_true(all(a$migrated))
  expect_equal(a$bearing[1], 0, tolerance = 0.5)
  expect_false(any(s$migrated[s$individual == "b"]))
  expect_equal(s$night_date[1], as.Date("2024-09-10"))
})

test_that("segments match only nearby radars", {
  s <- bat_nights(fake_track()); r <- fake_radars()
  m <- match_segments_to_radars(s, r, max_km = 50)
  expect_true(all(m$radar == "xxtst"))
  expect_true(all(s$segment_id %in% m$segment_id))
  req <- radar_requests(s, m)
  expect_setequal(as.character(req$date), c("2024-09-10", "2024-09-11", "2024-09-12"))
})

test_that("nightly traffic integrates density x speed over darkness", {
  v <- vpts_standardise(fake_vpts())
  expect_equal(unique(v$bin_km), 0.2)
  expect_equal(unique(v$dt_h), 5 / 60)
  n <- vpts_nightly(v)
  expect_equal(nrow(n), 1)
  # 4 bins (200-800 m) x 10 ind/km3 x 10 m/s x 3.6 x 0.2 km = 288 ind/km/h
  expect_equal(n$mtr_mean, 288)
  expect_equal(n$mt, 288 * n$hours)
  expect_equal(n$dir_mean, 0, tolerance = 1e-6)
  expect_equal(n$height_mean, 600)  # bin centres 300..900 -> 600
  expect_true(n$hours > 10 && n$hours < 12)
})

test_that("antenna-referenced heights are shifted to sea level", {
  v <- fake_vpts(); v$height_reference <- "antenna"; v$radar_height <- 500
  s <- vpts_standardise(v)
  expect_equal(min(s$height_asl), 500)
  expect_equal(min(s$height_agl_radar), 0)
})

test_that("segment windows and pairing give AUC and direction agreement", {
  # Night 1 busy, night 2 empty; bat 'a' migrates both nights, 'b' never.
  busy <- function(h, t) ifelse(h < 1000 & t < as.POSIXct("2024-09-11 12:00", tz = "UTC"), 10, 0)
  v <- vpts_standardise(fake_vpts(t1 = as.POSIXct("2024-09-12 12:00", tz = "UTC"), dens_fun = busy))
  s <- bat_nights(fake_track()); r <- fake_radars()
  m <- match_segments_to_radars(s, r)
  sr <- vpts_for_segments(v, s, m)
  expect_equal(nrow(sr), 4)
  p <- pair_segments(s, m, sr)
  expect_equal(nrow(p), 4)
  na <- night_agreement(p)
  expect_equal(na$n_migr, 2)
  expect_true(is.finite(na$auc_mtr))
  expect_equal(p$dir_diff[p$migrated & p$log_mtr > 0], 0, tolerance = 1)
})

test_that("fixes are placed in the right radar bin", {
  v <- vpts_standardise(fake_vpts(dens_fun = function(h, t) ifelse(h == 600, 50, ifelse(h < 1000, 10, 0))))
  fx <- sf::st_as_sf(tibble::tibble(
    individual = "a", species = "N", tag_type = "Tinyfox",
    time = as.POSIXct("2024-09-10 22:02", tz = "UTC"),
    lon = 9.05, lat = 47.75, altitude_asl = 650, heading = 10, fix_id = 1L),
    coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  fm <- match_fixes_to_radars(fx, fake_radars())
  expect_equal(nrow(fm), 1)
  out <- match_fixes_to_vpts(fx, fm, v)
  expect_equal(out$bin_height_asl, 600)
  expect_equal(out$dens_at_bat, 50)
  expect_equal(out$dens_rank_pct, 1)
  expect_equal(out$heading_diff, 10)
  # bins kept above the 200 m clutter floor: 200, 400 (10 each), 600 (50), 800 (10)
  expect_equal(out$cum_below, 20 / 80, tolerance = 1e-8)
})
