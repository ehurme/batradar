test_that("sun altitude matches solstice noon geometry", {
  a <- sun_altitude(as.POSIXct("2024-06-21 12:00", tz = "UTC"), 0, 51.5)
  expect_equal(a, 90 - 51.5 + 23.44, tolerance = 0.3)
  expect_lt(sun_altitude(as.POSIXct("2024-06-21 00:00", tz = "UTC"), 0, 51.5), -10)
})

test_that("sunset time is sensible for Konstanz in September", {
  s <- sunset_time(as.Date("2024-09-10"), 9.17, 47.66)
  h <- as.numeric(format(s, "%H")) + as.numeric(format(s, "%M")) / 60
  expect_true(h > 17.2 && h < 18.2)  # ~19:40 CEST = ~17:40 UTC
})

test_that("bearing, distance and angle helpers", {
  expect_equal(bearing_deg(0, 0, 0, 1), 0, tolerance = 1e-6)
  expect_equal(bearing_deg(0, 0, 1, 0), 90, tolerance = 1e-6)
  expect_equal(haversine_km(0, 0, 0, 1), 111.2, tolerance = 0.01)
  expect_equal(angle_diff(350, 10), 20)
  expect_equal(auc(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(auc(c(1, 2, 3, 4), c(TRUE, TRUE, FALSE, FALSE)), 0)
})
