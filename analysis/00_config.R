# Shared settings for the analysis scripts. Run scripts from the repository
# root, e.g. source("analysis/01_migration_trends.R").

# Load the functions in R/ (works without installing the package).
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(batradar)
}
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})

cfg <- list(
  # Tracking data: an .robj holding `bats_daily` (one fix per day) and
  # `bats_loc` (all fixes) as move2 objects. Override with an env variable.
  tracks_file = Sys.getenv(
    "BATRADAR_TRACKS",
    "C:/Users/User/Dropbox/MPI/Noctule/Data/rdata/move_icarus_bats_noonroost_daily_20261102.robj"),
  out_dir     = "output",
  cache_dir   = "cache/vpts",

  # Tag thresholds
  min_km      = 20,     # migration night: displacement >= this ...
  err_mult    = 2,      # ... and >= 2 x combined Sigfox location error
  max_gap_h   = 36,     # consecutive fixes further apart are not one night
  sun_elev    = -6,     # darkness threshold (deg)
  fix_min_agl = 50,     # drop near-ground "flight" fixes (uncalibrated pressure)

  # Radar
  vpts_source   = "baltrad",
  seg_radius_km = 50,   # segment passes this close to a radar
  fix_radius_km = 35,   # vol2bird outer range
  min_agl_radar = 200   # clutter floor above antenna (m)
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(...) file.path(cfg$out_dir, ...)

load_tracks <- function(file = cfg$tracks_file) {
  e <- new.env()
  load(file, envir = e)
  list(daily = e$bats_daily, loc = e$bats_loc)
}

theme_set(theme_bw(base_size = 11))
