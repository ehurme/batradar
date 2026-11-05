# batradar pipeline -------------------------------------------------------
# Copy into a project folder with batradar::use_pipeline(), edit the CONFIG
# block, then run targets::tar_make(). Radar downloads are cached in
# cache/vpts so re-runs only fetch new radar-months.
library(targets)
tar_option_set(packages = c("batradar", "sf", "dplyr"), format = "rds")

# CONFIG ------------------------------------------------------------------
cfg <- list(
  tracks_file   = "C:/Users/User/Dropbox/MPI/Noctule/Data/rdata/move_icarus_bats_noonroost_daily_20261102.robj",
  daily_obj     = "bats_daily",   # one fix per day -> nightly segments
  loc_obj       = "bats_loc",     # all fixes -> in-flight fixes
  vpts_source   = "baltrad",      # getRad source
  seg_radius_km = 50,             # segment must pass this close to a radar
  fix_radius_km = 35,             # flight fix within vol2bird range
  min_km        = 20,             # migration threshold (also 2x location error)
  fix_min_agl   = 50,             # drop near-ground "flight" fixes
  sun_elev      = -6,
  min_agl_radar = 200,            # radar clutter floor above antenna
  cache_dir     = "cache/vpts"
)

load_obj <- function(file, name) { e <- new.env(); load(file, envir = e); e[[name]] }

list(
  tar_target(tracks_file, cfg$tracks_file, format = "file"),
  tar_target(daily, load_obj(tracks_file, cfg$daily_obj)),
  tar_target(loc,   load_obj(tracks_file, cfg$loc_obj)),

  # Tags
  tar_target(segments, bat_nights(daily, min_km = cfg$min_km)),
  tar_target(fixes, bat_flight_fixes(loc, sun_elev = cfg$sun_elev, min_agl = cfg$fix_min_agl)),
  tar_target(envelope, tag_envelope(fixes, by = "species")),

  # Radars and matching
  tar_target(radars, radar_sites("opera")),
  tar_target(seg_matches, match_segments_to_radars(segments, radars, cfg$seg_radius_km)),
  tar_target(fix_matches, match_fixes_to_radars(fixes, radars, cfg$fix_radius_km)),
  tar_target(coverage, getRad::get_vpts_coverage(source = cfg$vpts_source)),
  tar_target(requests, radar_requests(segments, seg_matches, fixes, fix_matches)),

  # Radar data (cached per radar-month)
  tar_target(vpts_raw, fetch_vpts(requests, cfg$vpts_source, cfg$cache_dir, coverage = coverage)),
  tar_target(vpts, vpts_standardise(vpts_raw)),
  tar_target(radar_nights, vpts_nightly(vpts, cfg$sun_elev, cfg$min_agl_radar, radar_lonlat = radars)),
  tar_target(seg_radar, vpts_for_segments(vpts, segments, seg_matches, cfg$sun_elev,
                                          cfg$min_agl_radar, radar_lonlat = radars)),
  tar_target(fix_vpts, match_fixes_to_vpts(fixes, fix_matches, vpts, min_agl = cfg$min_agl_radar)),

  # 1. Night-level agreement: does radar see migration when our bats move?
  tar_target(paired, pair_segments(segments, seg_matches, seg_radar)),
  tar_target(agree_all, night_agreement(paired)),
  tar_target(agree_species, night_agreement(paired, by = c("species", "season"))),
  tar_target(detection_gam, fit_detection_gam(paired)),

  # 2. Same models, both data sets (phenology now; add ERA5 covariates later)
  tar_target(paired_gams, {
    rn <- radar_nights[radar_nights$radar %in% unique(paired$radar) &
                         radar_nights$yday >= min(paired$yday) &
                         radar_nights$yday <= max(paired$yday), ]
    fit_paired_gams(paired, rn, covariates = "yday")
  }),
  tar_target(effects, compare_partial_effects(paired_gams)),

  # 3. Bat filter: tag envelope gain + learned cell model
  tar_target(cells, vpts_cells(vpts, sun_elev = cfg$sun_elev, min_agl = cfg$min_agl_radar,
                               radar_lonlat = radars)),
  tar_target(envelope_gain_tbl, {
    purrr::map_dfr(envelope$species, function(sp) {
      p <- paired[paired$species == sp, ]
      if (sum(p$migrated) < 5) return(NULL)
      cbind(species = sp, envelope_gain(cells, envelope[envelope$species == sp, ], p))
    })
  }),
  tar_target(labelled, label_cells(cells, fix_vpts)),
  tar_target(bat_filter, fit_bat_filter(labelled)),

  # Report tables
  tar_target(report, list(agree_all = agree_all, agree_species = agree_species,
                          effects = effects$summary, envelope = envelope,
                          envelope_gain = envelope_gain_tbl,
                          bat_filter_cv = bat_filter$cv,
                          fix_height = fix_vpts))
)
