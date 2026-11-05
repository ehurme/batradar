# Match tags to European weather radars and download VPTS ---------------------
# Needs internet access to the OPERA radar list and the Aloft VPTS archive
# (getRad, from https://aloftdata.r-universe.dev). Downloads are cached per
# radar-month in cfg$cache_dir, so this can be re-run safely.
source("analysis/00_config.R")

tr <- load_tracks()
segments <- if (file.exists(out("segments.rds"))) readRDS(out("segments.rds")) else
  bat_nights(tr$daily, min_km = cfg$min_km, err_mult = cfg$err_mult, max_gap_h = cfg$max_gap_h)
fixes <- bat_flight_fixes(tr$loc, sun_elev = cfg$sun_elev, min_agl = cfg$fix_min_agl)
saveRDS(fixes, out("flight_fixes.rds"))

radars <- radar_sites("opera")
saveRDS(radars, out("radars.rds"))

seg_matches <- match_segments_to_radars(segments, radars, cfg$seg_radius_km)
fix_matches <- match_fixes_to_radars(fixes, radars, cfg$fix_radius_km)
saveRDS(list(seg = seg_matches, fix = fix_matches), out("radar_matches.rds"))

coverage <- getRad::get_vpts_coverage(source = cfg$vpts_source)
requests <- radar_requests(segments, seg_matches, fixes, fix_matches)
message(nrow(requests), " radar-days requested across ",
        dplyr::n_distinct(requests$radar), " radars")

vpts_raw <- fetch_vpts(requests, cfg$vpts_source, cfg$cache_dir, coverage = coverage)
fails <- attr(vpts_raw, "failures")
if (nrow(fails)) write.csv(fails, out("vpts_failures.csv"), row.names = FALSE)

# Overview map: radars used and migration segments
used <- radars[radars$radar %in% unique(seg_matches$radar), ]
p <- ggplot() +
  geom_sf(data = segments[segments$migrated, ], colour = "grey50", linewidth = 0.2) +
  geom_sf(data = used, colour = "red", size = 1.5) +
  labs(title = "Radars within reach of tagged-bat migration segments")
ggsave(out("fig_map_radars.png"), p, width = 8, height = 8, dpi = 200)
