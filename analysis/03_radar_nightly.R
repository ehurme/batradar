# Summarise radar profiles for nights and tag windows --------------------------
source("analysis/00_config.R")

radars   <- readRDS(out("radars.rds"))
segments <- readRDS(out("segments.rds"))
fixes    <- readRDS(out("flight_fixes.rds"))
m        <- readRDS(out("radar_matches.rds"))

vpts <- vpts_standardise(read_vpts_cache(cfg$cache_dir))

radar_nights <- vpts_nightly(vpts, cfg$sun_elev, cfg$min_agl_radar, radar_lonlat = radars)
seg_radar    <- vpts_for_segments(vpts, segments, m$seg, cfg$sun_elev,
                                  cfg$min_agl_radar, radar_lonlat = radars)
fix_vpts     <- match_fixes_to_vpts(fixes, m$fix, vpts, min_agl = cfg$min_agl_radar)

saveRDS(radar_nights, out("radar_nights.rds"))
saveRDS(seg_radar,    out("segment_radar.rds"))
saveRDS(fix_vpts,     out("fix_vpts.rds"))

# Radar migration phenology at the matched radars
p <- radar_nights |>
  filter(coverage >= 0.5) |>
  mutate(year = format(night_date, "%Y")) |>
  ggplot(aes(yday, log1p(mtr_mean))) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10)) +
  facet_wrap(~ year) +
  labs(x = "Day of year", y = "log(1 + mean MTR) [animals/km/h]",
       title = "Radar nocturnal migration at radars near tagged bats")
ggsave(out("fig_radar_phenology.png"), p, width = 9, height = 6, dpi = 200)
