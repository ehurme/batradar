# Which part of the radar signal is bats? ------------------------------------
# a) Tag flight envelope: restrict radar traffic to the heights and times
#    tagged bats fly; does agreement with tag migration improve?
# b) Cell model: presence (cells holding a tagged bat) vs background cells
#    from the same radar-nights; does adding radar features beat a
#    height + timing null model?
source("analysis/00_config.R")

fixes    <- readRDS(out("flight_fixes.rds"))
radars   <- readRDS(out("radars.rds"))
paired   <- readRDS(out("paired.rds"))
fix_vpts <- readRDS(out("fix_vpts.rds"))

vpts  <- vpts_standardise(read_vpts_cache(cfg$cache_dir))
cells <- vpts_cells(vpts, sun_elev = cfg$sun_elev, min_agl = cfg$min_agl_radar,
                    radar_lonlat = radars)
# Optional: ERA5 winds on radar bins -> airspeed
# cells <- vpts_cells(vpts, wind = era5_bins, ...)

# Where in the radar column do tagged bats fly?
fv <- fix_vpts |> filter(in_profile)
write.csv(fv |> group_by(species) |>
            summarise(n = n(), dens_rank_median = median(dens_rank_pct, na.rm = TRUE),
                      cum_below_median = median(cum_below, na.rm = TRUE),
                      bat_minus_radar_height = median(altitude_asl - height_mean_radar, na.rm = TRUE),
                      heading_diff_median = median(heading_diff, na.rm = TRUE)),
          out("fixes_in_radar_column.csv"), row.names = FALSE)

p <- ggplot(fv, aes(cum_below)) +
  geom_histogram(breaks = seq(0, 1, 0.05)) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Share of radar column density below the tagged bat",
       y = "Fixes", title = "0 = bat below all radar targets, 1 = above all")
ggsave(out("fig_fix_position_in_column.png"), p, width = 9, height = 5, dpi = 200)

# a) Envelope
env <- tag_envelope(fixes, by = "species")
write.csv(env, out("tag_flight_envelope.csv"), row.names = FALSE)
gain <- purrr::map_dfr(env$species, function(sp) {
  p <- paired[paired$species == sp, ]
  if (sum(p$migrated) < 5) return(NULL)
  cbind(species = sp, envelope_gain(cells, env[env$species == sp, ], p))
})
write.csv(gain, out("envelope_gain.csv"), row.names = FALSE)
print(gain)

# b) Cell model
lab <- label_cells(cells, fix_vpts, n_bg = 20)
bf  <- fit_bat_filter(lab)
write.csv(bf$cv, out("bat_filter_cv.csv"), row.names = FALSE)
saveRDS(bf$model, out("bat_filter_gam.rds"))
cat(sprintf("Bat filter AUC: full %.3f, null %.3f, gain %.3f\n",
            bf$auc_full, bf$auc_null, bf$auc_gain))
