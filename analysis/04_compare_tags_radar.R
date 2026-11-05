# Do radar and tags agree? ----------------------------------------------------
# 1. Detection: radar traffic on nights bats migrated past vs stayed.
# 2. Direction: bat bearing vs radar mean track direction.
# 3. Same GAM on both data sets: compare partial effects.
source("analysis/00_config.R")

segments     <- readRDS(out("segments.rds"))
m            <- readRDS(out("radar_matches.rds"))
seg_radar    <- readRDS(out("segment_radar.rds"))
radar_nights <- readRDS(out("radar_nights.rds"))

paired <- pair_segments(segments, m$seg, seg_radar)
saveRDS(paired, out("paired.rds"))

agree <- bind_rows(
  night_agreement(paired) |> mutate(species = "all", season = "all"),
  night_agreement(paired, by = c("species", "season")))
write.csv(agree, out("agreement_detection_direction.csv"), row.names = FALSE)
print(agree, width = Inf)

det <- fit_detection_gam(paired)
print(summary(det))
saveRDS(det, out("detection_gam.rds"))

# Radar traffic on migration vs non-migration nights
p1 <- ggplot(paired, aes(migrated, log_mtr)) +
  geom_violin(fill = "grey85") + geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_grid(season ~ species) +
  labs(x = "Tagged bat migrated past the radar", y = "log(1 + radar MTR)")
ggsave(out("fig_detection.png"), p1, width = 10, height = 6, dpi = 200)

# Direction agreement
p2 <- paired |> filter(migrated) |>
  ggplot(aes(bearing, dir_mean, colour = season)) +
  geom_abline(linetype = 2) + geom_point(alpha = 0.6) +
  coord_equal(xlim = c(0, 360), ylim = c(0, 360)) +
  facet_wrap(~ species) +
  labs(x = "Bat bearing (deg)", y = "Radar mean direction (deg)")
ggsave(out("fig_direction_agreement.png"), p2, width = 9, height = 6, dpi = 200)

# Same model, both data sets: phenology now. Join ERA5 covariates to both
# `paired` and `rn` (by radar and night_date) and add them to `covariates`.
rn <- radar_nights |>
  filter(radar %in% unique(paired$radar), coverage >= 0.5,
         yday >= min(paired$yday), yday <= max(paired$yday))
fits <- fit_paired_gams(paired, rn, covariates = "yday")
eff  <- compare_partial_effects(fits)
write.csv(eff$summary, out("paired_gam_effects.csv"), row.names = FALSE)
print(eff$summary)

p3 <- ggplot(eff$curves, aes(x, z, colour = model)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ covariate, scales = "free_x") +
  labs(x = NULL, y = "Standardised partial effect", colour = NULL)
ggsave(out("fig_paired_effects.png"), p3, width = 8, height = 4, dpi = 200)
