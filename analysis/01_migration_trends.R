# Tag-based migration trends ----------------------------------------------
# Nightly movement segments from daily Sigfox fixes: how often, when, how far
# and in which direction each species migrates, by season and year.
source("analysis/00_config.R")

tr <- load_tracks()
segments <- bat_nights(tr$daily, min_km = cfg$min_km, err_mult = cfg$err_mult,
                       max_gap_h = cfg$max_gap_h)
saveRDS(segments, out("segments.rds"))

# Tables ---------------------------------------------------------------------
summ_sp  <- migration_summary(segments, by = c("species", "season"))
summ_yr  <- migration_summary(segments, by = c("species", "season", "year"))
ind      <- individual_migration(segments)
pheno    <- migration_phenology(segments, by = c("species", "season"))
write.csv(summ_sp, out("migration_summary_species_season.csv"), row.names = FALSE)
write.csv(summ_yr, out("migration_summary_species_season_year.csv"), row.names = FALSE)
write.csv(ind,     out("individual_migration.csv"), row.names = FALSE)
print(summ_sp, width = Inf)

# Phenology model --------------------------------------------------------------
pheno_fit <- fit_phenology_gam(segments)
saveRDS(pheno_fit, out("phenology_gam.rds"))
print(summary(pheno_fit)$s.table)
pheno_pred <- predict_phenology(pheno_fit)

doy_breaks <- c(60, 91, 121, 152, 182, 213, 244, 274, 305)
doy_labels <- c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")

p_pheno <- ggplot(pheno_pred, aes(yday, p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ group, scales = "free_x") +
  scale_x_continuous(breaks = doy_breaks, labels = doy_labels) +
  labs(x = NULL, y = "P(migration night)",
       title = "Nightly migration probability (population-level GAM)")
ggsave(out("fig_phenology_gam.png"), p_pheno, width = 9, height = 6, dpi = 200)

# Distances ------------------------------------------------------------------
mig <- st_drop_geometry(segments) |> filter(migrated)
p_dist <- ggplot(mig, aes(dist_km, fill = season)) +
  geom_histogram(binwidth = 20, position = "identity", alpha = 0.6) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Displacement on migration nights (km)", y = "Nights")
ggsave(out("fig_distance.png"), p_dist, width = 9, height = 6, dpi = 200)

# Directions -----------------------------------------------------------------
p_dir <- ggplot(mig, aes(bearing, fill = season)) +
  geom_histogram(breaks = seq(0, 360, 15), alpha = 0.7, position = "identity") +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270),
                     labels = c("N", "E", "S", "W")) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = NULL, y = "Migration nights")
ggsave(out("fig_direction.png"), p_dir, width = 9, height = 7, dpi = 200)

# Map of migration segments ----------------------------------------------------
p_map <- ggplot(segments[segments$migrated, ]) +
  geom_sf(aes(colour = season), linewidth = 0.3, alpha = 0.5) +
  facet_wrap(~ species) +
  labs(colour = NULL)
ggsave(out("fig_map_segments.png"), p_map, width = 10, height = 8, dpi = 200)

# Year-to-year timing ----------------------------------------------------------
p_year <- summ_yr |>
  filter(migr_nights >= 5) |>
  ggplot(aes(factor(year), yday_median, colour = season)) +
  geom_pointrange(aes(ymin = yday_first, ymax = yday_last),
                  position = position_dodge(width = 0.4)) +
  facet_wrap(~ species) +
  scale_y_continuous(breaks = doy_breaks, labels = doy_labels) +
  labs(x = NULL, y = "Migration nights (first, median, last)")
ggsave(out("fig_timing_by_year.png"), p_year, width = 9, height = 6, dpi = 200)
