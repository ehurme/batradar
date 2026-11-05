# Run every step in order from the repository root.
for (f in c("01_migration_trends.R", "02_radar_download.R", "03_radar_nightly.R",
            "04_compare_tags_radar.R", "05_bat_filter.R")) {
  message("== ", f)
  source(file.path("analysis", f), local = new.env())
}
