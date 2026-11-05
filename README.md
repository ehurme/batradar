# batradar

Code for analysing bat migration from Sigfox tracking data (ICARUS Bats) and
comparing it with European weather-radar vertical profiles (VPTS from the
[Aloft](https://aloftdata.eu) archive, via
[`getRad`](https://aloftdata.r-universe.dev/getRad)).

The repository is laid out as an R package so the functions in `R/` are
documented and tested, but you don't need to install it: the scripts in
`analysis/` load everything with `pkgload::load_all()`.

```
R/            functions (tracks, radar matching, VPTS summaries, models, trends)
analysis/     numbered scripts: run from the repository root
inst/pipeline targets version of the same workflow
tests/        unit tests on synthetic data (no real data in the repo)
dev/          renv setup
```

No tracking data, radar downloads or outputs are committed: `output/` and
`cache/` are git-ignored.

## Workflow

| Script | What it does | Needs internet |
|---|---|---|
| `analysis/01_migration_trends.R` | Nightly movement segments from daily fixes; per species/season/year summaries; individual migration metrics; phenology GAM; figures (phenology, distances, directions, map, timing by year) | no |
| `analysis/02_radar_download.R` | OPERA radar sites, tag-radar matching, VPTS download with a per radar-month cache | yes |
| `analysis/03_radar_nightly.R` | Radar migration per radar-night and per tag segment window; places in-flight fixes in radar height bins | no |
| `analysis/04_compare_tags_radar.R` | Detection (radar traffic on nights bats migrated vs stayed), direction agreement, same GAM fitted to both data sets | no |
| `analysis/05_bat_filter.R` | Where tagged bats sit in the radar column; tag flight envelope as a bat filter; presence/background cell model vs a height + timing null | no |

```r
# once, on a machine with internet
source("dev/setup_renv.R")

# point to the tracking file (an .robj with move2 objects bats_daily and bats_loc)
Sys.setenv(BATRADAR_TRACKS = "C:/Users/User/Dropbox/MPI/Noctule/Data/rdata/move_icarus_bats_noonroost_daily_20261102.robj")
source("analysis/run_all.R")
```

Thresholds (migration distance, radar radius, darkness, clutter floor) are in
`analysis/00_config.R`.

## Why the radar comparison is built this way

vol2bird VPTS separate biology from rain but **not bats from birds**. The tags
are used in two ways: as an independent record of when and where bats
migrated, and as ground truth for which part of the radar signal is likely to
be bats.

- **Night dates** are anchored at local solar noon for tags and radar alike.
- **Migration night** = displacement ≥ max(20 km, 2 × combined Sigfox location
  error).
- **Segment radar traffic** is the mean traffic rate (animals/km/h) over the
  dark part of the window between two fixes. A total would grow with window
  length on its own.
- **Heights**: tag pressure altitude (m a.s.l.) is matched to radar bins in
  m a.s.l.; antenna-referenced profiles are shifted by `radar_height`; radar
  bins < 200 m above the antenna are dropped.
- **Bat filter evaluation** always reports the gain over a height + timing null
  model. A high AUC alone only means the model learned where bats fly.
- **Weather**: join ERA5 variables to the segment table and to `radar_nights`
  and pass their names to `fit_paired_gams(covariates = ...)`; pass ERA5 winds
  on radar bins to `vpts_cells(wind = ...)` for airspeed.

## Known data caveat

Some "in flight" fixes have pressure altitudes at or below ground level
(notably Nanofox on *N. noctula*), which points to uncalibrated pressure
altitude or the flight flag catching non-flying fixes. `fix_min_agl` (default
50 m) drops these; calibrating against ERA5 surface pressure would be better.

## Tests

```r
devtools::test()
```
