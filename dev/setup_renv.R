# One-time setup on your machine (needs internet to CRAN + aloftdata r-universe).
# Run from the package root.
options(repos = c(aloftdata = "https://aloftdata.r-universe.dev",
                  CRAN = "https://cloud.r-project.org"))
install.packages("renv")
renv::init(bare = TRUE)
renv::install(c("getRad", "bioRad", "move2", "targets", "tarchetypes", "ggplot2",
                "sf", "dplyr", "tidyr", "mgcv", "purrr", "tibble", "cli", "rlang",
                "testthat", "devtools", "roxygen2"))
renv::install(".")      # batradar itself
renv::snapshot()
