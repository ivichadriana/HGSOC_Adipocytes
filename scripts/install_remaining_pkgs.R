## This runs *inside* env_deconv_R after conda‑env creation
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = parallel::detectCores(),
  install.packages.compile.from.source = "never"  # use binaries where possible
)

need <- function(p) !requireNamespace(p, quietly = TRUE)

if (need("BiocManager"))
  install.packages("BiocManager", dependencies = FALSE)

if (need("remotes"))
  install.packages("remotes", dependencies = FALSE)

## ── Bioconductor pkgs not on conda-forge / bioconda ───────────────────────
BiocManager::install(
  c("GSVA", "genefu", "impute", "SpatialExperiment", "amap"),
  ask = FALSE, update = FALSE
)

## ── GitHub‑only packages ──────────────────────────────────────────────────
remotes::install_github(
  "bhklab/consensusOV",
  build_vignettes = FALSE,
  upgrade         = "never"
)

remotes::install_github(
  "humengying0907/InstaPrism",
  build_vignettes = FALSE,
  upgrade         = "never"
)

cat("✓ All remaining R packages installed\n")
