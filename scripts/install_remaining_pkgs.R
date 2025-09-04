# ------------------------------------------------------------------------------
# Description:
#   Completes R package setup *inside* the conda environment `env_deconv_R`
#   by installing Bioconductor and GitHub-only packages that are not available
#   (or pinned) via conda-forge/bioconda. Prefers prebuilt binaries and sets
#   parallel install options.
#
# What it does:
#   - Ensures BiocManager and remotes are available
#   - Installs Bioconductor pkgs: GSVA, genefu, impute, SpatialExperiment, amap
#   - Installs GitHub pkgs: bhklab/consensusOV, humengying0907/InstaPrism
#   - Prints a success message on completion
#
# Inputs:
#   - None (runs online installs; requires network access)
#
# Outputs:
#   - Packages installed into the active R library of `env_deconv_R`
#
# Requirements:
#   - Activate conda env first: `conda activate env_deconv_R`
#   - Internet access; R can reach CRAN, Bioconductor, and GitHub
# ------------------------------------------------------------------------------

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
