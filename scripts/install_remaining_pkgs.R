## install_remaining_pkgs.R  – install what Conda can’t
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus  = parallel::detectCores(),
  # ← prevent CRAN/Bioc from upgrading Conda-installed binaries
  install.packages.compile.from.source = "never"
)

need  <- function(pkg) !requireNamespace(pkg, quietly = TRUE)

if (need("BiocManager"))
    install.packages("BiocManager", dependencies = FALSE)

## 1) consensusOV  (skip upgrades of existing deps)
if (need("consensusOV"))
    BiocManager::install("consensusOV",
                         update       = FALSE,     # ← KEY
                         ask          = FALSE,
                         dependencies = TRUE)

## 2) NMF  (skip upgrades too)
if (need("NMF"))
    install.packages("NMF", dependencies = TRUE, INSTALL_opts = "--no-test-load")

## 3) InstaPrism from GitHub
if (need("InstaPrism")) {
    if (need("remotes"))
        install.packages("remotes", dependencies = FALSE)
    remotes::install_github("humengying0907/InstaPrism",
                            upgrade = "never")      # ← don’t touch RcppArmadillo
}

message("✓ R-side package installation finished")
