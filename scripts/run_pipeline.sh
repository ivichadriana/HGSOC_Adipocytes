#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Pipeline – master driver script
# ------------------------------------------------------------------------------
# This script orchestrates the full analysis workflow, including:
#   1. Creating / updating the required conda environments.
#   2. Unzipping raw input data (Python).
#   3. Running the core deconvolution pipeline (deconvolution in R).
#   4. Executing a series of Jupyter notebooks with Papermill to generate
#      downstream analysis and visualizations.
#
# Usage (example):
#   sbatch run_pipeline.sh   # Run via SLURM (recommended)
#   ./run_pipeline.sh        # Run locally (make sure env paths exist)
#
# Requirements:
#   * Two environments defined by scripts/0_create_environment.py:
#       - env_deconv      : Python‑focused tools (papermill, etc.)
#       - env_deconv_R    : R runtime with required Bioconductor packages.
# ------------------------------------------------------------------------------

set -Eeuo pipefail      # -E  : inherit traps, -e : exit on error, -u : unset vars, -o pipefail : propagate failures
IFS=$'\n\t'             # safer default IFS
PROJECT_ROOT="$(realpath "$(dirname "$0")/..")"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH+:$PYTHONPATH}"

# ---------------------------- 0. Environment setup ----------------------------
# Build or update conda environments.
# bash scripts/0_create_environments.sh
source activate base 
# ------------------------------ 1. Python phase -------------------------------
echo "→ Activating Python environment: env_deconv"
conda activate env_deconv
sleep 5  # give conda a moment to fully initialize

echo "→ Unzipping input data"
python scripts/1_unzip_input_data.py

echo "→ Deactivating Python environment"
conda deactivate
sleep 5

# ------------------------------- 2. R phase -----------------------------------
echo "→ Activating R environment: env_deconv_R"
set +u                    # disable “unbound variable” checking
conda activate env_deconv_R
set -u                    # re‑enable after activation
sleep 5

# Run the R‑based deconvolution workflow (each step is idempotent).
echo "→ Running R deconvolution pipeline"
# Rscript --vanilla scripts/2_process_data_and_subtypes.R
Rscript --vanilla scripts/3_prepare_deconvolution.R
Rscript --vanilla scripts/4_run_deconvolution.R

echo "→ Deactivating R environment"
set +u               
conda deactivate
set -u              
sleep 5

# --------------------------- 3. Notebook reporting ----------------------------
Reactivate Python environment for Papermill‑based reporting.

echo "→ Reactivating Python environment for reporting"
set +u 
conda activate env_deconv
set -u   

# Resolve the project base directory. Works for SLURM and local execution.
BASE_DIR=$(realpath "${SLURM_SUBMIT_DIR:-$(pwd)}/..")
NB_DIR="${BASE_DIR}/notebooks"

# Array of notebook basenames (without .ipynb)
NOTEBOOKS=(
  analysis_proportions_visualizations
  analysis_proportions_vs_bmi_and_age
  analysis_proportions_vs_race
  analysis_proportions_vs_residual
  analysis_proportions_vs_stage
  analysis_proportions_vs_subtype
  analysis_proportions_vs_subtype_and_race
  analysis_proportions_vs_survival
  analysis_proportions_vs_survival_and_race
  analysis_proportions_vs_tissue
)

# Iterate over notebooks and execute each with Papermill, saving *_RUN.ipynb.
for nb in "${NOTEBOOKS[@]}"; do
  IN="${NB_DIR}/${nb}.ipynb"
  OUT="${NB_DIR}/${nb}_RUN.ipynb"

  echo "→ Executing ${nb}.ipynb"
  papermill "$IN" "$OUT"

done

echo "✓ Pipeline complete! All notebooks executed successfully."