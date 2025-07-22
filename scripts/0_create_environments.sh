#!/usr/bin/env bash
set -euo pipefail

###############################################################################
#  Helper: create-or-update an env from the right YAML                       #
###############################################################################
make_env () {
  local env_name="$1" mac_yaml="$2" linux_yaml="$3"

  case "$(uname -s)" in
    Darwin*) yaml="$mac_yaml"   ;;
    Linux*)  yaml="$linux_yaml" ;;
    *) echo "Unsupported OS: $(uname -s)"; exit 1 ;;
  esac
  echo "→ Using YAML for $env_name: $(basename "$yaml")"

  if command -v mamba >/dev/null 2>&1; then
    CONDA_EXEC=mamba
  else
    CONDA_EXEC=conda
  fi

  if $CONDA_EXEC env list | grep -qE "^${env_name}\s"; then
    echo "→ Updating existing env: $env_name"
    $CONDA_EXEC env update -n "$env_name" -f "$yaml" --prune
  else
    echo "→ Creating env: $env_name"
    $CONDA_EXEC env create -f "$yaml"
  fi
}

###############################################################################
#  Locate YAML directory and R post-install script                            #
###############################################################################
script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
env_dir="${script_dir}/../environments"
r_post="${script_dir}/install_remaining_pkgs.R"

###############################################################################
#  1) Build / update the *R* environment                                      #
###############################################################################
make_env \
  env_deconv_R \
  "${env_dir}/env_deconv_R_mac.yml" \
  "${env_dir}/env_deconv_R_linux.yml"

echo "→ Installing remaining R packages (if any)…"
conda run -n env_deconv_R Rscript "$r_post"
echo "✓ Finished bootstrap for env_deconv_R"

###############################################################################
#  2) Build / update the *Python* environment                                 #
###############################################################################
make_env \
  env_deconv_py \
  "${env_dir}/env_deconv_mac.yml" \
  "${env_dir}/env_deconv_linux.yml"

echo "✓ Finished bootstrap for env_deconv_py"

echo "✓  All environments created and up-to-date!"

