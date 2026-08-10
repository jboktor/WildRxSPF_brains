#!/bin/bash
# Render DSPIN Analysis with the broken-conda Quarto layout fixed via env vars.
# Usage: bash notebooks/shell_scripts/render_dspin_analysis.sh
set -euo pipefail
ENV="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics"
export QUARTO_SHARE_PATH="${ENV}/share/quarto"
export PATH="${ENV}/bin:${PATH}"
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
cd "${WKDIR}/notebooks"
echo "[$(date)] QUARTO_SHARE_PATH=${QUARTO_SHARE_PATH}"
"${ENV}/bin/quarto" render 0X_DSPIN_Analysis.qmd --to html
ls -lh 0X_DSPIN_Analysis.html
