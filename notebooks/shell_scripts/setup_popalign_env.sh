#!/bin/bash
# One-time: create the `popalign` conda env (popAlign is not importable in any
# existing env). Additive only — does not touch spatialomics.
# popAlign itself is NOT pip-installed: jobs put the lab fork on PYTHONPATH
#   /resnick/groups/mthomson/jboktor/spatial_genomics/popalign  (jboktor/popalign @ d9a0c1a)
#
# Usage: bash notebooks/shell_scripts/setup_popalign_env.sh
set -euo pipefail
MAMBA="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/bin/mamba"
ENV="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/popalign"
POPALIGN_SRC="/resnick/groups/mthomson/jboktor/spatial_genomics/popalign"

if [ ! -x "${ENV}/bin/python" ]; then
  "${MAMBA}" create -y -n popalign -c conda-forge \
    python=3.9 "numpy<2" scipy pandas scikit-learn matplotlib seaborn \
    plotly ipywidgets adjusttext fastcluster umap-learn
fi
PYTHONPATH="${POPALIGN_SRC}" "${ENV}/bin/python" -c "import popalign as PA; print('popalign OK', PA.__file__)"
