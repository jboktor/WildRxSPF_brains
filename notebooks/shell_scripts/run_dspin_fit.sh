#!/bin/bash
#SBATCH --job-name="dspin_fit"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_fit_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_fit_%j.err"
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Usage:
#   sbatch notebooks/shell_scripts/run_dspin_fit.sh          # full interactive (≤5)
#   sbatch notebooks/shell_scripts/run_dspin_fit.sh smoke    # one slice smoke test

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
MODE="${1:-full}"

mkdir -p "${WKDIR}/.cluster_runs/dspin" "${WKDIR}/figures/DSPIN/snRNASeq"
cd "${WKDIR}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

SMOKE=0
if [[ "${MODE}" == "smoke" ]]; then
  SMOKE=1
fi

echo "[$(date)] MODE=${MODE} SMOKE=${SMOKE}"

"${PY}" - << PY
from pathlib import Path
import json
import time
import traceback
from datetime import datetime, timezone

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from dspin.dspin import DSPIN

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
DSPIN_DIR = WKDIR / "data/interim/dspin"
FIG_DIR = WKDIR / "figures/DSPIN/snRNASeq"
FIG_DIR.mkdir(parents=True, exist_ok=True)
MANIFEST = DSPIN_DIR / "pilot_manifest.csv"

SMOKE_TEST = bool(int("${SMOKE}"))
NUM_SPIN_SMOKE, NUM_REPEAT_SMOKE = 15, 3
NUM_SPIN_FULL, NUM_REPEAT_FULL = 20, 10

manifest = pd.read_csv(MANIFEST)
manifest = manifest.query("prep_status == 'ok'" if "prep_status" in manifest else "True").copy()
if SMOKE_TEST:
    # smallest slice for smoke
    if "n_cells" in manifest.columns:
        run_df = manifest.sort_values("n_cells").iloc[[0]].copy()
    else:
        run_df = manifest.iloc[[0]].copy()
else:
    run_df = manifest.query("run_tier == 'interactive'").copy() if "run_tier" in manifest else manifest.head(5)

print(run_df[["tissue", "supertype_name"]].to_string(index=False))


def slice_dir(row) -> Path:
    if "slice_path" in row and pd.notna(row["slice_path"]):
        return Path(row["slice_path"])
    return DSPIN_DIR / f"{row['tissue']}__{row['supertype_slug']}"


def matrix_to_tidy(mat, row_names, col_names, value_name):
    arr = np.asarray(mat)
    df = pd.DataFrame(arr, index=list(row_names)[: arr.shape[0]], columns=list(col_names)[: arr.shape[1]])
    return (
        df.rename_axis("program_id")
        .reset_index()
        .melt(id_vars="program_id", var_name="sample_id", value_name=value_name)
    )


def save_network_preview(network, program_names, out_png: Path, title: str):
    fig, ax = plt.subplots(figsize=(8, 7))
    M = np.asarray(network, dtype=float)
    vmax = np.nanmax(np.abs(M)) if M.size else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    im = ax.imshow(M, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_title(title)
    labels = [str(p)[:18] for p in program_names]
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def ensure_leiden(adata: ad.AnnData) -> ad.AnnData:
    if "leiden" in adata.obs.columns:
        return adata
    n_comps = min(30, max(2, adata.n_obs - 1), max(2, adata.n_vars - 1))
    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)))
    try:
        sc.tl.leiden(adata, key_added="leiden", flavor="igraph", directed=False, n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, key_added="leiden")
    return adata


def fit_one(row, num_spin: int, num_repeat: int) -> dict:
    sdir = slice_dir(row)
    h5ad = sdir / "filtered.h5ad"
    if not h5ad.exists():
        h5ad = sdir / "raw_counts.h5ad"
    if not h5ad.exists():
        return {"status": "error", "error": f"missing h5ad under {sdir}"}

    t0 = time.time()
    adata = ad.read_h5ad(h5ad)
    xmax = float(np.asarray(adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)))
    if xmax > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata = ensure_leiden(adata)

    save_path = str(sdir / "dspin_run")
    Path(save_path).mkdir(parents=True, exist_ok=True)

    status = {
        "tissue": row["tissue"],
        "supertype_name": row["supertype_name"],
        "num_spin": num_spin,
        "num_repeat": num_repeat,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_SPF": int((adata.obs["microbiome"].astype(str) == "SPF").sum()),
        "n_WildR": int((adata.obs["microbiome"].astype(str) == "WildR").sum()),
        "started": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "method": None,
        "error": None,
    }

    print(f"DSPIN init {row['tissue']} {row['supertype_name']} n={adata.n_obs} g={adata.n_vars}", flush=True)
    model = DSPIN(adata, save_path, num_spin=num_spin)
    print("gene_program_discovery...", flush=True)
    model.gene_program_discovery(num_repeat=num_repeat, seed=0, cluster_key="leiden")

    try:
        print("network_inference auto...", flush=True)
        model.network_inference(sample_id_key="sample_id", method="auto")
        status["method"] = "auto"
    except Exception as e:
        print("auto failed, pseudo_likelihood:", e, flush=True)
        model.network_inference(sample_id_key="sample_id", method="pseudo_likelihood")
        status["method"] = "pseudo_likelihood"

    print("response_relative_to_control...", flush=True)
    model.response_relative_to_control(
        sample_id_key="sample_id",
        if_control_key="if_control",
        batch_key="batch",
    )

    program_names = list(getattr(model, "name_list_short", None) or [f"P{i}" for i in range(num_spin)])
    sample_list = list(model.sample_list)

    rel = np.asarray(model.relative_responses)
    if rel.ndim == 2 and rel.shape[0] == len(sample_list):
        tidy_rel = matrix_to_tidy(rel.T, program_names, sample_list, "relative_h")
    elif rel.ndim == 2 and rel.shape[1] == len(sample_list):
        tidy_rel = matrix_to_tidy(rel, program_names, sample_list, "relative_h")
    else:
        R = np.atleast_2d(rel)
        tidy_rel = matrix_to_tidy(R, program_names[: R.shape[0]], [f"c{j}" for j in range(R.shape[1])], "relative_h")
    tidy_rel.to_csv(sdir / "relative_responses.csv", index=False)

    resp = np.asarray(model.responses)
    if resp.ndim == 2 and resp.shape[0] == len(sample_list):
        matrix_to_tidy(resp.T, program_names, sample_list, "h").to_csv(sdir / "responses.csv", index=False)
    elif resp.ndim == 2:
        matrix_to_tidy(resp, program_names[: resp.shape[0]], sample_list[: resp.shape[1]], "h").to_csv(
            sdir / "responses.csv", index=False
        )

    net = np.asarray(model.network)
    net_df = pd.DataFrame(net, index=program_names[: net.shape[0]], columns=program_names[: net.shape[1]])
    net_df.to_csv(sdir / "network.csv")

    rows = []
    full_names = list(getattr(model, "name_list", program_names))
    for pname, full in zip(program_names, full_names):
        if "-" in str(full):
            genes = str(full).split("-", 1)[1].split(",")
            for g in genes:
                if g and g != "nan":
                    rows.append({"program_id": pname, "gene": g, "weight": np.nan})
    pd.DataFrame(rows).to_csv(sdir / "programs_genes.csv", index=False)

    sample_meta = (
        adata.obs.groupby(["sample_id", "microbiome", "batch", "if_control"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    sample_meta.to_csv(sdir / "sample_meta.csv", index=False)

    png = sdir / "network_preview.png"
    save_network_preview(net, list(net_df.index), png, title=f"{row['tissue']} | {row['supertype_name']}")
    fig_copy = FIG_DIR / f"{row['tissue']}__{row['supertype_slug']}_network_preview.png"
    fig_copy.write_bytes(png.read_bytes())

    status["status"] = "ok"
    status["elapsed_s"] = round(time.time() - t0, 1)
    status["finished"] = datetime.now(timezone.utc).isoformat()
    (sdir / "fit_status.json").write_text(json.dumps(status, indent=2))
    print("DONE", status, flush=True)
    return status


results = []
for _, row in run_df.iterrows():
    print("=" * 60, flush=True)
    num_spin = NUM_SPIN_SMOKE if SMOKE_TEST else NUM_SPIN_FULL
    num_repeat = NUM_REPEAT_SMOKE if SMOKE_TEST else NUM_REPEAT_FULL
    try:
        results.append(fit_one(row, num_spin=num_spin, num_repeat=num_repeat))
    except Exception as e:
        sdir = slice_dir(row)
        sdir.mkdir(parents=True, exist_ok=True)
        err = {
            "status": "error",
            "tissue": row.get("tissue"),
            "supertype_name": row.get("supertype_name"),
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
        }
        (sdir / "fit_status.json").write_text(json.dumps(err, indent=2))
        print("ERROR", err["error"], flush=True)
        results.append(err)

res_df = pd.DataFrame(results)
out = DSPIN_DIR / ("fit_summary_smoke.csv" if SMOKE_TEST else "fit_summary.csv")
res_df.to_csv(out, index=False)
print("Wrote", out, flush=True)
print(res_df.to_string(index=False), flush=True)
PY

echo "[$(date)] fit script finished"
