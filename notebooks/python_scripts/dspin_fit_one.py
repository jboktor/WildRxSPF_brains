#!/usr/bin/env python
"""Fit one DSPIN program-level model for a manifest row (SLURM array friendly)."""
from __future__ import annotations

import argparse
import json
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

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


def unit_dir(row) -> Path:
    """Path key is tissue__{supertype_label} (never slugify name — '/' hazard)."""
    for key in ("unit_path", "slice_path"):
        if key in row and pd.notna(row[key]):
            return Path(row[key])
    if "supertype_label" in row and pd.notna(row["supertype_label"]):
        return DSPIN_DIR / f"{row['tissue']}__{row['supertype_label']}"
    raise KeyError("manifest row needs unit_path or supertype_label")


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
    sdir = unit_dir(row)
    h5ad = sdir / "filtered.h5ad"
    if not h5ad.exists():
        h5ad = sdir / "raw_counts.h5ad"
    if not h5ad.exists():
        return {"status": "error", "error": f"missing h5ad under {sdir}"}

    t0 = time.time()
    adata = ad.read_h5ad(h5ad)
    xmax = float(adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X))
    if xmax > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata = ensure_leiden(adata)

    save_path = str(sdir / "dspin_run")
    Path(save_path).mkdir(parents=True, exist_ok=True)

    label = row.get("supertype_label", "")
    name = row.get("supertype_name", "")
    status = {
        "tissue": row["tissue"],
        "supertype_label": None if pd.isna(label) else str(label),
        "supertype_name": None if pd.isna(name) else str(name),
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

    print(
        f"DSPIN init {row['tissue']} {label} ({name}) n={adata.n_obs} g={adata.n_vars}",
        flush=True,
    )
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
        tidy_rel = matrix_to_tidy(
            R, program_names[: R.shape[0]], [f"c{j}" for j in range(R.shape[1])], "relative_h"
        )
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
    title_name = row.get("supertype_name", row.get("supertype_label", ""))
    save_network_preview(net, list(net_df.index), png, title=f"{row['tissue']} | {title_name}")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig_stem = row.get("supertype_label") or sdir.name.split("__", 1)[-1]
    fig_copy = FIG_DIR / f"{row['tissue']}__{fig_stem}_network_preview.png"
    fig_copy.write_bytes(png.read_bytes())

    status["status"] = "ok"
    status["elapsed_s"] = round(time.time() - t0, 1)
    status["finished"] = datetime.now(timezone.utc).isoformat()
    (sdir / "fit_status.json").write_text(json.dumps(status, indent=2))
    print("DONE", status, flush=True)
    return status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--num-spin", type=int, default=20)
    ap.add_argument("--num-repeat", type=int, default=10)
    ap.add_argument("--skip-ok", action="store_true", help="Skip if fit_status.json already ok")
    args = ap.parse_args()

    man = pd.read_csv(args.manifest)
    if args.index < 0 or args.index >= len(man):
        raise SystemExit(f"index {args.index} out of range 0..{len(man)-1}")
    row = man.iloc[args.index]
    sdir = unit_dir(row)
    status_path = sdir / "fit_status.json"
    if args.skip_ok and status_path.exists():
        try:
            prev = json.loads(status_path.read_text())
            if prev.get("status") == "ok" and int(prev.get("num_spin", 0)) >= args.num_spin:
                print(
                    "SKIP already ok:",
                    row["tissue"],
                    row.get("supertype_label"),
                    flush=True,
                )
                return
        except Exception:
            pass

    try:
        fit_one(row, num_spin=args.num_spin, num_repeat=args.num_repeat)
    except Exception as e:
        sdir.mkdir(parents=True, exist_ok=True)
        err = {
            "status": "error",
            "tissue": row.get("tissue"),
            "supertype_label": None if pd.isna(row.get("supertype_label")) else str(row.get("supertype_label")),
            "supertype_name": None if pd.isna(row.get("supertype_name")) else str(row.get("supertype_name")),
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
        }
        status_path.write_text(json.dumps(err, indent=2))
        print("ERROR", err["error"], flush=True)
        raise


if __name__ == "__main__":
    main()
