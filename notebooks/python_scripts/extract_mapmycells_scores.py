#!/usr/bin/env python
"""Extract per-cell MapMyCells mapping-quality scores (+ depth covariates).

Flattens the MapMyCells *extended* JSON output (the CSV only carries the
bootstrap probability) into one compact parquet per modality, joined with the
sample/ROI metadata and per-cell depth covariates needed for QC modelling.

Per taxonomy level (class / subclass / supertype / cluster) we keep:
  <lvl>_label, <lvl>_name
  <lvl>_boot_prob   bootstrapping_probability (fraction of 100 bootstrap
                    iterations voting for the assigned node, conditional on
                    the parent assignment)
  <lvl>_agg_prob    aggregate_probability (product of boot_prob down the tree)
  <lvl>_avg_corr    avg_correlation to the assigned node's reference profile
  <lvl>_ru1_prob    probability of the top runner-up
  <lvl>_ru1_corr    correlation of the top runner-up
  <lvl>_ru_max_corr max correlation among the (up to 5) runners-up
  (supertype was a --drop_level; its scores are back-filled from cluster votes)

Usage (spatialomics env):
  PY=/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python
  $PY notebooks/python_scripts/extract_mapmycells_scores.py --modality seqFISH
  $PY notebooks/python_scripts/extract_mapmycells_scores.py --modality SN   # ~25 GB RAM -> sbatch
"""
import argparse
import datetime as dt
import glob
import json
import os
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd

WKDIR = "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
OUT_DIR = f"{WKDIR}/data/interim/MapMyCells_QC"
SEQFISH_DIR = f"{WKDIR}/data/input/SpatialAnalysis_2024-07"
SN_JSON = f"{WKDIR}/data/interim/MapMyCells/output.json"
SN_H5AD = f"{WKDIR}/data/input/parse_bio_anndata/anndata_ensembl.h5ad"

LEVELS = {
    "CCN20230722_CLAS": "class",
    "CCN20230722_SUBC": "subclass",
    "CCN20230722_SUPT": "supertype",
    "CCN20230722_CLUS": "cluster",
}

# ROI -> microbiome group + bregma, copied from 01_QC_prep_SeqFISH.qmd
SEQFISH_GROUP = {
    "roi1_run1": "SPF", "roi2_run1": "WildR", "roi3_run1": "WildR", "roi4_run1": "SPF",
    "roi1_run2": "WildR", "roi2_run2": "SPF", "roi3_run2": "SPF", "roi4_run2": "WildR",
    "roi1_run3": "SPF", "roi2_run3": "WildR", "roi3_run3": "WildR", "roi4_run3": "SPF",
    "roi1_run4": "SPF", "roi2_run4": "WildR", "roi3_run4": "WildR", "roi4_run4": "SPF",
}
SEQFISH_BREGMA = {
    "roi1_run1": -1.85, "roi2_run1": -1.95, "roi3_run1": -2.45, "roi4_run1": -1.85,
    "roi1_run2": -1.65, "roi2_run2": -1.65, "roi3_run2": -2.4, "roi4_run2": -1.85,
    "roi1_run3": -1.85, "roi2_run3": -2.15, "roi3_run3": -1.95, "roi4_run3": -1.95,
    "roi1_run4": -1.35, "roi2_run4": -1.85, "roi3_run4": -1.85, "roi4_run4": -1.75,
}


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def flatten_results(results, name_mapper):
    """List of MapMyCells per-cell result dicts -> wide DataFrame."""
    cols = {"cell_id": [r["cell_id"] for r in results]}
    for key, lvl in LEVELS.items():
        names = name_mapper[key]
        lab, boot, agg, corr, ru1p, ru1c, rumax = [], [], [], [], [], [], []
        for r in results:
            x = r[key]
            lab.append(x["assignment"])
            boot.append(x["bootstrapping_probability"])
            agg.append(x["aggregate_probability"])
            corr.append(x["avg_correlation"])
            rp, rc = x.get("runner_up_probability") or [], x.get("runner_up_correlation") or []
            ru1p.append(rp[0] if rp else np.nan)
            ru1c.append(rc[0] if rc else np.nan)
            rumax.append(max(rc) if rc else np.nan)
        cols[f"{lvl}_label"] = lab
        cols[f"{lvl}_name"] = [names[l]["name"] if l in names else None for l in lab]
        cols[f"{lvl}_boot_prob"] = boot
        cols[f"{lvl}_agg_prob"] = agg
        cols[f"{lvl}_avg_corr"] = corr
        cols[f"{lvl}_ru1_prob"] = ru1p
        cols[f"{lvl}_ru1_corr"] = ru1c
        cols[f"{lvl}_ru_max_corr"] = rumax
    df = pd.DataFrame(cols)
    for c in df.columns:
        if c.endswith(("_prob", "_corr")):
            df[c] = df[c].astype("float32")
    return df


def load_json(path):
    log(f"reading {path} ({os.path.getsize(path) / 1e9:.2f} GB)")
    with open(path) as fh:
        j = json.load(fh)
    meta = {
        "mmc_version": j["metadata"].get("version"),
        "mmc_timestamp": j["metadata"].get("timestamp"),
        "n_unmapped_genes": j.get("n_unmapped_genes"),
    }
    return j["results"], j["taxonomy_tree"]["name_mapper"], meta


def run_seqfish():
    frames = []
    for d in sorted(glob.glob(f"{SEQFISH_DIR}/export_*_july2024/export_*_july2024")):
        m = re.search(r"export_(run\d)_(roi\d)_(block\d)_july2024$", d)
        run, roi, block = m.groups()
        sample = f"{roi}_{run}"
        js = f"{d}/mapmycells_output.json"
        if not os.path.exists(js):
            log(f"MISSING mapmycells json for {sample}: {js}")
            continue
        results, nm, meta = load_json(js)
        df = flatten_results(results, nm)
        # depth covariates from the raw (symbol) count matrix + cell attributes
        cbg = pd.read_csv(f"{d}/cell_by_gene_matrix_{roi}_{block}.csv", index_col=0)
        attr = pd.read_csv(f"{d}/cell_attributes_{roi}_{block}.csv")
        cov = pd.DataFrame({
            "label": cbg.index.astype(int),
            "n_counts": cbg.sum(axis=1).to_numpy(),
            "n_genes": (cbg > 0).sum(axis=1).to_numpy(),
        }).merge(attr[["label", "area"]], on="label", how="left")
        cov["cell_id"] = f"{roi}_{block}_" + cov["label"].astype(str)
        n0 = len(df)
        df = df.merge(cov.drop(columns="label"), on="cell_id", how="left")
        assert len(df) == n0, "covariate join duplicated rows"
        n_miss = df["n_counts"].isna().sum()
        if n_miss:
            log(f"WARNING {sample}: {n_miss} cells without count covariates")
        df.insert(0, "sample_id", sample)
        df.insert(1, "group", SEQFISH_GROUP[sample])
        df.insert(2, "run", run)
        df.insert(3, "roi", roi)
        df.insert(4, "block", block)
        df.insert(5, "bregma", SEQFISH_BREGMA[sample])
        df["n_panel_genes"] = cbg.shape[1]
        df["mmc_timestamp"] = meta["mmc_timestamp"]
        log(f"{sample}: {len(df):,} cells, {cbg.shape[1]} genes, group={SEQFISH_GROUP[sample]}")
        frames.append(df)
        del results, cbg
    return pd.concat(frames, ignore_index=True)


def h5_obs_col(obs, k):
    g = obs[k]
    if isinstance(g, h5py.Group):
        cats = g["categories"][:].astype(str)
        return cats[g["codes"][:]]
    v = g[:]
    return v.astype(str) if v.dtype.kind in "SO" else v


def run_sn():
    log(f"reading obs from {SN_H5AD}")
    with h5py.File(SN_H5AD, "r") as f:
        obs = f["obs"]
        obs_df = pd.DataFrame({
            "cell_id": h5_obs_col(obs, obs.attrs["_index"]),
            "sample_id": h5_obs_col(obs, "sample"),
            "n_counts": h5_obs_col(obs, "tscp_count"),
            "n_genes": h5_obs_col(obs, "gene_count"),
        })
    assert obs_df["cell_id"].is_unique, "SN barcodes are not unique"
    results, nm, meta = load_json(SN_JSON)
    df = flatten_results(results, nm)
    del results
    n0 = len(df)
    df = df.merge(obs_df, on="cell_id", how="left")
    assert len(df) == n0
    log(f"{df['sample_id'].isna().sum()} mapped cells without obs metadata")
    # keep only this project's WildR / SPF samples (h5ad also holds ASO/WT samples)
    df = df[df["sample_id"].str.match(r"^(WildR|SPF)_", na=False)].copy()
    parts = df["sample_id"].str.extract(r"^(WildR|SPF)_(.+)_(AMY|HYP)$")
    df.insert(1, "group", parts[0])
    df.insert(2, "animal", parts[0] + "_" + parts[1])
    df.insert(3, "region", parts[2])
    df["mmc_timestamp"] = meta["mmc_timestamp"]
    log(f"SN: kept {len(df):,} WildR/SPF nuclei from {df['sample_id'].nunique()} samples")
    cols = ["sample_id"] + [c for c in df.columns if c != "sample_id"]
    return df[cols]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modality", choices=["seqFISH", "SN"], required=True)
    ap.add_argument("--date", default=dt.date.today().isoformat())
    a = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    out = f"{OUT_DIR}/mapmycells_scores_{a.modality}_{a.date}.parquet"
    if os.path.exists(out):
        sys.exit(f"refusing to overwrite existing {out}")
    df = run_seqfish() if a.modality == "seqFISH" else run_sn()
    df.to_parquet(out, index=False)
    log(f"wrote {out}: {df.shape[0]:,} x {df.shape[1]}")


if __name__ == "__main__":
    main()
