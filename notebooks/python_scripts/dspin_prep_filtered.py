#!/usr/bin/env python
"""HVG + log1p prep for DSPIN slices listed in a manifest CSV."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

MIN_CELLS_GENE = 10
N_HVG = 2000


def prepare_slice(h5ad_raw: Path, h5ad_out: Path) -> dict:
    adata = ad.read_h5ad(h5ad_raw)
    for col in ["sample_id", "batch", "if_control"]:
        if col not in adata.obs:
            raise KeyError(f"{h5ad_raw}: missing obs[{col}]")
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_GENE)
    adata.layers["counts"] = adata.X.copy()
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts", subset=False
        )
    except Exception as e:
        print(f"  seurat_v3 HVG failed ({e}); using seurat flavor")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat", subset=False)
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    h5ad_out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(h5ad_out, compression="gzip")
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_SPF": int((adata.obs["microbiome"].astype(str) == "SPF").sum()),
        "n_WildR": int((adata.obs["microbiome"].astype(str) == "WildR").sum()),
        "filtered": str(h5ad_out),
        "prep_status": "ok",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    man = pd.read_csv(args.manifest)
    dspin_root = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/data/interim/dspin")
    rows = []
    for _, row in man.iterrows():
        unit_path = None
        for key in ("unit_path", "slice_path"):
            if key in row and pd.notna(row[key]):
                unit_path = Path(row[key])
                break
        if unit_path is None:
            if "supertype_label" not in row or pd.isna(row["supertype_label"]):
                raise KeyError("manifest needs unit_path or supertype_label")
            unit_path = dspin_root / f"{row['tissue']}__{row['supertype_label']}"
        raw = unit_path / "raw_counts.h5ad"
        out = unit_path / "filtered.h5ad"
        print(f"=== {unit_path.name} ===", flush=True)
        if out.exists() and not args.force:
            print("  exists, skip", flush=True)
            rows.append({**row.to_dict(), "prep_status": "exists", "filtered": str(out)})
            continue
        if not raw.exists():
            print("  MISSING raw", raw, flush=True)
            rows.append({**row.to_dict(), "prep_status": "missing_raw"})
            continue
        info = prepare_slice(raw, out)
        print(" ", info, flush=True)
        rows.append({**row.to_dict(), **info})

    out_man = Path(args.manifest)
    prep_path = out_man.with_name(out_man.stem + "_prep_status.csv")
    pd.DataFrame(rows).to_csv(prep_path, index=False)
    updated = pd.DataFrame(rows)
    join_keys = ["tissue", "supertype_label"] if "supertype_label" in man.columns else ["tissue", "supertype_name"]
    keep = [c for c in man.columns if c not in ("n_genes", "prep_status", "filtered")]
    merge_cols = join_keys + ["n_genes", "prep_status", "filtered"]
    merged = man[keep].merge(updated[merge_cols], on=join_keys, how="left")
    merged.to_csv(out_man, index=False)
    print("Updated", out_man, "and", prep_path)
    print(merged["prep_status"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
