#!/usr/bin/env python
"""Export one pySCENIC loom per tissue x Allen subclass (or supertype), sliced from the class looms.

The class looms (data/input/scenic/snRNASeq/snRNASeq_{tissue}_{class_label}.loom)
already hold raw UMI counts for every nucleus of the class, with QC-passing
nuclei (>= 200 genes) and bare gene symbols. Slicing them avoids reloading the
full Seurat object. For AMY CLAS_01 the full 87k-nucleus loom (*_full87k) is
used, not the 45k subsample that was needed to fit the class GRN in walltime.

Levels: --level subclass (stem snRNASeq_{tissue}_SUBC_{NNN}) or
--level supertype (stem snRNASeq_{tissue}_SUPT_{NNNN}).

Eligibility (per tissue x subclass/supertype):
  >= 100 nuclei, >= 30 nuclei per microbiome group, and >= 3 animals per group
  with >= 10 nuclei each (the floor of the animal-level test).
Genes are re-filtered to those detected in >= 10 nuclei of the split.
Subclasses larger than --max-cells are stratified-subsampled (microbiome x
supertype, largest remainder); none exceeded 45k at the time of writing.

Outputs (data/interim/scenic/looms_<level>/):
  snRNASeq_{tissue}_{SUBC|SUPT}_{id}.loom + _meta.csv
  manifest_<level>.csv
  looms_<level>_{small,large}.txt   (absolute loom paths, one per line; tiers
                                      for the SLURM arrays, split at --large-n)

Usage:
  python export_scenic_subclass_looms.py [--level subclass|supertype] [--max-cells 45000] [--large-n 5000]
"""
import argparse
import re
from pathlib import Path

import loompy
import numpy as np
import pandas as pd

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
CLASS_DIR = WKDIR / "data/input/scenic/snRNASeq"

p = argparse.ArgumentParser()
p.add_argument("--level", choices=["subclass", "supertype"], default="subclass")
p.add_argument("--min-cells", type=int, default=100)
p.add_argument("--min-per-group", type=int, default=30)
p.add_argument("--min-animals", type=int, default=3)
p.add_argument("--min-cells-animal", type=int, default=10)
p.add_argument("--min-cells-gene", type=int, default=10)
p.add_argument("--max-cells", type=int, default=45000)
p.add_argument("--large-n", type=int, default=5000)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--glob", default="snRNASeq_*_CLAS_??.loom", help="class looms to slice (testing)")
p.add_argument("--out-dir", default=None, help="override output dir (testing)")
args = p.parse_args()

LEVEL = args.level
COL = f"{LEVEL}_name"
TAG = {"subclass": "SUBC", "supertype": "SUPT"}[LEVEL]
OUT_DIR = Path(args.out_dir) if args.out_dir else WKDIR / f"data/interim/scenic/looms_{LEVEL}"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def log(*msg):
    print(pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"), "-", *msg, flush=True)


def stratified_subsample(meta, n_target, rng):
    strata = meta.groupby(["microbiome", "supertype_name"]).indices
    sizes = pd.Series({k: len(v) for k, v in strata.items()})
    quota = sizes / sizes.sum() * n_target
    take = np.floor(quota).astype(int)
    rem = n_target - take.sum()
    take[(quota - take).sort_values(ascending=False).index[:rem]] += 1
    keep = np.concatenate([rng.choice(strata[k], size=min(take[k], sizes[k]), replace=False)
                           for k in sizes.index])
    return np.sort(keep)


class_looms = sorted(CLASS_DIR.glob(args.glob))
rng = np.random.default_rng(args.seed)
rows = []

for cl in class_looms:
    full = cl.with_name(cl.stem + "_full87k.loom")
    src = full if full.exists() else cl
    with loompy.connect(str(src), "r") as ds:
        ca = pd.DataFrame({k: ds.ca[k].astype(str) for k in
                           ["CellID", "sample", "microbiome", "tissue", "class_label",
                            "class_name", "subclass_name", "supertype_name"]})
        ca["animal"] = ca["sample"].str.replace(r"_(AMY|HYP)$", "", regex=True)
        genes = ds.ra["Gene"].astype(str)
        for subc, idx in ca.groupby(COL).indices.items():
            m = ca.iloc[idx]
            tissue = m["tissue"].iloc[0]
            subc_id = re.match(r"^(\d+)", subc).group(1)
            stem = f"snRNASeq_{tissue}_{TAG}_{subc_id}"
            n_by_group = m["microbiome"].value_counts()
            per_animal = m.groupby(["microbiome", "animal"]).size()
            n_anim = {g: int((per_animal.get(g, pd.Series(dtype=int)) >= args.min_cells_animal).sum())
                      for g in ["WildR", "SPF"]}
            row = dict(split_level=LEVEL, tissue=tissue, class_label=m["class_label"].iloc[0],
                       class_name=m["class_name"].iloc[0], subclass_name=m["subclass_name"].iloc[0],
                       supertype_name=m["supertype_name"].iloc[0] if LEVEL == "supertype" else None,
                       stem=stem,
                       n_cells_total=len(m), n_SPF=int(n_by_group.get("SPF", 0)),
                       n_WildR=int(n_by_group.get("WildR", 0)),
                       n_animals_WildR_ok=n_anim["WildR"], n_animals_SPF_ok=n_anim["SPF"],
                       source_loom=src.name)
            reason = None
            if len(m) < args.min_cells:
                reason = f"n<{args.min_cells}"
            elif min(row["n_SPF"], row["n_WildR"]) < args.min_per_group:
                reason = f"<{args.min_per_group} per group"
            elif min(n_anim.values()) < args.min_animals:
                reason = f"<{args.min_animals} animals/group with >={args.min_cells_animal} nuclei"
            if reason:
                rows.append({**row, "status": f"skipped: {reason}"})
                continue

            if len(idx) > args.max_cells:
                sub = stratified_subsample(m.reset_index(drop=True), args.max_cells, rng)
                idx = idx[sub]
                m = ca.iloc[idx]
                log(f"{stem}: subsampled {row['n_cells_total']} -> {len(idx)}")

            loom_path = OUT_DIR / f"{stem}.loom"
            meta_path = OUT_DIR / f"{stem}_meta.csv"
            if loom_path.exists() and loom_path.stat().st_size > 1000:
                log(f"exists, skip: {loom_path.name}")
                status = "exists"
            else:
                X = ds.sparse(cols=np.sort(idx)).tocsr()       # genes x cells, raw UMI
                m = ca.iloc[np.sort(idx)]
                gkeep = np.asarray((X > 0).sum(axis=1)).ravel() >= args.min_cells_gene
                X = X[gkeep]
                meta = m.rename(columns={"CellID": "cell_id"})[
                    ["cell_id", "sample", "microbiome", "tissue", "class_label", "class_name",
                     "subclass_name", "supertype_name"]]
                meta.to_csv(meta_path, index=False)
                col_attrs = {"CellID": m["CellID"].to_numpy(dtype=object)}
                for c in ["sample", "microbiome", "tissue", "class_label", "class_name",
                          "subclass_name", "supertype_name"]:
                    col_attrs[c] = m[c].to_numpy(dtype=object)
                loompy.create(str(loom_path), X.astype(np.float32),
                              {"Gene": genes[gkeep].astype(object)}, col_attrs)
                log(f"wrote {loom_path.name}: {X.shape[0]} genes x {X.shape[1]} nuclei")
                status = "wrote"
            rows.append({**row, "n_cells": len(idx), "loom": str(loom_path),
                         "meta": str(meta_path), "status": status})

man = pd.DataFrame(rows).sort_values(["tissue", "subclass_name", "stem"])
man.to_csv(OUT_DIR / f"manifest_{LEVEL}.csv", index=False)
ok = man[man["status"].isin(["wrote", "exists"])]
(OUT_DIR / f"looms_{LEVEL}_large.txt").write_text(
    "\n".join(ok.loc[ok["n_cells"] >= args.large_n, "loom"]) + "\n")
(OUT_DIR / f"looms_{LEVEL}_small.txt").write_text(
    "\n".join(ok.loc[ok["n_cells"] < args.large_n, "loom"]) + "\n")
log(f"{LEVEL} splits seen: {len(man)}; exported: {len(ok)} "
    f"(large >= {args.large_n}: {(ok['n_cells'] >= args.large_n).sum()}, "
    f"small: {(ok['n_cells'] < args.large_n).sum()})")
