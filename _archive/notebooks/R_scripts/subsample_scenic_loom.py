#!/usr/bin/env python
"""Stratified cell subsampling for an existing pySCENIC loom + meta pair.

Stratifies by microbiome x subclass_name, proportionally allocating the
target cell count (largest-remainder method) so WildR/SPF balance and
subtype diversity are preserved. Keeps the full gene set as-is (no re-filtering).

Usage:
  python subsample_scenic_loom.py \
    --loom snRNASeq_AMY_CS20230722_CLAS_01.loom \
    --meta snRNASeq_AMY_CS20230722_CLAS_01_meta.csv \
    --target-n 45000 --seed 0 \
    --backup-suffix _full87k
"""
import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import loompy
from scipy import sparse

parser = argparse.ArgumentParser()
parser.add_argument("--loom", required=True, help="loom basename or path, in-place input")
parser.add_argument("--meta", required=True, help="meta csv basename or path, in-place input")
parser.add_argument("--input-dir", default=None, help="directory holding --loom/--meta if given as basenames")
parser.add_argument("--target-n", type=int, required=True)
parser.add_argument("--min-per-group", type=int, default=5)
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--backup-suffix", default="_full")
args = parser.parse_args()

input_dir = Path(args.input_dir) if args.input_dir else None


def resolve(p):
    p = Path(p)
    if p.is_absolute() or p.exists():
        return p
    return (input_dir / p) if input_dir else p


loom_path = resolve(args.loom)
meta_path = resolve(args.meta)
stem = loom_path.name[: -len(".loom")]
backup_loom = loom_path.with_name(f"{stem}{args.backup_suffix}.loom")
backup_meta = meta_path.with_name(meta_path.stem + f"{args.backup_suffix}.csv")

# Resumable: a prior interrupted run may have already moved originals to the backup path.
source_loom = backup_loom if backup_loom.exists() else loom_path
source_meta = backup_meta if backup_meta.exists() else meta_path
assert source_loom.exists(), source_loom
assert source_meta.exists(), source_meta

meta = pd.read_csv(source_meta)
meta["cell_id"] = meta["cell_id"].astype(str)
n_total = len(meta)
assert args.target_n < n_total, f"target_n {args.target_n} >= existing {n_total}, nothing to do"

strata_col = "microbiome" if "subclass_name" not in meta.columns else None
group_cols = ["microbiome", "subclass_name"] if "subclass_name" in meta.columns else ["microbiome"]

group_sizes = meta.groupby(group_cols).size()
rng = np.random.default_rng(args.seed)

raw_alloc = (group_sizes / n_total * args.target_n)
floor_alloc = np.floor(raw_alloc).astype(int)
floor_alloc = np.minimum(floor_alloc, group_sizes)
floor_alloc = np.maximum(floor_alloc, np.minimum(args.min_per_group, group_sizes))

deficit = args.target_n - floor_alloc.sum()
remainder = (raw_alloc - floor_alloc).sort_values(ascending=False)
alloc = floor_alloc.copy()
idx = 0
groups_cycle = list(remainder.index)
while deficit > 0 and groups_cycle:
    g = groups_cycle[idx % len(groups_cycle)]
    if alloc[g] < group_sizes[g]:
        alloc[g] += 1
        deficit -= 1
    idx += 1
    if idx > 10 * len(groups_cycle) and deficit > 0:
        break  # all groups saturated; accept slightly under target

print(f"Target n={args.target_n}, allocated n={alloc.sum()} across {len(alloc)} strata")

sampled_ids = []
for g, n_g in alloc.items():
    mask = (meta[group_cols] == pd.Series(g, index=group_cols)).all(axis=1) if len(group_cols) > 1 else (meta[group_cols[0]] == g)
    pool = meta.loc[mask, "cell_id"].to_numpy()
    chosen = rng.choice(pool, size=int(n_g), replace=False) if n_g < len(pool) else pool
    sampled_ids.extend(chosen.tolist())

sampled_ids = set(sampled_ids)
meta_sub = meta[meta["cell_id"].isin(sampled_ids)].reset_index(drop=True)
print(f"Subsampled meta rows: {len(meta_sub)} (from {n_total})")
print(meta_sub.groupby("microbiome").size())

# --- back up originals (if not already done by a prior interrupted run) ---
if source_loom is loom_path:
    shutil.move(str(loom_path), str(backup_loom))
    print("Backed up loom ->", backup_loom)
else:
    print("Backup loom already present (resuming):", backup_loom)

if source_meta is meta_path:
    shutil.move(str(meta_path), str(backup_meta))
    print("Backed up meta ->", backup_meta)
else:
    print("Backup meta already present (resuming):", backup_meta)

# Fancy-indexing columns directly on a gzip-chunked loom dataset is pathologically
# slow (touches ~every chunk individually); bulk-read then slice in memory instead.
with h5py.File(str(backup_loom), "r") as f:
    cell_ids = f["col_attrs/CellID"][:].astype(str)
    keep_mask = np.isin(cell_ids, list(sampled_ids))
    keep_idx = np.where(keep_mask)[0]
    print(f"Loom cells matched: {keep_mask.sum()} / {len(cell_ids)}")

    print("Bulk-reading full matrix into memory...")
    full_matrix = f["matrix"][:, :]

    def decode(arr):
        # h5py returns vlen-string datasets as raw bytes; loompy.create mis-handles
        # bytes (stringifies to literal "b'...'"), so decode to str up front.
        if arr.dtype == object:
            return np.array(
                [x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object
            )
        return arr

    row_attrs = {k: decode(f[f"row_attrs/{k}"][:]) for k in f["row_attrs"].keys()}
    col_attrs_full = {k: decode(f[f"col_attrs/{k}"][:]) for k in f["col_attrs"].keys()}

matrix_sub = sparse.csc_matrix(full_matrix[:, keep_idx])
col_attrs = {k: v[keep_idx] for k, v in col_attrs_full.items()}

loompy.create(str(loom_path), matrix_sub, row_attrs, col_attrs)
meta_sub.to_csv(meta_path, index=False)
print("Wrote subsampled loom ->", loom_path, "shape", matrix_sub.shape)
print("Wrote subsampled meta ->", meta_path)
