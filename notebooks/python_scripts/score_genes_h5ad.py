#!/usr/bin/env python
"""Score gene sets on a DSPIN unit h5ad; print TSV: cell_id, sample_id, score_pos, score_signed"""
import sys
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

def strip(g):
    return g[:-7] if g.endswith("-GRCm39") else g

def main():
    h5, meta_csv, pos_csv, neg_csv, out_csv = sys.argv[1:6]
    pos = [g.strip() for g in open(pos_csv) if g.strip()]
    neg = [g.strip() for g in open(neg_csv) if g.strip()] if neg_csv != "-" else []
    A = sc.read_h5ad(h5)
    genes = [strip(str(g)) for g in A.var_names]
    # uniquify
    seen = {}
    uniq = []
    for g in genes:
        if g in seen:
            seen[g] += 1
            uniq.append(f"{g}.{seen[g]}")
        else:
            seen[g] = 0
            uniq.append(g)
    A.var_names = uniq
    A.var_names_make_unique()
    pos_i = [g for g in pos if g in A.var_names]
    neg_i = [g for g in neg if g in A.var_names]
    if pos_i:
        sc.tl.score_genes(A, gene_list=pos_i, score_name="score_pos", use_raw=False)
    else:
        A.obs["score_pos"] = 0.0
    if neg_i:
        sc.tl.score_genes(A, gene_list=neg_i, score_name="score_neg", use_raw=False)
    else:
        A.obs["score_neg"] = 0.0
    meta = pd.read_csv(meta_csv)
    n = min(A.n_obs, len(meta))
    out = pd.DataFrame({
        "cell_id": meta["cell_id"].astype(str).values[:n],
        "sample_id": meta["sample_id"].astype(str).values[:n],
        "score_pos_only": np.asarray(A.obs["score_pos"])[:n],
        "score_signed": (np.asarray(A.obs["score_pos"]) - np.asarray(A.obs["score_neg"]))[:n],
    })
    out.to_csv(out_csv, index=False)
    print(f"OK n={n} pos={len(pos_i)} neg={len(neg_i)}", file=sys.stderr)

if __name__ == "__main__":
    main()
