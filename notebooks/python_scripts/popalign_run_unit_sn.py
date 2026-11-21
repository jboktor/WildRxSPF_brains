#!/usr/bin/env python
"""Run popAlign on one snRNA-seq analysis unit (tissue × class, subclass or supertype).

One manifest row = one unit = one SLURM array task. Samples are animals
(6 SPF + 6 WildR per tissue); SPF animals are the popAlign control replicates.

Pipeline
  1. popAlign core: load_multiplexed → normalize → gene filter → oNMF feature
     space → per-animal GMMs (+ bootstrap replicates).
  2. Unit "atlas" GMM: one GMM fit on a pool of up to N nuclei per animal from
     all 12 animals (6 vs 6, so groups stay balanced). K is chosen by held-out BIC
     starting at K = 1, so "no substructure" is a valid outcome. Stability follows
     clusterboot (Hennig 2007): 20 bootstraps of 80% of the pool, K re-selected in
     each, and each atlas component scored by the Jaccard overlap of its nuclei
     with the best-matching bootstrap component (>= 0.75 = stable).
  3. Native popAlign alignment of every animal GMM (and replicates) to the atlas
     (method='conservative'), plus native delta (Δw, Δμ, ΔΣ) and LLR ranking.
  4. Tidy per-animal tables for animal-level statistics in R
     (component fractions + per-animal component centroids), component marker
     genes, Allen-label composition, and the oNMF gene programs.

Usage
  popalign_run_unit_sn.py --manifest data/interim/popalign/snRNASeq/manifest.csv --index 0
  (index defaults to $SLURM_ARRAY_TASK_ID)
"""
import argparse
import json
import os
import pickle
import sys
import time
import traceback

# popAlign parallelises with multiprocessing.Pool (fork). Multithreaded BLAS in
# the parent deadlocks forked workers, so pin BLAS to 1 thread before numpy loads.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import matplotlib
matplotlib.use("Agg")
from matplotlib.legend import Legend
if not hasattr(Legend, "legendHandles"):  # popAlign plot_deltas uses the pre-3.9 name
    Legend.legendHandles = property(lambda self: self.legend_handles)
import numpy as np
import pandas as pd
from scipy import sparse as ss
from sklearn import mixture as smix

import multiprocessing as mp
from joblib import Parallel, delayed
from scipy.stats import rankdata
import popalign as PA

# popAlign's forked Pool workers hang inside build_gmms after the oNMF/GSEA steps
# (parent state not fork-safe). Workers from a clean forkserver avoid this; the
# installed popAlign source is left untouched.
PA.popalign.Pool = mp.get_context("forkserver").Pool

# ── parameters by level ──────────────────────────────────────────────────────
# K ranges are half-open (popAlign convention) and start at 1 (no forced split).
# nfeats spans the paper's 10-20 programs with headroom on both sides.
PARAMS = {
    "class": dict(onmf_ncells=5000, nfeats=list(range(6, 25, 2)), gmm_ks=(1, 16),
                  atlas_ks=(1, 21), atlas_cells_per_animal=1000),
    "subclass": dict(onmf_ncells=3000, nfeats=list(range(4, 23, 2)), gmm_ks=(1, 13),
                     atlas_ks=(1, 16), atlas_cells_per_animal=700),
    "supertype": dict(onmf_ncells=2000, nfeats=list(range(4, 21, 2)), gmm_ks=(1, 9),
                      atlas_ks=(1, 11), atlas_cells_per_animal=400),
}
# scaling_factor="median_umi": beta = unit median UMI per nucleus (paper: beta ~
# median total counts per cell, so log(beta*g/sum+1) stays on the raw-count scale).
COMMON = dict(scaling_factor="median_umi", gene_offset=1.0, onmf_nreps=3, onmf_niter=300,
              gmm_niters=3, gmm_training=0.7, gmm_nreplicates=10, criteria="bic",
              atlas_nboot=20, atlas_boot_frac=0.8, stable_jaccard=0.75,
              n_markers=200, marker_min_pct=0.10, marker_min_dpct=0.05,
              marker_min_log2fc=0.5, marker_min_auc=0.6, seed=42)


def log(*a):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]", *a, flush=True)


def _fit_one(Ct, Cv, k, it, reg_covar, seed):
    g = smix.GaussianMixture(n_components=k, covariance_type="full",
                             reg_covar=reg_covar, random_state=seed + 101 * it + k).fit(Ct)
    return g, dict(k=k, iter=it, bic=g.bic(Cv), aic=g.aic(Cv))


def fit_gmm(C, ks, niters, reg_covar, seed, valid_frac=0.3, n_jobs=1):
    """Pick K by held-out BIC (popAlign's criterion), return fitted GMM + IC table."""
    rng = np.random.default_rng(seed)
    idx = rng.permutation(C.shape[0])
    nv = int(len(idx) * valid_frac)
    Cv, Ct = C[idx[:nv]], C[idx[nv:]]
    grid = [(k, it) for k in range(ks[0], ks[1]) if k < Ct.shape[0] for it in range(niters)]
    if n_jobs > 1:
        res = Parallel(n_jobs=n_jobs, backend="loky")(delayed(_fit_one)(Ct, Cv, k, it, reg_covar, seed) for k, it in grid)
    else:
        res = [_fit_one(Ct, Cv, k, it, reg_covar, seed) for k, it in grid]
    rows = [r for _, r in res]
    best = res[int(np.argmin([r["bic"] for r in rows]))][0]
    return best, pd.DataFrame(rows)


def jaccard(a, b):
    u = np.logical_or(a, b).sum()
    return np.logical_and(a, b).sum() / u if u else np.nan


def boot_once(pool, atlas_lab, K, ks, niters, reg_covar, seed, frac):
    """clusterboot-style resample: refit with K re-selected, Jaccard per atlas component."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(pool.shape[0], int(pool.shape[0] * frac), replace=False)
    g, _ = fit_gmm(pool[idx], ks, niters, reg_covar, seed)
    a, b = atlas_lab[idx], g.predict(pool[idx])
    jac = [max(jaccard(a == j, b == c) for c in range(g.n_components)) if (a == j).any() else np.nan
           for j in range(K)]
    return g.n_components, jac


def auc_rows(X, rows, in_mask):
    """AUC (Mann-Whitney) of each gene row of csr X for in_mask vs rest."""
    n_in = in_mask.sum(); n_out = len(in_mask) - n_in
    out = np.empty(len(rows))
    for i, g in enumerate(rows):
        r = rankdata(X[g].toarray().ravel())
        out[i] = (r[in_mask].sum() - n_in * (n_in + 1) / 2) / (n_in * n_out)
    return out


def jd_matrix(g_test, g_ref):
    arr = np.zeros((g_test.n_components, g_ref.n_components))
    for i in range(g_test.n_components):
        for j in range(g_ref.n_components):
            arr[i, j] = PA.JeffreyDiv(g_test.means_[i], g_test.covariances_[i],
                                      g_ref.means_[j], g_ref.covariances_[j])
    return arr


def group_stats(X, labels, K):
    """Per-component mean and fraction-expressing for a genes × cells sparse matrix."""
    n = X.shape[1]
    ind = ss.csr_matrix((np.ones(n), (np.arange(n), labels)), shape=(n, K))
    counts = np.asarray(ind.sum(axis=0)).ravel()
    tot_sum = np.asarray(X.sum(axis=1)).ravel()
    tot_nz = np.asarray((X > 0).sum(axis=1)).ravel()
    s_in = np.asarray((X @ ind).todense())
    nz_in = np.asarray(((X > 0).astype(np.float32) @ ind).todense())
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_in = s_in / counts
        pct_in = nz_in / counts
        mean_out = (tot_sum[:, None] - s_in) / (n - counts)
        pct_out = (tot_nz[:, None] - nz_in) / (n - counts)
    return mean_in, pct_in, mean_out, pct_out, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--index", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--skip-done", action="store_true")
    ap.add_argument("--no-pickle", action="store_true")
    ap.add_argument("--quick", action="store_true", help="smoke test: tiny parameter grid")
    ap.add_argument("--output-dir", default=None, help="override manifest output_dir")
    a = ap.parse_args()

    man = pd.read_csv(a.manifest)
    row = man.iloc[a.index]
    lv = row["level"]
    P = {**COMMON, **PARAMS[lv]}
    if a.quick:
        P.update(nfeats=[3, 5], onmf_nreps=1, onmf_niter=50, gmm_ks=(1, 4), atlas_ks=(1, 5),
                 gmm_niters=1, gmm_nreplicates=2, atlas_nboot=3)
    out = a.output_dir or row["output_dir"]
    inp = row["input_dir"]
    os.makedirs(out, exist_ok=True)
    status_f = os.path.join(out, "status.json")
    if a.skip_done and os.path.exists(status_f):
        if json.load(open(status_f)).get("status") == "ok":
            log("Already done:", row["unit_id"]); return
    np.random.seed(P["seed"])
    t0 = time.time()
    status = dict(unit_id=row["unit_id"], level=lv, tissue=row["tissue"], label=row["label"],
                  name=row["name"], params=P, steps={})
    log("Unit", row["unit_id"], "|", row["name"])

    # ── 1. popAlign core ────────────────────────────────────────────────────
    pop = PA.load_multiplexed(matrix=os.path.join(inp, "matrix.mtx"),
                              barcodes=os.path.join(inp, "barcodes.tsv"),
                              metafile=os.path.join(inp, "metadata.csv"),
                              genes=os.path.join(inp, "genes.tsv"),
                              controlstring="SPF", outputfolder=out)
    pop["ncores"] = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    meta = pop["meta"].reset_index(drop=True)
    log("Samples:", pop["order"], "| cells:", len(meta))

    # raw depth per nucleus (same order as pop["meta"] / cat_data)
    raw_umi = np.concatenate([np.asarray(pop["samples"][x]["M"].sum(axis=0)).ravel() for x in pop["order"]])
    beta = float(np.round(np.median(raw_umi))) if P["scaling_factor"] == "median_umi" else float(P["scaling_factor"])
    status["steps"]["beta"] = beta
    log("Scaling factor beta =", beta, "(median UMI per nucleus)")
    PA.normalize(pop, scaling_factor=beta)
    PA.plot_gene_filter(pop, offset=P["gene_offset"])
    # mt-/Rpl/Rps already removed at export (popAlign's prefix filter is human-only)
    PA.filter(pop, remove_ribsomal=False, remove_mitochondrial=False)
    status["steps"]["n_filtered_genes"] = int(len(pop["filtered_genes"]))

    PA.onmf(pop, ncells=P["onmf_ncells"], nfeats=P["nfeats"], nreps=P["onmf_nreps"], niter=P["onmf_niter"])
    PA.choose_featureset(pop, alpha=3, multiplier=3)
    m = int(pop["nfeats"])
    status["steps"]["n_features"] = m
    status["steps"]["n_features_at_bound"] = bool(m in (min(P["nfeats"]), max(P["nfeats"])))
    log("oNMF features:", m)

    PA.build_gmms(pop, ks=P["gmm_ks"], niters=P["gmm_niters"], training=P["gmm_training"],
                  nreplicates=P["gmm_nreplicates"], reg_covar="auto", types=None,
                  featuretype="onmf", criteria=P["criteria"])
    reg = pop["reg_covar"]

    # ── 2. atlas GMM on a group-balanced pool ───────────────────────────────
    # cap per animal (no equalising to the smallest animal; 6 vs 6 keeps groups balanced)
    rng = np.random.default_rng(P["seed"])
    ncores = pop["ncores"]
    parts = []
    for x in pop["order"]:
        Cx = PA.get_coeff(pop, x)
        parts.append(Cx[rng.choice(Cx.shape[0], min(P["atlas_cells_per_animal"], Cx.shape[0]), replace=False)])
    pool = np.vstack(parts)
    atlas, ic = fit_gmm(pool, P["atlas_ks"], P["gmm_niters"], reg, P["seed"], n_jobs=ncores)
    K = atlas.n_components
    ic.to_csv(os.path.join(out, "atlas_model_selection.csv"), index=False)
    status["steps"].update(atlas_K=int(K), atlas_pool=int(pool.shape[0]),
                           atlas_K_at_upper_bound=bool(K == P["atlas_ks"][1] - 1))
    log("Atlas K =", K, "| pool", pool.shape)

    # clusterboot stability (K re-selected in every bootstrap)
    atlas_lab = atlas.predict(pool)
    boots = Parallel(n_jobs=ncores, backend="loky")(
        delayed(boot_once)(pool, atlas_lab, K, P["atlas_ks"], 1, reg, P["seed"] + 1000 + b, P["atlas_boot_frac"])
        for b in range(P["atlas_nboot"]))
    boot_K = np.array([bk for bk, _ in boots])
    J = np.array([jac for _, jac in boots], dtype=float)          # boots × K
    stab = np.nanmean(J, axis=0)
    pd.DataFrame([dict(boot=b, boot_K=boot_K[b], comp=j, jaccard=J[b, j])
                  for b in range(len(boots)) for j in range(K)]) \
        .to_csv(os.path.join(out, "atlas_bootstrap_stability.csv"), index=False)
    status["steps"].update(boot_K_agreement=float((boot_K == K).mean()),
                           boot_K_values=sorted(map(int, set(boot_K))))
    log("Bootstrap K agreement:", status["steps"]["boot_K_agreement"], "| stability:", np.round(stab, 2))

    # ── 3. assign every cell to atlas components ────────────────────────────
    C_all = PA.cat_data(pop, "C")
    resp = atlas.predict_proba(C_all)
    hard = resp.argmax(axis=1)
    cells = meta[["cell_barcode", "sample_id", "microbiome"]].copy()
    for col in ["class_label", "subclass_label", "supertype_label", "cluster_label", "cell_type"]:
        if col in meta.columns:
            cells[col] = meta[col].values
    cells["comp"] = hard
    cells["max_prob"] = resp.max(axis=1)
    cells["loglik_atlas"] = atlas.score_samples(C_all)
    cells["log10_umi"] = np.log10(raw_umi)
    coef = pd.DataFrame(C_all, columns=[f"F{i+1}" for i in range(m)])
    pd.concat([cells, coef], axis=1).to_csv(os.path.join(out, "cells.csv.gz"), index=False)

    # per-animal fractions + per-animal component centroids in feature space
    frac_rows, cen_rows = [], []
    for s, d in cells.groupby("sample_id"):
        ix = d.index.values
        R = resp[ix]; Cs = C_all[ix]
        for j in range(K):
            w = R[:, j]
            frac_rows.append(dict(sample_id=s, microbiome=d.microbiome.iloc[0], comp=j,
                                  n_cells=len(ix), n_hard=int((hard[ix] == j).sum()),
                                  soft_frac=w.mean(), hard_frac=(hard[ix] == j).mean()))
            if w.sum() > 1e-8:
                mu = (w[:, None] * Cs).sum(0) / w.sum()
                cen_rows.append(dict(sample_id=s, microbiome=d.microbiome.iloc[0], comp=j,
                                     eff_n=w.sum(), **{f"F{i+1}": mu[i] for i in range(m)}))
    pd.DataFrame(frac_rows).to_csv(os.path.join(out, "animal_component_fractions.csv"), index=False)
    pd.DataFrame(cen_rows).to_csv(os.path.join(out, "animal_component_centroids.csv"), index=False)

    # popAlign upper-cases symbols on load (load_genes); restore mouse case from
    # genes.tsv (same order) so GO / D-SPIN / MAST joins work downstream.
    orig = pd.read_csv(os.path.join(inp, "genes.tsv"), sep="\t", header=None)[0].astype(str).values
    to_orig = dict(zip(np.char.upper(orig.astype(str)), orig))
    genes_all = orig if len(orig) == len(pop["genes"]) else np.array([to_orig.get(g, g) for g in pop["genes"]])

    # atlas component table (feature-space means + gene-space projection)
    W = np.asarray(pop["W"])
    genes_f = np.array([to_orig.get(g, g) for g in pop["filtered_genes"]])
    comp_rows = []
    for j in range(K):
        gproj = W @ atlas.means_[j]
        top = genes_f[np.argsort(-gproj)[:25]]
        comp_rows.append(dict(comp=j, atlas_weight=atlas.weights_[j], stability=stab[j],
                              frac_boots_jaccard_ge_075=float(np.nanmean(J[:, j] >= P["stable_jaccard"])),
                              n_cells_hard=int((hard == j).sum()),
                              top_projected_genes=";".join(top),
                              **{f"mu_F{i+1}": atlas.means_[j][i] for i in range(m)}))
    pd.DataFrame(comp_rows).to_csv(os.path.join(out, "atlas_components.csv"), index=False)

    # oNMF gene programs
    pd.DataFrame(W, index=genes_f, columns=[f"F{i+1}" for i in range(m)]) \
        .rename_axis("gene").to_csv(os.path.join(out, "onmf_W.csv.gz"))
    pd.DataFrame([dict(feature=f"F{i+1}", rank=r + 1, gene=genes_f[g], weight=W[g, i])
                  for i in range(m) for r, g in enumerate(np.argsort(-W[:, i])[:50])]) \
        .to_csv(os.path.join(out, "onmf_top_genes.csv"), index=False)

    # component markers: specificity-ranked (comp vs rest of unit; all genes).
    # log2FC on the normalised linear scale, expressing-fraction difference, then AUC.
    mk = []
    if K > 1:
        M = PA.cat_data(pop, "M").tocsr()             # log(beta*g/sum + 1)
        Mlin = M.copy(); Mlin.data = np.expm1(Mlin.data)
        mean_in, pct_in, mean_out, pct_out, counts = group_stats(Mlin, hard, K)
        pseudo = 1.0 / beta * 100                       # ~0.01 normalised counts at median depth
        for j in range(K):
            if counts[j] == 0 or counts[j] == len(hard):
                continue
            lfc = np.log2((mean_in[:, j] + pseudo) / (mean_out[:, j] + pseudo))
            dpct = pct_in[:, j] - pct_out[:, j]
            cand = np.where((pct_in[:, j] >= P["marker_min_pct"]) & (dpct >= P["marker_min_dpct"]) &
                            (lfc >= P["marker_min_log2fc"]))[0]
            cand = cand[np.argsort(-lfc[cand])][:3 * P["n_markers"]]
            if not len(cand):
                continue
            auc = auc_rows(M, cand, hard == j)
            for g, a_ in sorted(zip(cand, auc), key=lambda t: -t[1])[:P["n_markers"]]:
                if a_ < P["marker_min_auc"]:
                    break
                mk.append(dict(comp=j, gene=genes_all[g], auc=a_, log2fc=lfc[g],
                               pct_in=pct_in[g, j], pct_out=pct_out[g, j],
                               mean_in=mean_in[g, j], mean_out=mean_out[g, j]))
        del M, Mlin
    pd.DataFrame(mk, columns=["comp", "gene", "auc", "log2fc", "pct_in", "pct_out", "mean_in", "mean_out"]) \
        .to_csv(os.path.join(out, "atlas_markers.csv"), index=False)

    # ── 4. native popAlign: align animals to atlas, deltas, ranking ─────────
    pop["samples"]["ATLAS"] = {"gmm": atlas, "gmm_types": [str(j) for j in range(K)],
                               "gmm_means": atlas.means_, "C": pool,
                               "replicates": {0: {"gmm": atlas, "gmm_means": atlas.means_,
                                                  "gmm_types": [str(j) for j in range(K)]}}}
    PA.align(pop, ref="ATLAS", method="conservative")
    al = []
    for x in pop["order"]:
        for rep in range(pop["nreplicates"]):
            g = pop["samples"][x]["replicates"][rep]["gmm"]
            arr = pop["samples"][x]["replicates"][rep]["alignments"]
            for r_ in np.atleast_2d(arr):
                ti, ri = int(r_[0]), int(r_[1])
                al.append(dict(sample_id=x, rep=rep, test_comp=ti, atlas_comp=ri, jd=r_[2],
                               test_weight=g.weights_[ti], atlas_weight=atlas.weights_[ri],
                               test_K=g.n_components))
    pd.DataFrame(al).to_csv(os.path.join(out, "alignments.csv"), index=False)

    # rank() silently drops animals with <= mincells and then fails on their
    # missing control means, so keep every animal in.
    n_min = min(PA.get_coeff(pop, x).shape[0] for x in pop["order"])
    rank_k = int(min(200, max(n_min - 5, 10)))
    for step, fn in [("plot_deltas", lambda: PA.plot_deltas(pop, sortby="w", pthresh=0.05)),
                     ("rank", lambda: PA.rank(pop, ref="ATLAS", k=rank_k, niter=100, method="LLR",
                                              mincells=n_min - 1))]:
        try:
            fn()
            status["steps"][step] = "ok"
        except Exception as e:  # native helpers are secondary; never fail the unit on them
            status["steps"][step] = f"failed: {e!r}"
            log(step, "failed:", e)
    if pop.get("deltas"):
        pd.concat([pop["deltas"][k]["combined"].T.assign(atlas_comp=k) for k in pop["deltas"]]) \
            .to_csv(os.path.join(out, "popalign_deltas.csv"), index=False)
    if isinstance(pop.get("rankings"), pd.DataFrame):
        pop["rankings"].rename_axis("sample_id").to_csv(os.path.join(out, "popalign_rankings.csv"))

    if not a.no_pickle:
        for x in pop["samples"]:
            for k_ in ("M", "M_norm"):
                pop["samples"][x].pop(k_, None)
        with open(os.path.join(out, "popalign_obj.pkl"), "wb") as fh:
            pickle.dump(pop, fh)

    status.update(status="ok", runtime_min=round((time.time() - t0) / 60, 2))
    json.dump(status, open(status_f, "w"), indent=2, default=str)
    log("Done in", status["runtime_min"], "min")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
