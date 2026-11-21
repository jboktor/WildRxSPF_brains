# Input inventory — olfactory–amygdalar chemosensory relay (VR-1…VR-4)

Generated for the four relay notebooks. Paths are absolute under
`/resnick/groups/mthomson/jboktor/WILDRxSPF_brains`. Unit of replication is the
**animal** (n = 6 WildR, n = 6 SPF). Do not re-run MAST, D-SPIN, or pySCENIC GRN.

Terminology in all outputs: **"the olfactory–amygdalar chemosensory relay,
including its vomeronasal-recipient nodes"** — not "the vomeronasal relay."

---

## Status summary

| Artifact | Status | Used by |
|---|---|---|
| MAST hit catalogue (padj < 0.1, constructed on load) | **Found** (source RDS; curated 243×27 CSV not on disk) | VR-1, VR-3, VR-4 |
| `hit_triage_table.csv` | **Missing** — notebooks flag and proceed without triage annotations | VR-1, VR-3 |
| MAST full results (all genes / unit) | **Found** (`supertype_label/*_2026-01-04.rds`) | VR-3, VR-4 |
| D-SPIN `star_programs.csv` | **Found** (245 rows) | VR-1, VR-2 |
| D-SPIN per-animal `h` | **Found** (`responses.csv` + `relative_responses.csv` per unit) | VR-1 |
| D-SPIN `programs_genes.csv` | **Found** (per unit folder) | VR-2 |
| pySCENIC regulons pre-filter | **Found** (`regulons.csv` per tissue×class) | VR-4 |
| pySCENIC AUCell matrix | **Found** (`auc_mtx.csv` per tissue×class) | VR-4 |
| snRNA counts + metadata | **Found** (Seurat + per-unit h5ad) | VR-1…4 |
| seqFISH / spatial object | **Found** (prenorm Seurat + MapMyCells + CCF rebuild) | VR-2, VR-3 |

**VR-1 can run as specified** — per-animal `h` exists in each unit’s
`responses.csv` / `relative_responses.csv` (keyed by `sample_id`).

---

## Artifact details

### MAST hit catalogue

| Field | Value |
|---|---|
| Path | `data/interim/mast/snRNASeq/supertype_mast_signif_hits_2026-01-04.rds` |
| Shape | 148,286 × 10 (all tested gene×supertype×tissue rows; filter `p_adjust < 0.1` to get hits) |
| Columns | `celltype_tissue`, `primerid`, `Pr(>Chisq)`, `coef`, `ci.hi`, `ci.lo`, `p_adjust`, `tissue`, `cell_type`, `supertype_name` |
| Gene IDs | `primerid` with `-GRCm39` suffix (strip before matching symbols) |
| Notes | No standalone `hit_catalogue_243.csv`. Notebooks build the relay catalogue from AMY hits in chain positions 1–4 after §0.4 exclusions. Expected core: `Eps8l1`, `Zbtb16`, `Hdac9`, `Trpm3`, `Rasa1`; signed-negative: `Itpkb`. |

Related aggregates: `class_*` / `subclass_*` RDS in the same folder. Optional CSV export
`mast_signif_hits_2026-08-05_p_adjust_0.1.csv` may be written by `0X_DEG_Analysis.qmd`.

### `hit_triage_table.csv`

**Not found** anywhere in the repo. Deviation: VR notebooks record absence and do
not invent triage labels.

### MAST full results

| Field | Value |
|---|---|
| Path pattern | `data/interim/mast/snRNASeq/supertype_label/{supertype_label}_{tissue}_2026-01-04.rds` |
| Example | `…/CS20230722_SUPT_0055_AMY_2026-01-04.rds` → 9,164 × 6 |
| Columns | `primerid`, `Pr(>Chisq)`, `coef`, `ci.hi`, `ci.lo`, `p_adjust` |
| Satisfies | VR-3 gradient regression; VR-4 MAST context |

### D-SPIN star programs

| Field | Value |
|---|---|
| Path | `tables/DSPIN/star_programs.csv` |
| Shape | 245 × 8 |
| Columns | `tissue`, `supertype_name`, `supertype_label`, `program_id`, `n_wildr_samples`, `mean_h`, `sd_h`, `same_sign` |
| Key rows for GABA limb | `0338` P10, `0355` P5, `0359` P6 (all AMY, negative mean_h) |

Annotated twin: `tables/DSPIN/star_programs_david.csv`.

### D-SPIN per-animal `h` (required for VR-1)

| Field | Value |
|---|---|
| Unit dirs | `data/interim/dspin/{tissue}__{supertype_label}/` (152 units) |
| Absolute `h` | `responses.csv` — cols `program_id`, `sample_id`, `h` |
| Relative `h` | `relative_responses.csv` — cols `program_id`, `sample_id`, `relative_h` |
| Sample meta | `sample_meta.csv` — `sample_id`, `microbiome`, `batch`, `if_control`, `n_cells` |
| Cell meta | `cells_meta.csv` — animal, tissue, class/subclass/supertype |

VR-1 uses **absolute** `h` from `responses.csv` for the three GABA programs, then
z-scores across the 12 animals. `relative_h` is retained as a sensitivity check.

Relay glut units with h5ad present: `0039–0042`, `0055–0058`, `0061`, `0507`.
**Missing dspin units:** `0059`, `0060`, `0076` (0076 is excluded from relay
interpretation anyway as hippocampal carry-over).

### D-SPIN program genes

| Field | Value |
|---|---|
| Path | `data/interim/dspin/{unit}/programs_genes.csv` |
| Columns | `program_id`, `gene`, `weight` (weight often blank) |
| Satisfies | VR-2 Jaccard / enrichment / MAST bridge |

### pySCENIC regulons (pre-significance-filter)

| Field | Value |
|---|---|
| Path | `data/interim/scenic/by_class/snRNASeq_{tissue}_{class_label}/regulons.csv` |
| Format | pySCENIC ctx multi-header CSV (TF, MotifID, AUC, NES, …) |
| Relevant AMY classes | `CLAS_09` (09 CNU-LGE GABA), `CLAS_11` (11 CNU-HYa GABA), `CLAS_13` (13 CNU-HYa Glut) |
| Zbtb16 check (inventory) | **Not called as TF** in CLAS_09/11/13; not present in `auc_mtx` columns → VR-4 case **(a)** pending notebook confirmation |

Diff tables (post-test): `diff_regulons_WildR_vs_SPF.csv` per split;
aggregates under `data/processed/scenic_by_class/`.

### pySCENIC AUCell

| Field | Value |
|---|---|
| Path | `…/by_class/…/auc_mtx.csv` |
| Shape example | CLAS_09 ≈ 13,262 cells × ~55–80 regulons |
| Columns | `Cell`, then `TF(+)` columns |
| Cell↔animal join | `data/interim/scenic/cell_meta.csv` (`cell_id`, `sample`, `microbiome`, `tissue`, `supertype_name`, …) |

**Do not pool raw AUCell across classes** — z-score within class first (spec §VR-4 failure modes).

### snRNA count matrix + metadata

| Field | Value |
|---|---|
| Primary Seurat | `data/interim/seurat/snRNASeq/seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds` (~39 GB) |
| Filter | `keep_cell == TRUE`; droplist `…/celltype_droplist_2025-12-30.rds` |
| Lightweight path | Per-unit `raw_counts.h5ad` + `cells_meta.csv` under `data/interim/dspin/` |
| Metadata fields | `sample` / `sample_id`, `microbiome`, `tissue`, `class_*`, `subclass_*`, `supertype_*` |

### seqFISH / spatial

| Field | Value |
|---|---|
| Prenorm Seurat | `data/interim/seurat/seqFISH/merged_roi_seurat_prenorm_2026-02-16.rds` (~4.4 GB) |
| Legacy meta | `data/interim/merged_roi_seurat_metadata_2024-07-15.rds` — coords `center_x/y`, `roi`, `bregma`, `group` (no Allen supertype) |
| MapMyCells (Allen hierarchy) | `data/interim/MapMyCells/output.csv` (~719k rows) — `class/subclass/supertype/cluster` labels + bootstrap probs |
| CCF rebuild | `data/interim/registration/Allen_CCF_rebuild/` (`ccf2d.csv`, `ccf3d.csv`, `rebuild_decisions.csv`) |
| Notes | Spatial cell-type reference may differ from snRNA MapMyCells; VR-2/VR-3 must report which reference was used. MeApd/MeApv partition must be **within section**, not a global threshold. |

---

## Relay chain map (§0.8)

| Position | Node | Supertypes (match by name prefix / label) |
|---|---|---|
| 1 | Cortical input | PIR-ENTl L2/3 IT (`0039`–`0042` and siblings) |
| 2 | Cortical amygdala | `0059 COAp Grxcr2 Glut 2`, other COAp (`0058`, …) |
| 3 | MeA glutamatergic | `0055`, `0056`, `0057` MEA Slc17a7 Glut; `0507 MEA-COA-BMA` |
| 3G | MeA GABAergic output (parallel limb) | `0338`, `0355`, `0359` MEA-BST |
| 4 | Basolateral complex | `0060`, `0061` LA-BLA-BMA-PA |

Never put 3G into ordinal regressions over 1–4. Exclude `0076 CA3 Glut 2` from relay interpretation.

## Permanent gene exclusions (§0.4)

`ENSMUSG00000095041`, `Cmss1`, `Lars2`, `Mid1` — never inside composite scores.

---

## Notebook outputs (expected)

| Experiment | Notebook | Output dir |
|---|---|---|
| VR-1 | `notebooks/0X_VR1_Input_Output_Coupling.qmd` | `tables/VR/VR1/`, `figures/VR/VR1/` |
| VR-2 | `notebooks/0X_VR2_MeA_BST_Program_Identity.qmd` | `tables/VR/VR2/`, `figures/VR/VR2/` |
| VR-3 | `notebooks/0X_VR3_Relay_Gradient_MeApd.qmd` | `tables/VR/VR3/`, `figures/VR/VR3/` |
| VR-4 | `notebooks/0X_VR4_Relay_Regulons.qmd` | `tables/VR/VR4/`, `figures/VR/VR4/` |

Shared helpers: `notebooks/R_scripts/vr_relay_helpers.R`.
