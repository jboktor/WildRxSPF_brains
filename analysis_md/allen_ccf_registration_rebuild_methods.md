# Allen CCF registration rebuild — methods and outputs

## How the latest workflow is executed

From `notebooks/python_scripts/`, using the `spatialomics` conda environment:

```bash
cd /resnick/groups/mthomson/jboktor/WILDRxSPF_brains/notebooks/python_scripts

/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
  -m allen_ccf_repro.cli --all --jobs 5 --promote
```

Useful variants:

- `--samples roi1_run2 roi3_run3` — run a subset of slices
- `--skip-fits --promote` — re-render and promote figures from cached fits (does **not** recompute Elastix)
- `--no-render` — write decisions and the provenance manifest only

Package: `notebooks/python_scripts/allen_ccf_repro/`  
Frozen region map: `notebooks/python_scripts/allen_ccf_repro/configs/region_registry.yaml` (version **4**)

The current `run_manifest.json` may be rewritten by a promote-only pass (`fits_computed_this_run: false`). The command to recompute fits from scratch is:

```bash
python -m allen_ccf_repro.cli --all --jobs 5 --promote
```

---

## Methods (computational)

### Overview

Thirteen hemibrain **seqFISH** spatial genomics sections (`roi*_run*`) were registered to the Allen Mouse Brain Common Coordinate Framework (CCF) at 25 µm resolution using a single automated pipeline (`allen_ccf_repro`). For each section the pipeline (i) chooses atlas plane geometry (anterior–posterior level, optional AP obliquity, medial crop margin), (ii) builds tissue/artifact masks and multi-region corresponding landmarks from cell subclass labels, (iii) runs a candidate ladder of Elastix similarity / affine / B-spline fits, (iv) accepts or refuses each candidate with quantitative gates against the historical baseline registration, and (v) writes final CCF coordinates and quality-control figures. Rejected slices retain the historical fit; status fields and figure titles state that the baseline was retained.

### Inputs

- Slice parameters / Elastix workbook:  
  `data/interim/registration/Allen_CCF_optimization/slice_positions_25um_final.xlsx`
- Historical cell CCF table (baseline):  
  `data/interim/registration/Allen_CCF_optimization/ccf3d.csv`
- Atlas volumes under `data/input/allen_registration_ref/`:  
  `annotation_25.nrrd`, `ara_nissl_25.nrrd`, `annotations_25_border_only.tif`, `annotation_10.nrrd`, `allen_name_to_annots.pkl`
- Region registry v4 mapping `subclass_label_transfer` values to Allen annotation IDs, with safety roles (`hard` / `soft` / `mask_only`), landmark strategies, and acceptance-gate flags

### Per-slice computational sequence

1. **Fixed image and cells**  
   DAPI and cell centroids are prepared in padded fixed-image space (`fixed_x`, `fixed_y`) using the workbook’s reflection, scale, and pad settings. Per-slice Elastix parameters (`spline_grid_size`, iterations, histogram bins) are taken from the workbook when present; otherwise defaults are grid spacing 13 voxels, 2000 B-spline iterations, and bending-energy weight 0.

2. **Baseline reference**  
   Historical CCF coordinates from `ccf3d.csv` are scored for tissue containment and per-region containment/distance against a coronal atlas plane at the workbook’s `allen_slice_num` with medial margin `MIN_MEDIAL_EXTRA = 8` voxels past midline column 228.

3. **AP level, obliquity, and medial-margin search**  
   One stable affine transform is computed at the baseline AP level. Candidate AP levels (±8 atlas slices) and optional mediolateral-dependent AP gradients are scored by re-indexing the annotation through that fixed transform (no re-registration noise between levels). Medial crop margins `{0, 4, 8, 14, 22}` each require their own registration. Usable margins must keep the boundary clamp ratio ≤ 3.0; among near-tied scores, margins that leave room past the dentate gyrus tip are preferred.

4. **Geometry ladder**  
   Candidates are tried in order: searched geometry → baseline AP with the chosen medial margin → historical medial margin if different. Selection prioritizes (1) absence of structural gate failures, (2) fewer regional failures, then (3) mean gated-region containment gain.

5. **Masks**  
   Tissue is the union of the DAPI silhouette and cell-supported territory. Bright cell-free blocks are marked as artifacts. Missing-tissue exclusion is estimated from an undistorted similarity transform (rigid plus uniform scale). Elastix runs unmasked unless dead zones exist; when they do, the fixed mask is the full field minus those zones (no moving mask).

6. **Landmarks** (registry-driven)  
   Density-filtered cell clouds are paired to atlas structures:
   - Dentate gyrus granule cells (`037 DG Glut` ↔ Allen ID 632): mediolateral quantile **ladder** of upper/lower blade edges; when the geometry-only fit is sound, dedicated **medial tip/crest anchors** are added
   - CA fields: mediolateral **ribbon** medians
   - Caudoputamen: dorsal/ventral **lattice** edge stations (up to 8 points)
   - Medial and lateral habenula: robust centroids (protected)
   - Ventricles and choroid plexus: component centroids  
   One eligible region (preferentially CA3) is held out of fitting for validation.

7. **Candidate registration ladder**  
   - `geometry_only` (landmark weight 0)
   - Landmark-weighted fits at weights `{0.03, 0.015, 0.01, 0.05}` via Elastix corresponding points  
   If a fit folds (spatial Jacobian ≤ 0), it is retried along a stiffening ladder of (bending weight, control-grid scale): `(0.15, 1.5)`, `(0.25, 2.0)`, `(0.35, 2.5)`. Default bending weight is otherwise 0.

8. **Acceptance gates** (relative to the historical baseline)  
   Structural failures include tissue-containment loss > 1 percentage point, median displacement of already-correct anchor cells > 8 voxels (when ≥ 200 anchors exist), non-positive Jacobian, medial clamp ratio > 3, and atlas compression into excluded tissue. Regional failures include gated or protected region losses only when median cell-to-region distance also worsens by more than 1 voxel (or gross loss exceeds 20 pp), holdout CA3 relative worsening > 20%, and net mean gated gain < −1 pp or protected gain < −5 pp. Among accepted candidates within 3 pp of the best gated gain, selection prefers lower dentate gyrus p90 distance, then higher DG and caudoputamen containment, then protected-region gain.

9. **Outputs and refusal policy**  
   If accepted, the candidate’s cell coordinates and annotation plane are written. If refused, historical baseline coordinates are written (joined with fixed-image coordinates for QC rendering) and figures use baseline geometry; status is `rejected by gates; baseline retained`.

10. **Rendering and promotion**  
    Figures are written to `figures/Allen_CCF_alignment_optimized/final_rebuild/`, then promoted into `figures/Allen_CCF_alignment_optimized/final/` (stale per-sample files cleared first). `FINAL_FITS_LIBRARY.jpg` is a three-column montage: registration diagnostics | cells + Allen borders | 3D CCF surface.

11. **Consolidated coordinate tables**  
    After all slices finish, per-slice `work/<sample_id>/ccf2d_rebuild.csv` files are concatenated into project-level `Allen_CCF_rebuild/ccf2d.csv` and `Allen_CCF_rebuild/ccf3d.csv` (886,883 cells across the 13 registered sections).

### Coordinate convention

Final CCF columns are micrometers in Allen space: `ccfx`, `ccfy`, `ccfz` (25 µm atlas). `cells_final.csv` also retains `fixed_x` / `fixed_y` and `atlas_x` / `atlas_y` for quality control.

---

## Exact output paths (latest run)

Absolute project root:

`/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/`

### Run-level provenance

| Artifact | Path |
|---|---|
| Decisions table | `data/interim/registration/Allen_CCF_rebuild/rebuild_decisions.csv` |
| Provenance manifest | `data/interim/registration/Allen_CCF_rebuild/run_manifest.json` |
| Per-slice work root | `data/interim/registration/Allen_CCF_rebuild/work/<sample_id>/` |

### Final usable coordinate CSVs (project-level)

These are the consolidated tables to use for downstream analysis. They are assembled automatically at the end of every CLI run from the thirteen per-slice rebuild outputs (886,883 cells):

| File | Absolute path | Contents |
|---|---|---|
| **`ccf2d.csv`** | `/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/data/interim/registration/Allen_CCF_rebuild/ccf2d.csv` | Primary CCF columns only (`ccfx`, `ccfy`, `ccfz`, `annotation`) |
| **`ccf3d.csv`** | `/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/data/interim/registration/Allen_CCF_rebuild/ccf3d.csv` | Same cells plus synced `_2` duplicates (`ccfx_2`, `ccfy_2`, `ccfz_2`, `annotation_2`) |

Relative to the project root:

- `data/interim/registration/Allen_CCF_rebuild/ccf2d.csv`
- `data/interim/registration/Allen_CCF_rebuild/ccf3d.csv`

Do **not** use the older copies under `data/interim/registration/ccf2d.csv`, `data/interim/registration/ccf3d.csv`, or `data/interim/registration/Allen_CCF_optimization/` for rebuild results; those are historical baselines.

### Per-slice coordinate CSVs

For each of the 13 `sample_id`s (`roi1_run1`, `roi1_run2`, `roi2_run1`, `roi2_run2`, `roi2_run4`, `roi3_run1`, `roi3_run2`, `roi3_run3`, `roi3_run4`, `roi4_run1`, `roi4_run2`, `roi4_run3`, `roi4_run4`):

| File | Path | Contents |
|---|---|---|
| Per-slice rebuild table | `data/interim/registration/Allen_CCF_rebuild/work/<sample_id>/ccf2d_rebuild.csv` | Full cell table with updated `ccfx`, `ccfy`, `ccfz` and `annotation` (inputs to the consolidated CSVs) |
| Slim registration table | `data/interim/registration/Allen_CCF_rebuild/work/<sample_id>/cells_final.csv` | `fixed_x/y`, `atlas_x/y`, `ccfx/y/z`, `annotation` |
| Decision detail | `data/interim/registration/Allen_CCF_rebuild/work/<sample_id>/detail.json` | AP search, landmarks, candidate ladder, gates |
| Render cache | `data/interim/registration/Allen_CCF_rebuild/work/<sample_id>/render_inputs.npz` | Fixed/moving images, masks, and annotation plane used for figures |

Example (accepted slice):

`/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/data/interim/registration/Allen_CCF_rebuild/work/roi1_run2/ccf2d_rebuild.csv`

### Figures (promoted finals)

Root:

`figures/Allen_CCF_alignment_optimized/final/`

| Artifact | Path |
|---|---|
| Library montage (3 columns) | `figures/Allen_CCF_alignment_optimized/final/FINAL_FITS_LIBRARY.jpg` |
| Figure manifest | `figures/Allen_CCF_alignment_optimized/final/FINAL_FIGURE_MANIFEST.csv` |
| README | `figures/Allen_CCF_alignment_optimized/final/README.md` |
| Per-slice QC | `figures/Allen_CCF_alignment_optimized/final/<sample_id> 2d reg.png` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id> 3d reg.png` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id> regional targets.png` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id> mask QC.png` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id> landmark QC.png` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id>_slice_###_allen_###_reg_5.jpg` |
| | `figures/Allen_CCF_alignment_optimized/final/<sample_id>_slice_###_allen_###_reg_5_cells_all.jpg` |

Staging copies (same content before/alongside promote):

`figures/Allen_CCF_alignment_optimized/final_rebuild/`

### Latest acceptance state

From `rebuild_decisions.csv`:

- **Accepted (11):** `roi1_run1`, `roi1_run2`, `roi2_run1`, `roi2_run2`, `roi3_run1`, `roi3_run2`, `roi3_run3`, `roi3_run4`, `roi4_run1`, `roi4_run3`, `roi4_run4`
- **Baseline retained (2):** `roi2_run4`, `roi4_run2`

### Related notebook narrative

`notebooks/0X_Allen_CCF_Registration_optimized.ipynb` (rebuild section near the end).
