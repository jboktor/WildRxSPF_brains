# Hemisphere-aware 2DG × seqFISH analysis plan

**Status:** steps A–F implemented (2026-11-14e); G–H still open  

**Scope:** Allen CCF co-embedding of group-level 2DG maps with hemibrain seqFISH NEUROeSTIMator sections  
**Unit of replication (animal):** unchanged — still unavailable for seqFISH until a verified section→animal key exists; hemisphere policy must not invent animal-level inference.

---

## 1. Problem

SeqFISH sections are **hemibrain** acquisitions. During Allen CCF registration, several sections were **left–right reflected** (`reflect=True` in `slice_positions_25um_final.xlsx`) so that tissue could be fit against a consistent medial atlas crop. That is convenient for registration QC, but it:

1. Places originally contralateral tissue onto the **same** CCF hemisphere in downstream 3D views.
2. Makes left/right biological structure look like a single-sided cloud.
3. Invites analyses that silently **collapse** both hemispheres or treat mirrored coordinates as anatomical truth.

Separately, the 2DG Template→CCF pipeline applies a **global** horizontal flip (`FLIP_SOURCE_LEFT_RIGHT = TRUE`) because the source images’ left side is the brain’s right (email convention). That flip is an orientation correction for whole coronal 2DG sections that already span both hemispheres — it is **not** the same operation as seqFISH per-section hemibrain reflection.

**Goal:** keep full-brain 2DG and the true anatomical side of every seqFISH section visible; annotate hemispheres explicitly; classify bilateral patterns in 2DG and NEUROeSTIMator as symmetric, unilateral, or opposing — without forcing every section onto one side.

---

## 2. Current-state audit (facts)

### 2.1 SeqFISH registration reflections

From `data/interim/registration/Allen_CCF_optimization/slice_positions_25um_final.xlsx` (`slice_parameters.reflect`):

| `reflect` | Sections |
|-----------|----------|
| **TRUE** (mirrored into atlas medial half) | `roi2_run1`, `roi2_run2`, `roi2_run4`, `roi4_run1`, `roi4_run2`, `roi4_run3`, `roi4_run4` |
| **FALSE** (kept as imaged) | `roi1_run1`, `roi1_run2`, `roi1_run4`, `roi3_run1`, `roi3_run2`, `roi3_run3`, `roi3_run4` |

Registered cells in the 2026-11-11 join sit almost entirely at **low ML** (`ccfz` ≈ 0.5–6.0 mm; atlas midline ≈ **5.7 mm**). That is expected after medial-half registration, including for `reflect=True` sections.

Allen CCFv3 annotation IDs are largely L/R-shared (~95% id match under mirror at a test plane), so reflecting `ccfz` about midline for **display** usually preserves structure names. Still treat contralateral ID remap as a QC check, not a free lunch.

### 2.2 2DG orientation

- Global source flip before Elastix: image-left → anatomical-right, matching CCF (`ccfz = 0` = anatomical left).
- After warp, Intruder significant voxels occupy **both** hemispheres (~47% ML < midline, ~53% ML > midline across slices).
- Per-slice median ML still wanders leftish/rightish with AP — anatomical, not a hemibrain reflect flag.

### 2.3 What must not be done

- Do **not** re-fit Elastix by forcing all seqFISH onto one side as the long-term analytical coordinate system.
- Do **not** average left and right 2DG parcels without a symmetry class.
- Do **not** treat mirrored display coordinates as a new registration truth without provenance columns.
- Do **not** claim animal-level hemisphere effects from seqFISH without a section→animal key.

---

## 3. Principles

1. **Full brain first for 2DG.** Keep both hemispheres in all map summaries; hemisphere is a factor, not a discard rule.
2. **Original side first for seqFISH display.** Registration may use a reflected moving image; analytical/display coordinates should carry an **un-reflected anatomical side** plus the registration-side copy.
3. **Annotate, don’t collapse.** Every cell/voxel row gets `hemisphere`, `reflect_applied`, and `side_for_analysis`.
4. **Symmetry is a result, not an assumption.** Bilateral patterns are labeled *symmetric / unilateral / opposing / indeterminate* with explicit tests.
5. **Provenance is cheap; silent flips are expensive.** Workbook `reflect`, 2DG `FLIP_SOURCE_LEFT_RIGHT`, midline definition, and display transforms are written to dated tables.

---

## 4. Coordinate and labeling policy

### 4.1 Midline

Use Allen 25 µm atlas midline column **228** → `ML_MIDLINE_UM = 5700`.  
Mirror transform for display/analysis un-reflect:

\[
\mathrm{ccfz\_anatomical} =
\begin{cases}
2\cdot 5700 - \mathrm{ccfz\_registered} & \text{if } \texttt{reflect}=\mathrm{TRUE}\\
\mathrm{ccfz\_registered} & \text{otherwise}
\end{cases}
\]

Keep `ccfz` (registered) unchanged in primary join tables; add `ccfz_anatomical` (and optionally `ccfy` unchanged).

### 4.2 Hemisphere labels

| Label | Rule |
|-------|------|
| `left` | `ccfz_anatomical < 5700` |
| `right` | `ccfz_anatomical > 5700` |
| `midline_band` | `|ccfz_anatomical - 5700| ≤ 250` µm (optional QC band; exclude from strict L/R contrasts) |

Also store:

- `reflect_applied` (from workbook)
- `registration_hemisphere` (side in registered coords)
- `anatomical_hemisphere` (side after un-reflect)
- `side_source` = `workbook_reflect | geometry_inferred | manual_override`

### 4.3 Dual coordinate tracks

| Track | Purpose |
|-------|---------|
| **Registered CCF** (`ccfx,ccfy,ccfz`) | Reproducible join to current Elastix outputs / QC figures |
| **Anatomical CCF** (`ccfz_anatomical`, same AP/DV) | Visuals, hemisphere annotation, bilateral symmetry tests |

Do not overwrite registered coordinates in `data/interim/registration/Allen_CCF_rebuild/`.

### 4.4 2DG

- Retain global orientation flip (convention correction).
- Label every overlay voxel and region aggregate with `hemisphere` from registered `ccfz` vs midline.
- For region summaries, produce **bilateral**, **left**, **right**, and **L−R difference** columns rather than a single pooled mean when testing symmetry classes.

---

## 5. Symmetry signatures (2DG and NEUROeSTIMator)

For each Allen structure (and optionally each AP slab), compute left and right signed effects:

- 2DG Baseline: \(t_L, t_R\) (WildR − SPF)
- 2DG Intruder: \(t_L, t_R\)
- seqFISH NE: \(\Delta_L, \Delta_R\) (WildR − SPF mean `predicted_activity`, descriptive until animal key exists)

### 5.1 Classification rules (pre-specify; tune thresholds only in a dated config)

Let \(a,b\) be the two hemispheres’ effects for one modality/map.

| Class | Definition (operational) | Biological reading |
|-------|--------------------------|--------------------|
| **Symmetric** | \(\mathrm{sign}(a)=\mathrm{sign}(b)\) and \(\min(\|a\|,\|b\|) \ge \tau\) and \(\|a-b\| / (\|a\|+\|b\|+\varepsilon) \le \delta\) | Same direction, similar magnitude both sides |
| **Unilateral** | One side \(\ge \tau\), other \(\le \tau_0\) (near zero / NS) | Effect largely confined to one hemisphere |
| **Opposing** | \(\mathrm{sign}(a)=-\mathrm{sign}(b)\) and both \(\|a\|,\|b\| \ge \tau\) | Contralateral push–pull |
| **Indeterminate** | Sparse voxels/cells, midline-heavy, or fails QC | Do not interpret |

Suggested starting thresholds (to lock in config before looking at H-tests):  
\(\tau\) = 1.725 for 2DG voxel-level; for region means use permutation or bootstrap CIs.  
For NE descriptive \(\Delta\), use a null from gate or section shuffles — **not** animal-level p-values.

### 5.2 Cross-modality signature concordance

Build a region×map table:

- `twodg_baseline_symmetry_class`
- `twodg_intruder_symmetry_class`
- `ne_symmetry_class`
- `baseline_vs_intruder_class_agree` (logical)
- `intruder_vs_ne_class_agree` (logical)
- `lateralization_index` = \((E_R - E_L) / (|E_R| + |E_L| + \varepsilon)\)

Questions this enables (descriptive until animal key):

1. Are Intruder WildR−SPF effects more lateralized than Baseline?
2. Do NE lateralization indices track Intruder more than Baseline (echo of H3, now with side)?
3. Which limbic/hypothalamic parcels are opposing in 2DG but symmetric in NE (or vice versa)?

### 5.3 What “full brain” means operationally

- **2DG tests:** primary = bilateral pool; secondary = L-only, R-only, and L−R interaction / lateralization index.
- **seqFISH:** never invent the missing hemisphere; analyze available anatomical side, and when stacking sections use `anatomical_hemisphere` as a covariate / facet — not a reason to drop reflected sections.
- **Co-embedding plots:** show anatomical-side seqFISH and full bilateral 2DG; legend encodes hemisphere and reflect provenance.

---

## 6. Visual remediation (near-term)

1. Load workbook `reflect` into a dated provenance CSV.
2. In Plotly 3D / transition / AP-overlap views, plot seqFISH with `ccfz_anatomical`.
3. Hover fields: `anatomical_hemisphere`, `reflect_applied`, registered ML.
4. Optional toggle: **Registered side** vs **Anatomical side** (default anatomical).
5. Z-stack browser: seqFISH is not the moving image there; for 2DG keep as-is (already bilateral). Add a midline guide and L/R tags on colorbars (“left ML ← midline → right ML”).
6. Do **not** silently rewrite `ccf3d.csv`; display-layer only until a formal `ccf3d_anatomical_*.csv` export is reviewed.

---

## 7. Analysis roadmap (ordered)

| Step | Deliverable | Notes |
|------|-------------|-------|
| A | `seqFISH_section_hemisphere_provenance_*.csv` | **Done** (2026-11-14) |
| B | Display un-reflect in interactive viewers | **Done** in `2026-11-14d` Plotly viewers |
| C | `*_with_hemisphere_*.rds/csv` join add-on | **Done** `2026-11-14e` |
| D | Region bilateral summaries for Baseline, Intruder, I−B maps | **Done** `2DG_region_bilateral_summaries_2026-11-14e.csv` |
| E | NE descriptive bilateral / available-side summaries | **Done** `SeqFISH_NE_region_*_2026-11-14e.csv` |
| F | Symmetry-class tables + figure cards | **Done** `2DG_NE_hemisphere_symmetry_classes_2026-11-14e.*` |
| G | Revisit H2–H4 with hemisphere factor | Open — pre-specified; gated `run_analysis` |
| H | Only after section→animal key | Open — mixed models with hemisphere nested in animal |

---

## 8. Risks and open questions

1. **Which physical hemisphere was imaged** for each `reflect=True` section — confirm against lab notes / DAPI montage, not only the flag.
2. Midline band cells and ventricle-crossing structures inflate “opposing” calls — use `midline_band` exclusion.
3. 2DG group maps are not animal-level; lateralization is a **group-map property**, not an individual laterality claim.
4. RGB sibling maps (Baseline / I−B) inherit Template warps; hemisphere labels follow warped `ccfz`.
5. If a section’s `reflect` flag is wrong, anatomical un-reflect will put it on the incorrect side — keep a manual override column.

---

## 9. Immediate decisions (locked for implementation)

1. Midline = **5700 µm**.
2. Default visuals use **anatomical** ML for seqFISH; registered ML retained in data.
3. 2DG global flip stays; 2DG hemisphere from registered ML.
4. Symmetry classes follow §5.1 with dated threshold config before hypothesis re-runs.
5. Full-brain 2DG is the default analysis domain; hemisphere splits are additional, not replacements.

---

## 10. Related paths

- SeqFISH rebuild methods: `analysis_md/allen_ccf_registration_rebuild_methods.md`
- Workbook: `data/interim/registration/Allen_CCF_optimization/slice_positions_25um_final.xlsx`
- 2DG flip: `notebooks/python_scripts/register_2dg_to_allen_ccf.py` (`FLIP_SOURCE_LEFT_RIGHT`)
- Integration notebook: `notebooks/0X_2DG-SeqFISH_Allen-CCF_Integration.qmd`
