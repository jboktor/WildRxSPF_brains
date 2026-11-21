# WILDRxSPF_brains

Spatial (seqFISH) + single-nucleus RNA-seq analysis of WildR vs. SPF mouse brains
(amygdala `AMY`, hypothalamus `HYP`). Unit of replication is the **animal**.

**Canonical project root:** `/resnick/groups/mthomson/jboktor/WILDRxSPF_brains`
(`wkdir`). The old `/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains` is a
*different, near-empty directory*. Some older scripts still hardcode it. Use the
`mthomson` path in all new code, and point out stale `MazmanianLab` project
paths when you come across them. Software and env paths under
`/resnick/groups/MazmanianLab/jboktor/software/` are correct and stay as they are.

---

## 1. Safety rules (highest priority)

1. **Never delete anything without asking first.** This covers files, directories,
   data, figures, logs, `.cluster_runs/` registries and git history. That means no
   `rm`, `rm -rf`, `git clean`, `git reset --hard`, `git checkout -- <file>`,
   `git stash drop`, or overwriting with `>`, unless I have explicitly approved
   that specific action in this conversation. Approving one deletion does not
   approve the next.
2. **Prefer not deleting at all.** Data is expensive to regenerate. If you find
   junk, report it to me in a short list (path, size, why you think it's junk)
   and let me handle it. Examples: stale scratch files, `test2.py`, `scratch.R`,
   `*_OLD` folders, duplicate RDS files, `figures/temp`.
3. **Never overwrite an existing data artifact.** Write a new date-stamped
   version (see §4) instead of clobbering `data/interim/*` or `data/processed/*`.
   If overwriting truly is the right call, ask first.
4. **`data/input/` is read-only.** Never modify, move or rename raw/external inputs.
5. **Be careful with moves and renames** of anything that existing code
   references. Grep for references first, and update them in the same change.
6. **Secrets:** never read, print, or commit `.env`, `.secrets/`, or
   API keys in `.Renviron`.
7. Don't commit or push unless I ask.

When unsure whether an action is destructive, treat it as destructive and ask.

---

## 2. Repository layout — put things where they belong

```
WILDRxSPF_brains/
├── data/
│   ├── input/        # raw + external inputs (READ-ONLY): references, atlases, raw exports
│   ├── interim/      # intermediate objects: per-step Seurat/h5ad, model fits, tool outputs
│   └── processed/    # final, analysis-ready objects that downstream notebooks load
├── notebooks/        # .qmd analysis records (+ legacy .ipynb)
│   ├── R_scripts/       # R scripts & function libraries (*.R only)
│   ├── python_scripts/  # Python scripts (*.py only)
│   └── shell_scripts/   # sbatch / shell wrappers (*.sh)
├── figures/<Analysis>/   # all figure images (png/pdf/svg), one subfolder per analysis
├── tables/<Analysis>/    # all result tables (csv/tsv/xlsx), one subfolder per analysis
├── analysis_md/          # markdown method notes, input inventories, design docs
├── interpretations/      # written interpretation reports of results
├── batchtools_templates/ # SLURM templates for batchtools (per conda env)
├── .cluster_runs/        # SLURM logs + batchtools registries (see §6)
└── quarto-aesthetic-kit/ # house notebook style (see §5)
```

Rules:
- **Inputs vs. intermediates vs. finals:** results of a computation go to
  `data/interim/<analysis>/`. An object only goes to `data/processed/` once it
  is the cleaned version that other analyses load.
- **Figures go to `figures/<Analysis>/`** and **tables to `tables/<Analysis>/`**,
  not into `notebooks/` or `notebooks/output/`. `notebooks/output/` and the
  top-level `output/` are legacy. Don't add new things there.
- **Scripts go in the folder for their language.** A `.py` file belongs in
  `python_scripts/`, even when an R workflow calls it.
- Large data, figures, `.cluster_runs/`, rendered `notebooks/*.html` and
  `notebooks/misc_*` are **gitignored** (see `.gitignore`). Don't try to commit
  them, and remember that git doesn't back them up. This is one more reason
  never to delete.
- Build paths from `wkdir` with `glue("{wkdir}/...")` / `file.path(wkdir, ...)`.
  Never use paths relative to wherever the code happens to be running.

---

## 3. How I work: QMD notebooks are historical records

My `.qmd` files are **rich, descriptive records of how a whole analysis was done**.
Someone should be able to read one months later and see what was run, why,
with what parameters, where the outputs live and what we concluded.
**They are not active pipelines.** Re-rendering a notebook must be cheap and
must never recompute heavy steps or submit jobs.

Pattern (follow the existing `0X_VR*`, `0X_DEG_Analysis-top-hits.qmd` and
`0X_FluxAnalysis.qmd` notebooks):

- A `setup-display` chunk at the top sets `run_analysis <- FALSE`. It loads only
  the libraries and small artifacts that the display needs.
- **Every compute chunk** is gated with `#| eval: !expr run_analysis`. This
  includes loading big objects, fitting models, writing outputs, and any
  `sbatch` / `batchtools::submitJobs` call. Never add an ungated `sbatch` to a `.qmd`.
- **Active code at render time** should only load saved tables and figures from
  disk for display, including interactive tables and plots (DT, plotly, etc.).
- Heavy logic lives in `notebooks/R_scripts/` / `python_scripts/`. The notebook
  `source()`s those files or shows the exact command that ran them, so the
  record is complete without the notebook doing the work.
- When I run something directly (not via the notebook), **still record it in the
  relevant `.qmd`**. Add a gated chunk with the exact command or code, the
  parameters, and where it wrote outputs.
- If an artifact the notebook expects is missing, the display chunk should print
  a clear "not yet computed — set `run_analysis <- TRUE`" message, not error.
- Surface key file paths with `.path` chips so readers know where things live.

---

## 4. Naming conventions

**Notebooks** (`notebooks/`): `<NN>_<Topic-Name>_<Modality>.qmd`
- `NN` = position in the main pipeline (`01_`, `02_`, … `07.2_` for sub-steps).
- `0X_` = analysis not yet placed in the numbered sequence (most current work).
- `misc_` = tutorials, downloads, one-offs (gitignored).
- Modality suffix: `_SeqFISH`, `_SN` / `_snRNASeq` (single-nucleus). Use
  `_SN` for new files, and add the suffix when the modality isn't obvious.
- Related series share a prefix and index, e.g. `0X_VR1_…`, `0X_VR2_…`.

**Data artifacts:** `<descriptive_name>_<YYYY-MM-DD>.<ext>`. The date goes at the
end, and the description says what the object is and its key processing state,
e.g. `merged_roi_seurat_filtered_staNMF-updated_2024-07-24.rds`. (Older
files put the date first. Use date-last for new files.)

**Scripts:** verb-first snake_case describing the action:
`run_<tool>_<scope>.sh` for sbatch wrappers, `export_*`, `aggregate_*`,
`*_functions.R` / `*_helpers.R` for function libraries.

**Figures / tables:** `<Analysis>_<what>_<comparison-or-subset>.png|csv`, inside
`figures/<Analysis>/` or `tables/<Analysis>/`.

**Renaming:** you may rename **new files you created** if the new name is clearly
better for clarity or structure. Ask before renaming files that already exist,
and update every reference when you do.

---

## 5. New QMD notebooks → always use the quarto-aesthetic-kit

Every new `.qmd` is built from `quarto-aesthetic-kit/`. **Read
`quarto-aesthetic-kit/CLAUDE.md` first and follow it exactly.** It covers the
section order (Overview → Background → Results → Next Steps → Code & Provenance),
figure cards, takeaways, provenance labels (`cell-fig-*`, `cell-tbl-*`) and
the quality checklist. Start from `quarto-aesthetic-kit/templates/analysis-template.qmd`.

Repo-specific adaptations:
- There is no root `_quarto.yml`, so use **single-file mode**. Copy the
  front matter from an existing `0X_VR*.qmd`, with theme paths relative to
  `notebooks/`: `../quarto-aesthetic-kit/theme/aesthetic.scss` and
  `../quarto-aesthetic-kit/theme/fonts.html`.
- `author: "Joe Boktor"`. Code chunks are R (`{r}`) unless the analysis is Python.
- Figures and tables are written to `{wkdir}/figures/<Analysis>/` and
  `{wkdir}/tables/<Analysis>/`, not to the kit's local `figures/` folder.
- Combine the kit's results-up/code-down pattern with the `run_analysis` gate (§3).
- Never hand-write CSS or colors. All styling comes from the theme.

---

## 6. Compute: SLURM vs. running directly

**Use SLURM** when a job is heavy enough that it shouldn't run on the login
node. As a rule of thumb, that means it will take more than ~10–15 min, needs
more than ~16 GB RAM or more than ~4 cores, needs a GPU, or fans out over many
units (cell types, tissues, splits). **Otherwise run it directly.** Quick
loads, summaries, small plots and test runs on a subset are fine to run directly.

**How to submit:**
- **Shell/sbatch jobs:** put the script in `notebooks/shell_scripts/run_<tool>_<scope>.sh`.
  Use `set -euo pipefail`, absolute env binaries, a header comment that shows
  how to submit/override, and `--mail-user=jboktor@caltech.edu --mail-type=FAIL,END`.
  Use array jobs (`--array=0-N%K`) driven by a manifest in `data/interim/` for fan-out.
- **R fan-out:** use `batchtools` with a template from `batchtools_templates/`
  (e.g. `batchtools.slurm_mamba-spatialomics.tmpl`). See `0X_MAST_SN.qmd`.
- Partitions in use: `expansion` (CPU) and `gpu` (e.g. `--gres=gpu:nvidia_h200:1`).
  Request resources sensibly, without padding them by 10×.
- Record the submission command in the relevant `.qmd` inside a
  `#| eval: !expr run_analysis` chunk, so rendering never resubmits.

**Logging runs in `.cluster_runs/`** (always use the `mthomson` path):
- **batchtools registries:** `.cluster_runs/<YYYY-MM-DD_HH:MM:SS>_<AnalysisName>/`,
  created with `glue("{get_time()}_<AnalysisName>")` (`get_time()` is in
  `R_scripts/_misc_functions.R`), e.g. `2026-01-04_14:57:54_MAST_SingleNuclei`.
- **sbatch scripts:** log into a per-tool folder `.cluster_runs/<tool>/`, with
  separate stdout/stderr:
  ```
  #SBATCH --output=".../WILDRxSPF_brains/.cluster_runs/<tool>/log_<job>_%A_%a.log"   # arrays
  #SBATCH --error=".../WILDRxSPF_brains/.cluster_runs/<tool>/log_<job>_%A_%a.err"
  #SBATCH --output=".../WILDRxSPF_brains/.cluster_runs/<tool>/log_<job>_%x_%j.log"   # single jobs
  ```
- Print timestamped progress lines in scripts (`log()` helper as in `run_scFEA.sh`).
- After submitting, tell me the job ID(s) and log path. Check logs before
  claiming success.

---

## 7. Environment

Never use bare `R`, `Rscript`, `python` or `quarto` from PATH. Use the
**spatialomics** env binaries:

```bash
/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript script.R
/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python script.py
/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/quarto render notebooks/foo.qmd
```

Tool-specific envs (scFEA, neuroestimator, pySCENIC) have their own
`batchtools_templates/*.tmpl` and sbatch activation lines. Reuse those.
Render `.qmd` files from `notebooks/`, and report the exact command and env used.

---

## 8. Working style

- Read existing notebooks and scripts for the analysis before writing new ones,
  and reuse helpers (`_misc_functions.R`, `helper_general.R`, `*_functions.R`).
- Scientific care matters more than speed. Respect the animal as the unit of
  replication, state statistical choices explicitly, and flag missing inputs
  instead of silently working around them (see `analysis_md/00_input_inventory.md`).
- Don't re-run expensive upstream steps (MAST, D-SPIN, pySCENIC GRN, scFEA)
  unless I ask. Load their saved outputs.
- Summaries back to me: what ran, where outputs landed, what's still pending,
  and any junk or inconsistencies you noticed (report them; don't fix them by deleting).

---

## 9. Skill routing

When the user's request matches an available skill, ALWAYS invoke it using the Skill
tool as your FIRST action. Do NOT answer directly, do NOT use other tools first.
The skill has specialized workflows that produce better results than ad-hoc answers.

Key routing rules:
- Product ideas, "is this worth building", brainstorming → invoke office-hours
- Bugs, errors, "why is this broken", 500 errors → invoke investigate
- Ship, deploy, push, create PR → invoke ship
- QA, test the site, find bugs → invoke qa
- Code review, check my diff → invoke review
- Update docs after shipping → invoke document-release
- Weekly retro → invoke retro
- Design system, brand → invoke design-consultation
- Visual audit, design polish → invoke design-review
- Architecture review → invoke plan-eng-review
- Save progress, checkpoint, resume → invoke checkpoint
- Code quality, health check → invoke health
