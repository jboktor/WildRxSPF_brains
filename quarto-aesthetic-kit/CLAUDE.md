# Analysis Notebook Style Guide — Quarto (.qmd)

> **For Claude Code / Cowork:** This file is the authoritative spec for every
> `.qmd` analysis notebook we produce. When generating or editing a `.qmd`, load
> this file into context and follow it exactly. Do not improvise structure,
> styling, or section order. The goal is that every notebook looks and reads as
> if one person made it, in our consistent house aesthetic.

---

## 0. How to use this kit in a session

This folder is self-contained and portable. Drop it into any analysis project:

```
your-project/
├── CLAUDE.md              ← this guide (Claude reads it automatically)
├── _quarto.yml            ← project styling defaults (do not edit per-notebook)
├── theme/
│   └── aesthetic.scss        ← the visual theme (the look & feel)
├── templates/
│   └── analysis-template.qmd   ← copy this to start a new notebook
└── _output/               ← rendered .html lands here
```

**To reproducibly apply it, tell Claude one of:**

- *"Start a new analysis notebook from the aesthetic-kit template."* → Claude copies
  `templates/analysis-template.qmd`, renames it, and fills in your analysis.
- *"Style this notebook to our standard."* → Claude restructures an existing
  `.qmd` into the section order in §2 and wires it to `theme/aesthetic.scss`.

**Two ways to attach the theme:**

1. **Project mode (preferred):** keep `_quarto.yml` at the project root. Every
   `.qmd` in the folder inherits the theme automatically — nothing to declare
   per file.
2. **Single-file mode:** if there is no `_quarto.yml`, paste the front-matter
   block from §3 into the top of the `.qmd` and ship `theme/aesthetic.scss`
   alongside it.

**Render with:** `quarto render notebook.qmd` (or `quarto preview` while editing).

Rule of thumb for Claude: **never hand-write CSS or colors into a notebook.**
All styling comes from `theme/aesthetic.scss`. Notebooks only use the documented
component classes (§5).

---

## 1. The non-negotiables (read first)

1. **One H1, stating exactly what this notebook does.** Title-case, specific,
   no fluff. e.g. `# Differential Expression of Tumor-Infiltrating T Cells`.
2. **Results live at the TOP. Code lives at the BOTTOM.** A reader sees the
   story — figures, takeaways, tables — before any code. The executable code is
   collected in a single appendix at the end (§4).
3. **Every figure and table is traceable to the cell that made it** (provenance,
   §6). A reader can always answer "what code produced this?".
4. **Fixed section order:** Overview → Background → Results → Next Steps →
   Code & Provenance. Always. (§2)
5. **The aesthetic is fixed.** Inter typeface, the Caltech orange palette, soft
   white cards. Use the component classes; don't invent new ones.

---

## 2. Mandatory document structure

Every notebook follows this exact spine. Section titles are verbatim.

### `# <Notebook Title>`
A single H1 — the clearest possible statement of what happened in this notebook.

### `## Overview`
**One short paragraph.** What we accomplished and what this notebook contains.
A busy reader should understand the whole notebook from this alone. Optionally
follow with a `.stat-row` of 2–4 headline numbers (n samples, key effect size,
etc.). No background, no methods here — just the "what."

### `## Background`
**One, at most two short paragraphs.** Why we ran this analysis — the
motivation, the question, and the value we expect to extract. This is the
"why." No results here.

### `## Results`
The heart of the notebook. A sequence of **figure cards**, optionally
interleaved with **data tables** and **key file paths**, ordered so the story
flows. For each figure:

- A **multi-panel figure** inside a `.figure-card` (panels labelled A, B, C…).
- A **caption** (`.fig-caption`) that describes the analysis behind the figure —
  what was computed and how, in 1–3 sentences.
- A **Key takeaways** block (`.takeaways`) directly under the caption: sharp,
  concise, scientific bullets stating what the figure *answers* or *resolves*
  about the question we asked. Not descriptions of the figure — conclusions.
- A **provenance tag** linking to the code cell that generated it (§6).

Interleave key **data tables** and **file-path chips** (`.path`) wherever they
make the story most interpretable — e.g. the processed-data path right after the
QC figure, the DE results table right after the volcano plot.

### `## Next Steps`
A short bullet list. First 1–2 bullets recap what we did; the rest state where
this should go next and why that's the sensible direction.

### `# Code & Provenance` {.code-appendix}
The full, runnable notebook. Every executable chunk lives here, in order, each
**labelled** so the results above can point back to it (§4, §6).

---

## 3. Front-matter (single-file mode)

If not using the project `_quarto.yml`, every notebook starts with:

`categories` render as the orange **eyebrow pills** above the title — use them
as the kicker (analysis type, data modality, project).

```yaml
---
title: "Differential Expression of Tumor-Infiltrating T Cells"
subtitle: "Bulk RNA-seq, responders vs. non-responders"
date: today
author: "Aisona Research"
categories: [RNA-seq, Differential Expression]   # -> eyebrow pills above title
format:
  html:
    theme: [cosmo, theme/aesthetic.scss]
    include-in-header: theme/fonts.html   # loads Inter + JetBrains Mono webfonts
    grid:
      body-width: 1280px                  # wide canvas (see §8)
      sidebar-width: 180px
      margin-width: 240px
    toc: true
    toc-location: right
    toc-title: "On this page"
    page-layout: article
    code-fold: true
    code-tools: true
    code-copy: true
    embed-resources: true
    fig-cap-location: bottom
execute:
  echo: true
  warning: false
  message: false
  freeze: auto
---
```

`embed-resources: true` produces a single portable `.html` you can email or drop
in Slack. `code-tools: true` adds the global show/hide-code toggle (top-right).

---

## 4. The "results-up, code-down" mechanism

This is the most important mechanical rule. We use the **narrative-first**
pattern: code runs in the appendix and writes its outputs to disk; the Key
Results section embeds those saved artifacts.

**In the appendix (`# Code & Provenance`):** each analysis chunk is labelled and
saves its figure/table to disk.

````markdown
```{r}
#| label: cell-fig-volcano
#| fig-show: hide
library(ggplot2)
p <- ggplot(de, aes(log2FC, -log10(padj))) + geom_point()
ggsave("figures/volcano.png", p, width = 7, height = 5, dpi = 200)
write.csv(de_top, "tables/de_top_genes.csv", row.names = FALSE)
```
````

**In Results (top):** embed the saved image inside a figure card and point
back to the cell.

```markdown
::: {.figure-card}
::: {.panel-grid}
[A]{.panel-label}
![](figures/volcano.png)

[B]{.panel-label}
![](figures/heatmap.png)
:::

::: {.fig-caption}
[**Figure 1.**]{.fig-num} Differential expression between responders and
non-responders. (A) Volcano plot of 18,402 genes; (B) z-scored expression of the
top 40 DE genes.
:::

::: {.takeaways}
**Key takeaways**

- 312 genes are significantly upregulated in responders (padj < 0.05, |log2FC| > 1).
- The signature is dominated by cytotoxicity programs (GZMB, PRF1, IFNG).
- Non-responders show no comparable enrichment — the effect is responder-specific.
:::

[Source: cell `cell-fig-volcano`](#cell-fig-volcano){.provenance}
:::
```

**Why this pattern (use it by default):**

- Code is unambiguously at the bottom; results are unambiguously at the top.
- The narrative renders even before code finishes (artifacts are on disk).
- Provenance is a literal hyperlink from result → generating cell.

**Alternative (compute-once):** for heavy compute, run it in one appendix chunk
with `#| output: false`, store results in variables, and add tiny display chunks.
Prefer the narrative-first pattern above unless compute cost forces otherwise.

Always create `figures/` and `tables/` at the top of the appendix:

````markdown
```{r}
#| label: setup
#| include: false
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
```
````

---

## 5. Component reference (copy-paste)

All classes are defined in `theme/aesthetic.scss`. These are the only building
blocks a notebook should use.

**Eyebrow label** (uppercase kicker above a heading):
```markdown
[METHOD]{.eyebrow}
```

**Stat row** (headline numbers under Overview):
```markdown
::: {.stat-row}
::: {.stat}
[1,284]{.stat-value}
[Cells passing QC]{.stat-label}
:::
::: {.stat}
[312]{.stat-value}
[DE genes (padj<0.05)]{.stat-label}
:::
:::
```

**Figure card + panel grid** — see §4. Single-panel: omit `.panel-grid` and
`.panel-label`, just put one image in the `.figure-card`. Add `.fig-wide`
(`::: {.figure-card .fig-wide}`) for a figure that should span the full content
width — wide bar charts, time series, heatmaps. Multi-panel grids already fill
the column.

**Eyebrow / kicker** — two options: the document-level kicker is the front-matter
`categories:` list (renders as orange pills above the title). For an inline
section kicker, put `[METHOD]{.eyebrow}` immediately above a heading.

**Key takeaways** — see §4. Always a `.takeaways` block with a bold
"Key takeaways" first line, then 2–5 bullets.

**Provenance tag**:
```markdown
[Source: cell `cell-fig-umap`](#cell-fig-umap){.provenance}
```

**File-path chip** (inline):
```markdown
Processed object written to [data/processed/cells.h5ad]{.path}.
```

**Data table** — author as a normal Markdown/`knitr::kable()`/`gt` table; the
theme styles it (orange header rule, hover rows). Give it a caption with
`: My caption {#tbl-id}` so it can be cross-referenced and traced.

**Section divider**:
```markdown
:::: {.section-rule}
::::
```

**Callouts** (use sparingly, for caveats): standard Quarto
`::: {.callout-note}` / `callout-warning` — recolored to brand automatically.

---

## 6. Provenance rules

Provenance is mandatory and has two halves:

1. **Every appendix code chunk is labelled** with a stable `#| label:` of the
   form `cell-<thing>` (e.g. `cell-fig-umap`, `cell-tbl-de`, `cell-qc`). The
   label becomes the chunk's HTML anchor (`#cell-fig-umap`).
2. **Every result up top carries a `.provenance` link** to its generating cell's
   anchor. Clicking it smooth-scrolls to the code and highlights it.

Conventions:

- Figures: cell label `cell-fig-<name>`, saved to `figures/<name>.png`.
- Tables: cell label `cell-tbl-<name>`, saved to `tables/<name>.csv`.
- One chunk = one primary artifact. Don't bury three figures in one cell; it
  breaks the 1:1 result→code mapping.
- If a figure path is referenced in the text, also surface it as a `.path` chip
  so a reader knows where the artifact lives on disk.

---

## 7. Readability & flow rules

- **Audience is mixed.** Overview, Background, captions, and takeaways must read
  cleanly for a non-technical reader. Code and methods detail stay in the
  appendix for the technical reader.
- **Prose, not bullet-dumps,** in Overview and Background. Bullets are for Key
  takeaways and Next Steps.
- **Takeaways are conclusions, not descriptions.** "Responders show a cytotoxic
  signature" — not "the plot shows red dots on the right."
- **One idea per figure card.** If a figure needs more than ~4 takeaway bullets,
  it's probably two figures.
- **Captions describe the analysis;** takeaways describe the finding. Keep them
  distinct.
- **Tables are short by default** — show the top N rows, link the full CSV via a
  `.path` chip.
- **Whitespace is part of the design.** Don't crowd cards; let sections breathe.
- **Let content fill the wide column.** Prose, figures, and tables all use the
  full ~1280px width; lean on whitespace and short paragraphs (rather than a
  narrow column) to keep things readable.

---

## 8. Aesthetic spec (the source of truth, mirrored in aesthetic.scss)

Claude should not need these numbers day-to-day (the theme encodes them), but
they define the look and must not drift.

**Typeface**

| Role | Font | Notes |
|---|---|---|
| All UI/body/headings | **Inter** (400/500/600/700) | from Google Fonts |
| Code / monospace | **JetBrains Mono** | code, paths, provenance tags |

**Color palette**

| Token | Hex | Use |
|---|---|---|
| `--bh-orange` | `#FF6C0C` | primary brand: links, accents, table rule, takeaway marker |
| `--bh-orange-light` | `#FF9240` | hover states, secondary accents |
| `--bh-orange-tint` | `#FFF3EB` | takeaway / highlight / inline-code backgrounds |
| ink | `#0A0A0A` | body text |
| `--bh-gray` | `#767676` | secondary text |
| muted | `#555555` | captions |
| `--bh-gray-light` | `#A5A5A5` | eyebrow labels |
| `--bh-border` | `#E5E5E5` | hairline borders |
| `--bh-surface` | `#F8F8F8` | card / code / panel fill |
| `--bh-slate` | `#27455C` | dark accent (paths), use sparingly |
| white | `#FFFFFF` | page + card background |

**Type scale** (document base = 17px / line-height 1.62): H1 2.55rem/700,
H2 1.7rem/600, H3 1.34rem/600, H4 1.12rem. Headings use slightly negative
letter-spacing (−0.015em); eyebrow labels use +0.06em uppercase.

**Shape & depth:** card radius 8px, small radius 4px, pills 20px. Shadows are
soft and low: resting `0 2px 10px -3px rgba(120,120,120,.35)`; figure-card hover
lifts 2px with a faint orange glow. Borders are 1px hairlines, never heavy.

**Motion / reactivity:** links and cards transition on hover (~0.15–0.2s ease);
figure cards lift on hover; provenance jumps use `scroll-behavior: smooth`; the
active TOC item is marked in orange with a left rule. Keep motion subtle — this
is a research document, not a marketing page.

**Layout (wide canvas, not a blog column):** the content column is widened to
**1280px** via Quarto's `grid.body-width` (default is ~800px), and **everything —
prose, headings, figure cards, panel grids, KPI strips, tables — fills that full
width**. Type is set at 17px / 1.72 line-height for comfort across the wider
measure. White background, generous whitespace, a right-hand sticky TOC titled
"On this page," figures centered, captions below.

The width levers live in `_quarto.yml`:

```yaml
format:
  html:
    page-layout: article
    grid:
      body-width: 1280px     # widen the center column (raise toward 1440 if needed)
      sidebar-width: 180px
      margin-width: 240px
      gutter-width: 1.8rem
```

Prose fills the full column by default (`--bh-measure: none` in `aesthetic.scss`).
If you ever want a calmer, capped reading measure for a text-heavy notebook, set
`--bh-measure` back to a value like `74ch`. For a *single* figure that should span
the full width (e.g. a wide time series or per-sample bar chart), add the
`.fig-wide` modifier to its card: `::: {.figure-card .fig-wide}`. Without it,
single figures center at a calmer max width; multi-panel grids always fill the
column.

---

## 9. Quality checklist (Claude self-checks before delivering)

- [ ] Exactly one H1, and it states what the notebook does.
- [ ] Sections present in order: Overview, Background, Results, Next Steps,
      Code & Provenance.
- [ ] Overview is one paragraph; Background ≤ 2 paragraphs.
- [ ] Every figure is in a `.figure-card` with a caption **and** a `.takeaways`
      block.
- [ ] Every figure/table has a `.provenance` link to a labelled appendix cell.
- [ ] All executable code is in the `# Code & Provenance` appendix, nowhere else.
- [ ] Key file paths surfaced as `.path` chips.
- [ ] Theme attached via `theme/aesthetic.scss`; no inline CSS/colors in the .qmd.
- [ ] Renders cleanly: `quarto render` with no errors; `_output/*.html` opens
      and is self-contained.
