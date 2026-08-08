# Quarto Aesthetic Kit

A reusable styling + structure kit for analysis notebooks (`.qmd`), in the
visual language of a clean, modern reference style. Results read at the top, code
lives in an appendix at the bottom, and every figure traces back to the cell
that made it.

## Contents

| File | What it is |
|------|------------|
| `CLAUDE.md` | **The guide.** Authoritative spec for structure, aesthetic, and provenance. Claude reads this when generating any `.qmd`. |
| `_quarto.yml` | Project styling defaults — every `.qmd` in the folder inherits the theme. |
| `theme/aesthetic.scss` | The visual theme (Inter, Caltech orange palette, cards, takeaways, provenance tags). |
| `templates/analysis-template.qmd` | A fully worked, runnable example to copy when starting a new notebook. |

## Quick start

1. Install [Quarto](https://quarto.org) (`quarto --version` to check).
2. Copy `templates/analysis-template.qmd` to a new name and edit it.
3. `quarto render your-notebook.qmd` → opens a self-contained `.html` in `_output/`.

The template runs in Python (needs `numpy`, `pandas`, `scipy`, `matplotlib`,
`jupyter`). Swap to R by changing the `{python}` cells to `{r}` — the styling and
structure are language-agnostic.

## Using it with Claude

Keep this folder in your project and just ask:

- *"Start a new analysis notebook from the aesthetic-kit template."*
- *"Style this notebook to our standard."*

Claude follows `CLAUDE.md` automatically.
