# VR-1 results

## Headline (decision rule)

- Per-animal separation on `D`: **0 / 6** WildR exceed all SPF.
- Cliff's δ (WildR vs SPF on `D`): **0.611** (95% bootstrap CI -0.056 – 1.000).
- Spearman(input, output) across 12 animals: **ρ = -0.084** (CI -0.677 – 0.527).
- Wilcoxon rank-sum p (not headline): 0.0927.

## Branch

**Decision:** no clear separation on `D` → retire the coupled-circuit sentence (valuable negative).

## Technical covariates

- Spearman(input, median UMI)=NA; Spearman(output, median UMI)=NA
- Spearman(input, n_nuclei)=-0.329; Spearman(output, n_nuclei)=-0.392

## Signature

- Positive genes (n=344): Eps8l1, Zbtb16, Hdac9, Rasa1, 9330159F19Rik, A830018L16Rik, Abcb4, Abcc5, Abcf1, Acsl4, Actr3b, Adam23, Adcy1, Adcy8, Agap2, Ahi1, Alcam, Alg10b, Angpt1, Ankfn1, Ankrd10, Ap3d1, Apool, Arap2, Arhgap21, Arhgef17, Arhgef2, Arpp21, Asic2, Atf7ip, Atp1a1, Atrx, Atxn2l, BC005537, Basp1, Bhlhe40, Brd9, Cacna2d1, Cadm3, Calm1, Calm2, Camk4, Camta2, Ccdc186, Ccdc88a, Cct5, Cdc42, Cdc42bpa, Cdc5l, Cdh10
- Signed-negative genes (n=1): Itpkb
- Scoring: `scanpy.tl.score_genes` on per-unit `raw_counts.h5ad`.
- Aggregation: equal weight across glutamatergic AMY relay supertypes.
- Comparison: **between-animal** WildR vs SPF on within-animal `D`.

