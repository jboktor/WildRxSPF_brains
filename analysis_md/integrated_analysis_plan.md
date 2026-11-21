# Integrated analysis plan: MAST × pySCENIC × D-SPIN × popAlign

*Written 2026-10-06. Motivated by the popAlign results (`notebooks/0X_popAlign_SN.qmd`, run
2026-10-04, synthesis 2026-10-05): composition is stable, state shifts concentrate in the
MeA → BST → medial hypothalamus network, and popAlign leads overlap MAST / D-SPIN only weakly at the
level of identity markers. Implemented in `notebooks/0X_Integrated-Analysis_SN.qmd`.*

## Plan

The next step is to integrate MAST, pySCENIC, D-SPIN and popAlign on a shared scaffold: the same
tissue × Allen-supertype units (which D-SPIN, popAlign and MAST already share), the same nuclei,
and the animal as the unit of replication. The goal is to move from four separate lists of hits to
one model of how the microbiome acts, running from regulator to program to cell state to gene.

1. **Gene-level bridge.** For every unit, take the oNMF program that moves most between WildR and
   SPF in popAlign, and test its top-weighted genes for overlap and direction-of-effect agreement
   with three gene sets from the same unit: MAST DEGs, the genes of the D-SPIN programs with the
   largest WildR–SPF response differences, and the target sets of pySCENIC regulons. This shows
   whether popAlign's state shifts are the same biology as the DEGs and network responses, and the
   regulon overlap nominates candidate upstream transcription factors.
2. **Cell-level bridge.** Because AUCell, D-SPIN program activity and popAlign component
   assignments are all defined per nucleus, score regulon activity and D-SPIN program activity
   within each popAlign subpopulation. Then re-run MAST-style pseudobulk tests per animal ×
   subpopulation. This asks whether a microbiome effect is uniform across a cell type or confined
   to one subpopulation, and lets us order the evidence: regulon activity (pySCENIC) → coordinated
   program response (D-SPIN) → subpopulation state (popAlign) → individual DEGs (MAST).
3. **Animal-level integration.** Summarise each of the 12 animals by its popAlign state score per
   unit (leave-one-animal-out axis, so no animal scores itself), D-SPIN response vector `h`, mean
   regulon AUC, a MAST signature score, and the VR-1/VR-2 relay scores. Use rank correlations and
   multi-block summaries with permutation tests to ask whether one shared microbiome axis explains
   variation across modalities, and whether the animals most shifted in one modality are also the
   most shifted in the others.
4. **Circuit test across modalities.** Ask whether MAST DEG burden, D-SPIN response magnitude and
   regulon shifts concentrate in the MeA → BST → medial hypothalamus network the way popAlign state
   shifts do, using the same unit-level enrichment tests and size/level/tissue adjustments.
   Agreement among methods built on different principles would make the circuit result far more
   credible than any single analysis.

**Output.** One integrated evidence table per unit and lead, ranked by how many modalities agree,
used to decide which populations and genes go forward to seqFISH validation.

**Caveat.** All four methods analyse the same 12 animals, so integration makes the results more
interpretable and coherent but is not independent replication. Confirmatory tests remain the
spatial seqFISH data and, ideally, a new cohort.
