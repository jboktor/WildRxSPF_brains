# Figures for the integrated analysis. Sourced by run_integrated_analysis_sn.R --stage figures
# (expects `P` and helpers in scope). Light: reads saved tables only.
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(patchwork) })

# Readability defaults: large base text, full-size legend text and keys, nothing squeezed.
th <- theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 13, colour = "#555555"),
        plot.title.position = "plot", legend.position = "top", legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"), legend.key.size = unit(1, "lines"),
        axis.text = element_text(size = 13, colour = "#0A0A0A"), axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"), strip.background = element_rect(fill = "#FFF3EB", colour = NA),
        panel.grid.minor = element_blank(), plot.margin = margin(12, 16, 12, 12))
sv <- function(p, f, w, h) ggsave(file.path(P$fig, f), p, width = w, height = h, dpi = 200, bg = "white")
rd <- function(f) read_csv(file.path(P$tbl, f), show_col_types = FALSE)
ORANGE <- "#FF6C0C"; GREY <- "#A5A5A5"; SLATE <- "#27455C"
lvl <- c("class", "subclass", "supertype")

# ── Fig 1 · gene bridge: does the program popAlign sees moving carry the MAST signal? ──
gb <- rd("integrated_gene_bridge.csv") %>% mutate(level = factor(str_to_title(level), str_to_title(lvl)))
gbu <- gb %>% distinct(unit_id, shifted_feature, .keep_all = TRUE)       # one row per unit × shifted program
long <- gbu %>% filter(!is.na(mast_dirz_shift), !is.na(mast_dirz_control)) %>%
  transmute(unit_id, level, shift = mast_dirz_shift, control = mast_dirz_control) %>%
  pivot_longer(c(control, shift), names_to = "program", values_to = "z") %>%
  mutate(program = factor(if_else(program == "shift", "Most-shifted", "Least-shifted\n(control)"), c("Least-shifted\n(control)", "Most-shifted")))
pw <- gbu %>% filter(!is.na(mast_dirz_shift), !is.na(mast_dirz_control)) %>% group_by(level) %>%
  summarise(p = wilcox.test(mast_dirz_shift, mast_dirz_control, paired = TRUE)$p.value, n = n(), .groups = "drop") %>%
  mutate(lab = glue("n = {n}\npaired p = {format(signif(p, 2), scientific = TRUE)}"))
ylim <- quantile(long$z, c(0.01, 0.99)); ylim[2] <- ylim[2] + 3.2   # headroom for the p-value labels
p1a <- ggplot(long, aes(program, z, fill = program)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.6) +
  geom_violin(alpha = 0.3, colour = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.6) +
  geom_text(data = pw, aes(x = 1.5, y = ylim[2], label = lab), inherit.aes = FALSE, size = 4.4, vjust = 1, lineheight = 0.9) +
  facet_wrap(~ level) + coord_cartesian(ylim = ylim) +
  scale_fill_manual(values = c(GREY, ORANGE), guide = "none") +
  labs(x = "popAlign program in the same component", y = "MAST direction agreement (z)",
       title = "A · popAlign's moving programs carry the MAST expression signal",
       subtitle = "z > 0: the program's top 50 genes change in MAST in the same direction the program moves in popAlign") + th
ag <- gbu %>% summarise(
  `MAST: genes shift in\nthe same direction (z > 0)` = c(mean(mast_dirz_shift > 0, na.rm = TRUE), mean(mast_dirz_control > 0, na.rm = TRUE)),
  `D-SPIN: matches a program\nthat responds (p < 0.05)` = c(mean(dspin_any_responsive_overlap_shift, na.rm = TRUE), mean(dspin_any_responsive_overlap_control, na.rm = TRUE)),
  `pySCENIC: ≥ 3 regulons\ntarget the program` = c(mean(scenic_n_enriched_shift >= 3, na.rm = TRUE), mean(scenic_n_enriched_control >= 3, na.rm = TRUE))) %>%
  mutate(arm = c("Most-shifted program", "Least-shifted program (control)")) %>%
  pivot_longer(-arm) %>% mutate(arm = factor(arm, c("Least-shifted program (control)", "Most-shifted program")))
p1b <- ggplot(ag, aes(value, name, fill = arm)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = scales::percent(value, 1)), position = position_dodge(0.75), hjust = -0.15, size = 4.6) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = c(GREY, ORANGE), name = NULL) + guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Share of popAlign components", y = NULL,
       title = "B · Which other method sees the same program move?") + th
sv(p1a / p1b + plot_layout(heights = c(1.2, 1)), "integrated_fig1_gene_bridge.png", 13, 12)

# ── Fig 2 · MAST at the animal level, and uniform vs subpopulation-confined effects ──
pg <- readRDS(file.path(P$out, "pseudobulk_unit_genes.rds"))
mast_all <- map_dfr(unique(pg$unit_id), function(u) { p <- str_split(u, "__")[[1]]
  m <- load_mast_unit(P, "supertype", p[2], p[3]); if (is.null(m)) NULL else mutate(m, unit_id = u) })
mj <- inner_join(mast_all, pg, by = c("unit_id", "gene"), suffix = c("_mast", "_pb"))
grp <- mj %>% mutate(g = case_when(padj_mast < 0.1 ~ "MAST DEGs\n(padj < 0.1)", padj_mast >= 0.5 ~ "Not DE in MAST\n(padj ≥ 0.5)", TRUE ~ NA_character_)) %>%
  filter(!is.na(g)) %>% group_by(g) %>%
  summarise(`Same direction` = mean(sign(logFC) == sign(coef)), `Animal-level p < 0.05` = mean(p_pb < 0.05),
            `Animal-level FDR < 0.1` = mean(padj_pb < 0.1), n = n(), .groups = "drop") %>%
  pivot_longer(-c(g, n)) %>% mutate(name = factor(name, c("Same direction", "Animal-level p < 0.05", "Animal-level FDR < 0.1")),
                                   g = factor(g, c("Not DE in MAST\n(padj ≥ 0.5)", "MAST DEGs\n(padj < 0.1)")))
p2a <- ggplot(grp, aes(value, name, fill = g)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = scales::percent(value, 0.1)), position = position_dodge(0.75), hjust = -0.12, size = 4.6) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c(GREY, SLATE), name = NULL) + guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = glue("Share of genes (MAST DEGs: n = {format(sum(mj$padj_mast < 0.1), big.mark = ',')}, 117 supertypes)"), y = NULL,
       title = "A · MAST DEGs re-tested with animals as replicates (pseudobulk, 6 vs 6)",
       subtitle = "Nominal significance is enriched ~9×, but nothing survives FDR. Direction agrees even for non-DE genes (shared data)") + th
pc <- rd("integrated_pseudobulk_component_summary.csv") %>%
  left_join(gb %>% dplyr::select(unit_id, comp, p_state), by = c("unit_id", "comp")) %>%
  mutate(grp = case_when(!robust ~ "Non-robust", p_state < 0.05 ~ "State-shifted\n(popAlign p < 0.05)", TRUE ~ "Other robust")) %>%
  mutate(grp = factor(grp, c("Non-robust", "Other robust", "State-shifted\n(popAlign p < 0.05)")))
pcs <- pc %>% group_by(grp) %>% summarise(m = mean(excess_p01), se = sd(excess_p01) / sqrt(n()), n = n(), .groups = "drop")
p2b <- ggplot(pcs, aes(grp, m, fill = grp)) + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.18, linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_text(aes(y = 0.05, label = glue("n = {n}")), vjust = 0, size = 4.4, colour = "white", fontface = "bold") +
  annotate("text", x = 0.55, y = 1.04, label = "chance", hjust = 0, vjust = 0, size = 4.2, colour = "grey30") +
  scale_fill_manual(values = c("#D0D0D0", GREY, ORANGE), guide = "none") +
  labs(x = "Subpopulation (popAlign component)", y = "Animal-level DE signal\n(observed / expected p < 0.01)",
       title = "B · Subpopulations popAlign flags carry more DE",
       subtitle = "Compare groups, not the line: the 6-vs-6 test is\nconservative here, so most units sit below chance") + th
het <- rd("integrated_heterogeneity_summary.csv")
p2c <- ggplot(het, aes(excess_het_p01)) + geom_histogram(bins = 25, fill = "#FF9240", colour = "white") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey30") +
  annotate("text", x = 1.05, y = Inf, label = "chance", hjust = 0, vjust = 1.5, size = 4.2, colour = "grey30") +
  labs(x = "Microbiome × subpopulation interaction\n(observed / expected p < 0.01)", y = "Supertype units",
       title = "C · Effects are uniform across subpopulations",
       subtitle = glue("Only {sum(het$n_het_fdr10 > 0)} of {nrow(het)} units have any gene whose\nmicrobiome effect differs between subpopulations (FDR < 0.1)")) + th
sv(p2a / (p2b | p2c) + plot_layout(heights = c(0.9, 1.1)), "integrated_fig2_animal_level_de.png", 15, 11)

# ── Fig 3 · does each readout concentrate in the MeA–BST–medial hypothalamus network? ──
ce <- rd("integrated_circuit_enrichment.csv") %>%
  mutate(lab = glue("{modality}"), sig = p_adj < 0.05,
         lab = fct_reorder(lab, ratio_adj),
         rates = glue("{scales::percent(rate_circuit, 0.1)} vs {scales::percent(rate_other, 0.1)}"))
p3 <- ggplot(ce, aes(ratio_adj, lab)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = ratio_adj_lo, xmax = pmin(ratio_adj_hi, 40), colour = inference), height = 0.22, linewidth = 1.1) +
  geom_point(aes(colour = inference, shape = sig), size = 4.5, stroke = 1.4) +
  geom_text(aes(x = 60, label = glue("{rates}   p = {signif(p_adj, 2)}")), hjust = 0, size = 4.5) +
  scale_x_log10(breaks = c(0.3, 1, 3, 10, 30), labels = c("0.3", "1", "3", "10", "30"), limits = c(0.25, 900)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = c(animals = ORANGE, nuclei = SLATE), labels = c(animals = "animals (6 vs 6)", nuclei = "nuclei"), name = "Replicates:") +
  annotate("text", x = 60, y = nrow(ce) + 0.6, label = "hit rate: network vs other", hjust = 0, size = 4.3, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(x = "Hit-rate ratio, network vs other supertypes (adjusted for unit size and tissue; log scale)", y = NULL,
       title = "Only popAlign's state shifts concentrate in the MeA → BST → medial hypothalamus network",
       subtitle = "Hit rate = share of a unit's tests with p < 0.05 (MAST: padj < 0.1). Quasi-binomial model on hit counts; 95% CI") +
  th + theme(plot.margin = margin(30, 16, 12, 12))
sv(p3, "integrated_fig3_circuit_by_modality.png", 14, 6)

# ── Fig 4 · animal level: generalisation and coupling ───────────────────────
pr <- rd("integrated_animal_loao_prediction.csv") %>%
  mutate(scope = if_else(scope == "network", "Network units", "Other units"),
         modality = recode(modality, Expression = "Pseudobulk expression", `D-SPIN` = "D-SPIN programs", pySCENIC = "pySCENIC regulons", popAlign = "popAlign state"),
         sig = p_one_sided < 0.05)
p4a <- ggplot(pr, aes(auc, modality, colour = scope)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = sig), size = 4.5, stroke = 1.3, position = position_dodge(width = 0.55)) +
  facet_wrap(~ tissue) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "one-sided p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = c(`Network units` = ORANGE, `Other units` = GREY), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Held-out animal placed by the direction learned from the other 11 (AUC; 0.5 = chance)", y = NULL,
       title = "A · Only popAlign's network state direction places held-out animals (significant in HYP)",
       subtitle = "Leave-one-animal-out AUC. Leave-one-out is biased BELOW 0.5 under the null, so values < 0.5 are expected noise;\nonly AUC > 0.5 is tested (one-sided, all 924 relabelings)") + th
cp <- rd("integrated_animal_coupling.csv") %>% filter(!str_detect(m1, "^VR"), !str_detect(m2, "^VR")) %>%
  mutate(across(c(m1, m2), ~ recode(.x, Expression = "Pseudobulk", `D-SPIN` = "D-SPIN", pySCENIC = "pySCENIC", popAlign = "popAlign")),
         scope = if_else(scope == "network", "Network units", "Other units"), pair = glue("{m1} – {m2}"))
p4b <- ggplot(cp, aes(rho_within, fct_rev(factor(pair)), colour = scope)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(aes(shape = p_within < 0.05), size = 4.3, stroke = 1.3, position = position_dodge(width = 0.55)) +
  facet_wrap(~ tissue) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "permutation p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = c(`Network units` = ORANGE, `Other units` = GREY), name = NULL) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Spearman ρ between modality scores, within group (group means removed)", y = NULL,
       title = "B · Modalities agree on which animals are extreme, everywhere",
       subtitle = "Coupling is as strong outside the network as inside: shared animal-to-animal variation, not a network signal") + th
sv(p4a / p4b + plot_layout(heights = c(0.8, 1.2)), "integrated_fig4_animal_level.png", 14, 13)
log_msg("Figures written to ", P$fig)
