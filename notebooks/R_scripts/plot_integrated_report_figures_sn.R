#!/usr/bin/env Rscript
# Report figures for notebooks/0X_Integrated-Analysis_SN.qmd, multi-level pySCENIC update.
# Reads saved tables only:
#   tables/SN_Integrated/run_2026-10-06/   base integrated run (gene bridge, pseudobulk, heterogeneity)
#   tables/SN_Integrated/run_2026-10-07/   run_integrated_scenic_multilevel_sn.R (native-GRN pySCENIC,
#                                          circuit, animal level, effect sizes, unit concordance)
# Writes figures/SN_Integrated/run_2026-10-07/integrated_fig*.png. The 2026-10-06 figures (drawn by
# plot_integrated_analysis_sn.R) are left untouched as the record of the base run.
#
# Style follows plot_scenic_report_figures.R / mapmycells_qc_functions.R: 12 in wide, 200 dpi,
# theme_minimal base 13, panel tags A/B/C, one fixed colour per method.
#
# Run:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_integrated_report_figures_sn.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(scales); library(patchwork) })
wk <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
base_tbl <- file.path(wk, "tables/SN_Integrated/run_2026-10-06")
base_out <- file.path(wk, "data/interim/sn_integrated/run_2026-10-06")
new_tbl <- file.path(wk, "tables/SN_Integrated/run_2026-10-07")
fig_dir <- file.path(wk, "figures/SN_Integrated/run_2026-10-07")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
rb <- function(f) read_csv(file.path(base_tbl, f), show_col_types = FALSE)
rn <- function(f) read_csv(file.path(new_tbl, f), show_col_types = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")

# ---- Style ----
fig_w <- 12; fig_dpi <- 200
ORANGE <- "#E8590C"; BLUE <- "#2E6FA8"; SLATE <- "#27455C"; LORANGE <- "#FF9240"; GREY <- "#A5A5A5"; LGREY <- "#D9D9D9"
method_cols <- c(popAlign = ORANGE, Pseudobulk = SLATE, `D-SPIN` = LORANGE, pySCENIC = BLUE, MAST = "#7A7A7A")
scope_cols <- c(`Network units` = ORANGE, `Other units` = GREY)
theme_int <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.3),
          axis.text = element_text(colour = "#555555"), axis.title = element_text(colour = "#555555"),
          plot.title = element_text(face = "bold", colour = "#0A0A0A"), plot.subtitle = element_text(colour = "#555555"),
          strip.text = element_text(face = "bold", colour = "#0A0A0A", hjust = 0),
          legend.position = "top", plot.title.position = "plot", plot.caption.position = "plot",
          plot.tag = element_text(face = "bold", size = base_size + 2),
          plot.background = element_rect(fill = "white", colour = NA))
}
save_fig <- function(p, file, h) { ggsave(file.path(fig_dir, file), p, width = fig_w, height = h, dpi = fig_dpi, bg = "white"); log_msg("  wrote ", file) }
p_lab <- function(p) if_else(p < 0.001, "p < 0.001", glue("p = {formatC(signif(p, 2), format = 'fg', digits = 2)}"))
method_of <- function(x) case_when(str_detect(x, "popAlign") ~ "popAlign", str_detect(x, "Pseudobulk|Expression") ~ "Pseudobulk",
                                   str_detect(x, "D-SPIN") ~ "D-SPIN", str_detect(x, "pySCENIC") ~ "pySCENIC", str_detect(x, "MAST") ~ "MAST")

# ════════════════════════════════════════════════════════════════════════════
# Fig 1 · gene bridge
# ════════════════════════════════════════════════════════════════════════════
gb <- rb("integrated_gene_bridge.csv")
gbu <- gb %>% distinct(unit_id, shifted_feature, .keep_all = TRUE)
ok <- gbu %>% filter(!is.na(mast_dirz_shift), !is.na(mast_dirz_control))
pw <- ok %>% group_by(level) %>% summarise(p = wilcox.test(mast_dirz_shift, mast_dirz_control, paired = TRUE)$p.value, n = n(), .groups = "drop") %>%
  mutate(strip = glue("{str_to_title(level)}\n{n} programs · paired {p_lab(p)}"))
long <- ok %>% left_join(pw %>% dplyr::select(level, strip), by = "level") %>%
  transmute(strip = factor(strip, pw$strip[match(c("class", "subclass", "supertype"), pw$level)]),
            `Most-shifted` = mast_dirz_shift, `Least-shifted (control)` = mast_dirz_control) %>%
  pivot_longer(-strip, names_to = "program", values_to = "z") %>%
  mutate(program = factor(program, c("Least-shifted (control)", "Most-shifted")))
p1a <- ggplot(long, aes(program, z, fill = program)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_violin(alpha = 0.3, colour = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.45) +
  facet_wrap(~ strip, ncol = 3) + coord_cartesian(ylim = quantile(long$z, c(0.01, 0.99))) +
  scale_fill_manual(values = c(GREY, ORANGE), guide = "none") +
  scale_x_discrete(labels = c("Least-shifted\n(control)", "Most-shifted")) +
  labs(x = NULL, y = "MAST direction agreement (z)",
       title = "popAlign's moving programs carry the MAST expression signal",
       subtitle = "z > 0: the program's top 50 genes change in MAST in the direction the program moves in popAlign") +
  theme_int() + theme(strip.text = element_text(size = 11))
nat <- rn("integrated_gene_bridge_scenic_native.csv") %>% distinct(unit_id, shifted_feature, .keep_all = TRUE)
share <- function(x) mean(x, na.rm = TRUE)
ag <- tribble(
  ~readout, ~arm, ~value,
  "MAST: top genes shift the\nsame way (z > 0)", "shift", share(gbu$mast_dirz_shift > 0),
  "MAST: top genes shift the\nsame way (z > 0)", "control", share(gbu$mast_dirz_control > 0),
  "D-SPIN: overlaps a program\nthat responds (p < 0.05)", "shift", share(gbu$dspin_any_responsive_overlap_shift),
  "D-SPIN: overlaps a program\nthat responds (p < 0.05)", "control", share(gbu$dspin_any_responsive_overlap_control),
  "pySCENIC, class GRN:\n≥ 3 regulons target it", "shift", share(gbu$scenic_n_enriched_shift >= 3),
  "pySCENIC, class GRN:\n≥ 3 regulons target it", "control", share(gbu$scenic_n_enriched_control >= 3),
  "pySCENIC, own GRN:\n≥ 3 regulons target it", "shift", share(nat$native_n_enriched_shift >= 3),
  "pySCENIC, own GRN:\n≥ 3 regulons target it", "control", share(nat$native_n_enriched_control >= 3),
  "pySCENIC, own GRN: a targeting\nregulon shifts per animal", "shift", share(nat$native_enriched_shifted_shift >= 1),
  "pySCENIC, own GRN: a targeting\nregulon shifts per animal", "control", share(nat$native_enriched_shifted_control >= 1)) %>%
  mutate(readout = fct_rev(fct_inorder(readout)),
         arm = factor(if_else(arm == "shift", "Most-shifted program", "Least-shifted program (control)"),
                      c("Least-shifted program (control)", "Most-shifted program")))
p1b <- ggplot(ag, aes(value, readout, fill = arm)) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  geom_text(aes(label = percent(value, 1)), position = position_dodge(0.78), hjust = -0.15, size = 3.6, colour = "#333333") +
  scale_x_continuous(labels = percent, limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c(GREY, ORANGE), name = NULL) + guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Share of popAlign components", y = NULL, title = "Which other method sees the same program move?") + theme_int()
save_fig(p1a / p1b + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "A"), "integrated_fig1_gene_bridge.png", 11.5)

# ════════════════════════════════════════════════════════════════════════════
# Fig 2 · MAST at the animal level; subpopulation-confined vs uniform effects
# ════════════════════════════════════════════════════════════════════════════
source(file.path(wk, "notebooks/R_scripts/integrated_sn_helpers.R"))
P <- int_paths("2026-10-06")
pg <- readRDS(file.path(base_out, "pseudobulk_unit_genes.rds"))
mast_all <- map_dfr(unique(pg$unit_id), function(u) { p <- str_split(u, "__")[[1]]
  m <- load_mast_unit(P, "supertype", p[2], p[3]); if (is.null(m)) NULL else mutate(m, unit_id = u) })
mj <- inner_join(mast_all, pg, by = c("unit_id", "gene"), suffix = c("_mast", "_pb"))
grp <- mj %>% mutate(g = case_when(padj_mast < 0.1 ~ "MAST DEGs (padj < 0.1)", padj_mast >= 0.5 ~ "Not DE in MAST (padj ≥ 0.5)")) %>%
  filter(!is.na(g)) %>% group_by(g) %>%
  summarise(`Same direction` = mean(sign(logFC) == sign(coef)), `Animal-level p < 0.05` = mean(p_pb < 0.05),
            `Animal-level FDR < 0.1` = mean(padj_pb < 0.1), .groups = "drop") %>%
  pivot_longer(-g) %>% mutate(name = factor(name, rev(c("Same direction", "Animal-level p < 0.05", "Animal-level FDR < 0.1"))),
                              g = factor(g, c("Not DE in MAST (padj ≥ 0.5)", "MAST DEGs (padj < 0.1)")))
p2a <- ggplot(grp, aes(value, name, fill = g)) + geom_col(position = position_dodge(0.78), width = 0.72) +
  geom_text(aes(label = percent(value, 0.1)), position = position_dodge(0.78), hjust = -0.12, size = 3.7, colour = "#333333") +
  scale_x_continuous(labels = percent, limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c(GREY, SLATE), name = NULL) + guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = glue("Share of genes ({comma(sum(mj$padj_mast < 0.1))} MAST DEGs in {n_distinct(mj$unit_id[mj$padj_mast < 0.1])} supertypes)"), y = NULL,
       title = "MAST DEGs re-tested with animals as replicates (pseudobulk, 6 vs 6)",
       subtitle = "Nominal significance is strongly enriched, but nothing survives FDR") + theme_int()
pc <- rb("integrated_pseudobulk_component_summary.csv") %>% left_join(gb %>% distinct(unit_id, comp, p_state), by = c("unit_id", "comp")) %>%
  mutate(grp = case_when(!robust ~ "Non-robust", coalesce(p_state < 0.05, FALSE) ~ "State-shifted", TRUE ~ "Other robust"),
         grp = factor(grp, c("Non-robust", "Other robust", "State-shifted")))
pcs <- pc %>% group_by(grp) %>% summarise(m = mean(excess_p01), se = sd(excess_p01) / sqrt(n()), n = n(), .groups = "drop")
p2b <- ggplot(pcs, aes(grp, m, fill = grp)) + geom_col(width = 0.62) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.15, linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_text(aes(y = 0.04, label = glue("n = {n}")), vjust = 0, size = 3.6, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c(LGREY, GREY, ORANGE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "popAlign subpopulation", y = "Observed / expected genes at p < 0.01",
       title = "Flagged subpopulations carry more DE",
       subtitle = "Pseudobulk per animal × subpopulation; mean ± SE") + theme_int()
het <- rb("integrated_heterogeneity_summary.csv")
p2c <- ggplot(het, aes(excess_het_p01)) + geom_histogram(bins = 25, fill = LORANGE, colour = "white") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35") +
  labs(x = "Microbiome × subpopulation interaction\n(observed / expected genes at p < 0.01)", y = "Supertype units",
       title = "Effects are shared across subpopulations",
       subtitle = glue("{sum(het$n_het_fdr10 > 0)} of {nrow(het)} supertypes have any gene with a\nsubpopulation-specific effect (FDR < 0.1)")) + theme_int()
save_fig(p2a / (p2b | p2c) + plot_layout(heights = c(0.75, 1.1)) + plot_annotation(tag_levels = "A"), "integrated_fig2_animal_level_de.png", 10.5)

# ════════════════════════════════════════════════════════════════════════════
# Fig 3 · network concentration per method (forest + aligned text table)
# ════════════════════════════════════════════════════════════════════════════
ce <- rn("integrated_circuit_enrichment.csv") %>%
  mutate(label = recode(modality, `MAST DEGs (nuclei)` = "MAST DEGs (nuclei)", `pySCENIC regulons (class GRN)` = "pySCENIC · class GRN (base run)",
                        `pySCENIC regulons (supertype GRN)` = "pySCENIC · supertype GRN", `pySCENIC regulons (own GRN, all levels)` = "pySCENIC · own GRN, all levels"),
         label = fct_reorder(label, ratio_adj), method = method_of(modality), sig = p_adj < 0.05,
         repl = if_else(inference == "nuclei", "nuclei", "animals (6 vs 6)"))
p3a <- ggplot(ce, aes(ratio_adj, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = ratio_adj_lo, xmax = pmin(ratio_adj_hi, 40), colour = method), height = 0.25, linewidth = 0.8) +
  geom_point(aes(colour = method, shape = sig), size = 3.6, stroke = 1.2, fill = "white") +
  scale_x_log10(breaks = c(0.3, 1, 3, 10, 30), labels = c("0.3", "1", "3", "10", "30")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = method_cols, guide = "none") +
  labs(x = "Hit-rate ratio, network vs other units\n(adjusted for unit size and tissue; log scale)", y = NULL) +
  theme_int()
tab <- ce %>% transmute(label, Network = percent(rate_circuit, 0.1), Other = percent(rate_other, 0.1), `p (adj.)` = p_lab(p_adj),
                        Replicates = repl) %>%
  pivot_longer(-label, names_to = "col") %>% mutate(col = factor(col, c("Network", "Other", "p (adj.)", "Replicates")))
p3b <- ggplot(tab, aes(col, label, label = value)) + geom_text(size = 3.6, colour = "#333333") +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) +
  theme_int() + theme(axis.text.y = element_blank(), panel.grid = element_blank(), axis.text.x.top = element_text(face = "bold", colour = "#0A0A0A"))
save_fig(p3a + p3b + plot_layout(widths = c(1.15, 1)) +
           plot_annotation(title = "Only popAlign's state shifts concentrate in the MeA → BST → medial hypothalamus network",
                           subtitle = "Hit rate = share of a unit's tests at p < 0.05 (MAST: padj < 0.1). Quasi-binomial model on hit counts, 95% CI",
                           theme = theme_int()), "integrated_fig3_circuit_by_modality.png", 5.4)

# ════════════════════════════════════════════════════════════════════════════
# Fig 4 · animal level: held-out prediction and cross-method coupling
# ════════════════════════════════════════════════════════════════════════════
pr <- rn("integrated_animal_loao_prediction.csv") %>%
  mutate(scope = if_else(scope == "network", "Network units", "Other units"),
         modality = recode(modality, Expression = "Pseudobulk", popAlign = "popAlign", pySCENIC = "pySCENIC (supertype GRN)"),
         modality = factor(modality, rev(c("popAlign", "D-SPIN", "Pseudobulk", "pySCENIC (supertype GRN)"))), sig = p_one_sided < 0.05)
p4a <- ggplot(pr, aes(auc, modality, colour = scope)) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "#F6F6F6") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = sig), size = 3.4, stroke = 1.1, position = position_dodge(width = 0.5)) +
  facet_wrap(~ tissue) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "one-sided p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = scope_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "Held-out animal placed by the direction learned from the other 11 (AUC)", y = NULL,
       title = "Only popAlign's network state direction places held-out animals",
       subtitle = "Leave-one-animal-out is biased below 0.5 under the null (shaded), so only AUC > 0.5 is tested (exact, 924 relabellings)") + theme_int()
cp <- rn("integrated_animal_coupling.csv") %>% filter(!str_detect(m1, "^VR"), !str_detect(m2, "^VR")) %>%
  mutate(across(c(m1, m2), ~ recode(.x, Expression = "Pseudobulk")), scope = if_else(scope == "network", "Network units", "Other units"),
         pair = fct_rev(factor(glue("{m1} – {m2}"))))
p4b <- ggplot(cp, aes(rho_within, pair, colour = scope)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_point(aes(shape = p_within < 0.05), size = 3.4, stroke = 1.1, position = position_dodge(width = 0.5)) +
  facet_wrap(~ tissue) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "permutation p < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = scope_cols, name = NULL) + scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Spearman ρ between method scores within group (group means removed)", y = NULL,
       title = "Methods agree on which animals are atypical, inside and outside the network",
       subtitle = "Shared animal-to-animal variation, not a microbiome or network signal") + theme_int()
save_fig(p4a / p4b + plot_layout(heights = c(0.7, 1.2)) + plot_annotation(tag_levels = "A"), "integrated_fig4_animal_level.png", 11.5)

# ════════════════════════════════════════════════════════════════════════════
# Fig 5 · native-GRN regulon activity per popAlign subpopulation
# ════════════════════════════════════════════════════════════════════════════
cs <- rn("integrated_scenic_native_component.csv") %>%
  mutate(grp = case_when(!robust ~ "Non-robust", coalesce(p_state < 0.05, FALSE) ~ "State-shifted", TRUE ~ "Other robust"),
         grp = factor(grp, c("Non-robust", "Other robust", "State-shifted")), level = factor(str_to_title(level), c("Class", "Subclass", "Supertype")))
set.seed(3)
ratio_ci <- function(d) { r <- sum(d$n_p05) / sum(d$perm_mean)
  b <- replicate(2000, { ix <- sample(nrow(d), replace = TRUE); sum(d$n_p05[ix]) / sum(d$perm_mean[ix]) })
  tibble(ratio = r, lo = quantile(b, 0.025), hi = quantile(b, 0.975), n = nrow(d)) }
rr <- bind_rows(cs %>% group_by(level, grp) %>% group_modify(~ ratio_ci(.x)) %>% ungroup(),
                cs %>% mutate(level = "All levels") %>% group_by(level, grp) %>% group_modify(~ ratio_ci(.x)) %>% ungroup()) %>%
  mutate(level = factor(level, c("All levels", "Class", "Subclass", "Supertype")))
p5a <- ggplot(rr, aes(grp, ratio, colour = grp)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.55, linewidth = 0.8) +
  geom_text(aes(y = lo, label = glue("n = {n}")), vjust = 1.6, size = 3.2, colour = "#555555") +
  facet_wrap(~ level, nrow = 1) +
  scale_colour_manual(values = c(GREY, SLATE, ORANGE), guide = "none") +
  scale_x_discrete(labels = c("Non-\nrobust", "Other\nrobust", "State-\nshifted")) +
  labs(x = NULL, y = "Regulons at p < 0.05:\nobserved / permutation expectation",
       title = "Regulon activity shifts where popAlign sees subpopulation states move",
       subtitle = "Per-animal regulon AUC in each subpopulation, from the GRN fitted on that unit's own nuclei; bootstrap 95% CI over subpopulations") +
  theme_int()
ec <- cs %>% filter(!is.na(p_perm_count))
p5b <- ggplot(ec, aes(p_perm_count, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  stat_ecdf(linewidth = 0.9) +
  scale_colour_manual(values = c(GREY, SLATE, ORANGE), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "Subpopulation permutation p (excess of regulon hits)", y = "Cumulative share",
       title = "Per-subpopulation excess",
       subtitle = glue("{sum(ec$p_perm_count < 0.05)} of {nrow(ec)} subpopulations at p < 0.05\n(≈ {round(0.05 * nrow(ec))} by chance); none at FDR < 0.1")) +
  theme_int() + theme(legend.position = "inside", legend.position.inside = c(0.75, 0.2), legend.background = element_rect(fill = "white", colour = NA))
nb <- rn("integrated_gene_bridge_scenic_native.csv") %>% filter(native_enriched_tested_shift > 0 | native_enriched_tested_control > 0) %>%
  distinct(unit_id, shifted_feature, .keep_all = TRUE) %>%
  summarise(across(c(native_enriched_shifted_shift, native_enriched_tested_shift, native_enriched_shifted_control, native_enriched_tested_control), ~ sum(.x, na.rm = TRUE))) %>%
  { tibble(arm = factor(c("Least-shifted (control)", "Most-shifted"), c("Least-shifted (control)", "Most-shifted")),
           k = c(.$native_enriched_shifted_control, .$native_enriched_shifted_shift), n = c(.$native_enriched_tested_control, .$native_enriched_tested_shift)) } %>%
  mutate(rate = k / n, lo = map2_dbl(k, n, ~ binom.test(.x, .y)$conf.int[1]), hi = map2_dbl(k, n, ~ binom.test(.x, .y)$conf.int[2]))
p5c <- ggplot(nb, aes(arm, rate, fill = arm)) + geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.5) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey45") +
  geom_text(aes(y = hi, label = glue("{k} / {comma(n)}")), vjust = -0.6, size = 3.4, colour = "#333333") +
  scale_fill_manual(values = c(GREY, ORANGE), guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(labels = c("Least-shifted\n(control)", "Most-shifted")) +
  labs(x = "popAlign program", y = "Targeting regulons shifted\nper animal (p < 0.05)", title = "Regulons that target the program",
       subtitle = "Own-GRN regulons enriched for the\nprogram's top genes (q < 0.05)") + theme_int()
save_fig(p5a / (p5b | p5c) + plot_layout(heights = c(0.9, 1)) + plot_annotation(tag_levels = "A"), "integrated_fig5_scenic_subpopulations.png", 10.5)

# ════════════════════════════════════════════════════════════════════════════
# Fig 6 · effect sizes of every method vs the exact animal-relabelling null; power
# ════════════════════════════════════════════════════════════════════════════
es <- rn("integrated_effect_size_summary.csv")
dens <- rn("integrated_effect_size_density.csv")
nl <- read_csv(file.path(wk, "data/interim/sn_integrated/run_2026-10-07/effect_size_permutation_null.csv"), show_col_types = FALSE)
g_crit <- qt(0.975, 10) * sqrt(1 / 3) * (1 - 3 / (4 * 12 - 9))
valid <- c("Pseudobulk: all genes", "popAlign programs", "D-SPIN programs", "pySCENIC regulons")
col6 <- c(`Pseudobulk: all genes` = SLATE, `popAlign programs` = ORANGE, `D-SPIN programs` = LORANGE, `pySCENIC regulons` = BLUE,
          `Pseudobulk: MAST DEGs` = "#7A7A7A")
dlv <- c(valid, "Pseudobulk: MAST DEGs")
# all-genes density = weighted mix of the two pseudobulk subsets
nn <- es %>% filter(tissue == "Both") %>% dplyr::select(method, n_tests)
w_deg <- nn$n_tests[nn$method == "Pseudobulk: MAST DEGs"]; w_oth <- nn$n_tests[nn$method == "Pseudobulk: other genes"]
allg <- dens %>% filter(method %in% c("Pseudobulk: MAST DEGs", "Pseudobulk: other genes")) %>%
  group_by(abs_g) %>% summarise(density = sum(density * if_else(method == "Pseudobulk: MAST DEGs", w_deg, w_oth)) / (w_deg + w_oth), .groups = "drop") %>%
  mutate(method = "Pseudobulk: all genes")
dn <- bind_rows(dens, allg) %>% filter(method %in% dlv) %>%
  mutate(method = factor(method, dlv, labels = c("Pseudobulk\n(all genes)", "popAlign\nprograms", "D-SPIN\nprograms", "pySCENIC\nregulons", "MAST DEGs\n(selected on these animals)")))
nul <- dens %>% filter(str_detect(method, "^Null")) %>% dplyr::select(-method)
col6b <- setNames(unname(col6[dlv]), levels(dn$method))
p6a <- ggplot(dn, aes(abs_g, density)) +
  geom_area(data = nul, fill = "#EFEFEF", colour = "grey45", linetype = 2, linewidth = 0.35) +
  geom_line(aes(colour = method), linewidth = 0.9) +
  geom_vline(xintercept = g_crit, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ method, nrow = 1) +
  scale_colour_manual(values = col6b, guide = "none") + scale_x_continuous(limits = c(0, 3), breaks = 0:3) +
  labs(x = "|Hedges g|, WildR vs SPF (6 vs 6 animals)", y = "Density",
       title = "Every method's animal-level effect sizes match the no-effect distribution",
       subtitle = glue("Grey: no-effect distribution for 6 vs 6. Line: |g| = {round(g_crit, 2)} (p = 0.05). MAST DEGs are shifted only because they were chosen using these animals")) +
  theme_int(12) + theme(strip.text = element_text(size = 10.5))
nlv <- nl %>% filter(method %in% valid) %>% mutate(method = factor(method, valid, labels = c("Pseudobulk (all genes)", "popAlign programs", "D-SPIN programs", "pySCENIC regulons")))
obs <- nlv %>% filter(observed) %>% left_join(es %>% filter(method %in% valid) %>%
                                                mutate(method = factor(method, valid, labels = levels(nlv$method))) %>% dplyr::select(method, tissue, p_perm), by = c("method", "tissue"))
p6b <- ggplot(nlv, aes(mean_g2)) +
  geom_histogram(aes(fill = method), bins = 40, colour = NA, alpha = 0.55) +
  geom_vline(data = obs, aes(xintercept = mean_g2, colour = method), linewidth = 0.9) +
  geom_text(data = obs, aes(x = Inf, y = Inf, label = glue("observed: {p_lab(p_perm)}")), hjust = 1.05, vjust = 1.4, size = 3, colour = "#333333") +
  facet_grid(tissue ~ method, scales = "free") +
  scale_fill_manual(values = setNames(unname(col6[valid]), levels(nlv$method)), guide = "none") +
  scale_colour_manual(values = setNames(unname(col6[valid]), levels(nlv$method)), guide = "none") +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 3) +
  labs(x = "Mean g² over all tests, for each of the 924 ways to split the 12 animals 6 / 6", y = "Splits",
       title = "The real WildR / SPF split separates animals no more than a random split",
       subtitle = "p = share of random splits with at least the observed mean g² (exact; one tissue = one row)") +
  theme_int(12) + theme(strip.text = element_text(size = 10.5), strip.text.y = element_text(angle = 0))
pwc <- rn("integrated_power_curves.csv") %>% mutate(alpha = factor(if_else(alpha == 0.05, "α = 0.05 (one test)", "α = 10⁻⁴ (genome-wide scale)"),
                                                                    c("α = 0.05 (one test)", "α = 10⁻⁴ (genome-wide scale)")),
                                                     g = factor(glue("g = {g}"), glue("g = {c(0.5, 0.8, 1.2, 1.8)}")))
p6c <- ggplot(pwc, aes(n, power, colour = g, linetype = alpha)) +
  geom_hline(yintercept = 0.8, colour = "grey65", linewidth = 0.3) +
  geom_vline(xintercept = 6, linetype = 2, colour = ORANGE) +
  annotate("label", x = 7, y = 0.55, label = "this study:\n6 / group", hjust = 0, size = 3.2, colour = ORANGE, label.size = 0, fill = "white") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("#C9D3DB", "#8FA3B3", "#4F6F87", SLATE), name = NULL) +
  scale_linetype_manual(values = c(1, 3), name = NULL) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) + scale_x_continuous(breaks = c(6, 20, 40, 60, 80, 100)) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(x = "Animals per group", y = "Power", title = "Animals needed per effect size",
       subtitle = "Two-sample t-test on per-animal values") + theme_int()
save_fig(p6a / p6b / p6c + plot_layout(heights = c(0.75, 1, 0.95)) + plot_annotation(tag_levels = "A"), "integrated_fig6_effect_sizes_power.png", 15)

# ════════════════════════════════════════════════════════════════════════════
# Fig 8 · hidden animal structure: which splits of the animals carry the most signal?
# ════════════════════════════════════════════════════════════════════════════
# Top 5% of the 924 splits by pseudobulk mean g² (all genes): how often is each pair of animals on
# the same side? Clusters = animal groups that differ systematically, whatever the microbiome.
coassign <- function(tis) {
  d <- nl %>% filter(method == "Pseudobulk: all genes", tissue == tis) %>% arrange(desc(mean_g2))
  top <- d %>% slice_head(n = ceiling(0.05 * nrow(d)))
  an <- sort(unique(unlist(str_split(d$lab[1:2], ","))))
  M <- matrix(0, length(an), length(an), dimnames = list(an, an))
  for (l in top$lab) { s <- str_split(l, ",")[[1]]; side <- an %in% s; M <- M + outer(side, side, "==") }
  M <- M / nrow(top); hc <- hclust(as.dist(1 - M), "average")
  as_tibble(M, rownames = "a1") %>% pivot_longer(-a1, names_to = "a2", values_to = "same") %>%
    mutate(tissue = tis, a1 = factor(a1, an[hc$order]), a2 = factor(a2, an[hc$order]), obs_rank = which(d$observed))
}
co <- bind_rows(coassign("AMY"), coassign("HYP")) %>%
  mutate(lab1 = str_remove(a1, "_(AMY|HYP)$"), lab2 = str_remove(a2, "_(AMY|HYP)$"))
p8a <- ggplot(co, aes(a2, fct_rev(a1), fill = same)) + geom_tile(colour = "white", linewidth = 0.4) +
  facet_wrap(~ tissue, scales = "free") +
  scale_x_discrete(labels = function(x) str_remove(x, "_(AMY|HYP)$")) + scale_y_discrete(labels = function(x) str_remove(x, "_(AMY|HYP)$")) +
  scale_fill_gradient(low = "#F4F6F8", high = SLATE, labels = percent, name = "same side in\ntop 5% splits") +
  labs(x = NULL, y = NULL, title = "Animals fall into groups that cut across WildR and SPF",
       subtitle = "Of the 46 random 6 / 6 splits with the largest expression differences, how often each pair of animals is on the same side") +
  theme_int(12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "right")
ag8 <- nl %>% filter(method %in% valid) %>% dplyr::select(method, tissue, lab, mean_g2) %>%
  pivot_wider(names_from = method, values_from = mean_g2) %>% group_by(tissue) %>%
  summarise(across(c(`popAlign programs`, `D-SPIN programs`, `pySCENIC regulons`), ~ cor(.x, `Pseudobulk: all genes`, method = "spearman")), .groups = "drop") %>%
  pivot_longer(-tissue, names_to = "method", values_to = "rho") %>% mutate(method = factor(method, rev(valid[-1])))
p8b <- ggplot(ag8, aes(rho, method, fill = method)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), hjust = -0.2, size = 3.4, colour = "#333333") +
  facet_wrap(~ tissue) + scale_fill_manual(values = col6, guide = "none") + scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Spearman ρ with pseudobulk across the 924 splits", y = NULL, title = "Every method sees the same hidden structure",
       subtitle = "Splits that separate animals in expression also separate them in programs and regulons") + theme_int()
save_fig(p8a / p8b + plot_layout(heights = c(1.5, 0.6)) + plot_annotation(tag_levels = "A"), "integrated_fig8_hidden_animal_structure.png", 10.5)

# ════════════════════════════════════════════════════════════════════════════
# Fig 7 · do methods agree on WHICH supertypes carry signal?
# ════════════════════════════════════════════════════════════════════════════
cc <- rn("integrated_unit_concordance.csv") %>% mutate(across(c(m1, m2), method_of))
ord <- c("popAlign", "Pseudobulk", "D-SPIN", "pySCENIC", "MAST")
cm <- bind_rows(cc, cc %>% rename(m1 = m2, m2 = m1)) %>%
  mutate(m1 = factor(m1, ord), m2 = factor(m2, ord)) %>% filter(as.integer(m1) > as.integer(m2)) %>%
  mutate(lab = glue("{sprintf('%.2f', rho_partial)}{case_when(p_partial < 0.001 ~ '***', p_partial < 0.01 ~ '**', p_partial < 0.05 ~ '*', TRUE ~ '')}"))
p7a <- ggplot(cm, aes(m2, fct_rev(m1), fill = rho_partial)) + geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = lab), size = 4, colour = "#0A0A0A") +
  scale_fill_gradient2(low = BLUE, mid = "#F0EFEC", high = ORANGE, midpoint = 0, limits = c(-0.5, 0.5), oob = squish, name = "partial ρ") +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "Same units, same signal?",
       subtitle = "Spearman ρ between supertypes' hit rates, adjusted\nfor nuclei and tissue. * p < 0.05, ** < 0.01, *** < 0.001") +
  theme_int() + theme(panel.grid = element_blank(), legend.position = "right")
sz <- rn("integrated_unit_size_dependence.csv") %>%
  mutate(method = method_of(modality), method = fct_reorder(method, rho_size), repl = if_else(inference == "nuclei", "nuclei", "animals"))
p7b <- ggplot(sz, aes(rho_size, method, fill = method)) + geom_col(width = 0.6) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_text(aes(x = rho_size + if_else(rho_size > 0, 0.04, -0.04), label = sprintf("%.2f", rho_size), hjust = if_else(rho_size > 0, 0, 1)), size = 3.4, colour = "#333333") +
  scale_fill_manual(values = method_cols, guide = "none") + scale_x_continuous(limits = c(-0.9, 0.95)) +
  labs(x = "Spearman ρ, hit rate vs log nuclei per supertype", y = NULL, title = "Hit rate tracks unit size",
       subtitle = "Nuclei-level MAST finds more in big units;\nanimal-level methods find more in small ones") + theme_int()
save_fig(p7a + p7b + plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = "A"), "integrated_fig7_unit_concordance.png", 5.6)
log_msg("Figures written to ", fig_dir)
