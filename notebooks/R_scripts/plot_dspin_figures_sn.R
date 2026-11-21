#!/usr/bin/env Rscript
# Report figures for 0X_DSPIN_Analysis.qmd, drawn from the animal-level summary written by
# summarize_dspin_responses_sn.R. Light (reads CSVs; ~30 s). Never touches fits.
#
# Usage
#   Rscript notebooks/R_scripts/plot_dspin_figures_sn.R --tag 2026-10-06
#
# Writes figures/DSPIN/snRNASeq/summary_<tag>/DSPIN_*.png and two small display tables
# (MAST concordance, curated DAVID term frequencies) into tables/DSPIN/summary_<tag>/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(patchwork)
  library(ggrepel)
  library(scico)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1]]
}
tag <- get_arg("--tag", "2026-10-06")

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
dspin_dir <- file.path(wkdir, "data/interim/dspin")
sum_dir <- file.path(dspin_dir, glue("summary_{tag}"))
tbl_dir <- file.path(wkdir, "tables/DSPIN", glue("summary_{tag}"))
fig_dir <- file.path(wkdir, "figures/DSPIN/snRNASeq", glue("summary_{tag}"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prog <- read_csv(file.path(tbl_dir, "DSPIN_program_tests.csv"), show_col_types = FALSE)
unit_t <- read_csv(file.path(tbl_dir, "DSPIN_unit_tests.csv"), show_col_types = FALSE)
cal <- read_csv(file.path(tbl_dir, "DSPIN_calibration.csv"), show_col_types = FALSE)
animals <- read_csv(file.path(sum_dir, glue("dspin_animal_relative_h_{tag}.csv")), show_col_types = FALSE)
david_all <- read_csv(file.path(sum_dir, glue("dspin_david_all_terms_{tag}.csv")), show_col_types = FALSE)

# ── house style ──────────────────────────────────────────────────────────────
mb_cols <- c(SPF = "#27455C", WildR = "#FF6C0C")
tis_cols <- c(AMY = "#3B7EA1", HYP = "#A36AAB")
theme_dspin <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#EDEDED"),
      panel.border = element_rect(colour = "#D0D0D0"),
      strip.background = element_rect(fill = "#FFF3EB", colour = NA),
      strip.text = element_text(face = "bold", size = rel(0.85), colour = "#222222"),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(colour = "#555555", size = rel(0.85)),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      axis.title = element_text(size = rel(0.9)),
      plot.tag = element_text(face = "bold", size = rel(1.2))
    )
}
theme_set(theme_dspin())
save_png <- function(p, name, w, h) {
  f <- file.path(fig_dir, name)
  ggsave(f, p, width = w, height = h, dpi = 200, bg = "white")
  message("wrote ", f)
}
short_unit <- function(x, n = 30) str_trunc(str_remove(x, "^\\d+ "), n)
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

unit_t <- unit_t %>% mutate(unit_short = paste(tissue, short_unit(supertype_name, 34), sep = " · "))
g <- read_csv(file.path(tbl_dir, "DSPIN_global_summary.csv"), show_col_types = FALSE)

# ── Fig 1. Cohort: nuclei per unit and per animal ────────────────────────────
f1a <- unit_t %>%
  ggplot(aes(n_cells, median_cells_per_animal, colour = tissue)) +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = tis_cols) +
  labs(x = "Nuclei in unit (log scale)", y = "Median nuclei per animal (log scale)", colour = NULL,
       title = "140 supertype × tissue units", subtitle = "Each unit has 6 SPF + 6 WildR animals")
per_animal <- animals %>%
  distinct(unit_id, sample_id, microbiome, n_cells) %>%
  mutate(tissue = str_extract(sample_id, "(AMY|HYP)$"),
         animal = str_remove(sample_id, "_(AMY|HYP)$")) %>%
  group_by(tissue, animal, microbiome) %>%
  summarise(nuclei = sum(n_cells), .groups = "drop") %>%
  mutate(animal = fct_reorder(animal, nuclei))
f1b <- per_animal %>%
  ggplot(aes(nuclei, animal, fill = microbiome)) +
  geom_col(width = 0.7) +
  facet_wrap(~tissue, scales = "free_y") +
  scale_fill_manual(values = mb_cols) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Nuclei across all fitted units", y = NULL, fill = NULL,
       title = "Nuclei contributed per animal", subtitle = "Animal = DSPIN condition (sample_id)")
save_png(f1a + f1b + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A"),
         "DSPIN_cohort_overview.png", 13, 5.6)

# ── Fig 2. Calibration: is there more signal than the permutation null? ──────
p_prog <- prog %>% mutate(bin = cut(p_perm, breaks = seq(0, 1, 0.05), include.lowest = TRUE))
exp_per_bin <- nrow(prog) * 0.05
f2a <- ggplot(prog, aes(p_perm)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "#9DB3C4", colour = "white") +
  geom_hline(yintercept = exp_per_bin, linetype = "dashed", colour = "#FF6C0C", linewidth = 0.8) +
  annotate("label", x = 0.5, y = exp_per_bin * 1.22, label = "dashed = expected if WildR had no effect",
           colour = "#FF6C0C", size = 3.8, label.size = 0, fill = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = 0.01)) +
  labs(x = "Exact permutation p (WildR vs SPF, 924 relabelings)", y = "Programs",
       title = glue("Program tests: {scales::comma(nrow(prog))} programs × 12 animals"),
       subtitle = "A real microbiome effect would pile up on the left")
qq_df <- function(p) {
  tibble(obs = sort(p)) %>% mutate(exp = ppoints(n()), lo = qbeta(0.025, row_number(), n() - row_number() + 1),
                                    hi = qbeta(0.975, row_number(), n() - row_number() + 1))
}
qq_plot <- function(df, title, subtitle, floor_p = NULL) {
  pl <- ggplot(df, aes(-log10(exp), -log10(obs))) +
    geom_ribbon(aes(ymin = -log10(hi), ymax = -log10(lo)), fill = "#EDEDED") +
    geom_abline(colour = "#767676", linetype = "dashed") +
    geom_point(size = 1.6, colour = "#27455C", alpha = 0.8) +
    coord_equal() +
    labs(x = "Expected −log10 p", y = "Observed −log10 p", title = title, subtitle = subtitle)
  if (!is.null(floor_p)) pl <- pl + geom_hline(yintercept = -log10(floor_p), colour = "#FF6C0C", linetype = "dotted")
  pl
}
f2b <- qq_plot(qq_df(prog$p_perm), "Program-level QQ", "Dotted = smallest attainable p (2/924)", 2 / 924)
f2c <- qq_plot(qq_df(unit_t$p_unit), "Unit-level omnibus QQ", "Σ t² over 20 programs, permuted")
save_png((f2a | f2b | f2c) + plot_layout(widths = c(1.5, 1, 1)) + plot_annotation(tag_levels = "A"),
         "DSPIN_calibration.png", 15, 5.4)

# ── Fig 3. Noise floor: |h| tracks sampling depth, not microbiome ────────────
null_ratio <- sqrt(7 * 5) / 6   # E[ratio_loo] under no effect, n = 6 SPF
f3a <- unit_t %>%
  ggplot(aes(mean_abs_h_SPF_loo, mean_abs_h_WildR, colour = tissue)) +
  geom_abline(slope = 1, colour = "#767676", linetype = "dashed") +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_colour_manual(values = tis_cols) +
  coord_equal() +
  labs(x = "SPF animals: mean |h| vs other SPF", y = "WildR animals: mean |h| vs SPF",
       colour = NULL, title = "WildR animals are no further from SPF\nthan SPF animals are from each other",
       subtitle = glue("WildR > SPF in {round(100 * g$frac_units_WildR_gt_SPF_loo)}% of units; ",
                       "median ratio {round(g$median_abs_h_ratio_loo, 2)} (null ≈ {round(null_ratio, 2)})"))
old_top <- unit_t %>% arrange(desc(mean_abs_h_WildR)) %>% slice_head(n = 6)
f3b <- unit_t %>%
  ggplot(aes(median_cells_per_animal, mean_abs_h_WildR, colour = tissue)) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, colour = "#767676", linewidth = 0.7, formula = y ~ x) +
  geom_text_repel(data = old_top, aes(label = short_unit(supertype_name, 24)), size = 3.3,
                  colour = "#222222", min.segment.length = 0, box.padding = 0.5, seed = 1) +
  scale_x_log10() +
  scale_colour_manual(values = tis_cols) +
  labs(x = "Median nuclei per animal (log scale)", y = "WildR mean |relative h|", colour = NULL,
       title = "The old 'strongest response' ranking is sampling noise",
       subtitle = glue("Spearman ρ = {round(g$spearman_rho_absh_vs_cells, 2)}; labelled = previous top-6 units"))
an <- animals %>%
  mutate(tissue = str_extract(sample_id, "(AMY|HYP)$"), animal = str_remove(sample_id, "_(AMY|HYP)$")) %>%
  group_by(tissue, animal, microbiome) %>%
  summarise(mean_abs_h = mean(abs(relative_h)), med_cells = median(n_cells), .groups = "drop")
f3c <- an %>%
  ggplot(aes(med_cells, mean_abs_h, colour = microbiome)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = animal), size = 3, show.legend = FALSE, seed = 2, max.overlaps = 30) +
  facet_wrap(~tissue, scales = "free_x") +
  scale_colour_manual(values = mb_cols) +
  labs(x = "Median nuclei per unit for this animal", y = "Mean |relative h|", colour = NULL,
       title = "Per-animal deviation", subtitle = "Animals with fewer nuclei deviate more; groups interleave")
save_png((f3a | f3b) / f3c + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "A"),
         "DSPIN_noise_floor.png", 14, 11.5)

# ── Fig 4. Unit- and program-level volcanoes ─────────────────────────────────
f4a <- unit_t %>%
  ggplot(aes(rms_hedges_g, -log10(p_unit), colour = tissue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#A5A5A5") +
  geom_point(size = 2.4, alpha = 0.85) +
  geom_text_repel(data = slice_min(unit_t, p_unit, n = 8, with_ties = FALSE),
                  aes(label = short_unit(supertype_name, 26)), size = 3.3, colour = "#222222",
                  min.segment.length = 0, seed = 3, box.padding = 0.4) +
  scale_colour_manual(values = tis_cols) +
  labs(x = "RMS Hedges' g across the unit's 20 programs", y = "−log10 omnibus permutation p", colour = NULL,
       title = "Units", subtitle = glue("{g$n_unit_p_lt_0.05} / {g$n_units} units p < 0.05 (≈ {round(0.05 * g$n_units)} expected by chance); none q < 0.1"))
leads <- prog %>% filter(lead)
f4b <- prog %>%
  mutate(class = if_else(lead, "lead", "other")) %>%
  ggplot(aes(hedges_g, -log10(p_perm))) +
  geom_jitter(data = ~ filter(.x, !lead), colour = "#BFC9D1", size = 1.1, alpha = 0.6, width = 0, height = 0.03) +
  geom_point(data = leads, aes(colour = tissue), size = 2.8) +
  geom_text_repel(data = leads, aes(label = paste0(short_unit(supertype_name, 20), " ", program_id), colour = tissue),
                  size = 3, min.segment.length = 0, seed = 4, max.overlaps = 40, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.0065), linetype = "dashed", colour = "#A5A5A5") +
  scale_colour_manual(values = tis_cols) +
  labs(x = "Hedges' g (WildR − SPF relative h)", y = "−log10 permutation p (jittered)", colour = NULL,
       title = "Programs", subtitle = glue("{g$n_leads} leads (p ≤ 6/924, ≥ 4/6 WildR beyond SPF range) vs ≈ {round(cal$expected[cal$level == 'program'][2])} expected by chance"))
save_png(f4a + f4b + plot_annotation(tag_levels = "A"), "DSPIN_volcanoes.png", 15, 6.4)

# ── Fig 5. Lead programs, animal by animal ───────────────────────────────────
lead_lab <- leads %>%
  arrange(p_perm, desc(abs(hedges_g))) %>%
  mutate(facet = glue("{tissue} · {short_unit(supertype_name, 26)}\n{program_id}: {word(top_genes, 1, 3, sep = fixed(', '))}\np = {fmt_p(p_perm)}, g = {round(hedges_g, 1)}") %>% as.character())
lead_an <- animals %>%
  inner_join(lead_lab %>% select(unit_id, program_id, facet), by = c("unit_id", "program_id")) %>%
  mutate(facet = factor(facet, levels = lead_lab$facet), microbiome = factor(microbiome, c("SPF", "WildR")))
f5 <- lead_an %>%
  ggplot(aes(microbiome, relative_h, colour = microbiome)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4, colour = "#555555") +
  geom_beeswarm(size = 2.6, cex = 3) +
  facet_wrap(~facet, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = mb_cols, guide = "none") +
  labs(x = NULL, y = "Relative program response h (vs SPF mean)",
       title = "The 13 pre-specified lead programs, one point per animal",
       subtitle = "Bar = group mean. Nominal only: the count matches chance (Figure 2), so treat these as hypotheses, not findings.") +
  theme(strip.text = element_text(size = 9, face = "plain", hjust = 0, lineheight = 1.05))
save_png(f5, "DSPIN_lead_programs_animals.png", 15, 10)

# ── Fig 6. Program heatmaps for the 4 best-scoring units ─────────────────────
top_units <- unit_t %>% slice_min(p_unit, n = 4, with_ties = FALSE)
heat_one <- function(uid) {
  u <- top_units %>% filter(unit_id == uid)
  pr <- prog %>% filter(unit_id == uid) %>% arrange(delta_h) %>%
    mutate(row = glue("{program_id} · {word(top_genes, 1, 3, sep = fixed(', '))}{if_else(lead, '  ★', '')}") %>% as.character())
  d <- animals %>% filter(unit_id == uid) %>%
    inner_join(pr %>% select(program_id, row), by = "program_id") %>%
    mutate(row = factor(row, levels = pr$row),
           animal = str_remove(sample_id, "_(AMY|HYP)$"),
           microbiome = factor(microbiome, c("SPF", "WildR")))
  lim <- quantile(abs(d$relative_h), 0.98)
  ggplot(d, aes(animal, row, fill = pmax(pmin(relative_h, lim), -lim))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(. ~ microbiome, scales = "free_x", space = "free_x") +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim), name = "relative h") +
    labs(x = NULL, y = NULL, title = u$unit_short,
         subtitle = glue("omnibus p = {fmt_p(u$p_unit)} · rows ordered by WildR − SPF")) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9), axis.text.y = element_text(size = 9.5),
          panel.grid = element_blank(), legend.position = "right",
          legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.35, "cm"))
}
f6 <- wrap_plots(map(top_units$unit_id, heat_one), ncol = 2) + plot_annotation(tag_levels = "A")
save_png(f6, "DSPIN_heatmaps_top_units.png", 16, 13)

# ── Fig 7. Coupling networks J for the same units ────────────────────────────
net_one <- function(uid) {
  u <- top_units %>% filter(unit_id == uid)
  J <- read_csv(file.path(dspin_dir, uid, "network.csv"), show_col_types = FALSE)
  M <- as.matrix(J[, -1]); rownames(M) <- J[[1]]
  diag(M) <- NA
  M0 <- M; M0[is.na(M0)] <- 0
  ord <- rownames(M)[hclust(as.dist(1 - cor(M0)), method = "average")$order]
  pr <- prog %>% filter(unit_id == uid) %>%
    mutate(lab = glue("{program_id} · {word(top_genes, 1, 2, sep = fixed(', '))}") %>% as.character())
  lab <- setNames(pr$lab, pr$program_id)
  d <- as_tibble(M, rownames = "i") %>% pivot_longer(-i, names_to = "j", values_to = "J") %>%
    mutate(i = factor(i, ord), j = factor(j, rev(ord)))
  lim <- quantile(abs(d$J), 0.98, na.rm = TRUE)
  ggplot(d, aes(i, j, fill = pmax(pmin(J, lim), -lim))) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim), na.value = "#F4F4F4", name = "J") +
    scale_x_discrete(labels = function(x) x) +
    scale_y_discrete(labels = function(x) lab[x]) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = u$unit_short, subtitle = "Off-diagonal couplings; clustered order") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8.5), axis.text.y = element_text(size = 9),
          panel.grid = element_blank(), legend.position = "right",
          legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.35, "cm"))
}
f7 <- wrap_plots(map(top_units$unit_id, net_one), ncol = 2) + plot_annotation(tag_levels = "A")
save_png(f7, "DSPIN_networks_top_units.png", 15, 12.5)

# ── Fig 8. MAST concordance + curated DAVID terms ────────────────────────────
mc <- prog %>% filter(n_mast_hits >= 3)
ct <- cor.test(mc$delta_h, mc$mast_mean_coef, method = "spearman", exact = FALSE)
rho_t <- mc %>% group_by(tissue) %>% summarise(rho = cor(delta_h, mast_mean_coef, method = "spearman"), n = n(), .groups = "drop")
f8a <- mc %>%
  ggplot(aes(delta_h, mast_mean_coef, colour = tissue, size = n_mast_hits)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD") + geom_vline(xintercept = 0, colour = "#BDBDBD") +
  geom_point(alpha = 0.85, stroke = 0) +
  geom_smooth(aes(group = 1, size = NULL), method = "lm", se = TRUE, colour = "#555555", fill = "#EDEDED",
              linewidth = 0.7, show.legend = FALSE, formula = y ~ x) +
  scale_colour_manual(values = tis_cols) +
  scale_size_area(max_size = 6, name = "MAST hits\nin program") +
  labs(x = "D-SPIN Δh (WildR − SPF)", y = "Mean MAST coef of member hit genes", colour = NULL,
       title = "Where MAST finds DE genes, D-SPIN h tends to move the same way",
       subtitle = glue("{nrow(mc)} programs with ≥ 3 MAST hits (padj < 0.1); Spearman ρ = {round(ct$estimate, 2)}, p = {fmt_p(ct$p.value)}\n",
                       "AMY ρ = {round(rho_t$rho[rho_t$tissue == 'AMY'], 2)} (n = {rho_t$n[rho_t$tissue == 'AMY']}), ",
                       "HYP ρ = {round(rho_t$rho[rho_t$tissue == 'HYP'], 2)} (n = {rho_t$n[rho_t$tissue == 'HYP']}); same sign in {round(100 * mean(sign(mc$delta_h) == sign(mc$mast_mean_coef)))}%")) +
  theme(legend.position = "right")
mast_conc <- tibble(n_programs = nrow(mc), spearman_rho = unname(ct$estimate), p = ct$p.value,
                    rho_AMY = rho_t$rho[rho_t$tissue == "AMY"], rho_HYP = rho_t$rho[rho_t$tissue == "HYP"],
                    frac_sign_agree = mean(sign(mc$delta_h) == sign(mc$mast_mean_coef)))
write_csv(mast_conc, file.path(tbl_dir, "DSPIN_mast_concordance.csv"))
write_csv(mc %>% arrange(desc(n_mast_hits)) %>%
            select(unit_short, program_id, top_genes, delta_h, p_perm, n_mast_hits, mast_mean_coef, mast_genes),
          file.path(tbl_dir, "DSPIN_mast_overlap_programs.csv"))

curated <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "KEGG_PATHWAY")
term_freq <- david_all %>%
  filter(category %in% curated, benjamini < 0.05, count >= 3) %>%
  mutate(term = str_remove(term, "^GO:\\d+~") %>% str_remove("^mmu\\d+:"),
         category = recode(category, GOTERM_BP_DIRECT = "GO BP", GOTERM_MF_DIRECT = "GO MF",
                           GOTERM_CC_DIRECT = "GO CC", KEGG_PATHWAY = "KEGG")) %>%
  distinct(unit_id, program_id, category, term, tissue) %>%
  count(category, term, tissue, name = "n_programs") %>%
  group_by(category, term) %>% mutate(total = sum(n_programs)) %>% ungroup()
write_csv(term_freq %>% arrange(desc(total)), file.path(tbl_dir, "DSPIN_david_curated_term_frequency.csv"))
top_terms <- term_freq %>% distinct(category, term, total) %>%
  filter(category != "GO CC") %>% slice_max(total, n = 22, with_ties = FALSE)
f8b <- term_freq %>% semi_join(top_terms, by = c("category", "term")) %>%
  mutate(term_lab = fct_reorder(str_trunc(glue("{term}  [{category}]"), 58), total)) %>%
  ggplot(aes(n_programs, term_lab, fill = tissue)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = tis_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Programs with this term (Benjamini < 0.05)", y = NULL, fill = NULL,
       title = "What the programs are made of",
       subtitle = "Most frequent GO BP / MF and KEGG terms across all 2,800 programs (GO CC omitted)") +
  theme(axis.text.y = element_text(size = 10))
save_png(f8a + f8b + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A"),
         "DSPIN_mast_and_annotation.png", 16, 7.2)

message("Done → ", fig_dir)
