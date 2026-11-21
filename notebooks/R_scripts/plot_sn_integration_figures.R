# plot_sn_integration_figures.R
# Figures for notebooks/0X_Integrated-Results_SN.qmd, drawn from the saved
# harmonized tables written by harmonize_sn_results.R (no recomputation).
#
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_sn_integration_figures.R [run_tag]

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
fig_dir <- file.path(wkdir, "figures/SN_Integration")
tbl_dir <- file.path(wkdir, "tables/SN_Integration")

method_levels <- c("MAST", "pySCENIC", "popAlign", "DSPIN")
tissue_cols <- c(AMY = "#2a78d6", HYP = "#eb6834")
ink <- "#0A0A0A"; ink2 <- "#555555"; hair <- "#E5E5E5"

theme_sn <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = hair, linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = "bold"),
          plot.title = element_text(colour = ink, face = "bold"),
          legend.position = "top", plot.background = element_rect(fill = "white", colour = NA))
}

# Panel A: strongest adjusted p per cell-type unit, by method
plot_min_padj <- function(summ, sig_q = 0.1) {
  df <- summ %>%
    select(unit_key, tissue, level, ends_with("_min_padj")) %>%
    pivot_longer(ends_with("_min_padj"), names_to = "method", values_to = "min_padj") %>%
    mutate(method = factor(str_remove(method, "_min_padj"), levels = rev(method_levels))) %>%
    filter(!is.na(min_padj)) %>%
    mutate(y = pmin(-log10(min_padj), 20))
  ggplot(df, aes(y, method, colour = tissue)) +
    geom_vline(xintercept = -log10(sig_q), linetype = "dashed", colour = ink2, linewidth = 0.4) +
    geom_jitter(height = 0.22, width = 0, size = 1.6, alpha = 0.75, stroke = 0) +
    scale_colour_manual(values = tissue_cols, name = NULL) +
    labs(x = expression(-log[10]~"(smallest adjusted p in unit), capped at 20"), y = NULL,
         title = "Strongest WildR vs SPF signal per cell-type unit") +
    theme_sn()
}

# Panel B: class rollup tile (every hit assigned to its parent Allen class)
plot_class_rollup <- function(roll) {
  grid <- roll %>% distinct(tissue, class_label, class_name) %>%
    crossing(method = method_levels) %>%
    left_join(roll, by = c("tissue", "class_label", "class_name", "method")) %>%
    mutate(method = factor(method, levels = method_levels),
           state = case_when(is.na(n_tested) ~ "not run", n_sig == 0 ~ "0", TRUE ~ "hit"),
           fill_val = if_else(state == "hit", log10(n_sig + 1), NA_real_),
           lab = case_when(state == "not run" ~ "–", TRUE ~ as.character(n_sig)),
           row = glue("{tissue} · {class_name}"))
  ggplot(grid, aes(method, fct_rev(row))) +
    geom_tile(aes(fill = fill_val), colour = "white", linewidth = 1) +
    geom_tile(data = filter(grid, state == "not run"), fill = "#F2F2F2", colour = "white", linewidth = 1) +
    geom_text(aes(label = lab, colour = fill_val > 2.3), size = 2.9, show.legend = FALSE) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", na.value = "white",
                        name = "log10(n significant + 1)") +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white"), na.value = ink2) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = "Significant features rolled up to Allen class") +
    theme_sn() + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8))
}

run_plot_sn_figures <- function(run_tag) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  summ <- read_csv(file.path(tbl_dir, glue("SN_Integration_celltype_summary_{run_tag}.csv")), show_col_types = FALSE)
  roll <- read_csv(file.path(tbl_dir, glue("SN_Integration_class_rollup_{run_tag}.csv")), show_col_types = FALSE)
  ggsave(file.path(fig_dir, glue("SN_Integration_min-padj_by_method_{run_tag}.png")),
         plot_min_padj(summ), width = 8, height = 4.5, dpi = 200, bg = "white")
  ggsave(file.path(fig_dir, glue("SN_Integration_class-rollup_tile_{run_tag}.png")),
         plot_class_rollup(roll), width = 7, height = 9.5, dpi = 200, bg = "white")
  invisible(TRUE)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  run_plot_sn_figures(if (length(args)) args[1] else as.character(Sys.Date()))
}
