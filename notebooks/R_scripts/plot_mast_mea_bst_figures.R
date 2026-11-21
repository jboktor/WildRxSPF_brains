# Paper-figure drafts from the raw single-nucleus MAST tables (fit 2026-01-04).
#
# Three options for the MeA / BST social-behavior thread:
#   1. Scope of testing across every amygdala supertype (all genes, including non-hits).
#   2. Paired volcanoes: MEA Slc17a7 Glut_3 vs MEA-BST Sox6 Gaba_2.
#   3. Coefficient intervals for genes that pass BH 0.1 in all three MEA Slc17a7 supertypes.
#
# Run:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_mast_mea_bst_figures.R
#
# Writes new dated files under figures/MAST_Raw_SN and tables/MAST_Raw_SN.
# Refuses to overwrite a path that already exists.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(data.table)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
mast_dir <- file.path(wkdir, "data/interim/mast/snRNASeq/supertype_label")
fig_dir <- file.path(wkdir, "figures/MAST_Raw_SN")
tbl_dir <- file.path(wkdir, "tables/MAST_Raw_SN")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- "2026-11-06"
fig_stamp <- "2026-11-06_r2"
padj_cut <- 0.1
raw_date <- "2026-01-04"

bh_orange <- "#FF6C0C"
bh_orange_lt <- "#FF9240"
bh_orange_dk <- "#9C3D00"
bh_slate <- "#27455C"
bh_gray <- "#767676"
bh_ns <- "#D4D4D4"
bh_ink <- "#0A0A0A"

theme_paper <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "#E5E5E5"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = base + 1, color = bh_ink),
      plot.subtitle = element_text(size = base - 1, color = "#555555"),
      axis.title = element_text(color = "#555555"),
      axis.text = element_text(color = bh_ink),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", color = bh_ink)
    )
}

write_new <- function(path, writer) {
  if (file.exists(path)) {
    message("Already exists, left untouched: ", path)
    return(invisible(FALSE))
  }
  writer(path)
  message("Wrote ", path)
  invisible(TRUE)
}

save_plot <- function(stem, plot, width, height) {
  png_path <- file.path(fig_dir, glue("{stem}_{fig_stamp}.png"))
  pdf_path <- file.path(fig_dir, glue("{stem}_{fig_stamp}.pdf"))
  write_new(png_path, function(p) ggsave(p, plot, width = width, height = height, dpi = 300, bg = "white"))
  write_new(pdf_path, function(p) ggsave(p, plot, width = width, height = height, device = cairo_pdf, bg = "white"))
}

load_unit <- function(label, tissue = "AMY") {
  f <- file.path(mast_dir, glue("{label}_{tissue}_{raw_date}.rds"))
  if (!file.exists(f)) stop("Missing raw table: ", f)
  as_tibble(readRDS(f)) %>%
    mutate(
      gene = sub("-GRCm39$", "", primerid),
      neglog10 = -log10(pmax(p_adjust, 1e-300)),
      label = label,
      tissue = tissue
    )
}

# --- names for every amygdala supertype raw table --------------------------------
hits <- read_csv(
  file.path(wkdir, "data/interim/mast/snRNASeq/mast_signif_hits_2026-08-17_p_adjust_0.1.csv"),
  show_col_types = FALSE
)
hit_names <- hits %>%
  filter(cell_type_level == "supertype") %>%
  distinct(label = cell_type, tissue, cell_type_name)

manifest <- read_csv(
  file.path(wkdir, "data/interim/popalign/snRNASeq/manifest_2026-10-04.csv"),
  show_col_types = FALSE
) %>%
  filter(level == "supertype") %>%
  distinct(label, tissue, manifest_name = name)

files <- list.files(mast_dir, pattern = glue("_AMY_{raw_date}\\.rds$"), full.names = TRUE)
message("Summarizing ", length(files), " amygdala supertype tables")

scope <- rbindlist(lapply(files, function(f) {
  x <- as.data.table(readRDS(f))
  data.table(
    label = str_match(basename(f), "^(.*)_AMY_\\d{4}-\\d{2}-\\d{2}\\.rds$")[, 2],
    tissue = "AMY",
    n_genes = nrow(x),
    n_padj_lt_0.1 = sum(x$p_adjust < padj_cut, na.rm = TRUE),
    n_up = sum(x$p_adjust < padj_cut & x$coef > 0, na.rm = TRUE),
    n_down = sum(x$p_adjust < padj_cut & x$coef < 0, na.rm = TRUE)
  )
})) %>%
  as_tibble() %>%
  left_join(hit_names, by = c("label", "tissue")) %>%
  left_join(manifest, by = c("label", "tissue")) %>%
  mutate(
    cell_type_name = coalesce(cell_type_name, manifest_name, label),
    territory = case_when(
      str_detect(cell_type_name, "MEA Slc17a7") ~ "MEA Slc17a7 Glut",
      str_detect(cell_type_name, "MEA-BST|CEA-BST") ~ "MEA/CEA–BST",
      TRUE ~ "Other amygdala supertype"
    ),
    territory = factor(territory, levels = c("MEA Slc17a7 Glut", "MEA/CEA–BST", "Other amygdala supertype"))
  ) %>%
  arrange(desc(n_padj_lt_0.1), cell_type_name) %>%
  mutate(rank = row_number())

scope_path <- file.path(tbl_dir, glue("MAST_AMY_supertype_scope_{stamp}.csv"))
write_new(scope_path, function(p) write_csv(scope %>% select(-manifest_name), p))

# Label the story cell types plus the largest non-story peak so the tail is readable.
label_scope <- scope %>%
  filter(
    territory == "MEA Slc17a7 Glut" |
      (territory == "MEA/CEA–BST" & (n_padj_lt_0.1 >= 10 | str_detect(cell_type_name, "0338|0368"))) |
      str_detect(cell_type_name, "0059 COAp")
  ) %>%
  mutate(
    label_short = str_replace(cell_type_name, "^[0-9]+ ", ""),
    nudge_x = 16,
    nudge_y = case_when(
      str_detect(label_short, "COAp") ~ 20,
      str_detect(label_short, "Glut_3") ~ 55,
      str_detect(label_short, "Glut_1") ~ 30,
      str_detect(label_short, "Glut_2") ~ -20,
      str_detect(label_short, "Cyp26b1") ~ 45,
      str_detect(label_short, "Sox6") ~ 30,
      str_detect(label_short, "Nfib") ~ 8,
      TRUE ~ 0
    )
  )

p_scope <- ggplot(scope, aes(rank, n_padj_lt_0.1, color = territory)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#E5E5E5") +
  geom_point(size = 2.1, alpha = 0.9) +
  geom_text_repel(
    data = label_scope,
    aes(label = label_short),
    size = 2.7,
    color = bh_ink,
    seed = 1,
    nudge_x = label_scope$nudge_x,
    nudge_y = label_scope$nudge_y,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = bh_gray,
    segment.size = 0.3,
    direction = "y",
    hjust = 0,
    xlim = c(6, 78),
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "MEA Slc17a7 Glut" = bh_orange,
    "MEA/CEA–BST" = bh_slate,
    "Other amygdala supertype" = bh_ns
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Where amygdala MAST testing produced BH hits",
    subtitle = glue(
      "{nrow(scope)} supertypes, {format(sum(scope$n_genes), big.mark = ',')} gene tests. ",
      "Positive coefficient = higher in WildR. BH is within each supertype."
    ),
    x = "Amygdala supertype, ranked by BH hits",
    y = "Genes with BH p_adjust < 0.1"
  ) +
  theme_paper() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot("MAST_AMY_supertype_scope", p_scope, width = 8.2, height = 5.2)

# --- three MEA Slc17a7 raw tables and the concordant gene set -------------------
mea_ids <- c(
  "0055 MEA Slc17a7 Glut_1" = "CS20230722_SUPT_0055",
  "0056 MEA Slc17a7 Glut_2" = "CS20230722_SUPT_0056",
  "0057 MEA Slc17a7 Glut_3" = "CS20230722_SUPT_0057"
)
mea <- bind_rows(lapply(names(mea_ids), function(nm) {
  load_unit(mea_ids[[nm]]) %>% mutate(supertype = nm)
}))

shared <- mea %>%
  filter(p_adjust < padj_cut) %>%
  group_by(gene) %>%
  summarise(
    n_supertypes = n(),
    n_up = sum(coef > 0),
    n_down = sum(coef < 0),
    .groups = "drop"
  ) %>%
  filter(n_supertypes == 3)

shared_long <- mea %>%
  inner_join(shared, by = "gene") %>%
  mutate(
    direction = if_else(n_up == 3, "Higher in WildR in all three", "Higher in SPF in all three"),
    gene = reorder(gene, coef)
  )

sig_path <- file.path(tbl_dir, glue("MAST_MEA_Slc17a7_shared_genes_{stamp}.csv"))
write_new(sig_path, function(p) {
  write_csv(
    shared_long %>%
      select(supertype, gene, primerid, coef, ci_lo = ci.lo, ci_hi = ci.hi, hurdle_p = `Pr(>Chisq)`, p_adjust, direction) %>%
      arrange(gene, supertype),
    p
  )
})

# --- paired volcanoes -----------------------------------------------------------
sox <- load_unit("CS20230722_SUPT_0338") %>%
  mutate(supertype = "0338 MEA-BST Sox6 Gaba_2")
glut3 <- mea %>% filter(supertype == "0057 MEA Slc17a7 Glut_3")

volcano_df <- bind_rows(
  glut3 %>% mutate(panel = "MEA Slc17a7 Glut_3", signature = gene %in% shared$gene),
  sox %>% mutate(panel = "MEA-BST Sox6 Gaba_2", signature = FALSE)
) %>%
  mutate(
    point_class = case_when(
      signature ~ "Shared across MEA Slc17a7",
      p_adjust < padj_cut & coef > 0 ~ "BH hit, higher in WildR",
      p_adjust < padj_cut & coef < 0 ~ "BH hit, higher in SPF",
      TRUE ~ "Not significant"
    ),
    point_class = factor(point_class, levels = c(
      "Shared across MEA Slc17a7",
      "BH hit, higher in WildR",
      "BH hit, higher in SPF",
      "Not significant"
    )),
    panel = factor(panel, levels = c("MEA Slc17a7 Glut_3", "MEA-BST Sox6 Gaba_2"))
  )

label_volc <- volcano_df %>%
  filter(signature | (panel == "MEA-BST Sox6 Gaba_2" & p_adjust < padj_cut))

x_lim <- quantile(abs(volcano_df$coef), 0.995, na.rm = TRUE)
y_lim <- max(volcano_df$neglog10, na.rm = TRUE) * 1.08

p_volc <- ggplot(volcano_df, aes(coef, neglog10)) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.35, color = bh_gray) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = bh_gray) +
  geom_point(aes(color = point_class, size = point_class), alpha = 0.85) +
  geom_text_repel(
    data = label_volc,
    aes(label = gene),
    size = 2.6,
    color = bh_ink,
    seed = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.25,
    segment.color = bh_gray,
    segment.size = 0.25,
    show.legend = FALSE
  ) +
  facet_wrap(~panel, nrow = 1) +
  coord_cartesian(xlim = c(-x_lim * 1.35, x_lim * 1.28), ylim = c(0, y_lim), clip = "off") +
  scale_color_manual(values = c(
    "Shared across MEA Slc17a7" = bh_orange_dk,
    "BH hit, higher in WildR" = bh_orange,
    "BH hit, higher in SPF" = bh_slate,
    "Not significant" = bh_ns
  )) +
  scale_size_manual(values = c(
    "Shared across MEA Slc17a7" = 2.4,
    "BH hit, higher in WildR" = 1.5,
    "BH hit, higher in SPF" = 1.5,
    "Not significant" = 0.9
  ), guide = "none") +
  labs(
    title = "Excitatory MeA is transcriptionally busy; MEA–BST Sox6 GABA is not",
    subtitle = glue(
      "Every tested gene. Dashed line: BH {padj_cut} within the supertype. ",
      "Dark points: {nrow(shared)} genes with BH < {padj_cut} in Glut_1, Glut_2, and Glut_3."
    ),
    x = "WildR coefficient (log2 CPM+1)",
    y = expression(-log[10] ~ "BH p_adjust")
  ) +
  theme_paper() +
  theme(
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 16, 8, 12)
  )

save_plot("MAST_MEA_BST_volcano_pair", p_volc, width = 10.4, height = 5.6)

# --- signature intervals --------------------------------------------------------
gene_order <- shared_long %>%
  group_by(gene) %>%
  summarise(m = mean(coef), .groups = "drop") %>%
  arrange(m) %>%
  pull(gene)

p_sig <- shared_long %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    supertype_short = str_extract(supertype, "Glut_[123]")
  ) %>%
  ggplot(aes(coef, gene, color = supertype_short, shape = supertype_short)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = bh_gray) +
  geom_errorbar(
    aes(xmin = ci.lo, xmax = ci.hi),
    orientation = "y",
    width = 0.25,
    linewidth = 0.4,
    position = position_dodge(width = 0.55),
    show.legend = FALSE
  ) +
  geom_point(size = 2.3, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c(
    "Glut_1" = bh_orange_lt,
    "Glut_2" = bh_orange,
    "Glut_3" = bh_orange_dk
  )) +
  labs(
    title = "Genes significant in every MEA Slc17a7 supertype",
    subtitle = "BH p_adjust < 0.1 in Glut_1, Glut_2, and Glut_3. Sign is shared. Bars are the model CI.",
    x = "WildR coefficient (log2 CPM+1)",
    y = NULL
  ) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

save_plot("MAST_MEA_Slc17a7_shared_signature", p_sig, width = 7.2, height = 6.2)

message("Done. Concordant genes: ", paste(sort(as.character(unique(shared_long$gene))), collapse = ", "))
message("Scope units: ", nrow(scope),
        " | MEA Slc17a7 BH genes: ",
        paste(scope$n_padj_lt_0.1[scope$territory == "MEA Slc17a7 Glut"], collapse = ", "))
