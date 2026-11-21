# harmonize_sn_results.R
# Harmonize the four primary single-nucleus (snRNA-seq) WildR vs SPF analyses
# (MAST, popAlign, D-SPIN, pySCENIC) into one long feature-level table plus
# cell-type and class-rollup summaries.
#
# Display notebook: notebooks/0X_Integrated-Results_SN.qmd
# Run directly:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/harmonize_sn_results.R
#
# Inputs are the *saved* outputs of each method (nothing upstream is recomputed).
# The only new statistic computed here is an animal-level Welch t-test on D-SPIN
# relative responses (6 WildR vs 6 SPF animals per program), because D-SPIN
# itself does not emit p-values.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
run_date <- as.character(Sys.Date())
sig_q <- 0.1  # harmonized significance threshold applied to every method's adjusted p

log <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), glue(..., .envir = parent.frame()))

# ---- Source registry: the exact artifact each method is loaded from ----------
sn_sources <- list(
  mast = list(
    hits     = file.path(wkdir, "data/interim/mast/snRNASeq/mast_signif_hits_2026-08-17_p_adjust_0.1.csv"),
    per_unit = file.path(wkdir, "data/interim/mast/snRNASeq/{class,subclass,supertype}_label/")
  ),
  popalign = list(
    catalog = file.path(wkdir, "tables/popAlign/snRNASeq/run_2026-10-04_synth-2026-10-05/popAlign_subpopulation_catalog.csv"),
    tbl_dir = file.path(wkdir, "tables/popAlign/snRNASeq/run_2026-10-04_synth-2026-10-05")
  ),
  dspin = list(
    fit_summary = file.path(wkdir, "data/interim/dspin/fit_summary.csv"),
    unit_dir    = file.path(wkdir, "data/interim/dspin"),
    stars       = file.path(wkdir, "tables/DSPIN/star_programs_david.csv")
  ),
  scenic = list(
    all_splits = file.path(wkdir, "data/processed/scenic_by_class/diff_regulons_all_splits.csv"),
    run_status = file.path(wkdir, "data/processed/scenic_by_class/run_status.csv")
  ),
  hierarchy_src = file.path(wkdir, "data/interim/dspin/cell_meta_keep_cell.csv")
)

out_tbl_dir <- file.path(wkdir, "tables/SN_Integration")
out_int_dir <- file.path(wkdir, "data/interim/sn_integration")

# Common column set every loader returns
harmonized_cols <- c(
  "method", "tissue", "level", "celltype_label", "celltype_name",
  "feature_type", "feature_id", "genes", "test", "effect", "effect_metric",
  "direction", "p_value", "p_adj", "padj_scope", "native_hit", "native_rule",
  "inference_unit", "n_SPF", "n_WildR", "source_file"
)

dir_from_sign <- function(x) case_when(is.na(x) ~ NA_character_, x > 0 ~ "WildR_up", x < 0 ~ "SPF_up", TRUE ~ "none")

# ---- Allen hierarchy lookup (supertype -> subclass -> class) -----------------
build_hierarchy <- function(src = sn_sources$hierarchy_src) {
  log("Building Allen hierarchy lookup from {basename(src)}")
  fread(src, select = c("class_label", "class_name", "subclass_label", "subclass_name",
                        "supertype_label", "supertype_name")) %>%
    as_tibble() %>% distinct()
}

# Long form: one row per (level, label) with its parent class/subclass
hierarchy_long <- function(hier) {
  bind_rows(
    hier %>% distinct(class_label, class_name) %>%
      transmute(level = "class", celltype_label = class_label, celltype_name_ref = class_name,
                class_label, class_name, subclass_label = NA_character_, subclass_name = NA_character_),
    hier %>% distinct(class_label, class_name, subclass_label, subclass_name) %>%
      transmute(level = "subclass", celltype_label = subclass_label, celltype_name_ref = subclass_name,
                class_label, class_name, subclass_label, subclass_name),
    hier %>% distinct(class_label, class_name, subclass_label, subclass_name, supertype_label, supertype_name) %>%
      transmute(level = "supertype", celltype_label = supertype_label, celltype_name_ref = supertype_name,
                class_label, class_name, subclass_label, subclass_name)
  ) %>% distinct(level, celltype_label, .keep_all = TRUE)
}

# ---- MAST --------------------------------------------------------------------
# Hurdle LRT on microbiomeWildR, ~ microbiome + log_std_nFeature_RNA, per
# tissue x Allen unit (R_scripts/mast_functions.R). Saved table holds hits with
# BH p_adjust < 0.1 only (full per-unit results live in *_label/ RDS files).
load_mast <- function(f = sn_sources$mast$hits) {
  log("Loading MAST: {basename(f)}")
  read_csv(f, show_col_types = FALSE) %>%
    transmute(
      method = "MAST", tissue, level = cell_type_level,
      celltype_label = cell_type, celltype_name = cell_type_name,
      feature_type = "gene", feature_id = str_remove(primerid, "-GRCm39$"),
      genes = feature_id, test = "hurdle_LRT",
      effect = coef, effect_metric = "MAST logFC coef (log2 CPM)",
      direction = dir_from_sign(coef), p_value = pr_chisq, p_adj = p_adjust,
      padj_scope = "BH within unit (all tested genes)",
      native_hit = p_adjust < 0.1, native_rule = "p_adjust < 0.1",
      inference_unit = "nuclei (<=2000/group; no animal term)",
      n_SPF = NA_real_, n_WildR = NA_real_, source_file = f
    )
}

# ---- popAlign ----------------------------------------------------------------
# One Gaussian-mixture component = one subpopulation. Two animal-level tests per
# component: abundance (permutation over animals) and state shift (permutation).
load_popalign <- function(f = sn_sources$popalign$catalog) {
  log("Loading popAlign catalog: {basename(f)}")
  cat <- read_csv(f, show_col_types = FALSE)
  base <- cat %>% transmute(
    method = "popAlign", tissue, level, celltype_label = label, celltype_name = name,
    feature_type = "subpopulation", feature_id = paste0("comp", comp),
    genes = top_markers, source_file = f,
    inference_unit = "animals (6 vs 6, permutation)",
    n_SPF = NA_real_, n_WildR = NA_real_,
    log2fc_frac, p_abund, q_abund, abundance_hit, shift, p_state, q_state, state_hit
  )
  bind_rows(
    base %>% transmute(
      method, tissue, level, celltype_label, celltype_name, feature_type, feature_id, genes,
      test = "abundance", effect = log2fc_frac,
      effect_metric = "log2 FC of mean component fraction (WildR/SPF)",
      direction = dir_from_sign(log2fc_frac), p_value = p_abund, p_adj = q_abund,
      padj_scope = "BH across all popAlign components",
      native_hit = as.logical(abundance_hit), native_rule = "synthesis abundance_hit",
      inference_unit, n_SPF, n_WildR, source_file),
    base %>% transmute(
      method, tissue, level, celltype_label, celltype_name, feature_type, feature_id, genes,
      test = "state_shift", effect = shift,
      effect_metric = "state-shift statistic (unsigned)",
      direction = NA_character_, p_value = p_state, p_adj = q_state,
      padj_scope = "BH across all popAlign components",
      native_hit = as.logical(state_hit), native_rule = "synthesis state_hit",
      inference_unit, n_SPF, n_WildR, source_file)
  )
}

# ---- D-SPIN ------------------------------------------------------------------
# Per supertype x tissue unit: program responses h per animal (sample_id),
# relative to SPF. Effect = mean WildR relative_h. Test added here: Welch t-test
# of relative_h WildR vs SPF across animals, BH across all D-SPIN programs.
load_dspin <- function(fs = sn_sources$dspin$fit_summary,
                       unit_dir = sn_sources$dspin$unit_dir,
                       stars_f = sn_sources$dspin$stars) {
  log("Loading D-SPIN units from {unit_dir}")
  fit <- read_csv(fs, show_col_types = FALSE) %>% filter(status == "ok")
  stars <- read_csv(stars_f, show_col_types = FALSE) %>%
    distinct(tissue, supertype_label, program_id) %>% mutate(star = TRUE)

  welch_p <- function(x, g) {
    a <- x[g == "WildR"]; b <- x[g == "SPF"]
    if (length(a) < 2 || length(b) < 2 || sd(c(a, b)) == 0) return(NA_real_)
    tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_)
  }

  per_unit <- map(seq_len(nrow(fit)), function(i) {
    ud <- file.path(unit_dir, fit$unit_id[i])
    rr <- file.path(ud, "relative_responses.csv"); sm <- file.path(ud, "sample_meta.csv")
    if (!file.exists(rr) || !file.exists(sm)) return(NULL)
    resp <- read_csv(rr, show_col_types = FALSE) %>%
      left_join(read_csv(sm, show_col_types = FALSE) %>% select(sample_id, microbiome), by = "sample_id")
    pg <- file.path(ud, "programs_genes.csv")
    genes <- if (file.exists(pg)) {
      read_csv(pg, show_col_types = FALSE) %>% group_by(program_id) %>%
        summarise(genes = paste(head(gene, 10), collapse = ", "), .groups = "drop")
    } else tibble(program_id = character(), genes = character())
    resp %>% group_by(program_id) %>%
      summarise(
        effect = mean(relative_h[microbiome == "WildR"], na.rm = TRUE),
        p_value = welch_p(relative_h, microbiome),
        n_SPF = sum(microbiome == "SPF"), n_WildR = sum(microbiome == "WildR"),
        .groups = "drop") %>%
      left_join(genes, by = "program_id") %>%
      # fit_summary leaves supertype_label blank for the 5 pilot units; unit_id is always "<tissue>__<label>"
      mutate(tissue = fit$tissue[i],
             celltype_label = coalesce(fit$supertype_label[i], str_remove(fit$unit_id[i], "^[A-Z]+__")),
             celltype_name = fit$supertype_name[i], source_file = rr)
  }) %>% bind_rows()

  per_unit %>%
    left_join(stars, by = c("tissue", "celltype_label" = "supertype_label", "program_id")) %>%
    mutate(p_adj = p.adjust(p_value, "BH")) %>%
    transmute(
      method = "DSPIN", tissue, level = "supertype", celltype_label, celltype_name,
      feature_type = "program", feature_id = program_id, genes,
      test = "welch_t_relative_h", effect,
      effect_metric = "mean WildR relative h (SPF-centred)",
      direction = dir_from_sign(effect), p_value, p_adj,
      padj_scope = "BH across all D-SPIN programs (computed here)",
      native_hit = coalesce(star, FALSE), native_rule = "star program (same sign in all WildR animals)",
      inference_unit = "animals (6 vs 6, Welch t)",
      n_SPF, n_WildR, source_file
    )
}

# ---- pySCENIC ----------------------------------------------------------------
# Production run: one GRN per tissue x Allen class (raw UMIs); WildR vs SPF =
# Wilcoxon on AUCell within split. Regulons are split-specific.
load_scenic <- function(f = sn_sources$scenic$all_splits) {
  log("Loading pySCENIC class splits: {basename(f)}")
  read_csv(f, show_col_types = FALSE) %>%
    transmute(
      method = "pySCENIC", tissue, level = split_level,
      celltype_label = class_label, celltype_name = class_name,
      feature_type = "regulon", feature_id = regulon, genes = tf,
      test = "wilcoxon_AUCell", effect = delta,
      effect_metric = "delta mean AUCell (WildR - SPF)",
      direction = dir_from_sign(delta), p_value = p, p_adj = padj,
      padj_scope = "BH within split (all regulons)",
      native_hit = as.logical(hit), native_rule = "padj < 0.05",
      inference_unit = "nuclei (Wilcoxon; no animal term)",
      n_SPF, n_WildR, source_file = f
    )
}

# ---- Harmonize ---------------------------------------------------------------
harmonize_sn <- function(hier) {
  hl <- hierarchy_long(hier)
  bind_rows(load_mast(), load_popalign(), load_dspin(), load_scenic()) %>%
    select(all_of(harmonized_cols)) %>%
    mutate(significant = !is.na(p_adj) & p_adj < sig_q) %>%
    left_join(hl, by = c("level", "celltype_label")) %>%
    mutate(celltype_name = coalesce(celltype_name, celltype_name_ref),
           class_label = if_else(level == "class", celltype_label, class_label),
           unit_key = glue("{level}__{tissue}__{celltype_label}")) %>%
    select(method, unit_key, tissue, level, celltype_label, celltype_name,
           class_label, class_name, subclass_label, subclass_name,
           feature_type, feature_id, genes, test, effect, effect_metric, direction,
           p_value, p_adj, padj_scope, significant, native_hit, native_rule,
           inference_unit, n_SPF, n_WildR, source_file)
}

summarise_celltypes <- function(long) {
  long %>%
    group_by(method, unit_key, tissue, level, celltype_label, celltype_name, class_label, class_name) %>%
    summarise(n_tested = n(), n_sig = sum(significant), n_native = sum(native_hit, na.rm = TRUE),
              min_padj = suppressWarnings(min(p_adj, na.rm = TRUE)), .groups = "drop") %>%
    mutate(min_padj = if_else(is.finite(min_padj), min_padj, NA_real_),
           # the saved MAST table holds hits only, so its row count is not the number of genes tested
           n_tested = if_else(method == "MAST", NA_integer_, n_tested)) %>%
    pivot_wider(names_from = method, values_from = c(n_tested, n_sig, n_native, min_padj),
                names_glue = "{method}_{.value}") %>%
    mutate(n_methods_sig = rowSums(across(ends_with("_n_sig"), ~ coalesce(.x, 0L) > 0))) %>%
    arrange(desc(n_methods_sig), tissue, level, celltype_label)
}

# Every hit rolled up to its parent Allen class, so class-only (pySCENIC) and
# supertype-only (D-SPIN) results can be compared on one axis.
rollup_to_class <- function(long) {
  long %>% filter(!is.na(class_label)) %>%
    group_by(tissue, class_label, class_name, method) %>%
    summarise(n_tested = n(), n_sig = sum(significant),
              n_units_sig = n_distinct(unit_key[significant]), .groups = "drop") %>%
    arrange(tissue, class_label, method)
}

run_harmonize_sn <- function() {
  dir.create(out_tbl_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_int_dir, recursive = TRUE, showWarnings = FALSE)
  # Never overwrite: if today's outputs already exist, tag this run <date>_r2, _r3, ...
  stems <- function(tag) c(file.path(out_tbl_dir, glue("SN_Integration_harmonized_long_{tag}.csv")),
                           file.path(out_int_dir, glue("sn_harmonized_results_{tag}.rds")))
  run_tag <- run_date; k <- 1
  while (any(file.exists(stems(run_tag)))) { k <- k + 1; run_tag <- glue("{run_date}_r{k}") }
  log("Output tag: {run_tag}")
  safe_out <- function(path) {
    if (file.exists(path)) stop(glue("Refusing to overwrite existing artifact: {path}"))
    path
  }

  hier <- build_hierarchy()
  write_csv(hier, safe_out(file.path(out_int_dir, glue("allen_hierarchy_lookup_{run_tag}.csv"))))

  long <- harmonize_sn(hier)
  summ <- summarise_celltypes(long)
  roll <- rollup_to_class(long)
  registry <- long %>% group_by(method, source_file) %>%
    summarise(n_rows = n(), n_sig = sum(significant), .groups = "drop") %>%
    mutate(file_mtime = as.character(file.mtime(source_file)))
  # D-SPIN source_file is per-unit; collapse to the directory for the registry
  registry <- registry %>%
    mutate(source_file = if_else(method == "DSPIN", sn_sources$dspin$unit_dir, source_file)) %>%
    group_by(method, source_file) %>%
    summarise(n_rows = sum(n_rows), n_sig = sum(n_sig), file_mtime = max(file_mtime), .groups = "drop")

  write_csv(long, safe_out(file.path(out_tbl_dir, glue("SN_Integration_harmonized_long_{run_tag}.csv"))))
  write_csv(summ, safe_out(file.path(out_tbl_dir, glue("SN_Integration_celltype_summary_{run_tag}.csv"))))
  write_csv(roll, safe_out(file.path(out_tbl_dir, glue("SN_Integration_class_rollup_{run_tag}.csv"))))
  write_csv(registry, safe_out(file.path(out_tbl_dir, glue("SN_Integration_source_registry_{run_tag}.csv"))))
  saveRDS(long, safe_out(file.path(out_int_dir, glue("sn_harmonized_results_{run_tag}.rds"))))

  log("Harmonized {nrow(long)} rows; {sum(long$significant)} at p_adj < {sig_q}")
  print(count(long, method, test, significant))
  invisible(list(long = long, summary = summ, rollup = roll, registry = registry))
}

if (sys.nframe() == 0) run_harmonize_sn()
