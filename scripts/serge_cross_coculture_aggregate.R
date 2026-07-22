#!/usr/bin/env Rscript
# scripts/serge_cross_coculture_aggregate.R
#
# Aggregate the five lean Serge co-culture multiomics runs into cross-co-culture
# comparison tables + figures for the summary deck. Every current pipeline figure
# is single-co-culture; this is the only place the five are compared side by side.
#
# Reuses:
#   - the config-driven run-dir discovery from scripts/build_multiomics_index.R
#   - the jointly-significant rule from R/domain/multiomics/report_template_multiomics.Rmd
#     (alpha_joint = 0.05: nominally significant in every omics layer AND surviving
#      FDR on the Fisher-combined p; direction from per-layer NES signs).
#
# Nothing here is wired into _targets.R. Run AFTER the lean per-contrast pipelines:
#   Rscript scripts/serge_cross_coculture_aggregate.R
#
# Writes to outputs/Serge_multiomics/cross_coculture_summary/{tables,figures}/.

suppressWarnings(suppressMessages({
  library(yaml)
  library(ggplot2)
}))

`%||%` <- function(a, b) if (is.null(a)) b else a

ALPHA_JOINT <- 0.05                 # matches report_template_multiomics.Rmd
OUT_DIR     <- "outputs/Serge_multiomics/cross_coculture_summary"
TAB_DIR     <- file.path(OUT_DIR, "tables")
FIG_DIR     <- file.path(OUT_DIR, "figures")
# Display order for the co-cultures across the deck (modest -> larger response).
COCULTURE_ORDER <- c("E.coli", "L.rhamnosus", "MixSpp", "E.faecalis", "B.subtilis")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Discover the five lean run dirs (reuse build_multiomics_index.R logic) ----
discover_runs <- function() {
  cfg_files <- list.files("config", pattern = "multiomics_serge_.*\\.yaml$", full.names = TRUE)
  cfg_files <- cfg_files[!grepl("_complete\\.yaml$|_lean\\.yaml$", cfg_files)]
  stopifnot(length(cfg_files) > 0)
  runs <- lapply(cfg_files, function(f) {
    cfg     <- yaml::read_yaml(f)
    token   <- sub("\\.yaml$", "", sub("^multiomics_serge_", "", basename(f)))
    out_abs <- file.path(cfg$project$dir %||% ".", cfg$paths$out)
    hits <- if (dir.exists(out_abs)) {
      list.files(out_abs, pattern = "cross_omics_pathways_meta_analysis\\.csv$",
                 recursive = TRUE, full.names = TRUE)
    } else character(0)
    # A contrast dir may hold superseded runs; keep the newest cross_enrichment dir.
    hits <- hits[!grepl("/per_contrast/", hits)]
    if (!length(hits)) return(NULL)
    hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
    list(token = token, enrichment_dir = dirname(hits[[1]]))
  })
  runs <- Filter(Negate(is.null), runs)
  stopifnot(length(runs) > 0)
  runs
}

# ---- 2. Jointly-significant derivation (verbatim rule from the report template) ----
#' Derive jointly-significant pathways for one run.
#'
#' @param enrichment_dir Path to a run's multiomics/cross_enrichment directory.
#' @return data.frame with pathway_name, per-layer p, combined p/padj, per-layer
#'   NES and a direction label; zero rows if none qualify.
joint_sig_for_run <- function(enrichment_dir) {
  meta_file <- file.path(enrichment_dir, "cross_omics_pathways_meta_analysis.csv")
  if (!file.exists(meta_file)) return(NULL)
  meta_tbl <- read.csv(meta_file, stringsAsFactors = FALSE)
  pcols <- grep("^pval_", names(meta_tbl), value = TRUE)
  layer_names <- sub("^pval_", "", pcols)
  if (length(pcols) < 2 || !("combined_padj" %in% names(meta_tbl))) {
    return(meta_tbl[0, , drop = FALSE])
  }

  keep <- rep(TRUE, nrow(meta_tbl))
  for (pc in pcols) keep <- keep & !is.na(meta_tbl[[pc]]) & meta_tbl[[pc]] < ALPHA_JOINT
  keep <- keep & !is.na(meta_tbl$combined_padj) & meta_tbl$combined_padj < ALPHA_JOINT
  joint <- meta_tbl[keep, , drop = FALSE]
  joint <- joint[order(joint$combined_padj), , drop = FALSE]

  nes_cols <- paste0("NES_", layer_names)
  if (nrow(joint) > 0) {
    for (ln in layer_names) {
      nes_file <- file.path(enrichment_dir, sprintf("%s_enriched_pathways.csv", ln))
      joint[[paste0("NES_", ln)]] <- NA_real_
      if (file.exists(nes_file)) {
        nes_tbl <- tryCatch(read.csv(nes_file, stringsAsFactors = FALSE),
                            error = function(e) NULL)
        if (!is.null(nes_tbl) && all(c("pathway_name", "NES") %in% names(nes_tbl))) {
          nes_tbl <- nes_tbl[!duplicated(nes_tbl$pathway_name), , drop = FALSE]
          idx <- match(joint$pathway_name, nes_tbl$pathway_name)
          joint[[paste0("NES_", ln)]] <- suppressWarnings(as.numeric(nes_tbl$NES[idx]))
        }
      }
    }
    nes_mat <- as.matrix(joint[, nes_cols, drop = FALSE])
    joint$direction <- apply(nes_mat, 1, function(v) {
      if (any(is.na(v))) return("undetermined")
      if (all(v > 0)) "same (up in all layers)"
      else if (all(v < 0)) "same (down in all layers)"
      else "opposite"
    })
  } else {
    joint$direction <- character(0)
  }
  show_cols <- intersect(
    c("pathway_name", pcols, "combined_pval", "combined_padj", nes_cols, "direction"),
    names(joint))
  joint[, show_cols, drop = FALSE]
}

# ---- 3. Run, collect, verify against the 3 existing TSVs ----
runs <- discover_runs()
message("Discovered ", length(runs), " lean runs: ",
        paste(vapply(runs, `[[`, character(1), "token"), collapse = ", "))

per_run <- lapply(runs, function(r) {
  js <- joint_sig_for_run(r$enrichment_dir)
  if (is.null(js) || !nrow(js)) return(NULL)
  js$coculture <- r$token
  js
})
names(per_run) <- vapply(runs, `[[`, character(1), "token")

# Verification: where the pipeline already wrote the TSV, our derived pathway set
# must match it exactly. This gives confidence in the E.coli/MixSpp re-derivation
# (those two runs have no TSV of their own).
verify_against_tsv <- function(r, derived) {
  tsv <- file.path(r$enrichment_dir, "cross_omics_jointly_significant_pathways.tsv")
  if (!file.exists(tsv)) return(sprintf("%-12s no TSV (re-derived only)", r$token))
  ref <- read.delim(tsv, stringsAsFactors = FALSE)
  a <- sort(unique(derived$pathway_name %||% character(0)))
  b <- sort(unique(ref$pathway_name))
  ok <- identical(a, b)
  sprintf("%-12s TSV match: %s (derived %d / tsv %d)", r$token,
          if (ok) "YES" else "NO", length(a), length(b))
}
for (r in runs) {
  cat(verify_against_tsv(r, per_run[[r$token]] %||% data.frame(pathway_name = character(0))), "\n")
}

combined <- do.call(rbind, per_run)
if (is.null(combined) || !nrow(combined)) {
  stop("No jointly-significant pathways in any run; nothing to aggregate.")
}
# One row per pathway x co-culture; keep the most significant if a name repeats.
combined <- combined[!duplicated(combined[, c("pathway_name", "coculture")]), , drop = FALSE]
present <- intersect(COCULTURE_ORDER, unique(combined$coculture))
combined$coculture <- factor(combined$coculture, levels = present)

# ---- 4. Comparison tables ----
# Long table (one row per pathway x co-culture that reached joint significance).
write.csv(combined, file.path(TAB_DIR, "jointly_significant_long.csv"), row.names = FALSE)

# Pathway x co-culture direction matrix.
dir_short <- c("same (up in all layers)" = "up",
               "same (down in all layers)" = "down",
               "opposite" = "opposite", "undetermined" = "undetermined")
combined$dir_short <- dir_short[combined$direction]
mat <- reshape(
  combined[, c("pathway_name", "coculture", "dir_short")],
  idvar = "pathway_name", timevar = "coculture", direction = "wide")
names(mat) <- sub("^dir_short\\.", "", names(mat))
for (cc in present) if (!cc %in% names(mat)) mat[[cc]] <- NA_character_
# Order rows by how many co-cultures share the pathway; reorder the count vector
# together with the rows so the two never drift apart.
n_by_path <- rowSums(!is.na(mat[, present, drop = FALSE]))
ord <- order(-n_by_path, mat$pathway_name)
mat <- mat[ord, , drop = FALSE]
n_by_path <- n_by_path[ord]
mat_out <- cbind(pathway_name = mat$pathway_name, n_cocultures = n_by_path,
                 mat[, present, drop = FALSE])
write.csv(mat_out, file.path(TAB_DIR, "pathway_by_coculture_matrix.csv"), row.names = FALSE)

# Per-co-culture summary counts.
summ <- do.call(rbind, lapply(present, function(cc) {
  d <- combined[combined$coculture == cc, , drop = FALSE]
  data.frame(
    coculture   = cc,
    n_joint     = nrow(d),
    n_same_up   = sum(d$direction == "same (up in all layers)"),
    n_same_down = sum(d$direction == "same (down in all layers)"),
    n_opposite  = sum(d$direction == "opposite"),
    pct_agree   = if (nrow(d)) round(100 * mean(grepl("^same", d$direction)), 1) else NA_real_)
}))
write.csv(summ, file.path(TAB_DIR, "per_coculture_summary.csv"), row.names = FALSE)

# ---- 5. Figures (ggplot2, matches the repo plotting family) ----
theme_deck <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom")
dir_fill <- c(up = "#B2564B", down = "#3E6D9C", opposite = "#8C8C8C",
              undetermined = "#D9D9D9")

# 5a. Pathway x co-culture direction tile matrix (shared themes read as full rows).
long_tile <- combined[, c("pathway_name", "coculture", "dir_short")]
long_tile$pathway_name <- factor(long_tile$pathway_name, levels = rev(mat$pathway_name))
p_mat <- ggplot(long_tile, aes(coculture, pathway_name, fill = dir_short)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = dir_fill, name = "RNA + protein direction",
                    na.value = "grey96") +
  labs(x = NULL, y = NULL,
       title = "Jointly significant pathways across co-cultures",
       subtitle = "Cells shown only where a pathway is significant in RNA and protein (fgsea meta-analysis)") +
  theme_deck
n_path <- length(unique(long_tile$pathway_name))
ggsave(file.path(FIG_DIR, "pathway_by_coculture_matrix.png"), p_mat,
       width = 9, height = max(4, 0.28 * n_path + 2), dpi = 300, limitsize = FALSE)

# 5b. Shared vs private: how many co-cultures each pathway is jointly significant in.
freq_df <- data.frame(pathway_name = mat$pathway_name, n = n_by_path)
shared_only <- freq_df[freq_df$n >= 2, , drop = FALSE]
if (nrow(shared_only)) {
  shared_only$pathway_name <- factor(shared_only$pathway_name,
                                     levels = rev(shared_only$pathway_name))
  p_shared <- ggplot(shared_only, aes(n, pathway_name)) +
    geom_col(fill = "#4F6D7A", width = 0.7) +
    scale_x_continuous(breaks = seq_len(length(present)), expand = expansion(c(0, 0.05))) +
    labs(x = "Number of co-cultures", y = NULL,
         title = "Pathways shared across co-cultures",
         subtitle = "Jointly significant (RNA + protein) in two or more co-cultures") +
    theme_deck + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, "shared_pathways_bar.png"), p_shared,
         width = 9, height = max(3, 0.3 * nrow(shared_only) + 1.5), dpi = 300, limitsize = FALSE)
}

# 5c. Per-co-culture jointly-significant counts, split by direction.
count_long <- reshape(
  summ[, c("coculture", "n_same_up", "n_same_down", "n_opposite")],
  varying = c("n_same_up", "n_same_down", "n_opposite"),
  v.names = "n", timevar = "direction",
  times = c("up", "down", "opposite"), direction = "long")
count_long$coculture <- factor(count_long$coculture, levels = rev(present))
count_long$direction <- factor(count_long$direction, levels = c("up", "down", "opposite"))
p_counts <- ggplot(count_long, aes(n, coculture, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = dir_fill, name = "RNA + protein direction") +
  labs(x = "Jointly significant pathways", y = NULL,
       title = "Cross-omics response per co-culture",
       subtitle = "Pathways significant in RNA and protein, by agreement direction") +
  theme_deck
ggsave(file.path(FIG_DIR, "per_coculture_counts.png"), p_counts,
       width = 9, height = 3.5, dpi = 300)

message("Wrote tables to ", TAB_DIR)
message("Wrote figures to ", FIG_DIR)
cat("\nPer-co-culture jointly-significant counts:\n")
print(summ, row.names = FALSE)
