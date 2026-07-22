#!/usr/bin/env Rscript
# scripts/serge_cross_coculture_markers.R
#
# Cross-co-culture "top markers" = the top robust features from each lean run
# (consensus/consensus/robust_features.csv: features ranked robust across MOFA +
# DIABLO). Ranks the top N per co-culture, annotates the EHI_ locus tags with a
# readable description via the ehi_xp_map -> pg_matrix bridge, and marks which
# features recur across co-cultures (shared) vs appear in only one (unique).
#
#   Rscript scripts/serge_cross_coculture_markers.R
#
# Writes to outputs/Serge_multiomics/cross_coculture_summary/{tables,figures}/.

suppressWarnings(suppressMessages(library(ggplot2)))

TOPN    <- 12                    # top robust features shown per co-culture
OUT_DIR <- "outputs/Serge_multiomics/cross_coculture_summary"
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
MULTI   <- "outputs/Serge_multiomics"
COCULTURE_ORDER <- c("E.coli", "L.rhamnosus", "MixSpp", "E.faecalis", "B.subtilis")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Annotation bridge: EHI_ gene -> protein description ----
data_dir <- "data/Serge_Ankri_June2026"
xp_map <- read.csv(file.path(data_dir, "RNAseq/ehi_xp_map.csv"), stringsAsFactors = FALSE)
pg <- read.delim(file.path(data_dir, "proteomics/98858-81_hystolytica_MBR.pg_matrix.v2.tsv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
# Explode semicolon-joined protein groups so a single XP still resolves.
pg_expl <- do.call(rbind, lapply(seq_len(nrow(pg)), function(i) {
  ids <- trimws(strsplit(pg$Protein.Group[i], ";", fixed = TRUE)[[1]])
  data.frame(id = ids, desc = pg$First.Protein.Description[i], stringsAsFactors = FALSE)
}))
pg_expl <- pg_expl[!duplicated(pg_expl$id), ]
# Return the XP accession (reliably mapped) and a best-effort description.
annotate <- function(ehi) {
  xp <- xp_map$protein_id[match(ehi, xp_map$gene_id)]
  desc <- pg_expl$desc[match(xp, pg_expl$id)]
  desc <- ifelse(is.na(desc) | desc == "", "uncharacterised / not in protein set", substr(desc, 1, 60))
  data.frame(xp = ifelse(is.na(xp), "", xp), description = desc, stringsAsFactors = FALSE)
}

# ---- Collect top robust features per co-culture ----
find_robust <- function(token) {
  f <- list.files(file.path(MULTI, token), pattern = "^robust_features\\.csv$",
                  recursive = TRUE, full.names = TRUE)
  f <- f[grepl("/consensus/", f)]
  if (!length(f)) return(NULL)
  f[order(file.info(f)$mtime, decreasing = TRUE)][1]
}
per_run <- list()
for (cc in COCULTURE_ORDER) {
  f <- find_robust(cc)
  if (is.null(f)) { message("no robust_features for ", cc); next }
  x <- read.csv(f, stringsAsFactors = FALSE)
  x <- x[order(x$combined_rank), , drop = FALSE]
  top <- head(x, TOPN)
  top$coculture <- cc
  top$rank <- seq_len(nrow(top))
  a <- annotate(top$feature)
  top$xp <- a$xp
  top$description <- a$description
  per_run[[cc]] <- top[, c("coculture", "rank", "feature", "xp", "description",
                           "mofa", "diablo", "combined_rank")]
}
present <- intersect(COCULTURE_ORDER, names(per_run))
long <- do.call(rbind, per_run)
write.csv(long, file.path(TAB_DIR, "robust_top_by_coculture.csv"), row.names = FALSE)

# ---- Shared vs unique across the top-N sets ----
feat_cc <- lapply(per_run, function(d) d$feature)
all_feats <- unique(unlist(feat_cc))
presence <- sapply(present, function(cc) all_feats %in% feat_cc[[cc]])
rownames(presence) <- all_feats
n_cc <- rowSums(presence)
overlap <- data.frame(feature = all_feats, description = annotate(all_feats)$description,
                      n_cocultures = as.integer(n_cc),
                      cocultures = apply(presence, 1, function(v) paste(present[v], collapse = ";")),
                      stringsAsFactors = FALSE)
overlap <- overlap[order(-overlap$n_cocultures, overlap$feature), , drop = FALSE]
write.csv(overlap, file.path(TAB_DIR, "robust_feature_overlap.csv"), row.names = FALSE)

n_shared <- sum(n_cc >= 2)
summ <- data.frame(
  coculture = present,
  n_robust_total = sapply(present, function(cc) {
    f <- find_robust(cc); if (is.null(f)) NA_integer_ else nrow(read.csv(f))
  }),
  n_top_shown = sapply(present, function(cc) nrow(per_run[[cc]])),
  n_unique_in_top = sapply(present, function(cc)
    sum(overlap$n_cocultures[match(feat_cc[[cc]], overlap$feature)] == 1)))
write.csv(summ, file.path(TAB_DIR, "robust_feature_summary.csv"), row.names = FALSE)

message(sprintf("Top-%d robust features per co-culture; %d of %d distinct features are shared by >=2 co-cultures.",
                TOPN, n_shared, length(all_feats)))
print(summ, row.names = FALSE)

# ---- Figure: shared robust features (features in >=2 co-cultures' top-N) ----
theme_deck <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text.y = element_text(size = 9))
shared_feats <- overlap$feature[overlap$n_cocultures >= 2]
if (length(shared_feats)) {
  tile <- long[long$feature %in% shared_feats, c("feature", "coculture", "rank")]
  # Label rows with feature + short description for readability.
  lab <- setNames(paste0(overlap$feature, " — ", substr(overlap$description, 1, 34)), overlap$feature)
  tile$row <- lab[tile$feature]
  ord <- overlap$feature[overlap$feature %in% shared_feats]
  tile$row <- factor(tile$row, levels = rev(lab[ord]))
  tile$coculture <- factor(tile$coculture, levels = present)
  p <- ggplot(tile, aes(coculture, row, fill = rank)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = rank), size = 3, colour = "white") +
    scale_fill_gradient(low = "#3E6D9C", high = "#B9C7D6", name = "rank in co-culture", trans = "reverse") +
    labs(x = NULL, y = NULL,
         title = "Robust integration features shared across co-cultures",
         subtitle = "Features in the top-12 robust set of two or more co-cultures (rank shown in cell)") +
    theme_deck
  ggsave(file.path(FIG_DIR, "robust_features_shared.png"), p,
         width = 9, height = max(2.5, 0.4 * length(shared_feats) + 1.6), dpi = 300, limitsize = FALSE)
} else {
  message("No robust feature is shared across >=2 co-cultures' top-N sets.")
}

# ---- Figure: how co-culture-specific the top robust features are ----
spec <- data.frame(coculture = factor(present, levels = rev(present)),
                   unique = summ$n_unique_in_top,
                   shared = summ$n_top_shown - summ$n_unique_in_top)
spec_long <- reshape(spec, varying = c("unique", "shared"), v.names = "n",
                     timevar = "kind", times = c("co-culture-specific", "shared with another"),
                     direction = "long")
p2 <- ggplot(spec_long, aes(n, coculture, fill = kind)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("co-culture-specific" = "#4F6D7A", "shared with another" = "#C6A15B"), name = NULL) +
  labs(x = "Top robust features (of 12 shown)", y = NULL,
       title = "Top robust features are mostly co-culture-specific",
       subtitle = "Split of each co-culture's top-12 robust integration features") +
  theme_deck
ggsave(file.path(FIG_DIR, "robust_features_specificity.png"), p2, width = 9, height = 3.2, dpi = 300)

message("Wrote tables + figures to ", OUT_DIR)
