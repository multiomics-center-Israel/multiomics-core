library(readxl)

# --- Load both datasets ---
elah_file <- "/home/multiomics-core/data/elah_pick_rna/Proteomics/DEP_results_with_annot.Aug2025.from_Elah.v2EP - Copy.xlsx"
elah <- read_excel(elah_file, sheet="Prot_results")

our_file <- "/home/multiomics-core/outputs/elah_pick_protein/Results_Elah_Pick_Protein_1/proteomics/Datasets/limma_multimp_summary_p0.05.tsv"
ours <- read.delim(our_file, stringsAsFactors=FALSE)

cat("=== ID FORMAT ANALYSIS ===\n")
cat("Elah IDs (first 10):\n")
print(head(elah$name, 10))
cat("\nOur IDs (first 10):\n")
print(head(ours$FeatureID, 10))

# Normalize our IDs: strip version numbers (.1, .2) and take first ID before semicolons
our_base_ids <- gsub("\\.[0-9]+$", "", sapply(strsplit(ours$FeatureID, ";"), `[`, 1))
cat("\nOur normalized IDs (first 10):\n")
print(head(our_base_ids, 10))

# Check for multi-protein groups in our data
multi <- grepl(";", ours$FeatureID)
cat("\nOur proteins with multiple IDs (protein groups):", sum(multi), "out of", nrow(ours), "\n")

# Now find overlap
common <- intersect(elah$name, our_base_ids)
cat("\nCommon proteins after normalization:", length(common), "\n")
cat("Only in Elah:", sum(!elah$name %in% our_base_ids), "\n")
cat("Only in ours:", sum(!our_base_ids %in% elah$name), "\n")

# --- Build matched dataset ---
ours$base_id <- our_base_ids

# Handle duplicates in our data (from protein groups)
ours_dedup <- ours[!duplicated(ours$base_id), ]
cat("\nOur proteins after dedup:", nrow(ours_dedup), "\n")

# Merge
merged <- merge(
  data.frame(id = elah$name,
             elah_ratio = as.numeric(elah$OXY_vs_ANA_ratio),
             elah_pval = as.numeric(elah$OXY_vs_ANA_p.val),
             elah_padj = as.numeric(elah$OXY_vs_ANA_p.adj),
             elah_sig = as.logical(elah$significant),
             stringsAsFactors = FALSE),
  data.frame(id = ours_dedup$base_id,
             our_linearFC = ours_dedup$linearFC.imputs.Oxidativestress_vs_Anaerobic,
             our_pval = ours_dedup$pvalue.imputs.Oxidativestress_vs_Anaerobic,
             our_padj = ours_dedup$padj.imputs.Oxidativestress_vs_Anaerobic,
             stringsAsFactors = FALSE),
  by = "id"
)
cat("\nMerged dataset:", nrow(merged), "proteins\n")

# === FOLD CHANGE ANALYSIS ===
cat("\n=== FOLD CHANGE ANALYSIS ===\n")
cat("Our linearFC range:", range(merged$our_linearFC, na.rm=TRUE), "\n")
cat("Our linearFC quantiles:\n")
print(quantile(merged$our_linearFC, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE))

cat("\nElah ratio (log2FC) range:", range(merged$elah_ratio, na.rm=TRUE), "\n")
cat("Elah ratio quantiles:\n")
print(quantile(merged$elah_ratio, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE))

# Convert our linearFC to log2FC
# linearFC is a signed linear fold change: negative = downregulated
# For negative: actual ratio = 1/|linearFC|, log2FC = -log2(|linearFC|)
# For positive: actual ratio = linearFC, log2FC = log2(linearFC)
merged$our_log2fc <- ifelse(merged$our_linearFC < 0,
                             -log2(abs(merged$our_linearFC)),
                             log2(abs(merged$our_linearFC)))

# Spot check
cat("\n=== SPOT-CHECK FC CONVERSION ===\n")
test_ids <- c("KAE8301209", "KAE8301222", "KAE8301224")
for (tid in test_ids) {
  idx <- which(merged$id == tid)
  if (length(idx) > 0) {
    cat("Protein", tid, ":\n")
    cat("  Elah ratio (log2FC):", merged$elah_ratio[idx], "\n")
    cat("  Our linearFC:", merged$our_linearFC[idx], "\n")
    cat("  Our converted log2FC:", merged$our_log2fc[idx], "\n\n")
  }
}

cat("Our log2FC range:", range(merged$our_log2fc, na.rm=TRUE), "\n")

# === EFFECT SIZE CORRELATION ===
cat("\n=== EFFECT SIZE CORRELATION ===\n")
valid <- complete.cases(merged$elah_ratio, merged$our_log2fc)
cat("Valid pairs:", sum(valid), "\n")

cor_p <- cor(merged$elah_ratio[valid], merged$our_log2fc[valid], method="pearson")
cor_s <- cor(merged$elah_ratio[valid], merged$our_log2fc[valid], method="spearman")
cat("Pearson correlation (log2FC):", round(cor_p, 4), "\n")
cat("Spearman correlation (log2FC):", round(cor_s, 4), "\n")

mad_val <- mean(abs(merged$elah_ratio[valid] - merged$our_log2fc[valid]))
cat("Mean absolute difference in log2FC:", round(mad_val, 4), "\n")

# === SIGNIFICANCE COMPARISON ===
cat("\n=== SIGNIFICANCE COMPARISON ===\n")

elah_sig_005 <- merged$elah_padj < 0.05
our_sig_005 <- merged$our_padj < 0.05

cat("\nAt padj < 0.05:\n")
cat("  Both significant:    ", sum(elah_sig_005 & our_sig_005, na.rm=TRUE), "\n")
cat("  Elah only:           ", sum(elah_sig_005 & !our_sig_005, na.rm=TRUE), "\n")
cat("  Ours only:           ", sum(!elah_sig_005 & our_sig_005, na.rm=TRUE), "\n")
cat("  Neither:             ", sum(!elah_sig_005 & !our_sig_005, na.rm=TRUE), "\n")
cat("  Total Elah sig:      ", sum(elah_sig_005, na.rm=TRUE), "\n")
cat("  Total our sig:       ", sum(our_sig_005, na.rm=TRUE), "\n")

elah_sig_01 <- merged$elah_padj < 0.1
our_sig_01 <- merged$our_padj < 0.1

cat("\nAt padj < 0.1:\n")
cat("  Both significant:    ", sum(elah_sig_01 & our_sig_01, na.rm=TRUE), "\n")
cat("  Elah only:           ", sum(elah_sig_01 & !our_sig_01, na.rm=TRUE), "\n")
cat("  Ours only:           ", sum(!elah_sig_01 & our_sig_01, na.rm=TRUE), "\n")
cat("  Neither:             ", sum(!elah_sig_01 & !our_sig_01, na.rm=TRUE), "\n")

# Direction agreement
both_sig <- which(elah_sig_005 & our_sig_005)
if (length(both_sig) > 0) {
  same_dir <- sum(sign(merged$elah_ratio[both_sig]) == sign(merged$our_log2fc[both_sig]), na.rm=TRUE)
  cat("\nAmong", length(both_sig), "jointly significant (padj<0.05):\n")
  cat("  Same direction:", same_dir, "(", round(100*same_dir/length(both_sig), 1), "%)\n")
  cat("  Opposite direction:", length(both_sig) - same_dir, "\n")
}

# === P-VALUE CORRELATION ===
cat("\n=== P-VALUE CORRELATION ===\n")
valid_p <- complete.cases(merged$elah_padj, merged$our_padj) &
           merged$elah_padj > 0 & merged$our_padj > 0
cat("Valid p-value pairs:", sum(valid_p), "\n")
cor_pval <- cor(-log10(merged$elah_padj[valid_p]), -log10(merged$our_padj[valid_p]), method="spearman")
cat("Spearman correlation of -log10(padj):", round(cor_pval, 4), "\n")
cor_pval_raw <- cor(-log10(merged$elah_pval[valid_p]), -log10(merged$our_pval[valid_p]), method="spearman")
cat("Spearman correlation of -log10(raw pval):", round(cor_pval_raw, 4), "\n")

# === TOP DISCORDANT PROTEINS ===
cat("\n=== TOP DISCORDANT PROTEINS ===\n")
cat("(Significant in one analysis but not the other, at padj < 0.05)\n\n")

elah_only <- merged[elah_sig_005 & !our_sig_005 & !is.na(elah_sig_005) & !is.na(our_sig_005), ]
elah_only <- elah_only[order(elah_only$elah_padj), ]
cat("--- Significant in Elah only (top 15 by Elah padj) ---\n")
print(head(elah_only[, c("id", "elah_ratio", "elah_padj", "our_log2fc", "our_padj")], 15))

ours_only <- merged[!elah_sig_005 & our_sig_005 & !is.na(elah_sig_005) & !is.na(our_sig_005), ]
ours_only <- ours_only[order(ours_only$our_padj), ]
cat("\n--- Significant in ours only (top 15 by our padj) ---\n")
print(head(ours_only[, c("id", "elah_ratio", "elah_padj", "our_log2fc", "our_padj")], 15))

# === OVERALL SUMMARY ===
cat("\n\n========================================\n")
cat("=== OVERALL SUMMARY ===\n")
cat("========================================\n")
cat("Total proteins: Elah =", nrow(elah), ", Ours =", nrow(ours), "\n")
cat("Matched proteins:", nrow(merged), "\n")
cat("Effect size (log2FC) Pearson r:", round(cor_p, 4), "\n")
cat("Effect size (log2FC) Spearman r:", round(cor_s, 4), "\n")
cat("-log10(padj) Spearman r:", round(cor_pval, 4), "\n")
cat("-log10(raw pval) Spearman r:", round(cor_pval_raw, 4), "\n")

cat("\n--- Method Summary ---\n")
cat("Elah: DEP R package, t-test, BH correction\n")
cat("Ours: limma moderated t-test, Perseus-style multi-imputation (10 rounds), consensus filtering\n")
cat("Key difference: Our pipeline uses moderated t-statistics (empirical Bayes) and\n")
cat("  multiple imputation with consensus voting, which tends to be more conservative\n")
cat("  for proteins with high missing-value rates\n")

dep_cols <- grep("DEP|dep|centered|imputed", colnames(elah), value=TRUE, ignore.case=TRUE)
if (length(dep_cols) > 0) cat("\nPossible DEP-related columns in Elah:", paste(dep_cols, collapse=", "), "\n")
