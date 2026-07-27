#!/usr/bin/env Rscript
#
# Candidate-gene group comparison: classical one-way ANOVA and the limma
# moderated F-test on the same genes, samples and BH correction, with a heatmap
# per test (PDF + SVG) and every computed value written to tab-delimited files.
#
# Usage
#   Rscript candidate_gene_anova_heatmaps.R \
#     --genes=gene_ids.txt \
#     --counts=../RNAseq/rna_counts.tsv \
#     --meta=../RNAseq/rna_metadata.csv \
#     --groups=E.histolytica,B.subtilis_E.histolytica,E.coli_E.histolytica,E.faecalis_E.histolytica,L.rhamnosus_E.histolytica \
#     --outdir=biofilm_candidate_clustering/candidate_anova_script \
#     --prefix=candidates_nomix
#
#   Run it a second time with the extra group appended to --groups (and a
#   different --prefix) to get the MixSpp-included counterpart.
#
# Outputs, all under --outdir and stamped with --prefix:
#   <prefix>_results.tsv             every gene: F, p and BH-FDR for both tests,
#                                    per-group means, peak group, annotation
#   <prefix>_anova_significant.tsv   genes with fdr_anova < --fdr
#   <prefix>_limma_significant.tsv   genes with fdr_limma < --fdr
#   <prefix>_matrix_log2cpm.tsv      the exact matrix the tests ran on
#   <prefix>_sample_groups.tsv       sample -> group map, in drawn column order
#   <prefix>_anova_heatmap.{pdf,svg} heatmap of the ANOVA-significant genes
#   <prefix>_limma_heatmap.{pdf,svg} heatmap of the limma-significant genes
#   <prefix>_*_heatmap_zvalues.tsv   the z-scores actually plotted
#
# SVG is written by the cairo device, which stores glyphs as outlines: the
# figure stays vector and scales cleanly, but the labels are not re-typeable.

suppressPackageStartupMessages({
  library(data.table); library(limma); library(ComplexHeatmap)
  library(circlize); library(grid)
})

# ---------------------------------------------------------------- arguments --

#' Parse `--key=value` command-line arguments with defaults
#'
#' @param defaults Named list of default values; names are the accepted keys.
#' @param argv Character vector of raw arguments.
#' @return The defaults list with supplied values overridden (as characters).
parse_args <- function(defaults, argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) && any(argv %in% c("-h", "--help"))) {
    cat("see the header of this file for usage\n"); quit(status = 0)
  }
  bad <- argv[!grepl("^--[^=]+=", argv)]
  if (length(bad)) stop("arguments must look like --key=value; got: ",
                        paste(bad, collapse = " "), call. = FALSE)
  keys <- sub("^--([^=]+)=.*$", "\\1", argv)
  vals <- sub("^--[^=]+=", "", argv)
  unknown <- setdiff(keys, names(defaults))
  if (length(unknown)) stop("unknown argument(s): ", paste(unknown, collapse = ", "),
                            "\nknown: ", paste(names(defaults), collapse = ", "),
                            call. = FALSE)
  defaults[keys] <- as.list(vals)
  defaults
}

opt <- parse_args(list(
  genes       = NULL,   # file of gene IDs (one per line, or a table whose first column holds them)
  counts      = NULL,   # gene-by-sample count table (TSV, first column = gene ID)
  meta        = NULL,   # sample metadata (CSV/TSV) with a Group column
  groups      = "",     # comma-separated groups to compare; "" = all groups present
  group_col   = "Group",
  sample_col  = "SampleName",
  annotation  = "",     # optional per-gene annotation table (TSV/CSV)
  ann_id_col  = "Gene ID",
  outdir      = ".",
  prefix      = "candidates",
  fdr         = "0.10",
  min_cpm     = "1",
  min_samples = "4",
  seed        = "42"))

for (req in c("genes", "counts", "meta")) {
  if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
}
fdr_cut     <- as.numeric(opt$fdr)
min_cpm     <- as.numeric(opt$min_cpm)
min_samples <- as.integer(opt$min_samples)
seed        <- as.integer(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------- inputs --

#' Read a delimited table, guessing the separator from the file extension
#'
#' @param path File path (.csv is comma-separated, anything else tab-separated).
#' @return A data.table.
read_table_any <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    pkg <- c("readxl", "openxlsx")[c(requireNamespace("readxl", quietly = TRUE),
                                     requireNamespace("openxlsx", quietly = TRUE))][1]
    if (is.na(pkg)) stop("reading ", basename(path), " needs readxl or openxlsx, ",
                         "neither is installed here — export the sheet to TSV first",
                         call. = FALSE)
    return(as.data.table(if (pkg == "readxl") readxl::read_excel(path)
                         else openxlsx::read.xlsx(path)))
  }
  fread(path, sep = if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t")
}

#' Read the gene list: a bare one-per-line file, or the first column of a table
#'
#' @param path File path.
#' @return Character vector of unique, non-empty gene IDs.
read_gene_list <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  ids <- if (grepl("[\t,]", first)) as.character(read_table_any(path)[[1]])
         else readLines(path, warn = FALSE)
  ids <- trimws(ids)
  unique(ids[nzchar(ids)])
}

gene_ids  <- read_gene_list(opt$genes)
counts_dt <- read_table_any(opt$counts)
counts    <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- as.character(counts_dt[[1]])
storage.mode(counts) <- "numeric"

meta <- read_table_any(opt$meta)
if (!opt$group_col %in% names(meta))
  stop("--meta has no column '", opt$group_col, "'; columns are: ",
       paste(names(meta), collapse = ", "), call. = FALSE)
# Sample IDs live either in a named column or in the table's first column.
samp_ids <- as.character(
  if (opt$sample_col %in% names(meta)) meta[[opt$sample_col]] else meta[[1]])
samp_grp <- as.character(meta[[opt$group_col]])

wanted <- if (nzchar(opt$groups)) {
  trimws(strsplit(opt$groups, ",")[[1]])
} else {
  unique(samp_grp[samp_ids %in% colnames(counts)])
}
missing_groups <- setdiff(wanted, unique(samp_grp))
if (length(missing_groups))
  stop("group(s) not in --meta: ", paste(missing_groups, collapse = ", "), call. = FALSE)

keep    <- samp_grp %in% wanted & samp_ids %in% colnames(counts)
ord     <- order(match(samp_grp[keep], wanted))          # control/reference group first
samples <- samp_ids[keep][ord]
grp     <- factor(samp_grp[keep][ord], levels = wanted)
if (any(table(grp) < 2))
  stop("every group needs at least 2 samples; got ",
       paste(sprintf("%s=%d", names(table(grp)), table(grp)), collapse = ", "),
       call. = FALSE)

# ------------------------------------------------------------ normalisation --

# Library sizes come from the full count table: scaling by the candidate-gene
# subtotal would not correct for sequencing depth.
lib_sizes <- colSums(counts, na.rm = TRUE)
found     <- intersect(gene_ids, rownames(counts))
if (!length(found)) stop("none of the ", length(gene_ids),
                         " requested genes are rows of --counts", call. = FALSE)

cpm <- t(t(counts[found, samples, drop = FALSE]) / lib_sizes[samples]) * 1e6
# Low-coverage filter: a gene must reach min_cpm in at least min_samples samples,
# so one deep library cannot carry a hit on its own.
expressed <- rowSums(cpm >= min_cpm, na.rm = TRUE) >= min_samples
mat <- log2(cpm[expressed, , drop = FALSE] + 1)

cat(sprintf("genes requested %d | in counts %d | passing CPM>=%g in >=%d samples %d\n",
            length(gene_ids), length(found), min_cpm, min_samples, nrow(mat)))
cat(sprintf("samples %d across %d groups: %s\n", ncol(mat), nlevels(grp),
            paste(sprintf("%s=%d", levels(grp), table(grp)), collapse = ", ")))
if (!nrow(mat)) stop("no gene survived the coverage filter", call. = FALSE)

# ------------------------------------------------------------------- tests ---

#' Per-gene one-way ANOVA F-test across a grouping factor
#'
#' @param mat Numeric gene-by-sample matrix.
#' @param grp Factor of length ncol(mat) giving each sample's group.
#' @return Matrix with columns F and p, one row per gene (NA for constant rows).
anova_ftest <- function(mat, grp) {
  grp <- droplevels(grp)
  res <- t(vapply(seq_len(nrow(mat)), function(i) {
    y <- mat[i, ]
    if (length(unique(y)) < 2) return(c(F = NA_real_, p = NA_real_))
    s <- summary(aov(y ~ grp))[[1]]
    c(F = s[["F value"]][1], p = s[["Pr(>F)"]][1])
  }, numeric(2)))
  dimnames(res) <- list(rownames(mat), c("F", "p"))
  res
}

#' Per-gene limma moderated F-test (empirical-Bayes shrinkage, limma-trend)
#'
#' @param mat Numeric gene-by-sample matrix of log2 CPM.
#' @param grp Factor of length ncol(mat) giving each sample's group.
#' @return Matrix with columns F and p, one row per gene, in the order of mat.
limma_ftest <- function(mat, grp) {
  design <- model.matrix(~ droplevels(grp))
  fit    <- eBayes(lmFit(mat, design), trend = TRUE)
  tt     <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = "none")
  cat(sprintf("limma: eBayes prior df %.2f on %d residual df\n",
              fit$df.prior, fit$df.residual[1]))
  matrix(c(tt$F, tt$P.Value), ncol = 2,
         dimnames = list(rownames(mat), c("F", "p")))
}

#' Per-group row means (gene x group matrix)
#'
#' @param mat Numeric gene-by-sample matrix.
#' @param grp Factor of length ncol(mat) giving each sample's group.
#' @return Numeric gene-by-group matrix of means.
group_means <- function(mat, grp) {
  gl <- levels(droplevels(grp))
  m  <- vapply(gl, function(g) rowMeans(mat[, grp == g, drop = FALSE]), numeric(nrow(mat)))
  matrix(m, nrow = nrow(mat), dimnames = list(rownames(mat), gl))
}

af <- anova_ftest(mat, grp)
lf <- limma_ftest(mat, grp)
gm <- group_means(mat, grp)

res <- data.table(
  gene_id    = rownames(mat),
  F_anova    = round(af[, "F"], 3), p_anova = af[, "p"],
  F_limma    = round(lf[, "F"], 3), p_limma = lf[, "p"],
  peak_group = colnames(gm)[max.col(gm, ties.method = "first")])
res[, `:=`(fdr_anova = p.adjust(p_anova, "BH"),
           fdr_limma = p.adjust(p_limma, "BH"))]

# Optional per-gene annotation: everything except the ID column is carried
# through to the tables, and the flag/product columns feed the heatmap.
ann <- NULL
if (nzchar(opt$annotation)) {
  ann <- read_table_any(opt$annotation)
  if (!opt$ann_id_col %in% names(ann))
    stop("--annotation has no column '", opt$ann_id_col, "'; columns are: ",
         paste(names(ann), collapse = ", "), call. = FALSE)
  ann <- ann[!duplicated(ann[[opt$ann_id_col]])]
  setnames(ann, opt$ann_id_col, "gene_id")
  ann[, gene_id := as.character(gene_id)]
  res <- merge(res, ann, by = "gene_id", all.x = TRUE, sort = FALSE)
}

setorder(res, fdr_limma, p_limma)
res <- cbind(res, as.data.table(round(gm[res$gene_id, , drop = FALSE], 3)))

out <- function(suffix) file.path(opt$outdir, paste0(opt$prefix, suffix))
fwrite(res, out("_results.tsv"), sep = "\t")
fwrite(data.table(gene_id = rownames(mat), as.data.table(round(mat, 4))),
       out("_matrix_log2cpm.tsv"), sep = "\t")
fwrite(data.table(sample = samples, group = as.character(grp)),
       out("_sample_groups.tsv"), sep = "\t")

sig_anova <- res[!is.na(fdr_anova) & fdr_anova < fdr_cut][order(fdr_anova, p_anova)]
sig_limma <- res[!is.na(fdr_limma) & fdr_limma < fdr_cut][order(fdr_limma, p_limma)]
fwrite(sig_anova, out("_anova_significant.tsv"), sep = "\t")
fwrite(sig_limma, out("_limma_significant.tsv"), sep = "\t")

cat(sprintf("BH-FDR < %.2f : ANOVA %d, limma %d, shared %d\n", fdr_cut,
            nrow(sig_anova), nrow(sig_limma),
            length(intersect(sig_anova$gene_id, sig_limma$gene_id))))

# ---------------------------------------------------------------- heatmaps ---

# Fixed colours for the groups this project uses; anything else gets a
# reproducible fallback so the script still works on another gene/group set.
known_cols <- c(E.histolytica = "#4D4D4D", B.subtilis_E.histolytica = "#1B9E77",
                E.coli_E.histolytica = "#D95F02", E.faecalis_E.histolytica = "#7570B3",
                L.rhamnosus_E.histolytica = "#E7298A", MixSpp_E.histolytica = "#66A61E")
fallback  <- grDevices::hcl.colors(max(1L, nlevels(grp)), "Dark 3")
grp_cols  <- stats::setNames(ifelse(levels(grp) %in% names(known_cols),
                                    known_cols[levels(grp)], fallback), levels(grp))
flag_pal  <- c("FALSE" = "grey92", "TRUE" = "grey25")
flag_cols <- grep("Flag|flag", names(res), value = TRUE)
label_col <- intersect(c("Product Name", "product", "Product", "description"), names(res))[1]

#' Draw one heatmap of z-scored expression to PDF and SVG
#'
#' @param sig data.table of significant genes (must contain gene_id, peak_group).
#' @param title Panel title.
#' @param stem Output file stem, appended to the --prefix path.
#' @return Invisibly TRUE if a figure was drawn, FALSE if there was nothing to draw.
draw_heatmap <- function(sig, title, stem) {
  if (!nrow(sig)) { cat("nothing to draw for", stem, "\n"); return(invisible(FALSE)) }
  z <- t(scale(t(mat[sig$gene_id, , drop = FALSE]))); z[is.na(z)] <- 0

  top_ann <- HeatmapAnnotation(Group = grp, col = list(Group = grp_cols),
                               annotation_name_gp = gpar(fontsize = 8))
  ann_list <- c(
    list(`peak group` = factor(sig$peak_group, levels = levels(grp))),
    stats::setNames(lapply(flag_cols, function(f) toupper(as.character(sig[[f]]))),
                    flag_cols))
  col_list <- c(list(`peak group` = grp_cols),
                stats::setNames(rep(list(flag_pal), length(flag_cols)), flag_cols))
  row_ann <- do.call(rowAnnotation, c(ann_list, list(col = col_list,
    annotation_name_gp = gpar(fontsize = 7),
    show_legend = c(TRUE, rep(FALSE, length(flag_cols))))))

  lab <- if (!is.na(label_col)) paste0(sig$gene_id, "  ", sig[[label_col]]) else sig$gene_id
  lab <- ifelse(nchar(lab) > 48, paste0(substr(lab, 1, 45), "..."), lab)

  ht <- Heatmap(z, name = "z_log2CPM",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_columns = FALSE, clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    row_labels = lab, row_names_gp = gpar(fontsize = 6.5), show_column_names = FALSE,
    top_annotation = top_ann, left_annotation = row_ann,
    column_title = title, column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)))

  fig_h <- max(5, nrow(z) * 0.18 + 2.5)
  for (dev in c("pdf", "svg")) {
    path <- out(paste0(stem, ".", dev))
    if (dev == "pdf") grDevices::pdf(path, width = 8.5, height = fig_h)
    else              grDevices::svg(path, width = 8.5, height = fig_h)
    # Row clustering is deterministic, but seed anyway so reruns are byte-stable.
    set.seed(seed); draw(ht, merge_legend = TRUE); invisible(dev.off())
    cat("saved:", path, "\n")
  }
  # Persist the plotted values so the figure can be checked against a table.
  fwrite(data.table(gene_id = rownames(z), as.data.table(round(z, 4))),
         out(paste0(stem, "_zvalues.tsv")), sep = "\t")
  invisible(TRUE)
}

draw_heatmap(sig_anova, sprintf("%s — one-way ANOVA, BH-FDR < %.2f (%d genes, %d groups)",
                                opt$prefix, fdr_cut, nrow(sig_anova), nlevels(grp)),
             "_anova_heatmap")
draw_heatmap(sig_limma, sprintf("%s — limma moderated F, BH-FDR < %.2f (%d genes, %d groups)",
                                opt$prefix, fdr_cut, nrow(sig_limma), nlevels(grp)),
             "_limma_heatmap")

cat("done:", normalizePath(opt$outdir), "\n")
