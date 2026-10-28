# R/domain/multiomics/06c_rna_protein_pathway_correlation.R
# Per-pathway agreement of RNA and protein log2 fold changes.
#
# Reads the gene-protein log2FC pairs written by run_rna_protein_correlation()
# (06b) and, for every pathway in one configured GMT, correlates the RNA and
# protein log2FC of the pairs whose gene is a member. The report draws the
# interactive scatter from the two tables written here and recomputes nothing.

#' Resolve the GMT used for the per-pathway RNA-protein correlation
#'
#' The one place that reads
#' \code{modes$multiomics$rna_protein_correlation$pathway_gmt}. The view is off
#' unless the key names a file, so a project that never set it is unchanged.
#'
#' @param config Full pipeline config.
#' @return Absolute path to an existing GMT, or NULL when the key is unset or
#'   the file does not exist (the latter with a warning).
resolve_rna_protein_pathway_gmt <- function(config) {
    gmt <- config$modes$multiomics$rna_protein_correlation$pathway_gmt
    if (is.null(gmt) || length(gmt) == 0 || !nzchar(gmt[[1]])) return(NULL)
    if (length(gmt) > 1) {
        stop("modes$multiomics$rna_protein_correlation$pathway_gmt takes one GMT ",
             "path, got ", length(gmt), ". Give the single pathway collection ",
             "(e.g. a KEGG GMT) the scatter should offer.", call. = FALSE)
    }
    path <- resolve_input_path(config, gmt[[1]])
    if (!file.exists(path)) {
        warning("RNA-protein pathway correlation: pathway_gmt not found, skipped: ",
                path, call. = FALSE)
        return(NULL)
    }
    path
}

#' Resolve the smallest pathway size for the per-pathway RNA-protein correlation
#'
#' The one place that reads
#' \code{modes$multiomics$rna_protein_correlation$min_pairs}: how many mapped
#' gene-protein pairs a pathway needs in a contrast before it is scored.
#'
#' @param config Full pipeline config.
#' @return A whole number of at least 3; 11 (pathways with more than 10 pairs)
#'   when the key is unset.
resolve_rna_protein_pathway_min_pairs <- function(config) {
    n <- config$modes$multiomics$rna_protein_correlation$min_pairs
    if (is.null(n)) return(11L)
    # cor.test needs three points; below that a pathway has no Pearson p-value.
    if (length(n) != 1 || !is.numeric(n) || is.na(n) || n < 3 || n != round(n)) {
        stop("modes$multiomics$rna_protein_correlation$min_pairs must be one whole ",
             "number of at least 3 (the pairs a pathway needs to be scored); got ",
             paste(format(n), collapse = ", "), ". Leave it unset for the default 11.",
             call. = FALSE)
    }
    as.integer(n)
}

#' Correlate RNA and protein log2FC within each pathway
#'
#' Pairs join pathways on gene symbol, compared upper-cased. Each gene-protein
#' pair counts once, as it appears in the scatter, so a gene mapped to two
#' protein groups contributes two points. Pathways with fewer than
#' \code{min_pairs} finite pairs in a contrast are left out of that contrast.
#' P-values are per pathway and not adjusted across pathways.
#'
#' @param pairs Data frame with \code{gene_symbol}, \code{rna_log2FC},
#'   \code{protein_log2FC} and \code{contrast}; \code{gene_id} and
#'   \code{protein_id} are carried into the membership table when present.
#' @param gene_sets Named list of character vectors (pathway id -> symbols), as
#'   returned by \code{read_gmt()}; its \code{descriptions} attribute, when
#'   present, supplies the pathway names.
#' @param min_pairs Smallest number of pairs a pathway needs; the default 11
#'   keeps pathways with more than 10.
#' @return List with \code{stats} (one row per contrast x pathway: contrast,
#'   pathway_id, pathway_name, n_pairs, pearson_r, pearson_p, spearman_rho,
#'   spearman_p, same_sign_frac, min_pairs; ordered by contrast, then
#'   pearson_r from highest to lowest) and \code{members} (contrast,
#'   pathway_id, gene_symbol, gene_id, protein_id for the kept pathways).
#'   Both are zero-row data frames when no pathway qualifies.
#' @examples
#' pairs <- data.frame(gene_symbol = paste0("G", 1:60), contrast = "A_vs_B",
#'                     rna_log2FC = seq(-3, 3, length.out = 60),
#'                     protein_log2FC = seq(-1, 1, length.out = 60))
#' compute_pathway_rna_protein_correlation(pairs, list(P1 = paste0("G", 1:60)))
compute_pathway_rna_protein_correlation <- function(pairs, gene_sets, min_pairs = 11) {
    empty <- list(
        stats = data.frame(contrast = character(), pathway_id = character(),
                           pathway_name = character(), n_pairs = integer(),
                           pearson_r = numeric(), pearson_p = numeric(),
                           spearman_rho = numeric(), spearman_p = numeric(),
                           same_sign_frac = numeric(), min_pairs = integer(),
                           stringsAsFactors = FALSE),
        members = data.frame(contrast = character(), pathway_id = character(),
                             gene_symbol = character(), gene_id = character(),
                             protein_id = character(), stringsAsFactors = FALSE)
    )
    need <- c("gene_symbol", "rna_log2FC", "protein_log2FC", "contrast")
    if (is.null(pairs) || nrow(pairs) == 0 || !all(need %in% names(pairs)) ||
        length(gene_sets) == 0) {
        return(empty)
    }

    pairs <- pairs[is.finite(pairs$rna_log2FC) & is.finite(pairs$protein_log2FC) &
                       !is.na(pairs$gene_symbol) & nzchar(pairs$gene_symbol), , drop = FALSE]
    for (col in c("gene_id", "protein_id")) {
        if (!col %in% names(pairs)) pairs[[col]] <- NA_character_
    }
    key <- toupper(pairs$gene_symbol)
    descriptions <- attr(gene_sets, "descriptions")

    stats_rows <- list()
    member_rows <- list()
    for (cname in sort(unique(pairs$contrast))) {
        in_contrast <- pairs$contrast == cname
        for (pid in names(gene_sets)) {
            idx <- which(in_contrast & key %in% toupper(gene_sets[[pid]]))
            if (length(idx) < min_pairs) next
            x <- pairs$rna_log2FC[idx]
            y <- pairs$protein_log2FC[idx]
            pe <- suppressWarnings(stats::cor.test(x, y, method = "pearson"))
            sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
            nonzero <- x != 0 & y != 0
            nm <- if (!is.null(descriptions) && !is.na(descriptions[pid])) descriptions[[pid]] else pid
            stats_rows[[length(stats_rows) + 1]] <- data.frame(
                contrast = cname, pathway_id = pid, pathway_name = nm,
                n_pairs = length(idx),
                pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
                spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
                same_sign_frac = if (any(nonzero)) mean(sign(x[nonzero]) == sign(y[nonzero])) else NA_real_,
                min_pairs = as.integer(min_pairs),
                stringsAsFactors = FALSE
            )
            member_rows[[length(member_rows) + 1]] <- data.frame(
                contrast = cname, pathway_id = pid,
                gene_symbol = pairs$gene_symbol[idx],
                gene_id = as.character(pairs$gene_id[idx]),
                protein_id = as.character(pairs$protein_id[idx]),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(stats_rows) == 0) return(empty)

    st <- do.call(rbind, stats_rows)
    # pathway_id breaks ties so the dropdown order is the same on every run.
    st <- st[order(st$contrast, -st$pearson_r, st$pathway_id), , drop = FALSE]
    rownames(st) <- NULL
    mem <- do.call(rbind, member_rows)
    mem <- mem[order(mem$contrast, mem$pathway_id, mem$gene_symbol, mem$protein_id), , drop = FALSE]
    rownames(mem) <- NULL
    list(stats = st, members = mem)
}

#' Write the per-pathway RNA-protein correlation tables
#'
#' Side-effect step called by \code{run_rna_protein_correlation()} once the
#' pairs table exists. Writes \code{tables/rna_protein_pathway_correlation.csv}
#' and \code{tables/rna_protein_pathway_members.csv}; writes nothing when
#' \code{pathway_gmt} is unset or no pathway reaches \code{min_pairs}.
#'
#' @param pairs Gene-protein log2FC pairs (see
#'   \code{compute_pathway_rna_protein_correlation()}).
#' @param config Full pipeline config.
#' @param out_dir RNA-protein correlation output directory.
#' @return Invisibly the result of
#'   \code{compute_pathway_rna_protein_correlation()}, or NULL when skipped.
write_pathway_rna_protein_correlation <- function(pairs, config, out_dir) {
    gmt <- resolve_rna_protein_pathway_gmt(config)
    if (is.null(gmt)) {
        message("  Pathway RNA-protein correlation: no pathway_gmt set, skipped")
        return(invisible(NULL))
    }
    res <- compute_pathway_rna_protein_correlation(
        pairs, read_gmt(gmt), min_pairs = resolve_rna_protein_pathway_min_pairs(config))
    if (nrow(res$stats) == 0) {
        message("  Pathway RNA-protein correlation: no pathway has enough ",
                "mapped pairs in ", basename(gmt), ", nothing written")
        return(invisible(NULL))
    }
    res$stats$gmt_file <- basename(gmt)
    dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    write.csv(res$stats, file.path(out_dir, "tables", "rna_protein_pathway_correlation.csv"),
              row.names = FALSE)
    write.csv(res$members, file.path(out_dir, "tables", "rna_protein_pathway_members.csv"),
              row.names = FALSE)
    message("  Pathway RNA-protein correlation: ", nrow(res$stats),
            " contrast x pathway rows from ", basename(gmt))
    invisible(res)
}
