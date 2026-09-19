#' MultiGSEA Correlation Plots
#'
#' Functions for visualizing cross-omics enrichment correlations
#' and generating pathway overlay plots using pathview.
#'
#' @name multigsea_plots
NULL

#' Run MultiGSEA Correlation Analysis
#'
#' Generates scatter plots comparing enrichment scores between pairs of omics.
#'
#' @param enrichment_results List containing enrichment results from run_multiomics_enrichment.
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for plots and tables.
#' @return A list of ggplot objects.
#' @export
run_multigsea_plots <- function(enrichment_results, config, out_dir = NULL) {
    message("=== Running MultiGSEA Correlation Analysis ===")

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available for MultiGSEA.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_multigsea %||% TRUE)) {
        message("MultiGSEA analysis disabled in config.")
        return(NULL)
    }

    p_thresh <- mg_config$pvalue_threshold %||% 0.05
    corr_method <- mg_config$correlation_method %||% "pearson"

    # Resolved once and passed down explicitly. Every consumer of a term id below
    # sits inside this function or is called from it, so nothing has to reach for
    # config again -- and the pairwise panels, the combined panel and the
    # per-contrast plots all key pathways the same way as a result.
    kegg_org <- resolve_kegg_org_code(config$global$organism)

    # Create output directory
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Extract results per omic
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    if (length(omics_names) < 2) {
        message("Need at least 2 omics with enrichment results for MultiGSEA.")
        return(NULL)
    }

    # Identify pairs
    pairs <- utils::combn(omics_names, 2, simplify = FALSE)
    plots <- list()

    for (pair in pairs) {
        omic1 <- pair[1]
        omic2 <- pair[2]

        res1 <- per_omics[[omic1]]
        res2 <- per_omics[[omic2]]

        if (is.null(res1) || is.null(res2)) next

        term1 <- .multigsea_term_ids(res1, kegg_org)
        term2 <- .multigsea_term_ids(res2, kegg_org)
        if (is.null(term1) || is.null(term2)) {
            message("No term ID column in ", omic1, " or ", omic2,
                    " enrichment results; skipping this pair")
            next
        }
        # Names first, from the original frames: the collapse below keeps one row
        # per pathway by p-value, and that row is not necessarily the one carrying
        # the readable name. Scoring picks a representative; display should still
        # see every row that could name the pathway. Built before `term` is
        # overwritten, too, so a frame whose own `term` carries its name still has
        # that text to offer.
        id_to_name <- .multigsea_term_names(list(res1, res2), kegg_org)

        # A row whose every identifier column is blank has no key at all. It
        # cannot be matched or labelled, and carrying it into common_terms hands
        # resolve_term() an NA, where nchar(NA_character_) is NA and the `if`
        # below it stops the whole pairwise loop.
        res1 <- res1[!is.na(term1), , drop = FALSE]; term1 <- term1[!is.na(term1)]
        res2 <- res2[!is.na(term2), , drop = FALSE]; term2 <- term2[!is.na(term2)]

        # An omic left with nothing keyed has no enrichment identity to correlate.
        # The union below cannot see that -- the other side alone can carry it past
        # the length check -- and the pair would be plotted against a column of
        # imputed zeros.
        if (length(term1) == 0L || length(term2) == 0L) {
            message("No usable pathway identity in ", omic1, " or ", omic2,
                    "; skipping this pair")
            next
        }

        # One row per pathway per omic before anything matches on it: match()
        # below would otherwise resolve a duplicated key by row order.
        keep1 <- .multigsea_collapse_duplicate_terms(res1, term1)
        keep2 <- .multigsea_collapse_duplicate_terms(res2, term2)
        res1 <- res1[keep1, , drop = FALSE]; term1 <- term1[keep1]
        res2 <- res2[keep2, , drop = FALSE]; term2 <- term2[keep2]

        res1$term <- term1
        res2$term <- term2

        # Union of terms
        common_terms <- union(res1$term, res2$term)

        if (length(common_terms) < 3) {
            message("Too few terms (union) between ", omic1, " and ", omic2)
            next
        }

        # Align data
        df1 <- res1[match(common_terms, res1$term), ]
        df2 <- res2[match(common_terms, res2$term), ]

        # Calculate -log10(padj) scores
        get_score <- function(df) {
            # Robustly find p-adj column
            padj_col <- NULL
            for (col in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
                if (col %in% colnames(df)) {
                    padj_col <- col
                    break
                }
            }

            if (is.null(padj_col)) {
                warning("No p-value column found in enrichment results")
                return(rep(0, nrow(df)))
            }

            padj <- df[[padj_col]]
            padj[is.na(padj)] <- 1

            # Handle zero p-values
            non_zeros <- padj[padj > 0 & padj < 1]
            if (length(non_zeros) > 0) {
                min_nz <- min(non_zeros, na.rm = TRUE)
            } else {
                min_nz <- 1e-10
            }

            padj[padj == 0] <- min_nz / 10
            -log10(padj)
        }

        score1 <- get_score(df1)
        score2 <- get_score(df2)

        # Extract count (setSize) and compute fold enrichment from GeneRatio
        get_count <- function(df) {
            if ("setSize" %in% colnames(df)) return(df$setSize)
            if ("Count" %in% colnames(df)) return(df$Count)
            return(rep(NA_real_, nrow(df)))
        }

        parse_gene_ratio <- function(df) {
            if ("Fold Enrichment" %in% colnames(df)) return(df[["Fold Enrichment"]])
            if ("fold_enrichment" %in% colnames(df)) return(df$fold_enrichment)
            if ("GeneRatio" %in% colnames(df)) {
                gr <- as.character(df$GeneRatio)
                parts <- strsplit(gr, "/")
                ratio <- vapply(parts, function(p) {
                    if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2])
                    else NA_real_
                }, numeric(1))
                return(ratio)
            }
            return(rep(NA_real_, nrow(df)))
        }

        count1 <- get_count(df1)
        count2 <- get_count(df2)
        fold1 <- parse_gene_ratio(df1)
        fold2 <- parse_gene_ratio(df2)

        plot_df <- data.frame(
            term = common_terms,
            x = score1,
            y = score2,
            count1 = count1,
            count2 = count2,
            fold1 = fold1,
            fold2 = fold2,
            stringsAsFactors = FALSE
        )

        # Impute NAs (terms present in one omic but not the other)
        plot_df$count1[is.na(plot_df$count1)] <- 0
        plot_df$count2[is.na(plot_df$count2)] <- 0
        plot_df$fold1[is.na(plot_df$fold1)] <- 0
        plot_df$fold2[is.na(plot_df$fold2)] <- 0

        plot_df$avg_count <- (plot_df$count1 + plot_df$count2) / 2
        plot_df$avg_fold <- (plot_df$fold1 + plot_df$fold2) / 2

        # Resolve term to readable name, truncate long names
        resolve_term <- function(term) {
            # Look up pathway name from ID
            if (term %in% names(id_to_name) && nzchar(id_to_name[[term]])) {
                t <- id_to_name[[term]]
            } else {
                # Same stripping rule as the name map synthesizes with. Two
                # implementations of it is what let the map shadow this fallback
                # with a worse label.
                t <- .multigsea_readable_from_identifier(term)
            }
            if (nchar(t) > 50) t <- paste0(substr(t, 1, 47), "...")
            t
        }
        plot_df$label <- vapply(plot_df$term, resolve_term, character(1))

        # Calculate correlation
        cor_res <- cor.test(plot_df$x, plot_df$y, method = corr_method)
        cor_val <- round(cor_res$estimate, 3)
        p_val <- signif(cor_res$p.value, 3)

        omic1_label <- gsub("_", " ", tools::toTitleCase(omic1))
        omic2_label <- gsub("_", " ", tools::toTitleCase(omic2))

        # Determine whether we have meaningful count/fold data
        has_count <- any(plot_df$avg_count > 0, na.rm = TRUE)
        has_fold <- any(plot_df$avg_fold > 0, na.rm = TRUE)

        # Create plot matching DAVID scatter style
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y))

        if (has_count && has_fold) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(size = avg_count, color = avg_fold), alpha = 0.7
            ) +
            ggplot2::scale_color_viridis_c(name = "Avg Fold Enrichment") +
            ggplot2::scale_size_continuous(name = "Avg Count")
        } else if (has_count) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(size = avg_count), color = "steelblue", alpha = 0.7
            ) +
            ggplot2::scale_size_continuous(name = "Avg Count")
        } else if (has_fold) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(color = avg_fold), size = 3, alpha = 0.7
            ) +
            ggplot2::scale_color_viridis_c(name = "Avg Fold Enrichment")
        } else {
            p <- p + ggplot2::geom_point(color = "steelblue", size = 3, alpha = 0.7)
        }

        p <- p +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::labs(
                title = paste0(omic1_label, " vs ", omic2_label),
                subtitle = paste0("Total Unique Terms: ", nrow(plot_df),
                                  "  |  ", corr_method, " r = ", cor_val, ", p = ", p_val),
                x = paste0(omic1_label, " [-log10(FDR)]"),
                y = paste0(omic2_label, " [-log10(FDR)]")
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.title = ggplot2::element_text(face = "bold"),
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5)
            )

        # Label terms with ggrepel — top terms from each axis + any with both
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            cut_score <- -log10(p_thresh)
            both_sig <- plot_df[plot_df$x > cut_score & plot_df$y > cut_score, ]
            top_x <- head(plot_df[order(plot_df$x, decreasing = TRUE), ], 10)
            top_y <- head(plot_df[order(plot_df$y, decreasing = TRUE), ], 5)
            label_df <- unique(rbind(both_sig, top_x, top_y))
            label_df <- label_df[label_df$x > 0 | label_df$y > 0, ]

            p <- p + ggrepel::geom_text_repel(
                data = label_df,
                ggplot2::aes(label = label),
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.5,
                force = 2,
                min.segment.length = 0
            )
        }

        # Save plot and data
        if (!is.null(out_dir)) {
            filename <- paste0("multigsea_", omic1, "_vs_", omic2)
            ggplot2::ggsave(
                file.path(out_dir, paste0(filename, ".png")),
                plot = p, width = 8, height = 8, dpi = 300
            )
            ggplot2::ggsave(
                file.path(out_dir, paste0(filename, ".pdf")),
                plot = p, width = 8, height = 8
            )
            write.csv(plot_df, file.path(out_dir, paste0(filename, ".csv")),
                      row.names = FALSE)
        }

        plots[[paste0(omic1, "_vs_", omic2)]] <- p
    }

    # Generate combined 2x2 plot
    if (length(plots) > 0 && !is.null(out_dir)) {
        combined <- plot_multigsea_combined(plots, per_omics, out_dir, kegg_org)
        if (!is.null(combined)) {
            plots[["combined"]] <- combined
        }
    }

    # --- Per-contrast MultiGSEA plots ---
    contrast_names_mg <- unique(unlist(lapply(per_omics, function(df) {
        if (is.data.frame(df) && "contrast" %in% colnames(df)) unique(df$contrast)
        else NULL
    })))

    if (length(contrast_names_mg) > 1 && !is.null(out_dir)) {
        message("  Generating per-contrast MultiGSEA plots for ",
                length(contrast_names_mg), " contrasts")
        per_contrast_dir <- file.path(out_dir, "per_contrast")

        for (cname in contrast_names_mg) {
            safe_dir <- gsub("[^a-zA-Z0-9._-]", "_", cname)
            contrast_out <- file.path(per_contrast_dir, safe_dir)
            dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

            per_omics_contrast <- list()
            for (om in omics_names) {
                df <- per_omics[[om]]
                if (is.data.frame(df) && "contrast" %in% colnames(df)) {
                    df_c <- df[df$contrast == cname, , drop = FALSE]
                } else {
                    df_c <- df
                }
                if (is.data.frame(df_c) && nrow(df_c) > 0) {
                    per_omics_contrast[[om]] <- df_c
                }
            }

            if (length(per_omics_contrast) < 2) next

            contrast_pairs <- utils::combn(names(per_omics_contrast), 2, simplify = FALSE)
            for (pair in contrast_pairs) {
                omic1 <- pair[1]
                omic2 <- pair[2]
                res1 <- per_omics_contrast[[omic1]]
                res2 <- per_omics_contrast[[omic2]]
                if (is.null(res1) || is.null(res2)) next

                tryCatch({
                    .save_multigsea_pair_plot(
                        res1, res2, omic1, omic2,
                        corr_method = corr_method,
                        p_thresh = p_thresh,
                        out_dir = contrast_out,
                        kegg_org = kegg_org
                    )
                }, error = function(e) {
                    message("    MultiGSEA plot failed for ", omic1, " vs ", omic2,
                            " (", cname, "): ", e$message)
                })
            }
        }
        message("  Per-contrast MultiGSEA output saved to: ", per_contrast_dir)
    }

    message("MultiGSEA plots generated: ", length(plots))
    return(plots)
}


#' Save a single MultiGSEA pairwise scatter plot
#'
#' Generates and saves the enrichment correlation scatter plot for one pair
#' of omics layers. Used by both the combined and per-contrast MultiGSEA code.
#'
#' @param res1 Enrichment data frame for omics 1
#' @param res2 Enrichment data frame for omics 2
#' @param omic1 Name of omics 1
#' @param omic2 Name of omics 2
#' @param corr_method Correlation method (default "pearson")
#' @param p_thresh P-value threshold for labeling (default 0.05)
#' @param out_dir Output directory for saved files
#' @param kegg_org Active KEGG organism code for the run, or NULL. Passed in
#'   rather than read from config so this function has no hidden dependency.
#' @return Invisible NULL
.save_multigsea_pair_plot <- function(res1, res2, omic1, omic2,
                                      corr_method = "pearson",
                                      p_thresh = 0.05, out_dir,
                                      kegg_org = NULL) {

    term1 <- .multigsea_term_ids(res1, kegg_org)
    term2 <- .multigsea_term_ids(res2, kegg_org)
    if (is.null(term1) || is.null(term2)) return(invisible(NULL))

    # Names first, from the original frames and before `term` is overwritten, as
    # in the pairwise loop above.
    id_to_name <- .multigsea_term_names(list(res1, res2), kegg_org)

    # Rows with no identity at all cannot be matched or labelled, as above.
    res1 <- res1[!is.na(term1), , drop = FALSE]; term1 <- term1[!is.na(term1)]
    res2 <- res2[!is.na(term2), , drop = FALSE]; term2 <- term2[!is.na(term2)]

    # An omic left with nothing keyed has no enrichment identity to correlate,
    # and the union below cannot see that -- the other side alone carries it
    # past the length check, against a column of imputed zeros.
    if (length(term1) == 0L || length(term2) == 0L) return(invisible(NULL))

    # One row per pathway per omic, as in the pairwise loop above.
    keep1 <- .multigsea_collapse_duplicate_terms(res1, term1)
    keep2 <- .multigsea_collapse_duplicate_terms(res2, term2)
    res1 <- res1[keep1, , drop = FALSE]; term1 <- term1[keep1]
    res2 <- res2[keep2, , drop = FALSE]; term2 <- term2[keep2]

    res1$term <- term1
    res2$term <- term2

    common_terms <- union(res1$term, res2$term)
    if (length(common_terms) < 3) return(invisible(NULL))

    df1 <- res1[match(common_terms, res1$term), ]
    df2 <- res2[match(common_terms, res2$term), ]

    get_score <- function(df) {
        padj_col <- NULL
        for (col in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
            if (col %in% colnames(df)) { padj_col <- col; break }
        }
        if (is.null(padj_col)) return(rep(0, nrow(df)))
        padj <- df[[padj_col]]
        padj[is.na(padj)] <- 1
        non_zeros <- padj[padj > 0 & padj < 1]
        min_nz <- if (length(non_zeros) > 0) min(non_zeros, na.rm = TRUE) else 1e-10
        padj[padj == 0] <- min_nz / 10
        -log10(padj)
    }

    score1 <- get_score(df1)
    score2 <- get_score(df2)

    resolve_term <- function(term) {
        if (term %in% names(id_to_name) && nzchar(id_to_name[[term]])) {
            t <- id_to_name[[term]]
        } else {
            # Same stripping rule as the name map synthesizes with, as above.
            t <- .multigsea_readable_from_identifier(term)
        }
        if (nchar(t) > 50) t <- paste0(substr(t, 1, 47), "...")
        t
    }

    plot_df <- data.frame(
        term = common_terms,
        x = score1, y = score2,
        stringsAsFactors = FALSE
    )
    plot_df$label <- vapply(plot_df$term, resolve_term, character(1))

    cor_res <- cor.test(plot_df$x, plot_df$y, method = corr_method)
    cor_val <- round(cor_res$estimate, 3)
    p_val <- signif(cor_res$p.value, 3)

    omic1_label <- gsub("_", " ", tools::toTitleCase(omic1))
    omic2_label <- gsub("_", " ", tools::toTitleCase(omic2))

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(color = "steelblue", size = 3, alpha = 0.7) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::labs(
            title = paste0(omic1_label, " vs ", omic2_label),
            subtitle = paste0("Terms: ", nrow(plot_df),
                              "  |  ", corr_method, " r = ", cor_val, ", p = ", p_val),
            x = paste0(omic1_label, " [-log10(FDR)]"),
            y = paste0(omic2_label, " [-log10(FDR)]")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(face = "bold"),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5)
        )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
        cut_score <- -log10(p_thresh)
        both_sig <- plot_df[plot_df$x > cut_score & plot_df$y > cut_score, ]
        top_x <- utils::head(plot_df[order(plot_df$x, decreasing = TRUE), ], 10)
        top_y <- utils::head(plot_df[order(plot_df$y, decreasing = TRUE), ], 5)
        label_df <- unique(rbind(both_sig, top_x, top_y))
        label_df <- label_df[label_df$x > 0 | label_df$y > 0, ]
        p <- p + ggrepel::geom_text_repel(
            data = label_df, ggplot2::aes(label = label),
            size = 3, max.overlaps = Inf, box.padding = 0.5,
            force = 2, min.segment.length = 0
        )
    }

    filename <- paste0("multigsea_", omic1, "_vs_", omic2)
    ggplot2::ggsave(file.path(out_dir, paste0(filename, ".png")),
                    plot = p, width = 8, height = 8, dpi = 300)
    write.csv(plot_df, file.path(out_dir, paste0(filename, ".csv")),
              row.names = FALSE)

    invisible(NULL)
}


#' Term IDs of a per-omic enrichment table
#'
#' The per-omic tables key each gene set on `pathway` (e.g. "GO:0000027") and
#' keep the readable label in `pathway_name`; clusterProfiler-style tables use
#' `ID` and `Description`. Matching two omics on the stable ID keeps them
#' aligned. There is no row-name fallback on purpose: bound tables carry
#' positional row names, which would pair unrelated pathways by row number.
#'
#' Normalization is delegated to \code{normalize_pathway_join_key()}, so the two
#' omics are matched on the same normalized key the cross-omics join uses: the
#' KEGG forms hsa00010, map00010, ko00010 and 00010 are one pathway, while GO,
#' PFAM, InterPro and custom gene-set names come through byte-identical.
#'
#' @param df Enrichment data frame for one omic.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Character vector of term IDs, one per row, or NULL when the table has
#'   no recognised ID column.
.multigsea_term_ids <- function(df, kegg_org = NULL) {
    id <- .multigsea_identity(df, kegg_org)
    if (is.null(id)) NULL else id$key
}


#' Where each row's identity came from, and what it was
#'
#' One object, resolved once, carrying the key **and** the decision that produced
#' it. Everything downstream reads this rather than inferring from the columns
#' again: a display fallback that re-derives "did this row key on its ID?" gets it
#' wrong whenever the answer is subtler than the column being present, and after
#' normalization the raw text it needed is no longer recoverable from the key.
#'
#' The ladder is `term`, then \code{pathway_join_key()}'s `ID` -> `pathway` ->
#' `Description`. `term` is MultiGSEA's own: the pairwise callers write their
#' resolved identity onto it, so it precedes the shared ladder by design.
#'
#' Trimming decides only whether a candidate is usable. `raw` is the original
#' selected value, whitespace and all, and `key` is that value passed through
#' \code{normalize_pathway_join_key()} -- the KEGG rule is not reimplemented here.
#'
#' @param df Enrichment data frame for one omic.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return List of three aligned character vectors -- `key`, `raw` and `source`
#'   (one of "term", "ID", "pathway", "Description", or NA for a row where none
#'   carried a value) -- or NULL when the table has no candidate column at all.
#' @keywords internal
.multigsea_identity <- function(df, kegg_org = NULL) {
    candidates <- c("term", "ID", "pathway", "Description")
    present <- intersect(candidates, colnames(df))
    if (length(present) == 0L) return(NULL)

    n <- nrow(df)
    raw <- rep(NA_character_, n)
    source <- rep(NA_character_, n)

    for (col in present) {
        v <- as.character(df[[col]])
        fill <- is.na(raw) & .multigsea_usable_label(v)
        raw[fill] <- v[fill]
        source[fill] <- col
    }

    list(key = normalize_pathway_join_key(raw, kegg_org), raw = raw, source = source)
}


#' The column MultiGSEA scores a table on
#'
#' The same preference order the scoring and summary paths already use, in one
#' place so a collapse cannot reduce on a different column than the one plotted.
#'
#' @param df Enrichment data frame for one omic.
#' @return Column name, or NULL when the table carries none of them.
#' @keywords internal
.multigsea_padj_col <- function(df) {
    for (col in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
        if (col %in% colnames(df)) return(col)
    }
    NULL
}


#' One row per normalized term within one omic
#'
#' Normalizing the identity makes hsa00010 and map00010 the same pathway, so a
#' bound per-omic table can now hold two rows for it. Left alone those duplicates
#' make \code{match()} pick whichever row came first -- row order standing in for
#' biology -- and let one omic satisfy a co-significance test on its own.
#'
#' The reduction rule is the one #201 settled on for exactly this collapse:
#' \code{merge_pathway_pvalues()} aggregates a normalized key with \code{FUN =
#' min}. Nothing is recomputed or re-adjusted; the most significant of rows
#' already declared the same pathway is kept.
#'
#' @param df Enrichment data frame for one omic.
#' @param terms Normalized term ids aligned to the rows of \code{df}.
#' @return Integer row indices to keep, in their original order. Rows whose term
#'   is NA are all kept: they have no key to be duplicates of.
#' @keywords internal
.multigsea_collapse_duplicate_terms <- function(df, terms) {
    n <- length(terms)
    if (n == 0L) return(integer(0))

    padj_col <- .multigsea_padj_col(df)
    p <- if (is.null(padj_col)) {
        rep(NA_real_, n)
    } else {
        v <- df[[padj_col]]
        if (!is.numeric(v)) v <- suppressWarnings(as.numeric(as.character(v)))
        v
    }

    # Sorting by term then p puts the most significant row of each term first;
    # na.last keeps a row with no p-value behind one that has one.
    ord <- order(terms, p, na.last = TRUE)
    keep <- ord[!duplicated(terms[ord]) | is.na(terms[ord])]
    sort(unique(keep))
}


#' Terms significant in at least two distinct omics
#'
#' Counted over distinct omics rather than rows. "Co-significant" has to mean two
#' layers agreed, and normalizing the identity lets one layer hold several rows
#' for one pathway -- counting rows would let a single omic satisfy it alone.
#'
#' @param summary_df Data frame with \code{term}, \code{omic} and \code{padj}.
#' @param alpha Significance threshold, unchanged from the panel's own.
#' @return Character vector of terms, possibly empty.
#' @keywords internal
.multigsea_cosignificant_terms <- function(summary_df, alpha = 0.05) {
    sig <- summary_df[!is.na(summary_df$padj) & summary_df$padj < alpha, ,
                      drop = FALSE]
    if (nrow(sig) == 0L) return(character(0))
    counts <- table(unique(sig[, c("term", "omic")])$term)
    names(counts[counts >= 2])
}


#' Does this label actually say anything?
#'
#' @param x Character vector of candidate labels.
#' @return Logical vector; FALSE for NA, empty and whitespace-only labels.
#' @keywords internal
.multigsea_usable_label <- function(x) {
    x <- as.character(x)
    !is.na(x) & nzchar(trimws(x))
}


#' Readable text of an identifier that carries its own name
#'
#' The non-model KEGG fallback names a gene set "<accession> <readable name>"
#' in a single string. \code{resolve_term()} used to recover the readable part by
#' stripping the accession off the term it was given; normalizing the key removes
#' the prefix it stripped, so the name map has to carry that remainder instead or
#' the plot label would regress to a bare map number.
#'
#' It reproduces \code{resolve_term()}'s fallback exactly -- both of its strips,
#' and its rule that a strip leaving nothing behind falls back to the identifier
#' itself. That is the point: whatever this returns is installed in the name map,
#' which shadows the fallback, so anything it does not strip is text
#' \code{resolve_term()} would have removed and the label would regress.
#'
#' Deliberately organism-agnostic, unlike the identity side. `resolve_term()`
#' strips any two- or three-letter prefix, so a `mmu#####` label on a human run
#' is still shortened for display even though `mmu00010` is correctly *not* the
#' same pathway as `hsa00010` for joining. Display and identity answer different
#' questions here.
#'
#' @param ids Character vector of identifiers.
#' @return Character vector the same length as \code{ids}.
#' @keywords internal
.multigsea_readable_from_identifier <- function(ids) {
    ids <- as.character(ids)
    stripped <- sub("^GO:[0-9]+~", "", trimws(ids))
    stripped <- sub("^[a-z]{2,3}[0-9]{5}[[:space:]]*", "", stripped)
    ifelse(!is.na(stripped) & nzchar(stripped), stripped, ids)
}


#' Readable names for enrichment term IDs
#'
#' Driven by \code{.multigsea_identity()} rather than by re-inspecting which
#' column probably supplied the key. Asking the columns again is what put wrong
#' labels on rows: it cannot see `term`, it cannot see the raw text once the key
#' is normalized, and a column being present is not the same as that row having
#' keyed on it.
#'
#' @param dfs List of per-omic enrichment data frames.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Named character vector mapping term ID to readable name, keeping the
#'   first name seen for each ID; empty when no table carries a name column.
.multigsea_term_names <- function(dfs, kegg_org = NULL) {
    all_keys <- character(0)
    all_vals <- character(0)
    all_explicit <- logical(0)

    for (df in dfs) {
        if (!is.data.frame(df)) next
        id <- .multigsea_identity(df, kegg_org)
        if (is.null(id)) next
        keys <- id$key
        n <- length(keys)

        # Resolved per row throughout: a column filled for only some rows must not
        # leave the rest permanently blank just because the column exists.

        # 1. An explicit name column wins.
        vals <- if ("pathway_name" %in% colnames(df)) {
            as.character(df$pathway_name)
        } else {
            rep(NA_character_, n)
        }

        col_or_na <- function(nm) {
            if (nm %in% colnames(df)) as.character(df[[nm]]) else rep(NA_character_, n)
        }

        # 2. A row that keyed on its ID has the readable label in a column beside
        #    it: the older shape puts it in `pathway`, clusterProfiler in
        #    `Description`. Only rows that actually keyed on ID take this.
        keyed_on_id <- !is.na(id$source) & id$source == "ID"
        beside_id <- col_or_na("pathway")
        desc_v <- col_or_na("Description")
        fill <- !.multigsea_usable_label(beside_id)
        beside_id[fill] <- desc_v[fill]
        beside_id[!keyed_on_id] <- NA_character_

        gap <- !.multigsea_usable_label(vals)
        vals[gap] <- beside_id[gap]

        # Rungs 1 and 2 are a name the table actually states. Rung 3 below only
        # re-reads the identifier, so it must never outrank a stated name.
        explicit <- .multigsea_usable_label(vals)

        # 3. Anything still unnamed falls back to the readable text of the very
        #    value that produced its key, whichever column that was. That covers a
        #    `term` carrying its own name, an ID-only table whose accession would
        #    otherwise be lost to normalization, and the pathway/Description rungs.
        gap <- !.multigsea_usable_label(vals)
        vals[gap] <- .multigsea_readable_from_identifier(id$raw)[gap]

        # A label with nothing in it must not reserve a key. hsa00010 and
        # map00010 collapse to one key now, so a blank name in the layer that
        # happens to be visited first would otherwise lock out a readable name in
        # the next one, and the panel would fall back to the bare map number.
        keep <- .multigsea_usable_label(vals) & !is.na(keys)
        all_keys <- c(all_keys, keys[keep])
        all_vals <- c(all_vals, vals[keep])
        all_explicit <- c(all_explicit, explicit[keep])
    }

    if (length(all_keys) == 0L) return(character(0))

    # Resolved across every frame at once rather than frame by frame, so one
    # normalized key gets one label: two omics holding bare hsa00010 and map00010
    # must not name the same pathway differently.
    #
    # A name the table states outranks one synthesized from an identifier, whether
    # or not that identifier looks like an accession -- a bare "GO:0006915" is not
    # a KEGG accession but it is not a name either, and before this map existed it
    # contributed nothing and the stated name won. Among synthesized labels an
    # accession still loses to anything else, and first seen breaks the remaining
    # ties -- seq_along does that explicitly rather than relying on a stable sort.
    readable <- !is_kegg_pathway_accession(all_vals, kegg_org)
    ord <- order(!all_explicit, !readable, seq_along(all_vals))

    id_to_name <- setNames(all_vals[ord], all_keys[ord])
    id_to_name[!duplicated(names(id_to_name))]
}


#' Combined MultiGSEA Enrichment Plot
#'
#' Creates a 2x2 grid combining pairwise scatter plots with a summary dot plot.
#'
#' @param pairwise_plots Named list of pairwise ggplot objects.
#' @param per_omics Per-omics enrichment results list.
#' @param out_dir Output directory.
#' @param kegg_org Active KEGG organism code for the run, or NULL. Passed in
#'   rather than read from config so this function has no hidden dependency.
#' @return The combined ggplot/patchwork object, or NULL on failure.
plot_multigsea_combined <- function(pairwise_plots, per_omics, out_dir,
                                     kegg_org = NULL) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        message("Package 'patchwork' not available. Skipping combined plot.")
        return(NULL)
    }

    # Build summary dot plot of top co-significant pathways across omics
    omics_names <- names(per_omics)
    summary_rows <- list()

    # Resolved once over every participating omic, not per omic: a normalized key
    # must carry one label, or two layers holding bare hsa00010 and map00010 would
    # label the same pathway differently and take a y-axis row each.
    id_to_name <- .multigsea_term_names(per_omics, kegg_org)

    for (om in omics_names) {
        res <- per_omics[[om]]
        if (is.null(res) || !is.data.frame(res)) next

        term_ids <- .multigsea_term_ids(res, kegg_org)
        if (is.null(term_ids)) next

        padj_col <- .multigsea_padj_col(res)
        if (is.null(padj_col)) next

        # Rows with no identity cannot be placed on the panel, as in the pairwise
        # paths above.
        res <- res[!is.na(term_ids), , drop = FALSE]
        term_ids <- term_ids[!is.na(term_ids)]
        if (length(term_ids) == 0L) next

        # Readable name where any omic has one, otherwise the ID itself
        term_names <- unname(id_to_name[term_ids])
        term_names[is.na(term_names)] <- term_ids[is.na(term_names)]

        df_tmp <- data.frame(
            term = term_ids,
            pathway_name = term_names,
            padj = res[[padj_col]],
            omic = om,
            stringsAsFactors = FALSE
        )
        df_tmp <- df_tmp[!is.na(df_tmp$padj), ]
        # One row per pathway per omic, so the co-significance count below cannot
        # read two prefix variants from this layer as two layers.
        df_tmp <- df_tmp[.multigsea_collapse_duplicate_terms(df_tmp, df_tmp$term), ,
                         drop = FALSE]
        summary_rows[[length(summary_rows) + 1]] <- df_tmp
    }

    if (length(summary_rows) == 0) {
        message("No data for summary panel.")
        return(NULL)
    }

    summary_df <- do.call(rbind, summary_rows)

    co_sig <- .multigsea_cosignificant_terms(summary_df)

    if (length(co_sig) == 0) {
        # Fall back to top terms by lowest p-value across any omic
        top_terms <- unique(summary_df$term[order(summary_df$padj)])[1:min(15, length(unique(summary_df$term)))]
        plot_data <- summary_df[summary_df$term %in% top_terms, ]
        panel_title <- "Top Enriched Pathways"
    } else {
        plot_data <- summary_df[summary_df$term %in% co_sig, ]
        # Keep top 15 by mean -log10(padj) across omics
        mean_scores <- tapply(-log10(pmax(plot_data$padj, 1e-300)), plot_data$term, mean, na.rm = TRUE)
        top_terms <- names(sort(mean_scores, decreasing = TRUE))[1:min(15, length(mean_scores))]
        plot_data <- plot_data[plot_data$term %in% top_terms, ]
        panel_title <- "Co-Significant Pathways"
    }

    plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, 1e-300))
    # Use pathway name, truncate long names
    plot_data$term_short <- ifelse(
        !is.na(plot_data$pathway_name) & nzchar(plot_data$pathway_name),
        plot_data$pathway_name,
        plot_data$term
    )
    plot_data$term_short <- ifelse(
        nchar(plot_data$term_short) > 40,
        paste0(substr(plot_data$term_short, 1, 37), "..."),
        plot_data$term_short
    )
    plot_data$omic_label <- gsub("_", " ", tools::toTitleCase(plot_data$omic))

    summary_panel <- ggplot2::ggplot(plot_data,
        ggplot2::aes(x = omic_label, y = stats::reorder(term_short, neg_log10_padj),
                     size = neg_log10_padj, color = neg_log10_padj)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_viridis_c(name = "-log10(padj)") +
        ggplot2::scale_size_continuous(name = "-log10(padj)", range = c(2, 8)) +
        ggplot2::labs(title = panel_title, x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.title = ggplot2::element_text(face = "bold")
        )

    # Assemble 2x2 grid: up to 3 pairwise + summary panel
    scatter_plots <- pairwise_plots[!names(pairwise_plots) %in% "combined"]
    panels <- scatter_plots[1:min(3, length(scatter_plots))]
    panels[[length(panels) + 1]] <- summary_panel

    combined <- patchwork::wrap_plots(panels, ncol = 2)

    ggplot2::ggsave(
        file.path(out_dir, "multigsea_combined_enrichment.png"),
        plot = combined, width = 16, height = 16, dpi = 300
    )
    message("Saved combined MultiGSEA plot.")

    return(combined)
}


#' Run Pathview Visualization for Agreed Pathways
#'
#' Generates KEGG pathway overlay plots for pathways significant in multiple omics.
#'
#' @param enrichment_results List containing enrichment results.
#' @param mae_data MultiAssayExperiment object (or list with harmonized_omics).
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for pathview plots.
#' @return List of generated plot paths.
#' @export
run_multigsea_pathview <- function(enrichment_results, mae_data, config, out_dir = NULL) {
    message("=== Running MultiGSEA Pathview Visualization ===")

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("Package 'pathview' not installed. Skipping.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_pathview %||% TRUE)) {
        message("Pathview analysis disabled in config.")
        return(NULL)
    }

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available.")
        return(NULL)
    }

    # Create output directory
    if (is.null(out_dir)) {
        out_dir <- tempdir()
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Identify Agreed KEGG Pathways
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    kegg_pathways <- list()

    for (omic in omics_names) {
        res <- per_omics[[omic]]
        if (!is.null(res) && nrow(res) > 0) {
            # Check for KEGG results
            if ("type" %in% colnames(res)) {
                kegg_res <- res[res$type == "KEGG", ]
            } else {
                # Heuristic: IDs start with "hsa" or "map" or numeric
                kegg_res <- res[grep("^hsa|^map|^[0-9]+$", res$ID), ]
            }

            if (nrow(kegg_res) > 0) {
                kegg_pathways[[omic]] <- kegg_res$ID
            }
        }
    }

    if (length(kegg_pathways) < 2) {
        message("Less than 2 omics have KEGG results. Skipping Pathview.")
        return(NULL)
    }

    # Find common pathways
    all_kegg <- unlist(kegg_pathways)
    if (length(all_kegg) == 0) {
        message("No KEGG pathways found.")
        return(NULL)
    }

    pathway_counts <- table(all_kegg)
    common_pathways <- names(pathway_counts)[pathway_counts >= 2]

    if (length(common_pathways) == 0) {
        message("No agreed KEGG pathways found between omics.")
        return(NULL)
    }

    message("Found ", length(common_pathways), " agreed KEGG pathways.")

    # Prepare Data for Pathview
    # Helper to get fold changes
    get_logfc <- function(omic_name, id_col = "entrez_id") {
        if (!omic_name %in% names(mae_data$harmonized_omics)) {
            return(NULL)
        }

        dat <- mae_data$harmonized_omics[[omic_name]]
        de <- dat$de_table %||% dat$da_table
        anno <- dat$feature_annotation

        if (is.null(de) || is.null(anno)) {
            return(NULL)
        }

        merged <- merge(de, anno, by = "feature_id")

        if (!id_col %in% colnames(merged)) {
            return(NULL)
        }

        # Get LogFC column
        fc_col <- grep("logFC|log2FoldChange", colnames(merged), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(fc_col)) {
            return(NULL)
        }

        # Create named vector (take mean for duplicates)
        vec <- tapply(merged[[fc_col]], merged[[id_col]], mean, na.rm = TRUE)
        return(vec)
    }

    # Gene data (Transcriptomics + Proteomics)
    gene_data <- NULL
    rna_fc <- get_logfc("transcriptomics", "entrez_id")
    prot_fc <- get_logfc("proteomics", "entrez_id")

    if (!is.null(rna_fc) && !is.null(prot_fc)) {
        all_genes <- unique(c(names(rna_fc), names(prot_fc)))
        gene_data <- matrix(NA,
            nrow = length(all_genes), ncol = 2,
            dimnames = list(all_genes, c("Transcriptomics", "Proteomics"))
        )
        idx_rna <- match(names(rna_fc), all_genes)
        gene_data[idx_rna, 1] <- rna_fc
        idx_prot <- match(names(prot_fc), all_genes)
        gene_data[idx_prot, 2] <- prot_fc
    } else if (!is.null(rna_fc)) {
        gene_data <- rna_fc
    } else if (!is.null(prot_fc)) {
        gene_data <- prot_fc
    }

    # Metabolomics Data (CPD Data)
    cpd_data <- NULL
    met_fc <- get_logfc("metabolomics", "kegg_id")
    if (!is.null(met_fc)) {
        cpd_data <- met_fc
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        message("No valid Entrez/KEGG IDs found in data for Pathview.")
        return(NULL)
    }

    # Run Pathview
    generated_plots <- list()

    cwd <- getwd()
    setwd(out_dir)
    on.exit(setwd(cwd))

    # Determine organism code
    organism <- config$global$organism %||% "human"
    if (organism == "c_elegans") {
        species_code <- "cel"
    } else {
        species_code <- mg_config$organism_code %||% "hsa"
    }

    for (pid in common_pathways) {
        # Clean ID
        clean_pid <- sub("^[a-z]+", "", pid)

        tryCatch({
            pv_out <- pathview::pathview(
                gene.data = gene_data,
                cpd.data = cpd_data,
                pathway.id = clean_pid,
                species = species_code,
                out.suffix = "multiomics",
                temp.file = TRUE,
                kegg.dir = out_dir,
                keys.align = "y",
                match.data = TRUE,
                multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                same.layer = FALSE
            )

            outfile <- paste0(species_code, clean_pid, ".multiomics.png")
            if (file.exists(outfile)) {
                generated_plots[[pid]] <- file.path(out_dir, outfile)
                message("Generated Pathview: ", outfile)
            }
        }, error = function(e) {
            message("Pathview failed for ", pid, ": ", e$message)
        })
    }

    return(generated_plots)
}


# =============================================================================
# Multi-ORA: Combined Over-Representation Analysis across omics
# =============================================================================

#' Run multi-omics ORA (multi-ORA)
#'
#' Pools significant features from all omics layers, maps them to KEGG
#' pathways, and runs a single combined hypergeometric test per pathway.
#' Gene-based omics (RNA, proteomics) are pooled via ENTREZID; metabolomics
#' uses compound-level enrichment separately. Results are combined using
#' Fisher's method.
#'
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result with MAE and pre-processing data
#' @param config Full config object
#' @param out_dir Output directory for results and plots
#' @param per_omics_enrichment This run's per-omics enrichment frames
#'   (`multiomics_cross_enrichment$per_omics`), used by the no-OrgDb pathview
#'   fallback to select pathways from the current run rather than from whatever
#'   enrichment CSVs an earlier run left on disk.
#' @return List with: results (data.frame), plots (list of paths)
run_multi_ora <- function(de_results, harmonization_res, config, out_dir,
                          per_omics_enrichment = NULL) {

    message("=== Running Multi-ORA (combined cross-omics ORA) ===")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Both pathview renderers below write into this directory and the report
    # finds their maps by globbing it, so a map this run no longer produces
    # would otherwise linger on the page as a current result. Cleared here,
    # once, rather than by either renderer: neither owns the other's output.
    #
    # Before the input check on purpose. A rerun with one omics layer fewer
    # produces no Multi-ORA at all, and that is exactly the rerun whose stale
    # maps would otherwise stay on the page looking current.
    clear_multi_ora_pathview_outputs(out_dir)

    # One gate for every pathview renderer below, resolved once. The cleanup
    # above deliberately precedes it: turning the maps off has to remove the
    # previous run's, or the report keeps showing maps nobody asked for.
    run_pathview <- isTRUE(
        (config$modes$multiomics$enrichment$pathview$run_pathview %||% TRUE))
    if (!run_pathview) {
        message("Multi-ORA: pathway maps disabled by enrichment.pathview.run_pathview")
    }

    if (is.null(de_results) || length(de_results) < 2) {
        message("Multi-ORA requires DE results from at least 2 omics layers")
        return(NULL)
    }

    # Track errors/warnings for reporting in HTML
    .ora_issues <- character(0)

    organism <- config$global$organism
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)

    if (is.null(kegg_org) || is.null(org_db)) {
        # No KEGG organism / OrgDb (e.g. non-model organisms): fall back to
        # GMT-based ORA when the omics blocks supply custom gene sets.
        gmt_res <- run_multi_ora_gmt(de_results, harmonization_res, config, out_dir)
        if (is.null(gmt_res)) {
            message("Multi-ORA: organism annotation not available for ", organism,
                    " and no usable per-omic GMT gene sets")
        }
        # Without an OrgDb, pathway maps are still reachable: KEGG's reference
        # maps are organism-independent, so a configured feature-to-KO map puts
        # this run's features onto them in KO space.
        if (run_pathview) {
            tryCatch(
                generate_per_omic_union_pathview(de_results, harmonization_res,
                                                 config, out_dir,
                                                 per_omics_enrichment = per_omics_enrichment),
                error = function(e) message("  Union pathview failed: ", conditionMessage(e))
            )
        }
        return(gmt_res)
    }

    # --- Collect per-omics significant KEGG gene IDs ---
    gene_omics <- c("transcriptomics", "proteomics")
    per_omics_sig <- list()
    per_omics_universe <- list()

    for (om in intersect(gene_omics, names(de_results))) {
        de_data <- de_results[[om]]
        de_tables <- extract_de_tables(de_data, om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next

        # Map to ENTREZ IDs
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
            error = function(e) {
                msg <- paste0(om, " ENTREZ ID mapping failed: ", e$message)
                message("  Multi-ORA: ", msg)
                .ora_issues <<- c(.ora_issues, msg)
                NULL
            }
        )
        if (is.null(id_map) || nrow(id_map) == 0) {
            msg <- paste0(om, " produced no ENTREZ mappings")
            message("  Multi-ORA: ", msg, ", skipping")
            .ora_issues <<- c(.ora_issues, msg)
            next
        }

        # Collect significant features (padj < 0.05 in any contrast)
        sig_ids <- character(0)
        all_ids <- character(0)
        for (nm in names(de_tables)) {
            df <- de_tables[[nm]]
            df_mapped <- merge(df, id_map, by = "feature_id")
            all_ids <- c(all_ids, df_mapped$ENTREZID)

            sig <- df_mapped$ENTREZID[!is.na(df_mapped$padj) & df_mapped$padj < 0.05]
            if (length(sig) < 5) {
                sig <- df_mapped$ENTREZID[!is.na(df_mapped$pvalue) & df_mapped$pvalue < 0.05]
            }
            sig_ids <- c(sig_ids, sig)
        }

        sig_ids <- unique(sig_ids[!is.na(sig_ids)])
        all_ids <- unique(all_ids[!is.na(all_ids)])

        if (length(sig_ids) > 0) {
            per_omics_sig[[om]] <- sig_ids
            per_omics_universe[[om]] <- all_ids
            message("  ", om, ": ", length(sig_ids), " sig / ", length(all_ids), " total ENTREZ IDs")
        }
    }

    if (length(per_omics_sig) == 0) {
        message("Multi-ORA: no gene-based omics produced significant features")
        .ora_issues <- c(.ora_issues, "No gene-based omics produced significant ENTREZ features")
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    # Pool significant genes across gene-based omics
    pooled_sig <- unique(unlist(per_omics_sig))
    pooled_universe <- unique(unlist(per_omics_universe))
    message("  Pooled: ", length(pooled_sig), " sig / ", length(pooled_universe),
            " universe genes across ", length(per_omics_sig), " omics")

    # Convert ENTREZID to KEGG gene IDs
    kegg_conv <- tryCatch(
        convert_entrez_to_kegg(pooled_universe, kegg_org),
        error = function(e) {
            msg <- paste0("KEGG API error during ENTREZ-to-KEGG conversion: ", e$message)
            message("  Multi-ORA: ", msg)
            .ora_issues <<- c(.ora_issues, msg)
            NULL
        }
    )
    if (is.null(kegg_conv)) {
        message("Multi-ORA: KEGG ID conversion failed (possible KEGG API timeout)")
        .ora_issues <- c(.ora_issues, "KEGG ID conversion failed — KEGG REST API may be unavailable")
        # Save issues and return
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    pooled_sig_kegg <- unique(kegg_conv[pooled_sig])
    pooled_sig_kegg <- pooled_sig_kegg[!is.na(pooled_sig_kegg)]
    pooled_univ_kegg <- unique(kegg_conv[pooled_universe])
    pooled_univ_kegg <- pooled_univ_kegg[!is.na(pooled_univ_kegg)]

    # Also convert per-omics sig lists to KEGG IDs (for per-omics ORA)
    per_omics_sig_kegg <- lapply(per_omics_sig, function(ids) {
        k <- kegg_conv[ids]
        unique(k[!is.na(k)])
    })

    # --- Run pooled ORA ---
    message("  Running pooled gene ORA...")
    pooled_ora <- run_multi_ora_kegg(
        sig_genes = pooled_sig_kegg,
        universe = pooled_univ_kegg,
        kegg_org = kegg_org,
        label = "pooled",
        exclude_classes = .excluded_pathway_classes(config)
    )

    # --- Run per-omics ORA (with same universe) ---
    per_omics_ora <- list()
    for (om in names(per_omics_sig_kegg)) {
        message("  Running per-omics ORA for ", om, "...")
        per_omics_ora[[om]] <- run_multi_ora_kegg(
            sig_genes = per_omics_sig_kegg[[om]],
            universe = pooled_univ_kegg,
            kegg_org = kegg_org,
            label = om,
            exclude_classes = .excluded_pathway_classes(config)
        )
    }

    # --- Metabolomics compound ORA (separate) ---
    metab_ora <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_ora <- tryCatch({
            de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics", harmonization_res)
            id_map <- map_metabolite_ids_to_kegg(de_tables, harmonization_res)
            if (!is.null(id_map) && nrow(id_map) > 0) {
                full_universe <- unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])
                # Merge DE stats with KEGG IDs
                de_df <- do.call(rbind, de_tables)
                de_mapped <- merge(de_df, id_map, by = "feature_id")
                de_mapped$KEGG_ID <- de_mapped$KEGG_CPD
                run_compound_ora(de_mapped, out_dir, 2, 500, 0.1,
                                 universe = full_universe,
                                 exclude_classes = .excluded_pathway_classes(config))
            } else NULL
        }, error = function(e) {
            message("  Metabolomics compound ORA failed: ", e$message)
            NULL
        })
        if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
            message("  metabolomics: ", nrow(metab_ora), " enriched compound pathways")
            # Save compound ORA results for dedicated report section
            write.csv(metab_ora, file.path(out_dir, "compound_ora_results.csv"),
                      row.names = FALSE)
            # Generate compound ORA barplot
            compound_bp_path <- file.path(out_dir, "compound_ora_barplot.png")
            png(compound_bp_path, width = 1000, height = 700, res = 120)
            tryCatch({
                plot_multi_ora_barplot(metab_ora, "Metabolomics Compound ORA (KEGG Pathways)")
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
            })
            dev.off()
            message("  Saved compound ORA barplot and results table")
        }
    }

    # --- Combine results into summary table ---
    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, metab_ora)

    if (is.null(combined) || nrow(combined) == 0) {
        message("Multi-ORA: no enriched pathways found")
        .ora_issues <- c(.ora_issues, "No enriched pathways found after ORA")
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
        return(NULL)
    }

    # Save any accumulated issues (even on success, partial issues are useful)
    if (length(.ora_issues) > 0) {
        writeLines(.ora_issues, file.path(out_dir, "multi_ora_issues.txt"))
    }

    # Write results
    write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)
    message("  Multi-ORA found ", nrow(combined), " enriched pathways (pooled)")

    # --- Generate plots ---
    plots <- list()

    # 1. Pooled ORA barplot
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        plots$pooled_barplot <- file.path(out_dir, "multi_ora_pooled_barplot.png")
        png(plots$pooled_barplot, width = 1000, height = 700, res = 120)
        tryCatch({
            plot_multi_ora_barplot(pooled_ora, "Pooled Multi-ORA (All Gene-Based Omics)")
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()
    }

    # 2. Multi-ORA dot plot comparing per-omics + pooled
    plots$dotplot <- file.path(out_dir, "multi_ora_dotplot.png")
    tryCatch({
        plot_multi_ora_dotplot(combined, per_omics_ora, metab_ora, out_dir)
    }, error = function(e) {
        message("  Multi-ORA dot plot failed: ", e$message)
    })

    # 3. Upset-style bar showing number of omics supporting each pathway
    plots$support_barplot <- file.path(out_dir, "multi_ora_support_barplot.png")
    tryCatch({
        plot_multi_ora_support(combined, out_dir)
    }, error = function(e) {
        message("  Multi-ORA support plot failed: ", e$message)
    })

    # 4. Pathview maps for pathways supported by >= 2 omics
    plots$pathview_pdf <- if (!run_pathview) NULL else tryCatch({
        generate_multi_ora_pathview(
            combined = combined,
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir,
            min_support = 2
        )
    }, error = function(e) {
        message("  Multi-ORA pathview failed: ", e$message)
        NULL
    })

    # 5. Per-omics pathview: top metabolomics pathways + proteomics overlay,
    #    and top proteomics pathways + metabolomics overlay
    per_omics_pv <- if (!run_pathview) NULL else tryCatch({
        generate_per_omics_pathview(
            per_omics_ora = per_omics_ora,
            metab_ora = metab_ora,
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir,
            top_n = 5
        )
    }, error = function(e) {
        message("  Per-omics pathview failed: ", e$message)
        NULL
    })
    if (!is.null(per_omics_pv)) {
        plots$pathview_metabolomics_pdf <- per_omics_pv$metabolomics_pdf
        plots$pathview_proteomics_pdf <- per_omics_pv$proteomics_pdf
    }

    # --- Per-contrast Multi-ORA ---
    # Re-extract DE tables to get per-contrast names, then run ORA per contrast
    all_de_tables <- list()
    for (om in intersect(gene_omics, names(de_results))) {
        de_tables_tmp <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (!is.null(de_tables_tmp)) all_de_tables[[om]] <- de_tables_tmp
    }
    # Also extract metabolomics DE tables for per-contrast compound ORA
    metab_de_tables <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                             harmonization_res)
    }
    ora_contrast_names <- unique(unlist(lapply(all_de_tables, names)))

    if (length(ora_contrast_names) > 1) {
        message("  Generating per-contrast Multi-ORA for ",
                length(ora_contrast_names), " contrasts")
        per_contrast_dir <- file.path(out_dir, "per_contrast")

        for (cname in ora_contrast_names) {
            safe_dir <- gsub("[^a-zA-Z0-9._-]", "_", cname)
            contrast_out <- file.path(per_contrast_dir, safe_dir)
            dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

            tryCatch({
                .run_multi_ora_contrast_group(
                    all_de_tables = all_de_tables,
                    contrast_name = cname,
                    harmonization_res = harmonization_res,
                    kegg_org = kegg_org,
                    org_db = org_db,
                    out_dir = contrast_out,
                    metab_de_tables = metab_de_tables,
                    exclude_classes = .excluded_pathway_classes(config)
                )
            }, error = function(e) {
                message("    Per-contrast Multi-ORA failed for ", cname, ": ", e$message)
            })
        }
        message("  Per-contrast Multi-ORA output saved to: ", per_contrast_dir)
    }

    message("Multi-ORA complete: ", nrow(combined), " pathways")

    list(
        pooled = pooled_ora,
        per_omics = per_omics_ora,
        metabolomics = metab_ora,
        combined = combined,
        plots = plots
    )
}


#' Run KEGG ORA with Fisher's exact test for multi-ORA
#'
#' Run Multi-ORA for a single contrast
#'
#' @param all_de_tables Named list (per omics) of named lists (per contrast) of DE data frames
#' @param contrast_name Contrast name to process
#' @param harmonization_res Harmonization result
#' @param kegg_org KEGG organism code
#' @param org_db Organism annotation database
#' @param out_dir Output directory for this contrast
#' @param metab_de_tables Metabolomics DE tables (named list per contrast), or NULL
#' @param exclude_classes BRITE classes this project excludes from its report,
#'   passed down so a per-contrast section cannot show a class the run-level
#'   sections removed.
#' @return Invisible NULL
.run_multi_ora_contrast_group <- function(all_de_tables, contrast_name,
                                           harmonization_res, kegg_org, org_db,
                                           out_dir, metab_de_tables = NULL,
                                           exclude_classes = NULL) {

    per_omics_sig <- list()
    per_omics_universe <- list()

    for (om in names(all_de_tables)) {
        de_tables <- all_de_tables[[om]]
        if (!contrast_name %in% names(de_tables)) next

        df <- de_tables[[contrast_name]]
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables[contrast_name], om, harmonization_res, org_db),
            error = function(e) NULL
        )
        if (is.null(id_map) || nrow(id_map) == 0) next

        df_mapped <- merge(df, id_map, by = "feature_id")
        all_ids <- unique(df_mapped$ENTREZID[!is.na(df_mapped$ENTREZID)])
        sig <- df_mapped$ENTREZID[!is.na(df_mapped$padj) & df_mapped$padj < 0.05]
        if (length(sig) < 5) {
            sig <- df_mapped$ENTREZID[!is.na(df_mapped$pvalue) & df_mapped$pvalue < 0.05]
        }
        sig <- unique(sig[!is.na(sig)])

        if (length(sig) > 0) {
            per_omics_sig[[om]] <- sig
            per_omics_universe[[om]] <- all_ids
        }
    }

    if (length(per_omics_sig) == 0) return(invisible(NULL))

    pooled_sig <- unique(unlist(per_omics_sig))
    pooled_universe <- unique(unlist(per_omics_universe))

    kegg_conv <- convert_entrez_to_kegg(pooled_universe, kegg_org)
    if (is.null(kegg_conv)) return(invisible(NULL))

    pooled_sig_kegg <- unique(kegg_conv[pooled_sig])
    pooled_sig_kegg <- pooled_sig_kegg[!is.na(pooled_sig_kegg)]
    pooled_univ_kegg <- unique(kegg_conv[pooled_universe])
    pooled_univ_kegg <- pooled_univ_kegg[!is.na(pooled_univ_kegg)]

    pooled_ora <- run_multi_ora_kegg(pooled_sig_kegg, pooled_univ_kegg, kegg_org,
                                     "pooled", exclude_classes = exclude_classes)

    per_omics_ora <- list()
    for (om in names(per_omics_sig)) {
        k <- kegg_conv[per_omics_sig[[om]]]
        k <- unique(k[!is.na(k)])
        per_omics_ora[[om]] <- run_multi_ora_kegg(k, pooled_univ_kegg, kegg_org, om,
                                                  exclude_classes = exclude_classes)
    }

    # Run per-contrast metabolomics compound ORA if data is available
    contrast_metab_ora <- NULL
    if (!is.null(metab_de_tables)) {
        # Find matching contrast name (try exact, then fuzzy match)
        metab_cname <- if (contrast_name %in% names(metab_de_tables)) {
            contrast_name
        } else {
            # Fuzzy: try normalized contrast names
            norm_cn <- tolower(gsub("[^a-z0-9]", "", contrast_name))
            metab_norms <- tolower(gsub("[^a-z0-9]", "", names(metab_de_tables)))
            idx <- match(norm_cn, metab_norms)
            if (!is.na(idx)) names(metab_de_tables)[idx] else NULL
        }
        if (!is.null(metab_cname)) {
            contrast_metab_ora <- tryCatch({
                de_list <- metab_de_tables[metab_cname]
                id_map <- map_metabolite_ids_to_kegg(de_list, harmonization_res)
                if (!is.null(id_map) && nrow(id_map) > 0) {
                    full_universe <- unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])
                    de_df <- de_list[[1]]
                    de_mapped <- merge(de_df, id_map, by = "feature_id")
                    de_mapped$KEGG_ID <- de_mapped$KEGG_CPD
                    run_compound_ora(de_mapped, out_dir, 2, 500, 0.1,
                                     universe = full_universe,
                                     exclude_classes = exclude_classes)
                } else NULL
            }, error = function(e) {
                message("    Per-contrast compound ORA failed for ", contrast_name,
                        ": ", e$message)
                NULL
            })
        }
    }

    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, contrast_metab_ora)
    if (!is.null(combined) && nrow(combined) > 0) {
        write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)

        if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
            bp_path <- file.path(out_dir, "multi_ora_pooled_barplot.png")
            png(bp_path, width = 1000, height = 700, res = 120)
            tryCatch({
                plot_multi_ora_barplot(pooled_ora,
                    paste0("Multi-ORA (", gsub("_", " ", contrast_name), ")"))
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
            })
            dev.off()
        }
    }

    invisible(NULL)
}


#' @param sig_genes Significant KEGG gene IDs
#' @param universe All KEGG gene IDs (shared universe)
#' @param kegg_org KEGG organism code
#' @param label Label for messages
#' @param pval_cutoff Adjusted p-value cutoff for the preferred branch.
#' @param exclude_classes BRITE classes this project leaves out of its report.
#'   Applied to the finished table on the way out, so everything downstream --
#'   the summary, the plots, the OrgDb pathview renderers -- inherits an already
#'   filtered input instead of filtering again.
#' @return data.frame with ORA results
run_multi_ora_kegg <- function(sig_genes, universe, kegg_org,
                                label = "pooled", pval_cutoff = 0.1,
                                exclude_classes = NULL) {

    if (length(sig_genes) < 3) {
        message("    ", label, ": too few significant genes (", length(sig_genes), ")")
        return(NULL)
    }

    # Try clusterProfiler first — use lenient cutoff, filter manually after
    ora_res <- tryCatch({
        res <- clusterProfiler::enrichKEGG(
            gene = sig_genes,
            universe = universe,
            organism = kegg_org,
            keyType = "kegg",
            minGSSize = 5,
            maxGSSize = 500,
            # Open both cutoffs and let this wrapper threshold below; the default
            # qvalueCutoff = 0.2 would otherwise drop rows with raw p < 0.05 but
            # q >= 0.2 before the raw-p fallback sees them (same as the enricher path).
            pvalueCutoff = 1.0,
            qvalueCutoff = 1.0
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway = df$Description,
                ID = df$ID,
                pvalue = df$pvalue,
                padj = df$p.adjust,
                GeneRatio = df$GeneRatio,
                Count = df$Count,
                geneID = df$geneID,
                stringsAsFactors = FALSE
            )
            # Filter: prefer padj, fall back to pvalue < 0.05.
            # Which branch is taken is decided on the unfiltered results, so a
            # project's class exclusion cannot move the run from adjusted hits
            # to the raw-p fallback. Exclusion applies to whichever table this
            # chose, on its way out.
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) {
                return(.exclude_kegg_classes(padj_hits, exclude_classes,
                                             kegg_org, label))
            }
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    ", label, ": padj too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                return(.exclude_kegg_classes(pval_hits, exclude_classes,
                                             kegg_org, label))
            }
        }
        NULL
    }, error = function(e) {
        message("    ", label, " clusterProfiler ORA failed: ", e$message)
        NULL
    })

    if (!is.null(ora_res)) {
        message("    ", label, ": ", nrow(ora_res), " enriched pathways")
        return(ora_res)
    }

    # Fallback: Fisher's exact test
    .exclude_kegg_classes(
        run_ora_kegg_fisher(sig_genes, universe, kegg_org, 5, 500, pval_cutoff),
        exclude_classes, kegg_org, label)
}


#' Drop excluded KEGG classes from a finished gene-ORA table
#'
#' The one place the gene-based multi-ORA applies the exclusion, so every return
#' path of \code{run_multi_ora_kegg()} filters identically and the tables that
#' feed the summary, the plots and the OrgDb pathview renderers arrive already
#' filtered -- rather than each of those growing a filter of its own.
#'
#' Applied to completed results: the tested universe, the p-values and the
#' adjustment behind them are exactly what they were.
#'
#' @param df Finished ORA table, or NULL.
#' @param exclude_classes BRITE classes to drop; NULL or empty is a no-op.
#' @param kegg_org Active KEGG organism code, for accession recognition.
#' @param label Short context word for the message naming what was dropped.
#' @param classification Resolved class table, defaulted lazily to the fetch so
#'   a call with nothing to exclude never reaches the network, and so a test can
#'   supply one without depending on whether a machine has any.
#' @return \code{df} with the excluded rows removed, or NULL when nothing is
#'   left -- the same "no pathways" shape every other path here returns.
#' @keywords internal
.exclude_kegg_classes <- function(df, exclude_classes = NULL, kegg_org = NULL,
                                  label = "pathways",
                                  classification = kegg_pathway_categories()) {
    if (is.null(df) || nrow(df) == 0) return(df)
    if (length(unlist(exclude_classes)) == 0) return(df)

    ids <- if ("ID" %in% names(df)) df$ID else df$pathway
    df <- df[keep_kegg_pathways(ids, exclude = exclude_classes,
                                kegg_org = kegg_org, label = label,
                                classification = classification), ,
             drop = FALSE]
    if (nrow(df) == 0) NULL else df
}


#' Convert a GMT file to clusterProfiler TERM2GENE / TERM2NAME frames
#'
#' Reuses \code{read_gmt()} (a named list of gene vectors carrying a
#' \code{descriptions} attribute) and reshapes it for
#' \code{clusterProfiler::enricher}.
#'
#' @param gmt_file Path to a GMT file.
#' @return A list with \code{t2g} (data.frame term, gene) and \code{t2n}
#'   (data.frame term, name), or NULL if the file has no usable gene sets.
gmt_to_term2gene <- function(gmt_file) {
    gs <- tryCatch(read_gmt(gmt_file), error = function(e) NULL)
    if (is.null(gs) || length(gs) == 0) return(NULL)
    descr <- attr(gs, "descriptions")
    t2g <- data.frame(
        term = rep(names(gs), lengths(gs)),
        gene = unlist(gs, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    t2n <- if (!is.null(descr)) {
        data.frame(term = names(descr), name = unname(descr), stringsAsFactors = FALSE)
    } else NULL
    list(t2g = t2g, t2n = t2n)
}


#' Run ORA with clusterProfiler::enricher against a custom TERM2GENE
#'
#' GMT-based counterpart to \code{run_multi_ora_kegg()}; returns the same
#' data.frame shape so the shared multi-ORA summary and plots consume it
#' unchanged.
#'
#' @param sig_genes Significant feature IDs (native namespace, e.g. EHI_ / XP_).
#' @param universe All tested feature IDs (shared background).
#' @param term2gene data.frame(term, gene) gene-set membership.
#' @param term2name Optional data.frame(term, name) for readable pathway labels.
#' @param label Label used in progress messages.
#' @param pval_cutoff Adjusted-p cutoff (falls back to raw p < 0.05).
#' @return data.frame(pathway, ID, pvalue, padj, GeneRatio, Count, geneID), or NULL.
run_multi_ora_enricher <- function(sig_genes, universe, term2gene, term2name = NULL,
                                   label = "pooled", pval_cutoff = 0.1) {

    if (length(sig_genes) < 3) {
        message("    ", label, ": too few significant genes (", length(sig_genes), ")")
        return(NULL)
    }

    ora_res <- tryCatch({
        res <- clusterProfiler::enricher(
            gene = sig_genes,
            universe = universe,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            minGSSize = 5,
            maxGSSize = 500,
            # Ask enricher for every tested set (p and q cutoffs both open); this
            # wrapper does its own padj / raw-p thresholding below. Leaving the
            # default qvalueCutoff = 0.2 would silently drop rows with raw
            # p < 0.05 but q >= 0.2 before the raw-p fallback ever sees them.
            pvalueCutoff = 1.0,
            qvalueCutoff = 1.0
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway   = if (!is.null(df$Description)) df$Description else df$ID,
                ID        = df$ID,
                pvalue    = df$pvalue,
                padj      = df$p.adjust,
                GeneRatio = df$GeneRatio,
                Count     = df$Count,
                geneID    = df$geneID,
                stringsAsFactors = FALSE
            )
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) return(padj_hits)
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    ", label, ": padj too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                return(pval_hits)
            }
        }
        NULL
    }, error = function(e) {
        message("    ", label, " enricher ORA failed: ", e$message)
        NULL
    })

    if (!is.null(ora_res)) {
        message("    ", label, ": ", nrow(ora_res), " enriched pathways")
    }
    ora_res
}


#' GMT-based multi-omics ORA (fallback when no KEGG/OrgDb organism)
#'
#' Runs over-representation analysis on the pooled and per-omic significant
#' gene-based DE features using each omic's custom GMT
#' (\code{modes.<omic>.pathway.gmt_file}), via \code{clusterProfiler::enricher}.
#' The per-omic GMTs share pathway IDs across namespaces (e.g. EHI_ for RNA and
#' XP_ for proteomics), so their TERM2GENE tables are row-bound into one
#' collection and a pooled mixed-namespace hit list enriches against it — each
#' omic's hits match that omic's members under the same term. A gene measured in
#' two omics contributes once per omic (pooled multi-omic evidence).
#'
#' Per-contrast ORA is intentionally not run here: the multi-omics runs are
#' single-contrast, so pooled + per-omic already covers the whole comparison.
#'
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result (for \code{extract_de_tables}).
#' @param config Full config object.
#' @param out_dir Output directory for results and plots.
#' @return list(pooled, per_omics, metabolomics, combined, plots), or NULL.
run_multi_ora_gmt <- function(de_results, harmonization_res, config, out_dir) {

    gene_omics   <- c("transcriptomics", "proteomics")
    omic_cfg_key <- c(transcriptomics = "rna", proteomics = "proteomics")

    per_omics_sig  <- list()
    per_omics_univ <- list()
    per_omics_t2g  <- list()
    per_omics_t2n  <- list()

    for (om in intersect(gene_omics, names(de_results))) {
        # gmt_file may be a single path or a YAML list of paths (GO + KEGG);
        # read_gmt() already merges several files, so only these guards needed
        # to vectorise — with a list they were comparing length-2 vectors and
        # aborting the whole Multi-ORA step.
        gmt_path <- unlist(config$modes[[omic_cfg_key[[om]]]]$pathway$gmt_file,
                           use.names = FALSE)
        if (length(gmt_path) == 0 || !any(nzchar(gmt_path))) next
        # Resolve like every other user-supplied input (metabolomics enrichment,
        # data files): absolute paths pass through, relative ones resolve under
        # the raw/ data dir. resolve_raw_path() would mangle an absolute path.
        gmt_abs <- resolve_input_path(config, gmt_path)
        if (any(!file.exists(gmt_abs))) {
            message("  Multi-ORA (GMT): ", om, " gmt_file not found: ",
                    paste(gmt_abs[!file.exists(gmt_abs)], collapse = ", "))
            next
        }
        gs <- gmt_to_term2gene(gmt_abs)
        if (is.null(gs) || nrow(gs$t2g) == 0) next

        de_tables <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next

        sig <- character(0); univ <- character(0)
        for (nm in names(de_tables)) {
            df  <- de_tables[[nm]]
            # Match GMT members by bare ID: RNA feature_ids may carry a "Gene:"
            # prefix (e.g. "Gene:EHI_012345") while custom GMTs list the bare ID.
            # The single-omics RNA pathway module strips it the same way.
            fid <- sub("^Gene:", "", as.character(df$feature_id))
            univ <- c(univ, fid)
            s <- fid[!is.na(df$padj) & df$padj < 0.05]
            if (length(s) < 5) s <- fid[!is.na(df$pvalue) & df$pvalue < 0.05]
            sig <- c(sig, s)
        }
        sig  <- unique(sig[!is.na(sig)])
        univ <- unique(univ[!is.na(univ)])
        if (length(sig) == 0) next

        per_omics_sig[[om]]  <- sig
        per_omics_univ[[om]] <- univ
        per_omics_t2g[[om]]  <- gs$t2g
        per_omics_t2n[[om]]  <- gs$t2n
        message("  Multi-ORA (GMT) ", om, ": ", length(sig), " sig / ",
                length(univ), " tested features")
    }

    # A cross-omics ORA needs at least two gene-based GMT omics. With only one
    # (e.g. a non-model RNA+metabolomics run, where metabolomics is not covered
    # by this gene-based fallback) the pooled result would be single-omic while
    # presenting as multi-omics, so decline instead.
    if (length(per_omics_sig) < 2) {
        message("Multi-ORA (GMT): need >= 2 gene-based omics with a usable GMT + ",
                "sig features; found ", length(per_omics_sig), " -- skipping")
        return(NULL)
    }

    comb_t2g    <- unique(do.call(rbind, per_omics_t2g))
    comb_t2n    <- unique(do.call(rbind, Filter(Negate(is.null), per_omics_t2n)))
    pooled_sig  <- unique(unlist(per_omics_sig))
    pooled_univ <- unique(unlist(per_omics_univ))

    message("  Running pooled GMT ORA...")
    pooled_ora <- run_multi_ora_enricher(pooled_sig, pooled_univ, comb_t2g, comb_t2n, "pooled")

    per_omics_ora <- list()
    for (om in names(per_omics_sig)) {
        message("  Running per-omics GMT ORA for ", om, "...")
        per_omics_ora[[om]] <- run_multi_ora_enricher(
            per_omics_sig[[om]], per_omics_univ[[om]],
            per_omics_t2g[[om]], per_omics_t2n[[om]], om)
    }

    combined <- build_multi_ora_summary(pooled_ora, per_omics_ora, NULL)
    if (is.null(combined) || nrow(combined) == 0) {
        message("Multi-ORA (GMT): no enriched pathways found")
        return(NULL)
    }
    write.csv(combined, file.path(out_dir, "multi_ora_results.csv"), row.names = FALSE)
    message("  Multi-ORA (GMT) found ", nrow(combined), " enriched pathways (pooled)")

    plots <- list()
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        plots$pooled_barplot <- file.path(out_dir, "multi_ora_pooled_barplot.png")
        png(plots$pooled_barplot, width = 1000, height = 700, res = 120)
        tryCatch(
            # The one figure in this file whose table really mixes namespaces:
            # comb_t2g above row-binds each omic's GMTs, so GO, KEGG and Pfam
            # sets land in one pool and the largest collection takes every slot
            # on set count alone. resolve_kegg_org_code() is NULL for an
            # organism KEGG does not know -- the usual case on this fallback --
            # and then only the map/ko/bare spellings count as KEGG, which is
            # the identity contract answering honestly rather than a gap.
            plot_multi_ora_barplot(pooled_ora, "Pooled Multi-ORA (GMT gene sets)",
                                   by_collection = TRUE,
                                   kegg_org = resolve_kegg_org_code(
                                       config$global$organism)),
            error = function(e) { plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2) }
        )
        dev.off()
    }
    plots$dotplot <- file.path(out_dir, "multi_ora_dotplot.png")
    tryCatch(plot_multi_ora_dotplot(combined, per_omics_ora, NULL, out_dir),
             error = function(e) message("  Multi-ORA dot plot failed: ", e$message))
    plots$support_barplot <- file.path(out_dir, "multi_ora_support_barplot.png")
    tryCatch(plot_multi_ora_support(combined, out_dir),
             error = function(e) message("  Multi-ORA support plot failed: ", e$message))

    list(pooled = pooled_ora, per_omics = per_omics_ora, metabolomics = NULL,
         combined = combined, plots = plots)
}


#' Build combined multi-ORA summary table
#'
#' Merges pooled ORA with per-omics results into a single table showing
#' which omics support each pathway. Uses full outer join so compound-only
#' pathways from metabolomics are also included.
build_multi_ora_summary <- function(pooled_ora, per_omics_ora, metab_ora) {

    if ((is.null(pooled_ora) || nrow(pooled_ora) == 0) &&
        (is.null(metab_ora) || nrow(metab_ora) == 0)) return(NULL)

    # Start from gene-based pooled ORA if available
    if (!is.null(pooled_ora) && nrow(pooled_ora) > 0) {
        summary <- pooled_ora[, c("pathway", "ID", "pvalue", "padj", "Count"), drop = FALSE]
        colnames(summary)[colnames(summary) == "pvalue"] <- "pooled_pvalue"
        colnames(summary)[colnames(summary) == "padj"] <- "pooled_padj"
        colnames(summary)[colnames(summary) == "Count"] <- "pooled_count"
    } else {
        summary <- data.frame(
            pathway = character(0), ID = character(0),
            pooled_pvalue = numeric(0), pooled_padj = numeric(0),
            pooled_count = integer(0), stringsAsFactors = FALSE
        )
    }

    # Add normalized ID column for cross-prefix matching
    summary$norm_id <- normalize_kegg_pathway_id(summary$ID)

    # Add per-omics support columns (padj and pvalue)
    omics_names <- names(per_omics_ora)
    for (om in omics_names) {
        om_res <- per_omics_ora[[om]]
        if (is.null(om_res) || nrow(om_res) == 0) {
            summary[[paste0(om, "_padj")]] <- NA_real_
            summary[[paste0(om, "_pvalue")]] <- NA_real_
            next
        }
        om_padj <- setNames(om_res$padj, om_res$ID)
        summary[[paste0(om, "_padj")]] <- om_padj[summary$ID]
        if ("pvalue" %in% colnames(om_res)) {
            om_pval <- setNames(om_res$pvalue, om_res$ID)
            summary[[paste0(om, "_pvalue")]] <- om_pval[summary$ID]
        }
    }

    # Add metabolomics via full outer join on normalized pathway IDs
    if (!is.null(metab_ora) && nrow(metab_ora) > 0 && "ID" %in% colnames(metab_ora)) {
        met_norm <- normalize_kegg_pathway_id(metab_ora$ID)
        met_padj <- setNames(metab_ora$padj, met_norm)
        met_pval <- setNames(metab_ora$pvalue, met_norm)

        # Match existing rows
        summary$metabolomics_padj <- met_padj[summary$norm_id]
        summary$metabolomics_pvalue <- met_pval[summary$norm_id]

        # Find compound-only pathways (not in gene ORA results)
        new_ids <- setdiff(met_norm, summary$norm_id)
        if (length(new_ids) > 0) {
            new_rows <- data.frame(
                pathway = metab_ora$pathway[match(new_ids, met_norm)],
                ID = metab_ora$ID[match(new_ids, met_norm)],
                pooled_pvalue = NA_real_,
                pooled_padj = NA_real_,
                pooled_count = NA_integer_,
                norm_id = new_ids,
                stringsAsFactors = FALSE
            )
            # Add per-omics columns as NA
            for (om in omics_names) {
                new_rows[[paste0(om, "_padj")]] <- NA_real_
                new_rows[[paste0(om, "_pvalue")]] <- NA_real_
            }
            new_rows$metabolomics_padj <- met_padj[new_ids]
            new_rows$metabolomics_pvalue <- met_pval[new_ids]

            summary <- rbind(summary, new_rows)
        }
    }

    # Count supporting omics (padj < 0.05)
    padj_cols <- grep("_padj$", colnames(summary), value = TRUE)
    padj_cols <- setdiff(padj_cols, "pooled_padj")
    padj_mat <- sapply(padj_cols, function(col) !is.na(summary[[col]]) & summary[[col]] < 0.05)
    if (!is.matrix(padj_mat)) padj_mat <- matrix(padj_mat, ncol = length(padj_cols))
    summary$n_omics_support <- rowSums(padj_mat)

    # Fallback: count support using pvalue < 0.05 when padj is too strict
    pval_cols <- grep("_pvalue$", colnames(summary), value = TRUE)
    pval_cols <- setdiff(pval_cols, "pooled_pvalue")
    if (length(pval_cols) > 0) {
        pval_mat <- sapply(pval_cols, function(col) !is.na(summary[[col]]) & summary[[col]] < 0.05)
        if (!is.matrix(pval_mat)) pval_mat <- matrix(pval_mat, ncol = length(pval_cols))
        summary$n_omics_support_pval <- rowSums(pval_mat)
    } else {
        summary$n_omics_support_pval <- summary$n_omics_support
    }

    # Sort: gene ORA pathways first (by pooled_pvalue), then compound-only.
    # metabolomics_padj only exists when metab_ora was supplied; order() rejects a
    # NULL key ("argument N is not a vector"), so add it only when present.
    sort_keys <- list(is.na(summary$pooled_pvalue), summary$pooled_pvalue)
    if (!is.null(summary$metabolomics_padj)) {
        sort_keys <- c(sort_keys, list(summary$metabolomics_padj))
    }
    summary <- summary[do.call(order, sort_keys), ]

    # Drop internal helper column
    summary$norm_id <- NULL
    summary
}


#' Bar colour for each gene-set collection
#'
#' Fixed by name rather than assigned in the order collections happen to appear,
#' so a collection keeps its colour between runs and between contrasts. KEGG
#' keeps the purple every one of these bar plots used before there were
#' collections, so the KEGG-only figures look as they did.
#'
#' @keywords internal
.ORA_COLLECTION_COLOURS <- c(KEGG     = "#7B2D8E",
                             GO       = "#2C7FB8",
                             Pfam     = "#41AB5D",
                             InterPro = "#D95F0E",
                             Other    = "#9E9E9E")


#' Which significance column a pathway figure should show
#'
#' One choice for the whole table, never row by row: a figure whose bars are
#' part adjusted and part raw p-values has no readable axis and no honest
#' caption. Adjusted p is preferred wherever the table actually carries usable
#' values -- \code{any(is.finite())} rather than column presence, since an
#' all-NA `padj` column is common and would otherwise silence every bar.
#'
#' @param df Enrichment table, expected to carry `padj` and/or `pvalue`.
#' @return A list with `column` (the column name to read, or NA_character_ when
#'   the table carries neither), `values` (that column as numeric, or all NA),
#'   and `adjusted` (TRUE when the choice fell on an adjusted p-value).
#' @keywords internal
.ora_display_score <- function(df) {
    as_num <- function(col) suppressWarnings(as.numeric(df[[col]]))

    if ("padj" %in% names(df)) {
        padj <- as_num("padj")
        if (any(is.finite(padj))) {
            return(list(column = "padj", values = padj, adjusted = TRUE))
        }
    }
    if ("pvalue" %in% names(df)) {
        return(list(column = "pvalue", values = as_num("pvalue"),
                    adjusted = FALSE))
    }
    list(column = NA_character_, values = rep(NA_real_, nrow(df)),
         adjusted = FALSE)
}


#' Order pathway rows strongest-first, deterministically
#'
#' Every key is derived from the data: the significance score, then the
#' normalized pathway identity, then the readable label. Arrival index is
#' deliberately not a key -- a figure whose contents depend on the order rows
#' happened to be bound in is not reproducible, and the tables reaching here are
#' assembled from several sources. Rows still equal on all three are
#' indistinguishable on every field this function can see, so their relative
#' order carries no meaning.
#'
#' @param df Enrichment table.
#' @param score Numeric score for each row, smaller is better, from
#'   \code{.ora_display_score()}.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return Integer permutation of \code{seq_len(nrow(df))}.
#' @keywords internal
.order_ora_rows <- function(df, score, kegg_org = NULL) {
    order(score,
          pathway_join_key(df, kegg_org),
          pathway_display_label(df),
          na.last = TRUE)
}


#' Classify pathway identifiers by the gene-set collection they come from
#'
#' A pooled ORA over several GMTs mixes namespaces that nothing downstream can
#' tell apart once the tables are row-bound, so the collection has to be read
#' back off the identifier.
#'
#' KEGG is decided by the identity contract and nothing else:
#' \code{pathway_join_key()} resolves identity per row and normalizes the
#' spellings this pipeline produces, and \code{is_kegg_pathway_accession()}
#' decides which of those are KEGG accessions. A shape-only rule such as "two to
#' four letters then five digits" would claim a custom gene set named
#' `abcd12345`, which is exactly the silent misclassification the contract
#' exists to prevent.
#'
#' The other three patterns are fully anchored for the same reason: an
#' unanchored `^GO:?[0-9]+` claims `GO12345_signalling`, a custom set name that
#' has nothing to do with the Gene Ontology.
#'
#' Nothing is dropped and nothing is rewritten. A row whose identity is missing
#' entirely is `Other`, not an error and not a gap.
#'
#' @param df Enrichment table.
#' @param kegg_org Active KEGG organism code for the run, or NULL when the
#'   organism has no KEGG code.
#' @param keys Pre-computed join keys, defaulting to deriving them. The default
#'   is lazy, so a caller that already has them -- as
#'   \code{select_top_ora_per_collection()} does -- does not pay for them twice.
#' @return Character vector, one collection name per row of \code{df}: one of
#'   "KEGG", "GO", "Pfam", "InterPro" or "Other".
#' @examples
#' df <- data.frame(ID = c("hsa04110", "GO:0006915", "PF00069", "myset"))
#' classify_pathway_collection(df, kegg_org = "hsa")
#' # "KEGG" "GO" "Pfam" "Other"
classify_pathway_collection <- function(df, kegg_org = NULL,
                                        keys = pathway_join_key(df, kegg_org)) {
    out <- rep("Other", length(keys))

    # KEGG first and by contract. The remaining patterns are only ever offered
    # keys KEGG has already declined, so no identifier can match two rules.
    is_kegg <- is_kegg_pathway_accession(keys, kegg_org)
    out[is_kegg] <- "KEGG"

    rest <- !is_kegg & !is.na(keys)
    out[rest & grepl("^GO:[0-9]{7}$", keys)]  <- "GO"
    out[rest & grepl("^PF[0-9]{5}$", keys)]   <- "Pfam"
    out[rest & grepl("^IPR[0-9]{6}$", keys)]  <- "InterPro"

    out
}


#' Pick top ORA terms while keeping every collection represented
#'
#' A pooled ORA over GO + KEGG + Pfam is dominated by GO on set count alone --
#' thousands of GO terms against a few hundred KEGG maps -- so a plain top-n
#' leaves a figure that looks like a GO-only analysis however much evidence the
#' other collections hold.
#'
#' Terms are therefore drawn round-robin: collections enter the rotation
#' best-first, and each contributes its next-best term per pass until the figure
#' is full. No quota is computed, which is the point -- a collection with fewer
#' terms than the others simply stops appearing in later passes and the
#' remaining slots go to whoever still has terms, so there is nothing to
#' redistribute and no allocation to get wrong.
#'
#' Round-robin decides *membership* only. The rows come back in global evidence
#' order, so the figure still reads strongest-first and a reader is never told
#' that the second bar outranks the third when it does not.
#'
#' With a single collection this is exactly the plain top-n. With fewer slots
#' than collections the strongest collections are still the ones shown.
#'
#' Selection runs on whatever table it is handed: the KEGG class exclusion has
#' already been applied to these results by \code{run_multi_ora_kegg()} on the
#' way out, and nothing here reclassifies or re-filters a pathway.
#'
#' @param ora_df ORA table carrying `pvalue` and/or `padj`, plus the identity
#'   columns \code{pathway_join_key()} reads.
#' @param top_n Maximum number of rows to keep.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return A subset of \code{ora_df}, in global evidence order.
#' @examples
#' ora <- data.frame(ID = c("GO:0000001", "GO:0000002", "hsa04110"),
#'                   pathway = c("a", "b", "c"), padj = c(0.001, 0.002, 0.04))
#' select_top_ora_per_collection(ora, top_n = 2, kegg_org = "hsa")$ID
#' # "GO:0000001" "hsa04110"  -- KEGG is not displaced by the second GO term
select_top_ora_per_collection <- function(ora_df, top_n = 20, kegg_org = NULL) {
    if (is.null(ora_df) || nrow(ora_df) == 0) return(ora_df)

    score <- .ora_display_score(ora_df)$values
    keys  <- pathway_join_key(ora_df, kegg_org)

    ord   <- .order_ora_rows(ora_df, score, kegg_org)
    df    <- ora_df[ord, , drop = FALSE]
    score <- score[ord]
    keys  <- keys[ord]

    # Nothing to balance: every row is shown, so the collections cannot crowd
    # each other out and the round-robin below would only reorder them.
    if (nrow(df) <= top_n) return(df)

    coll <- classify_pathway_collection(df, kegg_org, keys = keys)
    by_coll <- split(seq_len(nrow(df)), coll)

    # Best-first, so fewer slots than collections still shows the strongest
    # ones. The collection name breaks a tie on the best score, because
    # split() names the groups and nothing else here distinguishes them.
    best <- vapply(by_coll, function(i) score[i[1]], numeric(1))
    by_coll <- by_coll[order(best, names(by_coll), na.last = TRUE)]

    keep <- integer(0)
    rank <- 1L
    while (length(keep) < top_n && any(lengths(by_coll) >= rank)) {
        for (idx in by_coll) {
            if (length(keep) >= top_n) break
            if (length(idx) >= rank) keep <- c(keep, idx[rank])
        }
        rank <- rank + 1L
    }

    # sort() restores global evidence order: df is already in it, so the row
    # numbers carry it.
    df[sort(keep), , drop = FALSE]
}


#' Plot multi-ORA pooled barplot
#'
#' Horizontal bar plot of the strongest enriched terms.
#'
#' Bars show the adjusted p-value wherever the table carries usable ones and the
#' raw p-value otherwise, and that one choice drives the ranking, the bar
#' length, the threshold line and the axis label together. It was previously
#' ranked and drawn on raw p while the report legend described a padj threshold,
#' so the figure and its caption disagreed about which statistic was on screen.
#'
#' Collection-aware selection is opt-in, because only the pooled GMT figure
#' actually mixes namespaces -- the KEGG bar plots hold KEGG maps alone, where
#' round-robin would be an elaborate way of writing plain top-n. Where it is on,
#' each bar is tagged and coloured by its collection so GO terms, KEGG maps and
#' Pfam domains are never pooled into one anonymous ranking.
#'
#' @param ora_df ORA table carrying `pvalue` and/or `padj`.
#' @param title Plot title.
#' @param top_n Maximum number of bars.
#' @param by_collection Draw terms round-robin across gene-set collections
#'   rather than by plain top-n, tagging and colouring each bar with its
#'   collection. FALSE by default: the caller opts in.
#' @param kegg_org Active KEGG organism code, or NULL. Read only when
#'   \code{by_collection} is TRUE, which is the only branch that classifies.
#' @return Invisibly, a data frame describing the bars drawn -- `key`, `label`,
#'   `collection` (NA where nothing was classified), `score` -- in the order
#'   they were selected, so which terms a figure leads with can be checked
#'   without reading pixels. NULL when there was nothing to draw.
plot_multi_ora_barplot <- function(ora_df, title, top_n = 20,
                                    by_collection = FALSE, kegg_org = NULL) {
    if (is.null(ora_df) || nrow(ora_df) == 0) return(invisible(NULL))

    # Resolved on the whole table, then read from the selected rows. Choosing
    # again on the subset could land on a different column than the ranking
    # used, which is how a figure ends up ordered by one statistic and drawn
    # with another.
    score_col <- .ora_display_score(ora_df)
    if (is.na(score_col$column)) return(invisible(NULL))

    if (isTRUE(by_collection)) {
        df <- select_top_ora_per_collection(ora_df, top_n, kegg_org)
        collection <- classify_pathway_collection(df, kegg_org)
    } else {
        ord <- .order_ora_rows(ora_df, score_col$values, kegg_org)
        df <- ora_df[ord[seq_len(min(top_n, nrow(ora_df)))], , drop = FALSE]
        collection <- rep(NA_character_, nrow(df))
    }

    keys <- pathway_join_key(df, kegg_org)
    label <- pathway_display_label(df)
    # A row with no readable text still gets a bar; the key is a poorer label
    # than a name but a better one than "NA".
    label[is.na(label)] <- keys[is.na(label)]
    label[is.na(label)] <- "(unnamed pathway)"
    label <- ifelse(nchar(label) > 50, paste0(substr(label, 1, 47), "..."),
                    label)

    mixed <- by_collection && length(unique(collection)) > 1
    if (mixed) label <- paste0("[", collection, "] ", label)

    score <- suppressWarnings(as.numeric(df[[score_col$column]]))
    neg_log_p <- pmin(-log10(score + 1e-300), 15)

    # Coloured by collection whenever anything was classified, not only when
    # several were found: a GMT run that happens to yield GO terms alone should
    # not be drawn in the colour this file reserves for KEGG. The tag and the
    # legend are what a single collection does not need, since there is nothing
    # to tell apart.
    bar_col <- if (isTRUE(by_collection)) {
        unname(.ORA_COLLECTION_COLOURS[collection])
    } else {
        .ORA_COLLECTION_COLOURS[["KEGG"]]
    }

    p_label <- if (score_col$adjusted) "adjusted p-value" else "p-value"

    # with_par rather than a bare par(): these are called straight from tests
    # and from renderers that draw more than one figure to a device, and a
    # 17-line left margin left behind is not this function's to leave.
    withr::with_par(list(mar = c(5, 17, 3, 2)), {
        barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(label),
                las = 1, cex.names = 0.65, col = rev(bar_col),
                xlab = paste0("-log10(", p_label, ")"),
                main = title)
        abline(v = -log10(0.05), col = "red", lty = 2)
        if (mixed) {
            drawn <- unique(collection)
            legend("bottomright", legend = drawn, bty = "n", cex = 0.7,
                   fill = unname(.ORA_COLLECTION_COLOURS[drawn]))
        }
    })

    invisible(data.frame(key = keys, label = label, collection = collection,
                         score = score, stringsAsFactors = FALSE))
}


#' Plot multi-ORA dot plot comparing per-omics and pooled results
plot_multi_ora_dotplot <- function(combined, per_omics_ora, metab_ora, out_dir, top_n = 20) {

    if (is.null(combined) || nrow(combined) == 0) return(invisible(NULL))

    top <- combined[seq_len(min(top_n, nrow(combined))), ]

    # Collect per-omics data for these pathways
    plot_data <- list()

    # Pooled
    plot_data[["Pooled"]] <- data.frame(
        pathway = top$pathway,
        source = "Pooled",
        neg_log10_padj = -log10(top$pooled_padj + 1e-300),
        count = top$pooled_count,
        stringsAsFactors = FALSE
    )

    # Per-omics
    for (om in names(per_omics_ora)) {
        om_res <- per_omics_ora[[om]]
        if (is.null(om_res) || nrow(om_res) == 0) next
        padj_map <- setNames(om_res$padj, om_res$pathway)
        count_map <- setNames(om_res$Count, om_res$pathway)

        # unname: a lookup miss yields an NA-named element, which data.frame()
        # would otherwise promote to a (missing) row name and abort.
        matched_padj <- unname(padj_map[top$pathway])
        matched_count <- unname(count_map[top$pathway])

        plot_data[[om]] <- data.frame(
            pathway = top$pathway,
            source = gsub("_", " ", tools::toTitleCase(om)),
            neg_log10_padj = -log10(ifelse(is.na(matched_padj), 1, matched_padj) + 1e-300),
            count = ifelse(is.na(matched_count), 0, matched_count),
            stringsAsFactors = FALSE
        )
    }

    # Metabolomics compound ORA (join by normalized KEGG pathway ID)
    if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
        met_padj_map <- setNames(metab_ora$padj, normalize_kegg_pathway_id(metab_ora$ID))
        met_count_map <- setNames(metab_ora$setSize, normalize_kegg_pathway_id(metab_ora$ID))
        norm_top_ids <- normalize_kegg_pathway_id(top$ID)

        matched_padj <- unname(met_padj_map[norm_top_ids])
        matched_count <- unname(met_count_map[norm_top_ids])

        plot_data[["Metabolomics"]] <- data.frame(
            pathway = top$pathway,
            source = "Metabolomics",
            neg_log10_padj = -log10(ifelse(is.na(matched_padj), 1, matched_padj) + 1e-300),
            count = ifelse(is.na(matched_count), 0, matched_count),
            stringsAsFactors = FALSE
        )
    }

    plot_df <- do.call(rbind, plot_data)
    plot_df$neg_log10_padj <- pmin(plot_df$neg_log10_padj, 15)
    # Remove entries where pathway was not found in the per-omics results. Guard
    # against NA scores (pathways with no pooled/per-omic padj) — an NA in the
    # logical index yields NA row names and aborts the subset.
    plot_df <- plot_df[!is.na(plot_df$neg_log10_padj) & plot_df$neg_log10_padj > 0, ]

    if (nrow(plot_df) == 0) return(invisible(NULL))

    # Truncate long pathway names
    plot_df$pathway_short <- ifelse(
        nchar(plot_df$pathway) > 45,
        paste0(substr(plot_df$pathway, 1, 42), "..."),
        plot_df$pathway
    )

    # ggplot dot plot
    source_colors <- c(
        "Pooled" = "#7B2D8E",
        "Transcriptomics" = "#E41A1C",
        "Proteomics" = "#377EB8",
        "Metabolomics" = "#4DAF4A"
    )

    p <- ggplot2::ggplot(plot_df,
        ggplot2::aes(
            x = neg_log10_padj,
            y = stats::reorder(pathway_short, neg_log10_padj),
            color = source,
            size = count
        )) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::scale_color_manual(values = source_colors, name = "Source") +
        ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
        ggplot2::labs(
            title = "Multi-ORA: Pooled vs Per-Omics Enrichment",
            x = "-log10(padj)",
            y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    ggplot2::ggsave(
        file.path(out_dir, "multi_ora_dotplot.png"),
        plot = p, width = 12, height = max(6, 2 + top_n * 0.3), dpi = 300
    )
    message("  Saved multi-ORA dot plot")
}


#' Plot multi-ORA omics support barplot
#'
#' Shows how many omics layers support each enriched pathway.
plot_multi_ora_support <- function(combined, out_dir, top_n = 25) {
    if (is.null(combined) || nrow(combined) == 0) return(invisible(NULL))

    top <- combined[seq_len(min(top_n, nrow(combined))), ]
    top$label <- ifelse(nchar(top$pathway) > 45,
                         paste0(substr(top$pathway, 1, 42), "..."),
                         top$pathway)

    support_colors <- c("0" = "grey70", "1" = "#FDB863", "2" = "#E66101", "3" = "#B35806")

    p <- ggplot2::ggplot(top,
        ggplot2::aes(
            x = -log10(pooled_pvalue + 1e-300),
            y = stats::reorder(label, -pooled_pvalue),
            fill = factor(n_omics_support)
        )) +
        ggplot2::geom_col(alpha = 0.85) +
        ggplot2::scale_fill_manual(values = support_colors, name = "# Omics\nSupporting") +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
        ggplot2::labs(
            title = "Multi-ORA: Pathway Enrichment with Omics Support",
            x = "-log10(p-value)",
            y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    ggplot2::ggsave(
        file.path(out_dir, "multi_ora_support_barplot.png"),
        plot = p, width = 12, height = max(6, 2 + top_n * 0.25), dpi = 300
    )
    message("  Saved multi-ORA support barplot")
}


#' Generate pathview overlays for multi-omics-supported KEGG pathways
#'
#' For pathways enriched in >= min_support omics, renders KEGG pathway maps
#' with per-omics logFC coloring and compiles them into a single PDF.
#'
#' @param combined Multi-ORA summary table (from build_multi_ora_summary)
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result
#' @param config Full config
#' @param out_dir Output directory for pathview files
#' @param min_support Minimum n_omics_support to include (default 2)
#' @return Character path to compiled PDF, or NULL
generate_multi_ora_pathview <- function(combined, de_results, harmonization_res,
                                        config, out_dir, min_support = 2,
                                        top_n = 5) {

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("  Package 'pathview' not installed. Skipping pathway maps.")
        return(NULL)
    }

    if (is.null(combined) || !"n_omics_support" %in% colnames(combined)) return(NULL)

    # Read top_n from config if available
    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    # --- Prioritize pathways supported by both metabolomics AND proteomics ---
    has_met_padj <- "metabolomics_padj" %in% colnames(combined)
    has_prot_padj <- "proteomics_padj" %in% colnames(combined)
    has_met_pval <- "metabolomics_pvalue" %in% colnames(combined)
    has_prot_pval <- "proteomics_pvalue" %in% colnames(combined)

    supported <- NULL

    # Tier 1: both metabolomics + proteomics padj < 0.05
    if (has_met_padj && has_prot_padj) {
        tier1 <- combined[
            !is.na(combined$metabolomics_padj) & combined$metabolomics_padj < 0.05 &
            !is.na(combined$proteomics_padj) & combined$proteomics_padj < 0.05, ]
        if (nrow(tier1) > 0) {
            message("  Pathview: ", nrow(tier1),
                    " pathways supported by both metabolomics & proteomics (padj < 0.05)")
            supported <- tier1
        }
    }

    # Tier 2: both metabolomics + proteomics pvalue < 0.05
    if (is.null(supported) || nrow(supported) == 0) {
        if (has_met_pval && has_prot_pval) {
            tier2 <- combined[
                !is.na(combined$metabolomics_pvalue) & combined$metabolomics_pvalue < 0.05 &
                !is.na(combined$proteomics_pvalue) & combined$proteomics_pvalue < 0.05, ]
            if (nrow(tier2) > 0) {
                message("  Pathview: padj too strict; ", nrow(tier2),
                        " pathways with metabolomics & proteomics pvalue < 0.05")
                supported <- tier2
            }
        }
    }

    # Tier 3: general n_omics_support >= min_support (padj-based)
    if (is.null(supported) || nrow(supported) == 0) {
        supported <- combined[combined$n_omics_support >= min_support, ]
    }

    # Tier 4: general pvalue-based support
    if (nrow(supported) == 0 && "n_omics_support_pval" %in% colnames(combined)) {
        message("  No pathways with padj-based support >= ", min_support,
                ", falling back to pvalue < 0.05")
        supported <- combined[combined$n_omics_support_pval >= min_support, ]
    }

    # Tier 5: relax to single-omics support with pvalue
    if (nrow(supported) == 0 && min_support > 1) {
        message("  No pathways with pvalue-based support >= ", min_support,
                ", relaxing to single-omics support")
        if ("n_omics_support_pval" %in% colnames(combined)) {
            supported <- combined[combined$n_omics_support_pval >= 1, ]
        } else {
            supported <- combined[combined$n_omics_support >= 1, ]
        }
    }

    if (nrow(supported) == 0) {
        message("  No pathways with any omics support for pathview")
        return(NULL)
    }

    # Limit to top_n pathways sorted by pooled_pvalue
    supported <- supported[order(supported$pooled_pvalue, na.last = TRUE), ]
    supported <- head(supported, top_n)

    message("  Generating pathview for ", nrow(supported),
            " pathways (>= ", min_support, " omics support)")

    organism <- config$global$organism %||% "human"
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)

    if (is.null(kegg_org)) return(NULL)

    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Build per-contrast gene logFC and compound logFC ---
    # Extract DE tables and ID maps once per omics
    gene_de_tables <- list()
    gene_id_maps <- list()
    for (om in intersect(c("transcriptomics", "proteomics"), names(de_results))) {
        de_tables <- extract_de_tables(de_results[[om]], om, harmonization_res)
        if (is.null(de_tables) || length(de_tables) == 0) next
        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
            error = function(e) NULL
        )
        if (is.null(id_map) || nrow(id_map) == 0) next
        gene_de_tables[[om]] <- de_tables
        gene_id_maps[[om]] <- id_map
    }

    metab_de_tables <- NULL
    metab_id_map <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                              harmonization_res)
        metab_id_map <- tryCatch(
            map_metabolite_ids_to_kegg(metab_de_tables, harmonization_res),
            error = function(e) NULL
        )
    }

    # Determine contrast names (from first omics with multiple contrasts)
    contrast_names <- NULL
    for (om in names(gene_de_tables)) {
        nms <- names(gene_de_tables[[om]])
        if (!is.null(nms) && length(nms) > 0) {
            contrast_names <- nms
            break
        }
    }
    if (is.null(contrast_names)) {
        # Fallback to metabolomics contrast names
        if (!is.null(metab_de_tables)) contrast_names <- names(metab_de_tables)
    }
    if (is.null(contrast_names)) contrast_names <- "contrast_1"

    # pathview requires its 'bods' dataset in the global environment
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }

    cwd <- getwd()
    setwd(pv_dir)
    on.exit(setwd(cwd), add = TRUE)

    all_generated_pngs <- list()

    for (ci in seq_along(contrast_names)) {
        contrast <- contrast_names[ci]
        safe_contrast <- make.names(contrast)
        message("  Pathview for contrast: ", contrast)

        # Build gene logFC for this contrast
        gene_fc_list <- list()
        for (om in names(gene_de_tables)) {
            de_tbl <- gene_de_tables[[om]]
            idx <- min(ci, length(de_tbl))
            df <- de_tbl[[idx]]
            df_mapped <- merge(df, gene_id_maps[[om]], by = "feature_id")
            fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID, mean, na.rm = TRUE)
            fc_vec <- as.numeric(fc_arr)
            names(fc_vec) <- names(fc_arr)
            gene_fc_list[[om]] <- fc_vec
        }

        gene_data <- NULL
        if (length(gene_fc_list) == 2) {
            all_genes <- unique(c(names(gene_fc_list[[1]]), names(gene_fc_list[[2]])))
            gene_data <- matrix(NA, nrow = length(all_genes), ncol = 2,
                                dimnames = list(all_genes, names(gene_fc_list)))
            for (i in seq_along(gene_fc_list)) {
                fc <- gene_fc_list[[i]]
                gene_data[names(fc), i] <- fc
            }
        } else if (length(gene_fc_list) == 1) {
            gene_data <- gene_fc_list[[1]]
        }

        # Build compound logFC for this contrast
        cpd_data <- NULL
        if (!is.null(metab_id_map) && nrow(metab_id_map) > 0 &&
            !is.null(metab_de_tables) && length(metab_de_tables) > 0) {
            idx <- min(ci, length(metab_de_tables))
            df <- metab_de_tables[[idx]]
            df_mapped <- merge(df, metab_id_map, by = "feature_id")
            cpd_fc <- tapply(df_mapped$log2fc, df_mapped$KEGG_CPD, mean, na.rm = TRUE)
            cpd_data <- as.numeric(cpd_fc)
            names(cpd_data) <- names(cpd_fc)
        }

        if (is.null(gene_data) && is.null(cpd_data)) next

        # Run pathview per pathway for this contrast
        out_suffix <- paste0("multi_ora_", safe_contrast)
        contrast_pngs <- character(0)

        for (i in seq_len(nrow(supported))) {
            pid <- supported$ID[i]
            pw_name <- supported$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = out_suffix,
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                    same.layer = FALSE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".", out_suffix, ".multi.png"),
                    paste0(kegg_org, clean_pid, ".", out_suffix, ".png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    contrast_pngs <- c(contrast_pngs, file.path(pv_dir, found[1]))
                    message("    Pathview: ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(contrast_pngs) > 0) {
            all_generated_pngs[[contrast]] <- contrast_pngs
        }
    }

    if (length(all_generated_pngs) == 0) {
        message("  No pathview plots generated")
        return(NULL)
    }

    # --- Compile into a single PDF with contrast labels ---
    # "supported" is a claim: these pathways are enriched in two or more omics
    # layers. generate_per_omic_union_pathview(), the no-OrgDb fallback, unions
    # single-layer hits and so writes its own file -- sharing this name would
    # have presented one layer's evidence under this one's promise.
    pdf_path <- file.path(out_dir, "multi_ora_pathview_supported.pdf")
    tryCatch({
        grDevices::pdf(pdf_path, width = 12, height = 8)
        for (contrast in names(all_generated_pngs)) {
            pngs <- all_generated_pngs[[contrast]]
            # Title page for this contrast
            grid::grid.newpage()
            grid::grid.text(
                paste0("Contrast: ", contrast),
                gp = grid::gpar(fontsize = 24, fontface = "bold")
            )
            for (png_file in pngs) {
                img <- png::readPNG(png_file)
                grid::grid.newpage()
                # Add contrast label at top
                grid::grid.text(
                    contrast, x = 0.5, y = 0.98,
                    gp = grid::gpar(fontsize = 10, col = "grey40")
                )
                grid::grid.raster(img, y = 0.48, height = 0.92)
            }
        }
        grDevices::dev.off()
        total <- sum(lengths(all_generated_pngs))
        message("  Compiled ", total, " pathview maps (",
                length(all_generated_pngs), " contrasts) into: ", basename(pdf_path))
        pdf_path
    }, error = function(e) {
        message("  PDF compilation failed: ", e$message)
        tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        NULL
    })
}


#' Generate pathview plots for top pathways from each omics, overlaying both layers
#'
#' Two sets of plots:
#' 1. Top N metabolomics-enriched pathways with proteomics enzyme logFC overlaid
#' 2. Top N proteomics-enriched pathways with metabolomics compound logFC overlaid
#'
#' @param per_omics_ora Named list of per-omics ORA results
#' @param metab_ora Metabolomics compound ORA results
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result
#' @param config Full config
#' @param out_dir Output directory
#' @param top_n Number of top pathways per omics (default 5)
#' @return List with paths to compiled PDFs
generate_per_omics_pathview <- function(per_omics_ora, metab_ora, de_results,
                                        harmonization_res, config, out_dir,
                                        top_n = 5) {

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("  Package 'pathview' not installed. Skipping per-omics pathview.")
        return(NULL)
    }

    organism <- config$global$organism %||% "human"
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)
    if (is.null(kegg_org)) return(NULL)

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Build proteomics gene logFC vector ---
    gene_data <- NULL
    if ("proteomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$proteomics, "proteomics",
                                        harmonization_res)
        if (!is.null(de_tables) && length(de_tables) > 0) {
            id_map <- tryCatch(
                map_feature_ids_to_entrez(de_tables, "proteomics",
                                           harmonization_res, org_db),
                error = function(e) NULL
            )
            if (!is.null(id_map) && nrow(id_map) > 0) {
                df <- de_tables[[1]]
                df_mapped <- merge(df, id_map, by = "feature_id")
                fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID,
                                  mean, na.rm = TRUE)
                gene_data <- as.numeric(fc_arr)
                names(gene_data) <- names(fc_arr)
                message("  Pathview: ", length(gene_data),
                        " proteomics enzymes with logFC")
            }
        }
    }

    # Also try transcriptomics if available
    if (is.null(gene_data) && "transcriptomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$transcriptomics, "transcriptomics",
                                        harmonization_res)
        if (!is.null(de_tables) && length(de_tables) > 0) {
            id_map <- tryCatch(
                map_feature_ids_to_entrez(de_tables, "transcriptomics",
                                           harmonization_res, org_db),
                error = function(e) NULL
            )
            if (!is.null(id_map) && nrow(id_map) > 0) {
                df <- de_tables[[1]]
                df_mapped <- merge(df, id_map, by = "feature_id")
                fc_arr <- tapply(df_mapped$log2fc, df_mapped$ENTREZID,
                                  mean, na.rm = TRUE)
                gene_data <- as.numeric(fc_arr)
                names(gene_data) <- names(fc_arr)
                message("  Pathview: ", length(gene_data),
                        " transcriptomics genes with logFC")
            }
        }
    }

    # --- Build metabolomics compound logFC vector ---
    cpd_data <- NULL
    if ("metabolomics" %in% names(de_results)) {
        de_tables <- extract_de_tables(de_results$metabolomics, "metabolomics",
                                        harmonization_res)
        id_map <- tryCatch(
            map_metabolite_ids_to_kegg(de_tables, harmonization_res),
            error = function(e) NULL
        )
        if (!is.null(id_map) && nrow(id_map) > 0 && length(de_tables) > 0) {
            df <- de_tables[[1]]
            df_mapped <- merge(df, id_map, by = "feature_id")
            cpd_fc <- tapply(df_mapped$log2fc, df_mapped$KEGG_CPD,
                              mean, na.rm = TRUE)
            cpd_data <- as.numeric(cpd_fc)
            names(cpd_data) <- names(cpd_fc)
            message("  Pathview: ", length(cpd_data),
                    " metabolites with logFC")
        }
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        message("  No logFC data available for per-omics pathview")
        return(NULL)
    }

    # Ensure pathview bods dataset is loaded
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }

    cwd <- getwd()
    setwd(pv_dir)
    on.exit(setwd(cwd), add = TRUE)

    results <- list()

    # ------------------------------------------------------------------
    # Set 1: Top metabolomics pathways with proteomics enzyme overlay
    # ------------------------------------------------------------------
    if (!is.null(metab_ora) && nrow(metab_ora) > 0) {
        # Sort by pvalue, take top_n
        met_top <- metab_ora[order(metab_ora$pvalue), ]
        met_top <- head(met_top, top_n)

        message("  Pathview set 1: top ", nrow(met_top),
                " metabolomics pathways + proteomics overlay")

        met_pngs <- character(0)
        for (i in seq_len(nrow(met_top))) {
            pid <- met_top$ID[i]
            pw_name <- met_top$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = "metab_top",
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = FALSE,
                    same.layer  = TRUE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".metab_top.png"),
                    paste0(kegg_org, clean_pid, ".metab_top.multi.png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    met_pngs <- c(met_pngs, file.path(pv_dir, found[1]))
                    message("    ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(met_pngs) > 0) {
            pdf_path <- file.path(out_dir,
                                   "pathview_top_metabolomics_pathways.pdf")
            tryCatch({
                grDevices::pdf(pdf_path, width = 12, height = 8)
                for (png_file in met_pngs) {
                    img <- png::readPNG(png_file)
                    grid::grid.newpage()
                    grid::grid.raster(img)
                }
                grDevices::dev.off()
                message("  Compiled ", length(met_pngs),
                        " metabolomics pathview maps into: ",
                        basename(pdf_path))
                results$metabolomics_pdf <- pdf_path
            }, error = function(e) {
                message("  PDF compilation failed: ", e$message)
                tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            })
        }
    }

    # ------------------------------------------------------------------
    # Set 2: Top proteomics pathways with metabolomics compound overlay
    # ------------------------------------------------------------------
    prot_ora <- per_omics_ora[["proteomics"]]
    if (is.null(prot_ora) || nrow(prot_ora) == 0) {
        prot_ora <- per_omics_ora[["transcriptomics"]]
    }

    if (!is.null(prot_ora) && nrow(prot_ora) > 0) {
        prot_top <- prot_ora[order(prot_ora$pvalue), ]
        prot_top <- head(prot_top, top_n)

        message("  Pathview set 2: top ", nrow(prot_top),
                " proteomics pathways + metabolomics overlay")

        prot_pngs <- character(0)
        for (i in seq_len(nrow(prot_top))) {
            pid <- prot_top$ID[i]
            pw_name <- prot_top$pathway[i]
            clean_pid <- normalize_kegg_pathway_id(pid)

            tryCatch({
                pathview::pathview(
                    gene.data  = gene_data,
                    cpd.data   = cpd_data,
                    pathway.id = clean_pid,
                    species    = kegg_org,
                    out.suffix = "prot_top",
                    kegg.dir   = pv_dir,
                    keys.align = "y",
                    match.data = TRUE,
                    multi.state = FALSE,
                    same.layer  = TRUE
                )

                candidates <- c(
                    paste0(kegg_org, clean_pid, ".prot_top.png"),
                    paste0(kegg_org, clean_pid, ".prot_top.multi.png")
                )
                found <- candidates[file.exists(candidates)]
                if (length(found) > 0) {
                    prot_pngs <- c(prot_pngs, file.path(pv_dir, found[1]))
                    message("    ", pw_name, " (", pid, ")")
                }
            }, error = function(e) {
                message("    Pathview failed for ", pid, ": ", e$message)
            })
        }

        if (length(prot_pngs) > 0) {
            pdf_path <- file.path(out_dir,
                                   "pathview_top_proteomics_pathways.pdf")
            tryCatch({
                grDevices::pdf(pdf_path, width = 12, height = 8)
                for (png_file in prot_pngs) {
                    img <- png::readPNG(png_file)
                    grid::grid.newpage()
                    grid::grid.raster(img)
                }
                grDevices::dev.off()
                message("  Compiled ", length(prot_pngs),
                        " proteomics pathview maps into: ",
                        basename(pdf_path))
                results$proteomics_pdf <- pdf_path
            }, error = function(e) {
                message("  PDF compilation failed: ", e$message)
                tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            })
        }
    }

    if (length(results) == 0) {
        message("  No per-omics pathview plots generated")
        return(NULL)
    }

    results
}
