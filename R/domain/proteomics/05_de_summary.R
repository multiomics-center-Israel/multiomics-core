#' Summarize differential expression across multiple imputations (legacy-compatible)
summarize_limma_mult_imputation <- function(runs_de_tables, config) {
    de_cfg <- config$modes$proteomics$de
    imp_cfg <- config$modes$proteomics$imputation

    NO_REPETITIONS <- as.integer(imp_cfg$no_repetitions)
    MIN_NO_PASSED <- as.integer(imp_cfg$min_no_passed)

    use_adj_for_pass1 <- isTRUE(de_cfg$use_adj_for_pass1)
    p_cutoff <- as.numeric(de_cfg$p_cutoff)

    linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff)
    lfc_cutoff <- log2(linear_fc_cutoff)

    stopifnot(length(runs_de_tables) == NO_REPETITIONS)
    stopifnot(MIN_NO_PASSED >= 1, MIN_NO_PASSED <= NO_REPETITIONS)

    q <- MIN_NO_PASSED / NO_REPETITIONS

    contrasts <- names(runs_de_tables[[1]])
    stopifnot(length(contrasts) > 0)

    de_table_cfg <- config$modes$proteomics$de_table %||% list()
    id_col <- de_table_cfg$id_col %||% "FeatureID"

    # Optional columns that might exist
    extra_cols <- c("Protein.Names", "Genes", "First.Protein.Description")

    ref_df <- runs_de_tables[[1]][[contrasts[1]]]

    # Check if primary ID column exists
    if (!id_col %in% colnames(ref_df)) {
        stop(sprintf(
            "Configured ID column '%s' not found in DE results. Available: %s",
            id_col, paste(colnames(ref_df), collapse = ", ")
        ))
    }

    # Collect available columns
    available_cols <- c(id_col, intersect(extra_cols, colnames(ref_df)))

    out <- ref_df[, available_cols, drop = FALSE]
    ref_ids <- ref_df[[id_col]]

    # Validate alignment
    for (n in seq_len(NO_REPETITIONS)) {
        for (cn in contrasts) {
            cur <- runs_de_tables[[n]][[cn]]
            stopifnot(nrow(cur) == length(ref_ids))

            # Use dynamic ID column for check
            if (!all(cur[[id_col]] == ref_ids)) {
                stop(sprintf("Mismatch in ID column '%s' between imputations/contrasts", id_col))
            }
        }
    }

    for (cn in contrasts) {
        contrast_print <- gsub(" ", "", cn)

        logfc_mat <- sapply(seq_len(NO_REPETITIONS), function(n) runs_de_tables[[n]][[cn]][["logFC"]])
        p_mat <- sapply(seq_len(NO_REPETITIONS), function(n) runs_de_tables[[n]][[cn]][["P.Value"]])
        padj_mat <- sapply(seq_len(NO_REPETITIONS), function(n) runs_de_tables[[n]][[cn]][["adj.P.Val"]])

        pass1_mat <- sapply(seq_len(NO_REPETITIONS), function(n) {
            mark_pass1(runs_de_tables[[n]][[cn]],
                p_cutoff    = p_cutoff,
                lfc_cutoff  = lfc_cutoff,
                use_adj     = use_adj_for_pass1
            )
        })

        sum_pass <- rowSums(pass1_mat, na.rm = TRUE)

        pass_imputs <- ifelse(sum_pass >= MIN_NO_PASSED, 1, NA)

        linearRatio_imputs <- rowMeans(2^logfc_mat, na.rm = TRUE)
        linearFC_imputs <- ifelse(linearRatio_imputs >= 1, linearRatio_imputs, -1 / linearRatio_imputs)

        pvalue_imputs <- apply(p_mat, 1, quantile, probs = q, na.rm = TRUE)
        padj_imputs <- apply(padj_mat, 1, quantile, probs = q, na.rm = TRUE)

        out[[paste0("sum.pass.", contrast_print)]] <- sum_pass
        out[[paste0("pass.imputs.", contrast_print)]] <- pass_imputs
        out[[paste0("linearRatio.imputs.", contrast_print)]] <- linearRatio_imputs
        out[[paste0("linearFC.imputs.", contrast_print)]] <- signif(linearFC_imputs, 4)
        out[[paste0("pvalue.imputs.", contrast_print)]] <- pvalue_imputs
        out[[paste0("padj.imputs.", contrast_print)]] <- padj_imputs
    }

    out <- add_pass_any_contrast(out)
    return(out)
}

#' Build DE contrast summary statistics
#'
#' Creates a summary data.frame showing for each contrast:
#' - total_genes: Total number of genes/proteins analyzed
#' - n_significant: Number passing DE thresholds
#' - n_up: Number of upregulated (passed + positive FC)
#' - n_down: Number of downregulated (passed + negative FC)
#'
#' @param summary_df The summary dataframe from summarize_limma_mult_imputation()
#' @return A data.frame with one row per contrast
#' @export
build_de_contrast_summary <- function(summary_df) {
    if (is.null(summary_df) || nrow(summary_df) == 0) {
        return(NULL)
    }

    # Find all contrast names from pass columns
    pass_cols <- grep("^pass\\.imputs\\.", names(summary_df), value = TRUE)
    if (length(pass_cols) == 0) {
        return(NULL)
    }

    # Extract contrast names
    contrasts <- sub("^pass\\.imputs\\.", "", pass_cols)

    total_genes <- nrow(summary_df)

    results <- lapply(contrasts, function(cn) {
        pass_col <- paste0("pass.imputs.", cn)
        fc_col <- paste0("linearFC.imputs.", cn)

        # Check columns exist
        if (!pass_col %in% names(summary_df)) {
            return(data.frame(
                Contrast = cn,
                total_genes = total_genes,
                n_significant = NA_integer_,
                n_up = NA_integer_,
                n_down = NA_integer_,
                stringsAsFactors = FALSE
            ))
        }

        passed <- !is.na(summary_df[[pass_col]]) & summary_df[[pass_col]] == 1
        n_significant <- sum(passed, na.rm = TRUE)

        if (fc_col %in% names(summary_df)) {
            fc_vals <- summary_df[[fc_col]]
            n_up <- sum(passed & !is.na(fc_vals) & fc_vals > 0, na.rm = TRUE)
            n_down <- sum(passed & !is.na(fc_vals) & fc_vals < 0, na.rm = TRUE)
        } else {
            n_up <- NA_integer_
            n_down <- NA_integer_
        }

        data.frame(
            Contrast = cn,
            total_genes = total_genes,
            n_significant = n_significant,
            n_up = n_up,
            n_down = n_down,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, results)
}

#' Mark differential expression pass for a single result table
mark_pass1 <- function(de_table, p_cutoff, lfc_cutoff, use_adj = TRUE) {
    pcol <- if (isTRUE(use_adj)) "adj.P.Val" else "P.Value"
    stopifnot(all(c("logFC", pcol) %in% colnames(de_table)))

    pass <- (de_table[[pcol]] <= p_cutoff) & (abs(de_table[["logFC"]]) >= lfc_cutoff)
    ifelse(pass, 1, NA)
}

#' Helper to add pass_any_contrast column
add_pass_any_contrast <- function(summary_df, pass_prefix = "^pass\\.imputs\\.", out_col = "pass_any_contrast", out_n_col = "n_pass_contrasts") {
    pass_cols <- grep(pass_prefix, names(summary_df), value = TRUE)

    if (length(pass_cols) == 0) {
        summary_df[[out_n_col]] <- 0L
        summary_df[[out_col]] <- NA
        return(summary_df)
    }

    pass_mat <- as.matrix(summary_df[, pass_cols, drop = FALSE])
    n_pass <- rowSums(!is.na(pass_mat) & pass_mat == 1, na.rm = TRUE)

    summary_df[[out_n_col]] <- as.integer(n_pass)
    summary_df[[out_col]] <- ifelse(n_pass > 0, 1, NA)

    summary_df
}

#' Run limma differential analysis for proteomics with contrast support
#'
#' Fits a limma linear model on an imputed proteomics expression matrix and
#' returns per-contrast result tables with feature annotations.
#'
#' @return A list with aligned metadata, design matrix, contrasts, fitted model and per-contrast DE tables.
#' @export
run_limma_proteomics <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg) {
    p_cfg <- cfg$modes$proteomics

    sample_col <- p_cfg$effects$samples %||% "SampleID"
    group_col <- p_cfg$effects$color %||% "Condition"

    protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
    default_annot <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
    annot_cols <- unique(c(protein_id_col, p_cfg$id_columns$protein_annot %||% default_annot))

    assert_numeric_matrix(expr_imp, "expr_imp")
    meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)

    meta_aligned[[group_col]] <- factor(meta_aligned[[group_col]])
    orig_levels <- levels(meta_aligned[[group_col]])
    safe_levels <- make.names(orig_levels)
    level_map <- setNames(orig_levels, safe_levels)  # safe -> original
    levels(meta_aligned[[group_col]]) <- safe_levels

    design <- model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)), data = meta_aligned)
    colnames(design) <- safe_levels

    stopifnot(all(contrasts_df$Factor == group_col))
    # Map contrast numerator/denominator to safe level names
    safe_num <- make.names(contrasts_df$Numerator)
    safe_den <- make.names(contrasts_df$Denominator)
    contrast_formulas <- setNames(
        paste(safe_num, safe_den, sep = " - "),
        contrasts_df$Contrast_name
    )

    contrast_matrix <- limma::makeContrasts(contrasts = contrast_formulas, levels = design)
    colnames(contrast_matrix) <- names(contrast_formulas)

    fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_imp, design), contrast_matrix))

    ann <- align_annotations_to_expr(expr_imp, prot_tbl, protein_id_col, annot_cols)
    feature_id <- ann[[protein_id_col]]
    annot_out <- setdiff(annot_cols, protein_id_col)

    de_table_cfg <- p_cfg$de_table %||% list()
    target_id_col <- de_table_cfg$id_col %||% "FeatureID"

    de_tables <- lapply(colnames(contrast_matrix), function(cn) {
        de <- limma::topTable(fit2, coef = cn, adjust.method = "BH", sort.by = "none", number = Inf)
        de <- align_de_to_expr(de, expr_imp, contrast_name = cn)
        df_out <- data.frame(
            TEMP_ID_COL = feature_id,
            Contrast = cn,
            ann[, annot_out, drop = FALSE],
            de[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
            check.names = FALSE
        )
        # Rename the temporal ID column to the actual configured ID column (target for DE table)
        names(df_out)[names(df_out) == "TEMP_ID_COL"] <- target_id_col
        df_out
    })

    # IMPORTANT: names(de_tables) are used by summarize_limma_mult_imputation()
    names(de_tables) <- colnames(contrast_matrix)

    list(
        expr_imp = expr_imp,
        meta_aligned = meta_aligned,
        design = design,
        contrast_formulas = contrast_formulas,
        contrast_matrix = contrast_matrix,
        fit2 = fit2,
        de_tables = de_tables
    )
}

#' Run limma proteomics DE on a precomputed list of imputed datasets
run_limma_multimp <- function(imputations, meta, contrasts_df, prot_tbl, cfg, verbose = FALSE) {
    validate_proteomics_imputations(imputations = imputations, meta = meta, cfg = cfg)

    runs <- vector("list", length(imputations))

    for (i in seq_along(imputations)) {
        if (isTRUE(verbose)) message(sprintf("Limma on imputations: %d / %d", i, length(imputations)))

        runs[[i]] <- run_limma_proteomics(
            expr_imp     = imputations[[i]],
            meta         = meta,
            contrasts_df = contrasts_df,
            prot_tbl     = prot_tbl,
            cfg          = cfg
        )
    }
    runs
}

#' Extract DE feature IDs from a DE result object
#'
#' Omics-agnostic helper that returns feature identifiers passing
#' differential expression thresholds as defined in config.
#'
#' @param de_res DE result object containing a summary table
#' @param cfg Global config
#'
#' @return character vector of feature IDs
get_de_features <- function(de_res, cfg) {
    stopifnot(is.list(de_res), !is.null(de_res$summary_df))
    df <- de_res$summary_df
    stopifnot(is.data.frame(df))

    # Accept legacy + new column names (with safe defaults)
    id_col <- cfg$de_table$id_col %||% "FeatureID"
    pass_any_col <- cfg$de_table$pass_any_col %||% "pass_any_contrast"

    if (is.na(id_col) || is.null(id_col)) {
        stop("DE summary missing id column. Expected : ", id_col)
    }

    if (!id_col %in% colnames(df)) {
        # Fallback or strict error?
        # Given summary_df creation, FeatureID is standard.
        # If strictly missing, we stop.
        stop("Column '", id_col, "' not found in DE summary dataframe.")
    }

    # Standardize internally
    df$feature_id <- df[[id_col]]

    if (pass_any_col %in% colnames(df)) {
        feats <- df$feature_id[!is.na(df[[pass_any_col]]) & df[[pass_any_col]] == 1]
        return(unique(as.character(feats)))
    } else {
        warning("pass_any_col '", pass_any_col, "' not found. Returning empty feature set.")
        return(character(0))
    }
}
