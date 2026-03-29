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
        contrast_print <- normalize_contrast_name(cn)

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
        out[[paste0("linearFC.imputs.", contrast_print)]] <- signif(linearFC_imputs, 3)
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
    group_col <- p_cfg$de_table$group_col %||% p_cfg$effects$color %||% "Condition"
    

    protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
    default_annot <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
    annot_cols <- unique(c(protein_id_col, p_cfg$id_columns$protein_annot %||% default_annot))

    assert_numeric_matrix(expr_imp, "expr_imp")
    expr_imp <- align_matrix_to_meta(expr_imp, meta, sample_col)
    meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)

    meta_aligned[[group_col]] <- factor(meta_aligned[[group_col]])
    orig_levels <- levels(meta_aligned[[group_col]])
    safe_levels <- make.names(orig_levels)
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


#' Load pre-computed proteomics DE tables from config$modes$proteomics$files$de_table
#'
#' Reads CSV files and builds a summary_df + runs_de_tables structure
#' compatible with the proteomics DE contract used by downstream modules.
#'
#' @param config Full pipeline config.
#' @param contrasts_df Optional contrasts data frame (with Contrast_name column).
#'   When provided, contrast names from this table are used instead of filenames.
#' @return List matching mod_proteomics_de() output: method, summary_df,
#'   runs_de_tables, de_model, imputations, runs.
load_precomputed_proteomics_de <- function(config, contrasts_df = NULL) {
    cfg <- config$modes$proteomics
    de_cfg <- cfg$de %||% list()

    de_files <- cfg$files$de_table
    if (is.list(de_files)) de_files <- unlist(de_files)

    padj_cutoff <- as.numeric(de_cfg$p_cutoff %||% 0.05)
    linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff %||% 1.5)
    lfc_cutoff <- log2(linear_fc_cutoff)

    de_table_cfg <- cfg$de_table %||% list()
    id_col <- de_table_cfg$id_col %||% "FeatureID"

    # Use contrast names from contrasts_df when available (must match file count)
    if (!is.null(contrasts_df) && "Contrast_name" %in% colnames(contrasts_df) &&
        nrow(contrasts_df) == length(de_files)) {
        contrast_labels <- as.character(contrasts_df$Contrast_name)
    } else {
        contrast_labels <- vapply(de_files, function(f) {
            bn <- tools::file_path_sans_ext(basename(f))
            sub("^de_", "", bn)
        }, character(1), USE.NAMES = FALSE)
    }

    # Load per-contrast tables
    per_contrast <- list()
    for (i in seq_along(de_files)) {
        abs_path <- resolve_raw_path(config, de_files[i])
        if (!file.exists(abs_path)) {
            stop("Pre-computed proteomics DE table not found: ", abs_path)
        }

        raw <- read_table_auto(abs_path)
        cn <- colnames(raw)

        # Feature IDs
        feat_col <- cn[cn %in% c(id_col, "FeatureID", "Protein.Group",
                                  "protein_id", "feature_id")][1]
        if (is.na(feat_col)) {
            unnamed_idx <- match(TRUE, cn %in% c("...1", "", "X", "V1"))
            feat_ids <- if (!is.na(unnamed_idx)) as.character(raw[[unnamed_idx]]) else rownames(raw)
        } else {
            feat_ids <- as.character(raw[[feat_col]])
        }

        # logFC
        lfc_col <- cn[cn %in% c("logFC", "log2FoldChange", "log2FC",
                                 "log2(FC)", "log2.FC.")][1]
        lfc_vals <- if (!is.na(lfc_col)) as.numeric(raw[[lfc_col]]) else NA_real_

        # P.Value
        pval_col <- cn[cn %in% c("P.Value", "pvalue", "PValue", "p.value",
                                  "raw.pval")][1]
        pval_vals <- if (!is.na(pval_col)) as.numeric(raw[[pval_col]]) else NA_real_

        # adj.P.Val
        padj_col_name <- cn[cn %in% c("adj.P.Val", "padj", "FDR", "q.value",
                                       "p.adjust", "qvalue")][1]
        padj_vals <- if (!is.na(padj_col_name)) {
            as.numeric(raw[[padj_col_name]])
        } else {
            stats::p.adjust(pval_vals, method = "BH")
        }

        # AveExpr
        ave_col <- cn[cn %in% c("AveExpr", "baseMean", "logCPM")][1]
        ave_vals <- if (!is.na(ave_col)) as.numeric(raw[[ave_col]]) else NA_real_

        tbl <- data.frame(
            logFC     = lfc_vals,
            AveExpr   = ave_vals,
            t         = NA_real_,
            P.Value   = pval_vals,
            adj.P.Val = padj_vals,
            B         = NA_real_,
            stringsAsFactors = FALSE
        )
        tbl[[id_col]] <- feat_ids

        # Carry annotation columns if present
        for (a in c("Protein.Names", "Genes", "First.Protein.Description")) {
            if (a %in% cn) tbl[[a]] <- as.character(raw[[a]])
        }

        per_contrast[[contrast_labels[i]]] <- tbl
        message("  Loaded ", nrow(tbl), " features from ", basename(de_files[i]),
                " (label: ", contrast_labels[i], ")")
    }

    # Build summary_df matching the imputation-based naming convention
    all_features <- unique(unlist(lapply(per_contrast, function(t) t[[id_col]])))
    out <- data.frame(x__ = all_features, stringsAsFactors = FALSE)
    names(out)[1] <- id_col

    # Carry annotation columns from first table
    ref <- per_contrast[[1]]
    for (a in intersect(c("Protein.Names", "Genes", "First.Protein.Description"),
                        colnames(ref))) {
        idx <- match(all_features, ref[[id_col]])
        out[[a]] <- ref[[a]][idx]
    }

    for (ctr in names(per_contrast)) {
        tbl <- per_contrast[[ctr]]
        idx <- match(all_features, tbl[[id_col]])
        contrast_print <- normalize_contrast_name(ctr)

        lfc <- tbl$logFC[idx]
        linear_ratio <- 2^abs(lfc)
        linear_fc <- ifelse(lfc >= 0, linear_ratio, -linear_ratio)

        is_sig <- !is.na(tbl$adj.P.Val[idx]) &
                  tbl$adj.P.Val[idx] < padj_cutoff &
                  abs(lfc) >= lfc_cutoff
        pass <- ifelse(is_sig, 1, NA)

        out[[paste0("sum.pass.", contrast_print)]]          <- as.integer(!is.na(pass) & pass == 1)
        out[[paste0("pass.imputs.", contrast_print)]]       <- pass
        out[[paste0("linearRatio.imputs.", contrast_print)]] <- linear_ratio
        out[[paste0("linearFC.imputs.", contrast_print)]]   <- signif(linear_fc, 3)
        out[[paste0("pvalue.imputs.", contrast_print)]]     <- tbl$P.Value[idx]
        out[[paste0("padj.imputs.", contrast_print)]]       <- tbl$adj.P.Val[idx]
    }

    out <- add_pass_any_contrast(out)

    # Build runs_de_tables (single "imputation" with our loaded tables)
    runs_de_tables <- list(per_contrast)

    message("Proteomics DE: loaded ", length(per_contrast), " pre-computed contrast(s), ",
            sum(out$pass_any_contrast == 1, na.rm = TRUE), " significant features")

    list(
        method         = "precomputed",
        imputations    = NULL,
        runs           = NULL,
        runs_de_tables = runs_de_tables,
        summary_df     = out,
        de_model       = NULL
    )
}
