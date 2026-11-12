#' Batch correction of the RNA-seq display matrices
#'
#' Optional `modes.rna.batch_correction` block. It is independent of
#' `batch_in_model` ([get_rna_batch_model_config()]), which changes the DE
#' statistics only. This block changes expression values only, and only in a
#' separate copy:
#'
#'   * The input is the continuous log-scale analysis matrix (`pre$expr_work`:
#'     VST, TMM logCPM or preprocessed input) and the QC matrix
#'     ([rna_qc_expr()]). Raw counts are never corrected.
#'   * `pre` is not modified. The corrected matrices are returned as a separate
#'     object and are never DE input.
#'   * The biological group the DE compares is kept through the design of both
#'     backends ([correct_batch_limma()], [correct_batch_combat()]).
#'
#' The corrected values are for exploratory display and clustering. They do not
#' show that the variation that remains is biological.
#'
#' What ran is recorded by [compute_rna_batch_correction()] and written to
#' `Datasets/batch_correction.tsv` by [write_rna_batch_correction()], so text
#' that describes the correction reads that record, not the config.

#' Resolve the batch-correction config
#'
#' The one reading of `modes.rna.batch_correction`, used by the correction and
#' the config validation.
#'
#' @param rna_cfg The `config$modes$rna` list.
#' @return A list with `enabled` (logical), `method` (`"limma"` or `"combat"`;
#'   `NA` when disabled), `batch_cols` (character; empty when disabled) and
#'   `mean_only` (logical; `NA` unless the method is ComBat).
#' @examples
#' get_rna_batch_correction_config(list())  # disabled
#' get_rna_batch_correction_config(list(batch_correction = list(
#'     enabled = TRUE, batch_cols = "Batch")))  # limma
get_rna_batch_correction_config <- function(rna_cfg) {
    bc <- rna_cfg$batch_correction %||% list()
    if (!isTRUE(bc$enabled)) {
        return(list(enabled = FALSE, method = NA_character_,
                    batch_cols = character(0), mean_only = NA))
    }

    method <- bc$method %||% "limma"
    if (!is.character(method) || length(method) != 1 ||
        !method %in% c("limma", "combat")) {
        stop("[config] batch_correction$method must be \"limma\" or \"combat\" (got: ",
             paste(format(method), collapse = ", "), ").", call. = FALSE)
    }

    batch_cols <- as.character(unlist(bc$batch_cols %||% character(0)))
    if (!length(batch_cols) %in% 1:2 || !all(nzchar(batch_cols)) ||
        anyDuplicated(batch_cols)) {
        stop("[config] batch_correction$batch_cols must name one or two distinct metadata ",
             "columns when batch_correction$enabled is true (got: ",
             if (length(batch_cols)) paste(batch_cols, collapse = ", ") else "none", ").",
             call. = FALSE)
    }
    if (identical(method, "combat") && length(batch_cols) != 1) {
        stop("[config] batch_correction$method \"combat\" takes exactly one batch column ",
             "(got: ", paste(batch_cols, collapse = ", "), "). Use method \"limma\" for two ",
             "batch columns, or name one.", call. = FALSE)
    }

    list(
        enabled    = TRUE,
        method     = method,
        batch_cols = batch_cols,
        mean_only  = if (identical(method, "combat")) isTRUE(bc$combat_mean_only) else NA
    )
}

#' The group column the RNA DE compares
#'
#' The contrasts' `Factor` when a contrasts table exists. Without one,
#' `mod_rnaseq_de()` generates the contrasts from `effects$color`, then
#' `filtering$group_col`, then `"Condition"`; this follows the same order, so
#' the correction keeps the group that the DE tests.
#'
#' @param rna_cfg The `config$modes$rna` list.
#' @param contrasts_df The contrasts table (`pre$contrasts`), or `NULL`.
#' @return A single column name.
#' @examples
#' resolve_rna_de_factor_col(list(), data.frame(Factor = c("Condition", "Condition")))
#' resolve_rna_de_factor_col(list(effects = list(color = "Group")), NULL)
resolve_rna_de_factor_col <- function(rna_cfg, contrasts_df = NULL) {
    if (!is.null(contrasts_df) && nrow(contrasts_df) > 0) {
        factor_col <- unique(as.character(contrasts_df$Factor))
        if (length(factor_col) != 1) {
            stop("[batch_correction] All contrasts must use the same factor column (found: ",
                 paste(factor_col, collapse = ", "), ").", call. = FALSE)
        }
        return(factor_col)
    }
    group_col <- rna_cfg$effects$color %||% rna_cfg$filtering$group_col %||% "Condition"
    as.character(group_col)[1]
}

#' Check the batch and group columns and align them to the matrix
#'
#' Run on the samples of the matrix to correct, before any correction. Stops
#' on a missing, non-syntactic or NA-holding batch column, a batch column that
#' is the group column, a batch with fewer than two levels, a sample without a
#' group, and batch fully confounded with the group
#' ([assert_batch_not_confounded()]).
#'
#' Batch columns are converted with `factor(as.character(x))`, so numeric
#' labels (e.g. pair numbers 1, 2, 3) are categories.
#'
#' @param meta Sample metadata.
#' @param sample_ids Column names of the matrix to correct.
#' @param sample_col Metadata column holding sample IDs; the row names are used
#'   when the column is absent.
#' @param group_col The group column to keep ([resolve_rna_de_factor_col()]).
#' @param batch_cols One or two batch column names.
#' @return A list with `group` (factor) and `batch` (a list of one or two
#'   factors), each in the order of `sample_ids`.
prepare_rna_batch_correction <- function(meta, sample_ids, sample_col, group_col,
                                         batch_cols) {
    missing_cols <- setdiff(c(batch_cols, group_col), colnames(meta))
    if (length(missing_cols) > 0) {
        stop("[batch_correction] Column(s) not found in metadata: ",
             paste(missing_cols, collapse = ", "), ". Available columns: ",
             paste(colnames(meta), collapse = ", "), ".", call. = FALSE)
    }
    if (group_col %in% batch_cols) {
        stop("[batch_correction] Batch column '", group_col, "' is the group column the ",
             "contrasts compare; removing it would remove the biology.", call. = FALSE)
    }
    bad_names <- batch_cols[make.names(batch_cols) != batch_cols]
    if (length(bad_names) > 0) {
        stop("[batch_correction] Batch column name(s) ", paste(bad_names, collapse = ", "),
             " are not syntactic. Rename the column(s) in the metadata using letters, ",
             "digits, '.' and '_' only.", call. = FALSE)
    }

    meta_ids <- if (sample_col %in% colnames(meta)) {
        as.character(meta[[sample_col]])
    } else {
        rownames(meta)
    }
    rows <- match(sample_ids, meta_ids)
    if (anyNA(rows)) {
        stop("[batch_correction] ", sum(is.na(rows)), " sample(s) of the expression matrix ",
             "are not in the metadata: ",
             paste(utils::head(sample_ids[is.na(rows)], 10), collapse = ", "), ".",
             call. = FALSE)
    }

    values_of <- function(col, what) {
        vals <- as.character(meta[[col]])[rows]
        bad <- is.na(vals) | !nzchar(vals)
        if (any(bad)) {
            stop("[batch_correction] ", what, " column '", col, "' has no value for ",
                 sum(bad), " sample(s): ",
                 paste(utils::head(sample_ids[bad], 10), collapse = ", "),
                 ". Every sample needs a value.", call. = FALSE)
        }
        factor(vals)
    }

    group <- values_of(group_col, "Group")
    batch <- lapply(batch_cols, function(b) {
        f <- values_of(b, "Batch")
        if (nlevels(f) < 2) {
            stop("[batch_correction] Batch column '", b, "' has a single level among the ",
                 "samples, so there is no batch effect to remove. Turn batch_correction ",
                 "off or choose another column.", call. = FALSE)
        }
        f
    })
    names(batch) <- batch_cols

    assert_batch_not_confounded(
        batch_var = batch[[1]],
        group_var = group,
        batch2    = if (length(batch) == 2) batch[[2]]
    )

    list(group = group, batch = batch)
}

#' Features that are constant within a batch
#'
#' ComBat standardizes each feature by its within-batch spread, which is
#' undefined for a feature with the same value in every sample of a batch
#' (e.g. a gene with no reads in that batch).
#'
#' @param expr Numeric matrix (features x samples).
#' @param batch Factor of batch labels, one per column.
#' @return Logical vector, `TRUE` for a feature that is constant within at
#'   least one batch of two or more samples.
#' @examples
#' m <- rbind(a = c(1, 1, 2, 3), b = c(1, 2, 3, 4))
#' rna_constant_within_batch(m, factor(c("B1", "B1", "B2", "B2")))  # a: TRUE
rna_constant_within_batch <- function(expr, batch) {
    constant <- rep(FALSE, nrow(expr))
    for (lv in levels(batch)) {
        idx <- which(batch == lv)
        if (length(idx) < 2) next
        sub <- expr[, idx, drop = FALSE]
        constant <- constant | (rowSums(sub != sub[, 1]) == 0)
    }
    unname(constant)
}

#' Batch-correct one RNA matrix
#'
#' limma corrects every feature. ComBat corrects the features that vary within
#' every batch and leaves the others as they were
#' ([rna_constant_within_batch()]); their number is returned so that it can be
#' reported. This is done here rather than left to `sva`, so that the result
#' does not depend on how the installed version treats such features.
#'
#' @param expr Numeric log-scale matrix (features x samples).
#' @param design Output of [prepare_rna_batch_correction()] for the columns of
#'   `expr`.
#' @param bc Output of [get_rna_batch_correction_config()], enabled.
#' @param seed Integer seed for ComBat.
#' @return A list with `expr` (the corrected matrix, same dim and dimnames as
#'   the input) and `n_unadjusted` (features left as they were).
correct_rna_batch_matrix <- function(expr, design, bc, seed = 42) {
    expr <- as.matrix(expr)
    attr(expr, "method") <- NULL
    n_unadjusted <- 0L

    if (identical(bc$method, "limma")) {
        corrected <- correct_batch_limma(
            expr,
            batch_var = design$batch[[1]],
            batch2    = if (length(design$batch) == 2) design$batch[[2]],
            group_var = design$group
        )
    } else {
        constant <- rna_constant_within_batch(expr, design$batch[[1]])
        n_unadjusted <- sum(constant)
        if (sum(!constant) < 2) {
            stop("[batch_correction] Fewer than two features vary within every batch, so ",
                 "ComBat has nothing to estimate from. Use method \"limma\".", call. = FALSE)
        }
        corrected <- expr
        corrected[!constant, ] <- correct_batch_combat(
            expr[!constant, , drop = FALSE],
            batch_var = design$batch[[1]],
            group_var = design$group,
            mean_only = isTRUE(bc$mean_only),
            seed      = seed
        )
    }

    if (!identical(dim(corrected), dim(expr)) || !all(is.finite(corrected))) {
        stop("[batch_correction] The ", bc$method, " correction returned a matrix of another ",
             "size or with non-finite values, so it is not used. Try method \"limma\"",
             if (identical(bc$method, "combat")) " or combat_mean_only: true" else "",
             ".", call. = FALSE)
    }
    dimnames(corrected) <- dimnames(expr)
    list(expr = corrected, n_unadjusted = as.integer(n_unadjusted))
}

#' File name of the exported batch-corrected matrix
#'
#' Deliberately not `rna_norm_*`: the report and the fact sheet identify the
#' normalized table by that prefix.
#'
#' @param method Correction method that ran.
#' @param input_matrix What produced the corrected matrix's input
#'   (`pre$info$expr_work_method`).
#' @return A file name.
rna_batch_corr_file_name <- function(method, input_matrix) {
    sprintf("rna_batch_corr_%s_%s.tsv", method, input_matrix)
}

#' Compute the batch-corrected RNA matrices
#'
#' Corrects the analysis matrix (`pre$expr_work`) and the QC matrix
#' ([rna_qc_expr()]) with the same backend and the same design. For VST and
#' preprocessed input the two are one matrix and it is corrected once. `pre` is
#' not modified.
#'
#' @param pre List returned by `preprocess_rna()`.
#' @param config Full config list.
#' @return A list with `ran` (logical). When `ran` is `TRUE` also `expr_work`
#'   and `qc_expr` (corrected matrices) and `record` (a one-row data.frame:
#'   `method`, `batch_cols` (`;`-separated), `group_col`, `mean_only`, `seed`,
#'   `input_matrix`, `qc_matrix`, `n_samples`, `n_features`,
#'   `n_features_unadjusted`, `n_features_unadjusted_qc`, `output_file`).
compute_rna_batch_correction <- function(pre, config) {
    cfg <- config$modes$rna
    bc <- get_rna_batch_correction_config(cfg)
    if (!bc$enabled) return(list(ran = FALSE))

    assert_pre_contract(pre, stage = "rna")
    expr_work <- as.matrix(pre$expr_work)
    qc_transform <- resolve_rna_qc_transform(pre$info, cfg$qc)

    group_col <- resolve_rna_de_factor_col(cfg, pre$contrasts)
    design <- prepare_rna_batch_correction(
        meta       = pre$meta,
        sample_ids = colnames(expr_work),
        sample_col = cfg$id_columns$sample_col %||% "SampleID",
        group_col  = group_col,
        batch_cols = bc$batch_cols
    )
    seed <- as.integer(config$params$seed %||% 42)

    work <- correct_rna_batch_matrix(expr_work, design, bc, seed = seed)
    # Only TMMlogCPM has a QC matrix of its own; otherwise it is expr_work.
    qc <- if (isTRUE(qc_transform$applied)) {
        correct_rna_batch_matrix(rna_qc_expr(pre, cfg$qc), design, bc, seed = seed)
    } else {
        work
    }

    input_matrix <- pre$info$expr_work_method %||% "norm"
    message(sprintf(
        "[batch_correction] %s on %s, batch column(s) %s, keeping '%s'. DE input is unchanged.",
        bc$method, input_matrix, paste(bc$batch_cols, collapse = " + "), group_col))
    if (work$n_unadjusted > 0 || qc$n_unadjusted > 0) {
        warning(sprintf(paste0(
            "[batch_correction] ComBat left %d feature(s) of the analysis matrix and %d of ",
            "the QC matrix as they were: they are constant within a batch."),
            work$n_unadjusted, qc$n_unadjusted), call. = FALSE)
    }

    list(
        ran       = TRUE,
        expr_work = work$expr,
        qc_expr   = qc$expr,
        record    = data.frame(
            method                   = bc$method,
            batch_cols               = paste(bc$batch_cols, collapse = ";"),
            group_col                = group_col,
            mean_only                = bc$mean_only,
            seed                     = if (identical(bc$method, "combat")) seed else NA_integer_,
            input_matrix             = input_matrix,
            qc_matrix                = qc_transform$label,
            n_samples                = ncol(expr_work),
            n_features               = nrow(expr_work),
            n_features_unadjusted    = work$n_unadjusted,
            n_features_unadjusted_qc = qc$n_unadjusted,
            output_file              = rna_batch_corr_file_name(bc$method, input_matrix),
            stringsAsFactors = FALSE
        )
    )
}

#' Write the batch-corrected matrix and the record of the correction
#'
#' Writes the corrected analysis matrix to
#' `rna_batch_corr_<method>_<input>.tsv` and the record to
#' `batch_correction.tsv`. Files of both kinds left in a reused output
#' directory by an earlier run are removed first, also when no correction ran:
#' readers take the record's existence as evidence of a correction.
#'
#' @param bc_res Output of [compute_rna_batch_correction()].
#' @param datasets_dir The run's `Datasets` directory.
#' @return The written paths, or `character(0)` when no correction ran.
write_rna_batch_correction <- function(bc_res, datasets_dir) {
    stale <- c(
        file.path(datasets_dir, "batch_correction.tsv"),
        list.files(datasets_dir, pattern = "^rna_batch_corr_.*\\.tsv$", full.names = TRUE)
    )
    stale <- stale[file.exists(stale)]
    if (length(stale) > 0) file.remove(stale)

    if (!isTRUE(bc_res$ran)) return(character(0))

    corrected_df <- as.data.frame(bc_res$expr_work, check.names = FALSE)
    corrected_df <- cbind(gene = rownames(corrected_df), corrected_df, stringsAsFactors = FALSE)
    rownames(corrected_df) <- NULL

    c(
        save_tsv(corrected_df, datasets_dir, bc_res$record$output_file),
        save_tsv(bc_res$record, datasets_dir, "batch_correction.tsv")
    )
}

#' Methods sentence for the batch correction, from the record of what ran
#'
#' @param record_tbl The `batch_correction.tsv` table
#'   ([compute_rna_batch_correction()]), or `NULL` when no correction ran.
#' @return A paragraph, or `NULL` when there is no record.
#' @examples
#' rna_batch_correction_sentence(data.frame(
#'     method = "limma", batch_cols = "Batch", group_col = "Condition",
#'     mean_only = NA, input_matrix = "VST", n_features_unadjusted = 0,
#'     n_features_unadjusted_qc = 0))
rna_batch_correction_sentence <- function(record_tbl) {
    if (is.null(record_tbl) || nrow(record_tbl) == 0) return(NULL)
    rec <- record_tbl[1, , drop = FALSE]
    cols <- strsplit(as.character(rec$batch_cols), ";", fixed = TRUE)[[1]]

    input <- switch(as.character(rec$input_matrix),
        VST          = "variance-stabilized (VST) expression matrix",
        TMMlogCPM    = "TMM log-CPM expression matrix",
        preprocessed = "preprocessed expression matrix, as supplied",
        "normalized expression matrix"
    )
    backend <- if (identical(as.character(rec$method), "combat")) {
        paste0("sva::ComBat (parametric empirical Bayes",
               if (isTRUE(as.logical(rec$mean_only))) ", batch means only" else "", ")")
    } else {
        "limma::removeBatchEffect"
    }

    sentence <- sprintf(paste0(
        "For exploratory display, a batch-corrected copy of the %s was computed with %s ",
        "(batch column%s %s), with the group column %s kept in the design."),
        input, backend, if (length(cols) > 1) "s" else "",
        paste(cols, collapse = " and "), rec$group_col)

    n_unadjusted <- max(as.numeric(rec$n_features_unadjusted),
                        as.numeric(rec$n_features_unadjusted_qc), na.rm = TRUE)
    if (is.finite(n_unadjusted) && n_unadjusted > 0) {
        sentence <- paste(sentence, sprintf(paste0(
            "Up to %s gene(s) with a constant value within a batch were left uncorrected."),
            format(n_unadjusted, big.mark = ",")))
    }

    paste(sentence, paste0(
        "The corrected values were used only for the outputs labelled batch-corrected; ",
        "differential expression and every other result used the uncorrected data. ",
        "Batch and group can be separated only as far as the experimental design allows, ",
        "and the corrected values do not show that the remaining variation is biological."))
}

#' How the clustering group column relates to what the correction kept
#'
#' The correction keeps one column in its design, the group the DE compares.
#' Partition clustering and binary patterns then summarize the corrected values
#' by `clustering$group_col`, which may be another column. Whether its signal
#' survived the correction depends on how the two relate:
#'
#'   * `"same"`: the clustering column is the kept column.
#'   * `"determined_by_kept_group"`: every kept group falls in one clustering
#'     group, so the clustering column's effect lies inside the kept design.
#'   * `"batch_column"`: the clustering column is a column the correction
#'     removed.
#'   * `"confounded_with_batch"`: within the kept groups, the clustering column
#'     cannot be told apart from batch (the design of kept group, clustering
#'     group and batch is not identifiable), so the correction may have removed
#'     it.
#'   * `"not_protected"`: the clustering column is identifiable next to batch,
#'     but it was not in the correction's model. Identifiable is not the same
#'     as preserved: its signal was not explicitly protected and may be altered
#'     by the correction, also when it is balanced across batches.
#'
#' The first two and the last allow the corrected clustering to run; the last
#' one with a warning.
#'
#' @param kept Kept-group labels per sample.
#' @param clustering Clustering-group labels per sample.
#' @param batch List of one or two vectors of batch labels per sample.
#' @param kept_col,clustering_col,batch_cols The column names behind them.
#' @return A list with `relation` (one of the values above) and `proceed`
#'   (logical).
#' @examples
#' classify_rna_clustering_group(
#'     kept = c("A", "A", "B", "B"), clustering = c("x", "y", "x", "y"),
#'     batch = list(c("b1", "b2", "b2", "b1")),
#'     kept_col = "Condition", clustering_col = "Line", batch_cols = "Batch")
classify_rna_clustering_group <- function(kept, clustering, batch, kept_col,
                                          clustering_col, batch_cols) {
    result <- function(relation) {
        list(relation = relation,
             proceed = relation %in% c("same", "determined_by_kept_group", "not_protected"))
    }
    if (clustering_col %in% batch_cols) return(result("batch_column"))
    if (identical(clustering_col, kept_col)) return(result("same"))

    as_levels <- function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "NA"
        factor(x)
    }
    kept <- as_levels(kept)
    clustering <- as_levels(clustering)
    if (all(rowSums(table(kept, clustering) > 0) == 1)) {
        return(result("determined_by_kept_group"))
    }

    dummies <- function(f) {
        if (nlevels(f) < 2) return(NULL)
        stats::model.matrix(~f)[, -1, drop = FALSE]
    }
    keep_design <- cbind(stats::model.matrix(~kept), dummies(clustering))
    batch_design <- do.call(cbind, lapply(batch, function(b) dummies(as_levels(b))))
    # The correction already required kept group + batch to be identifiable.
    # Adding the clustering group must leave every batch column identifiable.
    n_batch <- if (is.null(batch_design)) 0L else ncol(batch_design)
    identifiable <- qr(cbind(keep_design, batch_design))$rank ==
        qr(keep_design)$rank + n_batch
    result(if (identifiable) "not_protected" else "confounded_with_batch")
}

#' Record of the clustering of the batch-corrected matrix
#'
#' @param status `"completed"` or `"skipped"`.
#' @param reason Short code for a skip (`"batch_column"`,
#'   `"confounded_with_batch"`, `"no_output"`), or `NA` when completed.
#' @param method Correction method that ran.
#' @param kept_col The group column the correction kept.
#' @param clustering_col `clustering$group_col`.
#' @param relation Output `relation` of [classify_rna_clustering_group()].
#' @param batch_cols Batch column(s) that were corrected for.
#' @return A one-row data.frame: `status`, `reason`, `method`,
#'   `kept_group_col`, `clustering_group_col`, `relation`, `batch_cols`
#'   (`;`-separated), `output_dir`.
rna_batch_clustering_record <- function(status, reason, method, kept_col, clustering_col,
                                        relation, batch_cols) {
    data.frame(
        status               = status,
        reason               = reason,
        method               = method,
        kept_group_col       = kept_col,
        clustering_group_col = clustering_col,
        relation             = relation,
        batch_cols           = paste(batch_cols, collapse = ";"),
        output_dir           = file.path(sprintf("Batch_corrected_%s", method), "Clustering"),
        stringsAsFactors = FALSE
    )
}

#' Reader-facing note on the clustering of the batch-corrected matrix
#'
#' Says why it was skipped, or, when it ran on a clustering group that the
#' correction did not protect, what that means. Used for the run log and the
#' report, from the record of what happened.
#'
#' @param record The `clustering_batch_correction.tsv` table
#'   ([rna_batch_clustering_record()]), or `NULL`.
#' @return A sentence or two, or `NULL` when there is nothing to say.
#' @examples
#' rna_batch_clustering_note(rna_batch_clustering_record(
#'     "skipped", "batch_column", "limma", "Condition", "Batch", "batch_column", "Batch"))
rna_batch_clustering_note <- function(record) {
    if (is.null(record) || nrow(record) == 0) return(NULL)
    rec <- record[1, , drop = FALSE]
    clustering_col <- rec$clustering_group_col
    kept_col <- rec$kept_group_col

    if (identical(as.character(rec$status), "skipped")) {
        why <- switch(as.character(rec$reason),
            batch_column = sprintf(paste0(
                "the clustering group column '%s' is a batch column that the correction ",
                "removed, so clustering the corrected values by it would not be meaningful"),
                clustering_col),
            confounded_with_batch = sprintf(paste0(
                "within the groups of '%s' (the column the correction kept), the clustering ",
                "group column '%s' cannot be told apart from batch, so the correction may ",
                "have removed its signal"), kept_col, clustering_col),
            no_output = paste0(
                "the clustering produced no output for the corrected values (for example, ",
                "too few differentially expressed genes)"),
            "of a condition recorded as unknown"
        )
        return(sprintf(paste0(
            "The clustering of the batch-corrected values was skipped: %s. The clustering ",
            "of the uncorrected values is not affected."), why))
    }

    if (identical(as.character(rec$relation), "not_protected")) {
        note <- sprintf(paste0(
            "The clustering groups genes by '%s', which is not the column the correction ",
            "kept ('%s') and is not determined by it. It can be told apart from batch in ",
            "this design, but it was not part of the correction's model: its signal was not ",
            "explicitly protected and may be altered by the correction, also where it is ",
            "balanced across batches."), clustering_col, kept_col)
        if (identical(as.character(rec$method), "combat")) {
            note <- paste(note, paste0(
                "ComBat also adjusts each batch's spread, which can reduce that signal ",
                "further."))
        }
        return(paste(note, "Read the batch-corrected clusters against the uncorrected ones."))
    }
    NULL
}

#' Partition-clustering directory that holds this run's heatmap
#'
#' The partition step writes into `Partition_clustering_<k>_clusters`. A reused
#' output directory can hold the directory of an earlier run with another `k`;
#' this returns the one that contains the heatmap, not simply the first.
#'
#' @param clustering_root A `Clustering` directory.
#' @return The directory path, or `NA_character_` when none holds a heatmap.
rna_partition_heatmap_dir <- function(clustering_root) {
    dirs <- list.dirs(file.path(clustering_root, "Partition_clustering"), recursive = FALSE)
    dirs <- dirs[file.exists(file.path(dirs, "Partition_clustering_heatmap.png"))]
    if (length(dirs) == 0) NA_character_ else sort(dirs)[1]
}

#' What the report shows for the clustering of the batch-corrected matrix
#'
#' Built from the record and the files on disk, independently of whether the
#' uncorrected clustering produced figures.
#'
#' @param record The `clustering_batch_correction.tsv` table, or `NULL` when
#'   no batch-corrected clustering was attempted.
#' @param rna_dir The RNA output directory.
#' @return A list with `show` (logical: the report has something to say),
#'   `status`, `method`, `note` ([rna_batch_clustering_note()]),
#'   `partition_heatmap` (path or `NA`) and `binary_heatmaps` (paths).
rna_batch_clustering_view <- function(record, rna_dir) {
    empty <- list(show = FALSE, status = NA_character_, method = NA_character_, note = NULL,
                  partition_heatmap = NA_character_, binary_heatmaps = character(0))
    if (is.null(record) || nrow(record) == 0) return(empty)

    status <- as.character(record$status[1])
    view <- empty
    view$status <- status
    view$method <- as.character(record$method[1])
    view$note <- rna_batch_clustering_note(record)
    if (identical(status, "completed")) {
        root <- file.path(rna_dir, as.character(record$output_dir[1]))
        part_dir <- rna_partition_heatmap_dir(root)
        if (!is.na(part_dir)) {
            view$partition_heatmap <- file.path(part_dir, "Partition_clustering_heatmap.png")
        }
        view$binary_heatmaps <- list.files(file.path(root, "Binary_patterns"),
                                           pattern = "^Heatmap_[01]+\\.png$", full.names = TRUE)
    }
    view$show <- identical(status, "skipped") || !is.na(view$partition_heatmap) ||
        length(view$binary_heatmaps) > 0
    view
}
