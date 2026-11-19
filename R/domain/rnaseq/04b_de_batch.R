#' Batch in the RNA-seq DE model
#'
#' Optional `modes.rna.batch_in_model` block: one or two metadata columns enter
#' the DE design as additive categorical covariates, so each contrast is
#' estimated after adjusting for batch.
#'
#'   * DESeq2: `~ Batch + <group>` (group last, see [run_deseq2_de()]).
#'   * limma (preprocessed input): `~ 0 + <group> + Batch`, treatment-coded.
#'
#' This changes the statistics only. It never modifies an expression matrix;
#' QC, normalization and clustering inputs are untouched. The contrasts and
#' their directions are the same as without batch.
#'
#' The fitted design is recorded by [rna_de_design_record()], written to
#' `Datasets/de_design.tsv`, and the report's Methods sentence is built from
#' that file by [rna_de_design_sentence()], so the text describes the model
#' that ran rather than the config that asked for it.

#' Resolve the batch-in-model config
#'
#' The one reading of `modes.rna.batch_in_model`, used by the DE, the config
#' validation and the exports.
#'
#' @param rna_cfg The `config$modes$rna` list.
#' @return A list with `enabled` (logical) and `batch_cols` (character; empty
#'   when disabled).
#' @examples
#' get_rna_batch_model_config(list())  # disabled
#' get_rna_batch_model_config(list(batch_in_model = list(enabled = TRUE,
#'                                                       batch_cols = "Batch")))
get_rna_batch_model_config <- function(rna_cfg) {
    bm <- rna_cfg$batch_in_model %||% list()
    enabled <- isTRUE(bm$enabled)
    batch_cols <- as.character(unlist(bm$batch_cols %||% character(0)))
    if (enabled) {
        if (!length(batch_cols) %in% 1:2 || !all(nzchar(batch_cols)) ||
            anyDuplicated(batch_cols)) {
            stop("[config] batch_in_model$batch_cols must name one or two distinct metadata ",
                 "columns when batch_in_model$enabled is true (got: ",
                 if (length(batch_cols)) paste(batch_cols, collapse = ", ") else "none", ").",
                 call. = FALSE)
        }
    }
    list(enabled = enabled, batch_cols = if (enabled) batch_cols else character(0))
}

#' Check the batch columns and the combined design, and make batch categorical
#'
#' Run on the samples that enter DE, before any model is built. Stops on a
#' missing, non-syntactic or NA-holding batch column, a batch column that is
#' the group column, a batch with fewer than two levels among these samples, a
#' combined design (group plus every batch column) that is not of full rank,
#' and a design that leaves no residual degrees of freedom.
#'
#' Batch columns are converted with `factor(as.character(x))`, so numeric
#' labels (e.g. pair numbers 1, 2, 3) are categories, not a continuous slope.
#'
#' @param meta Sample metadata (may hold more samples than enter DE).
#' @param sample_ids Character vector of the samples that enter DE.
#' @param sample_col Metadata column holding sample IDs.
#' @param factor_col The group column the contrasts compare.
#' @param batch_cols One or two batch column names.
#' @return `meta` with each batch column converted to a factor.
prepare_rna_batch_design <- function(meta, sample_ids, sample_col, factor_col, batch_cols) {
    missing_cols <- setdiff(batch_cols, colnames(meta))
    if (length(missing_cols) > 0) {
        stop("[batch_in_model] Batch column(s) not found in metadata: ",
             paste(missing_cols, collapse = ", "), ". Available columns: ",
             paste(colnames(meta), collapse = ", "), ".", call. = FALSE)
    }
    if (factor_col %in% batch_cols) {
        stop("[batch_in_model] Batch column '", factor_col, "' is the group column the ",
             "contrasts compare; batch and group cannot be the same variable.", call. = FALSE)
    }
    bad_names <- batch_cols[make.names(batch_cols) != batch_cols]
    if (length(bad_names) > 0) {
        stop("[batch_in_model] Batch column name(s) ", paste(bad_names, collapse = ", "),
             " cannot be used in a model formula. Rename the column(s) in the metadata ",
             "using letters, digits, '.' and '_' only.", call. = FALSE)
    }

    rows <- match(sample_ids, as.character(meta[[sample_col]]))
    rows <- rows[!is.na(rows)]

    for (b in batch_cols) {
        vals <- as.character(meta[[b]])
        bad <- rows[is.na(vals[rows]) | !nzchar(vals[rows])]
        if (length(bad) > 0) {
            stop("[batch_in_model] Batch column '", b, "' has no value for ", length(bad),
                 " sample(s): ", paste(head(as.character(meta[[sample_col]][bad]), 10),
                                       collapse = ", "),
                 ". Assign a batch to every sample that enters DE.", call. = FALSE)
        }
        meta[[b]] <- factor(vals)
        if (nlevels(droplevels(meta[[b]][rows])) < 2) {
            stop("[batch_in_model] Batch column '", b, "' has a single level among the ",
                 "samples that enter DE, so there is no batch effect to model. Turn ",
                 "batch_in_model off or choose another column.", call. = FALSE)
        }
    }

    sub <- droplevels(as.data.frame(meta[rows, , drop = FALSE]))
    assert_batch_not_confounded(
        batch_var = sub[[batch_cols[1]]],
        group_var = sub[[factor_col]],
        batch2    = if (length(batch_cols) == 2) sub[[batch_cols[2]]]
    )

    x <- stats::model.matrix(rna_batch_design_formula(factor_col, batch_cols), data = sub)
    if (nrow(x) <= ncol(x)) {
        stop("[batch_in_model] The design with batch has ", ncol(x), " coefficients for ",
             nrow(x), " samples, leaving no residual degrees of freedom to estimate ",
             "variability. Model fewer batch columns or turn batch_in_model off.",
             call. = FALSE)
    }
    meta
}

#' The DESeq2 design formula with batch first and the group last
#'
#' The group goes last so that `DESeq2::results(dds)` without a contrast still
#' reports a group comparison: in default mode it takes the last coefficient,
#' and in legacy (expanded) mode the last variable of the formula.
#'
#' @param factor_col The group column.
#' @param batch_cols Batch column names, possibly empty.
#' @return A formula.
rna_batch_design_formula <- function(factor_col, batch_cols = character(0)) {
    stats::as.formula(paste("~", paste(c(batch_cols, factor_col), collapse = " + ")))
}

#' Build the limma design for preprocessed RNA data with batch covariates
#'
#' One column per group level (no intercept), so the existing contrasts can
#' name the levels directly, followed by treatment-coded columns for each batch
#' factor, named `<batch_col>_<level>`.
#'
#' @param group Factor of group labels, levels already made syntactic.
#' @param meta Metadata aligned to the expression columns, batch columns
#'   already factors ([prepare_rna_batch_design()]).
#' @param batch_cols One or two batch column names.
#' @return A numeric design matrix, one row per sample.
build_rna_limma_design <- function(group, meta, batch_cols) {
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    for (b in batch_cols) {
        bv <- droplevels(meta[[b]])
        bm <- stats::model.matrix(~ bv)[, -1, drop = FALSE]
        colnames(bm) <- paste0(b, "_", levels(bv)[-1])
        design <- cbind(design, bm)
    }
    # unique = TRUE only renames later (batch) columns on a clash, so the group
    # level names the contrasts refer to stay as they are.
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    design
}

#' Record the DE design that was actually fitted
#'
#' @param engine `"DESeq2"` or `"limma"`.
#' @param mode DESeq2 mode (`"default"` / `"legacy"`), or `NA` for limma.
#' @param formula Character form of the fitted design.
#' @param batch_cols Batch columns in the model (empty when none).
#' @param design_matrix The fitted (non-expanded) design matrix.
#' @return A one-row data.frame: `engine`, `mode`, `formula`, `batch_cols`
#'   (`;`-separated, empty when none), `n_samples`, `n_coefficients`,
#'   `residual_df`.
rna_de_design_record <- function(engine, mode, formula, batch_cols, design_matrix) {
    data.frame(
        engine         = engine,
        mode           = mode,
        formula        = formula,
        batch_cols     = paste(batch_cols, collapse = ";"),
        n_samples      = nrow(design_matrix),
        n_coefficients = ncol(design_matrix),
        residual_df    = nrow(design_matrix) - qr(design_matrix)$rank,
        stringsAsFactors = FALSE
    )
}

#' Methods sentence for batch in the DE model, from the recorded design
#'
#' @param design_tbl The `de_design.tsv` table ([rna_de_design_record()]), or
#'   `NULL` when no DE model was fitted (e.g. precomputed DE).
#' @return A sentence, or `NULL` when the fitted model had no batch term, so
#'   the report keeps its existing wording.
#' @examples
#' rna_de_design_sentence(data.frame(engine = "DESeq2", mode = "default",
#'                                   formula = "~ Batch + Condition",
#'                                   batch_cols = "Batch"))
rna_de_design_sentence <- function(design_tbl) {
    if (is.null(design_tbl) || nrow(design_tbl) == 0) return(NULL)
    cols <- as.character(design_tbl$batch_cols[1])
    if (is.na(cols) || !nzchar(cols)) return(NULL)
    cols <- strsplit(cols, ";", fixed = TRUE)[[1]]
    sentence <- sprintf(paste0(
        "The %s model included batch as a covariate (design `%s`; batch column%s %s), ",
        "so each contrast was estimated after adjusting for batch."),
        design_tbl$engine[1], design_tbl$formula[1],
        if (length(cols) > 1) "s" else "", paste(cols, collapse = " and "))
    if (identical(design_tbl$engine[1], "DESeq2") && identical(design_tbl$mode[1], "legacy")) {
        sentence <- paste(sentence, paste0(
            "In legacy mode (betaPrior), the batch coefficients were shrunk by the same ",
            "kind of prior as the group coefficients."))
    }
    sentence
}
