
#' Build Final Results Table for Metabolomics
#'
#' Mirrors the proteomics structure via build_final_results_generic():
#'   feature_id | annotations | expression | DE stats | pass_any_contrast
#'
#' Required contract columns: feature_id, linearFC.<cn>, pvalue.<cn>,
#' padj.<cn>, upDown.<cn>, manual_cutoffs.<cn>, pass_any_contrast.
#'
#' Metabolomics-specific extras (from row_data): compound name, RT, m/z,
#' HMDB ID, molecular formula, etc.
#'
#' @param pre Preprocessing results (requires expr_work, row_data).
#' @param summary_df DE summary with aligned column naming (linearFC., padj., pass.).
#' @param contrast_labels Character vector of contrast labels.
#' @param row_data Optional feature annotation data.frame. Defaults to pre$row_data.
#' @param feature_id_col Name of ID column (default "feature_id").
#' @return A data.frame matching the proteomics final_results structure.
build_final_results_metabolomics <- function(
    pre,
    summary_df,
    contrast_labels,
    row_data       = NULL,
    feature_id_col = "feature_id",
    cv_contrasts_df = NULL,
    config          = NULL
) {
    # Wrap contrast labels into contrasts_df for the generic API
    contrasts_df <- data.frame(
        Contrast_name = contrast_labels,
        stringsAsFactors = FALSE
    )

    # Build annotation column mapping from row_data
    rd <- row_data %||% pre$row_data
    annot_cols <- NULL
    if (!is.null(rd)) {
        available <- setdiff(colnames(rd), feature_id_col)
        if (length(available) > 0) {
            annot_cols <- setNames(available, available)
        }
    }

    # Per-group CV requires the structured contrasts (Factor/Numerator/Denominator),
    # which contrast_labels alone do not carry.
    cv_cols <- build_group_cv_metabolomics(pre, cv_contrasts_df, config)

    build_final_results_generic(
        summary_df     = summary_df,
        expr_df        = pre$expr_work,
        contrasts_df   = contrasts_df,
        feature_id_col = feature_id_col,
        annot_cols     = annot_cols,
        row_data       = rd,
        fc_is_signed   = TRUE,
        mode           = "metabolomics",
        cv_cols        = cv_cols
    )
}

#' Build per-group CV columns for metabolomics final results
#'
#' CV is computed on the post-normalization LINEAR matrix, reconstructed from
#' \code{expr_work} as \code{2^expr_work - pseudocount}. This inversion is only
#' valid when the workspace is genuinely log2-based, so it is gated on a
#' normalization whitelist (see below). This keeps metabolomics CV on normalized
#' data, comparable to the RNA-seq and proteomics CV columns.
#'
#' Guards (return \code{NULL}, never approximate an inverse):
#' \itemize{
#'   \item If feature scaling is enabled (\code{normalization$scaling != "none"}),
#'     the back-transform is invalid (centering breaks positivity), so we fall
#'     back to the filtered PRE-normalization linear matrix (\code{expr_filt});
#'     the CV then includes technical variation that normalization would have
#'     removed (see the contract doc / column note).
#'   \item If the workspace is not provably log2 — i.e. \code{chosen_norm ==
#'     "median"} with a non-log2 \code{normalization$transform} (log10/glog10/
#'     none), or an unrecognized \code{chosen_norm} — we skip CV with a warning
#'     rather than \code{2^}-inverting a non-log2 matrix into wrong numbers.
#' }
#'
#' @param pre Metabolomics preprocessing results (uses \code{expr_work},
#'   \code{expr_filt}, \code{meta}).
#' @param contrasts_df Structured contrasts (Factor, Numerator, Denominator).
#' @param config Full pipeline config (feature flag, sample-ID column, scale).
#' @return Feature-indexed data.frame of \code{CV.<group>} columns, or NULL.
build_group_cv_metabolomics <- function(pre, contrasts_df, config = NULL) {
    if (is.null(config) || is.null(contrasts_df)) return(NULL)
    enabled <- config$modes$metabolomics$excel$group_cv %||% TRUE
    if (!isTRUE(enabled)) return(NULL)
    if (is.null(pre$expr_work) || is.null(pre$meta)) return(NULL)

    cfg_mode      <- config$modes$metabolomics %||% list()
    sample_id_col <- cfg_mode$effects$samples %||% "sample_id"
    norm_cfg      <- cfg_mode$preprocessing %||% list()
    scaling       <- tolower(norm_cfg$scaling %||% "none")
    pseudocount   <- norm_cfg$pseudocount %||% 1

    if (identical(scaling, "none")) {
        # The 2^x inverse is only valid on a genuinely log2 workspace. These
        # normalization methods hardcode log2(x + pseudocount) regardless of
        # config, so they are always invertible.
        # MAINTAINER NOTE: any NEW normalization method that always log2-
        # transforms MUST be added to this whitelist; otherwise CV will safely
        # (but silently) skip for it via the guard below.
        always_log2_norms <- c("tss", "pqn", "eigenms", "eigenms_forced")
        chosen_norm <- tolower(cfg_mode$preprocessing$chosen_norm %||% "")
        transform   <- tolower(norm_cfg$transform %||% "log2")
        # The median path honors the configurable transform, so it is only log2
        # when transform == "log2". Everything outside the whitelist (incl.
        # unknown chosen_norm) is not provably log2 -> skip rather than guess.
        provably_log2 <- chosen_norm %in% always_log2_norms ||
            (identical(chosen_norm, "median") && identical(transform, "log2"))
        if (!isTRUE(provably_log2)) {
            warning(sprintf(
                "metabolomics group CV: chosen_norm='%s' with transform='%s' is not a provably-log2 workspace; skipping CV (a 2^x back-transform would produce wrong values on a non-log2 matrix).",
                chosen_norm, transform
            ))
            return(NULL)
        }
        # Post-normalization linear matrix (exact inverse on the log2 workspace).
        expr_linear <- 2^as.matrix(pre$expr_work) - pseudocount
        # Floor tiny negatives from approximate back-transform / drift (matches
        # compute_linear_rsd()); keeps the value count for CV.
        neg <- is.finite(expr_linear) & expr_linear < 0
        expr_linear[neg] <- .Machine$double.eps
    } else {
        # Scaling enabled -> reconstruction invalid; use filtered pre-norm linear.
        warning(sprintf(
            "metabolomics group CV: feature scaling = '%s'; using pre-normalization 'expr_filt' (CV includes technical variation, not directly comparable to other omics).",
            scaling
        ))
        if (is.null(pre$expr_filt)) return(NULL)
        expr_linear <- as.matrix(pre$expr_filt)
        # Align to the columns/features actually present in expr_work
        common_cols <- intersect(colnames(pre$expr_work), colnames(expr_linear))
        common_rows <- intersect(rownames(pre$expr_work), rownames(expr_linear))
        expr_linear <- expr_linear[common_rows, common_cols, drop = FALSE]
    }

    compute_group_cv_columns(
        expr_linear   = expr_linear,
        sample_meta   = pre$meta,
        sample_id_col = sample_id_col,
        contrasts_df  = contrasts_df
    )
}


#' Write standardized metabolomics outputs (Stage 1)
#'
#' Saves the internal representation as TSV files for downstream use.
#' @return Character vector of written file paths.
write_metabolomics_outputs <- function(pre, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  out_ds <- dirs$datasets
  
  written <- character(0)
  
  # Expression matrix (normalized)
  expr_df <- as.data.frame(pre$expr_work)
  expr_df$feature_id <- rownames(pre$expr_work)
  expr_df <- expr_df[, c("feature_id", setdiff(colnames(expr_df), "feature_id")),
                     drop = FALSE]
  written <- c(written, save_tsv(expr_df, out_ds, "expr_normalized.tsv"))
  
  # Expression matrix (raw)
  raw_df <- as.data.frame(pre$expr_raw)
  raw_df$feature_id <- rownames(pre$expr_raw)
  raw_df <- raw_df[, c("feature_id", setdiff(colnames(raw_df), "feature_id")),
                   drop = FALSE]
  written <- c(written, save_tsv(raw_df, out_ds, "expr_raw.tsv"))
  
  # Metadata
  written <- c(written, save_tsv(pre$meta, out_ds, "metadata_aligned.tsv"))
  
  # Row data (feature annotations)
  written <- c(written, save_tsv(pre$row_data, out_ds, "feature_annotations.tsv"))
  
  # Sample map (if available, from CD raw parsing)
  if (!is.null(pre$sample_map)) {
    written <- c(written, save_tsv(pre$sample_map, out_ds, "sample_map_cd.tsv"))
  }

  written
}


#' Write metabolomics final results (TSV + Excel)
#'
#' Builds the final_results table from DE summary and exports:
#' - Final_results_metabolomics.tsv
#' - Final_results_ALL_P_*.xlsx (all features)
#' - Final_results_DE_P_*.xlsx (DE-only features with clustering order + z-scores)
#'
#' @param pre Preprocessing results.
#' @param de_res DE results (must contain summary_df and de_tables or contrasts).
#' @param config Full config.
#' @param out_dir Output directory.
#' @param clustering_res Optional clustering results for Excel row ordering.
#' @param inputs Optional loaded inputs; \code{inputs$contrasts} (Factor,
#'   Numerator, Denominator) is used for per-group CV column resolution.
#' @return Character vector of written file paths.
write_metabolomics_final_results <- function(pre, de_res, config, out_dir,
                                              clustering_res = NULL,
                                              inputs = NULL) {
    if (is.null(de_res) || is.null(de_res$summary_df)) {
        message("metabolomics final_results: skipped (no DE results)")
        return(character(0))
    }

    dirs <- create_legacy_output_dirs(out_dir)
    files <- character(0)

    # Extract contrast labels from DE results
    contrast_labels <- names(de_res$de_tables)
    if (is.null(contrast_labels) || length(contrast_labels) == 0) {
        stop("metabolomics final_results: no contrast labels found in de_res$de_tables")
    }

    # Build final_results table
    final_results <- build_final_results_metabolomics(
        pre             = pre,
        summary_df      = de_res$summary_df,
        contrast_labels = contrast_labels,
        row_data        = pre$row_data,
        feature_id_col  = "feature_id",
        cv_contrasts_df = inputs$contrasts,
        config          = config
    )

    # Write TSV
    files <- c(files, save_tsv(final_results, dirs$datasets,
                                "Final_results_metabolomics.tsv"))

    # Extract Excel config
    excel_cfg <- config$modes$metabolomics$excel %||% list()

    # Resolve sample ID column
    sample_id_col <- config$modes$metabolomics$effects$samples %||% "sample_id"

    # Write Excel files (ALL + DE-only)
    excel_files <- tryCatch(
        write_final_results_excels_legacy_generic(
            final_results   = final_results,
            config          = config,
            out_dir         = out_dir,
            mode            = "metabolomics",
            id_col          = "feature_id",
            expr_for_de     = pre$expr_work,
            with_cutoffs    = TRUE,
            clustering_res  = clustering_res,
            sample_meta     = pre$meta,
            sample_id_col   = sample_id_col,
            annotation_rows = excel_cfg$annotation_rows,
            sample_label_cols = excel_cfg$sample_label_cols
        ),
        error = function(e) {
            warning("Excel export failed: ", e$message)
            character(0)
        }
    )
    files <- c(files, excel_files)

    message("metabolomics final_results: wrote ", length(files), " files (",
            nrow(final_results), " features)")
    unique(files)
}