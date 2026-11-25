# 03_mod_rnaseq_de
#' Run (or load) RNA-seq differential expression
#'
#' @param pre Preprocessed RNA object (\code{rna_pre}).
#' @param inputs Loaded RNA inputs (\code{rna_inputs}); \code{inputs$contrasts}
#'   is the supplied contrasts table, or NULL.
#' @param config Validated config.
#' @param verbose Unused; kept for the module signature.
#' @return The DE result of the route that ran (\code{tables}, \code{design},
#'   and \code{dds} for DESeq2), plus \code{contrasts}: the contrast table the
#'   DE step used -- the supplied one, or the auto-generated one when none was
#'   supplied. On the pre-computed route it is set by
#'   \code{load_precomputed_rna_de()}, and is NULL when the table names came
#'   from the file names rather than from a contrasts table. Match it to
#'   \code{tables} by \code{Contrast_name}: a contrast DESeq2 skipped (a level
#'   absent from the metadata) has a row here but no table.
mod_rnaseq_de <- function(pre, inputs, config, verbose = FALSE) {
    assert_pre_contract(pre, stage = "rna")

    # Check for pre-computed DE tables (skip limma/DESeq2 if provided)
    de_table_files <- config$modes$rna$files$de_table
    has_precomputed <- !is.null(de_table_files) && length(de_table_files) > 0 &&
                       all(nzchar(unlist(de_table_files)))

    if (has_precomputed) {
        message("RNA DE: loading pre-computed DE tables (skipping limma/DESeq2)")
        # The loader sets $contrasts itself: it alone knows whether the
        # supplied names labelled the tables or it fell back to file names.
        return(load_precomputed_rna_de(config, contrasts_df = inputs$contrasts))
    }

    # Extract DE configuration
    de_cfg <- config$modes$rna$de
    batch_cfg <- get_rna_batch_model_config(config$modes$rna)

    # Get contrasts (already loaded by load_rna_inputs), auto-generate if missing
    contrasts_df <- inputs$contrasts
    if (is.null(contrasts_df)) {
        group_col <- config$modes$rna$effects$color %||%
                     config$modes$rna$filtering$group_col %||% "Condition"
        contrasts_df <- auto_generate_contrasts(pre$meta, group_col)
    }

    # Dispatch DE method based on source type
    source_type <- attr(pre, "source_type") %||%
                   attr(pre$de_input, "source_type") %||%
                   "matrix"

    de_res <- if (source_type == "preprocessed") {
        run_limma_rna_de(
            expr = pre$expr_work,
            meta = pre$meta,
            contrasts_df = contrasts_df,
            de_cfg = de_cfg,
            batch_cfg = batch_cfg
        )
    } else {
        run_deseq2_de(
            counts = pre$de_input,
            meta = pre$meta,
            contrasts_df = contrasts_df,
            de_cfg = de_cfg,
            batch_cfg = batch_cfg
        )
    }
    # Keep the contrast table this step resolved. Auto-generated contrasts
    # exist nowhere else (pre$contrasts is NULL then), and the table names
    # alone cannot say which level is the numerator; the omics bundle carries
    # this so integration need not resolve the contrasts a second time.
    de_res$contrasts <- contrasts_df
    de_res
}

#' Auto-generate all pairwise contrasts from metadata
#' @param meta    data.frame with sample metadata
#' @param group_col character, column name defining biological groups
#' @return data.frame with Contrast_name, Factor, Numerator, Denominator
auto_generate_contrasts <- function(meta, group_col) {
    lvls <- sort(unique(as.character(meta[[group_col]])))
    lvls <- lvls[!is.na(lvls) & nzchar(lvls)]
    if (length(lvls) < 2) stop("Cannot auto-generate contrasts: fewer than 2 groups in '", group_col, "'.")
    pairs <- combn(lvls, 2)
    df <- data.frame(
        Contrast_name = apply(pairs, 2, function(p) paste0(p[1], "_vs_", p[2])),
        Factor        = group_col,
        Numerator     = pairs[1, ],
        Denominator   = pairs[2, ],
        stringsAsFactors = FALSE
    )
    message(sprintf("Auto-generated %d pairwise contrast(s) from '%s': %s",
                    nrow(df), group_col, paste(df$Contrast_name, collapse = ", ")))
    df
}
