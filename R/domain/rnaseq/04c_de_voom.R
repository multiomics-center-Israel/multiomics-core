#' Run limma-voom differential expression on RNA-seq counts
#'
#' TMM-normalises the (already gene-filtered) count matrix with edgeR, converts
#' it to log2-CPM with precision weights via [limma::voom()], and fits a
#' moderated linear model with [limma::eBayes()].
#'
#' When `de$block_col` names a metadata column, `de$block_model` (see
#' [resolve_rna_block_model()]) says how the block enters. "random" (the
#' default here): samples sharing a level are treated as correlated, the
#' within-block correlation is estimated with [limma::duplicateCorrelation()]
#' and passed to [limma::lmFit()]; following the limma user's guide, voom and
#' duplicateCorrelation are each run twice. A non-finite consensus refits
#' without blocking, and `info` records that outcome. "fixed": the block is a
#' term of the design (`~ block + factor`), as in the DESeq2 path.
#'
#' Returns the same structure as [run_deseq2_de()] so downstream code works
#' unchanged. `baseMean` carries voom's AveExpr (mean log2-CPM), not DESeq2's
#' mean of normalised counts.
#'
#' @param counts Integer count matrix (genes x samples). tximport input is not
#'   supported here.
#' @param meta Data frame with sample metadata.
#' @param contrasts_df Data frame with columns: Contrast_name, Factor, Numerator, Denominator.
#' @param de_cfg DE configuration list (sample_col, sample_alignment, block_col).
#' @return List with 'dds' (NULL), 'tables' (named list of DE result data
#'   frames) and 'info' (method, design, factor_col, block_col, block_model, n_blocks,
#'   block_correlation, block_requested). Each contrast is fitted as one
#'   coefficient, with its denominator as the reference level.
run_voom_rna_de <- function(counts, meta, contrasts_df, de_cfg) {
    prep <- prepare_count_de_inputs(counts, meta, contrasts_df, de_cfg,
                                    method = "limma_voom", pkgs = c("limma", "edgeR"))
    counts <- prep$counts; meta <- prep$meta
    factor_col <- prep$factor_col; safe_levels <- prep$safe_levels; sample_col <- prep$sample_col

    design_formula <- stats::as.formula(paste0("~ 0 + ", factor_col))
    design <- stats::model.matrix(design_formula, data = meta)
    colnames(design) <- safe_levels

    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")

    block <- resolve_de_block(meta, list(de = de_cfg), sample_col)
    block_model <- if (is.null(block)) NULL else resolve_rna_block_model(de_cfg, method = "limma_voom")
    block_correlation <- NA_real_
    if (identical(block_model, "fixed")) {
        # Block as a fixed term: one coefficient per extra block level, no
        # correlation to estimate. Same refusal rules as the DESeq2 path.
        meta$.block <- block
        mm <- check_fixed_block_design(stats::as.formula(paste0("~ .block + ", factor_col)),
                                       meta, de_cfg$block_col, factor_col)
        v <- limma::voom(dge, mm)
        message(sprintf("[run_voom_rna_de] Fixed block '%s' (%d levels) in the design",
                        de_cfg$block_col, nlevels(block)))
    } else if (!is.null(block)) {
        v <- limma::voom(dge, design)
        corfit <- limma::duplicateCorrelation(v, design, block = block)
        v <- limma::voom(dge, design, block = block, correlation = corfit$consensus)
        corfit <- limma::duplicateCorrelation(v, design, block = block)
        block_correlation <- corfit$consensus
        if (is.finite(block_correlation)) {
            message(sprintf("[run_voom_rna_de] Blocking on '%s' (%d levels), consensus correlation %.3f",
                            de_cfg$block_col, nlevels(block), block_correlation))
        } else {
            warning("duplicateCorrelation returned a non-finite consensus for de$block_col '",
                    de_cfg$block_col, "'; fitting without blocking.")
            v <- limma::voom(dge, design)
            block <- NULL
            block_model <- NULL
        }
    } else {
        v <- limma::voom(dge, design)
    }

    # One fit per contrast, with the denominator as the reference level, so the
    # contrast is a single coefficient. With voom weights and a block
    # correlation the ~ 0 + factor coefficients are correlated gene by gene, and
    # contrasts.fit() is then only approximate. A reparametrisation spans the
    # same column space, so voom weights and the consensus correlation above
    # still apply unchanged.
    tables <- list()
    for (i in seq_len(nrow(contrasts_df))) {
        cn  <- contrasts_df$Contrast_name[i]
        num <- make.names(contrasts_df$Numerator[i])
        den <- make.names(contrasts_df$Denominator[i])
        if (!num %in% safe_levels || !den %in% safe_levels) {
            warning(sprintf("Contrast '%s': level '%s' or '%s' not in '%s'. Skipping.",
                            cn, contrasts_df$Numerator[i], contrasts_df$Denominator[i], factor_col))
            next
        }
        lev <- stats::relevel(meta[[factor_col]], ref = den)
        design_cn <- if (identical(block_model, "fixed")) {
            stats::model.matrix(~ block + lev)
        } else {
            stats::model.matrix(~ lev)
        }
        colnames(design_cn) <- sub("^lev", "", colnames(design_cn))
        fit <- if (identical(block_model, "random")) {
            limma::lmFit(v, design_cn, block = block, correlation = block_correlation)
        } else {
            limma::lmFit(v, design_cn)
        }
        fit <- limma::eBayes(fit)
        de <- limma::topTable(fit, coef = num, adjust.method = "BH",
                              sort.by = "none", number = Inf)
        tab <- data.frame(
            baseMean       = de$AveExpr,
            log2FoldChange = de$logFC,
            lfcSE          = (fit$stdev.unscaled[, num] * sqrt(fit$s2.post))[rownames(de)],
            stat           = de$t,
            pvalue         = de$P.Value,
            padj           = de$adj.P.Val,
            FeatureID      = rownames(de),
            stringsAsFactors = FALSE
        )
        rownames(tab) <- rownames(de)
        tables[[cn]] <- tab
    }

    list(
        dds    = NULL,
        tables = tables,
        info   = list(
            method            = "limma_voom",
            design            = if (identical(block_model, "fixed")) {
                                    paste0("~", de_cfg$block_col, " + ", factor_col)
                                } else paste0("~", factor_col),
            factor_col        = factor_col,
            block_col         = if (is.null(block)) NULL else de_cfg$block_col,
            block_model       = block_model,
            n_blocks          = if (is.null(block)) NA_integer_ else nlevels(block),
            block_correlation = block_correlation,
            block_requested   = de_cfg$block_col
        )
    )
}

#' How a configured block enters the RNA-seq DE model
#'
#' The one reading of \code{de$block_model}, shared by the validator and the
#' limma-voom fit. DESeq2 has no random effect, so its default is "fixed";
#' limma-voom defaults to "random" (duplicateCorrelation).
#'
#' @param de_cfg The \code{modes$rna$de} config list.
#' @param method The DE method that will run; defaults to \code{de$method}.
#' @return "fixed", "random", or NULL when no \code{block_col} is set.
resolve_rna_block_model <- function(de_cfg, method = de_cfg$method %||% "deseq2") {
    bc <- de_cfg$block_col
    if (is.null(bc) || !nzchar(bc)) return(NULL)
    default <- if (identical(tolower(method), "limma_voom")) "random" else "fixed"
    tolower(de_cfg$block_model %||% default)
}

#' Validate and align inputs for a count-based limma/edgeR DE fit
#'
#' Shared by [run_voom_rna_de()] and [run_edger_rna_de()]: both take a raw
#' count matrix, one contrast factor and the same sample alignment rules.
#'
#' @param counts Integer count matrix (genes x samples); tximport is refused.
#' @param meta Sample metadata.
#' @param contrasts_df Contrasts (Contrast_name, Factor, Numerator, Denominator).
#' @param de_cfg DE config list (sample_col, sample_alignment).
#' @param method Method name, used in messages.
#' @param pkgs Packages that must be installed.
#' @return A list with \code{counts} (matrix, columns in metadata order),
#'   \code{meta} (aligned, factor column with syntactic levels),
#'   \code{factor_col}, \code{safe_levels} and \code{sample_col}.
prepare_count_de_inputs <- function(counts, meta, contrasts_df, de_cfg, method, pkgs) {
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Package '", pkg, "' is required for de$method: ", method, ".", call. = FALSE)
        }
    }
    if (identical(detect_source_type(counts), "tximport")) {
        stop("de$method ", method, " takes a count matrix, but this run has tximport input.\n",
             "  Use de$method deseq2 for tximport, or supply a gene-level count matrix.",
             call. = FALSE)
    }
    validate_count_matrix(counts)
    counts <- as.matrix(counts)

    if (nrow(contrasts_df) == 0) stop("contrasts_df is empty")
    factor_col <- unique(contrasts_df$Factor)
    if (length(factor_col) != 1) stop("All contrasts must use the same factor column")
    factor_col <- factor_col[[1]]
    if (!factor_col %in% colnames(meta)) {
        stop(sprintf("Factor column '%s' not found in metadata", factor_col))
    }

    sample_col <- de_cfg$sample_col %||% de_cfg$id_columns$sample_col %||% "SampleID"
    if (!sample_col %in% colnames(meta)) {
        if (!is.null(rownames(meta)) && !any(rownames(meta) == "")) {
            meta[[sample_col]] <- rownames(meta)
        } else {
            stop(sprintf("Sample column '%s' not found in metadata", sample_col))
        }
    }
    meta <- align_samples_strict(colnames(counts), meta, sample_col,
                                 lenient = identical(de_cfg$sample_alignment, "lenient"))
    counts <- counts[, rownames(meta), drop = FALSE]

    meta[[factor_col]] <- factor(meta[[factor_col]])
    safe_levels <- make.names(levels(meta[[factor_col]]))
    levels(meta[[factor_col]]) <- safe_levels

    list(counts = counts, meta = meta, factor_col = factor_col,
         safe_levels = safe_levels, sample_col = sample_col)
}
