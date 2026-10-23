#' Run edgeR quasi-likelihood differential expression on RNA-seq counts
#'
#' TMM-normalises the (already gene-filtered) count matrix, estimates
#' negative-binomial dispersions with [edgeR::estimateDisp()], fits a
#' quasi-likelihood GLM with [edgeR::glmQLFit()] and tests each contrast with
#' [edgeR::glmQLFTest()].
#'
#' edgeR has no random effect, so a configured `de$block_col` always enters as
#' a fixed term (`~ block + factor`), as in the DESeq2 path. Each contrast is
#' fitted with its denominator as the reference level and tested as a single
#' coefficient; the dispersions are estimated once, since every such
#' reparametrisation spans the same column space.
#'
#' Returns the same structure as [run_deseq2_de()]. `baseMean` carries edgeR's
#' average log2-CPM, `lfcSE` is NA (the QL F-test reports no standard error),
#' and `stat` is the signed square root of the 1-df QL F statistic, a t-like
#' value used to rank genes for GSEA.
#'
#' @param counts Integer count matrix (genes x samples). tximport input is not
#'   supported here.
#' @param meta Data frame with sample metadata.
#' @param contrasts_df Data frame with columns: Contrast_name, Factor, Numerator, Denominator.
#' @param de_cfg DE configuration list (sample_col, sample_alignment, block_col).
#' @return List with 'dds' (NULL), 'tables' (named list of DE result data
#'   frames) and 'info' (method, design, factor_col, block_col, block_model,
#'   n_blocks, edger_version).
run_edger_rna_de <- function(counts, meta, contrasts_df, de_cfg) {
    prep <- prepare_count_de_inputs(counts, meta, contrasts_df, de_cfg,
                                    method = "edger_ql", pkgs = "edgeR")
    counts <- prep$counts; meta <- prep$meta
    factor_col <- prep$factor_col; safe_levels <- prep$safe_levels

    block <- resolve_de_block(meta, list(de = de_cfg), prep$sample_col)
    if (!is.null(block)) {
        meta$.block <- block
        design <- check_fixed_block_design(stats::as.formula(paste0("~ .block + ", factor_col)),
                                           meta, de_cfg$block_col, factor_col)
        message(sprintf("[run_edger_rna_de] Fixed block '%s' (%d levels) in the design",
                        de_cfg$block_col, nlevels(block)))
    } else {
        design <- stats::model.matrix(stats::as.formula(paste0("~ ", factor_col)), data = meta)
    }

    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
    dge <- edgeR::estimateDisp(dge, design)

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
        if (identical(num, den)) {
            warning(sprintf("Contrast '%s' compares level '%s' with itself. Skipping.",
                            cn, contrasts_df$Numerator[i]))
            next
        }
        lev <- stats::relevel(meta[[factor_col]], ref = den)
        design_cn <- if (is.null(block)) stats::model.matrix(~ lev) else stats::model.matrix(~ block + lev)
        # Full "lev" column name: a level such as "block2" must not match a block column.
        fit <- edgeR::glmQLFit(dge, design_cn)
        qlf <- edgeR::glmQLFTest(fit, coef = paste0("lev", num))
        de <- edgeR::topTags(qlf, n = Inf, sort.by = "none", adjust.method = "BH")$table
        tab <- data.frame(
            baseMean       = de$logCPM,
            log2FoldChange = de$logFC,
            lfcSE          = NA_real_,
            stat           = sign(de$logFC) * sqrt(de$F),
            pvalue         = de$PValue,
            padj           = de$FDR,
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
            method        = "edger_ql",
            design        = if (is.null(block)) paste0("~", factor_col)
                            else paste0("~", de_cfg$block_col, " + ", factor_col),
            factor_col    = factor_col,
            sample_col    = prep$sample_col,
            block_col     = if (is.null(block)) NULL else de_cfg$block_col,
            block_model   = if (is.null(block)) NULL else "fixed",
            n_blocks      = if (is.null(block)) NA_integer_ else nlevels(block),
            edger_version = as.character(utils::packageVersion("edgeR"))
        )
    )
}
