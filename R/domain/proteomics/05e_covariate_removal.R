# R/domain/proteomics/05e_covariate_removal.R
#
# The DE matrix with the fitted de$covariates term taken out and the condition
# effect kept. One place builds it, and every consumer reads it from here: the
# whole-matrix sample heatmaps, the PCA before/after panels, the clustering
# heatmap and the Mean.adj.<group> columns of the results tables. Each of them
# then shows the matrix the pooled limma model describes, rather than values
# that still carry the term the model adjusted for.

#' The matrix the proteomics DE model was fitted on
#'
#' The first imputation draw, which is also the \code{.norm} block of
#' final_results.tsv; the preprocessing matrix when the DE step made no draws
#' (precomputed DE tables). Same rule as \code{mod_proteomics_exports()}.
#'
#' @param pre Preprocessed proteomics object.
#' @param de_res DE result from \code{mod_proteomics_de()}.
#' @return Numeric matrix, features x samples.
resolve_proteomics_de_matrix <- function(pre, de_res) {
    if (length(de_res$imputations) > 0) {
        as.matrix(de_res$imputations[[1]])
    } else {
        as.matrix(pre$expr_imp_single)
    }
}

#' Remove continuous covariate effects, keeping the group effect
#'
#' A thin wrapper around \code{limma::removeBatchEffect()} with
#' \code{design = model.matrix(~ 0 + group)}: each feature is fitted on
#' group + covariates and the covariate part of the fit is subtracted. With
#' covariates centred on their mean (as [resolve_de_covariates()] returns them),
#' the group means of the result equal the group coefficients of that fit, so
#' differences between them reproduce the model's contrasts.
#'
#' @param expr_mat Numeric matrix, features x samples, without missing values.
#' @param covariates Numeric matrix, one row per sample (in column order of
#'   \code{expr_mat}), one column per covariate.
#' @param group Vector of condition labels, one per sample.
#' @param block Optional blocking factor, passed to \code{limma::lmFit()}.
#' @param correlation Within-block correlation for \code{block}.
#' @return Matrix with the same dimensions and names as \code{expr_mat}.
#' @examples
#' x <- matrix(rnorm(40), nrow = 4, dimnames = list(paste0("p", 1:4), paste0("s", 1:10)))
#' remove_covariate_effects(x, cbind(depth = scale(1:10, scale = FALSE)), rep(c("a", "b"), 5))
remove_covariate_effects <- function(expr_mat, covariates, group,
                                     block = NULL, correlation = NULL) {
    expr_mat <- as.matrix(expr_mat)
    cov_mat <- matrix(as.numeric(covariates), nrow = NROW(covariates),
                      dimnames = list(NULL, colnames(covariates)))
    if (nrow(cov_mat) != ncol(expr_mat) || length(group) != ncol(expr_mat)) {
        stop("remove_covariate_effects(): covariates and group need one entry per sample ",
             "(", ncol(expr_mat), " columns in the matrix); got ", nrow(cov_mat),
             " covariate rows and ", length(group), " group labels.", call. = FALSE)
    }
    if (anyNA(expr_mat)) {
        stop("remove_covariate_effects(): the matrix has missing values. It expects the ",
             "imputed DE matrix, which has none.", call. = FALSE)
    }
    design <- stats::model.matrix(~ 0 + factor(group))
    if (is.null(block)) {
        limma::removeBatchEffect(expr_mat, covariates = cov_mat, design = design)
    } else {
        limma::removeBatchEffect(expr_mat, covariates = cov_mat, design = design,
                                 block = block, correlation = correlation)
    }
}

#' Remove the configured DE covariates from a proteomics matrix
#'
#' Reads the covariates through [resolve_de_covariates()] and the condition
#' column through [get_proteomics_de_group_col()], as the limma fit does. When
#' \code{de$block_col} is also set, the within-block correlation is estimated
#' the way [run_limma_proteomics()] estimates it, so the removed slope is the
#' fitted one.
#'
#' @param expr_mat Numeric matrix, features x samples (an imputation draw).
#' @param meta Sample metadata.
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @return List: \code{mat} (covariates removed, or the input when none are
#'   configured), \code{before} (the input aligned to the metadata),
#'   \code{meta} (metadata aligned to the columns) and \code{covariates}
#'   (the centred covariate matrix, or NULL).
remove_proteomics_de_covariates <- function(expr_mat, meta, p_cfg) {
    sample_col <- p_cfg$effects$samples %||% "SampleID"
    expr_mat <- align_matrix_to_meta(as.matrix(expr_mat), meta, sample_col)
    meta_aligned <- align_meta_to_expr(expr_mat, meta, p_cfg)

    covariates <- resolve_de_covariates(meta_aligned, p_cfg)
    if (is.null(covariates)) {
        return(list(mat = expr_mat, before = expr_mat, meta = meta_aligned, covariates = NULL))
    }

    group <- factor(meta_aligned[[get_proteomics_de_group_col(p_cfg)]])
    block <- resolve_de_block(meta_aligned, p_cfg, sample_col)
    correlation <- NULL
    if (!is.null(block)) {
        design <- build_proteomics_de_design(group, covariates = covariates)
        correlation <- limma::duplicateCorrelation(expr_mat, design, block = block)$consensus
        # The fit drops the block on a non-finite consensus; do the same here.
        if (!is.finite(correlation)) {
            block <- NULL
            correlation <- NULL
        }
    }

    adjusted <- remove_covariate_effects(expr_mat, covariates, group,
                                         block = block, correlation = correlation)
    list(mat = adjusted, before = expr_mat, meta = meta_aligned, covariates = covariates)
}

#' The matrix the report's sample-level figures are drawn from
#'
#' The DE matrix ([resolve_proteomics_de_matrix()]), with the de$covariates
#' term removed when the DE model fitted one (evidence: \code{de_res$covariate_design},
#' which [run_limma_proteomics()] fills only then). Without covariates the DE
#' matrix is returned unchanged.
#'
#' @param pre Preprocessed proteomics object.
#' @param de_res DE result from \code{mod_proteomics_de()}.
#' @param config Full pipeline config.
#' @return The list from [remove_proteomics_de_covariates()] plus
#'   \code{observed} (logical matrix, TRUE where the value was measured before
#'   imputation, aligned to \code{mat}) and \code{label} (a short description of
#'   the matrix for figure captions).
resolve_proteomics_display_matrix <- function(pre, de_res, config) {
    p_cfg <- config$modes$proteomics
    de_mat <- resolve_proteomics_de_matrix(pre, de_res)

    res <- if (is.null(de_res$covariate_design)) {
        sample_col <- p_cfg$effects$samples %||% "SampleID"
        m <- align_matrix_to_meta(de_mat, pre$meta, sample_col)
        list(mat = m, before = m, meta = align_meta_to_expr(m, pre$meta, p_cfg),
             covariates = NULL)
    } else {
        remove_proteomics_de_covariates(de_mat, pre$meta, p_cfg)
    }

    # pre$expr_filt is the matrix every DE imputation draw starts from, so its
    # NAs mark exactly the cells the DE matrix filled in.
    filt <- as.matrix(pre$expr_filt)
    observed <- !is.na(filt[rownames(res$mat), colnames(res$mat), drop = FALSE])

    res$observed <- observed
    res$label <- if (is.null(res$covariates)) {
        "DE matrix"
    } else {
        sprintf("DE matrix, covariate%s %s removed (limma::removeBatchEffect, group kept)",
                if (ncol(res$covariates) > 1) "s" else "",
                paste(colnames(res$covariates), collapse = ", "))
    }
    res
}

#' Covariate-adjusted per-group means for the proteomics results table
#'
#' \code{Mean.adj.<group>}: the group means of each imputation draw after
#' [remove_proteomics_de_covariates()], averaged over the draws. Within one
#' draw the difference of two such means is that draw's fitted contrast, so
#' with one draw (\code{imputation$multi_imputation: false}) the differences
#' reproduce \code{log2FC.imputs}. With several draws they reproduce the mean
#' of the per-draw log2FC, which differs from \code{log2FC.imputs} (the log2 of
#' the mean linear ratio) by the jensen_gap the DE_reconciliation sheet reports.
#'
#' @param pre Preprocessed proteomics object (its \code{meta} is used).
#' @param de_res DE result from \code{mod_proteomics_de()}.
#' @param contrasts_df Contrasts table.
#' @param config Full pipeline config.
#' @return Feature-indexed data.frame of \code{Mean.adj.<group>} columns, or
#'   NULL when no covariate was fitted, no draw exists, or the group summary
#'   block is switched off (\code{excel$group_cv: false}, which also drops
#'   \code{Mean.<group>}).
build_group_mean_adjusted_proteomics <- function(pre, de_res, contrasts_df, config) {
    p_cfg <- config$modes$proteomics %||% list()
    if (is.null(de_res$covariate_design) || length(de_res$imputations) == 0) return(NULL)
    if (!isTRUE(p_cfg$excel$group_cv %||% TRUE)) return(NULL)

    sample_id_col <- p_cfg$effects$samples %||% p_cfg$id_columns$sample_col %||% "SampleID"
    per_draw <- lapply(de_res$imputations, function(x) {
        adj <- remove_proteomics_de_covariates(x, pre$meta, p_cfg)$mat
        compute_group_mean_columns(adj, pre$meta, sample_id_col, contrasts_df,
                                   prefix = "Mean.adj.")
    })
    if (is.null(per_draw[[1]])) return(NULL)
    if (length(per_draw) == 1) return(per_draw[[1]])

    pooled <- Reduce(`+`, lapply(per_draw, as.matrix)) / length(per_draw)
    as.data.frame(pooled, stringsAsFactors = FALSE)
}

#' One-line description of the proteomics limma design
#'
#' Used in the Excel "How to read" notes when the model carried covariates.
#'
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @param covariates Covariate column names that entered the model.
#' @return A single string, e.g. \code{"~ 0 + group + depth"}.
describe_proteomics_de_design <- function(p_cfg, covariates) {
    terms <- c("0", get_proteomics_de_group_col(p_cfg), covariates)
    out <- paste("~", paste(terms, collapse = " + "))
    block_col <- p_cfg$de[["block_col"]]
    if (!is.null(block_col) && nzchar(block_col)) {
        out <- paste0(out, ", with '", block_col, "' as a blocking factor (duplicateCorrelation)")
    }
    out
}
