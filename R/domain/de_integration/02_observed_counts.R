# R/domain/de_integration/02_observed_counts.R
#
# How many values each group actually measured, per feature and contrast, so
# later steps can set aside fold changes that rest mostly on imputation.

#' Numerator and denominator groups of a contrast
#'
#' From the layer's contrasts file when it has one, matched on
#' \code{normalize_contrast_key()}; otherwise from a label of the form
#' "<A>_vs_<B>".
#'
#' @param contrast Contrast label.
#' @param contrasts_df Contrasts table (\code{Contrast_name}, \code{Factor},
#'   \code{Numerator}, \code{Denominator}), or NULL.
#' @return List with \code{numerator}, \code{denominator} and \code{factor}
#'   (NA where unknown), or NULL.
contrast_groups <- function(contrast, contrasts_df = NULL) {
    if (is.data.frame(contrasts_df) && "Contrast_name" %in% names(contrasts_df)) {
        i <- which(normalize_contrast_key(contrasts_df$Contrast_name) ==
                   normalize_contrast_key(contrast))
        if (length(i) == 1) {
            return(list(numerator = as.character(contrasts_df$Numerator[i]),
                        denominator = as.character(contrasts_df$Denominator[i]),
                        factor = as.character(contrasts_df$Factor[i] %||% NA)))
        }
    }
    parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) == 2 && all(nzchar(parts))) {
        return(list(numerator = parts[1], denominator = parts[2], factor = NA_character_))
    }
    NULL
}


#' Observed (non-imputed) values per group for one contrast
#'
#' A fold change resting on one measured value per group is mostly imputation,
#' and later steps flag it through \code{well_observed}. Counts come from the
#' table's own \code{N.observed.<group>} / \code{n_obs.<group>} columns when it
#' has them, otherwise from the layer's unimputed matrix and sample sheet
#' (\code{observed} block), using the grouping \code{compute_group_stat_columns()}
#' applies to the exported means. Without either, they are NA.
#'
#' @param ly The layer's config.
#' @param df The layer's table.
#' @param ids Feature ids, one per row of \code{df}.
#' @param contrast Contrast label.
#' @param config Full config, for path resolution.
#' @param observed Pre-read observed inputs from \code{read_observed_inputs()},
#'   or NULL.
#' @return List with \code{num}, \code{den} (integer vectors, one per row) and
#'   \code{source}.
layer_observed_counts <- function(ly, df, ids, contrast, config, observed = NULL) {
    none <- list(num = rep(NA_integer_, nrow(df)), den = rep(NA_integer_, nrow(df)),
                 source = "not available")

    if (identical(ly$format, "generic")) {
        cols <- ly$columns %||% list()
        if (!is.null(cols$n_obs_num) && !is.null(cols$n_obs_den) &&
            all(c(cols$n_obs_num, cols$n_obs_den) %in% names(df))) {
            return(list(num = as.integer(df[[cols$n_obs_num]]),
                        den = as.integer(df[[cols$n_obs_den]]),
                        source = paste(cols$n_obs_num, cols$n_obs_den, sep = ", ")))
        }
        return(none)
    }

    grp <- contrast_groups(contrast, observed$contrasts)
    if (is.null(grp)) return(none)

    cols <- de_layer_columns(ly$format, contrast)
    for (prefix in cols$n_obs_prefixes) {
        num_col <- paste0(prefix, grp$numerator)
        den_col <- paste0(prefix, grp$denominator)
        if (all(c(num_col, den_col) %in% names(df))) {
            return(list(num = as.integer(df[[num_col]]), den = as.integer(df[[den_col]]),
                        source = paste(num_col, den_col, sep = ", ")))
        }
    }

    if (is.null(observed)) return(none)
    counts <- compute_group_stat_columns(
        expr = observed$matrix, sample_meta = observed$samplesheet,
        sample_id_col = ly$observed$sample_col,
        contrasts_df = observed$contrasts[
            normalize_contrast_key(observed$contrasts$Contrast_name) ==
                normalize_contrast_key(contrast), , drop = FALSE],
        stat_fn = function(m) rowSums(!is.na(m)),
        prefix = "N.observed.")
    num_col <- paste0("N.observed.", grp$numerator)
    den_col <- paste0("N.observed.", grp$denominator)
    if (is.null(counts) || !all(c(num_col, den_col) %in% names(counts))) return(none)
    idx <- match(ids, rownames(counts))
    list(num = as.integer(counts[[num_col]][idx]), den = as.integer(counts[[den_col]][idx]),
         source = paste0("unimputed matrix ", basename(ly$observed$matrix)))
}


#' Read a layer's observed-count inputs once
#'
#' @param ly The layer's config (with an \code{observed} block).
#' @param config Full config, for path resolution.
#' @return List with \code{matrix} (features x samples, NA where unobserved),
#'   \code{samplesheet} and \code{contrasts}; NULL when the layer has no
#'   \code{observed} block.
read_observed_inputs <- function(ly, config) {
    obs <- ly$observed
    if (is.null(obs)) return(NULL)
    mat <- read_table_auto(resolve_input_path(config, obs$matrix))
    id_col <- obs$id_col %||% names(mat)[1]
    ids <- as.character(mat[[id_col]])
    mat <- as.matrix(mat[, setdiff(names(mat), id_col), drop = FALSE])
    storage.mode(mat) <- "double"
    rownames(mat) <- ids
    list(
        matrix = mat,
        samplesheet = read_samplesheet(resolve_input_path(config, obs$samplesheet)),
        contrasts = read_table_auto(resolve_input_path(config, obs$contrasts_file))
    )
}
