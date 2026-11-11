# R/domain/proteomics/07g_sample_figures.R
#
# Inputs shared by the whole-matrix sample figures (correlation and z-score
# heatmaps, PCA before/after covariate removal) and the DE clustering heatmap:
# which metadata columns annotate the samples, their colours, and which
# proteins enter a top-N-by-variance set.

#' Metadata columns that annotate samples in the proteomics heatmaps and PCA
#'
#' The single reader of \code{effects$heatmap_annotations}, the key the RNA-seq
#' and metabolomics modes already use for the same purpose. Every proteomics
#' heatmap with a sample annotation goes through it.
#'
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @return Character vector of column names, or NULL when the key is unset
#'   (callers then keep their condition-column annotation unchanged).
get_proteomics_sample_annotation_cols <- function(p_cfg) {
    cols <- (p_cfg$effects %||% list())[["heatmap_annotations"]]
    if (is.null(cols) || length(cols) == 0) return(NULL)
    as.character(unlist(cols))
}

#' Order the levels of one sample-annotation column
#'
#' One rule for every annotation bar:
#' \enumerate{
#'   \item values that all start with a number ("0", "4", "24", or "0h", "4h",
#'     "24h") are ordered by that number, so 4h comes before 24h;
#'   \item otherwise a negative/positive pair (no/yes, absent/present,
#'     false/true, any case) puts the negative first;
#'   \item otherwise the order of first appearance in the metadata.
#' }
#'
#' @param x A metadata column.
#' @return A factor with the levels in that order.
#' @examples
#' order_annotation_levels(c("24h", "0h", "4h", "0h"))
order_annotation_levels <- function(x) {
    chr <- as.character(x)
    vals <- unique(chr[!is.na(chr)])
    lead <- suppressWarnings(as.numeric(sub("^\\s*([-+]?[0-9]*\\.?[0-9]+).*$", "\\1", vals)))
    starts_numeric <- grepl("^\\s*[-+]?[0-9]*\\.?[0-9]+", vals)
    lv <- if (length(vals) > 0 && all(starts_numeric) && !anyNA(lead)) {
        vals[order(lead, seq_along(vals))]
    } else {
        neg <- c("no", "absent", "false")
        pos <- c("yes", "present", "true")
        low <- tolower(vals)
        pair <- which(neg %in% low & pos %in% low)
        if (length(vals) == 2 && length(pair) == 1) {
            vals[order(match(low, c(neg[pair], pos[pair])))]
        } else {
            vals
        }
    }
    factor(chr, levels = lv)
}

#' Sample annotation data.frame for pheatmap and PCA
#'
#' With \code{effects$heatmap_annotations} set, those columns in that order,
#' each a factor ordered by [order_annotation_levels()]. Columns missing from
#' the metadata are dropped with a warning; if none is left, the colour column
#' is used, as the RNA-seq mode does. Without the key it is
#' \code{build_heatmap_annotation_col()}, the annotation the clustering heatmap
#' has always used, so nothing changes for other projects.
#'
#' @param meta Sample metadata.
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @return data.frame with sample IDs as rownames, or NULL.
build_proteomics_sample_annotation <- function(meta, p_cfg) {
    cols <- get_proteomics_sample_annotation_cols(p_cfg)
    if (is.null(cols)) return(build_heatmap_annotation_col(meta, p_cfg))

    missing <- setdiff(cols, colnames(meta))
    if (length(missing) > 0) {
        warning("effects$heatmap_annotations names column(s) not in the sample metadata: ",
                paste(missing, collapse = ", "), ". ",
                if (length(missing) == length(cols)) "Using the colour column instead." else
                    "They are left out.", call. = FALSE)
        cols <- setdiff(cols, missing)
        if (length(cols) == 0) return(build_heatmap_annotation_col(meta, p_cfg))
    }
    sample_col <- p_cfg$effects$samples %||% p_cfg$id_columns$sample_col %||% "SampleID"
    annot <- as.data.frame(lapply(meta[, cols, drop = FALSE], order_annotation_levels),
                           stringsAsFactors = FALSE)
    colnames(annot) <- cols
    rownames(annot) <- as.character(meta[[sample_col]])
    annot
}

#' Colours for the sample annotation bars
#'
#' The palettes of the exploratory sample heatmaps: a light-to-dark purple
#' ramp for a column whose values start with numbers (time), grey then brown
#' for a no/yes-type pair (in level order). Any other column is left to
#' pheatmap's default colours.
#'
#' @param annot Output of [build_proteomics_sample_annotation()].
#' @return Named list of named colour vectors, as pheatmap's
#'   \code{annotation_colors} expects, or NULL when no column has a palette.
build_sample_annotation_colors <- function(annot) {
    if (is.null(annot) || ncol(annot) == 0) return(NULL)
    out <- lapply(names(annot), function(col) {
        lv <- levels(factor(annot[[col]]))
        n <- length(lv)
        if (n >= 2 && all(grepl("^\\s*[-+]?[0-9]*\\.?[0-9]+", lv))) {
            return(stats::setNames(
                grDevices::colorRampPalette(c("#DADAEB", "#9E9AC8", "#54278F"))(n), lv))
        }
        if (n == 2 && any(tolower(lv[1]) == c("no", "absent", "false"))) {
            return(stats::setNames(c("#BDBDBD", "#B8864B"), lv))
        }
        NULL
    })
    names(out) <- names(annot)
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0) NULL else out
}

#' Put the annotation columns in pheatmap's drawing order
#'
#' pheatmap draws the first annotation column at the bottom of the bar stack
#' and lists it last in the legend. Reversing the columns at the pheatmap
#' boundary puts the first configured column (e.g. time) on top and first in
#' the legend. Only when \code{effects$heatmap_annotations} is set, so the
#' default annotation of other projects is drawn exactly as before; the PCA
#' keeps the configured order (colour = first column, shape = second).
#'
#' @param annot Annotation data.frame from [build_proteomics_sample_annotation()].
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @return \code{annot} with its columns reversed, or unchanged.
annotation_for_pheatmap <- function(annot, p_cfg) {
    if (is.null(annot) || ncol(annot) < 2 ||
        is.null(get_proteomics_sample_annotation_cols(p_cfg))) {
        return(annot)
    }
    annot[, rev(seq_len(ncol(annot))), drop = FALSE]
}

#' Sample annotation and colours for one proteomics heatmap
#'
#' What every pre-existing proteomics QC heatmap calls: NULL when
#' \code{effects$heatmap_annotations} is unset, so the caller keeps the
#' annotation it has always drawn; otherwise the resolver's annotation aligned
#' to \code{sample_ids} and in pheatmap's drawing order
#' ([annotation_for_pheatmap()]), with its colours.
#'
#' @param meta Sample metadata.
#' @param p_cfg The \code{modes$proteomics} config branch.
#' @param sample_ids Column names of the matrix to be drawn.
#' @return NULL, or list(annot, colors).
proteomics_heatmap_annotation <- function(meta, p_cfg, sample_ids) {
    if (is.null(get_proteomics_sample_annotation_cols(p_cfg))) return(NULL)
    annot <- build_proteomics_sample_annotation(meta, p_cfg)
    if (is.null(annot)) return(NULL)
    annot <- annot[intersect(sample_ids, rownames(annot)), , drop = FALSE]
    annot <- annotation_for_pheatmap(annot, p_cfg)
    list(annot = annot, colors = build_sample_annotation_colors(annot))
}

#' Proteins for a whole-matrix sample figure
#'
#' Only fully measured proteins are eligible (observed in every sample, so no
#' imputed value enters the set or its ranking). With \code{n} given, the
#' \code{n} with the largest variance across samples, ranked on
#' \code{expr_mat}; ties keep their input order.
#'
#' @param expr_mat Numeric matrix, features x samples (the matrix drawn).
#' @param observed Logical matrix of the same shape, TRUE where measured.
#' @param n Number of proteins, or NULL for every fully measured protein.
#' @return Character vector of feature IDs: in decreasing variance for a
#'   top-N set, in matrix order for the full set.
#' @examples
#' x <- matrix(c(1, 2, 3, 1, 5, 9, 2, 2, 3), nrow = 3, byrow = TRUE,
#'             dimnames = list(c("a", "b", "c"), c("s1", "s2", "s3")))
#' obs <- x > 0
#' obs["b", "s2"] <- FALSE
#' select_top_variable_measured(x, obs, n = 1)  # "a": "b" is not fully measured
select_top_variable_measured <- function(expr_mat, observed, n = NULL) {
    expr_mat <- as.matrix(expr_mat)
    observed <- as.matrix(observed)
    if (!identical(dim(expr_mat), dim(observed))) {
        stop("select_top_variable_measured(): expr_mat and observed must have the same shape.",
             call. = FALSE)
    }
    complete <- select_complete_case_features(!observed)
    ids <- rownames(expr_mat)[complete]
    if (is.null(n)) return(ids)
    v <- apply(expr_mat[complete, , drop = FALSE], 1, stats::var)
    ids[order(-v, seq_along(v))][seq_len(min(n, length(ids)))]
}
