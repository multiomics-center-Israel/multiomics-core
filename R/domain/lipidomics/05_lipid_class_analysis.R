# R/domain/lipidomics/05_lipid_class_analysis.R
#
# Lipidomics-specific analyses:
#   1. Lipid class composition (absolute and normalized)
#   2. Chain length and saturation profiling
#   3. Lipid class enrichment (ORA on DE significant lipids)
#
# These analyses leverage the structural hierarchy unique to lipidomics.


# ==== LIPID CLASS COMPOSITION ================================================

#' Compute lipid class composition per sample group
#'
#' Aggregates expression by lipid class and computes both absolute and
#' normalized (within-sample percentage) compositions.
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(class_abs, class_norm, class_summary, per_sample)
compute_lipid_class_composition <- function(pre, config) {
    cfg <- config$modes$lipidomics

    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    # Use raw filtered intensities (not log-transformed)
    mat <- pre$expr_filt
    meta <- pre$meta
    row_data <- pre$row_data

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])

    classes <- row_data$lipid_class

    # Per-sample class totals
    unique_classes <- sort(unique(classes))
    class_mat <- matrix(0, nrow = length(unique_classes), ncol = ncol(mat),
                        dimnames = list(unique_classes, colnames(mat)))

    for (cl in unique_classes) {
        idx <- which(classes == cl)
        if (length(idx) == 1) {
            class_mat[cl, ] <- mat[idx, ]
        } else {
            class_mat[cl, ] <- colSums(mat[idx, , drop = FALSE], na.rm = TRUE)
        }
    }

    # Normalized: each sample sums to 100%
    col_totals <- colSums(class_mat, na.rm = TRUE)
    col_totals[col_totals == 0] <- 1
    class_mat_norm <- sweep(class_mat, 2, col_totals, FUN = "/") * 100

    # Per-group mean
    groups <- unique(condition)
    class_abs_mean <- sapply(groups, function(g) {
        idx <- which(condition == g)
        if (length(idx) == 1) class_mat[, idx] else rowMeans(class_mat[, idx, drop = FALSE], na.rm = TRUE)
    })
    colnames(class_abs_mean) <- groups

    class_norm_mean <- sapply(groups, function(g) {
        idx <- which(condition == g)
        if (length(idx) == 1) class_mat_norm[, idx] else rowMeans(class_mat_norm[, idx, drop = FALSE], na.rm = TRUE)
    })
    colnames(class_norm_mean) <- groups

    # Feature count per class
    class_count <- table(classes)

    list(
        class_abs       = class_abs_mean,
        class_norm      = class_norm_mean,
        class_count     = class_count,
        per_sample_abs  = class_mat,
        per_sample_norm = class_mat_norm,
        classes         = unique_classes
    )
}


#' Plot lipid class composition barplot
#'
#' @param class_data List from compute_lipid_class_composition().
#' @param type       "absolute" or "normalized".
#' @param title      Plot title.
#' @return ggplot object.
plot_lipid_class_barplot <- function(class_data, type = "normalized",
                                     title = NULL) {
    mat <- if (type == "normalized") class_data$class_norm else class_data$class_abs

    df_long <- data.frame(
        lipid_class = rep(rownames(mat), ncol(mat)),
        group       = rep(colnames(mat), each = nrow(mat)),
        value       = as.vector(mat),
        stringsAsFactors = FALSE
    )

    if (is.null(title)) {
        title <- if (type == "normalized") {
            "Lipid Class Composition (Normalized %)"
        } else {
            "Lipid Class Composition (Absolute Intensity)"
        }
    }

    y_lab <- if (type == "normalized") "Percentage (%)" else "Mean Intensity"

    ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = value, fill = lipid_class)) +
        ggplot2::geom_col(position = "stack", width = 0.7) +
        ggplot2::labs(title = title, x = NULL, y = y_lab, fill = "Lipid Class") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            axis.text.x = ggplot2::element_text(size = 11),
            legend.position = "right"
        )
}


# ==== CHAIN PROFILING =========================================================

#' Compute chain saturation summary per group
#'
#' Categorizes lipids by saturation level (SFA, MUFA, PUFA) and computes
#' abundance per group.
#'
#' @param pre    Pre-processed lipidomics list.
#' @param config Full config.
#' @return data.frame with columns: saturation, group, mean_intensity, count
compute_chain_saturation <- function(pre, config) {
    cfg <- config$modes$lipidomics
    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat <- pre$expr_filt
    meta <- pre$meta
    row_data <- pre$row_data

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])

    dbs <- row_data$total_double_bonds
    saturation <- ifelse(is.na(dbs), "Unknown",
                  ifelse(dbs == 0, "SFA (0 DB)",
                  ifelse(dbs == 1, "MUFA (1 DB)", "PUFA (2+ DB)")))

    groups <- unique(condition)
    results <- list()

    for (sat in c("SFA (0 DB)", "MUFA (1 DB)", "PUFA (2+ DB)", "Unknown")) {
        idx <- which(saturation == sat)
        if (length(idx) == 0) next
        for (g in groups) {
            g_idx <- which(condition == g)
            if (length(idx) == 1) {
                mean_int <- mean(mat[idx, g_idx], na.rm = TRUE)
            } else {
                mean_int <- mean(colSums(mat[idx, g_idx, drop = FALSE], na.rm = TRUE))
            }
            results[[length(results) + 1]] <- data.frame(
                saturation = sat, group = g,
                mean_intensity = mean_int, count = length(idx),
                stringsAsFactors = FALSE
            )
        }
    }

    do.call(rbind, results)
}


#' Compute chain length distribution per group
#'
#' @param pre    Pre-processed lipidomics list.
#' @param config Full config.
#' @return data.frame with columns: chain_length_bin, group, mean_intensity, count
compute_chain_length_distribution <- function(pre, config) {
    cfg <- config$modes$lipidomics
    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat <- pre$expr_filt
    meta <- pre$meta
    row_data <- pre$row_data

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])

    carbons <- row_data$total_carbons
    # Bin into categories
    bins <- ifelse(is.na(carbons), "Unknown",
            ifelse(carbons <= 20, "Short (<=20C)",
            ifelse(carbons <= 40, "Medium (21-40C)",
            ifelse(carbons <= 60, "Long (41-60C)", "Very Long (>60C)"))))

    groups <- unique(condition)
    results <- list()

    for (b in c("Short (<=20C)", "Medium (21-40C)", "Long (41-60C)",
                "Very Long (>60C)", "Unknown")) {
        idx <- which(bins == b)
        if (length(idx) == 0) next
        for (g in groups) {
            g_idx <- which(condition == g)
            if (length(idx) == 1) {
                mean_int <- mean(mat[idx, g_idx], na.rm = TRUE)
            } else {
                mean_int <- mean(colSums(mat[idx, g_idx, drop = FALSE], na.rm = TRUE))
            }
            results[[length(results) + 1]] <- data.frame(
                chain_length_bin = b, group = g,
                mean_intensity = mean_int, count = length(idx),
                stringsAsFactors = FALSE
            )
        }
    }

    do.call(rbind, results)
}


#' Plot chain saturation barplot
#'
#' @param sat_data data.frame from compute_chain_saturation().
#' @return ggplot object.
plot_chain_saturation <- function(sat_data) {
    sat_data <- sat_data[sat_data$saturation != "Unknown", ]

    ggplot2::ggplot(sat_data,
                    ggplot2::aes(x = group, y = mean_intensity, fill = saturation)) +
        ggplot2::geom_col(position = "dodge", width = 0.7) +
        ggplot2::scale_fill_manual(values = c(
            "SFA (0 DB)"   = "#4DAF4A",
            "MUFA (1 DB)"  = "#377EB8",
            "PUFA (2+ DB)" = "#E41A1C"
        )) +
        ggplot2::labs(
            title = "Chain Saturation Profile",
            x = NULL, y = "Mean Total Intensity", fill = "Saturation"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )
}


#' Plot chain length distribution barplot
#'
#' @param len_data data.frame from compute_chain_length_distribution().
#' @return ggplot object.
plot_chain_length <- function(len_data) {
    len_data <- len_data[len_data$chain_length_bin != "Unknown", ]
    len_data$chain_length_bin <- factor(len_data$chain_length_bin,
        levels = c("Short (<=20C)", "Medium (21-40C)",
                   "Long (41-60C)", "Very Long (>60C)"))

    ggplot2::ggplot(len_data,
                    ggplot2::aes(x = group, y = mean_intensity, fill = chain_length_bin)) +
        ggplot2::geom_col(position = "dodge", width = 0.7) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
            title = "Chain Length Distribution",
            x = NULL, y = "Mean Total Intensity", fill = "Chain Length"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )
}


# ==== LIPID CLASS ENRICHMENT =================================================

#' Lipid class over-representation analysis
#'
#' Tests whether DE-significant lipids are enriched in specific classes
#' using Fisher's exact test (one per class).
#'
#' @param de_res   DE result list (must have summary_df with pass_any_contrast).
#' @param row_data Row data with lipid_class column.
#' @return data.frame with class, n_class, n_sig, pct_sig, odds_ratio,
#'   p_value, p_adjusted
lipid_class_ora <- function(de_res, row_data) {
    if (is.null(de_res) || is.null(de_res$summary_df)) return(NULL)

    summary_df <- de_res$summary_df
    is_sig <- summary_df$pass_any_contrast == 1

    # Map feature_id -> lipid_class
    class_map <- stats::setNames(
        as.character(row_data$lipid_class),
        as.character(row_data$feature_id)
    )

    feat_classes <- class_map[summary_df$feature_id]
    unique_classes <- sort(unique(feat_classes[!is.na(feat_classes)]))

    n_total <- length(feat_classes)
    n_sig_total <- sum(is_sig, na.rm = TRUE)

    if (n_sig_total == 0) {
        message("lipid class ORA: no significant features — skipping")
        return(NULL)
    }

    results <- lapply(unique_classes, function(cl) {
        in_class <- !is.na(feat_classes) & feat_classes == cl
        # 2x2 table: sig/not-sig x in-class/not-in-class
        a <- sum(is_sig & in_class)
        b <- sum(is_sig & !in_class)
        c_val <- sum(!is_sig & in_class)
        d <- sum(!is_sig & !in_class)

        ft <- stats::fisher.test(matrix(c(a, b, c_val, d), nrow = 2))

        data.frame(
            lipid_class = cl,
            n_class     = sum(in_class),
            n_sig       = a,
            pct_sig     = round(100 * a / max(sum(in_class), 1), 1),
            odds_ratio  = round(ft$estimate, 3),
            p_value     = ft$p.value,
            stringsAsFactors = FALSE
        )
    })

    result_df <- do.call(rbind, results)
    result_df$p_adjusted <- stats::p.adjust(result_df$p_value, method = "BH")
    result_df <- result_df[order(result_df$p_value), ]
    rownames(result_df) <- NULL

    result_df
}


#' Plot lipid class ORA results
#'
#' @param ora_df data.frame from lipid_class_ora().
#' @param top_n  Number of classes to show.
#' @return ggplot object or NULL.
plot_lipid_class_ora <- function(ora_df, top_n = 20) {
    if (is.null(ora_df) || nrow(ora_df) == 0) return(NULL)

    top_df <- utils::head(ora_df, top_n)
    top_df$lipid_class <- factor(top_df$lipid_class,
                                  levels = rev(top_df$lipid_class))
    top_df$neglog10p <- -log10(pmax(top_df$p_adjusted, 1e-300))

    ggplot2::ggplot(top_df,
                    ggplot2::aes(x = neglog10p, y = lipid_class, size = n_sig)) +
        ggplot2::geom_point(color = "steelblue", alpha = 0.8) +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                             color = "grey50") +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(
            title = "Lipid Class Over-Representation (Fisher's Exact)",
            x = "-log10(FDR)", y = NULL, size = "# Significant"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )
}
