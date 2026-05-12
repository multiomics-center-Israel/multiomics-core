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

#' Convert p-value to significance stars
#'
#' @param p_value Numeric p-value.
#' @return Character string: "***", "**", "*", or "".
add_significance_stars <- function(p_value) {
    ifelse(p_value < 0.001, "***",
    ifelse(p_value < 0.01, "**",
    ifelse(p_value < 0.05, "*", "")))
}


#' Compute per-sample intensities grouped by a category and bond_type
#'
#' Internal helper for chain saturation and chain length computations.
#'
#' @param mat        Expression matrix (features x samples).
#' @param condition  Character vector of group labels per sample.
#' @param bond_type  Character vector of bond types per feature.
#' @param category   Character vector of category labels per feature.
#' @param cat_col    Name for the category column in output.
#' @return data.frame with columns: bond_type, <cat_col>, group, sample_id,
#'         intensity, count
.compute_chain_per_sample <- function(mat, condition, bond_type, category,
                                      cat_col = "category") {
    results <- list()
    for (bt in unique(bond_type)) {
        for (cat in unique(category)) {
            idx <- which(bond_type == bt & category == cat)
            if (length(idx) == 0) next
            for (g in unique(condition)) {
                g_idx <- which(condition == g)
                for (si in seq_along(g_idx)) {
                    s <- g_idx[si]
                    if (length(idx) == 1) {
                        int_val <- mat[idx, s]
                    } else {
                        int_val <- sum(mat[idx, s], na.rm = TRUE)
                    }
                    row <- data.frame(
                        bond_type  = bt,
                        group      = g,
                        sample_id  = colnames(mat)[s],
                        intensity  = int_val,
                        count      = length(idx),
                        stringsAsFactors = FALSE
                    )
                    row[[cat_col]] <- cat
                    results[[length(results) + 1]] <- row
                }
            }
        }
    }
    do.call(rbind, results)
}


#' Compute significance annotations for faceted chain plots
#'
#' Runs a t-test per (bond_type, category) comparing groups.
#'
#' @param per_sample  data.frame from .compute_chain_per_sample().
#' @param cat_col     Name of the category column.
#' @return data.frame with bond_type, <cat_col>, label, y_position
.compute_chain_significance <- function(per_sample, cat_col = "category") {
    groups <- unique(per_sample$group)
    if (length(groups) != 2) return(NULL)

    combos <- unique(per_sample[, c("bond_type", cat_col)])
    annot <- list()

    for (i in seq_len(nrow(combos))) {
        bt  <- combos$bond_type[i]
        cat <- combos[[cat_col]][i]
        sub <- per_sample[per_sample$bond_type == bt &
                          per_sample[[cat_col]] == cat, ]
        vals_a <- sub$intensity[sub$group == groups[1]]
        vals_b <- sub$intensity[sub$group == groups[2]]
        if (length(vals_a) < 2 || length(vals_b) < 2) next
        if (sd(vals_a, na.rm = TRUE) == 0 && sd(vals_b, na.rm = TRUE) == 0) next
        pval <- tryCatch(
            stats::t.test(vals_a, vals_b)$p.value,
            error = function(e) NA_real_
        )
        if (is.na(pval)) next
        star <- add_significance_stars(pval)
        if (nzchar(star)) {
            y_pos <- max(c(vals_a, vals_b), na.rm = TRUE) * 1.05
            row <- data.frame(
                bond_type  = bt,
                label      = star,
                y_position = y_pos,
                stringsAsFactors = FALSE
            )
            row[[cat_col]] <- cat
            annot[[length(annot) + 1]] <- row
        }
    }
    if (length(annot) == 0) return(NULL)
    do.call(rbind, annot)
}


#' Compute chain saturation summary per group, split by bond type
#'
#' Categorizes lipids by saturation level (SFA, MUFA, PUFA) and computes
#' abundance per group, split by bond type (Acyl, Alkyl, Alkenyl).
#'
#' @param pre    Pre-processed lipidomics list.
#' @param config Full config.
#' @return data.frame with columns: bond_type, saturation, group,
#'         mean_intensity, count
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

    # Bond type: fall back to "Acyl" if column missing
    bond_type <- if ("bond_type" %in% colnames(row_data)) {
        as.character(row_data$bond_type)
    } else {
        rep("Acyl", nrow(row_data))
    }

    per_sample <- .compute_chain_per_sample(
        mat, condition, bond_type, saturation, cat_col = "saturation"
    )

    # Aggregate to mean per (bond_type, saturation, group)
    agg <- stats::aggregate(
        intensity ~ bond_type + saturation + group,
        data = per_sample, FUN = mean, na.rm = TRUE
    )
    names(agg)[names(agg) == "intensity"] <- "mean_intensity"

    # Add count (number of features per combo, take first occurrence)
    cnt <- stats::aggregate(
        count ~ bond_type + saturation + group,
        data = per_sample[!duplicated(paste(per_sample$bond_type,
                                             per_sample$saturation,
                                             per_sample$group)), ],
        FUN = max
    )
    agg <- merge(agg, cnt, by = c("bond_type", "saturation", "group"),
                 all.x = TRUE)

    # Store per-sample data as attribute for significance testing in plots
    attr(agg, "per_sample") <- per_sample
    agg
}


#' Compute chain length distribution per group, split by bond type
#'
#' @param pre    Pre-processed lipidomics list.
#' @param config Full config.
#' @return data.frame with columns: bond_type, chain_length_bin, group,
#'         mean_intensity, count
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
    bins <- ifelse(is.na(carbons), "Unknown",
            ifelse(carbons <= 20, "Short (<=20C)",
            ifelse(carbons <= 40, "Medium (21-40C)",
            ifelse(carbons <= 60, "Long (41-60C)", "Very Long (>60C)"))))

    # Bond type: fall back to "Acyl" if column missing
    bond_type <- if ("bond_type" %in% colnames(row_data)) {
        as.character(row_data$bond_type)
    } else {
        rep("Acyl", nrow(row_data))
    }

    per_sample <- .compute_chain_per_sample(
        mat, condition, bond_type, bins, cat_col = "chain_length_bin"
    )

    # Aggregate to mean per (bond_type, chain_length_bin, group)
    agg <- stats::aggregate(
        intensity ~ bond_type + chain_length_bin + group,
        data = per_sample, FUN = mean, na.rm = TRUE
    )
    names(agg)[names(agg) == "intensity"] <- "mean_intensity"

    cnt <- stats::aggregate(
        count ~ bond_type + chain_length_bin + group,
        data = per_sample[!duplicated(paste(per_sample$bond_type,
                                             per_sample$chain_length_bin,
                                             per_sample$group)), ],
        FUN = max
    )
    agg <- merge(agg, cnt, by = c("bond_type", "chain_length_bin", "group"),
                 all.x = TRUE)

    attr(agg, "per_sample") <- per_sample
    agg
}


#' Plot chain saturation barplot, faceted by bond type
#'
#' @param sat_data data.frame from compute_chain_saturation().
#' @return ggplot object.
plot_chain_saturation <- function(sat_data) {
    per_sample <- attr(sat_data, "per_sample")
    sat_data <- sat_data[sat_data$saturation != "Unknown", ]

    # Filter to bond types that have data
    keep_bt <- unique(sat_data$bond_type[sat_data$mean_intensity > 0])
    sat_data <- sat_data[sat_data$bond_type %in% keep_bt, ]

    # Restructure: x = saturation category, fill = group (so brackets go between groups)
    p <- ggplot2::ggplot(sat_data,
                    ggplot2::aes(x = saturation, y = mean_intensity, fill = group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                          width = 0.65) +
        ggplot2::labs(
            title = "Chain Saturation Profile",
            x = NULL, y = "Mean Total Intensity", fill = "Group"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            strip.text = ggplot2::element_text(face = "bold", size = 11)
        )

    # Facet if more than one bond type
    if (length(keep_bt) > 1) {
        p <- p + ggplot2::facet_wrap(~ bond_type, scales = "free_y")
    }

    # Add significance brackets between the two group bars
    if (!is.null(per_sample)) {
        per_sample <- per_sample[per_sample$saturation != "Unknown" &
                                  per_sample$bond_type %in% keep_bt, ]
        sig_annot <- .compute_chain_significance(per_sample,
                                                  cat_col = "saturation")
        if (!is.null(sig_annot) && nrow(sig_annot) > 0) {
            sig_annot <- sig_annot[sig_annot$bond_type %in% keep_bt, ]
            if (nrow(sig_annot) > 0) {
                # Position bracket between the two dodged bars
                dodge_w <- 0.7
                cats <- levels(factor(sat_data$saturation))
                if (is.null(cats)) cats <- unique(sat_data$saturation)
                sig_annot$x_cat <- match(sig_annot$saturation, cats)
                sig_annot$xmin <- sig_annot$x_cat - dodge_w / 2 + 0.05
                sig_annot$xmax <- sig_annot$x_cat + dodge_w / 2 - 0.05
                sig_annot$x <- sig_annot$x_cat
                # Bracket line
                p <- p + ggplot2::geom_segment(
                    data = sig_annot,
                    ggplot2::aes(x = xmin, xend = xmax,
                                 y = y_position, yend = y_position),
                    inherit.aes = FALSE, linewidth = 0.4
                ) +
                ggplot2::geom_text(
                    data = sig_annot,
                    ggplot2::aes(x = x, y = y_position, label = label),
                    inherit.aes = FALSE, size = 5, vjust = -0.3
                )
            }
        }
    }

    p
}


#' Plot chain length distribution barplot, faceted by bond type
#'
#' @param len_data data.frame from compute_chain_length_distribution().
#' @return ggplot object.
plot_chain_length <- function(len_data) {
    per_sample <- attr(len_data, "per_sample")
    len_data <- len_data[len_data$chain_length_bin != "Unknown", ]
    len_data$chain_length_bin <- factor(len_data$chain_length_bin,
        levels = c("Short (<=20C)", "Medium (21-40C)",
                   "Long (41-60C)", "Very Long (>60C)"))

    # Filter to bond types that have data
    keep_bt <- unique(len_data$bond_type[len_data$mean_intensity > 0])
    len_data <- len_data[len_data$bond_type %in% keep_bt, ]

    # Restructure: x = chain length category, fill = group (so brackets go between groups)
    p <- ggplot2::ggplot(len_data,
                    ggplot2::aes(x = chain_length_bin, y = mean_intensity,
                                 fill = group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                          width = 0.65) +
        ggplot2::labs(
            title = "Chain Length Distribution",
            x = NULL, y = "Mean Total Intensity", fill = "Group"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            strip.text = ggplot2::element_text(face = "bold", size = 11)
        )

    # Facet if more than one bond type
    if (length(keep_bt) > 1) {
        p <- p + ggplot2::facet_wrap(~ bond_type, scales = "free_y")
    }

    # Add significance brackets between the two group bars
    if (!is.null(per_sample)) {
        per_sample <- per_sample[per_sample$chain_length_bin != "Unknown" &
                                  per_sample$bond_type %in% keep_bt, ]
        sig_annot <- .compute_chain_significance(per_sample,
                                                  cat_col = "chain_length_bin")
        if (!is.null(sig_annot) && nrow(sig_annot) > 0) {
            sig_annot <- sig_annot[sig_annot$bond_type %in% keep_bt, ]
            if (nrow(sig_annot) > 0) {
                dodge_w <- 0.7
                cats <- levels(len_data$chain_length_bin)
                if (is.null(cats)) cats <- unique(len_data$chain_length_bin)
                sig_annot$x_cat <- match(sig_annot$chain_length_bin, cats)
                sig_annot$xmin <- sig_annot$x_cat - dodge_w / 2 + 0.05
                sig_annot$xmax <- sig_annot$x_cat + dodge_w / 2 - 0.05
                sig_annot$x <- sig_annot$x_cat
                p <- p + ggplot2::geom_segment(
                    data = sig_annot,
                    ggplot2::aes(x = xmin, xend = xmax,
                                 y = y_position, yend = y_position),
                    inherit.aes = FALSE, linewidth = 0.4
                ) +
                ggplot2::geom_text(
                    data = sig_annot,
                    ggplot2::aes(x = x, y = y_position, label = label),
                    inherit.aes = FALSE, size = 5, vjust = -0.3
                )
            }
        }
    }

    p
}


# ==== PER-CLASS CHAIN PROFILES ================================================

#' Extract individual fatty acid chains from a lipid name
#'
#' Parses names like "PE(18:1_18:2)" or "TG(16:0/18:1/18:2)" into
#' individual chain species vectors: c("18:1", "18:2"), etc.
#'
#' @param lipid_name Character scalar.
#' @return Character vector of chain species (e.g. c("16:0", "18:1")).
extract_individual_chains <- function(lipid_name) {
    # Remove class prefix: "PE(18:1_18:2)" -> "18:1_18:2"
    chains_str <- sub("^[A-Za-z()\\-]+\\(", "", lipid_name)
    chains_str <- sub("\\)$", "", chains_str)
    # Split on _ or /
    chains <- unlist(strsplit(chains_str, "[_/]"))
    # Keep only valid C:DB patterns
    chains <- grep("^\\d+:\\d+$", chains, value = TRUE)
    chains
}


#' Compute per-class chain profiles with group comparison
#'
#' For each lipid class, extracts individual chain species from lipid names,
#' aggregates chain abundance per experimental group, and tests for
#' significance (Wilcoxon rank-sum or t-test).
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return Named list keyed by lipid class.
#'   Each element is a data.frame with columns:
#'   chain, group, mean_intensity, se, p_value.
compute_per_class_chain_profiles <- function(pre, config) {
    cfg <- config$modes$lipidomics
    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    mat      <- pre$expr_filt
    meta     <- pre$meta
    row_data <- pre$row_data

    meta      <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])
    groups    <- unique(condition)

    classes   <- as.character(row_data$lipid_class)
    lipid_ids <- as.character(row_data$feature_id)
    unique_classes <- sort(unique(classes[!is.na(classes)]))

    # Map feature_id -> lipid name (prefer display_name if available)
    name_col <- if ("display_name" %in% colnames(row_data)) "display_name" else
                if ("Name" %in% colnames(row_data)) "Name" else "feature_id"
    lipid_names <- as.character(row_data[[name_col]])

    results <- list()

    for (cl in unique_classes) {
        cl_idx <- which(classes == cl)
        if (length(cl_idx) == 0) next

        # Build per-chain, per-sample intensity matrix
        chain_sample_list <- list()

        for (i in cl_idx) {
            chains <- extract_individual_chains(lipid_names[i])
            if (length(chains) == 0) next
            for (ch in chains) {
                if (is.null(chain_sample_list[[ch]])) {
                    chain_sample_list[[ch]] <- rep(0, ncol(mat))
                }
                chain_sample_list[[ch]] <- chain_sample_list[[ch]] + mat[i, ]
            }
        }

        if (length(chain_sample_list) == 0) next

        chain_names <- names(chain_sample_list)
        rows <- list()

        for (ch in chain_names) {
            vals <- chain_sample_list[[ch]]

            # Test between first two groups (primary comparison)
            p_val <- NA_real_
            if (length(groups) == 2) {
                g1_vals <- vals[condition == groups[1]]
                g2_vals <- vals[condition == groups[2]]
                p_val <- tryCatch({
                    if (length(g1_vals) >= 3 && length(g2_vals) >= 3) {
                        stats::wilcox.test(g1_vals, g2_vals)$p.value
                    } else {
                        stats::t.test(g1_vals, g2_vals)$p.value
                    }
                }, error = function(e) NA_real_)
            }

            for (g in groups) {
                g_idx <- which(condition == g)
                g_vals <- vals[g_idx]
                rows[[length(rows) + 1]] <- data.frame(
                    chain          = ch,
                    group          = g,
                    mean_intensity = mean(g_vals, na.rm = TRUE),
                    se             = if (length(g_vals) > 1) {
                        stats::sd(g_vals, na.rm = TRUE) / sqrt(length(g_vals))
                    } else 0,
                    p_value        = p_val,
                    stringsAsFactors = FALSE
                )
            }
        }

        if (length(rows) > 0) {
            results[[cl]] <- do.call(rbind, rows)
        }
    }

    results
}


#' Plot per-class chain profile as a grouped barplot
#'
#' Creates a dodged barplot of chain abundances per group with error bars
#' and significance stars.
#'
#' @param chain_df  data.frame from one element of compute_per_class_chain_profiles().
#' @param class_name Character scalar – the lipid class name for the title.
#' @return ggplot object.
plot_per_class_chain_profile <- function(chain_df, class_name) {
    # Sort chains by total carbon:double-bond numerically
    chain_order <- unique(chain_df$chain)
    chain_parts <- do.call(rbind, lapply(chain_order, function(ch) {
        parts <- as.integer(unlist(strsplit(ch, ":")))
        if (length(parts) == 2) parts else c(NA_integer_, NA_integer_)
    }))
    chain_order <- chain_order[order(chain_parts[, 1], chain_parts[, 2])]
    chain_df$chain <- factor(chain_df$chain, levels = chain_order)

    # Significance stars
    sig_df <- unique(chain_df[, c("chain", "p_value")])
    sig_df$star <- ifelse(is.na(sig_df$p_value), "",
                   ifelse(sig_df$p_value < 0.001, "***",
                   ifelse(sig_df$p_value < 0.01, "**",
                   ifelse(sig_df$p_value < 0.05, "*", ""))))

    # Y position for stars: max of (mean + se) per chain
    y_max <- stats::aggregate(
        cbind(ymax = mean_intensity + se) ~ chain,
        data = chain_df, FUN = max
    )
    sig_df <- merge(sig_df, y_max, by = "chain", all.x = TRUE)
    sig_df$ypos <- sig_df$ymax * 1.08

    p <- ggplot2::ggplot(chain_df,
                          ggplot2::aes(x = chain, y = mean_intensity, fill = group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = mean_intensity - se, ymax = mean_intensity + se),
            position = ggplot2::position_dodge(width = 0.8),
            width = 0.25
        ) +
        ggplot2::labs(
            title = paste0("Chain Profile \u2014 ", class_name),
            x = "Fatty Acid Chain", y = "Mean Intensity", fill = "Group"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title   = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 10)
        )

    # Add significance stars
    sig_df <- sig_df[sig_df$star != "", , drop = FALSE]
    if (nrow(sig_df) > 0) {
        p <- p + ggplot2::geom_text(
            data    = sig_df,
            mapping = ggplot2::aes(x = chain, y = ypos, label = star),
            inherit.aes = FALSE,
            size = 5, vjust = 0
        )
    }

    p
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


#' DE-aware structural bubble: chain length vs. double bonds
#'
#' One point per significant lipid: x = total carbons, y = total double bonds,
#' size = -log10(p), colour = direction (Up/Down by sign of logFC).
#' Significance: `p_col` <= `sig_cutoff` AND |logFC| >= `logfc_cutoff`.
#'
#' @param de_tbl       Per-contrast DE table with columns feature_id, logFC,
#'                     P.Value, adj.P.Val (output of extract_contrast_table()).
#' @param row_data     Feature annotations with feature_id, total_carbons,
#'                     total_double_bonds.
#' @param sig_cutoff   p-value threshold (default 0.05).
#' @param logfc_cutoff |logFC| threshold (default 0 = no FC filter).
#' @param p_col        Which p column to use ("adj.P.Val" or "P.Value").
#' @param title        Optional plot title.
#' @return ggplot object, or NULL when no significant annotated lipids remain.
plot_de_chain_bubble <- function(de_tbl, row_data,
                                  sig_cutoff   = 0.05,
                                  logfc_cutoff = 0,
                                  p_col        = "adj.P.Val",
                                  title        = NULL) {
    if (is.null(de_tbl) || nrow(de_tbl) == 0L) return(NULL)
    if (!"feature_id" %in% colnames(de_tbl) ||
        !"logFC"      %in% colnames(de_tbl)) {
        return(NULL)
    }
    if (!p_col %in% colnames(de_tbl)) {
        warning("plot_de_chain_bubble: '", p_col, "' not in de_tbl; falling back to P.Value.")
        p_col <- "P.Value"
        if (!p_col %in% colnames(de_tbl)) return(NULL)
    }
    if (is.null(row_data) ||
        !all(c("feature_id", "total_carbons", "total_double_bonds")
             %in% colnames(row_data))) {
        return(NULL)
    }

    keep_cols <- c("feature_id", "total_carbons", "total_double_bonds",
                   intersect(c("lipid_class", "bond_type"), colnames(row_data)))
    rd <- row_data[, keep_cols, drop = FALSE]
    df <- merge(de_tbl, rd, by = "feature_id", all.x = TRUE)

    pv <- df[[p_col]]
    sig <- !is.na(pv) & pv <= sig_cutoff &
        !is.na(df$logFC) & abs(df$logFC) >= logfc_cutoff &
        !is.na(df$total_carbons) & !is.na(df$total_double_bonds)
    df <- df[sig, , drop = FALSE]
    if (nrow(df) == 0L) return(NULL)

    df$direction <- ifelse(df$logFC > 0, "Up", "Down")
    df$nlp <- -log10(pmax(df[[p_col]], .Machine$double.eps))

    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    title <- title %||% sprintf(
        "Significant lipid structures (%s ≤ %.3g, |logFC| ≥ %.2f)",
        p_col, sig_cutoff, logfc_cutoff
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$total_carbons,
                                       y = .data$total_double_bonds,
                                       size = .data$nlp,
                                       colour = .data$direction)) +
        ggplot2::geom_jitter(alpha = 0.75, width = 0.2, height = 0.2) +
        ggplot2::scale_colour_manual(
            values = c(Up = "#d73027", Down = "#4575b4"),
            name   = "Direction"
        ) +
        ggplot2::scale_size_continuous(
            range = c(2, 8),
            name  = expression(-log[10] * "(p)")
        ) +
        ggplot2::labs(
            title = title,
            x = "Total carbon atoms",
            y = "Total double bonds"
        ) +
        ggplot2::theme_bw()
}
