# =============================================================================
# Advanced Statistical Methods for Proteomics (multiomics-core)
# =============================================================================
# Ported from multiomics/R/07c_advanced_stats.R with config adaptations:
#   log_message(...) → message(...)
#   config$design$condition_column → cfg$effects$color
#   config$design$random_effects → cfg$advanced_stats$random_effects
#   config$output$output_dir → output_dir parameter
#   da_results is a named list of per-contrast DE tables (from extract_de_table_for_pathway)

# -----------------------------------------------------------------------------
# Main Orchestrator
# -----------------------------------------------------------------------------

#' Run Advanced Statistical Analysis
#'
#' @param da_results      Named list of per-contrast DE tables, or single data.frame
#' @param normalized_data Numeric matrix (proteins x samples)
#' @param metadata        Sample metadata data.frame
#' @param config          Full pipeline config list
#' @param output_dir      Output directory
#' @return List of advanced stats results
#' @export
run_advanced_stats <- function(da_results, normalized_data, metadata, config,
                                output_dir) {
    message("=== Running Advanced Statistical Analysis ===")

    cfg <- config$modes$proteomics
    stats_cfg <- cfg$advanced_stats %||% list()

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plots_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    # Flatten da_results into a single data.frame (use first contrast if list)
    da_table <- NULL
    if (is.data.frame(da_results)) {
        da_table <- da_results
    } else if (is.list(da_results) && length(da_results) > 0) {
        da_table <- da_results[[1]]
    }

    if (!is.null(da_table) && !"protein_id" %in% colnames(da_table)) {
        if ("FeatureID" %in% colnames(da_table))
            da_table$protein_id <- da_table$FeatureID
        else if (!is.null(rownames(da_table)))
            da_table$protein_id <- rownames(da_table)
    }

    # Map column names
    if (!is.null(da_table)) {
        if (!"log2FoldChange" %in% colnames(da_table) && "logFC" %in% colnames(da_table))
            da_table$log2FoldChange <- da_table$logFC
        if (!"pvalue" %in% colnames(da_table) && "P.Value" %in% colnames(da_table))
            da_table$pvalue <- da_table$P.Value
        if (!"padj" %in% colnames(da_table) && "adj.P.Val" %in% colnames(da_table))
            da_table$padj <- da_table$adj.P.Val
    }

    results <- list()

    # 1. Effect sizes with bootstrap CIs
    if (isTRUE(stats_cfg$compute_effect_size_ci %||% TRUE)) {
        results$effect_sizes <- tryCatch(
            compute_effect_sizes_with_ci(normalized_data, metadata, da_table, config, output_dir),
            error = function(e) { message("  Effect size computation failed: ", e$message); NULL }
        )
    }

    # 2. Mixed effects models
    if (has_random_effects(metadata, config)) {
        results$mixed_effects <- tryCatch(
            run_mixed_effects_analysis(normalized_data, metadata, config, output_dir),
            error = function(e) { message("  Mixed effects failed: ", e$message); NULL }
        )
    }

    # 3. Robust regression
    if (isTRUE(stats_cfg$run_robust_regression)) {
        results$robust_regression <- tryCatch(
            run_robust_regression(normalized_data, metadata, config, output_dir),
            error = function(e) { message("  Robust regression failed: ", e$message); NULL }
        )
    }

    # 4. Power analysis
    if (isTRUE(stats_cfg$run_power_analysis %||% TRUE)) {
        results$power_analysis <- tryCatch(
            run_power_analysis(da_table, metadata, config, output_dir),
            error = function(e) { message("  Power analysis failed: ", e$message); NULL }
        )
    }

    # 5. Multiple testing comparison
    results$mtc_comparison <- tryCatch(
        compare_multiple_testing_methods(da_table, output_dir),
        error = function(e) { message("  MTC comparison failed: ", e$message); NULL }
    )

    # Generate plots
    figures <- generate_advanced_stats_plots(results, da_table, config, plots_dir)
    results$figures <- figures

    # Summary
    results$summary <- create_advanced_stats_summary(results)

    # Save
    save_advanced_stats_outputs(results, output_dir)
    message("Advanced statistical analysis complete!")
    results
}

# -----------------------------------------------------------------------------
# Effect Size with Bootstrap CIs
# -----------------------------------------------------------------------------

compute_effect_sizes_with_ci <- function(normalized_data, metadata, da_results,
                                          config, output_dir) {
    message("Computing effect sizes with bootstrap CIs...")

    cfg <- config$modes$proteomics
    condition_col <- cfg$effects$color %||% "Condition"
    ref_level <- cfg$de$reference_level %||% NULL

    if (!condition_col %in% colnames(metadata)) {
        message("WARNING: Condition column '", condition_col, "' not found.")
        return(NULL)
    }

    conditions <- metadata[[condition_col]]
    if (is.null(ref_level)) ref_level <- levels(factor(conditions))[1]

    n_boot   <- cfg$advanced_stats$bootstrap_n %||% 1000
    ci_level <- cfg$advanced_stats$ci_level %||% 0.95

    proteins <- if (!is.null(da_results)) da_results$protein_id else rownames(normalized_data)
    proteins <- intersect(proteins, rownames(normalized_data))

    if (length(proteins) == 0) return(NULL)

    message("  Processing ", length(proteins), " proteins...")

    effect_results <- lapply(proteins, function(protein) {
        tryCatch(
            compute_protein_effect_size(protein, normalized_data, conditions,
                                        ref_level, n_boot, ci_level),
            error = function(e) {
                data.frame(protein_id = protein, effect_size = NA,
                           effect_size_type = NA, ci_lower = NA, ci_upper = NA,
                           se = NA, stringsAsFactors = FALSE)
            }
        )
    })

    effect_df <- do.call(rbind, effect_results)

    # Merge DA results
    if (!is.null(da_results)) {
        merge_cols <- intersect(c("protein_id", "log2FoldChange", "pvalue", "padj"),
                                colnames(da_results))
        if (length(merge_cols) > 1) {
            effect_df <- merge(effect_df, da_results[, merge_cols, drop = FALSE],
                               by = "protein_id", all.x = TRUE)
        }
    }

    effect_df$significant <- effect_df$ci_lower > 0 | effect_df$ci_upper < 0

    write.csv(effect_df, file.path(output_dir, "effect_sizes_with_ci.csv"), row.names = FALSE)
    message("  Computed effect sizes for ", sum(!is.na(effect_df$effect_size)), " proteins")
    effect_df
}

compute_protein_effect_size <- function(protein, expr_data, conditions,
                                         ref_level, n_boot, ci_level) {
    if (!(protein %in% rownames(expr_data))) {
        return(data.frame(protein_id = protein, effect_size = NA,
                          effect_size_type = NA, ci_lower = NA, ci_upper = NA,
                          se = NA, stringsAsFactors = FALSE))
    }

    values <- as.numeric(expr_data[protein, ])
    unique_conds <- unique(conditions)

    if (length(unique_conds) == 2) {
        group1 <- values[conditions == ref_level]
        group2 <- values[conditions != ref_level]
        group1 <- group1[!is.na(group1)]
        group2 <- group2[!is.na(group2)]
        if (length(group1) < 2 || length(group2) < 2)
            return(data.frame(protein_id = protein, effect_size = NA,
                              effect_size_type = NA, ci_lower = NA, ci_upper = NA,
                              se = NA, stringsAsFactors = FALSE))

        effect_size <- compute_hedges_g(group1, group2)
        boot_effects <- replicate(n_boot, {
            compute_hedges_g(sample(group1, replace = TRUE), sample(group2, replace = TRUE))
        })
        boot_effects <- boot_effects[is.finite(boot_effects)]

        if (length(boot_effects) < 10) {
            ci_lower <- ci_upper <- se <- NA
        } else {
            alpha <- 1 - ci_level
            ci_lower <- quantile(boot_effects, alpha / 2, na.rm = TRUE)
            ci_upper <- quantile(boot_effects, 1 - alpha / 2, na.rm = TRUE)
            se <- sd(boot_effects, na.rm = TRUE)
        }

        data.frame(protein_id = protein, effect_size = effect_size,
                   effect_size_type = "Hedges_g", ci_lower = ci_lower,
                   ci_upper = ci_upper, se = se, stringsAsFactors = FALSE)

    } else {
        df <- data.frame(y = values, group = factor(conditions))
        df <- df[complete.cases(df), ]
        if (nrow(df) < 5)
            return(data.frame(protein_id = protein, effect_size = NA,
                              effect_size_type = NA, ci_lower = NA, ci_upper = NA,
                              se = NA, stringsAsFactors = FALSE))

        effect_size <- compute_eta_squared(df$y, df$group)
        boot_effects <- replicate(n_boot, {
            idx <- sample(nrow(df), replace = TRUE)
            compute_eta_squared(df$y[idx], df$group[idx])
        })
        boot_effects <- boot_effects[is.finite(boot_effects)]

        if (length(boot_effects) < 10) {
            ci_lower <- ci_upper <- se <- NA
        } else {
            alpha <- 1 - ci_level
            ci_lower <- quantile(boot_effects, alpha / 2, na.rm = TRUE)
            ci_upper <- quantile(boot_effects, 1 - alpha / 2, na.rm = TRUE)
            se <- sd(boot_effects, na.rm = TRUE)
        }

        data.frame(protein_id = protein, effect_size = effect_size,
                   effect_size_type = "eta_squared", ci_lower = ci_lower,
                   ci_upper = ci_upper, se = se, stringsAsFactors = FALSE)
    }
}

compute_hedges_g <- function(group1, group2) {
    group1 <- group1[!is.na(group1)]
    group2 <- group2[!is.na(group2)]
    n1 <- length(group1); n2 <- length(group2)
    if (n1 < 2 || n2 < 2) return(NA_real_)

    var1 <- var(group1); var2 <- var(group2)
    if (is.na(var1) || is.na(var2) || !is.finite(var1) || !is.finite(var2)) return(NA_real_)

    pooled_var <- ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    if (is.na(pooled_var) || !is.finite(pooled_var) || pooled_var < 0) return(NA_real_)

    pooled_sd <- sqrt(pooled_var)
    if (pooled_sd == 0 || !is.finite(pooled_sd)) return(NA_real_)

    d <- (mean(group2) - mean(group1)) / pooled_sd
    correction <- 1 - (3 / (4 * (n1 + n2) - 9))
    d * correction
}

compute_eta_squared <- function(y, group) {
    fit <- tryCatch(aov(y ~ group), error = function(e) NULL)
    if (is.null(fit)) return(NA)
    ss <- summary(fit)[[1]]
    ss["group", "Sum Sq"] / sum(ss[, "Sum Sq"])
}

# -----------------------------------------------------------------------------
# Mixed Effects Models
# -----------------------------------------------------------------------------

has_random_effects <- function(metadata, config) {
    cfg <- config$modes$proteomics
    random_cols <- cfg$advanced_stats$random_effects %||% NULL

    if (is.null(random_cols)) {
        potential <- c("subject", "patient", "donor", "replicate", "batch", "block")
        random_cols <- intersect(potential, colnames(metadata))
    }

    for (col in random_cols) {
        if (col %in% colnames(metadata) && any(table(metadata[[col]]) > 1))
            return(TRUE)
    }
    FALSE
}

run_mixed_effects_analysis <- function(normalized_data, metadata, config, output_dir) {
    message("Running mixed effects models...")

    if (!requireNamespace("lme4", quietly = TRUE)) {
        message("WARNING: lme4 not available. Skipping.")
        return(NULL)
    }

    cfg <- config$modes$proteomics
    condition_col <- cfg$effects$color %||% "Condition"
    random_cols <- cfg$advanced_stats$random_effects %||% c("subject", "patient", "donor")
    random_cols <- intersect(random_cols, colnames(metadata))
    if (length(random_cols) == 0) return(NULL)

    random_col <- random_cols[1]
    message("  Fixed: ", condition_col, " | Random: ", random_col)

    proteins <- rownames(normalized_data)

    mixed_results <- lapply(proteins, function(protein) {
        tryCatch(
            fit_mixed_model(protein, normalized_data, metadata, condition_col, random_col),
            error = function(e) {
                data.frame(protein_id = protein, fixed_effect = NA, fixed_se = NA,
                           fixed_pvalue = NA, random_variance = NA,
                           residual_variance = NA, icc = NA, converged = FALSE,
                           stringsAsFactors = FALSE)
            }
        )
    })

    mixed_df <- do.call(rbind, mixed_results)
    mixed_df$padj <- p.adjust(mixed_df$fixed_pvalue, method = "BH")
    mixed_df <- mixed_df[order(mixed_df$padj), ]

    write.csv(mixed_df, file.path(output_dir, "mixed_effects_results.csv"), row.names = FALSE)
    message("  Converged: ", sum(mixed_df$converged, na.rm = TRUE), "/", nrow(mixed_df))
    mixed_df
}

fit_mixed_model <- function(protein, expr_data, metadata, condition_col, random_col) {
    y <- as.numeric(expr_data[protein, ])
    df <- data.frame(y = y, condition = factor(metadata[[condition_col]]),
                     random = factor(metadata[[random_col]]))
    df <- df[complete.cases(df), ]
    if (nrow(df) < 5 || length(unique(df$random)) < 2)
        return(data.frame(protein_id = protein, fixed_effect = NA, fixed_se = NA,
                          fixed_pvalue = NA, random_variance = NA,
                          residual_variance = NA, icc = NA, converged = FALSE,
                          stringsAsFactors = FALSE))

    model <- lme4::lmer(y ~ condition + (1 | random), data = df)
    fixed_coef <- lme4::fixef(model)
    fixed_se <- sqrt(diag(as.matrix(vcov(model))))
    var_comp <- as.data.frame(lme4::VarCorr(model))
    random_var <- var_comp$vcov[var_comp$grp == "random"]
    resid_var <- var_comp$vcov[var_comp$grp == "Residual"]
    icc <- random_var / (random_var + resid_var)

    pvalue <- tryCatch({
        if (requireNamespace("lmerTest", quietly = TRUE)) {
            model_test <- lmerTest::lmer(y ~ condition + (1 | random), data = df)
            summary(model_test)$coefficients[2, "Pr(>|t|)"]
        } else {
            z <- fixed_coef[2] / fixed_se[2]
            2 * pnorm(-abs(z))
        }
    }, error = function(e) NA)

    data.frame(protein_id = protein, fixed_effect = fixed_coef[2], fixed_se = fixed_se[2],
               fixed_pvalue = pvalue, random_variance = random_var,
               residual_variance = resid_var, icc = icc, converged = TRUE,
               stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Robust Regression
# -----------------------------------------------------------------------------

run_robust_regression <- function(normalized_data, metadata, config, output_dir) {
    message("Running robust regression...")
    if (!requireNamespace("MASS", quietly = TRUE)) {
        message("WARNING: MASS not available. Skipping.")
        return(NULL)
    }
    cfg <- config$modes$proteomics
    condition_col <- cfg$effects$color %||% "Condition"
    proteins <- rownames(normalized_data)

    robust_results <- lapply(proteins, function(protein) {
        tryCatch(
            fit_robust_regression(protein, normalized_data, metadata, condition_col),
            error = function(e) {
                data.frame(protein_id = protein, robust_coef = NA, robust_se = NA,
                           robust_pvalue = NA, n_outliers = NA, stringsAsFactors = FALSE)
            }
        )
    })

    robust_df <- do.call(rbind, robust_results)
    robust_df$padj <- p.adjust(robust_df$robust_pvalue, method = "BH")
    write.csv(robust_df, file.path(output_dir, "robust_regression_results.csv"), row.names = FALSE)
    robust_df
}

fit_robust_regression <- function(protein, expr_data, metadata, condition_col) {
    y <- as.numeric(expr_data[protein, ])
    df <- data.frame(y = y, condition = factor(metadata[[condition_col]]))
    df <- df[complete.cases(df), ]
    if (nrow(df) < 5)
        return(data.frame(protein_id = protein, robust_coef = NA, robust_se = NA,
                          robust_pvalue = NA, n_outliers = NA, stringsAsFactors = FALSE))

    model <- MASS::rlm(y ~ condition, data = df)
    coefs <- summary(model)$coefficients
    if (nrow(coefs) < 2)
        return(data.frame(protein_id = protein, robust_coef = NA, robust_se = NA,
                          robust_pvalue = NA, n_outliers = NA, stringsAsFactors = FALSE))

    t_value <- coefs[2, 3]
    pvalue <- 2 * pt(-abs(t_value), nrow(df) - 2)

    data.frame(protein_id = protein, robust_coef = coefs[2, 1], robust_se = coefs[2, 2],
               robust_pvalue = pvalue, n_outliers = sum(model$w < 0.9), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Power Analysis
# -----------------------------------------------------------------------------

run_power_analysis <- function(da_results, metadata, config, output_dir) {
    message("Running power analysis...")

    cfg <- config$modes$proteomics
    condition_col <- cfg$effects$color %||% "Condition"

    if (!condition_col %in% colnames(metadata)) {
        message("  Condition column not found for power analysis.")
        return(NULL)
    }

    group_sizes <- table(metadata[[condition_col]])
    current_n <- min(group_sizes)

    if (is.null(da_results) || !is.data.frame(da_results)) {
        message("  No DA table for power analysis.")
        return(NULL)
    }

    padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(da_results))[1]
    fc_col   <- intersect(c("log2FoldChange", "logFC"), colnames(da_results))[1]
    if (is.na(padj_col) || is.na(fc_col)) {
        message("  Required columns not found for power analysis.")
        return(NULL)
    }

    significant <- !is.na(da_results[[padj_col]]) & da_results[[padj_col]] < 0.05
    effect_sizes <- abs(da_results[[fc_col]][significant])
    effect_sizes <- effect_sizes[!is.na(effect_sizes) & is.finite(effect_sizes)]

    if (length(effect_sizes) == 0) {
        message("  No significant proteins for power estimation.")
        return(NULL)
    }

    median_effect <- median(effect_sizes, na.rm = TRUE)
    if (is.na(median_effect) || !is.finite(median_effect)) return(NULL)

    message("  Current n per group: ", current_n, " | Median |log2FC|: ", round(median_effect, 2))

    sample_sizes <- c(3, 5, 8, 10, 15, 20, 30, 50)
    power_results <- lapply(sample_sizes, function(n) {
        effect_tests <- c(0.5, 1.0, 1.5, 2.0, median_effect)
        powers <- sapply(effect_tests, function(d) calculate_power(n, d))
        data.frame(n_per_group = n, effect_0.5 = powers[1], effect_1.0 = powers[2],
                   effect_1.5 = powers[3], effect_2.0 = powers[4],
                   effect_median = powers[5], stringsAsFactors = FALSE)
    })
    power_df <- do.call(rbind, power_results)

    required_n <- sapply(c(0.5, 1.0, 1.5, 2.0, median_effect), calculate_required_n)
    required_df <- data.frame(effect_size = c(0.5, 1.0, 1.5, 2.0, median_effect),
                              required_n = required_n, stringsAsFactors = FALSE)

    write.csv(power_df, file.path(output_dir, "power_by_sample_size.csv"), row.names = FALSE)
    write.csv(required_df, file.path(output_dir, "required_sample_sizes.csv"), row.names = FALSE)

    list(power_table = power_df, required_n = required_df,
         current_n = current_n, median_effect = median_effect)
}

calculate_power <- function(n, effect_size, alpha = 0.05) {
    if (!is.numeric(n) || !is.numeric(effect_size) || is.na(n) || is.na(effect_size))
        return(NA_real_)
    if (n < 2 || effect_size <= 0) return(NA_real_)
    df <- 2 * n - 2
    ncp <- effect_size * sqrt(n / 2)
    t_crit <- qt(1 - alpha / 2, df)
    1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)
}

calculate_required_n <- function(effect_size, power = 0.80, alpha = 0.05) {
    if (!is.numeric(effect_size) || is.na(effect_size) || effect_size <= 0) return(NA_integer_)
    for (n in seq(3, 500)) {
        p <- calculate_power(n, effect_size, alpha)
        if (!is.na(p) && p >= power) return(n)
    }
    NA_integer_
}

# -----------------------------------------------------------------------------
# Multiple Testing Comparison
# -----------------------------------------------------------------------------

compare_multiple_testing_methods <- function(da_results, output_dir) {
    message("Comparing multiple testing correction methods...")

    if (is.null(da_results) || !is.data.frame(da_results)) return(NULL)

    pval_col <- intersect(c("pvalue", "P.Value", "PValue", "p.value"), colnames(da_results))[1]
    if (is.na(pval_col)) return(NULL)

    pvalues <- da_results[[pval_col]][!is.na(da_results[[pval_col]])]
    if (length(pvalues) == 0) return(NULL)

    methods <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")
    mtc_comparison <- lapply(methods, function(method) {
        padj <- p.adjust(pvalues, method = method)
        data.frame(method = method,
                   n_significant_0.01 = sum(padj < 0.01, na.rm = TRUE),
                   n_significant_0.05 = sum(padj < 0.05, na.rm = TRUE),
                   n_significant_0.10 = sum(padj < 0.10, na.rm = TRUE),
                   median_padj = median(padj, na.rm = TRUE),
                   stringsAsFactors = FALSE)
    })
    mtc_df <- do.call(rbind, mtc_comparison)
    write.csv(mtc_df, file.path(output_dir, "mtc_comparison.csv"), row.names = FALSE)
    mtc_df
}

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

generate_advanced_stats_plots <- function(results, da_results, config, plots_dir) {
    message("Generating advanced statistics plots...")
    figures <- list()

    if (!is.null(results$effect_sizes))
        figures$forest <- tryCatch(plot_effect_size_forest(results$effect_sizes, plots_dir),
                                   error = function(e) NULL)
    if (!is.null(results$power_analysis))
        figures$power <- tryCatch(plot_power_curves(results$power_analysis, plots_dir),
                                  error = function(e) NULL)
    if (!is.null(results$mtc_comparison))
        figures$mtc <- tryCatch(plot_mtc_comparison(results$mtc_comparison, plots_dir),
                                error = function(e) NULL)
    if (!is.null(results$mixed_effects))
        figures$icc <- tryCatch(plot_icc_distribution(results$mixed_effects, plots_dir),
                                error = function(e) NULL)
    figures
}

plot_effect_size_forest <- function(effect_df, plots_dir) {
    sig_df <- effect_df[!is.na(effect_df$significant) & effect_df$significant &
                         !is.na(effect_df$effect_size), ]
    sig_df <- sig_df[order(-abs(sig_df$effect_size)), ]
    sig_df <- head(sig_df, 30)
    if (nrow(sig_df) == 0) return(NULL)

    sig_df$protein_id <- factor(sig_df$protein_id,
                                levels = sig_df$protein_id[order(sig_df$effect_size)])

    p <- ggplot2::ggplot(sig_df, ggplot2::aes(x = effect_size, y = protein_id)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
                                height = 0.3, color = "grey30") +
        ggplot2::geom_point(ggplot2::aes(color = effect_size > 0), size = 2) +
        ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                                    labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                                    name = "Direction") +
        ggplot2::labs(title = "Effect Sizes with 95% Bootstrap CI",
                      subtitle = paste0("Top ", nrow(sig_df), " significant proteins"),
                      x = "Effect Size", y = "Protein") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

    fig_path <- file.path(plots_dir, "effect_size_forest_plot.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 12)
    fig_path
}

plot_power_curves <- function(power_results, plots_dir) {
    power_df <- power_results$power_table

    power_long <- tidyr::pivot_longer(power_df, cols = tidyselect::starts_with("effect_"),
                                       names_to = "effect_size", values_to = "power")
    power_long$effect_size <- gsub("effect_", "", power_long$effect_size)
    power_long$effect_size <- ifelse(
        power_long$effect_size == "median",
        paste0("median (", round(power_results$median_effect, 2), ")"),
        power_long$effect_size
    )

    p1 <- ggplot2::ggplot(power_long,
               ggplot2::aes(x = n_per_group, y = power, color = effect_size)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
        ggplot2::geom_vline(xintercept = power_results$current_n,
                            linetype = "dotted", color = "blue") +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(title = "Statistical Power by Sample Size",
                      x = "Sample Size per Group", y = "Power",
                      color = "Effect Size (log2FC)") +
        ggplot2::theme_minimal()

    fig_path <- file.path(plots_dir, "power_analysis.png")
    ggplot2::ggsave(fig_path, p1, width = 10, height = 6)
    fig_path
}

plot_mtc_comparison <- function(mtc_df, plots_dir) {
    mtc_long <- tidyr::pivot_longer(mtc_df, cols = tidyselect::starts_with("n_significant"),
                                     names_to = "threshold", values_to = "n_significant")
    mtc_long$threshold <- gsub("n_significant_", "FDR < ", mtc_long$threshold)

    p <- ggplot2::ggplot(mtc_long,
             ggplot2::aes(x = method, y = n_significant, fill = threshold)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(title = "Multiple Testing Correction Comparison",
                      x = "Method", y = "Significant Proteins", fill = "Threshold") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    fig_path <- file.path(plots_dir, "mtc_comparison.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 6)
    fig_path
}

plot_icc_distribution <- function(mixed_df, plots_dir) {
    icc_values <- mixed_df$icc[mixed_df$converged & !is.na(mixed_df$icc)]
    if (length(icc_values) == 0) return(NULL)

    p <- ggplot2::ggplot(data.frame(icc = icc_values), ggplot2::aes(x = icc)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
        ggplot2::geom_vline(xintercept = median(icc_values), linetype = "dashed", color = "red") +
        ggplot2::labs(title = "Distribution of Intraclass Correlation Coefficients",
                      subtitle = "From mixed effects models", x = "ICC", y = "Count") +
        ggplot2::theme_minimal()

    fig_path <- file.path(plots_dir, "icc_distribution.png")
    ggplot2::ggsave(fig_path, p, width = 8, height = 6)
    fig_path
}

# -----------------------------------------------------------------------------
# Summary & Output
# -----------------------------------------------------------------------------

create_advanced_stats_summary <- function(results) {
    summary <- list()
    if (!is.null(results$effect_sizes)) {
        summary$n_significant_effects <- sum(results$effect_sizes$significant, na.rm = TRUE)
        summary$median_effect_size <- median(abs(results$effect_sizes$effect_size), na.rm = TRUE)
    }
    if (!is.null(results$mixed_effects)) {
        summary$n_converged_mixed <- sum(results$mixed_effects$converged, na.rm = TRUE)
        summary$median_icc <- median(results$mixed_effects$icc, na.rm = TRUE)
    }
    if (!is.null(results$power_analysis)) {
        pt <- results$power_analysis$power_table
        cn <- results$power_analysis$current_n
        row <- pt[pt$n_per_group == cn, ]
        if (nrow(row) > 0) summary$current_power <- row$effect_median[1]
    }
    summary
}

save_advanced_stats_outputs <- function(results, output_dir) {
    message("Saving advanced stats outputs...")
    saveRDS(results, file.path(output_dir, "advanced_stats_results.rds"))
}
