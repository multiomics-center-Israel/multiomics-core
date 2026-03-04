# =============================================================================
# Fallback Commentary Functions for Multi-Omics Pipeline Figures
# =============================================================================
# This module provides deterministic, data-driven commentary for pipeline figures
# when AI-powered commentary (Claude/OpenAI) is not available or fails
#
# All fallback functions follow a consistent structure:
#   - figure_id: Unique identifier
#   - title: Figure title
#   - what_is_this: Educational explanation
#   - observations: List of data-driven observations
#   - issues_checks: List of quality checks/potential issues
#   - next_steps: List of recommended follow-up actions
#   - confidence: Assessment confidence level
#   - limitations: Known limitations of the commentary


#' Generate fallback data-driven commentary (no LLM)
#'
#' Routes to appropriate fallback function based on figure type
#'
#' @param fig Figure row from figures_tbl
#' @param context Context information
#' @param mae_data Multi-assay experiment data (legacy format)
#' @param integration_results Integration results
#' @param concordance_results Concordance results
#' @param config Configuration
#' @return List with commentary structure
generate_fallback_commentary <- function(fig, context, mae_data, integration_results, concordance_results, config) {
    figure_id <- fig$figure_id

    # Route to specific fallback function based on figure type
    if (grepl("^mofa_variance", figure_id)) {
        return(fallback_mofa_variance(fig, context, integration_results))
    } else if (grepl("^mofa_factors", figure_id)) {
        return(fallback_mofa_factors(fig, context, integration_results))
    } else if (grepl("^mofa_weights", figure_id)) {
        return(fallback_mofa_weights(fig, context, integration_results))
    } else if (grepl("^diablo_sample", figure_id)) {
        return(fallback_diablo_sample(fig, context, integration_results))
    } else if (grepl("^diablo_loadings", figure_id)) {
        return(fallback_diablo_loadings(fig, context, integration_results))
    } else if (grepl("^diablo_circos", figure_id)) {
        return(fallback_diablo_circos(fig, context, integration_results))
    } else if (grepl("^diablo_cv", figure_id)) {
        return(fallback_diablo_cv(fig, context, integration_results))
    } else if (grepl("^snf_mds", figure_id)) {
        return(fallback_snf_mds(fig, context, integration_results))
    } else if (grepl("^snf_heatmap", figure_id)) {
        return(fallback_snf_heatmap(fig, context, integration_results))
    } else if (grepl("concordance", figure_id)) {
        return(fallback_concordance(fig, context, concordance_results))
    } else if (grepl("de_scatter", figure_id)) {
        return(fallback_de_scatter(fig, context, concordance_results))
    } else if (grepl("^ora_|^combined_enrichment", figure_id)) {
        return(fallback_enrichment(fig, context))
    } else {
        return(fallback_generic(fig, context))
    }
}


# =============================================================================
# MOFA2 Fallback Functions
# =============================================================================

fallback_mofa_variance <- function(fig, context, integration_results) {
    mofa <- integration_results$mofa

    observations <- list(
        "This heatmap displays the percentage of variance explained by each MOFA factor in each omics view.",
        paste0("The analysis included ", context$omics_present, " omics layers.")
    )

    if (!is.null(mofa$variance_explained)) {
        total_var <- sum(mofa$variance_explained$total)
        observations <- c(
            observations,
            paste0("Factors collectively explain approximately ", round(total_var, 1), "% of total variance.")
        )
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "MOFA2 variance explained heatmap showing how much variance each latent factor captures from each omics view. Factors capturing variance from multiple views represent shared biological signal.",
        observations = observations,
        issues_checks = list(
            "Check if some factors are dominated by a single omics view (may indicate technical variation)",
            "Factors explaining very little variance may be noise and could be excluded"
        ),
        next_steps = list(
            "Examine factor weights to identify key features",
            "Test factor associations with clinical variables",
            "Focus on factors with multi-view contributions for biological interpretation"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on summary statistics"
    )
}


fallback_mofa_factors <- function(fig, context, integration_results) {
    mofa <- integration_results$mofa

    observations <- list(
        "This scatter plot shows samples positioned according to their MOFA factor values.",
        paste0("Samples are colored by ", context$condition_col, "."),
        paste0("Analysis includes ", context$n_samples, " samples from ", context$omics_present, ".")
    )

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "MOFA2 factor plot showing samples in the latent factor space. Samples that cluster together share similar multi-omics profiles. Factor 1 typically captures the most variance.",
        observations = observations,
        issues_checks = list(
            "Check if sample groups separate along factor axes",
            "Look for outlier samples that don't cluster with their group",
            "Consider if separation is driven by biological signal or batch effects"
        ),
        next_steps = list(
            "Identify which factors correlate with condition",
            "Examine factor loadings to understand what drives separation",
            "Test factors for association with other clinical covariates"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary; visual patterns require manual inspection"
    )
}


fallback_mofa_weights <- function(fig, context, integration_results) {
    omic <- sub("mofa_weights_", "", fig$figure_id)

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = paste0("MOFA2 feature weights showing the contribution of individual ", omic, " features to each latent factor. High absolute weights indicate features that strongly define the factor."),
        observations = list(
            paste0("This plot shows top contributing features from ", omic, " to MOFA factors."),
            "Features with high positive or negative weights drive the factor pattern.",
            "Features with weights near zero contribute little to that factor."
        ),
        issues_checks = list(
            "Verify top features are biologically meaningful",
            "Check if housekeeping or technical features dominate (potential QC issue)"
        ),
        next_steps = list(
            "Perform pathway enrichment on high-weight features",
            "Cross-reference top features with differential analysis results",
            "Investigate features appearing in multiple important factors"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on figure metadata"
    )
}


# =============================================================================
# DIABLO Fallback Functions
# =============================================================================

fallback_diablo_sample <- function(fig, context, integration_results) {
    diablo <- integration_results$diablo

    observations <- list(
        "This plot shows samples in DIABLO component space, integrating all omics views.",
        paste0("Samples are colored by ", context$condition_col, "."),
        "DIABLO is a supervised method that maximizes separation between groups."
    )

    if (!is.null(context$diablo_cv_error)) {
        observations <- c(
            observations,
            paste0("Cross-validation error rate: ", round(context$diablo_cv_error * 100, 1), "%")
        )
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "DIABLO sample plot showing multi-omics integration results. Unlike PCA, DIABLO explicitly optimizes for group separation using features from all omics layers.",
        observations = observations,
        issues_checks = list(
            "Good separation suggests reliable multi-omics biomarker panel",
            "Overlapping groups may indicate weak signal or need for more samples",
            "Check for overfitting by examining CV error rates"
        ),
        next_steps = list(
            "Review selected features from each omics layer",
            "Examine the circos plot for cross-omics correlations",
            "Validate key features in an independent cohort"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary; classification performance requires CV metrics"
    )
}


fallback_diablo_loadings <- function(fig, context, integration_results) {
    omic <- sub("diablo_loadings_", "", fig$figure_id)

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = paste0("DIABLO loading plot showing selected ", omic, " features and their contribution to discriminant components. These are the features most important for multi-omics classification."),
        observations = list(
            paste0("Features from ", omic, " selected by DIABLO for optimal classification."),
            "Loading magnitude indicates feature importance for group discrimination.",
            "Features are selected through sparse regularization (LASSO-like)."
        ),
        issues_checks = list(
            "Verify selected features are biologically interpretable",
            "Check if key known markers are included in selection"
        ),
        next_steps = list(
            "Validate selected features with univariate tests",
            "Perform pathway analysis on selected feature set",
            "Consider these features for biomarker panel development"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on figure metadata"
    )
}


fallback_diablo_circos <- function(fig, context, integration_results) {
    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "DIABLO circos plot visualizing correlations between selected features across different omics layers. Lines connect features that co-vary, revealing cross-omics regulatory relationships.",
        observations = list(
            "Each arc represents a different omics layer.",
            "Lines connect features with strong correlations across omics.",
            "Dense connections suggest coordinated multi-omics patterns."
        ),
        issues_checks = list(
            "Look for biologically meaningful cross-omics connections",
            "Sparse connections may indicate independent omics patterns"
        ),
        next_steps = list(
            "Investigate the most connected features (hub features)",
            "Map connections to known biological pathways",
            "Consider connected feature pairs for mechanistic studies"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary; connection details require table inspection"
    )
}


fallback_diablo_cv <- function(fig, context, integration_results) {
    diablo <- integration_results$diablo

    observations <- list(
        "This plot shows classification error rate across cross-validation folds.",
        "Lower error rates indicate better classification performance.",
        "The optimal number of components balances accuracy and model complexity."
    )

    if (!is.null(context$diablo_cv_error)) {
        observations <- c(
            observations,
            paste0("Minimum CV error: ", round(context$diablo_cv_error * 100, 1), "%")
        )
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "DIABLO cross-validation error plot showing classification performance. This assesses how well the multi-omics signature generalizes to unseen samples.",
        observations = observations,
        issues_checks = list(
            "Error rate > 30% suggests weak discriminative signal",
            "Large variance across folds may indicate sample heterogeneity",
            "Overfitting suspected if training error << CV error"
        ),
        next_steps = list(
            "Select optimal component number based on error plateau",
            "Consider external validation on independent samples",
            "If error is high, revisit feature selection strategy"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on figure metadata"
    )
}


# =============================================================================
# SNF Fallback Functions
# =============================================================================

fallback_snf_mds <- function(fig, context, integration_results) {
    snf <- integration_results$snf

    observations <- list(
        "This MDS plot visualizes the fused similarity network from SNF.",
        "Samples are colored by their SNF cluster assignment.",
        "SNF integrates omics by fusing sample similarity networks."
    )

    if (!is.null(context$snf_nmi)) {
        observations <- c(
            observations,
            paste0("NMI with condition: ", round(context$snf_nmi, 3))
        )
    }

    if (!is.null(context$snf_n_clusters)) {
        observations <- c(
            observations,
            paste0("Number of clusters: ", context$snf_n_clusters)
        )
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "SNF MDS plot showing sample clustering based on fused multi-omics similarity. Unlike DIABLO, SNF is unsupervised and finds natural patient subgroups.",
        observations = observations,
        issues_checks = list(
            "Check if SNF clusters align with known condition groups",
            "Low NMI suggests clusters capture different structure than condition",
            "Look for potential novel patient subgroups"
        ),
        next_steps = list(
            "Compare SNF clusters with clinical outcomes",
            "Identify features that distinguish clusters",
            "Consider SNF for patient stratification"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary; clustering quality depends on NMI and silhouette"
    )
}


fallback_snf_heatmap <- function(fig, context, integration_results) {
    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "SNF fused network heatmap showing pairwise sample similarities after network fusion. High values (red) indicate similar multi-omics profiles. Block structure reveals sample clusters.",
        observations = list(
            "The heatmap shows the fused similarity matrix from SNF.",
            "Clear block structure indicates well-defined sample clusters.",
            "Off-diagonal patterns may reveal sub-structure within clusters."
        ),
        issues_checks = list(
            "Weak block structure suggests heterogeneous data",
            "Check if cluster blocks align with sample annotations"
        ),
        next_steps = list(
            "Examine cluster composition relative to condition",
            "Identify samples that bridge clusters",
            "Consider alternative number of clusters if structure is unclear"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on figure metadata"
    )
}


# =============================================================================
# Concordance Fallback Functions
# =============================================================================

fallback_concordance <- function(fig, context, concordance_results) {
    rp <- concordance_results$rna_protein

    observations <- list(
        "This histogram shows the distribution of RNA-protein correlation coefficients.",
        "Positive correlations indicate concordant expression patterns."
    )

    if (!is.null(context$concordance_mean_cor)) {
        observations <- c(
            observations,
            paste0("Mean correlation: ", round(context$concordance_mean_cor, 3)),
            paste0("Number of genes compared: ", context$concordance_n_genes),
            paste0("Percentage with positive correlation: ", round(context$concordance_pct_positive, 1), "%")
        )
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "RNA-protein concordance distribution showing how well mRNA expression predicts protein abundance. Positive correlations are expected, but post-transcriptional regulation causes deviations.",
        observations = observations,
        issues_checks = list(
            "Mean correlation < 0.3 suggests substantial post-transcriptional regulation",
            "Bimodal distribution may indicate two classes of genes",
            "Very low correlations could indicate sample mismatch or technical issues"
        ),
        next_steps = list(
            "Investigate genes with discordant RNA-protein patterns",
            "Check if discordant genes are enriched for specific functions",
            "Consider protein half-life and mRNA stability in interpretation"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on summary statistics"
    )
}


fallback_de_scatter <- function(fig, context, concordance_results) {
    observations <- list(
        "This scatter plot compares log2 fold changes between RNA and protein.",
        "Points on the diagonal show concordant direction of change.",
        "Off-diagonal points indicate discordant regulation."
    )

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = "RNA vs protein fold change scatter plot comparing differential expression results. Concordant features change in the same direction at both levels, suggesting transcriptional control.",
        observations = observations,
        issues_checks = list(
            "Low correlation suggests post-transcriptional regulation dominates",
            "Quadrant analysis: concordant (Q1, Q3) vs discordant (Q2, Q4)",
            "Outliers may represent interesting regulatory cases"
        ),
        next_steps = list(
            "Focus on concordant features for reliable biomarkers",
            "Investigate discordant features for post-transcriptional mechanisms",
            "Perform pathway analysis on concordant vs discordant sets"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary; requires visual inspection for patterns"
    )
}


# =============================================================================
# Enrichment Fallback Functions
# =============================================================================

fallback_enrichment <- function(fig, context) {
    omic <- sub("^ora_", "", fig$figure_id)

    if (fig$figure_id == "combined_enrichment") {
        return(list(
            figure_id = fig$figure_id,
            title = fig$title,
            what_is_this = "Combined multi-omics enrichment showing pathways significant across multiple omics layers. These represent convergent biological signals supported by independent molecular evidence.",
            observations = list(
                "Pathways enriched in multiple omics provide stronger biological evidence.",
                "Combined p-values integrate evidence across omics layers.",
                "Top pathways represent consensus biological signals."
            ),
            issues_checks = list(
                "Check overlap between pathway gene sets and measured features",
                "Pathways enriched in only one omics may still be important"
            ),
            next_steps = list(
                "Prioritize pathways with multi-omics support for follow-up",
                "Examine individual omics contributions to top pathways",
                "Integrate pathway findings with integration method results"
            ),
            confidence = "medium",
            limitations = "Data-driven commentary based on figure metadata"
        ))
    }

    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = paste0("Over-representation analysis (ORA) results for ", omic, " showing enriched biological pathways. Dot size indicates gene count, position shows significance."),
        observations = list(
            paste0("Enrichment analysis for ", omic, " differential features."),
            "Pathways are ranked by adjusted p-value significance.",
            "Larger dots indicate more features in the pathway."
        ),
        issues_checks = list(
            "Check for redundant pathway terms (GO hierarchy)",
            "Very broad terms may be less informative than specific ones"
        ),
        next_steps = list(
            "Compare enriched pathways across omics layers",
            "Focus on specific, interpretable pathways",
            "Examine leading-edge genes in top pathways"
        ),
        confidence = "medium",
        limitations = "Data-driven commentary based on figure metadata"
    )
}


# =============================================================================
# Generic Fallback Function
# =============================================================================

fallback_generic <- function(fig, context) {
    list(
        figure_id = fig$figure_id,
        title = fig$title,
        what_is_this = fig$description,
        observations = list(
            paste0("This is a ", fig$plot_type, " plot from the ", fig$section, " section."),
            paste0("Analysis includes ", context$n_samples, " samples and ", context$omics_present, " omics.")
        ),
        issues_checks = list(
            "Review the figure for expected patterns",
            "Check consistency with other analysis results"
        ),
        next_steps = list(
            "Cross-reference with related figures and tables",
            "Consider statistical significance of observed patterns"
        ),
        confidence = "low",
        limitations = "Generic data-driven commentary"
    )
}
