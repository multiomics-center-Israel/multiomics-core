# R/modules/metabolomics/04b_mod_network.R
#
# Metabolomics DE metabolite network module.
# Wraps domain network functions for pipeline integration.


#' Build DE metabolite network with KEGG reaction pair edges
#'
#' Extracts DE results and feature annotations, builds an interactive network
#' of DE metabolites connected by known KEGG reaction pairs, and saves
#' the visualization as a self-contained HTML file.
#'
#' @param de_res  DE results from mod_metabolomics_de().
#' @param pre     Preprocessing results (needs row_data for KEGG IDs).
#' @param config  Full pipeline config.
#' @param out_dir Output directory for this mode.
#' @return list with graph, nodes, edges, counts, and output file path.
#'   Returns NULL if network cannot be built (no KEGG IDs).
mod_metabolomics_network <- function(de_res, pre, config, out_dir) {
    metab_cfg <- config$modes$metabolomics

    # Extract DE table from mod_metabolomics_de() output
    # Structure: de_res$de_tables (named list of per-contrast data.frames)
    #            de_res$summary_df (full summary)
    de_table <- NULL
    if (!is.null(de_res$de_tables) && length(de_res$de_tables) > 0) {
        # Use first contrast table
        de_table <- de_res$de_tables[[1]]
    } else if (!is.null(de_res$summary_df)) {
        de_table <- de_res$summary_df
    } else if (is.data.frame(de_res)) {
        de_table <- de_res
    }

    if (is.null(de_table)) {
        message("metabolomics network: no DE results available — skipping")
        return(NULL)
    }

    # Feature annotations from row_data
    feature_ann <- pre$row_data
    if (is.null(feature_ann) || !"KEGG" %in% colnames(feature_ann)) {
        message("metabolomics network: no KEGG IDs in feature annotations — skipping")
        return(NULL)
    }

    p_cutoff <- metab_cfg$de$p_cutoff %||% 0.05

    # Build network
    network <- tryCatch(
        build_de_metabolite_network(
            de_res              = de_table,
            feature_annotations = feature_ann,
            p_cutoff            = p_cutoff,
            remove_isolated     = TRUE
        ),
        error = function(e) {
            warning("metabolomics network failed: ", e$message)
            NULL
        }
    )

    if (is.null(network)) return(NULL)

    # Save interactive HTML into Network/ subfolder
    net_dir <- file.path(out_dir, "Network")
    ensure_dir(net_dir)
    out_html <- file.path(net_dir, "de_metabolite_network.html")

    tryCatch(
        plot_metabolite_network_interactive(
            network_result = network,
            output_file    = out_html,
            title          = "DE Metabolite Network (KEGG Reaction Pairs)"
        ),
        error = function(e) {
            warning("Network HTML generation failed: ", e$message)
        }
    )

    # Save node/edge tables into Network/ subfolder
    if (!is.null(network$nodes) && nrow(network$nodes) > 0) {
        save_tsv(network$nodes, net_dir, "network_nodes.tsv")
    }
    if (!is.null(network$edges) && nrow(network$edges) > 0) {
        save_tsv(network$edges, net_dir, "network_edges.tsv")
    }

    network$html_file <- out_html
    network
}
