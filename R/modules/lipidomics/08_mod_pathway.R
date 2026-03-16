# R/modules/lipidomics/08_mod_pathway.R
#
# Lipidomics pathway and network analysis module.
# Wraps lipid pathway scoring, enzyme prediction, and STRING network.


#' Lipidomics pathway and network analysis module
#'
#' @param pre     List from preprocess_lipidomics().
#' @param de_res  List from mod_lipidomics_de().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(pathway_scores, enzyme_predictions, string_network, plots, files)
mod_lipidomics_pathway <- function(pre, de_res, config, out_dir) {
    assert_pre_contract(pre, stage = "lipidomics")

    cfg <- config$modes$lipidomics
    pw_cfg <- cfg$pathway %||% list()

    if (identical(pw_cfg$enabled, FALSE)) {
        message("lipidomics pathway: disabled in config -- skipping")
        return(NULL)
    }

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plots <- list()
    files <- character(0)

    p_threshold <- pw_cfg$p_threshold %||% 0.05
    organism    <- pw_cfg$organism    %||% "human"

    # ---- Lipid pathway scoring ----
    pathway_scores <- tryCatch(
        compute_lipid_pathway_scores(pre, config, de_res = de_res,
                                     organism = organism),
        error = function(e) {
            warning("Lipid pathway scoring failed: ", e$message)
            NULL
        }
    )

    if (!is.null(pathway_scores) && nrow(pathway_scores) > 0) {
        f_pw <- save_tsv(pathway_scores, out_ds, "lipid_pathway_scores.tsv")
        files <- c(files, f_pw)

        # Pathway interaction table
        pw_table <- format_pathway_table(pathway_scores)
        f_pw_table <- save_tsv(pw_table, out_ds, "lipid_pathway_interactions.tsv")
        files <- c(files, f_pw_table)

        # Pathway network plot
        f_net <- file.path(out_qc, "lipid_pathway_network.png")
        p_net <- tryCatch({
            pn <- plot_lipid_pathway_network(pathway_scores,
                                              p_threshold = p_threshold)
            ggplot2::ggsave(f_net, pn, width = 10, height = 9, dpi = 300)
            pn
        }, error = function(e) {
            warning("Pathway network plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_net)) {
            files <- c(files, f_net)
            plots$pathway_network <- p_net
        }

        # Pathway barplot
        f_bar <- file.path(out_qc, "lipid_pathway_barplot.png")
        p_bar <- tryCatch({
            pb <- plot_lipid_pathway_barplot(pathway_scores,
                                              p_threshold = p_threshold)
            ggplot2::ggsave(f_bar, pb, width = 9, height = 7, dpi = 300)
            pb
        }, error = function(e) {
            warning("Pathway barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_bar)) {
            files <- c(files, f_bar)
            plots$pathway_barplot <- p_bar
        }
    }

    # ---- Enzyme prediction ----
    enzyme_predictions <- NULL
    if (!is.null(pathway_scores)) {
        enzyme_predictions <- tryCatch(
            predict_lipid_enzymes(pathway_scores, p_threshold = p_threshold),
            error = function(e) {
                warning("Enzyme prediction failed: ", e$message)
                NULL
            }
        )

        if (!is.null(enzyme_predictions) && nrow(enzyme_predictions) > 0) {
            f_enz <- save_tsv(enzyme_predictions, out_ds,
                              "lipid_predicted_enzymes.tsv")
            files <- c(files, f_enz)
        }
    }

    # ---- STRING network (optional) ----
    string_network <- NULL
    if (!is.null(enzyme_predictions) &&
        nrow(enzyme_predictions) > 0 &&
        !identical(pw_cfg$string_network, FALSE)) {

        # Map organism name to STRING taxonomy ID if not explicitly set
        organism_taxid <- c(human = 9606L, mouse = 10090L, rat = 10116L)
        organism_lower <- tolower(organism)
        organism_key <- if (organism_lower %in% c("human", "homo sapiens")) {
            "human"
        } else if (organism_lower %in% c("mouse", "mus musculus")) {
            "mouse"
        } else if (organism_lower %in% c("rat", "rattus norvegicus")) {
            "rat"
        } else {
            "human"
        }
        default_taxid <- organism_taxid[[organism_key]]
        string_organism <- pw_cfg$string_organism %||% default_taxid
        score_threshold <- pw_cfg$string_score    %||% 400

        string_network <- tryCatch(
            run_enzyme_string_network(
                enzyme_list     = enzyme_predictions,
                organism        = string_organism,
                score_threshold = score_threshold
            ),
            error = function(e) {
                warning("STRING network failed: ", e$message)
                NULL
            }
        )

        if (!is.null(string_network)) {
            if (!is.null(string_network$network_df) &&
                nrow(string_network$network_df) > 0) {
                f_str_net <- save_tsv(string_network$network_df, out_ds,
                                      "string_protein_interactions.tsv")
                files <- c(files, f_str_net)
            }
            if (!is.null(string_network$enrichment_df) &&
                nrow(string_network$enrichment_df) > 0) {
                f_str_enr <- save_tsv(string_network$enrichment_df, out_ds,
                                      "string_protein_enrichment.tsv")
                files <- c(files, f_str_enr)
            }

            # Enzyme network plot
            f_enz_net <- file.path(out_qc, "enzyme_interaction_network.png")
            p_enz_net <- tryCatch({
                pe <- plot_enzyme_network(string_network, enzyme_predictions)
                ggplot2::ggsave(f_enz_net, pe, width = 10, height = 9, dpi = 300)
                pe
            }, error = function(e) {
                warning("Enzyme network plot failed: ", e$message)
                NULL
            })
            if (!is.null(p_enz_net)) {
                files <- c(files, f_enz_net)
                plots$enzyme_network <- p_enz_net
            }
        }
    }

    list(
        pathway_scores     = pathway_scores,
        enzyme_predictions = enzyme_predictions,
        string_network     = string_network,
        plots              = plots,
        files              = unique(files)
    )
}
