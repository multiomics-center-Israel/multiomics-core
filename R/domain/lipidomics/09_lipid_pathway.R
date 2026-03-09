# R/domain/lipidomics/09_lipid_pathway.R
#
# Lipid pathway and network analysis inspired by LipidOne's
# "Lipid System Biology" module (Nguyen et al. 2017).
#
# Part 1: Lipid metabolic pathway scoring (product/reactant ratios, Z-scores)
# Part 2: Gene/enzyme prediction and STRING network queries
# Part 3: Pathway interaction table formatting


# ==== PART 1: LIPID METABOLIC PATHWAY ANALYSIS ===============================

#' Lipid metabolic reaction database
#'
#' Returns a data.frame of known lipid metabolic transformations in human lipid
#' metabolism.  Each row describes a substrate -> product conversion catalyzed
#' by a specific enzyme or enzyme family.
#'
#' @param organism Character; currently only "human" is supported.
#' @return data.frame with columns: reactant_class, product_class, enzyme,
#'   reaction_type, description, gene_symbol
lipid_reaction_database <- function(organism = "human") {

    if (organism != "human") {
        warning("lipid_reaction_database: only 'human' currently supported; ",
                "returning human reactions")
    }

    reactions <- data.frame(
        reactant_class = c(
            # --- Glycerophospholipid metabolism ---
            "PC",   "PE",   "PS",   "PI",   "PG",          # PLA2 hydrolysis
            "PC",   "PE",                                   # PLD hydrolysis
            "PC",   "PE",                                   # PLC hydrolysis -> DG
            "PA",                                           # PAP/Lipin
            "DG",                                           # DGAT
            "DG",                                           # DGK
            "MG",                                           # MGAT
            "LPC",                                          # LPCAT re-acylation
            "LPE",                                          # LPEAT re-acylation
            "PE",                                           # PEMT
            "PC",                                           # PEMT reverse
            "DG",                                           # CPT
            "DG",                                           # EPT
            # --- Sphingolipid metabolism ---
            "Cer",  "SM",   "Cer",  "HexCer", "Hex2Cer",
            "Sph",  "Cer",  "Cer",
            "dhCer",                                        # DES
            "So",                                           # dhCer synthase (CerS on sphinganine)
            # --- PI signaling ---
            "PI",   "PIP",
            # --- Fatty acid / carnitine ---
            "FA",
            # --- Ether lipids ---
            "PC O", "PE O",
            # --- Cholesterol ester metabolism ---
            "CE"
        ),
        product_class = c(
            # PLA2 hydrolysis products
            "LPC",  "LPE",  "LPS",  "LPI",  "LPG",
            # PLD products
            "PA",   "PA",
            # PLC products
            "DG",   "DG",
            # PAP/Lipin
            "DG",
            # DGAT
            "TG",
            # DGK
            "PA",
            # MGAT
            "DG",
            # LPCAT
            "PC",
            # LPEAT
            "PE",
            # PEMT
            "PC",
            # PEMT reverse
            "PE",
            # CPT
            "PC",
            # EPT
            "PE",
            # Sphingolipid
            "SM",   "Cer",  "HexCer", "Hex2Cer", "GM3",
            "Cer",  "Sph",  "S1P",
            # DES
            "Cer",
            # dhCerS
            "dhCer",
            # PI signaling
            "PIP",  "PIP2",
            # CPT1
            "CAR",
            # Ether lipid PLA2
            "LPC O", "LPE O",
            # LCAT / CEH
            "FA"
        ),
        enzyme = c(
            "PLA2", "PLA2", "PLA2", "PLA2", "PLA2",
            "PLD",  "PLD",
            "PLC",  "PLC",
            "PAP/Lipin",
            "DGAT",
            "DGK",
            "MGAT",
            "LPCAT",
            "LPEAT",
            "PEMT",
            "PEMT_rev",
            "CPT",
            "EPT",
            "SMS",  "SMase", "GCS",  "LacCerS", "GM3S",
            "CerS", "CDase", "SphK",
            "DES",
            "CerS_dh",
            "PI4K", "PIP5K",
            "CPT1",
            "PLA2_ether", "PLA2_ether",
            "CEH"
        ),
        reaction_type = c(
            "hydrolysis", "hydrolysis", "hydrolysis", "hydrolysis", "hydrolysis",
            "hydrolysis", "hydrolysis",
            "hydrolysis", "hydrolysis",
            "dephosphorylation",
            "acylation",
            "phosphorylation",
            "acylation",
            "acylation",
            "acylation",
            "methylation",
            "headgroup_exchange",
            "transfer",
            "transfer",
            "transfer",   "hydrolysis", "glycosylation", "glycosylation", "glycosylation",
            "acylation",  "hydrolysis", "phosphorylation",
            "desaturation",
            "acylation",
            "phosphorylation", "phosphorylation",
            "conjugation",
            "hydrolysis", "hydrolysis",
            "hydrolysis"
        ),
        description = c(
            "Phospholipase A2 cleaves sn-2 acyl chain from PC",
            "Phospholipase A2 cleaves sn-2 acyl chain from PE",
            "Phospholipase A2 cleaves sn-2 acyl chain from PS",
            "Phospholipase A2 cleaves sn-2 acyl chain from PI",
            "Phospholipase A2 cleaves sn-2 acyl chain from PG",
            "Phospholipase D removes headgroup from PC",
            "Phospholipase D removes headgroup from PE",
            "Phospholipase C cleaves PC headgroup to DG",
            "Phospholipase C cleaves PE headgroup to DG",
            "PA phosphatase / Lipin dephosphorylates PA to DG",
            "Diacylglycerol acyltransferase adds acyl chain",
            "Diacylglycerol kinase phosphorylates DG",
            "Monoacylglycerol acyltransferase",
            "Lysophosphatidylcholine acyltransferase (Lands cycle)",
            "Lysophosphatidylethanolamine acyltransferase",
            "PE N-methyltransferase converts PE to PC",
            "Reverse of PEMT pathway",
            "CDP-choline:DAG cholinephosphotransferase",
            "CDP-ethanolamine:DAG ethanolaminephosphotransferase",
            "Sphingomyelin synthase transfers phosphocholine to Cer",
            "Sphingomyelinase hydrolyzes SM to Cer",
            "Glucosylceramide synthase",
            "Lactosylceramide synthase",
            "GM3 synthase (sialyltransferase)",
            "Ceramide synthase acylates sphingosine",
            "Ceramidase hydrolyzes Cer to sphingosine",
            "Sphingosine kinase phosphorylates Sph to S1P",
            "Dihydroceramide desaturase",
            "Ceramide synthase on sphinganine (de novo)",
            "PI 4-kinase",
            "PIP 5-kinase",
            "Carnitine palmitoyltransferase 1",
            "PLA2 on ether-linked PC",
            "PLA2 on ether-linked PE",
            "Cholesterol ester hydrolase"
        ),
        gene_symbol = c(
            "PLA2G2A",  "PLA2G2A",  "PLA2G2A",  "PLA2G2A",  "PLA2G2A",
            "PLD1",     "PLD1",
            "PLCG1",    "PLCG1",
            "LPIN1",
            "DGAT1",
            "DGKA",
            "MOGAT1",
            "LPCAT1",
            "LPEAT1",
            "PEMT",
            "PEMT",
            "CHPT1",
            "SELENOI",
            "SGMS1",    "SMPD1",    "UGCG",     "B4GALT5",  "ST3GAL5",
            "CERS2",    "ASAH1",    "SPHK1",
            "DEGS1",
            "CERS2",
            "PI4KA",    "PIP5K1A",
            "CPT1A",
            "PLA2G7",   "PLA2G7",
            "LIPE"
        ),
        stringsAsFactors = FALSE
    )

    reactions
}


#' Compute lipid pathway scores via product/reactant ratios
#'
#' For each reaction in the lipid reaction database, finds matching species in
#' the data, computes the mean product/reactant concentration ratio per sample,
#' then compares ratios between experimental groups with a t-test.  P-values are
#' converted to Z-scores (Nguyen et al. 2017 approach).
#'
#' @param pre     List from preprocess_lipidomics() (pre contract).
#' @param config  Full pipeline config.
#' @param de_res  Optional DE result list (not used directly but reserved).
#' @return data.frame with columns: reaction, reactant_class, product_class,
#'   enzyme, z_score, p_value, direction, n_reactant_species, n_product_species,
#'   gene_symbol
compute_lipid_pathway_scores <- function(pre, config, de_res = NULL) {

    cfg <- config$modes$lipidomics
    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples     %||% "sample_id"

    # Use raw filtered intensities (not log-transformed)
    mat      <- pre$expr_filt
    meta     <- pre$meta
    row_data <- pre$row_data

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- as.character(meta[[condition_col]])

    groups <- unique(condition)
    if (length(groups) != 2) {
        message("lipid pathway scores: requires exactly 2 groups, found ",
                length(groups), " — skipping")
        return(NULL)
    }

    # Identify group indices
    g1_idx <- which(condition == groups[1])
    g2_idx <- which(condition == groups[2])

    # Lipid classes present in data
    classes_in_data <- as.character(row_data$lipid_class)

    # Reaction database
    rxn_db <- lipid_reaction_database()

    results <- list()

    for (i in seq_len(nrow(rxn_db))) {
        rc <- rxn_db$reactant_class[i]
        pc <- rxn_db$product_class[i]

        reactant_idx <- which(classes_in_data == rc)
        product_idx  <- which(classes_in_data == pc)

        # Need at least 1 species on each side
        if (length(reactant_idx) == 0 || length(product_idx) == 0) next

        # Mean reactant intensity per sample
        if (length(reactant_idx) == 1) {
            reactant_means <- mat[reactant_idx, ]
        } else {
            reactant_means <- colMeans(mat[reactant_idx, , drop = FALSE],
                                       na.rm = TRUE)
        }

        # Mean product intensity per sample
        if (length(product_idx) == 1) {
            product_means <- mat[product_idx, ]
        } else {
            product_means <- colMeans(mat[product_idx, , drop = FALSE],
                                      na.rm = TRUE)
        }

        # Avoid division by zero
        reactant_means[reactant_means == 0] <- NA

        # Product/reactant ratio per sample
        ratios <- product_means / reactant_means

        # Split by group
        r_g1 <- ratios[g1_idx]
        r_g2 <- ratios[g2_idx]

        # Remove NAs
        r_g1 <- r_g1[!is.na(r_g1) & is.finite(r_g1)]
        r_g2 <- r_g2[!is.na(r_g2) & is.finite(r_g2)]

        # Need at least 2 samples per group for t-test
        if (length(r_g1) < 2 || length(r_g2) < 2) next

        tt <- tryCatch(
            stats::t.test(r_g2, r_g1),
            error = function(e) NULL
        )
        if (is.null(tt)) next

        p_val <- tt$p.value
        # Direction: positive z = reaction activated in group2 vs group1
        mean_diff <- mean(r_g2, na.rm = TRUE) - mean(r_g1, na.rm = TRUE)
        sign_dir <- sign(mean_diff)
        if (sign_dir == 0) sign_dir <- 1

        # Convert p-value to Z-score (two-sided)
        z_score <- sign_dir * stats::qnorm(1 - p_val / 2)
        if (!is.finite(z_score)) z_score <- sign_dir * 6  # cap at 6

        direction <- if (mean_diff > 0) "activated" else "deactivated"

        results[[length(results) + 1]] <- data.frame(
            reaction           = paste0(rc, " -> ", pc, " (", rxn_db$enzyme[i], ")"),
            reactant_class     = rc,
            product_class      = pc,
            enzyme             = rxn_db$enzyme[i],
            reaction_type      = rxn_db$reaction_type[i],
            z_score            = round(z_score, 3),
            p_value            = p_val,
            direction          = direction,
            n_reactant_species = length(reactant_idx),
            n_product_species  = length(product_idx),
            gene_symbol        = rxn_db$gene_symbol[i],
            group1             = groups[1],
            group2             = groups[2],
            stringsAsFactors   = FALSE
        )
    }

    if (length(results) == 0) {
        message("lipid pathway scores: no testable reactions found")
        return(NULL)
    }

    result_df <- do.call(rbind, results)
    result_df$p_adjusted <- stats::p.adjust(result_df$p_value, method = "BH")
    result_df <- result_df[order(result_df$p_value), ]
    rownames(result_df) <- NULL

    result_df
}


#' Plot lipid pathway network
#'
#' Network visualization where nodes are lipid classes and edges are metabolic
#' reactions.  Edge color encodes direction (activated vs deactivated) and edge
#' thickness encodes significance.
#'
#' @param pathway_scores data.frame from compute_lipid_pathway_scores().
#' @param p_threshold    Significance threshold for displaying edges.
#' @return ggplot object or NULL.
plot_lipid_pathway_network <- function(pathway_scores, p_threshold = 0.05) {

    if (is.null(pathway_scores) || nrow(pathway_scores) == 0) return(NULL)

    sig <- pathway_scores[pathway_scores$p_adjusted <= p_threshold, , drop = FALSE]
    if (nrow(sig) == 0) {
        # Fall back to top 10 by p-value
        sig <- utils::head(pathway_scores[order(pathway_scores$p_value), ], 10)
    }

    # Lipid class categories for node coloring
    class_category <- c(
        PC = "GP", PE = "GP", PS = "GP", PI = "GP", PG = "GP", PA = "GP",
        LPC = "GP", LPE = "GP", LPS = "GP", LPI = "GP", LPG = "GP",
        `PC O` = "GP", `PE O` = "GP", `LPC O` = "GP", `LPE O` = "GP",
        PIP = "GP", PIP2 = "GP",
        DG = "GL", TG = "GL", MG = "GL",
        Cer = "SP", SM = "SP", HexCer = "SP", Hex2Cer = "SP", GM3 = "SP",
        Sph = "SP", S1P = "SP", dhCer = "SP", So = "SP",
        FA = "FA", CAR = "FA",
        CE = "ST"
    )

    category_colors <- c(
        GP = "#3182BD", GL = "#31A354", SP = "#E6550D",
        FA = "#756BB1", ST = "#636363"
    )

    # Build node list
    all_nodes <- unique(c(sig$reactant_class, sig$product_class))
    node_cats <- vapply(all_nodes, function(n) {
        if (n %in% names(class_category)) class_category[[n]] else "Other"
    }, character(1))

    # Map abbreviations to full names for display
    dict <- lipid_class_dictionary()
    abbr_to_full <- stats::setNames(dict$full_name, dict$abbreviation)
    node_labels <- vapply(all_nodes, function(n) {
        full <- abbr_to_full[n]
        if (!is.na(full)) full else n
    }, character(1))

    use_igraph <- requireNamespace("igraph", quietly = TRUE)

    if (use_igraph) {
        # Build igraph for layout
        edge_df <- data.frame(
            from = sig$reactant_class,
            to   = sig$product_class,
            stringsAsFactors = FALSE
        )
        g <- igraph::graph_from_data_frame(edge_df, directed = TRUE,
                                            vertices = data.frame(name = all_nodes))
        layout <- igraph::layout_with_fr(g)
        node_x <- layout[, 1]
        node_y <- layout[, 2]
    } else {
        # Circular layout fallback
        n_nodes <- length(all_nodes)
        angles  <- seq(0, 2 * pi, length.out = n_nodes + 1)[seq_len(n_nodes)]
        node_x  <- cos(angles)
        node_y  <- sin(angles)
    }

    node_df <- data.frame(
        label    = all_nodes,
        display  = unname(node_labels),
        x        = node_x,
        y        = node_y,
        category = node_cats,
        stringsAsFactors = FALSE
    )

    # Build edge segments
    edge_list <- lapply(seq_len(nrow(sig)), function(i) {
        from_idx <- which(node_df$label == sig$reactant_class[i])
        to_idx   <- which(node_df$label == sig$product_class[i])
        data.frame(
            x     = node_df$x[from_idx],
            y     = node_df$y[from_idx],
            xend  = node_df$x[to_idx],
            yend  = node_df$y[to_idx],
            direction = sig$direction[i],
            thickness = -log10(pmax(sig$p_value[i], 1e-300)),
            enzyme    = sig$enzyme[i],
            stringsAsFactors = FALSE
        )
    })
    edge_df <- do.call(rbind, edge_list)

    # Compute label position offset (slightly outside node)
    if (use_igraph) {
        center_x <- mean(node_df$x)
        center_y <- mean(node_df$y)
        dx <- node_df$x - center_x
        dy <- node_df$y - center_y
        dist <- sqrt(dx^2 + dy^2)
        dist[dist == 0] <- 1
        label_offset <- 0.12 * max(dist)
        node_df$label_x <- node_df$x + label_offset * dx / dist
        node_df$label_y <- node_df$y + label_offset * dy / dist
    } else {
        node_df$label_x <- node_df$x * 1.15
        node_df$label_y <- node_df$y * 1.15
    }

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = edge_df,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                         color = direction, linewidth = thickness),
            arrow = grid::arrow(length = grid::unit(0.15, "inches"),
                                type = "closed"),
            alpha = 0.7
        ) +
        ggplot2::geom_point(
            data = node_df,
            ggplot2::aes(x = x, y = y, fill = category),
            shape = 21, size = 7, color = "black", stroke = 0.8
        ) +
        ggplot2::geom_text(
            data = node_df,
            ggplot2::aes(x = label_x, y = label_y, label = display),
            size = 3.0, fontface = "bold"
        ) +
        ggplot2::scale_color_manual(
            values = c(activated = "#2166AC", deactivated = "#B2182B"),
            name   = "Reaction Direction"
        ) +
        ggplot2::scale_fill_manual(
            values = category_colors,
            name   = "Lipid Category"
        ) +
        ggplot2::scale_linewidth_continuous(range = c(0.5, 3), name = "-log10(p)") +
        ggplot2::labs(
            title = "Lipid Metabolic Pathway Network"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 14,
                                                   hjust = 0.5),
            legend.position = "right"
        )

    p
}


#' Plot lipid pathway Z-score barplot
#'
#' Horizontal barplot of Z-scores for significant metabolic reactions.
#'
#' @param pathway_scores data.frame from compute_lipid_pathway_scores().
#' @param p_threshold    Significance threshold for filtering reactions.
#' @return ggplot object or NULL.
plot_lipid_pathway_barplot <- function(pathway_scores, p_threshold = 0.05) {

    if (is.null(pathway_scores) || nrow(pathway_scores) == 0) return(NULL)

    sig <- pathway_scores[pathway_scores$p_adjusted <= p_threshold, , drop = FALSE]
    if (nrow(sig) == 0) {
        message("lipid pathway barplot: no significant reactions at p < ",
                p_threshold, " — showing top 10")
        sig <- utils::head(pathway_scores[order(pathway_scores$p_value), ], 10)
    }

    sig <- sig[order(sig$z_score), ]

    # Create display labels with full class names
    dict <- lipid_class_dictionary()
    abbr_to_full <- stats::setNames(dict$full_name, dict$abbreviation)
    sig$reaction_display <- vapply(seq_len(nrow(sig)), function(i) {
        rc_full <- abbr_to_full[sig$reactant_class[i]]
        pc_full <- abbr_to_full[sig$product_class[i]]
        if (is.na(rc_full)) rc_full <- sig$reactant_class[i]
        if (is.na(pc_full)) pc_full <- sig$product_class[i]
        paste0(rc_full, " -> ", pc_full, " (", sig$enzyme[i], ")")
    }, character(1))
    sig$reaction_display <- factor(sig$reaction_display, levels = sig$reaction_display)

    ggplot2::ggplot(sig, ggplot2::aes(x = z_score, y = reaction_display,
                                       fill = direction)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
        ggplot2::scale_fill_manual(
            values = c(activated = "#2166AC", deactivated = "#B2182B"),
            name   = "Direction"
        ) +
        ggplot2::labs(
            title = "Lipid Metabolic Pathway Activity (Z-scores)",
            x     = "Z-score",
            y     = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13,
                                                hjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 9)
        )
}


# ==== PART 2: GENE / ENZYME PREDICTION =======================================

#' Predict lipid-metabolizing enzymes from pathway scores
#'
#' Extracts significantly altered enzymatic reactions and returns a gene-level
#' summary of predicted enzyme activity changes.
#'
#' @param pathway_scores data.frame from compute_lipid_pathway_scores().
#' @param p_threshold    Significance threshold.
#' @return data.frame with columns: gene_symbol, enzyme, reaction, z_score,
#'   p_value, p_adjusted, activation_state
predict_lipid_enzymes <- function(pathway_scores, p_threshold = 0.05) {

    if (is.null(pathway_scores) || nrow(pathway_scores) == 0) return(NULL)

    sig <- pathway_scores[pathway_scores$p_adjusted <= p_threshold, , drop = FALSE]
    if (nrow(sig) == 0) {
        message("predict_lipid_enzymes: no significant reactions at p < ",
                p_threshold)
        return(NULL)
    }

    enzyme_df <- data.frame(
        gene_symbol      = sig$gene_symbol,
        enzyme           = sig$enzyme,
        reaction         = sig$reaction,
        z_score          = sig$z_score,
        p_value          = sig$p_value,
        p_adjusted       = sig$p_adjusted,
        activation_state = sig$direction,
        stringsAsFactors = FALSE
    )

    # If a gene appears in multiple reactions, keep all but sort by significance
    enzyme_df <- enzyme_df[order(enzyme_df$p_value), ]
    rownames(enzyme_df) <- NULL

    enzyme_df
}


#' Query STRING database for enzyme interactions
#'
#' If the STRINGdb package is available, queries the STRING protein-protein
#' interaction network for the predicted enzymes and retrieves functional
#' enrichment.
#'
#' @param enzyme_list    Character vector of gene symbols.
#' @param organism       NCBI taxonomy ID (default 9606 = human).
#' @param score_threshold Minimum interaction score (0-1000, default 400).
#' @return list(network_df, enrichment_df) or NULL if STRINGdb not available.
run_enzyme_string_network <- function(enzyme_list,
                                      organism = 9606,
                                      score_threshold = 400) {

    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        message("run_enzyme_string_network: STRINGdb not available — skipping")
        return(NULL)
    }

    if (length(enzyme_list) < 2) {
        message("run_enzyme_string_network: need >= 2 enzymes — skipping")
        return(NULL)
    }

    result <- tryCatch({
        string_db <- STRINGdb::STRINGdb$new(
            version    = "11.5",
            species    = organism,
            score_threshold = score_threshold
        )

        # Map gene symbols to STRING IDs
        gene_df <- data.frame(gene = unique(enzyme_list),
                              stringsAsFactors = FALSE)
        mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

        if (is.null(mapped) || nrow(mapped) == 0) {
            message("run_enzyme_string_network: no genes mapped to STRING")
            return(NULL)
        }

        # Get interactions
        interactions <- string_db$get_interactions(mapped$STRING_id)

        network_df <- NULL
        if (!is.null(interactions) && nrow(interactions) > 0) {
            network_df <- data.frame(
                protein1 = interactions$from,
                protein2 = interactions$to,
                score    = interactions$combined_score,
                stringsAsFactors = FALSE
            )
        }

        # Functional enrichment
        enrichment_df <- tryCatch({
            enrich <- string_db$get_enrichment(mapped$STRING_id)
            if (!is.null(enrich) && nrow(enrich) > 0) enrich else NULL
        }, error = function(e) {
            message("STRING enrichment failed: ", conditionMessage(e))
            NULL
        })

        list(
            network_df    = network_df,
            enrichment_df = enrichment_df,
            mapped_genes  = mapped
        )
    }, error = function(e) {
        message("run_enzyme_string_network failed: ", conditionMessage(e))
        NULL
    })

    result
}


#' Plot enzyme interaction network
#'
#' Network visualization of predicted lipid-metabolizing enzymes with STRING
#' interactions (if available).  Nodes are colored by activation state.
#'
#' @param network_result    List from run_enzyme_string_network() (may be NULL).
#' @param enzyme_predictions data.frame from predict_lipid_enzymes().
#' @return ggplot object or NULL.
plot_enzyme_network <- function(network_result, enzyme_predictions) {

    if (is.null(enzyme_predictions) || nrow(enzyme_predictions) == 0) {
        return(NULL)
    }

    # Unique enzymes with their activation state (take strongest z-score)
    genes <- unique(enzyme_predictions$gene_symbol)
    gene_state <- vapply(genes, function(g) {
        sub_df <- enzyme_predictions[enzyme_predictions$gene_symbol == g, ]
        best <- sub_df[which.min(sub_df$p_value), ]
        best$activation_state
    }, character(1))

    gene_z <- vapply(genes, function(g) {
        sub_df <- enzyme_predictions[enzyme_predictions$gene_symbol == g, ]
        best <- sub_df[which.min(sub_df$p_value), ]
        best$z_score
    }, numeric(1))

    # Attempt igraph-based layout if we have network edges
    has_igraph <- requireNamespace("igraph", quietly = TRUE)
    has_edges  <- !is.null(network_result) && !is.null(network_result$network_df) &&
                  nrow(network_result$network_df) > 0

    if (has_igraph && has_edges) {
        # Map STRING IDs back to gene symbols for edges
        mapped <- network_result$mapped_genes
        id_to_gene <- stats::setNames(mapped$gene, mapped$STRING_id)

        net_df <- network_result$network_df
        net_df$gene1 <- id_to_gene[net_df$protein1]
        net_df$gene2 <- id_to_gene[net_df$protein2]
        net_df <- net_df[!is.na(net_df$gene1) & !is.na(net_df$gene2), ]

        if (nrow(net_df) > 0) {
            g <- igraph::graph_from_data_frame(
                net_df[, c("gene1", "gene2")],
                directed = FALSE,
                vertices = data.frame(name = genes)
            )
            layout <- igraph::layout_with_fr(g)
        } else {
            n <- length(genes)
            angles <- seq(0, 2 * pi, length.out = n + 1)[seq_len(n)]
            layout <- cbind(cos(angles), sin(angles))
            net_df <- NULL
        }
    } else {
        n <- length(genes)
        angles <- seq(0, 2 * pi, length.out = n + 1)[seq_len(n)]
        layout <- cbind(cos(angles), sin(angles))
        net_df <- NULL
    }

    node_df <- data.frame(
        gene     = genes,
        x        = layout[, 1],
        y        = layout[, 2],
        state    = gene_state,
        z_score  = gene_z,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot()

    # Draw edges if available
    if (!is.null(net_df) && nrow(net_df) > 0) {
        edge_segments <- lapply(seq_len(nrow(net_df)), function(i) {
            i1 <- which(node_df$gene == net_df$gene1[i])
            i2 <- which(node_df$gene == net_df$gene2[i])
            if (length(i1) == 0 || length(i2) == 0) return(NULL)
            data.frame(
                x = node_df$x[i1], y = node_df$y[i1],
                xend = node_df$x[i2], yend = node_df$y[i2],
                score = net_df$score[i],
                stringsAsFactors = FALSE
            )
        })
        edge_df <- do.call(rbind, edge_segments)

        if (!is.null(edge_df) && nrow(edge_df) > 0) {
            p <- p + ggplot2::geom_segment(
                data = edge_df,
                ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                             linewidth = score / 1000),
                color = "grey60", alpha = 0.5
            ) +
            ggplot2::scale_linewidth_continuous(range = c(0.3, 2),
                                                name = "Interaction\nScore")
        }
    }

    p <- p +
        ggplot2::geom_point(
            data = node_df,
            ggplot2::aes(x = x, y = y, fill = state),
            shape = 21, size = 8, color = "black", stroke = 0.8
        ) +
        ggplot2::geom_text(
            data = node_df,
            ggplot2::aes(x = x, y = y, label = gene),
            size = 2.8, fontface = "bold"
        ) +
        ggplot2::scale_fill_manual(
            values = c(activated = "#2166AC", deactivated = "#B2182B"),
            name   = "Activation"
        ) +
        ggplot2::labs(title = "Predicted Lipid Enzyme Network") +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14,
                                                hjust = 0.5),
            legend.position = "right"
        )

    p
}


# ==== PART 3: PATHWAY INTERACTION TABLE ======================================

#' Format pathway scores into a display-ready table
#'
#' @param pathway_scores data.frame from compute_lipid_pathway_scores().
#' @return data.frame with columns: From, To, Enzyme, Gene, Z_score, P_value,
#'   FDR, Direction, N_Reactants, N_Products
format_pathway_table <- function(pathway_scores) {

    if (is.null(pathway_scores) || nrow(pathway_scores) == 0) return(NULL)

    dict <- lipid_class_dictionary()
    abbr_to_full <- stats::setNames(dict$full_name, dict$abbreviation)

    from_full <- unname(abbr_to_full[pathway_scores$reactant_class])
    from_full[is.na(from_full)] <- pathway_scores$reactant_class[is.na(from_full)]
    to_full <- unname(abbr_to_full[pathway_scores$product_class])
    to_full[is.na(to_full)] <- pathway_scores$product_class[is.na(to_full)]

    tbl <- data.frame(
        From        = from_full,
        To          = to_full,
        Enzyme      = pathway_scores$enzyme,
        Gene        = pathway_scores$gene_symbol,
        Z_score     = round(pathway_scores$z_score, 2),
        P_value     = signif(pathway_scores$p_value, 3),
        FDR         = signif(pathway_scores$p_adjusted, 3),
        Direction   = pathway_scores$direction,
        N_Reactants = pathway_scores$n_reactant_species,
        N_Products  = pathway_scores$n_product_species,
        stringsAsFactors = FALSE
    )

    tbl <- tbl[order(tbl$P_value), ]
    rownames(tbl) <- NULL

    tbl
}
