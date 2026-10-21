# KEGG pathway coverage per omics layer: why a cell in a KEGG figure is empty
#
# A KEGG pathway can be missing from one layer's results for quite different
# reasons, and the figures used to draw all of them as the same grey:
#
#   - the pathway has no members of that omics type at all in KEGG (a
#     signalling map has no compounds; a map with no genes for this organism);
#   - it has members, but none -- or fewer than the producer's minimum -- were
#     measured and mapped in this run, so the layer could not test it;
#   - the layer could test it and simply reported no result.
#
# The first two are "not expected to give a result" and are drawn as one mark,
# a grey cell with a cross; the third is white. The class and the counts behind
# each cell are kept in the tables.
#
# Metabolomics reads the cached KEGG compound-pathway table and the local
# HMDB-to-KEGG file, never the network. Gene and protein layers need the KEGG
# gene-pathway links of the organism: read from a cache under the output
# directory, downloaded once when that cache is cold. Nothing here re-tests a
# pathway or changes a p-value; it only labels what each layer could see.


#' Layer name of the metabolomics block
#' @keywords internal
.METABOLITE_LAYER <- "metabolomics"

#' Coverage classes meaning "this layer could not test the pathway"
#' @keywords internal
.NOT_TESTABLE_CLASSES <- c("no_kegg_members", "no_measured_members",
                           "below_min_measured", "above_max_measured",
                           "no_id_mapping")


#' Minimum measured compounds a KEGG compound pathway needs to be tested
#'
#' \code{min_set_size} is configured on a gene-set scale, where 10 is modest,
#' but a KEGG compound pathway rarely has ten measured members in one
#' experiment, so the compound producers do not take it as given. The floor is
#' three measured compounds, in every producer: compound ORA and compound GSEA
#' in the cross-omics enrichment and compound ORA in the Multi-ORA. Written
#' once so they and the coverage table cannot drift apart.
#'
#' @param min_gs Configured minimum set size. Accepted so callers can pass what
#'   they hold, and deliberately not used: the compound floor does not follow
#'   the gene-set scale.
#' @return Integer 3.
#' @examples
#' compound_min_measured(10)  # 3
#' compound_min_measured(2)   # 3
compound_min_measured <- function(min_gs = NULL) {
    3L
}


#' Minimum measured genes a KEGG pathway needs in the KEGG-from-DE enrichment
#'
#' \code{run_kegg_enrichment_for_omics()} relaxes the configured minimum for a
#' small mapped list (fewer than 500 features) to between 3 and 5. This is that
#' rule, shared by the producer and the coverage table.
#'
#' @param min_gs Configured minimum set size.
#' @param n_mapped Number of mapped features in the list being tested.
#' @return Integer minimum.
#' @examples
#' gene_min_measured(10, 200)   # 5
#' gene_min_measured(10, 2000)  # 10
gene_min_measured <- function(min_gs, n_mapped) {
    if (n_mapped < 500) as.integer(max(3, min(min_gs, 5))) else as.integer(min_gs)
}


#' Minimum measured members per layer in the Multi-ORA
#'
#' Multi-ORA tests each gene layer's significant genes against the pooled gene
#' universe (minimum set size 5 over that universe), so for one layer a result
#' needs only one of its own genes in the pathway. Compound ORA there uses the
#' shared compound floor, \code{compound_min_measured()}.
#'
#' @return Named list of integer minima.
#' @keywords internal
.multi_ora_min_measured <- function() {
    list(transcriptomics = 1L, proteomics = 1L,
         metabolomics = compound_min_measured())
}


#' Classify a pathway's coverage from its counts
#'
#' @param n_measured Integer vector of measured members.
#' @param min_measured Minimum for a test.
#' @param max_measured Maximum for a test.
#' @return Character vector of classes.
#' @keywords internal
.coverage_class <- function(n_measured, min_measured, max_measured = Inf) {
    ifelse(n_measured == 0, "no_measured_members",
    ifelse(n_measured < min_measured, "below_min_measured",
    ifelse(n_measured > max_measured, "above_max_measured", "testable")))
}


#' Per-pathway coverage of one layer's measured features
#'
#' One row per KEGG pathway in \code{memberships}, counting its distinct
#' members and how many were measured, and classifying it as
#' \code{testable}, \code{no_measured_members}, \code{below_min_measured} or
#' \code{above_max_measured}. A pathway with no member at all is not in
#' \code{memberships}; it is classed \code{no_kegg_members} where a result key
#' is looked up (\code{annotate_pathway_coverage()}).
#'
#' @param memberships Data frame with a \code{pathway} column, a member column
#'   and optionally \code{name}.
#' @param measured Character vector of measured member ids (same id space).
#' @param min_measured Minimum measured members for a test.
#' @param max_measured Maximum measured members for a test; \code{Inf} for none.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @param member_col Name of the member column.
#' @return Data frame: \code{norm_id}, \code{pathway_id}, \code{pathway_name},
#'   \code{n_in_pathway}, \code{n_measured}, \code{n_min_required},
#'   \code{coverage}, ordered by \code{pathway_id}.
#' @examples
#' cpd <- data.frame(pathway = c("map00010", "map00010", "map00020"),
#'                   compound = c("C00001", "C00002", "C00003"))
#' pathway_feature_coverage(cpd, "C00001", min_measured = 2)
pathway_feature_coverage <- function(memberships, measured, min_measured,
                                     max_measured = Inf, kegg_org = NULL,
                                     member_col = "compound") {
    mem <- data.frame(pathway = as.character(memberships$pathway),
                      member = as.character(memberships[[member_col]]),
                      stringsAsFactors = FALSE)
    if ("name" %in% names(memberships)) mem$name <- as.character(memberships$name)
    mem <- mem[!is.na(mem$pathway) & !is.na(mem$member), , drop = FALSE]
    mem <- mem[!duplicated(mem[, c("pathway", "member")]), , drop = FALSE]
    measured <- unique(as.character(measured[!is.na(measured)]))

    ids <- sort(unique(mem$pathway))
    n_in <- as.integer(table(factor(mem$pathway, levels = ids)))
    hit <- mem[mem$member %in% measured, , drop = FALSE]
    n_meas <- as.integer(table(factor(hit$pathway, levels = ids)))
    pw_name <- if ("name" %in% names(mem)) {
        unname(stats::setNames(mem$name, mem$pathway)[ids])
    } else {
        rep(NA_character_, length(ids))
    }

    data.frame(
        norm_id        = normalize_pathway_join_key(ids, kegg_org),
        pathway_id     = ids,
        pathway_name   = pw_name,
        n_in_pathway   = n_in,
        n_measured     = n_meas,
        n_min_required = rep(as.integer(min_measured), length(ids)),
        coverage       = .coverage_class(n_meas, as.integer(min_measured), max_measured),
        stringsAsFactors = FALSE
    )
}


#' KEGG compound coverage of the measured metabolites
#'
#' The metabolite case of \code{pathway_feature_coverage()}, on the cached
#' compound-pathway table.
#'
#' @param cpd_pathways Compound-pathway associations (\code{pathway},
#'   \code{compound}, optionally \code{name}).
#' @param measured_cpds KEGG compound ids measured and mapped in this run.
#' @param min_measured Minimum, from \code{compound_min_measured()}.
#' @param max_measured Maximum measured compounds; \code{Inf} for none.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return See \code{pathway_feature_coverage()}.
compound_pathway_coverage <- function(cpd_pathways, measured_cpds, min_measured,
                                      max_measured = Inf, kegg_org = NULL) {
    pathway_feature_coverage(cpd_pathways, measured_cpds, min_measured,
                             max_measured, kegg_org, member_col = "compound")
}


#' KEGG gene-pathway links of an organism, cached
#'
#' Read from \code{<cache_dir>/kegg_gene_pathways_<org>.rds}; on a cold cache,
#' downloaded once from the KEGG REST API (\code{link/pathway/<org>}, the same
#' resource \code{run_ora_kegg_fisher()} reads) and written there.
#'
#' @param kegg_org KEGG organism code.
#' @param cache_dir Cache directory, or NULL for no caching.
#' @param download Whether a cold cache may be filled from the network.
#' @return Data frame with \code{pathway} and \code{gene} (organism prefix
#'   removed), or NULL.
get_kegg_gene_pathways <- function(kegg_org, cache_dir = NULL, download = TRUE) {
    if (is.null(kegg_org) || !nzchar(kegg_org)) return(NULL)
    cache_file <- if (is.null(cache_dir)) NULL
                  else file.path(cache_dir, paste0("kegg_gene_pathways_", kegg_org, ".rds"))
    if (!is.null(cache_file) && file.exists(cache_file)) {
        cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
        if (is.data.frame(cached) && all(c("pathway", "gene") %in% names(cached))) {
            return(cached)
        }
    }
    if (!download) return(NULL)

    lines <- tryCatch(
        readLines(paste0("https://rest.kegg.jp/link/pathway/", kegg_org), warn = FALSE),
        error = function(e) {
            message("    KEGG gene-pathway link download failed: ", e$message)
            NULL
        })
    if (is.null(lines) || length(lines) == 0) return(NULL)
    parts <- strsplit(lines, "\t")
    out <- data.frame(
        pathway = sub("^path:", "", vapply(parts, `[`, character(1), 2)),
        gene    = sub("^[^:]+:", "", vapply(parts, `[`, character(1), 1)),
        stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$pathway) & !is.na(out$gene), , drop = FALSE]
    if (!is.null(cache_file)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        tryCatch(saveRDS(out, cache_file), error = function(e) NULL)
    }
    out
}


#' Pathway coverage settings from the config, with defaults
#'
#' On unless switched off. The config validator fills the same default.
#'
#' @param config Full config object.
#' @return List with \code{enabled} (logical).
#' @keywords internal
.pathway_coverage_config <- function(config) {
    cfg <- config$modes$multiomics$enrichment$pathway_coverage %||% list()
    list(enabled = !isFALSE(cfg$enabled %||% TRUE))
}


#' Coverage of one gene-based layer
#'
#' @keywords internal
#' @return List with \code{pathways} and \code{summary}, or NULL.
.gene_layer_coverage <- function(om, de_data, harmonization_res, config,
                                 gene_pathways, kegg_org, min_measured = NULL) {
    de_tables <- extract_de_tables(de_data, om, harmonization_res)
    if (length(de_tables) == 0) return(NULL)
    feature_ids <- unique(unlist(lapply(de_tables, function(df) df$feature_id)))
    org_db <- get_organism_db(config$global$organism)
    id_map <- if (is.null(org_db)) NULL else tryCatch(
        map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
        error = function(e) NULL)
    entrez <- if (is.null(id_map)) character(0)
              else unique(as.character(id_map$ENTREZID[!is.na(id_map$ENTREZID)]))

    # Most KEGG organism codes key genes on the NCBI Gene ID; for the others
    # the ids are converted, as the enrichment step converts them.
    numeric_keys <- mean(grepl("^[0-9]+$", utils::head(gene_pathways$gene, 2000))) > 0.9
    measured <- if (length(entrez) == 0) character(0)
                else if (numeric_keys) entrez
                else {
                    conv <- tryCatch(convert_entrez_to_kegg(entrez, kegg_org),
                                     error = function(e) NULL)
                    if (is.null(conv)) character(0)
                    else unique(sub("^[^:]+:", "", as.character(conv)))
                }

    enrich_cfg <- config$modes$multiomics$enrichment
    min_req <- min_measured %||% gene_min_measured(enrich_cfg$min_set_size %||% 10,
                                                   length(measured))
    cov <- pathway_feature_coverage(gene_pathways, measured, min_req,
                                    enrich_cfg$max_set_size %||% 500, kegg_org,
                                    member_col = "gene")
    .finish_layer_coverage(om, cov, feature_ids,
                           n_mapped = if (is.null(id_map)) 0L else nrow(id_map),
                           measured = measured, min_req = min_req)
}


#' Coverage of the metabolomics layer
#'
#' @keywords internal
#' @return List with \code{pathways} and \code{summary}, or NULL.
.metabolite_layer_coverage <- function(de_data, harmonization_res, config,
                                       cpd_pathways, kegg_org, min_measured = NULL) {
    de_tables <- extract_de_tables(de_data, .METABOLITE_LAYER, harmonization_res)
    feature_ids <- unique(unlist(lapply(de_tables, function(df) df$feature_id)))
    id_map <- if (length(de_tables) == 0) NULL else tryCatch(
        map_metabolite_ids_to_kegg(de_tables, harmonization_res),
        error = function(e) NULL)
    measured <- if (is.null(id_map)) character(0)
                else unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])

    enrich_cfg <- config$modes$multiomics$enrichment
    min_req <- min_measured %||% compound_min_measured(enrich_cfg$min_set_size %||% 10)
    cov <- compound_pathway_coverage(cpd_pathways, measured, min_req,
                                     enrich_cfg$max_set_size %||% 500, kegg_org)
    .finish_layer_coverage(.METABOLITE_LAYER, cov, feature_ids,
                           n_mapped = if (is.null(id_map)) 0L else nrow(id_map),
                           measured = measured, min_req = min_req)
}


#' Attach the summary row and the no-mapping class to a layer's coverage
#' @keywords internal
.finish_layer_coverage <- function(om, cov, feature_ids, n_mapped, measured, min_req) {
    status <- if (length(measured) == 0) "no_id_mapping" else "mapped"
    if (identical(status, "no_id_mapping")) cov$coverage <- rep("no_id_mapping", nrow(cov))
    summary <- data.frame(
        layer                    = om,
        status                   = status,
        n_features               = length(feature_ids),
        n_mapped_features        = as.integer(n_mapped),
        n_distinct_ids           = length(measured),
        n_pathways_with_members  = nrow(cov),
        n_pathways_testable      = sum(cov$coverage == "testable"),
        n_min_required           = as.integer(min_req),
        stringsAsFactors         = FALSE
    )
    list(pathways = cov, summary = summary)
}


#' KEGG pathway coverage of every layer in this run
#'
#' For each layer with DE results: which KEGG pathways it could test, from its
#' own measured features. Metabolomics uses the compound-pathway cache the
#' compound ORA step filled (never the network); gene and protein layers use
#' the organism's KEGG gene-pathway links (\code{get_kegg_gene_pathways()}),
#' and are skipped when the organism has no KEGG code.
#'
#' The minimum per layer is the producer's own: \code{compound_min_measured()}
#' for metabolomics and \code{gene_min_measured()} for the gene layers, unless
#' \code{min_measured} names one.
#'
#' @param de_results Named list of per-omics DE results.
#' @param harmonization_res Harmonization result, for the ID mappings.
#' @param config Full config object.
#' @param compound_cache_dir Directory holding \code{kegg_compound_pathways.rds}.
#' @param gene_cache_dir Directory for the gene-pathway cache.
#' @param min_measured Named list of per-layer minima overriding the default.
#' @return List with \code{layers} (named list of per-layer tables from
#'   \code{pathway_feature_coverage()}) and \code{summary} (one row per layer);
#'   NULL, with a message, when switched off or no layer could be covered.
pathway_coverage_for_run <- function(de_results, harmonization_res, config,
                                     compound_cache_dir, gene_cache_dir,
                                     min_measured = list()) {
    if (!.pathway_coverage_config(config)$enabled) {
        message("  Pathway coverage switched off in config")
        return(NULL)
    }
    kegg_org <- resolve_kegg_org_code(config$global$organism)
    layers <- list()

    if (!is.null(de_results[[.METABOLITE_LAYER]])) {
        cpd <- .cached_compound_pathways(compound_cache_dir)
        if (is.null(cpd)) {
            message("  Pathway coverage: no cached KEGG compound-pathway table at ",
                    compound_cache_dir %||% "<none>", "; metabolomics skipped")
        } else {
            layers[[.METABOLITE_LAYER]] <- tryCatch(
                .metabolite_layer_coverage(de_results[[.METABOLITE_LAYER]],
                                           harmonization_res, config, cpd, kegg_org,
                                           min_measured[[.METABOLITE_LAYER]]),
                error = function(e) {
                    message("  Pathway coverage (metabolomics) failed: ", e$message)
                    NULL
                })
        }
    }

    gene_layers <- intersect(c("transcriptomics", "proteomics"), names(de_results))
    if (length(gene_layers) > 0) {
        gene_pw <- if (is.null(kegg_org)) NULL
                   else get_kegg_gene_pathways(kegg_org, gene_cache_dir)
        if (is.null(gene_pw)) {
            message("  Pathway coverage: no KEGG gene-pathway links for this ",
                    "organism; gene layers skipped")
        } else {
            for (om in gene_layers) {
                layers[[om]] <- tryCatch(
                    .gene_layer_coverage(om, de_results[[om]], harmonization_res,
                                         config, gene_pw, kegg_org, min_measured[[om]]),
                    error = function(e) {
                        message("  Pathway coverage (", om, ") failed: ", e$message)
                        NULL
                    })
            }
        }
    }

    layers <- Filter(Negate(is.null), layers)
    if (length(layers) == 0) return(NULL)
    for (om in names(layers)) {
        s <- layers[[om]]$summary
        message("  Pathway coverage, ", om, ": ", s$n_pathways_testable, " of ",
                s$n_pathways_with_members, " KEGG pathways have at least ",
                s$n_min_required, " measured members")
    }
    list(layers = lapply(layers, `[[`, "pathways"),
         summary = do.call(rbind, lapply(layers, `[[`, "summary")))
}


#' The same coverage with other per-layer minima
#'
#' Counts are unchanged; only the class follows the new minimum. Used where a
#' producer applies a different floor (Multi-ORA).
#'
#' @param coverage Result of \code{pathway_coverage_for_run()}, or NULL.
#' @param min_measured Named list of per-layer minima.
#' @param max_measured Maximum measured members; \code{Inf} for none.
#' @return The coverage with classes and \code{n_min_required} recomputed for
#'   the named layers.
reclassify_pathway_coverage <- function(coverage, min_measured, max_measured = Inf) {
    if (is.null(coverage)) return(NULL)
    for (om in intersect(names(min_measured), names(coverage$layers))) {
        tab <- coverage$layers[[om]]
        k <- as.integer(min_measured[[om]])
        unmapped <- tab$coverage == "no_id_mapping"
        tab$n_min_required <- rep(k, nrow(tab))
        tab$coverage[!unmapped] <- .coverage_class(tab$n_measured[!unmapped], k,
                                                   max_measured)
        coverage$layers[[om]] <- tab
        if (is.data.frame(coverage$summary)) {
            i <- coverage$summary$layer == om
            coverage$summary$n_min_required[i] <- k
            coverage$summary$n_pathways_testable[i] <- sum(tab$coverage == "testable")
        }
    }
    coverage
}


#' Coverage columns of one layer for a vector of pathway keys
#'
#' A KEGG key absent from the layer's table is a pathway with no members of
#' that omics type (\code{no_kegg_members}, counts 0); a key that is not a KEGG
#' accession -- a GO term, a custom set -- gets NA, so non-KEGG collections are
#' left as they were.
#'
#' @param keys Pathway keys.
#' @param tab One layer's coverage table.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return Data frame with \code{n_in_pathway}, \code{n_measured},
#'   \code{n_min_required}, \code{coverage}, one row per key.
layer_coverage_for_keys <- function(keys, tab, kegg_org = NULL) {
    keys <- normalize_pathway_join_key(as.character(keys), kegg_org)
    idx <- match(keys, tab$norm_id)
    absent <- is.na(idx) & is_kegg_pathway_accession(keys, kegg_org)
    unmapped <- nrow(tab) > 0 && all(tab$coverage == "no_id_mapping")
    out <- data.frame(n_in_pathway = tab$n_in_pathway[idx],
                      n_measured = tab$n_measured[idx],
                      n_min_required = tab$n_min_required[idx],
                      coverage = tab$coverage[idx],
                      stringsAsFactors = FALSE)
    out$n_in_pathway[absent] <- 0L
    out$n_measured[absent] <- 0L
    out$n_min_required[absent] <- if (nrow(tab) > 0) tab$n_min_required[1] else NA_integer_
    out$coverage[absent] <- if (unmapped) "no_id_mapping" else "no_kegg_members"
    out
}


#' Attach per-layer KEGG coverage to a table keyed on pathway identity
#'
#' Adds \code{n_in_pathway_<layer>}, \code{n_measured_<layer>},
#' \code{n_min_required_<layer>} and \code{coverage_<layer>} for each layer the
#' coverage holds (or the subset in \code{layers}).
#'
#' @param df Data frame carrying the key column.
#' @param coverage Result of \code{pathway_coverage_for_run()}; NULL returns
#'   \code{df} unchanged.
#' @param layers Layers to annotate, or NULL for all covered layers.
#' @param key_col Name of the key column.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return \code{df} with the columns added (replacing any already there).
annotate_pathway_coverage <- function(df, coverage, layers = NULL,
                                      key_col = "norm_id", kegg_org = NULL) {
    if (is.null(coverage) || is.null(df) || !key_col %in% names(df)) return(df)
    layers <- intersect(layers %||% names(coverage$layers), names(coverage$layers))
    for (om in layers) {
        cols <- layer_coverage_for_keys(df[[key_col]], coverage$layers[[om]], kegg_org)
        for (cn in names(cols)) df[[paste0(cn, "_", om)]] <- cols[[cn]]
    }
    df
}


#' Display state of a layer's cell in a KEGG figure
#'
#' The one rule every figure uses. A cell with a value is \code{scored} and
#' drawn on the figure's own scale. A cell without one is \code{neutral} when
#' the layer could test the pathway (white), and \code{not_testable} when it
#' could not -- no KEGG members of that omics type, none measured, fewer than
#' the minimum, or no ids mapped -- drawn grey with a cross. Unknown coverage
#' (non-KEGG pathways, no coverage) is NA, drawn as the figures always did.
#'
#' @param values Values shown (NA where there is no result); vector or matrix.
#' @param coverage_class Coverage class per cell, same length.
#' @return Character vector (or matrix shaped like \code{values}) of states.
#' @examples
#' pathway_cell_state(c(2, NA, NA, NA), c("testable", "testable",
#'                                         "no_kegg_members", NA))
#' # "scored" "neutral" "not_testable" NA
pathway_cell_state <- function(values, coverage_class) {
    v <- as.vector(values)
    cls <- as.character(coverage_class)
    if (length(cls) != length(v)) cls <- rep_len(cls, length(v))
    state <- rep(NA_character_, length(v))
    empty <- is.na(v)
    state[!empty] <- "scored"
    state[empty & cls %in% "testable"] <- "neutral"
    state[empty & cls %in% .NOT_TESTABLE_CLASSES] <- "not_testable"
    if (!is.null(dim(values))) {
        dim(state) <- dim(values)
        dimnames(state) <- dimnames(values)
    }
    state
}


#' Cell states for a heatmap whose columns are layers
#'
#' @param m Value matrix with layer names as columns.
#' @param classes Named list: per layer, the coverage class of each row of
#'   \code{m}. Layers without an entry are left NA.
#' @return Character matrix shaped like \code{m}, or NULL when no column has
#'   coverage.
#' @keywords internal
.coverage_state_matrix <- function(m, classes) {
    cols <- intersect(colnames(m), names(Filter(Negate(is.null), classes)))
    cols <- cols[vapply(cols, function(cn) any(!is.na(classes[[cn]])), logical(1))]
    if (length(cols) == 0) return(NULL)
    states <- matrix(NA_character_, nrow(m), ncol(m), dimnames = dimnames(m))
    for (cn in cols) states[, cn] <- pathway_cell_state(m[, cn], classes[[cn]])
    states
}


#' Colours, breaks and marks for a heatmap with coverage states
#'
#' The figure's own scale is kept exactly: the same palette over the same
#' breaks. Neutral cells get a value below the scale and one extra white
#' interval underneath it, with a spacer interval in the palette's first colour
#' so the scale's lowest real value keeps that colour (the bins are closed on
#' the right). Not-testable cells stay NA, drawn in the figure's NA grey, and
#' carry a cross.
#'
#' @param m Value matrix as drawn.
#' @param states State matrix from \code{.coverage_state_matrix()}, or NULL.
#' @param palette Colour vector of the figure's scale.
#' @param zlim Value range of the scale, from \code{.nondegenerate_range()}.
#' @return List: \code{values}, \code{colors}, \code{breaks},
#'   \code{legend_breaks} (NA for the default legend), \code{marks} (character
#'   matrix, or FALSE) and \code{has_states}.
#' @keywords internal
.coverage_heatmap_scale <- function(m, states, palette, zlim) {
    breaks <- seq(zlim[1], zlim[2], length.out = length(palette) + 1)
    out <- list(values = m, colors = palette, breaks = breaks,
                legend_breaks = NA, marks = FALSE, has_states = FALSE)
    if (is.null(states)) return(out)

    neutral <- !is.na(states) & states == "neutral"
    untestable <- !is.na(states) & states == "not_testable"
    if (any(neutral)) {
        d <- diff(zlim) / length(palette)
        out$values[neutral] <- zlim[1] - 1.5 * d
        out$colors <- c("white", palette[1], palette)
        out$breaks <- c(zlim[1] - 2 * d, zlim[1] - d, breaks)
        out$legend_breaks <- grid::grid.pretty(zlim)
    }
    if (any(untestable)) {
        marks <- matrix("", nrow(m), ncol(m), dimnames = dimnames(m))
        marks[untestable] <- "×"
        out$marks <- marks
    }
    out$has_states <- any(neutral) || any(untestable)
    out
}


#' Lookup status of a target cell, from its coverage
#'
#' Used only where the target has no result for the pathway. A testable
#' pathway becomes "tested, not reported" only when the target has results in
#' that contrast at all; otherwise nothing says it was tested there and the
#' original status stands.
#'
#' @param status Current \code{target_status}.
#' @param coverage_class Coverage class per row.
#' @param target_has_contrast Logical per row.
#' @return Character vector of statuses.
#' @keywords internal
.coverage_lookup_status <- function(status, coverage_class, target_has_contrast) {
    open <- status == "not in target results"
    cls <- as.character(coverage_class)
    labels <- c(no_kegg_members     = "not testable: no KEGG members",
                no_measured_members = "not testable: none measured",
                below_min_measured  = "not testable: fewer than the minimum measured",
                above_max_measured  = "not testable: above the maximum set size",
                no_id_mapping       = "not testable: no ids mapped")
    hit <- open & cls %in% names(labels)
    status[hit] <- unname(labels[cls[hit]])
    status[open & cls %in% "testable" & target_has_contrast] <- "tested, not reported"
    status
}


#' Remove the coverage exports a previous run left here
#'
#' @param out_dir Output directory.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_pathway_coverage <- function(out_dir) {
    if (is.null(out_dir) || !dir.exists(out_dir)) return(invisible(character(0)))
    stale <- file.path(out_dir, c("kegg_pathway_coverage.csv",
                                  "kegg_pathway_coverage_summary.tsv"))
    stale <- stale[file.exists(stale)]
    if (length(stale) > 0) unlink(stale)
    invisible(stale)
}


#' Write the per-layer coverage table and its summary
#'
#' Previous copies are removed first, so the files exist only if this call had
#' coverage to write. \code{kegg_pathway_coverage.csv} is long, one row per
#' layer and pathway: every KEGG pathway with at least one member of that
#' layer's type, plus every KEGG pathway in \code{extra_keys} that has none
#' (\code{no_kegg_members}).
#'
#' @param coverage Result of \code{pathway_coverage_for_run()}, or NULL.
#' @param out_dir Output directory.
#' @param extra_keys Pathway keys appearing in this run's results, or NULL.
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return Invisibly, the paths written, or NULL.
write_pathway_coverage <- function(coverage, out_dir, extra_keys = NULL,
                                   kegg_org = NULL) {
    .clear_pathway_coverage(out_dir)
    if (is.null(coverage) || is.null(out_dir)) return(invisible(NULL))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    keys <- unique(normalize_pathway_join_key(as.character(extra_keys), kegg_org))
    keys <- sort(keys[is_kegg_pathway_accession(keys, kegg_org)])
    rows <- lapply(names(coverage$layers), function(om) {
        tab <- coverage$layers[[om]]
        miss <- setdiff(keys, tab$norm_id)
        if (length(miss) > 0) {
            extra <- cbind(data.frame(norm_id = miss, pathway_id = NA_character_,
                                      pathway_name = NA_character_,
                                      stringsAsFactors = FALSE),
                           layer_coverage_for_keys(miss, tab, kegg_org))
            tab <- rbind(tab, extra[, names(tab)])
        }
        cbind(layer = rep(om, nrow(tab)), tab, stringsAsFactors = FALSE)
    })
    paths <- c(csv = file.path(out_dir, "kegg_pathway_coverage.csv"))
    utils::write.csv(do.call(rbind, rows), paths[["csv"]], row.names = FALSE)
    if (is.data.frame(coverage$summary)) {
        paths[["summary"]] <- file.path(out_dir, "kegg_pathway_coverage_summary.tsv")
        utils::write.table(coverage$summary, paths[["summary"]], sep = "\t",
                           quote = FALSE, row.names = FALSE)
    }
    invisible(paths)
}


#' KEGG pathway keys appearing in a list of enrichment tables
#'
#' @param tables List of enrichment data frames (NULL entries allowed).
#' @param kegg_org Active KEGG organism code, or NULL.
#' @return Character vector of unique join keys.
#' @keywords internal
.result_pathway_keys <- function(tables, kegg_org = NULL) {
    keys <- unlist(lapply(tables, function(df) {
        if (is.data.frame(df) && nrow(df) > 0) pathway_join_key(df, kegg_org) else NULL
    }), use.names = FALSE)
    unique(keys[!is.na(keys)])
}


#' Whether a result table carries pathway coverage, for the report
#'
#' Reads what the producer wrote; decides nothing.
#'
#' @param result_csv Path of the table a figure was drawn from.
#' @return TRUE when the table has at least one \code{coverage_<layer>} column
#'   with a classified row.
has_pathway_coverage <- function(result_csv) {
    if (is.null(result_csv) || !file.exists(result_csv)) return(FALSE)
    tbl <- tryCatch(utils::read.csv(result_csv, stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(tbl)) return(FALSE)
    cols <- grep("^coverage_", names(tbl), value = TRUE)
    any(vapply(cols, function(cn) any(!is.na(tbl[[cn]]) & nzchar(tbl[[cn]])),
               logical(1)))
}


#' Legend text explaining the coverage marks of one figure
#'
#' Plain and complete: every state the figure can show. Built here, beside the
#' rule it describes, rather than in the report.
#'
#' @param figure One of "heatmap" (meta-analysis and ORA heatmaps), "lookup",
#'   "multi_ora" or "multigsea".
#' @return A single string.
pathway_coverage_legend <- function(figure = c("heatmap", "lookup", "multi_ora",
                                               "multigsea")) {
    figure <- match.arg(figure)
    why <- paste("the layer is not expected to give a result there: KEGG lists no",
                 "features of that omics type in the pathway, none were measured,",
                 "or fewer than the minimum needed for a test were measured")
    tail <- paste("The class and counts per layer are in the table",
                  "(coverage_, n_measured_ and n_in_pathway_ columns) and in",
                  "kegg_pathway_coverage.csv.")
    switch(figure,
        heatmap = paste0(
            "KEGG figures, per layer: a coloured cell is a result on the scale ",
            "above. White: the layer could test the pathway but reported no ",
            "result. Grey with a cross: ", why, ". Plain grey: no result and no ",
            "coverage information. ", tail),
        lookup = paste0(
            "Looked-up cells without a result: white, the layer could test the ",
            "pathway but has no result in that contrast; grey with a cross, ", why,
            ". target_status and target_coverage give the case."),
        multi_ora = paste0(
            "Markers at zero stand for a layer without a result: an open circle, ",
            "the layer could test the pathway but did not report it; a cross, ",
            why, ". ", tail),
        multigsea = paste0(
            "A grey cross over a point: for at least one of the two layers, ",
            why, "; that layer's value is plotted at zero. Per-term coverage is ",
            "in the CSV beside the figure."))
}


#' Pathway coverage for the Multi-ORA, on its own minima
#'
#' Reuses the cross-omics coverage when it is supplied, reclassified on the
#' Multi-ORA floors (\code{.multi_ora_min_measured()}). A layer it lacks --
#' metabolomics when the cross-omics step had no compound cache -- is computed
#' here from the cache the Multi-ORA compound ORA itself filled.
#'
#' @param pathway_coverage Cross-omics coverage, or NULL.
#' @param de_results Named list of per-omics DE results.
#' @param harmonization_res Harmonization result.
#' @param config Full config object.
#' @param out_dir Multi-ORA output directory (its compound cache lives there).
#' @return Coverage list, or NULL.
multi_ora_pathway_coverage <- function(pathway_coverage, de_results,
                                       harmonization_res, config, out_dir) {
    if (!.pathway_coverage_config(config)$enabled) return(NULL)
    mins <- .multi_ora_min_measured()
    cov <- pathway_coverage
    if (!is.null(de_results[[.METABOLITE_LAYER]]) &&
        is.null(cov$layers[[.METABOLITE_LAYER]])) {
        own <- tryCatch(
            pathway_coverage_for_run(de_results[.METABOLITE_LAYER], harmonization_res,
                                     config, compound_cache_dir = out_dir,
                                     gene_cache_dir = NULL),
            error = function(e) NULL)
        if (!is.null(own)) {
            cov$layers[[.METABOLITE_LAYER]] <- own$layers[[.METABOLITE_LAYER]]
            cov$summary <- rbind(cov$summary, own$summary)
        }
    }
    if (is.null(cov) || length(cov$layers) == 0) return(NULL)
    reclassify_pathway_coverage(cov, mins,
        max_measured = config$modes$multiomics$enrichment$max_set_size %||% 500)
}


#' Whether any MultiGSEA pair table carries a not-testable mark
#'
#' Reads the CSVs the pair plots wrote beside themselves; decides nothing.
#'
#' @param csv_paths Paths of multigsea_<a>_vs_<b>.csv files.
#' @return TRUE when any row has \code{state_x} or \code{state_y} equal to
#'   "not_testable".
multigsea_has_coverage_marks <- function(csv_paths) {
    for (f in csv_paths[file.exists(csv_paths)]) {
        tbl <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE),
                        error = function(e) NULL)
        if (is.null(tbl)) next
        for (cn in intersect(c("state_x", "state_y"), names(tbl))) {
            if (any(tbl[[cn]] %in% "not_testable")) return(TRUE)
        }
    }
    FALSE
}
