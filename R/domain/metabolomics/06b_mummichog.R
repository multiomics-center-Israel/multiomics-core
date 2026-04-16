# R/domain/metabolomics/06b_mummichog.R
#
# Mummichog pathway enrichment for metabolomics via reticulate.
# Builds organism-specific metabolic models from KEGG (Azimuth JSON format),
# then calls the mummichog Python CLI for pathway/module analysis.
#
# Supports: honey bee (ame), E. coli (eco), combined (ame+eco).
#
# Dependencies: KEGGREST, reticulate, jsonlite, ggplot2


# ==== PYTHON ENVIRONMENT =====================================================

#' Ensure mummichog Python virtualenv exists
#'
#' Creates a project-local virtualenv and installs mummichog on first run.
#'
#' @param project_dir Project root directory.
#' @return Path to the virtualenv's Python executable.
ensure_mummichog_env <- function(project_dir = ".") {
    venv_path <- file.path(normalizePath(project_dir), ".mummichog_venv")
    python_bin <- file.path(venv_path, "bin", "python")

    if (!file.exists(python_bin)) {
        message("Creating mummichog virtualenv at: ", venv_path)
        if (!requireNamespace("reticulate", quietly = TRUE))
            stop("reticulate package required for mummichog")
        reticulate::virtualenv_create(venv_path)
        reticulate::virtualenv_install(venv_path, "mummichog")
    }

    python_bin
}


# ==== KEGG MODEL BUILDER =====================================================

#' Build a mummichog-native JSON metabolic model from KEGG
#'
#' Fetches pathway and compound data from KEGG for a given organism
#' and writes a JSON file in mummichog's native model format.
#' This format includes: Compounds (with mw), metabolic_pathways,
#' cpd_edges, cpd2pathways, dict_cpds_def, metabolic_rxns.
#'
#' @param organism  KEGG organism code (e.g., "ame", "eco").
#' @param output_dir Directory to save the JSON model file.
#' @param verbose    Print progress messages.
#' @return Path to the generated JSON model file.
build_kegg_organism_model <- function(organism, output_dir, verbose = TRUE) {
    if (!requireNamespace("KEGGREST", quietly = TRUE))
        stop("KEGGREST package required")
    if (!requireNamespace("jsonlite", quietly = TRUE))
        stop("jsonlite package required")

    out_file <- file.path(output_dir, paste0(organism, "_metabolic_model.json"))
    if (file.exists(out_file)) {
        if (verbose) message("Using cached model: ", out_file)
        return(out_file)
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) message("Building KEGG metabolic model for organism: ", organism)

    # Step 1: Get all metabolic pathways
    if (verbose) message("  Fetching pathway list...")
    pathway_list <- tryCatch(
        KEGGREST::keggList("pathway", organism),
        error = function(e) {
            warning("Failed to fetch KEGG pathways for ", organism, ": ", e$message)
            return(character(0))
        }
    )
    if (length(pathway_list) == 0) {
        warning("No pathways found for organism: ", organism)
        return(NULL)
    }

    pathway_ids <- names(pathway_list)
    pathway_names <- unname(pathway_list)
    # Strip species suffix: "Fatty acid metabolism - Apis mellifera (honey bee)" -> "Fatty acid metabolism"
    pathway_names <- sub(" - [A-Z].*$", "", pathway_names)
    if (verbose) message("  Found ", length(pathway_ids), " pathways")

    # Step 2: Fetch pathway details — get compounds per pathway
    if (verbose) message("  Fetching pathway details...")
    batch_size <- 10
    pathway_compounds <- list()  # pathway_id -> vector of compound IDs
    pathway_rxn_ids <- list()    # pathway_id -> vector of reaction IDs
    all_compound_ids <- character(0)

    for (i in seq(1, length(pathway_ids), by = batch_size)) {
        batch <- pathway_ids[i:min(i + batch_size - 1, length(pathway_ids))]
        tryCatch({
            pw_data <- KEGGREST::keggGet(batch)
            for (pw in pw_data) {
                pw_id <- pw$ENTRY
                cpd_ids <- names(pw$COMPOUND)
                if (!is.null(cpd_ids)) {
                    pathway_compounds[[pw_id]] <- cpd_ids
                    all_compound_ids <- c(all_compound_ids, cpd_ids)
                }
                rxn_ids <- names(pw$REACTION)
                if (!is.null(rxn_ids)) pathway_rxn_ids[[pw_id]] <- rxn_ids
            }
            Sys.sleep(0.5)
        }, error = function(e) {
            if (verbose) message("  Warning: batch fetch failed: ", e$message)
        })
        if (verbose && i %% 50 == 1)
            message("    Processed ", min(i + batch_size - 1, length(pathway_ids)),
                    "/", length(pathway_ids), " pathways")
    }

    all_compound_ids <- unique(all_compound_ids)
    if (verbose) message("  Found ", length(all_compound_ids), " unique compounds across pathways")

    # Step 3: Fetch reaction details for compound edges
    all_rxn_ids <- unique(unlist(pathway_rxn_ids))
    all_reactions <- list()
    cpd_edges <- list()

    if (length(all_rxn_ids) > 0) {
        if (verbose) message("  Fetching ", length(all_rxn_ids), " reaction details...")
        for (i in seq(1, length(all_rxn_ids), by = batch_size)) {
            batch <- all_rxn_ids[i:min(i + batch_size - 1, length(all_rxn_ids))]
            tryCatch({
                rxn_data <- KEGGREST::keggGet(paste0("rn:", batch))
                for (rxn in rxn_data) {
                    rxn_id <- rxn$ENTRY
                    eq <- rxn$EQUATION
                    if (!is.null(eq)) {
                        parts <- strsplit(eq, "<=>")[[1]]
                        if (length(parts) == 2) {
                            subs <- regmatches(parts[1], gregexpr("C[0-9]{5}", parts[1]))[[1]]
                            prods <- regmatches(parts[2], gregexpr("C[0-9]{5}", parts[2]))[[1]]
                            all_reactions[[rxn_id]] <- list(
                                id = rxn_id, reactants = subs, products = prods,
                                cpds = unique(c(subs, prods)), source = paste0("kegg_", organism)
                            )
                            for (s in subs) for (p in prods) cpd_edges <- c(cpd_edges, list(c(s, p)))
                            all_compound_ids <- c(all_compound_ids, subs, prods)
                        }
                    }
                }
                Sys.sleep(0.5)
            }, error = function(e) {
                if (verbose) message("  Warning: reaction fetch failed: ", e$message)
            })
            if (verbose && i %% 50 == 1)
                message("    Processed ", min(i + batch_size - 1, length(all_rxn_ids)),
                        "/", length(all_rxn_ids), " reactions")
        }
    }

    # If no reactions found, build edges from pathway co-membership
    if (length(cpd_edges) == 0) {
        if (verbose) message("  No reactions — building edges from pathway co-membership")
        for (pw_cpds in pathway_compounds) {
            if (length(pw_cpds) >= 2) {
                pairs <- utils::combn(pw_cpds, 2, simplify = FALSE)
                cpd_edges <- c(cpd_edges, pairs)
            }
        }
    }

    # Deduplicate edges
    if (length(cpd_edges) > 0) {
        edge_strs <- vapply(cpd_edges, function(e) paste(sort(e), collapse = ","), character(1))
        cpd_edges <- cpd_edges[!duplicated(edge_strs)]
    }
    if (verbose) message("  Built ", length(cpd_edges), " compound edges")

    all_compound_ids <- unique(all_compound_ids)

    # Step 4: Fetch compound details (mass, formula, name)
    if (verbose) message("  Fetching ", length(all_compound_ids), " compound details...")
    compounds <- list()      # mummichog format: id -> {formula, mw, name}
    dict_cpds_def <- list()  # id -> name

    for (i in seq(1, length(all_compound_ids), by = batch_size)) {
        batch <- all_compound_ids[i:min(i + batch_size - 1, length(all_compound_ids))]
        tryCatch({
            cpd_data <- KEGGREST::keggGet(batch)
            for (cpd in cpd_data) {
                cpd_id <- cpd$ENTRY
                mass <- if (!is.null(cpd$EXACT_MASS)) as.numeric(cpd$EXACT_MASS) else 0
                formula <- if (!is.null(cpd$FORMULA)) cpd$FORMULA else ""
                name <- if (!is.null(cpd$NAME)) cpd$NAME[1] else cpd_id
                name <- sub(";$", "", name)
                compounds[[cpd_id]] <- list(formula = formula, mw = mass, name = name, adducts = list())
                dict_cpds_def[[cpd_id]] <- name
            }
            Sys.sleep(0.5)
        }, error = function(e) {
            if (verbose) message("  Warning: compound fetch failed: ", e$message)
        })
        if (verbose && i %% 100 == 1)
            message("    Processed ", min(i + batch_size - 1, length(all_compound_ids)),
                    "/", length(all_compound_ids), " compounds")
    }

    # Step 5: If no reactions, create synthetic reactions from pathway compounds
    # Each pathway's compounds become a chain of reactions for the Azimuth format
    if (length(all_reactions) == 0 && length(pathway_compounds) > 0) {
        if (verbose) message("  Creating synthetic reactions from pathway co-membership...")
        rxn_counter <- 0
        for (pw_id in names(pathway_compounds)) {
            cpds <- pathway_compounds[[pw_id]]
            if (length(cpds) >= 2) {
                # Create pairwise reactions between consecutive compounds
                for (k in seq(1, length(cpds) - 1)) {
                    rxn_counter <- rxn_counter + 1
                    rxn_id <- sprintf("SYN_%s_%04d", pw_id, rxn_counter)
                    all_reactions[[rxn_id]] <- list(
                        id = rxn_id,
                        reactants = list(cpds[k]),
                        products = list(cpds[k + 1]),
                        cpds = c(cpds[k], cpds[k + 1]),
                        source = paste0("kegg_", organism)
                    )
                    cpd_edges <- c(cpd_edges, list(c(cpds[k], cpds[k + 1])))
                }
                # Update pathway_rxn_ids
                pw_rxns <- sprintf("SYN_%s_%04d", pw_id, seq_len(length(cpds) - 1) +
                                   rxn_counter - length(cpds) + 1)
                pathway_rxn_ids[[pw_id]] <- pw_rxns
            }
        }
        # Deduplicate edges
        if (length(cpd_edges) > 0) {
            edge_strs <- vapply(cpd_edges, function(e) paste(sort(e), collapse = ","), character(1))
            cpd_edges <- cpd_edges[!duplicated(edge_strs)]
        }
        if (verbose) message("  Created ", length(all_reactions), " synthetic reactions, ",
                             length(cpd_edges), " edges")
    }

    # Step 6: Build Azimuth JSON format (required by mummichog's read_user_json_model)
    list_of_compounds <- lapply(names(compounds), function(cpd_id) {
        cpd <- compounds[[cpd_id]]
        list(
            id = cpd_id,
            name = cpd$name,
            identifiers = list(`kegg.compound` = cpd_id),
            neutral_formula = cpd$formula,
            charge = 0,
            charged_formula = cpd$formula,
            neutral_mono_mass = cpd$mw,
            SMILES = "",
            inchi = ""
        )
    })

    list_of_reactions <- lapply(unname(all_reactions), function(rxn) {
        list(id = rxn$id, reactants = rxn$reactants, products = rxn$products)
    })

    list_of_pathways <- lapply(seq_along(pathway_ids), function(j) {
        pw_id <- pathway_ids[j]
        rxns <- pathway_rxn_ids[[pw_id]]
        if (is.null(rxns)) rxns <- character(0)
        list(id = pw_id, name = pathway_names[j], list_of_reactions = as.list(rxns))
    })

    model <- list(
        id = paste0(organism, "_kegg_model"),
        meta_data = list(
            species = organism,
            version = "1.0",
            source = "KEGG via KEGGREST",
            date = as.character(Sys.Date())
        ),
        list_of_compounds = list_of_compounds,
        list_of_reactions = list_of_reactions,
        list_of_pathways = list_of_pathways
    )

    jsonlite::write_json(model, out_file, auto_unbox = TRUE, pretty = FALSE)
    if (verbose) message("  Model saved: ", out_file)
    if (verbose) message("  Summary: ", length(compounds), " compounds, ",
                         length(all_reactions), " reactions, ",
                         length(list_of_pathways), " pathways")
    out_file
}


#' Build a combined model by merging two organism models (Azimuth format)
#'
#' @param model_files Character vector of model JSON file paths.
#' @param output_dir  Directory for the combined model.
#' @param combined_id Identifier for the combined model.
#' @return Path to combined model JSON.
build_combined_model <- function(model_files, output_dir, combined_id = "ame_eco") {
    out_file <- file.path(output_dir, paste0(combined_id, "_metabolic_model.json"))
    if (file.exists(out_file)) {
        message("Using cached combined model: ", out_file)
        return(out_file)
    }

    models <- lapply(model_files, function(f) jsonlite::read_json(f))

    # Merge compounds (deduplicate by id)
    all_cpds <- do.call(c, lapply(models, `[[`, "list_of_compounds"))
    cpd_ids_seen <- character(0)
    unique_cpds <- list()
    for (cpd in all_cpds) {
        if (!(cpd$id %in% cpd_ids_seen)) {
            cpd_ids_seen <- c(cpd_ids_seen, cpd$id)
            unique_cpds <- c(unique_cpds, list(cpd))
        }
    }

    # Merge reactions (deduplicate by id)
    all_rxns <- do.call(c, lapply(models, `[[`, "list_of_reactions"))
    rxn_ids_seen <- character(0)
    unique_rxns <- list()
    for (rxn in all_rxns) {
        if (!(rxn$id %in% rxn_ids_seen)) {
            rxn_ids_seen <- c(rxn_ids_seen, rxn$id)
            unique_rxns <- c(unique_rxns, list(rxn))
        }
    }

    # Merge pathways (prefix names with organism to avoid confusion)
    all_pws <- list()
    for (m in models) {
        org <- m$meta_data$species
        for (pw in m$list_of_pathways) {
            pw$name <- paste0("[", toupper(org), "] ", pw$name)
            pw$id <- paste0(org, "_", pw$id)
            all_pws <- c(all_pws, list(pw))
        }
    }

    combined <- list(
        id = combined_id,
        meta_data = list(
            species = combined_id,
            version = "1.0",
            source = "Merged KEGG models",
            date = as.character(Sys.Date())
        ),
        list_of_compounds = unique_cpds,
        list_of_reactions = unique_rxns,
        list_of_pathways = all_pws
    )

    jsonlite::write_json(combined, out_file, auto_unbox = TRUE, pretty = FALSE)
    message("Combined model saved: ", out_file, " (",
            length(unique_cpds), " compounds, ", length(unique_rxns),
            " reactions, ", length(all_pws), " pathways)")
    out_file
}


# ==== MUMMICHOG INPUT PREPARATION ============================================

#' Prepare mummichog input file from DE results + annotations
#'
#' @param de_res   data.frame with columns: feature_id, logFC, P.Value.
#' @param row_data data.frame with feature annotations (must have m/z and RT columns).
#' @param config   Pipeline config list.
#' @param output_dir Directory for input file.
#' @return Path to the mummichog input TSV file.
prepare_mummichog_input <- function(de_res, row_data, config, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Find m/z column (handle check.names=TRUE mangling)
    mz_col <- NULL
    for (candidate in c("m/z", "mz", "m.z", "MZ", "Mass", "m.z.")) {
        if (candidate %in% colnames(row_data)) {
            mz_col <- candidate
            break
        }
    }
    if (is.null(mz_col)) {
        # Try grep as fallback
        mz_idx <- grep("^m[./]z", colnames(row_data), ignore.case = TRUE)
        if (length(mz_idx) > 0) mz_col <- colnames(row_data)[mz_idx[1]]
    }
    if (is.null(mz_col)) stop("No m/z column found in row_data. Columns: ",
                               paste(colnames(row_data), collapse = ", "))

    # Find RT column (handle check.names=TRUE mangling)
    rt_col <- NULL
    for (candidate in c("RT [min]", "rt", "RT", "retention_time", "RetentionTime",
                         "RT..min.", "RT.min.")) {
        if (candidate %in% colnames(row_data)) {
            rt_col <- candidate
            break
        }
    }
    if (is.null(rt_col)) {
        rt_idx <- grep("^RT|retention", colnames(row_data), ignore.case = TRUE)
        if (length(rt_idx) > 0) rt_col <- colnames(row_data)[rt_idx[1]]
    }
    if (is.null(rt_col)) stop("No retention time column found in row_data. Columns: ",
                               paste(colnames(row_data), collapse = ", "))

    # Find feature_id column for joining
    id_col <- "feature_id"
    if (!id_col %in% colnames(row_data)) {
        # Try Metabolite column
        if ("Metabolite" %in% colnames(row_data)) {
            row_data$feature_id <- row_data$Metabolite
        } else {
            row_data$feature_id <- rownames(row_data)
        }
    }

    # Merge DE results with annotations
    merged <- merge(de_res, row_data[, c("feature_id", mz_col, rt_col)],
                    by = "feature_id")

    # Create mummichog input format
    input_df <- data.frame(
        `m.z`     = as.numeric(merged[[mz_col]]),
        `r.t`     = as.numeric(merged[[rt_col]]),
        `p.value` = merged$P.Value,
        `t.score` = merged$logFC,  # logFC as proxy for t-statistic
        check.names = FALSE
    )

    # Remove rows with missing values
    input_df <- input_df[complete.cases(input_df), ]

    out_file <- file.path(output_dir, "mummichog_input.tsv")
    write.table(input_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Mummichog input file: ", out_file, " (", nrow(input_df), " features)")
    out_file
}


# ==== MUMMICHOG EXECUTION ====================================================

#' Run mummichog for one organism model
#'
#' @param input_file   Path to mummichog input TSV.
#' @param model_file   Path to organism JSON model.
#' @param output_dir   Directory for mummichog outputs.
#' @param organism     Organism label for output naming.
#' @param config       Pipeline config (for p_cutoff, permutations, etc.).
#' @param python_bin   Path to Python binary in mummichog virtualenv.
#' @return list(table = data.frame, method = "mummichog", organism = organism)
run_mummichog_organism <- function(input_file, model_file, output_dir,
                                   organism, config, python_bin) {
    org_label <- gsub("[+]", "_", organism)
    work_dir  <- file.path(output_dir, paste0("run_", org_label))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()
    p_cutoff     <- mummi_cfg$p_cutoff %||% 0.05
    n_perm       <- mummi_cfg$n_permutations %||% 100
    ppm          <- mummi_cfg$tolerance_ppm %||% 10
    ion_mode     <- mummi_cfg$ionization_mode %||% "pos_default"

    # Map ionization mode to mummichog flag
    mode_flag <- switch(
        tolower(as.character(ion_mode)),
        "positive" = "positive",
        "negative" = "negative",
        "pos_default"
    )

    out_id <- paste0("mcg_", org_label)

    # Build command using the mummichog CLI entry point
    mummi_bin <- file.path(dirname(python_bin), "mummichog")
    cmd <- sprintf(
        "%s -f %s -n %s -o %s -k %s --mode %s -u %d -p %d -c %s",
        shQuote(mummi_bin),
        shQuote(normalizePath(input_file)),
        shQuote(normalizePath(model_file)),
        shQuote(out_id),
        shQuote(normalizePath(work_dir)),
        mode_flag,
        ppm,
        n_perm,
        as.character(p_cutoff)
    )

    message("Running mummichog for ", organism, "...")
    message("  Command: ", cmd)

    result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    message("  mummichog completed")

    # Parse output — mummichog creates a timestamped directory: {timestamp}.{out_id}
    # Find the most recent output directory matching the out_id pattern
    result_dirs <- list.dirs(work_dir, recursive = FALSE)
    result_dirs <- result_dirs[grepl(paste0("\\.", out_id, "$"), result_dirs)]
    if (length(result_dirs) == 0) {
        # Fallback: look for any directory with tables/
        result_dirs <- list.dirs(work_dir, recursive = FALSE)
        result_dirs <- result_dirs[file.exists(file.path(result_dirs, "tables"))]
    }
    if (length(result_dirs) == 0) {
        warning("No mummichog output directory found for ", organism, " in ", work_dir)
        return(NULL)
    }
    # Use the most recent one
    results_dir <- result_dirs[order(file.info(result_dirs)$mtime, decreasing = TRUE)][1]

    pathway_tsv <- list.files(file.path(results_dir, "tables"),
                              pattern = "mcg_pathwayanalysis.*\\.tsv$",
                              full.names = TRUE)

    if (length(pathway_tsv) == 0) {
        warning("No mummichog pathway results found for ", organism)
        return(NULL)
    }

    pathway_table <- read.delim(pathway_tsv[1], stringsAsFactors = FALSE,
                                check.names = FALSE)

    # Standardize column names
    colnames(pathway_table) <- gsub("^p-value$", "p_value",
                                    gsub(" ", "_", colnames(pathway_table)))
    if (!"p_value" %in% colnames(pathway_table) && "p.value" %in% colnames(pathway_table))
        colnames(pathway_table)[colnames(pathway_table) == "p.value"] <- "p_value"

    # Add FDR
    if ("p_value" %in% colnames(pathway_table)) {
        pathway_table$FDR <- p.adjust(as.numeric(pathway_table$p_value), method = "BH")
    }

    # Add organism column
    pathway_table$organism <- organism

    # Sort by p-value
    if ("p_value" %in% colnames(pathway_table)) {
        pathway_table <- pathway_table[order(as.numeric(pathway_table$p_value)), ]
    }

    list(
        table    = pathway_table,
        method   = "mummichog",
        organism = organism,
        out_dir  = results_dir
    )
}


# ==== ORCHESTRATOR ============================================================

#' Run mummichog enrichment for all configured organisms
#'
#' @param pre     Preprocessing results (expr_work, meta, row_data).
#' @param de_res  DE results list (from mod_metabolomics_de).
#' @param config  Full pipeline config.
#' @param out_dir Output directory for the metabolomics mode.
#' @return list(results = named list per organism, organisms = character vector)
run_mummichog_all <- function(pre, de_res, config, out_dir) {
    mummi_cfg <- config$modes$metabolomics$enrichment$mummichog %||% list()

    if (!isTRUE(mummi_cfg$enabled)) {
        message("mummichog: disabled — skipping")
        return(NULL)
    }

    organisms <- mummi_cfg$organisms
    # Fallback: derive KEGG code from metabolomics.organism if mummichog organisms not set
    if (is.null(organisms) || length(organisms) == 0) {
        metab_organism <- config$modes$metabolomics$organism %||% NULL
        kegg_map <- c("Homo sapiens" = "hsa", "Mus musculus" = "mmu",
                      "Rattus norvegicus" = "rno", "Danio rerio" = "dre",
                      "Caenorhabditis elegans" = "cel",
                      "Drosophila melanogaster" = "dme",
                      "Saccharomyces cerevisiae" = "sce",
                      "Apis mellifera" = "ame",
                      "Escherichia coli" = "eco")
        if (!is.null(metab_organism) && metab_organism %in% names(kegg_map)) {
            organisms <- kegg_map[[metab_organism]]
            message("mummichog: using KEGG code '", organisms, "' from organism '", metab_organism, "'")
        } else {
            organisms <- c("hsa")
            message("mummichog: no organisms configured, defaulting to 'hsa' (human)")
        }
    }
    if (is.list(organisms)) organisms <- unlist(organisms)

    mummi_dir <- file.path(out_dir, "mummichog_enrichment")
    dir.create(mummi_dir, recursive = TRUE, showWarnings = FALSE)

    # Ensure Python env
    project_dir <- config$project$dir %||% "."
    python_bin <- ensure_mummichog_env(project_dir)

    # Extract DE table from mod_metabolomics_de() output
    # Structure: de_res$de_tables (named list), de_res$summary_df
    de_table <- NULL
    if (!is.null(de_res$de_tables) && length(de_res$de_tables) > 0) {
        de_table <- de_res$de_tables[[1]]
    } else if (!is.null(de_res$summary_df)) {
        de_table <- de_res$summary_df
    } else if (is.data.frame(de_res)) {
        de_table <- de_res
    }

    if (is.null(de_table)) {
        warning("mummichog: could not extract DE results table")
        return(NULL)
    }

    # Prepare input file (shared across organisms)
    input_file <- prepare_mummichog_input(
        de_res    = de_table,
        row_data  = pre$row_data,
        config    = config,
        output_dir = mummi_dir
    )

    # Build organism models
    model_dir <- file.path(mummi_dir, "models")
    model_files <- list()

    # Identify base organisms (not combined)
    base_orgs <- organisms[!grepl("\\+", organisms)]
    for (org in base_orgs) {
        model_files[[org]] <- build_kegg_organism_model(
            organism   = org,
            output_dir = model_dir,
            verbose    = TRUE
        )
    }

    # Build combined models
    combined_orgs <- organisms[grepl("\\+", organisms)]
    for (combo in combined_orgs) {
        parts <- trimws(strsplit(combo, "\\+")[[1]])
        combo_label <- gsub("\\+", "_", combo)
        # Ensure base models exist
        missing <- parts[!parts %in% names(model_files)]
        for (org in missing) {
            model_files[[org]] <- build_kegg_organism_model(
                organism = org, output_dir = model_dir, verbose = TRUE
            )
        }
        model_files[[combo]] <- build_combined_model(
            model_files = unname(model_files[parts]),
            output_dir  = model_dir,
            combined_id = combo_label
        )
    }

    # Run mummichog per organism
    results <- list()
    for (org in organisms) {
        org_label <- gsub("\\+", "_", org)
        mf <- model_files[[org]]
        if (is.null(mf) || !file.exists(mf)) {
            warning("mummichog: model not available for ", org, " — skipping")
            next
        }

        res <- tryCatch(
            run_mummichog_organism(
                input_file = input_file,
                model_file = mf,
                output_dir = mummi_dir,
                organism   = org,
                config     = config,
                python_bin = python_bin
            ),
            error = function(e) {
                warning("mummichog (", org, ") failed: ", e$message)
                NULL
            }
        )
        results[[org_label]] <- res
    }

    list(
        results   = results,
        organisms = organisms
    )
}


# ==== PLOTTING ================================================================

#' Barplot for mummichog pathway results
#'
#' @param table   data.frame from mummichog pathway results.
#' @param top_n   Number of top pathways to show.
#' @param title   Plot title.
#' @return ggplot object.
plot_mummichog_barplot <- function(table, top_n = 20, title = "Mummichog Pathway Enrichment") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
    if (is.null(table) || nrow(table) == 0) return(NULL)

    # Use p_value or p-value column
    pval_col <- intersect(c("p_value", "p-value"), colnames(table))
    if (length(pval_col) == 0) return(NULL)

    df <- table
    df$pval <- as.numeric(df[[pval_col[1]]])
    df$neg_log10p <- -log10(pmax(df$pval, 1e-10))
    df <- utils::head(df[order(df$pval), ], top_n)

    df$pathway <- factor(df$pathway, levels = rev(df$pathway))

    ggplot2::ggplot(df, ggplot2::aes(x = neg_log10p, y = pathway)) +
        ggplot2::geom_col(fill = "#2980B9", alpha = 0.85) +
        ggplot2::labs(x = expression(-log[10](p-value)), y = NULL,
                      title = title) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


#' Lollipop plot for mummichog pathway results
#'
#' @param table   data.frame from mummichog pathway results.
#' @param top_n   Number of top pathways.
#' @param title   Plot title.
#' @return ggplot object.
plot_mummichog_lollipop <- function(table, top_n = 20,
                                    title = "Mummichog Pathway Enrichment") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
    if (is.null(table) || nrow(table) == 0) return(NULL)

    pval_col <- intersect(c("p_value", "p-value"), colnames(table))
    if (length(pval_col) == 0) return(NULL)

    size_col <- intersect(c("overlap_size"), colnames(table))

    df <- table
    df$pval <- as.numeric(df[[pval_col[1]]])
    df$neg_log10p <- -log10(pmax(df$pval, 1e-10))
    if (length(size_col) > 0) {
        df$hits <- as.numeric(df[[size_col[1]]])
    } else {
        df$hits <- 1
    }
    df <- utils::head(df[order(df$pval), ], top_n)
    df$pathway <- factor(df$pathway, levels = rev(df$pathway))

    ggplot2::ggplot(df, ggplot2::aes(x = neg_log10p, y = pathway)) +
        ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = pathway),
                              color = "#7F8C8D") +
        ggplot2::geom_point(ggplot2::aes(size = hits), color = "#E74C3C",
                            alpha = 0.8) +
        ggplot2::scale_size_continuous(range = c(3, 10), name = "Overlap size") +
        ggplot2::labs(x = expression(-log[10](p-value)), y = NULL,
                      title = title) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


#' Cross-organism comparison heatmap
#'
#' Shows which pathways are significant across different organism models.
#'
#' @param results_list Named list of mummichog result lists (one per organism).
#' @param fdr_cutoff   FDR threshold for marking significance.
#' @return ggplot object.
plot_mummichog_organism_heatmap <- function(results_list, fdr_cutoff = 0.1) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    # Collect all pathways and their p-values per organism
    all_data <- list()
    for (org in names(results_list)) {
        res <- results_list[[org]]
        if (is.null(res) || is.null(res$table)) next
        tbl <- res$table
        pval_col <- intersect(c("p_value", "p-value"), colnames(tbl))
        if (length(pval_col) == 0) next

        # Strip organism prefix from combined model pathway names
        pathway_names <- tbl$pathway
        pathway_names <- sub("^\\[(AME|ECO)\\] ", "", pathway_names)

        all_data[[org]] <- data.frame(
            pathway  = pathway_names,
            organism = org,
            pval     = as.numeric(tbl[[pval_col[1]]]),
            stringsAsFactors = FALSE
        )
    }

    if (length(all_data) == 0) return(NULL)
    combined <- do.call(rbind, all_data)
    combined$neg_log10p <- -log10(pmax(combined$pval, 1e-10))

    # Keep only pathways significant in at least one organism
    sig_pathways <- unique(combined$pathway[combined$pval < fdr_cutoff])
    if (length(sig_pathways) == 0) {
        # Fall back to top 20 by min p-value
        pw_min_p <- tapply(combined$pval, combined$pathway, min)
        sig_pathways <- names(sort(pw_min_p))[seq_len(min(20, length(pw_min_p)))]
    }

    combined <- combined[combined$pathway %in% sig_pathways, ]

    ggplot2::ggplot(combined,
                    ggplot2::aes(x = organism, y = pathway, fill = neg_log10p)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient(low = "#EBF5FB", high = "#1A5276",
                                     name = expression(-log[10](p))) +
        ggplot2::labs(x = "Organism model", y = NULL,
                      title = "Mummichog: Cross-Organism Pathway Comparison") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title = ggplot2::element_text(face = "bold")
        )
}
