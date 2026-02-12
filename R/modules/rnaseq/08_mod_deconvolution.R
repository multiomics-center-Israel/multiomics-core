#' Module wrapper: RNA-seq Cell Type Deconvolution
#'
#' Organism guard + delegates to domain.
#' Skips silently for non-human/mouse organisms.
#'
#' @param pre Preprocessed data list
#' @param de_res DE results list
#' @param config Full config
#' @param out_dir Output directory
#' @return List with deconvolution results, or NULL-like skip result
mod_rnaseq_deconvolution <- function(pre, de_res, config, out_dir) {
    rna_cfg <- config$modes$rna

    dc_cfg <- rna_cfg$deconvolution %||% list(enabled = FALSE)

    if (!isTRUE(dc_cfg$enabled)) {
        message("[mod_deconvolution] Deconvolution disabled in config. Skipping.")
        return(list(status = "skipped", reason = "Disabled in config"))
    }

    # Get organism from annotation config
    organism <- rna_cfg$annotation$organism %||% "auto"

    # If auto, try to detect
    if (organism == "auto") {
        gene_ids <- rownames(pre$expr_filt)
        if (!is.null(gene_ids)) {
            org_detect <- tryCatch(detect_organism_from_ids(gene_ids), error = function(e) NULL)
            if (!is.null(org_detect) && org_detect$organism != "unknown") {
                organism <- org_detect$organism
            }
        }
    }

    organism_norm <- tryCatch(normalize_organism_name(organism), error = function(e) organism)

    if (!is_deconvolution_supported(organism_norm)) {
        message("[mod_deconvolution] Organism '", organism_norm, "' not supported. Skipping.")
        return(list(
            status = "skipped",
            reason = paste0("Organism '", organism_norm, "' not supported (human/mouse only)")
        ))
    }

    message("[mod_deconvolution] Running deconvolution for: ", organism_norm)
    run_deconvolution(pre, de_res, rna_cfg, out_dir)
}
