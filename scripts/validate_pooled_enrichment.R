#!/usr/bin/env Rscript
# =============================================================================
# Validation of the Rubin-pooled ranking and its enrichment (proteomics)
#
# Standalone; not sourced by R/. Reads a finished proteomics run from its
# targets store and writes comparison tables and plots to a directory you name.
# Nothing in the pipeline reads these outputs.
#
# Usage (from the repository root, after a proteomics run with de.method limma):
#   Rscript scripts/validate_pooled_enrichment.R --store <targets store> --out <dir> [--k 100]
#
# What it compares, per contrast:
#   1. Protein rankings: legacy fgsea statistic (sign(log2FC) * -log10 p),
#      log2FC, the naive mean of per-draw moderated t, pooled z option A (the
#      pipeline method) and two option B variants. Spearman correlation and
#      top/bottom-k overlap.
#   2. Leave-one-imputation-out stability of option A.
#   3. Sensitivity to proteins with frac_imputed > 0.5.
#   4. Pathway overlap: legacy fgsea, pooled fgsea, cameraPR (FDR 0.05 / 0.10),
#      with sign agreement on shared sets.
#   5. cameraPR at inter.gene.cor 0, 0.01, 0.05, and the per-set inter-gene
#      correlation estimated by limma::interGeneCorrelation, averaged over draws.
#   6. p-value histograms and significance against set size.
#
# Option B mirrors mi4p (Chion, Carapito & Bertrand 2022, PLoS Comput Biol
# 18:e1010420): pool the unmoderated variances, sigma_d^2 * stdev.unscaled_d^2,
# with Rubin's rules, then moderate the pooled variance with limma::squeezeVar
# on the per-observation scale using the complete-data residual df, and refer t
# to df.residual + df.prior (no Barnard-Rubin df). mi4p also projects the K x K
# coefficient covariance by taking a maximum over groups; here the contrast
# variance is used directly, which is a deliberate difference. A second variant
# squeezes with the Barnard-Rubin df instead, so the two can be compared.
# =============================================================================

suppressPackageStartupMessages({
    library(targets)
    library(limma)
    library(ggplot2)
})

parse_args <- function(args) {
    out <- list(store = NULL, out = NULL, k = 100L)
    i <- 1L
    while (i <= length(args)) {
        key <- sub("^--", "", args[i])
        if (!key %in% names(out) || i == length(args)) stop("Unknown or incomplete argument: ", args[i])
        out[[key]] <- args[i + 1L]
        i <- i + 2L
    }
    if (is.null(out$store) || is.null(out$out)) {
        stop("Usage: Rscript scripts/validate_pooled_enrichment.R --store <targets store> --out <dir> [--k 100]")
    }
    out$k <- as.integer(out$k)
    out
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opts$out, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# Same load order as _targets.R, without the pipeline layer.
for (layer in c("core", "services", "domain", "modules")) {
    for (f in sort(list.files(file.path("R", layer), pattern = "\\.R$", full.names = TRUE,
                              recursive = TRUE))) {
        source(f)
    }
}

config <- tar_read(config, store = opts$store)
de_res <- tar_read(prot_de_res, store = opts$store)
pre    <- tar_read(prot_pre, store = opts$store)
if (!identical(de_res$method, "limma")) stop("This run used de.method '", de_res$method, "'; pooling needs 'limma'.")

runs <- de_res$runs
m <- length(runs)
observed <- !is.na(pre$expr_filt)
pw_cfg <- config$modes$proteomics$pathway %||% list()
min_size <- pw_cfg$min_size %||% 10
max_size <- pw_cfg$max_size %||% 500
seed <- config$params$seed %||% 1L
k <- opts$k

with_ranking <- function(cfg, ranking) {
    cfg$modes$proteomics$pathway$gsea_ranking <- ranking
    cfg
}

# ---- Option B: pool unmoderated variances, moderate after pooling ----------
option_b <- function(runs, cn, squeeze_df = c("residual", "barnard_rubin")) {
    squeeze_df <- match.arg(squeeze_df)
    fits <- lapply(runs, function(r) r$fit2)
    Q <- sapply(fits, function(f) f$coefficients[, cn])
    su2 <- sapply(fits, function(f) f$stdev.unscaled[, cn]^2)
    U <- sapply(fits, function(f) f$sigma^2) * su2
    df_res <- fits[[1]]$df.residual
    rp <- rubin_pool_estimates(Q, U, df_com = df_res)
    su2_bar <- rowMeans(su2)
    v <- rp$T / su2_bar
    df_sq <- if (squeeze_df == "residual") df_res else rp$df
    sq <- limma::squeezeVar(v, df = df_sq)
    t_b <- rp$estimate / sqrt(sq$var.post * su2_bar)
    limma::zscoreT(t_b, df = df_sq + sq$df.prior)
}

# Overlap of the top and bottom k of two rankings, on the features finite in
# both. k is cut to the number of such features and returned, so a row says
# how many it compared (top-100 of 20 features is all 20, a trivial 1.0).
topk_overlap <- function(a, b, k) {
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    k <- min(k, length(a))
    if (k == 0) return(c(top = NA_real_, bottom = NA_real_, k = 0))
    top <- length(intersect(order(-a)[1:k], order(-b)[1:k])) / k
    bot <- length(intersect(order(a)[1:k], order(b)[1:k])) / k
    c(top = top, bottom = bot, k = k)
}

# Spearman correlation of two rankings, NA when fewer than two features are
# finite in both (an empty sensitivity subset would otherwise make cor() stop
# the script before it writes anything).
spearman_finite <- function(a, b) {
    if (sum(is.finite(a) & is.finite(b)) < 2) return(NA_real_)
    stats::cor(a, b, method = "spearman", use = "complete.obs")
}

pairwise_compare <- function(df, k) {
    cols <- setdiff(names(df), "FeatureID")
    pairs <- utils::combn(cols, 2)
    do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
        a <- df[[pairs[1, j]]]; b <- df[[pairs[2, j]]]
        ov <- topk_overlap(a, b, k)
        data.frame(ranking_a = pairs[1, j], ranking_b = pairs[2, j],
                   spearman = spearman_finite(a, b),
                   top_k_overlap = ov[["top"]], bottom_k_overlap = ov[["bottom"]],
                   k_requested = k, k = ov[["k"]])
    }))
}

sig_sets <- function(tab, fdr) {
    if (is.null(tab) || nrow(tab) == 0) return(character(0))
    p <- if ("padj" %in% names(tab)) tab$padj else tab$FDR
    tab$pathway[!is.na(p) & p < fdr]
}

set_sign <- function(tab) {
    if (is.null(tab) || nrow(tab) == 0) return(setNames(numeric(0), character(0)))
    s <- if ("NES" %in% names(tab)) sign(tab$NES) else ifelse(tab$Direction == "Up", 1, -1)
    setNames(s, tab$pathway)
}

# ---- Rankings ---------------------------------------------------------------
pooled_tables <- build_pooled_ranking_tables(de_res, pre, config)
contrasts <- colnames(runs[[1]]$fit2$coefficients)
rank_cmp <- loo <- sens <- list()

for (cn in contrasts) {
    key <- normalize_contrast_name(cn)
    pa <- pool_limma_imputations(runs, cn, observed = observed)
    leg <- extract_de_table_for_pathway(de_res$summary_df, cn, with_ranking(config, "stat"))
    lfc <- extract_de_table_for_pathway(de_res$summary_df, cn, with_ranking(config, "lfc"))
    rk <- data.frame(
        FeatureID   = pa$FeatureID,
        legacy_stat = leg$stat[match(pa$FeatureID, leg$FeatureID)],
        lfc         = lfc$stat[match(pa$FeatureID, lfc$FeatureID)],
        mean_t      = rowMeans(sapply(runs, function(r) r$fit2$t[, cn])),
        z_A         = pa$z_pooled,
        z_B_mi4p    = option_b(runs, cn, "residual"),
        z_B_br      = option_b(runs, cn, "barnard_rubin")
    )
    utils::write.csv(cbind(rk, frac_imputed = pa$frac_imputed, fmi = pa$fmi),
                     file.path(opts$out, sprintf("rankings_%s.csv", key)), row.names = FALSE)
    rank_cmp[[key]] <- cbind(contrast = key, pairwise_compare(rk, k))

    if (m > 1) {
        loo[[key]] <- do.call(rbind, lapply(seq_len(m), function(i) {
            z_i <- pool_limma_imputations(runs[-i], cn)$z_pooled
            ov <- topk_overlap(pa$z_pooled, z_i, k)
            data.frame(contrast = key, left_out = i,
                       spearman = spearman_finite(pa$z_pooled, z_i),
                       top_k_overlap = ov[["top"]], bottom_k_overlap = ov[["bottom"]],
                   k_requested = k, k = ov[["k"]])
        }))
    }

    keep <- is.finite(pa$frac_imputed) & pa$frac_imputed <= 0.5
    sens[[key]] <- cbind(contrast = key, n_kept = sum(keep), n_total = nrow(pa),
                         pairwise_compare(rk[keep, , drop = FALSE], k))
}
utils::write.csv(do.call(rbind, rank_cmp), file.path(opts$out, "ranking_comparison.csv"), row.names = FALSE)
if (length(loo)) utils::write.csv(do.call(rbind, loo), file.path(opts$out, "loo_stability.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, sens), file.path(opts$out, "sensitivity_frac_imputed_le_0.5.csv"), row.names = FALSE)

# ---- Pathways ---------------------------------------------------------------
cfg_p <- config$modes$proteomics
skip_remap <- isTRUE((cfg_p$annotation %||% list())$skip_annotation)
legacy_tables <- lapply(setNames(names(pooled_tables), names(pooled_tables)), function(key) {
    tbl <- extract_de_table_for_pathway(de_res$summary_df, key, config)
    if (skip_remap) tbl else map_proteins_to_gene_symbols(tbl, de_res$summary_df, config)
})
gs_info <- prepare_proteomics_gene_sets(legacy_tables, config)
if (is.null(gs_info)) stop("No gene sets loaded; cannot compare pathway results.")
gene_sets <- gs_info$gene_sets

run_pa <- function(tables, method, rho = 0.01) {
    run_pathway_analysis(tables, gene_sets, method = method, min_size = min_size,
                         max_size = max_size, seed = seed, camerapr_inter_gene_cor = rho)
}
legacy_res <- run_pa(legacy_tables, "fgsea")
pooled_res <- run_pa(pooled_tables, c("fgsea", "camerapr"))
sub_tables <- lapply(pooled_tables, function(t) t[is.finite(t$frac_imputed) & t$frac_imputed <= 0.5, , drop = FALSE])
sub_res <- run_pa(sub_tables, c("fgsea", "camerapr"))
rho_res <- lapply(c(`0` = 0, `0.01` = 0.01, `0.05` = 0.05), function(r) run_pa(pooled_tables, "camerapr", r))

overlap <- rho_tab <- igc <- sens_pw <- pvals <- list()
for (key in names(pooled_tables)) {
    for (db in names(gene_sets)) {
        lf <- legacy_res[[key]][[paste0(db, "_fgsea")]]
        pf <- pooled_res[[key]][[paste0(db, "_fgsea")]]
        pc <- pooled_res[[key]][[paste0(db, "_camerapr")]]
        for (fdr in c(0.05, 0.10)) {
            s <- list(legacy_fgsea = sig_sets(lf, fdr), pooled_fgsea = sig_sets(pf, fdr),
                      camerapr = sig_sets(pc, fdr))
            signs <- list(legacy_fgsea = set_sign(lf), pooled_fgsea = set_sign(pf), camerapr = set_sign(pc))
            pr <- utils::combn(names(s), 2)
            overlap[[length(overlap) + 1]] <- do.call(rbind, lapply(seq_len(ncol(pr)), function(j) {
                a <- pr[1, j]; b <- pr[2, j]
                shared <- intersect(s[[a]], s[[b]])
                data.frame(contrast = key, database = db, fdr = fdr, method_a = a, method_b = b,
                           n_a = length(s[[a]]), n_b = length(s[[b]]), n_shared = length(shared),
                           sign_agreement = if (length(shared)) mean(signs[[a]][shared] == signs[[b]][shared]) else NA_real_)
            }))
            sp <- list(full_fgsea = sig_sets(pf, fdr), sub_fgsea = sig_sets(sub_res[[key]][[paste0(db, "_fgsea")]], fdr),
                       full_camerapr = sig_sets(pc, fdr), sub_camerapr = sig_sets(sub_res[[key]][[paste0(db, "_camerapr")]], fdr))
            sens_pw[[length(sens_pw) + 1]] <- data.frame(
                contrast = key, database = db, fdr = fdr,
                n_full_fgsea = length(sp$full_fgsea), n_sub_fgsea = length(sp$sub_fgsea),
                n_shared_fgsea = length(intersect(sp$full_fgsea, sp$sub_fgsea)),
                n_full_camerapr = length(sp$full_camerapr), n_sub_camerapr = length(sp$sub_camerapr),
                n_shared_camerapr = length(intersect(sp$full_camerapr, sp$sub_camerapr)))
        }
        for (r in names(rho_res)) {
            t <- rho_res[[r]][[key]][[paste0(db, "_camerapr")]]
            if (!is.null(t) && nrow(t)) rho_tab[[length(rho_tab) + 1]] <- data.frame(
                contrast = key, database = db, inter_gene_cor = as.numeric(r),
                pathway = t$pathway, PValue = t$PValue, FDR = t$FDR, Direction = t$Direction)
        }

        # Per-set inter-gene correlation on the proteins that were ranked,
        # averaged over the imputation draws.
        tbl <- pooled_tables[[key]]
        ranks <- .build_fgsea_ranks(tbl)
        sets <- .filter_sets_to_ranks(gene_sets[[db]], ranks, min_size, max_size)
        if (length(sets)) {
            igc[[length(igc) + 1]] <- do.call(rbind, lapply(names(sets), function(sn) {
                prots <- tbl$ProteinID[match(sets[[sn]], tbl$FeatureID)]
                cors <- vapply(runs, function(r) {
                    limma::interGeneCorrelation(r$expr_imp[prots, , drop = FALSE], r$design)$correlation
                }, numeric(1))
                data.frame(contrast = key, database = db, pathway = sn, NGenes = length(prots),
                           mean_inter_gene_cor = mean(cors), sd_over_draws = if (m > 1) stats::sd(cors) else NA_real_)
            }))
        }

        for (nm in c("legacy_fgsea", "pooled_fgsea", "camerapr")) {
            t <- switch(nm, legacy_fgsea = lf, pooled_fgsea = pf, camerapr = pc)
            if (is.null(t) || !nrow(t)) next
            pvals[[length(pvals) + 1]] <- data.frame(
                contrast = key, database = db, method = nm,
                p = if ("pval" %in% names(t)) t$pval else t$PValue,
                size = if ("size" %in% names(t)) t$size else t$NGenes)
        }
    }
}
utils::write.csv(do.call(rbind, overlap), file.path(opts$out, "pathway_overlap.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, sens_pw), file.path(opts$out, "pathway_sensitivity_frac_imputed_le_0.5.csv"), row.names = FALSE)
if (length(rho_tab)) utils::write.csv(do.call(rbind, rho_tab), file.path(opts$out, "camerapr_inter_gene_cor.csv"), row.names = FALSE)
if (length(igc)) utils::write.csv(do.call(rbind, igc), file.path(opts$out, "estimated_inter_gene_cor.csv"), row.names = FALSE)

if (length(pvals)) {
    pv <- do.call(rbind, pvals)
    for (key in unique(pv$contrast)) {
        d <- pv[pv$contrast == key, ]
        h <- ggplot(d, aes(x = p)) + geom_histogram(bins = 20, boundary = 0) +
            facet_grid(database ~ method, scales = "free_y") +
            labs(title = paste("Set-level p-values:", key), x = "p-value", y = "Sets") +
            theme_minimal()
        ggsave(file.path(plot_dir, sprintf("pvalue_hist_%s.png", key)), h, width = 10, height = 7)
        s <- ggplot(d, aes(x = size, y = -log10(p))) + geom_point(alpha = 0.4, size = 1) +
            scale_x_log10() + facet_grid(database ~ method) +
            labs(title = paste("Significance against set size:", key),
                 x = "Set size (log scale)", y = "-log10(p-value)") +
            theme_minimal()
        ggsave(file.path(plot_dir, sprintf("pvalue_vs_size_%s.png", key)), s, width = 10, height = 7)
    }
}

message("Validation outputs written to ", normalizePath(opts$out))
