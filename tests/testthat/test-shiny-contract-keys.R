# tests/testthat/test-shiny-contract-keys.R
#
# Canonical payload keys must stay present (as NULL) in every omics builder
# even when the value to store is NULL: `payload$key <- NULL` would delete the
# key from the list, so the builders use `payload["key"] <- list(value)`.
# Also pins that contrasts, clust_patterns_list and clust_heatmap_partition_fig
# are registered in the contract. Fixtures are synthetic and deterministic.

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

build_min_payload <- function(omics, pca_res = NULL, clustering_res = NULL,
                              contrasts = NULL) {
    samples <- paste0("s", 1:4)
    feats   <- paste0("f", 1:6)
    mat <- matrix(as.numeric(1:24), nrow = 6, ncol = 4,
                  dimnames = list(feats, samples))
    meta <- data.frame(sample_id = samples, grp = c("A", "A", "B", "B"),
                       row.names = samples, stringsAsFactors = FALSE)
    inputs <- list(contrasts = contrasts)

    switch(omics,
        rnaseq = {
            pre <- list(expr_filt = mat, expr_work = mat, meta = meta,
                        row_data = data.frame(gene_id = feats,
                                              stringsAsFactors = FALSE))
            config <- list(modes = list(rna = list(
                de = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
                normalization = list(method = "TMMlogCPM"),
                effects = list(color = list("grp")),
                id_columns = list(gene_id = "gene_id"))))
            quiet(build_shiny_payload_rnaseq(
                pre = pre, inputs = inputs, config = config,
                pca_res = pca_res, clustering_res = clustering_res))
        },
        proteomics = {
            pre <- list(expr_filt = mat, expr_imp_single = mat, meta = meta,
                        row_data = data.frame(protein_id = feats,
                                              stringsAsFactors = FALSE))
            config <- list(modes = list(proteomics = list(
                de = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
                normalization = list(method = "none"),
                effects = list(color = list("grp")),
                id_columns = list(protein_id = "protein_id"))))
            quiet(build_shiny_payload_proteomics(
                pre = pre, inputs = inputs, config = config,
                pca_res = pca_res, clustering_res = clustering_res))
        },
        metabolomics = {
            rd <- data.frame(feature_id = feats, row.names = feats,
                             stringsAsFactors = FALSE)
            pre <- list(expr_raw = mat, expr_filt = mat, expr_work = mat,
                        meta = meta, row_data = rd)
            config <- list(modes = list(metabolomics = list(
                de = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
                preprocessing = list(chosen_norm = "pqn", transform = "log2",
                                     scaling = "none"),
                effects = list(samples = "sample_id", color = "grp"))))
            quiet(build_shiny_payload_metabolomics(
                pre = pre, inputs = inputs, config = config,
                pca_res = pca_res, clustering_res = clustering_res,
                include_legacy = FALSE))
        }
    )
}

omics_types <- c("rnaseq", "proteomics", "metabolomics")

test_that("contrasts, clust_patterns_list and clust_heatmap_partition_fig are registered and optional", {
    new_keys <- c("contrasts", "clust_patterns_list", "clust_heatmap_partition_fig")
    expect_true(all(new_keys %in% get_canonical_keys()))
    expect_false(any(new_keys %in% get_required_keys()))
    expect_true(all(new_keys %in% names(init_shiny_payload("rnaseq"))))
})

test_that("every builder keeps every canonical key present when optional inputs are absent", {
    for (omics in omics_types) {
        p <- build_min_payload(omics)
        expect_true(all(get_canonical_keys() %in% names(p)), info = omics)
        expect_null(p$contrasts)
    }
})

test_that("canonical keys stay present when pca/clustering results carry none of the optional objects", {
    pca_res        <- list(plots = list())
    clustering_res <- list(objects = list(patterns = data.frame(pattern = "p1")))

    for (omics in omics_types) {
        p <- build_min_payload(omics, pca_res = pca_res,
                               clustering_res = clustering_res)
        expect_true(all(get_canonical_keys() %in% names(p)), info = omics)
        expect_null(p$clust_patterns_list)
    }
})

test_that("contrasts, clust_patterns_list and clust_heatmap_partition_fig are carried through", {
    contrasts      <- data.frame(Contrast_name = "A_vs_B", stringsAsFactors = FALSE)
    clustering_res <- list(
        objects = list(patterns      = data.frame(pattern = "p1"),
                       patterns_list = c("010", "101")),
        plots   = list(partition_heatmap = list(gtable = "G"))
    )

    for (omics in omics_types) {
        p <- build_min_payload(omics, clustering_res = clustering_res,
                               contrasts = contrasts)
        expect_identical(p$contrasts, contrasts, info = omics)
        expect_identical(p$clust_patterns_list, c("010", "101"), info = omics)
        expect_identical(p$clust_heatmap_partition_fig, "G", info = omics)
    }
})
