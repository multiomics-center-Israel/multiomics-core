# What happens when one omics layer's contrasts take different enrichment
# branches.
#
# run_kegg_enrichment_for_omics() picks its method per contrast, on the number of
# mapped features in that contrast:
#
#     use_method <- method
#     if (nrow(de_mapped) < 500) use_method <- "ora"
#
# So a layer whose contrasts straddle that line sends one through GSEA and the
# other through ORA, and the two return legitimately different columns -- ORA
# carries Fold_enrichment and Count, GSEA carries NES and core_enrichment.
#
# Combining those with rbind() aborts. The failure is quiet rather than loud:
# the abort happens outside the per-contrast tryCatch, the caller catches it and
# downgrades it to a warning, and the layer comes back NULL -- which the branch
# below reports as "No enriched pathways found", indistinguishable from a layer
# that genuinely enriched nothing. The per-omics CSV is written after the
# combine, so it is never written either.
#
# The contract these tests hold is therefore about the layer surviving, not
# about the combining helper: rows from both contrasts present, each method's
# own columns intact and NA on the other method's rows, and the layer not
# collapsing to the no-enrichment path.
#
# Everything here is synthetic. No KEGG, fgsea or clusterProfiler code runs --
# the two enrichment entry points are stubbed, because what is under test is
# what the function does with their results, not what they compute.

# The enrichment functions are sourced into the global environment by
# tar_source(), not exported from a package, so with_mocked_bindings() has no
# namespace to rebind and errors with "No packages loaded with pkgload". Assign
# into the function's own environment and restore on exit instead.
local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(run_kegg_enrichment_for_omics)
    nms <- names(stubs)

    had <- vapply(nms, exists, logical(1), envir = target, inherits = FALSE)
    old <- lapply(nms[had], get, envir = target, inherits = FALSE)
    names(old) <- nms[had]

    withr::defer({
        for (nm in nms) {
            if (nm %in% names(old)) {
                assign(nm, old[[nm]], envir = target)
            } else if (exists(nm, envir = target, inherits = FALSE)) {
                rm(list = nm, envir = target)
            }
        }
    }, envir = env)

    for (nm in nms) assign(nm, stubs[[nm]], envir = target)
    invisible(NULL)
}

# Contrast A clears the 500-feature threshold and takes the GSEA branch;
# contrast B sits under it and is forced to ORA. That is the whole point of the
# fixture, so the counts are chosen to straddle the threshold rather than to
# resemble any real experiment.
mixed_de_fixture <- function(n_a = 520L, n_b = 40L) {
    mk <- function(prefix, n) {
        data.frame(feature_id = paste0(prefix, "_", seq_len(n)),
                   log2FC = rep(1, n), pvalue = rep(0.01, n),
                   stringsAsFactors = FALSE)
    }
    list(contrast_A = mk("A", n_a), contrast_B = mk("B", n_b))
}

id_map_fixture <- function(de_tables) {
    ids <- unlist(lapply(de_tables, function(d) d$feature_id), use.names = FALSE)
    data.frame(feature_id = ids,
               ENTREZID = paste0("E", seq_along(ids)),
               stringsAsFactors = FALSE)
}

gsea_shaped <- function(n = 3L) {
    data.frame(ID = paste0("tst0000", seq_len(n)),
               Description = paste("GSEA pathway", seq_len(n)),
               NES = seq(2, by = -0.1, length.out = n),
               core_enrichment = rep("g1/g2", n),
               pvalue = rep(0.001, n), padj = rep(0.01, n),
               stringsAsFactors = FALSE)
}

ora_shaped <- function(n = 2L) {
    data.frame(ID = paste0("tst1000", seq_len(n)),
               Description = paste("ORA pathway", seq_len(n)),
               Fold_enrichment = seq(3, by = -0.5, length.out = n),
               Count = rep(7L, n),
               pvalue = rep(0.002, n), padj = rep(0.02, n),
               stringsAsFactors = FALSE)
}

# Drives the real function with everything upstream of the branch stubbed out.
run_mixed_layer <- function(out_dir, env = parent.frame()) {
    de_tables <- mixed_de_fixture()
    id_map <- id_map_fixture(de_tables)
    conv <- setNames(paste0("tst:", id_map$ENTREZID), id_map$ENTREZID)

    local_stubs(list(
        get_kegg_organism = function(...) "tst",
        get_organism_db   = function(...) "OrgDb.stub",
        extract_de_tables = function(...) de_tables,
        map_feature_ids_to_entrez = function(...) id_map,
        run_gsea_kegg = function(...) gsea_shaped(),
        run_ora_kegg  = function(...) ora_shaped()
    ), env = env)

    run_kegg_enrichment_for_omics(
        de_data = list(tables = de_tables),
        omics_type = "transcriptomics",
        harmonization_res = list(),
        organism = "test organism",
        config = list(modes = list(multiomics = list(
            enrichment = list(methods = "gsea")))),
        out_dir = out_dir,
        kegg_conv_cache = conv
    )
}


test_that("a layer whose contrasts take different branches is not lost", {
    out_dir <- withr::local_tempdir()

    res <- run_mixed_layer(out_dir)

    # The failure this guards against is silent: NULL here reaches the caller as
    # a warning and prints as "No enriched pathways found".
    expect_false(is.null(res))
    expect_gt(nrow(res), 0)
})

test_that("both contrasts contribute their rows", {
    out_dir <- withr::local_tempdir()

    res <- run_mixed_layer(out_dir)

    expect_setequal(unique(res$contrast), c("contrast_A", "contrast_B"))
    expect_equal(nrow(res), nrow(gsea_shaped()) + nrow(ora_shaped()))
    expect_setequal(unique(res$omics), "transcriptomics")
})

test_that("each method's own columns survive and are NA on the other's rows", {
    out_dir <- withr::local_tempdir()

    res <- run_mixed_layer(out_dir)

    # Both schemas are present in the combined table rather than one being
    # dropped to make the shapes agree.
    expect_true(all(c("NES", "core_enrichment", "Fold_enrichment", "Count")
                    %in% names(res)))

    gsea_rows <- res$contrast == "contrast_A"
    ora_rows  <- res$contrast == "contrast_B"

    # Which branch each contrast took, read off the columns it filled: the
    # larger contrast went to GSEA, the smaller was forced to ORA.
    expect_true(all(!is.na(res$NES[gsea_rows])))
    expect_true(all(is.na(res$NES[ora_rows])))
    expect_true(all(!is.na(res$Fold_enrichment[ora_rows])))
    expect_true(all(is.na(res$Fold_enrichment[gsea_rows])))

    # Columns both methods share are populated throughout.
    expect_true(all(!is.na(res$pvalue)))
    expect_true(all(!is.na(res$padj)))
})

test_that("the per-omics CSV is written, not skipped by an aborted combine", {
    out_dir <- withr::local_tempdir()

    res <- run_mixed_layer(out_dir)

    csv <- file.path(out_dir, "transcriptomics_kegg_enrichment.csv")
    expect_true(file.exists(csv))
    expect_equal(nrow(read.csv(csv, stringsAsFactors = FALSE)), nrow(res))
})

test_that("a layer whose contrasts agree on method is unchanged", {
    # The mixed case is the regression; this is the case that already worked,
    # asserted so the fix cannot quietly alter it.
    out_dir <- withr::local_tempdir()
    de_tables <- mixed_de_fixture(n_a = 520L, n_b = 600L)
    id_map <- id_map_fixture(de_tables)
    conv <- setNames(paste0("tst:", id_map$ENTREZID), id_map$ENTREZID)

    local_stubs(list(
        get_kegg_organism = function(...) "tst",
        get_organism_db   = function(...) "OrgDb.stub",
        extract_de_tables = function(...) de_tables,
        map_feature_ids_to_entrez = function(...) id_map,
        run_gsea_kegg = function(...) gsea_shaped(),
        run_ora_kegg  = function(...) stop("ORA branch should not be reached here")
    ))

    res <- run_kegg_enrichment_for_omics(
        de_data = list(tables = de_tables),
        omics_type = "transcriptomics",
        harmonization_res = list(),
        organism = "test organism",
        config = list(modes = list(multiomics = list(
            enrichment = list(methods = "gsea")))),
        out_dir = out_dir,
        kegg_conv_cache = conv
    )

    expect_equal(nrow(res), 2L * nrow(gsea_shaped()))
    expect_false(any(c("Fold_enrichment", "Count") %in% names(res)))
    expect_true(all(!is.na(res$NES)))
})


# ---- the helper the fix now depends on -------------------------------------

test_that(".rbind_fill keeps both schemas instead of aborting", {
    # No coverage existed for this helper, and the combine above now rests on
    # it, so the NA-filling behaviour is worth pinning directly.
    combined <- .rbind_fill(list(gsea_shaped(1L), ora_shaped(1L)))

    expect_equal(nrow(combined), 2L)
    expect_true(all(c("NES", "core_enrichment", "Fold_enrichment", "Count")
                    %in% names(combined)))
    expect_true(is.na(combined$NES[2]))
    expect_true(is.na(combined$Fold_enrichment[1]))
})
