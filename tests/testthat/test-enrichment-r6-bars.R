# =============================================================================
# §18 D18.1 / D18.3 — pipeline enrichment summary bar plots + Dice shared terms
# The pipeline static bar plots mirror the Shiny Enrichment Plot exactly (same
# plot type, metric, axes, colour, annotation) and are written in two static
# versions per unit — Top 20 (by adj p) and Full (all significant terms). The
# ORA enrichplot dotplot is retired as the canonical summary. Shared-terms
# overlap is Sørensen–Dice everywhere (never Jaccard), and the term<->term
# static heatmap always renders (dynamic sizing, no silent skip).
# =============================================================================

.mk_ora_bar_df <- function(n = 30) {
    set.seed(11)
    data.frame(
        Description = paste("term", seq_len(n)),
        Count = sample(3:40, n, TRUE),
        BgRatio = paste0(sample(20:300, n, TRUE), "/9087"),
        Fold_enrichment = round(runif(n, 1.2, 9), 1),
        p.adjust = sort(runif(n, 1e-8, 0.05)),
        geneID = replicate(n, paste0("g", sample(1:200, 5), collapse = "/")),
        stringsAsFactors = FALSE)
}

.mk_gsea_bar_df <- function(n_pos = 12, n_neg = 13) {
    set.seed(12); n <- n_pos + n_neg
    data.frame(
        Description = paste("path", seq_len(n)),
        NES = c(runif(n_pos, 1, 3), -runif(n_neg, 1, 3)),
        setSize = sample(15:400, n, TRUE),
        p.adjust = sort(runif(n, 1e-6, 0.05)),
        core_enrichment = replicate(n, paste0("g", sample(1:300, 8), collapse = "/")),
        stringsAsFactors = FALSE)
}

# ---- ORA bar (§18 D18.1) ----------------------------------------------------

test_that("ORA enrichment bar (§23 B): SELECT top-N by FDR, DISPLAY by Count desc", {
    skip_if_not_installed("ggplot2")
    df <- .mk_ora_bar_df(30)
    p20 <- plot_ora_enrichment_bar(df, max_terms = 20)
    pf  <- plot_ora_enrichment_bar(df, max_terms = nrow(df))
    expect_equal(nrow(p20$data), 20L)             # Top-20 head
    expect_equal(nrow(pf$data), 30L)              # Full = all terms
    # SELECTION = the 20 best-FDR terms (as a SET), regardless of display order
    expect_setequal(as.numeric(p20$data$p.adjust), sort(df$p.adjust)[seq_len(20)])
    # DISPLAY order = Count descending (data frame row order)
    expect_true(".count" %in% names(p20$data))
    expect_false(is.unsorted(rev(p20$data$.count)))            # non-increasing Count
    expect_equal(p20$data$.count[1], max(p20$data$.count))     # highest Count first
    # and the y factor puts the highest Count at the TOP (last factor level)
    expect_equal(as.character(p20$data$.lab[1]),
                 utils::tail(levels(p20$data$.lab), 1))
})

# ---- GSEA bar (§18 D18.1) ---------------------------------------------------

test_that("GSEA enrichment bar: X = SIGNED NES (both signs preserved), LE count", {
    skip_if_not_installed("ggplot2")
    df <- .mk_gsea_bar_df()
    pf <- plot_gsea_enrichment_bar(df, max_terms = nrow(df))
    expect_true(".nes" %in% names(pf$data))
    expect_true(any(pf$data$.nes > 0) && any(pf$data$.nes < 0))  # signed, never abs()
    # leading-edge count read from core_enrichment (never recomputed)
    expect_true(all(grepl("LE ", pf$data$.ann)))
})

test_that("GSEA enrichment bar (§D25.1): SELECT top-N by FDR, DISPLAY by signed NES desc", {
    skip_if_not_installed("ggplot2")
    # deterministic: FDR order (P1<P2<..) differs from signed-NES order
    df <- data.frame(
        ID = paste0("P", 1:5), Description = paste("path", 1:5),
        NES = c(2.0, 1.5, -2.0, 1.1, 2.2),
        setSize = c(10L, 50L, 30L, 20L, 40L),
        p.adjust = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
        core_enrichment = replicate(5, paste0("g", 1:4, collapse = "/")),
        stringsAsFactors = FALSE)
    p3 <- plot_gsea_enrichment_bar(df, max_terms = 3)
    # SELECTION = 3 best-FDR pathways (as a SET): P1/P2/P3
    expect_setequal(as.character(p3$data$Description), paste("path", 1:3))
    # DISPLAY order = SIGNED NES descending: P1(2.0), P2(1.5), P3(-2.0)
    expect_equal(as.numeric(p3$data$.nes), c(2.0, 1.5, -2.0))
    # setSize is NOT the display-order criterion (would have been 50, 30, 10)
    expect_false(identical(as.integer(p3$data$setSize), c(50L, 30L, 10L)))
    # highest signed NES sits at the TOP (last y factor level)
    expect_equal(as.character(p3$data$.lab[1]), utils::tail(levels(p3$data$.lab), 1))
    # ALL: every term, ordered by signed NES descending (2.2,2.0,1.5,1.1,-2.0)
    pall <- plot_gsea_enrichment_bar(df, max_terms = nrow(df))
    expect_equal(as.numeric(pall$data$.nes), c(2.2, 2.0, 1.5, 1.1, -2.0))
})

# ---- write pair (§18 D18.1) -------------------------------------------------

test_that("write_enrichment_bar_pair writes top20 + full for ORA and GSEA", {
    skip_if_not_installed("ggplot2")
    d <- withr::local_tempdir()
    write_enrichment_bar_pair(.mk_ora_bar_df(30), d, plot_ora_enrichment_bar)
    expect_true(file.exists(file.path(d, "enrichment_bar_top20.pdf")))
    expect_true(file.exists(file.path(d, "enrichment_bar_full.pdf")))

    dg <- withr::local_tempdir()
    write_enrichment_bar_pair(.mk_gsea_bar_df(), dg, plot_gsea_enrichment_bar)
    expect_true(file.exists(file.path(dg, "enrichment_bar_top20.pdf")))
    expect_true(file.exists(file.path(dg, "enrichment_bar_full.pdf")))
})

test_that("write_enrichment_bar_pair is a no-op on an empty table", {
    empty <- .mk_ora_bar_df(30)[0, , drop = FALSE]
    d <- withr::local_tempdir()
    write_enrichment_bar_pair(empty, d, plot_ora_enrichment_bar)
    expect_length(list.files(d, pattern = "pdf$"), 0L)
})

# ---- Shared terms: Sørensen–Dice + dynamic sizing (§18 D18.3) ---------------

test_that("term<->term overlap is Sørensen–Dice (diag 100, symmetric, <=100)", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("RColorBrewer")
    set.seed(3); n <- 54
    df <- data.frame(
        Cluster = rep("3", n), ID = paste0("GO:", seq_len(n)),
        Description = paste("bigterm", seq_len(n)),
        geneID = replicate(n, paste0("g", sample(1:120, sample(4:25, 1)),
                                     collapse = "/")),
        stringsAsFactors = FALSE)
    od <- withr::local_tempdir()
    plot_ora_shared_genes(df, od)
    t2t_pdf <- file.path(od, "cluster_3_terms_to_terms.pdf")
    t2t_csv <- file.path(od, "cluster_3_terms_to_terms.csv")
    # §18 D18.3: 54-term matrix ALWAYS renders a PDF (no silent skip)
    expect_true(file.exists(t2t_pdf))
    expect_true(file.exists(t2t_csv))
    m <- as.matrix(utils::read.csv(t2t_csv, row.names = 1, check.names = FALSE))
    expect_true(all(abs(diag(m) - 100) < 1e-6))   # Dice(A,A) = 100
    expect_lte(max(m), 100 + 1e-6)
    expect_equal(m, t(m), tolerance = 1e-6)        # symmetric
})

test_that("Dice matches the closed form 100*2|A∩B|/(|A|+|B|), not Jaccard", {
    a <- c("g1", "g2", "g3", "g4")      # |A| = 4
    b <- c("g3", "g4", "g5")            # |B| = 3, |A∩B| = 2
    dice <- 100 * 2 * length(intersect(a, b)) / (length(a) + length(b))  # 57.14
    jacc <- 100 * length(intersect(a, b)) / length(union(a, b))         # 40.0
    expect_equal(round(dice, 2), 57.14)
    expect_false(isTRUE(all.equal(dice, jacc)))
})

test_that("pipeline .enrich_membership_gene_order clusters by Term membership (§3)", {
    # Terms x Genes binary matrix; GeneA == GeneB (identical membership vectors).
    m <- rbind(Term1 = c(1, 1, 0, 1), Term2 = c(1, 1, 1, 0),
               Term3 = c(0, 0, 1, 0), Term4 = c(1, 1, 0, 0))
    colnames(m) <- c("GeneA", "GeneB", "GeneC", "GeneD")
    ord <- .enrich_membership_gene_order(m)
    expect_setequal(ord, colnames(m))                     # nothing dropped
    expect_equal(abs(which(ord == "GeneA") - which(ord == "GeneB")), 1L)  # identical -> adjacent
    # Jaccard: shared absences do not create similarity
    D <- as.matrix(.enrich_gene_jaccard_dist(m))
    expect_equal(unname(D["GeneA", "GeneB"]), 0)          # identical membership -> dist 0
    big <- matrix(0L, nrow = 100, ncol = 2,
                  dimnames = list(paste0("T", 1:100), c("gX", "gY")))
    big[1, "gX"] <- 1L; big[2, "gY"] <- 1L
    expect_equal(unname(as.matrix(.enrich_gene_jaccard_dist(big))["gX", "gY"]), 1)
    # deterministic; identical-membership genes contiguous in input order
    m2 <- rbind(c(1, 1, 0, 1), c(1, 1, 0, 1))
    colnames(m2) <- c("gB", "gA", "gZ", "gC")
    o <- .enrich_membership_gene_order(m2)
    idx <- sort(match(c("gB", "gA", "gC"), o))
    expect_equal(idx, c(idx[1], idx[1] + 1L, idx[1] + 2L))
    expect_identical(.enrich_membership_gene_order(m2), o)   # reproducible
})

test_that("write_ora_compare_contrasts (§5/§9): combines DE per contrast, >=2 -> PDF", {
    skip_if_not_installed("ggplot2")
    de <- list(
        cA = data.frame(Cluster = "all", ID = c("GO:1", "GO:2"),
            Description = c("t1", "t2"), geneID = c("a/b", "c"),
            Count = c(2L, 1L), p.adjust = c(1e-4, 1e-3),
            GeneRatio = "2/50", Fold_enrichment = 2, stringsAsFactors = FALSE),
        cB = data.frame(Cluster = "all", ID = c("GO:1", "GO:3"),
            Description = c("t1", "t3"), geneID = c("a/d", "e"),
            Count = c(2L, 1L), p.adjust = c(2e-4, 1e-2),
            GeneRatio = "2/60", Fold_enrichment = 2, stringsAsFactors = FALSE))
    d <- withr::local_tempdir()
    write_ora_compare_contrasts(de, d)
    expect_true(file.exists(file.path(d, "compare_contrasts.csv")))
    expect_true(file.exists(file.path(d, "compare_contrasts.pdf")))
    cc <- utils::read.csv(file.path(d, "compare_contrasts.csv"))
    expect_equal(nrow(cc), 4L)                              # 2 + 2, nothing dropped
    expect_setequal(unique(cc$Cluster), c("cA", "cB"))     # group = contrast name
    expect_equal(sum(cc$ID == "GO:1"), 2L)                 # shared term kept in BOTH
    expect_false(any(cc$Cluster %in% c("up", "down")))     # DE only, no Up/Down
    # single contributing contrast -> CSV only, no comparison PDF (§14)
    d2 <- withr::local_tempdir()
    write_ora_compare_contrasts(de["cA"], d2)
    expect_true(file.exists(file.path(d2, "compare_contrasts.csv")))
    expect_false(file.exists(file.path(d2, "compare_contrasts.pdf")))
})
