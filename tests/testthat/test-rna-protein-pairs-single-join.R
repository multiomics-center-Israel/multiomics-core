make_rna_de <- function() {
  data.frame(
    feature_id = c("g1", "g2", "g3", "g4"),
    log2FC     = c(1.0, -2.0, 0.5, 3.0),
    padj       = c(0.01, 0.20, 0.60, 0.001),
    stringsAsFactors = FALSE
  )
}

make_prot_de <- function() {
  data.frame(
    feature_id = c("p1", "p2a", "p2b", "pC;pD"),
    log2FC     = c(0.8, -1.5, -1.1, 2.2),
    padj       = c(0.02, 0.30, 0.04, 0.005),
    stringsAsFactors = FALSE
  )
}

make_mapping <- function() {
  data.frame(
    gene_id        = c("g1", "g2", "g2", "g4", "g3"),
    protein_id     = c("p1", "p2a", "p2b", "pC;pD", "absent"),
    mapping_source = "custom_file",
    stringsAsFactors = FALSE
  )
}

test_that("build_rna_protein_pairs joins through the mapping and keeps both ids", {
  out <- build_rna_protein_pairs(make_rna_de(), make_prot_de(), make_mapping())

  expect_true(all(c("gene_id", "protein_id", "logFC_rna", "logFC_prot",
                    "padj_rna", "padj_prot") %in% names(out)))
  # g3 maps to a protein that is absent from the DE table, so it drops out.
  expect_equal(nrow(out), 4)
  expect_false("g3" %in% out$gene_id)
})

test_that("a gene pairing with several protein groups keeps every pair", {
  # g2 -> p2a and p2b. Collapsing this would silently pick one protein per gene.
  out <- build_rna_protein_pairs(make_rna_de(), make_prot_de(), make_mapping())
  expect_equal(sum(out$gene_id == "g2"), 2)
  expect_setequal(out$protein_id[out$gene_id == "g2"], c("p2a", "p2b"))
})

test_that("a semicolon-separated protein group still joins", {
  # This is the case the second code path used to drop: 60 real pairs in the
  # wasp run, which moved the reported correlation from 0.0615 to 0.0636.
  out <- build_rna_protein_pairs(make_rna_de(), make_prot_de(), make_mapping())
  expect_true("pC;pD" %in% out$protein_id)
  expect_equal(out$logFC_prot[out$protein_id == "pC;pD"], 2.2)
})

test_that("values are carried through unchanged", {
  out <- build_rna_protein_pairs(make_rna_de(), make_prot_de(), make_mapping())
  row <- out[out$gene_id == "g1", ]
  expect_equal(row$logFC_rna, 1.0)
  expect_equal(row$logFC_prot, 0.8)
  expect_equal(row$padj_rna, 0.01)
  expect_equal(row$padj_prot, 0.02)
})

test_that("build_rna_protein_pairs returns zero rows rather than erroring", {
  expect_equal(nrow(build_rna_protein_pairs(make_rna_de(), make_prot_de(), NULL)), 0)
  empty <- data.frame(gene_id = character(0), protein_id = character(0),
                      stringsAsFactors = FALSE)
  expect_equal(nrow(build_rna_protein_pairs(make_rna_de(), make_prot_de(), empty)), 0)

  no_overlap <- data.frame(gene_id = "gX", protein_id = "pX",
                           stringsAsFactors = FALSE)
  expect_equal(nrow(build_rna_protein_pairs(make_rna_de(), make_prot_de(), no_overlap)), 0)
})

test_that("compute_de_concordance reports the same pairs as the shared join", {
  # The regression this guards: two paths joining independently and disagreeing.
  pairs <- build_rna_protein_pairs(make_rna_de(), make_prot_de(), make_mapping())
  conc <- compute_de_concordance(make_rna_de(), make_prot_de(), make_mapping(),
                                 out_dir = NULL, contrast_label = "c1")
  skip_if(is.null(conc), "compute_de_concordance applies a minimum-row guard")

  expect_setequal(paste(conc$gene_id, conc$protein_id),
                  paste(pairs$gene_id, pairs$protein_id))
  expect_equal(sort(conc$rna_log2FC), sort(pairs$logFC_rna))
  expect_equal(sort(conc$protein_log2FC), sort(pairs$logFC_prot))
})
