# tests/testthat/test-metabolite-network.R
#
# Unit tests for the pinned-reference metabolite network (domain 08 / 08b):
# reference read+validation, pure reaction-pair filtering with edge provenance,
# graph construction, and the guarantee that no live-KEGG path remains.

# Write a small reference .tsv.gz with the canonical "# key: value" header.
write_ref_gz <- function(path, rows,
                         schema_version = 1L,
                         method = "equation_side_cartesian_product",
                         extra_header = character(0),
                         columns = c("reaction_id", "substrate_id", "product_id",
                                     "equation", "equation_arrow")) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  writeLines(c(
    paste0("# schema_version: ", schema_version),
    paste0("# pair_definition_method: ", method),
    "# source: KEGG REST",
    extra_header,
    paste(columns, collapse = "\t")
  ), con)
  if (nrow(rows) > 0) {
    for (i in seq_len(nrow(rows))) {
      writeLines(paste(as.character(unlist(rows[i, ])), collapse = "\t"), con)
    }
  }
}

valid_rows <- function() {
  data.frame(
    reaction_id    = c("R00001", "R00002"),
    substrate_id   = c("C00001", "C00002"),
    product_id     = c("C00002", "C00003"),
    equation       = c("C00001 <=> C00002", "C00002 <=> C00003"),
    equation_arrow = c("<=>", "<=>"),
    stringsAsFactors = FALSE)
}

# ---- read + validate ------------------------------------------------------

test_that("read_kegg_reaction_pairs reads a valid reference + header", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, valid_rows())
  df <- read_kegg_reaction_pairs(f)

  expect_equal(nrow(df), 2L)
  expect_true(all(c("reaction_id", "substrate_id", "product_id", "equation",
                    "equation_arrow") %in% colnames(df)))
  hdr <- attr(df, "reference_header")
  expect_identical(hdr$pair_definition_method, "equation_side_cartesian_product")
})

test_that("validation rejects an unsupported schema_version", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, valid_rows(), schema_version = 2L)
  expect_error(read_kegg_reaction_pairs(f, expected_schema_version = 1L),
               "schema_version")
})

test_that("validation rejects a wrong pair_definition_method", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, valid_rows(), method = "kegg_rclass")
  expect_error(read_kegg_reaction_pairs(f), "pair_definition_method")
})

test_that("validation rejects a self-pair (substrate == product)", {
  rows <- valid_rows()
  rows$product_id[1] <- rows$substrate_id[1]   # C00001 <=> C00001
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, rows)
  expect_error(read_kegg_reaction_pairs(f), "self-pair")
})

test_that("validation rejects duplicate rows under the key", {
  rows <- valid_rows()[c(1, 1), ]                # duplicate R00001/C00001/C00002
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, rows)
  expect_error(read_kegg_reaction_pairs(f), "duplicate")
})

test_that("validation rejects a missing required column", {
  f <- withr::local_tempfile(fileext = ".tsv.gz")
  write_ref_gz(f, valid_rows()[, 1:3],
               columns = c("reaction_id", "substrate_id", "product_id"))
  expect_error(read_kegg_reaction_pairs(f), "missing required column")
})

test_that("read_kegg_reaction_pairs errors on a missing file (no fallback)", {
  expect_error(read_kegg_reaction_pairs("/no/such/ref.tsv.gz"), "not found")
})

# ---- pure filtering + edge provenance -------------------------------------

test_that("filter_reaction_pairs_to_features keeps in-set pairs with provenance", {
  rp <- data.frame(
    reaction_id  = c("R00001", "R00002", "R00003", "R00004"),
    substrate_id = c("C00001", "C00002", "C00001", "C00001"),
    product_id   = c("C00002", "C00003", "C00002", "C00009"),  # R00004 out of set
    stringsAsFactors = FALSE)
  edges <- filter_reaction_pairs_to_features(rp, c("C00001", "C00002", "C00003"))

  expect_equal(nrow(edges), 2L)
  e1 <- edges[edges$from == "C00001" & edges$to == "C00002", ]
  expect_identical(e1$reaction_ids, "R00001;R00003")  # sorted, unique, joined
  expect_equal(e1$n_reactions, 2L)
  e2 <- edges[edges$from == "C00002" & edges$to == "C00003", ]
  expect_equal(e2$n_reactions, 1L)
})

test_that("filter collapses reverse-direction pairs into one undirected edge", {
  rp <- data.frame(
    reaction_id  = c("R00001", "R00005"),
    substrate_id = c("C00001", "C00002"),   # R00005 is the reverse direction
    product_id   = c("C00002", "C00001"),
    stringsAsFactors = FALSE)
  edges <- filter_reaction_pairs_to_features(rp, c("C00001", "C00002"))
  expect_equal(nrow(edges), 1L)
  expect_identical(edges$reaction_ids, "R00001;R00005")
  expect_equal(edges$n_reactions, 2L)
})

test_that("filter excludes self-pairs and returns an empty frame when nothing matches", {
  rp <- data.frame(reaction_id = "R00001", substrate_id = "C00001",
                   product_id = "C00001", stringsAsFactors = FALSE)
  expect_equal(nrow(filter_reaction_pairs_to_features(rp, "C00001")), 0L)
  expect_equal(nrow(filter_reaction_pairs_to_features(rp, "C99999")), 0L)
})

# ---- graph construction from a reference (no KEGG) ------------------------

test_that("build_de_metabolite_network builds a graph carrying edge provenance", {
  skip_if_not_installed("igraph")
  de <- data.frame(feature_id = c("F1", "F2", "F3"),
                   logFC = c(1, -1, 2), P.Value = c(0.001, 0.02, 0.4),
                   stringsAsFactors = FALSE)
  ann <- data.frame(feature_id = c("F1", "F2", "F3"),
                    KEGG = c("C00001", "C00002", "C00003"),
                    stringsAsFactors = FALSE)
  rp <- valid_rows()

  net <- build_de_metabolite_network(de, ann, reaction_pairs = rp,
                                     p_cutoff = 0.05)
  # F1 (p .001) and F2 (p .02) are significant; F3 (p .4) is not -> the only
  # in-set edge is C00001-C00002 (R00001).
  expect_equal(igraph::ecount(net$graph), 1L)
  expect_identical(igraph::edge_attr(net$graph, "reaction_ids"), "R00001")
  expect_equal(as.integer(igraph::edge_attr(net$graph, "n_reactions")), 1L)
  # provenance is exported too (network_edges.tsv is written from net$edges)
  expect_true(all(c("reaction_ids", "n_reactions") %in% colnames(net$edges)))
})

# ---- the live-KEGG path is gone -------------------------------------------

test_that("no live-KEGG fetch function remains in the network path", {
  expect_false(exists("fetch_kegg_reaction_pairs", mode = "function"))
})
