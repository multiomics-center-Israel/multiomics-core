# tests/testthat/test-mummichog-ma-library.R
#
# MetaboAnalyst-compatible mummichog, component 1: the compound library built
# from our model JSON (06j). Unit tests on synthetic inputs, plus a parity test
# against the pinned MetaboAnalystR library builder
# (fixtures/mummichog_ma_parity/generate_library_reference.R). All synthetic.

ma_fixture_dir <- file.path(
  normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE),
  "tests", "testthat", "fixtures", "mummichog_ma_parity")

# Write a minimal Azimuth-shape model. `pathways` is a named list
# (pathway id -> compound ids); each pathway gets one reaction holding all of
# its compounds as reactants, so the reaction-derived membership is exactly the
# list given.
ma_write_model <- function(cpd, pathways, path = tempfile(fileext = ".json"),
                           pathway_names = names(pathways)) {
  rxns <- lapply(seq_along(pathways), function(k) {
    list(id = sprintf("R9%04d", k), reactants = I(pathways[[k]]),
         products = I(character()))
  })
  pws <- lapply(seq_along(pathways), function(k) {
    list(id = names(pathways)[[k]], name = pathway_names[[k]],
         list_of_reactions = I(sprintf("R9%04d", k)))
  })
  cpds <- lapply(seq_len(nrow(cpd)), function(i) {
    list(id = cpd$id[[i]], name = cpd$name[[i]],
         neutral_mono_mass = cpd$mw[[i]])
  })
  jsonlite::write_json(
    list(id = "syn_model", version = "20000101", meta_data = list(),
         list_of_compounds = cpds, list_of_reactions = rxns,
         list_of_pathways = pws),
    path, auto_unbox = TRUE, digits = NA)
  path
}

ma_cpd <- data.frame(id = c("C90001", "C90002", "C90003"),
                     name = c("syn A", "syn B", "syn C"),
                     mw = c(180.063388, 146.105528, 300.9995),
                     stringsAsFactors = FALSE)

# ---- adduct matrix ------------------------------------------------------------

test_that("adduct matrices carry upstream's 25 positive and 16 negative forms", {
  pos <- mmc_ma_adduct_matrix(c(180.063388, 146.105528), "positive")
  neg <- mmc_ma_adduct_matrix(c(180.063388, 146.105528), "negative")
  expect_identical(dim(pos), c(2L, 25L))
  expect_identical(dim(neg), c(2L, 16L))
  expect_identical(colnames(pos)[c(1, 2, 10, 14, 25)],
                   c("M[1+]", "M+H[1+]", "M+Na[1+]", "M-H2O+H[1+]",
                     "M-HCOOK+H[1+]"))
  expect_identical(colnames(neg)[c(1, 2, 8, 16)],
                   c("M-H[-]", "M-2H[2-]", "M-H2O-H[-]", "M-H+O[-]"))
  expect_null(rownames(pos))
})

test_that("adduct m/z follow upstream's formulas", {
  mw <- 180.063388
  p <- 1.00727646677
  pos <- mmc_ma_adduct_matrix(mw, "positive")
  neg <- mmc_ma_adduct_matrix(mw, "negative")
  # A single cell extracted from the matrix keeps its column name; the
  # expectations compare values only.
  expect_identical(unname(pos[1, "M[1+]"]), mw)
  expect_identical(unname(pos[1, "M+H[1+]"]), mw + p)
  expect_identical(unname(pos[1, "M+2H[2+]"]), mw/2 + p)
  expect_identical(unname(pos[1, "M+Na[1+]"]), mw + 21.9820 + p)
  expect_identical(unname(neg[1, "M-H[-]"]), mw - p)
  expect_identical(unname(neg[1, "M+Na-2H[-]"]), mw + 21.9820 - 2*p)
  expect_identical(unname(neg[1, "M+Cl[-]"]), mw + 34.9689)
  # textbook check, glucose: [M+H]+ 181.0707, [M-H]- 179.0561
  expect_equal(unname(pos[1, "M+H[1+]"]), 181.0707, tolerance = 1e-4)
  expect_equal(unname(neg[1, "M-H[-]"]), 179.0561, tolerance = 1e-4)
})

test_that("adduct matrix refuses missing or non-numeric masses", {
  expect_error(mmc_ma_adduct_matrix(c(100, NA), "positive"), "no missing")
  expect_error(mmc_ma_adduct_matrix("100", "negative"), "numeric")
  expect_error(mmc_ma_adduct_matrix(100, "mixed"))
})

# ---- m/z tree -----------------------------------------------------------------

test_that("the m/z tree bins row indices by floored adduct m/z, 50 < n < 2000", {
  # row 1: M[1+] = 100.2 -> bin 100; row 2: M[1+] = 49.95 -> excluded
  mz <- cbind(c(100.2, 49.95), c(101.9, 50.5), c(100.7, 2000.1))
  tr <- mmc_ma_cpd_tree(mz)
  expect_length(tr, 2001L)
  expect_identical(tr[[49]], "")
  expect_identical(tr[[2001]], "")
  expect_identical(tr[[100]], 1L)          # two adducts of row 1, deduplicated
  expect_identical(tr[[101]], 1L)
  expect_null(tr[[50]])                    # floor 50 is excluded (n > 50)
  expect_null(tr[[2000]])                  # floor 2000 is excluded (n < 2000)
  expect_false(2L %in% unlist(tr))
})

test_that("the m/z tree holds several rows in one bin in row order", {
  mz <- cbind(c(300.1, 300.9, 300.5))
  expect_identical(mmc_ma_cpd_tree(mz)[[300]], 1:3)
})

# ---- converter ------------------------------------------------------------------

test_that("a model JSON converts into the MetaboAnalystR library shape", {
  path <- ma_write_model(ma_cpd, list(syn01 = c("C90001", "C90002"),
                                      syn02 = c("C90002", "C90003")),
                         pathway_names = c("Pathway one", "Pathway two"))
  lib <- mmc_ma_library_from_model(path)

  expect_named(lib, c("pathways", "cpd.tree", "cpd.lib", "model"))
  expect_identical(lib$cpd.lib$id, ma_cpd$id)
  expect_identical(lib$cpd.lib$mw, ma_cpd$mw)
  expect_identical(lib$cpd.lib$adducts$positive,
                   mmc_ma_adduct_matrix(ma_cpd$mw, "positive"))
  expect_identical(lib$cpd.tree$negative,
                   mmc_ma_cpd_tree(mmc_ma_adduct_matrix(ma_cpd$mw, "negative")))
  expect_identical(lib$pathways$name, c("Pathway one", "Pathway two"))
  expect_identical(lib$pathways$id, c("syn01", "syn02"))
  expect_null(names(lib$pathways$cpds))
  expect_setequal(lib$pathways$cpds[[1]], c("C90001", "C90002"))
  expect_setequal(lib$pathways$cpds[[2]], c("C90002", "C90003"))
  expect_identical(lib$model$id, "syn_model")
})

test_that("pathway compounds absent from the compound table are dropped", {
  path <- ma_write_model(ma_cpd, list(syn01 = c("C90001", "C99999")))
  lib <- mmc_ma_library_from_model(path)
  expect_identical(lib$pathways$cpds[[1]], "C90001")
})

test_that("the converter refuses models it cannot build a library from", {
  expect_error(mmc_ma_library_from_model(tempfile()), "not found")

  no_mass <- ma_cpd
  no_mass$mw[2] <- NA
  expect_error(mmc_ma_library_from_model(ma_write_model(no_mass, list(a = "C90001"))),
               "neutral_mono_mass")

  dup <- ma_cpd
  dup$id[3] <- "C90001"
  expect_error(mmc_ma_library_from_model(ma_write_model(dup, list(a = "C90001"))),
               "unique")

  mcg2 <- tempfile(fileext = ".json")
  jsonlite::write_json(list(metabolic_pathways = list()), mcg2, auto_unbox = TRUE)
  expect_error(mmc_ma_library_from_model(mcg2), "list_of_compounds")
})

# ---- parity with the pinned MetaboAnalystR builder ------------------------------

test_that("the library matches the pinned MetaboAnalystR builder", {
  ref_file <- file.path(ma_fixture_dir, "library_reference.rds")
  skip_if_not(file.exists(ref_file), "library parity reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")

  inp <- ref$inputs
  cpd <- data.frame(id = inp$cpd.lib$id, name = inp$cpd.lib$name,
                    mw = inp$cpd.lib$mw, stringsAsFactors = FALSE)
  pws <- stats::setNames(inp$pathways$cpds, inp$pathways$id)
  lib <- mmc_ma_library_from_model(
    ma_write_model(cpd, pws, pathway_names = inp$pathways$name))

  # the JSON round trip must not have touched a mass
  expect_identical(lib$cpd.lib$mw, inp$cpd.lib$mw)
  for (m in c("positive", "negative")) {
    expect_identical(lib$cpd.lib$adducts[[m]], ref$adducts[[m]])
    expect_identical(lib$cpd.tree[[m]], ref$cpd.tree[[m]])
  }
})
