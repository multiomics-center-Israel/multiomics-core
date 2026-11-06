# tests/testthat/helper-mummichog-ma.R
#
# Shared helpers for the MetaboAnalyst-compatible mummichog tests (06j-06l).
# Everything is SYNTHETIC.

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

