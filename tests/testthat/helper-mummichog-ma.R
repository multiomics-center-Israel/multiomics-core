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

# One contrast's MetaboAnalyst-compatible result files, written by hand so the
# Pathway Evidence can be checked against values worked out in advance.
# SYNTHETIC throughout; KEGG-style ids only exercise the id comparison.
#
# Model:  Glycolysis          (syn00010): C00031, C00267, C00022
#         Pyruvate metabolism (syn00620): C00022, C00186
#         TCA cycle           (syn00020): C00158, C00022, C00042
#         C99999 is in no pathway.
# ECs (as the engine grouped them):
#   EC0001  F01 (+, C00031)  F02 (-, C00031)  F04 (+, C99999)   both polarities;
#           F04 matched only a compound outside every pathway
#   EC0002  "0012" (+, C00022 and C00186)  F06 (+, C00022)
#   EC0003  F01 (+, C00267)                                      F01 in two ECs
# F07 (-) has F06's m/z and a larger p: the engine dropped it, so it is in
# ma_input.tsv but in no match.
# Pathway table: Glycolysis ec_hits EC0001;EC0002;EC0003 (Sig.EC 3, FET 0.01);
#                Pyruvate metabolism ec_hits EC0002 (Sig.EC 1, FET 0.2).
# TCA cycle is never in the ORA table (no significant EC). With
# `gsea_only = TRUE` the GSEA table also tests it (p 0.005, FDR 0.03), giving a
# pathway with a GSEA signal and no significant EC overlap.
# `contrast`, `fet`, `holm`, `gsea_p`, `gsea_fdr` vary the contrast so that
# several fixtures in one `dir` form a multi-contrast run on one model.
#
# Returns the files, the pathway table as the report reads it, the model index
# (with the sha256 attribute mmc_load_model_pathways() attaches), row_data, the
# DE table, the config and the paths.
ma_evidence_fixture <- function(dir = tempfile("ma_ev_"), sig_ec = c(3, 1),
                                with_gsea = TRUE, record_sha = TRUE,
                                contrast = "HL_vs_LL",
                                fet = c(0.01, 0.2), holm = c(0.02, 0.2),
                                gsea_p = c(0.03, 0.6), gsea_fdr = c(0.06, 0.6),
                                gsea_only = FALSE, run_p_cutoff = 0.05) {
  ma_dir <- file.path(dir, "mummichog_pinned", contrast, "metaboanalyst")
  dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)
  cpd <- data.frame(
    id   = c("C00031", "C00267", "C00022", "C00186", "C99999", "C00158",
             "C00042"),
    name = c("D-Glucose", "alpha-D-Glucose", "Pyruvate", "(S)-Lactate",
             "Foo compound", "Citrate", "Succinate"),
    mw   = c(180.0634, 180.0634, 88.0160, 90.0317, 146.0691, 192.0270,
             118.0266),
    stringsAsFactors = FALSE)
  model_path <- ma_write_model(
    cpd, list(syn00010 = c("C00031", "C00267", "C00022"),
              syn00620 = c("C00022", "C00186"),
              syn00020 = c("C00158", "C00022", "C00042")),
    path = file.path(dir, "model.json"),
    pathway_names = c("Glycolysis", "Pyruvate metabolism", "TCA cycle"))
  sha <- .mmc_sha256_file(model_path)

  res <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    ma_rank = 1:2, pathway_id = c("syn00010", "syn00620"),
    pathway = c("Glycolysis", "Pyruvate metabolism"),
    "Pathway total" = c(3, 2), Hits.total = c(3, 2), Hits.sig = c(3, 2),
    Expected = c(0.5, 0.3), FET = fet, EASE = c(0.02, 0.3),
    Total.EC = c(4, 1), Sig.EC = sig_ec, Gamma = c(0.03, 0.25),
    Emp.Hits = c(1, 20), Empirical = c(0.01, 0.2),
    Holm_adj_FET = holm, Holm_adj_EASE = c(0.04, 0.3),
    Holm_adj_Gamma = c(0.06, 0.25),
    ec_hits = c("EC0001;EC0002;EC0003", "EC0002"))
  readr::write_tsv(res, file.path(ma_dir, "ma_pathway_results.tsv"), na = "NA")

  matched <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    feature_id = c("F01", "F02", "F04", "0012", "0012", "F06", "F01"),
    polarity   = c("positive", "negative", "positive", "positive", "positive",
                   "positive", "positive"),
    Query.Mass = c("181.0707", "179.0561", "147.0764", "89.0233", "89.0233",
                   "120.5", "181.0707"),
    Matched.Compound = c("C00031", "C00031", "C99999", "C00022", "C00186",
                         "C00022", "C00267"),
    Matched.Form = c("M+H[1+]", "M-H[-]", "M+Na[1+]", "M+H[1+]", "M+H[1+]",
                     "M+K[1+]", "M+H[1+]"),
    Retention.Time = c("2.1", "2.15", "2.2", "1", "1", "1.05", "2.1"),
    Mass.Diff = c("1e-04", "2e-04", "1e-04", "3e-04", "5e-04", "2e-04", "1e-04"),
    Empirical.Compound = c("EC0001", "EC0001", "EC0001", "EC0002", "EC0002",
                           "EC0002", "EC0003"))
  readr::write_tsv(matched, file.path(ma_dir, "ma_matched_features.tsv"), na = "NA")

  input <- data.frame(
    feature_id = c("F01", "F02", "F04", "0012", "F06", "F07"),
    mz = c(181.0707, 179.0561, 147.0764, 89.0233, 120.5, 120.5),
    rt = c(2.10, 2.15, 2.20, 1.00, 1.05, 1.06),
    p_value = c(0.001, 0.020, 0.300, 0.004, 0.005, 0.200),
    statistic = c(4.0, 2.5, -0.4, -3.0, 3.2, -1.1),
    polarity = c("positive", "negative", "positive", "positive", "positive",
                 "negative"),
    stringsAsFactors = FALSE)
  readr::write_tsv(input, file.path(ma_dir, "ma_input.tsv"), na = "NA")

  info <- list(engine = "metaboanalyst", contrast = contrast,
               params = list(ionization_mode = "mixed",
                             polarity_column = "Reference Ion",
                             tolerance_ppm = 10, p_cutoff = run_p_cutoff,
                             n_permutations = 100, force_primary_ion = TRUE,
                             rt_frac = 0.02, min_pathway_size = 3),
               model = list(id = "syn_model", version = "20000101",
                            sha256 = if (record_sha) sha else NULL),
               reference = list(repo = "xia-lab/MetaboAnalystR",
                                commit = "398476ae2a0c996e390925fa62adc46b5ecde334"),
               fitdistrplus_version = "0.0.0",
               gamma = list(fitted = FALSE, error = "synthetic fixture"),
               # As the engine records them: 6 DE features, all passed; F07
               # removed as a duplicate m/z; 4 of the 5 retained have p <= 0.05.
               counts = list(features_given = 6, features_kept = 5,
                             removed_out_of_range = 0, removed_duplicate_mz = 1,
                             significant_features = 4, empirical_compounds = 3,
                             significant_ecs = 3, pathways_tested = 3,
                             pathways_reported = 2))
  if (with_gsea) {
    gsea <- data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      Pathway = c("Glycolysis", "Pyruvate metabolism"),
      ES = c(0.7, -0.3), NES = c(1.6, -0.8), P.Value = gsea_p,
      padj = gsea_fdr,
      "Leading-edge EmpiricalCompounds" = c("EC0001; EC0002", "EC0002"),
      EmpiricalCompounds = c("EC0001; EC0002; EC0003", "EC0002"),
      "ES defaulted by reference implementation" = c(FALSE, FALSE),
      "ES fallback reason" = c(NA_character_, NA_character_))
    if (gsea_only) {
      # In the model and tested by GSEA, but with no significant EC, so the
      # engine leaves it out of the ORA table.
      gsea <- rbind(gsea, data.frame(
        check.names = FALSE, stringsAsFactors = FALSE,
        Pathway = "TCA cycle", ES = 0.8, NES = 1.9, P.Value = 0.005, padj = 0.03,
        "Leading-edge EmpiricalCompounds" = "EC0002",
        EmpiricalCompounds = "EC0002",
        "ES defaulted by reference implementation" = FALSE,
        "ES fallback reason" = NA_character_))
    }
    readr::write_tsv(gsea, file.path(ma_dir, "ma_gsea_results.tsv"), na = "NA")
    # As the engine writes it for the default method. A record without
    # `method` would be read as the legacy method (see 06n).
    info$gsea <- list(ran = TRUE, method = "fgsea",
                      method_label = mmc_gsea_method_label("fgsea", "0.0.0"),
                      metric_label = "signed t statistic (synthetic)",
                      log2err = list(max = 0.1, n_unavailable = 0),
                      params = list(nperm_simple = 10000, min_size = 3,
                                    max_size = 5000, seed = 123,
                                    fgsea_version = "0.0.0"),
                      deviations = list(),
                      reference = list(method = "fgsea::fgsea (fgseaMultilevel)",
                                       package = "fgsea", version = "0.0.0"),
                      adjustment = "BH (p.adjust 'fdr')", n_ec = 3, n_positions = 3,
                      n_pathways = 2)
  } else {
    info$gsea <- list(ran = FALSE, error = NULL)
  }
  jsonlite::write_json(info, file.path(ma_dir, "ma_run_info.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")

  files <- sort(list.files(ma_dir, full.names = TRUE))
  model <- read_mummichog_model_pathways(model_path)
  attr(model, "sha256") <- sha

  row_data <- data.frame(
    feature_id = c("F01", "F02", "F04", "0012", "F06", "F07"),
    Name = c("D-Glucose", "Citrate", "Foo compound", NA, "Pyruvate", "(S)-Lactate"),
    KEGG = c("C00031", "C00158", "C99999", NA, "C00022", "C00186"),
    stringsAsFactors = FALSE)
  de <- data.frame(feature_id = input$feature_id,
                   logFC = c(1.5, 0.9, -0.1, -1.2, 1.1, -0.4),
                   P.Value = input$p_value, stringsAsFactors = FALSE)
  cfg <- list(modes = list(metabolomics = list(
    organism = "Synthetic organism",
    enrichment = list(mummichog = list(
      enabled = TRUE, engine = "metaboanalyst", ionization_mode = "mixed",
      polarity_column = "Reference Ion", p_cutoff = 0.05,
      model_json = model_path)))),
    params = list(seed = 1))
  list(dir = dir, ma_dir = ma_dir, files = files, model = model,
       model_path = model_path, sha = sha,
       pathways = read_mummichog_ma_pathways(files),
       row_data = row_data, annot = normalize_metab_annotation(row_data),
       de = de, de_res = list(de_tables = stats::setNames(list(de), contrast)),
       contrast = contrast, config = cfg)
}

# Two contrasts on one model, built into real report sections. Worked out by
# hand (SYNTHETIC):
#   HL_vs_LL  Glycolysis  ORA FET 0.01, Holm 0.02 (adjusted); GSEA p 0.03, FDR 0.06 (nominal)
#             Pyruvate    ORA FET 0.2 (none);                  GSEA p 0.6 (none)
#             TCA cycle   not in the ORA table (no significant EC); GSEA p 0.005, FDR 0.03
#   LL_vs_DD  Glycolysis  ORA FET 0.3 (none);  GSEA p 0.02, FDR 0.04 -> GSEA-only, adjusted
#             Pyruvate    ORA FET 0.004, Holm 0.01 (adjusted); GSEA p 0.7 (none)
#             TCA cycle   no result
ma_two_contrast_sections <- function(dir = tempfile("ma_two_")) {
  a <- ma_evidence_fixture(dir, contrast = "HL_vs_LL", gsea_only = TRUE)
  b <- ma_evidence_fixture(dir, contrast = "LL_vs_DD",
                           fet = c(0.3, 0.004), holm = c(0.3, 0.01),
                           gsea_p = c(0.02, 0.7), gsea_fdr = c(0.04, 0.7))
  secs <- suppressWarnings(suppressMessages(build_mummichog_report_sections(
    list(HL_vs_LL = a$pathways, LL_vs_DD = b$pathways), a$config,
    files = c(a$files, b$files), row_data = a$row_data,
    de_res = list(de_tables = list(HL_vs_LL = a$de, LL_vs_DD = b$de)),
    model = b$model)))
  list(sections = secs, a = a, b = b, dir = dir)
}
