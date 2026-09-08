# tests/testthat/test-mummichog-evidence.R
#
# Unit tests for the mummichog pathway supporting-evidence layer (06f):
# the metabolic-model pathway loader, the EmpiricalCompound readers, the
# schema-adaptive annotation normaliser, the annotation-agreement rules and the
# end-to-end Pathway -> EC -> candidate -> feature -> annotation -> agreement
# trace.
#
# All fixtures are SYNTHETIC. The compound ids are real KEGG accessions only so
# that the ID-comparison path is exercised realistically; the pathway membership
# below is invented for the test and is not a biological claim.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# A minimal Azimuth-format metabolic model, the shape model_ref / model_json
# publish. Pathways are derived from their reactions' reactants + products,
# exactly as mummichog's json_convert_azmuth_mummichog() does.
#
#   Purine metabolism      : C00020 (AMP), C00362 (dGMP), C01367 (3'-AMP)
#   Fatty acid degradation : C00020 (AMP), C00186 (L-Lactate)
#   Glycolysis             : C00031 (D-Glucose), C00095 (D-Fructose)
write_model_json <- function(dir) {
  model <- list(
    id = "test_model",
    meta_data = list(version = "test-1"),
    list_of_compounds = list(
      list(id = "C00020", name = "AMP;Adenosine 5'-monophosphate",
           identifiers = list(kegg.compound = "C00020")),
      list(id = "C00362", name = "dGMP",
           identifiers = list(kegg.compound = "C00362")),
      list(id = "C01367", name = "3'-AMP",
           identifiers = list(kegg.compound = "C01367")),
      list(id = "C00031", name = "D-Glucose;Glucose",
           identifiers = list(kegg.compound = "C00031")),
      list(id = "C00186", name = "(S)-Lactate;L-Lactate",
           identifiers = list(kegg.compound = "C00186")),
      list(id = "C00095", name = "D-Fructose",
           identifiers = list(kegg.compound = "C00095"))
    ),
    list_of_reactions = list(
      list(id = "R1", reactants = list("C00020"), products = list("C00362")),
      list(id = "R2", reactants = list("C00020"), products = list("C01367")),
      list(id = "R3", reactants = list("C00020"), products = list("C00186")),
      list(id = "R4", reactants = list("C00031"), products = list("C00095"))
    ),
    list_of_pathways = list(
      list(id = "p1", name = "Purine metabolism",
           list_of_reactions = list("R1", "R2")),
      list(id = "p2", name = "Fatty acid degradation",
           list_of_reactions = list("R3")),
      list(id = "p3", name = "Glycolysis",
           list_of_reactions = list("R4"))
    )
  )
  path <- file.path(dir, "model.json")
  jsonlite::write_json(model, path, auto_unbox = TRUE)
  path
}

# A faithful (if tiny) mummichog v2 result tree for ONE contrast, laid out the
# way the pinned stage writes it (mummichog_pinned/<contrast>/v2/<stamp>/tables).
#
#   E196 -> candidates AMP / dGMP / 3'-AMP, built from TWO features (feat_1, feat_2)
#   E2   -> candidate D-Glucose,            one feature (feat_3)
#   E3   -> candidate L-Lactate,            one feature (feat_4)
#   E4   -> candidate D-Fructose,           one feature (feat_5)
build_evidence_tree <- function(root, contrast_dir = "HL_vs_LL") {
  tables <- file.path(root, "mummichog_pinned", contrast_dir, "v2",
                      "1700000000.1.run", "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "pathway\toverlap_size\tpathway_size\tp-value\toverlap_EmpiricalCompounds (id)\toverlap_features (id)\toverlap_features (name)",
    "Fatty acid degradation\t2\t2\t0.01\tE196,E3\tC00020\tAMP",
    "Purine metabolism\t1\t3\t0.04\tE196\tC00020\tAMP",
    "Glycolysis\t2\t2\t0.30\tE2,E4\t\t"
  ), file.path(tables, "mcg_pathwayanalysis_HL_vs_LL.tsv"))

  writeLines(c(
    "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names",
    "E196\trow1;row2\trow1_M+H[1+];row2_M+Na[1+]\tC00020;C00362;C01367\tAMP;Adenosine 5'-monophosphate$dGMP$3'-AMP",
    "E2\trow3\trow3_M+H[1+]\tC00031\tD-Glucose;Glucose",
    "E3\trow4\trow4_M-H[-]\tC00186\t(S)-Lactate;L-Lactate",
    "E4\trow5\trow5_M+H[1+]\tC00095\tD-Fructose"
  ), file.path(tables, "ListOfEmpiricalCompounds.tsv"))

  writeLines(c(
    "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    "row1\tE196\trow1_M+H[1+];row2_M+Na[1+]\tC00020;C00362;C01367\tAMP$dGMP$3'-AMP\trow1\t348.07\t3.1\t0.002\t4.1\tfeat_1",
    "row2\tE196\trow1_M+H[1+];row2_M+Na[1+]\tC00020;C00362;C01367\tAMP$dGMP$3'-AMP\trow2\t370.05\t3.1\t0.010\t2.6\tfeat_2",
    "row3\tE2\trow3_M+H[1+]\tC00031\tD-Glucose;Glucose\trow3\t181.07\t1.2\t0.030\t-2.0\tfeat_3",
    "row4\tE3\trow4_M-H[-]\tC00186\t(S)-Lactate;L-Lactate\trow4\t89.02\t0.8\t0.001\t-5.5\tfeat_4",
    "row5\tE4\trow5_M+H[1+]\tC00095\tD-Fructose\trow5\t181.07\t2.4\t0.400\t0.7\tfeat_5"
  ), file.path(tables, "userInput_to_EmpiricalCompounds.tsv"))

  readr::write_tsv(
    data.frame(dir = contrast_dir, contrast = "HL_vs_LL",
               stringsAsFactors = FALSE),
    file.path(root, "mummichog_pinned", "contrasts.tsv"))

  list_mummichog_files(root)
}

# Dataset A: KEGG ids AND MSI identification levels.
#   feat_1  AMP, KEGG C00020, Level 1  -> matches the AMP candidate by ID
#   feat_2  unannotated                -> Not assessed on its own
#   feat_3  D-Glucose, no id, Level 2  -> matches by normalised name
#   feat_4  D-Glucose + KEGG C00031    -> conflicts with the L-Lactate candidate
#   feat_5  unannotated                -> Not assessed
row_data_with_levels <- function() {
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    feature_id = paste0("feat_", 1:5),
    Name       = c("AMP", NA, "D-glucose", "D-Glucose", NA),
    KEGG       = c("C00020", NA, NA, "C00031", NA),
    identification_level = c(1L, NA, 2L, 1L, NA)
  )
}

# Dataset B: names only — no ids, no level system at all.
row_data_no_levels <- function() {
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    feature_id = paste0("feat_", 1:5),
    Name       = c("AMP", NA, "D-glucose", "D-Glucose", NA)
  )
}


# ---------------------------------------------------------------------------
# metabolic model loader
# ---------------------------------------------------------------------------

test_that("the Azimuth model loader derives pathway compounds from reactions", {
  dir <- withr::local_tempdir()
  m   <- read_mummichog_model_pathways(write_model_json(dir))

  expect_identical(m$format, "azimuth")
  expect_setequal(names(m$pathways),
                  c("Purine metabolism", "Fatty acid degradation", "Glycolysis"))
  expect_setequal(m$pathways[["Purine metabolism"]],
                  c("C00020", "C00362", "C01367"))
  expect_setequal(m$pathways[["Fatty acid degradation"]], c("C00020", "C00186"))
  expect_identical(unname(m$compound_kegg[["C00020"]]), "C00020")
  expect_match(unname(m$compound_names[["C00020"]]), "^AMP")
})

test_that("the mummichog-2 native model shape is read too", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "native.json")
  jsonlite::write_json(list(
    metabolic_pathways = list(
      list(id = "x1", name = "Glycolysis", cpds = list("C00031", "C00095"))),
    dict_cpds_def = list(C00031 = "D-Glucose", C00095 = "D-Fructose")
  ), path, auto_unbox = TRUE)

  m <- read_mummichog_model_pathways(path)
  expect_identical(m$format, "mummichog2")
  expect_setequal(m$pathways[["Glycolysis"]], c("C00031", "C00095"))
  expect_identical(unname(m$compound_names[["C00031"]]), "D-Glucose")
})

test_that("a malformed model fails loudly rather than degrading silently", {
  dir <- withr::local_tempdir()
  bad <- file.path(dir, "bad.json")
  jsonlite::write_json(list(something_else = 1), bad, auto_unbox = TRUE)

  expect_error(read_mummichog_model_pathways(bad), "neither 'list_of_pathways'")
  expect_error(read_mummichog_model_pathways(file.path(dir, "nope.json")),
               "not found")
})

test_that("the built-in human_mfn model reports itself unavailable, not broken", {
  cfg <- list(modes = list(metabolomics = list(
    organism = "Homo sapiens",
    enrichment = list(mummichog = list(enabled = TRUE)))))
  expect_message(res <- mmc_load_model_pathways(cfg), "human_mfn")
  expect_null(res)
})


# ---------------------------------------------------------------------------
# EmpiricalCompound readers
# ---------------------------------------------------------------------------

test_that("EC candidates are read in full, with names paired positionally", {
  root  <- withr::local_tempdir()
  files <- build_evidence_tree(root)
  cand  <- read_mummichog_ec_candidates(files)

  # every candidate of the multi-candidate EC survives (no "best guess" pick)
  e196 <- cand[cand$EID == "E196", ]
  expect_setequal(e196$compound_id, c("C00020", "C00362", "C01367"))
  expect_identical(e196$compound_name[e196$compound_id == "C00362"], "dGMP")
  expect_equal(nrow(cand), 6L)
})

test_that("EC features keep every underlying signal and recover the adduct", {
  root  <- withr::local_tempdir()
  files <- build_evidence_tree(root)
  feat  <- read_mummichog_ec_features(files)

  e196 <- feat[feat$EID == "E196", ]
  expect_equal(nrow(e196), 2L)                       # both features preserved
  expect_setequal(e196$feature_id, c("feat_1", "feat_2"))
  expect_equal(e196$mz[e196$feature_id == "feat_1"], 348.07)
  # per-feature ion parsed out of the EC's str_row_ion string
  expect_identical(e196$adduct[e196$feature_id == "feat_1"], "M+H[1+]")
  expect_identical(e196$adduct[e196$feature_id == "feat_2"], "M+Na[1+]")
  expect_identical(feat$adduct[feat$feature_id == "feat_4"], "M-H[-]")
})

test_that("EC readers return zero rows (not an error) when the files are absent", {
  empty <- withr::local_tempdir()
  expect_equal(nrow(read_mummichog_ec_candidates(list_mummichog_files(empty))), 0L)
  expect_equal(nrow(read_mummichog_ec_features(list_mummichog_files(empty))), 0L)
  expect_equal(nrow(read_mummichog_ec_candidates(character(0))), 0L)
})


# ---------------------------------------------------------------------------
# annotation normalisation (schema-adaptive)
# ---------------------------------------------------------------------------

test_that("annotations normalise with an identification-level system", {
  a <- normalize_metab_annotation(row_data_with_levels())

  expect_identical(a$original_annotation_name[1], "AMP")
  expect_identical(a$original_annotation_id[1], "C00020")
  expect_identical(a$original_annotation_id_type[1], "KEGG")
  expect_identical(a$original_annotation_confidence[1], "Level 1")
  expect_identical(a$original_annotation_confidence[3], "Level 2")
  # unannotated feature stays NA everywhere — nothing is inferred
  expect_true(is.na(a$original_annotation_name[2]))
  expect_true(is.na(a$original_annotation_id[2]))
})

test_that("annotations normalise without any level system", {
  a <- normalize_metab_annotation(row_data_no_levels())

  expect_identical(a$original_annotation_name[3], "D-glucose")
  expect_true(all(is.na(a$original_annotation_confidence)))
  expect_true(all(is.na(a$original_annotation_id)))
})

test_that("a textual confidence column is passed through verbatim", {
  rd <- data.frame(feature_id = c("f1", "f2"),
                   Name = c("AMP", "dGMP"),
                   Confidence = c("high (MS2)", "putative"),
                   stringsAsFactors = FALSE)
  a <- normalize_metab_annotation(rd)
  expect_identical(a$original_annotation_confidence, c("high (MS2)", "putative"))
})

test_that("annotation columns are found regardless of capitalisation", {
  rd <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                   feature_id = "f1", metabolite = "AMP", `kegg id` = "C00020",
                   `MSI level` = 2L)
  a <- normalize_metab_annotation(rd)
  expect_identical(a$original_annotation_name, "AMP")
  expect_identical(a$original_annotation_id, "C00020")
  expect_identical(a$original_annotation_confidence, "Level 2")
})

test_that("row_data with no annotation columns yields no usable annotation", {
  rd <- data.frame(feature_id = c("f1", "f2"), `m/z` = c(100, 200),
                   check.names = FALSE, stringsAsFactors = FALSE)
  a  <- normalize_metab_annotation(rd)
  expect_equal(nrow(a), 2L)
  expect_true(all(is.na(a$original_annotation_name)))
  expect_equal(nrow(normalize_metab_annotation(NULL)), 0L)
})


# ---------------------------------------------------------------------------
# annotation agreement rules
# ---------------------------------------------------------------------------

test_that("agreement matches on KEGG compound id when both sides have one", {
  expect_identical(
    mmc_annotation_agreement("C00020", "AMP",
                             candidate_ids = "C00020",
                             candidate_kegg = "C00020",
                             candidate_names = "AMP"),
    "Match")
})

test_that("agreement matches on a normalised name when ids are unavailable", {
  # different case and punctuation, same compound; no ids on either side
  expect_identical(
    mmc_annotation_agreement(NA_character_, "d-glucose",
                             candidate_ids = "glc",
                             candidate_kegg = NA_character_,
                             candidate_names = "D-Glucose;Glucose"),
    "Match")
  # synonyms inside the model's ";"-separated name list also count
  expect_identical(
    mmc_annotation_agreement(NA_character_, "Glucose",
                             candidate_ids = "glc",
                             candidate_kegg = NA_character_,
                             candidate_names = "D-Glucose;Glucose"),
    "Match")
})

test_that("agreement reports a conflict for a different metabolite", {
  # comparable KEGG ids that disagree
  expect_identical(
    mmc_annotation_agreement("C00031", "D-Glucose",
                             candidate_ids = "C00186",
                             candidate_kegg = "C00186",
                             candidate_names = "(S)-Lactate;L-Lactate"),
    "Conflict")
  # comparable names that disagree
  expect_identical(
    mmc_annotation_agreement(NA_character_, "D-Glucose",
                             candidate_ids = "lac",
                             candidate_kegg = NA_character_,
                             candidate_names = "L-Lactate"),
    "Conflict")
})

test_that("agreement is 'Not assessed' when there is nothing to compare", {
  # no original annotation at all
  expect_identical(
    mmc_annotation_agreement(NA_character_, NA_character_,
                             candidate_ids = "C00020",
                             candidate_kegg = "C00020",
                             candidate_names = "AMP"),
    "Not assessed")
  # an original name but no candidate ids or names to compare against
  expect_identical(
    mmc_annotation_agreement(NA_character_, "AMP",
                             candidate_ids = "x1",
                             candidate_kegg = NA_character_,
                             candidate_names = NA_character_),
    "Not assessed")
  # no candidates at all
  expect_identical(
    mmc_annotation_agreement("C00020", "AMP", character(0), character(0),
                             character(0)),
    "Not assessed")
})

test_that("agreement never matches on mass-like evidence alone", {
  # identical formula/mass is irrelevant: different compound, different name/id
  expect_identical(
    mmc_annotation_agreement("C00031", "D-Glucose",
                             candidate_ids = "C00095",
                             candidate_kegg = "C00095",
                             candidate_names = "D-Fructose"),
    "Conflict")
})


# ---------------------------------------------------------------------------
# end-to-end evidence trace
# ---------------------------------------------------------------------------

evidence_fixture <- function(row_data = row_data_with_levels()) {
  root  <- withr::local_tempdir(.local_envir = parent.frame())
  files <- build_evidence_tree(root)
  model <- read_mummichog_model_pathways(write_model_json(root))
  pw    <- read_mummichog_pathways(files)
  annot <- normalize_metab_annotation(row_data)
  list(files = files, model = model, pathways = pw, annot = annot,
       evidence = build_mummichog_pathway_evidence(pw, files, model, annot,
                                                   p_cutoff = 0.05))
}

test_that("pathways resolve to their overlapping EmpiricalCompounds", {
  f  <- evidence_fixture()
  ev <- f$evidence

  expect_setequal(unique(ev$ec_table$Pathway),
                  c("Fatty acid degradation", "Purine metabolism", "Glycolysis"))
  fad <- ev$ec_table[ev$ec_table$Pathway == "Fatty acid degradation", ]
  expect_setequal(fad$EmpiricalCompound, c("E196", "E3"))
  gly <- ev$ec_table[ev$ec_table$Pathway == "Glycolysis", ]
  expect_setequal(gly$EmpiricalCompound, c("E2", "E4"))
})

test_that("an EC with one candidate reports that candidate", {
  ev <- evidence_fixture()$evidence
  row <- ev$ec_table[ev$ec_table$Pathway == "Glycolysis" &
                       ev$ec_table$EmpiricalCompound == "E2", ]
  expect_identical(row[["Pathway-matching candidate(s)"]], "C00031")
  expect_identical(row[["Candidate KEGG ID(s)"]], "C00031")
  expect_equal(row[["# Features"]], 1)
})

test_that("only the candidate belonging to the selected pathway is reported", {
  ev <- evidence_fixture()$evidence
  # E196 carries AMP / dGMP / 3'-AMP, but only AMP is in Fatty acid degradation
  row <- ev$ec_table[ev$ec_table$Pathway == "Fatty acid degradation" &
                       ev$ec_table$EmpiricalCompound == "E196", ]
  expect_identical(row[["Pathway-matching candidate(s)"]], "C00020")
  expect_false(grepl("C00362|C01367", row[["Pathway-matching candidate(s)"]]))
})

test_that("several candidates of one EC in one pathway are all kept, EC counted once", {
  ev <- evidence_fixture()$evidence
  rows <- ev$ec_table[ev$ec_table$Pathway == "Purine metabolism", ]

  expect_equal(nrow(rows), 1L)                 # the EC is counted once
  cands <- strsplit(rows[["Pathway-matching candidate(s)"]], "; ")[[1]]
  expect_setequal(cands, c("C00020", "C00362", "C01367"))   # all three kept
  keggs <- strsplit(rows[["Candidate KEGG ID(s)"]], "; ")[[1]]
  expect_setequal(keggs, c("C00020", "C00362", "C01367"))
})

test_that("every measured feature behind an EC is preserved", {
  ev <- evidence_fixture()$evidence
  ft <- ev$feature_table[ev$feature_table$Pathway == "Fatty acid degradation" &
                           ev$feature_table$EmpiricalCompound == "E196", ]

  expect_equal(nrow(ft), 2L)
  expect_setequal(ft$Feature, c("feat_1", "feat_2"))
  # per-feature statistics carried through, nothing summed or averaged
  expect_equal(sort(ft[["m/z"]]), c(348.07, 370.05))
  expect_setequal(ft[["Adduct/ion"]], c("M+H[1+]", "M+Na[1+]"))
  expect_equal(ft[["Feature p-value"]][ft$Feature == "feat_1"], 0.002)
  expect_equal(ft[["Feature statistic"]][ft$Feature == "feat_1"], 4.1)
  expect_true(all(ft$Significant))                    # both p < 0.05
})

test_that("agreement is traced per feature and rolled up per EC", {
  ev <- evidence_fixture()$evidence
  ft <- ev$feature_table

  # KEGG-based match (feat_1 annotated as AMP / C00020)
  expect_identical(ft$Agreement[ft$Feature == "feat_1" &
                                  ft$Pathway == "Fatty acid degradation"],
                   "Match")
  # unannotated feature of the same EC
  expect_identical(ft$Agreement[ft$Feature == "feat_2" &
                                  ft$Pathway == "Fatty acid degradation"],
                   "Not assessed")
  # EC-level roll-up: one matching feature is enough
  expect_identical(
    ev$ec_table$Agreement[ev$ec_table$EmpiricalCompound == "E196" &
                            ev$ec_table$Pathway == "Fatty acid degradation"],
    "Match")

  # name-based match with no ids available (feat_3 "D-glucose" vs D-Glucose)
  expect_identical(ft$Agreement[ft$Feature == "feat_3"], "Match")
  # conflict: feat_4 is annotated as glucose, the candidate is lactate
  expect_identical(ft$Agreement[ft$Feature == "feat_4"], "Conflict")
  # missing annotation -> Not assessed
  expect_identical(ft$Agreement[ft$Feature == "feat_5"], "Not assessed")

  # conflicts are reported, never dropped from the mummichog result
  expect_true("E3" %in% ev$ec_table$EmpiricalCompound)
  fad <- ev$pathway_summary$Pathway == "Fatty acid degradation"
  expect_equal(ev$pathway_summary[["ECs Conflict"]][fad], 1)
  expect_equal(ev$pathway_summary[["ECs Match"]][fad], 1)
  expect_equal(ev$pathway_summary[["ECs Mixed"]][fad], 0)

  # EC-level feature counts describe the features behind the roll-up
  e196 <- ev$ec_table$EmpiricalCompound == "E196" &
    ev$ec_table$Pathway == "Fatty acid degradation"
  expect_equal(ev$ec_table$n_match[e196], 1)          # feat_1
  expect_equal(ev$ec_table$n_conflict[e196], 0)
  expect_equal(ev$ec_table$n_not_assessed[e196], 1)   # feat_2
  expect_equal(ev$ec_table[["# Features"]][e196],
               ev$ec_table$n_match[e196] + ev$ec_table$n_conflict[e196] +
                 ev$ec_table$n_not_assessed[e196])
})

test_that("the pathway summary passes the ORA statistics through untouched", {
  f  <- evidence_fixture()
  s  <- f$evidence$pathway_summary
  fad <- s[s$Pathway == "Fatty acid degradation", ]

  expect_equal(fad$Overlap, 2)
  expect_equal(fad[["Detected pathway size"]], 2)
  expect_equal(fad[["Enrichment ratio"]], 1)
  expect_equal(fad[["p.value"]], 0.01)
  expect_equal(fad[["Supporting ECs"]], 2)
  expect_equal(fad[["Supporting features"]], 3)   # feat_1, feat_2, feat_4
  # the two grains are separately reported and never conflated
  expect_equal(fad[["ECs Match"]] + fad[["ECs Conflict"]] +
                 fad[["ECs Mixed"]] + fad[["ECs Not assessed"]],
               fad[["Supporting ECs"]])
  expect_equal(fad[["features Match"]] + fad[["features Conflict"]] +
                 fad[["features Not assessed"]],
               fad[["Supporting features"]])
  # sorted by ascending empirical p-value
  expect_equal(s[["p.value"]], sort(s[["p.value"]]))
})

test_that("a dataset with no annotation levels still traces evidence", {
  ev <- evidence_fixture(row_data_no_levels())$evidence

  expect_true(all(is.na(ev$ec_table[["Annotation confidence"]])))
  # feat_1 has no KEGG here, so the name comparison decides
  expect_identical(
    ev$feature_table$Agreement[ev$feature_table$Feature == "feat_1" &
                                 ev$feature_table$Pathway == "Purine metabolism"],
    "Match")
})

test_that("evidence returns NULL when its inputs are genuinely unavailable", {
  f <- evidence_fixture()

  # no model -> no pathway membership -> no evidence (not an error)
  expect_null(build_mummichog_pathway_evidence(f$pathways, f$files, NULL, f$annot))
  # no EC tables among the files
  empty_dir <- withr::local_tempdir()
  expect_message(
    res <- build_mummichog_pathway_evidence(f$pathways,
                                            list_mummichog_files(empty_dir),
                                            f$model, f$annot),
    "no EmpiricalCompound tables")
  expect_null(res)
  # empty pathway table
  expect_null(build_mummichog_pathway_evidence(data.frame(), f$files, f$model,
                                               f$annot))
  expect_null(build_mummichog_pathway_evidence(NULL, f$files, f$model, f$annot))
})

test_that("evidence works with a zero-row annotation table", {
  f  <- evidence_fixture()
  ev <- build_mummichog_pathway_evidence(f$pathways, f$files, f$model,
                                         normalize_metab_annotation(NULL))
  expect_false(is.null(ev))
  expect_true(all(ev$feature_table$Agreement == "Not assessed"))
  expect_true(all(ev$ec_table$Agreement == "Not assessed"))
})


# ---------------------------------------------------------------------------
# file grouping shared with the pathway read-back
# ---------------------------------------------------------------------------

test_that("group_mummichog_files_by_contrast recovers the original label", {
  root  <- withr::local_tempdir()
  files <- build_evidence_tree(root, contrast_dir = "HL_vs_LL")
  groups <- group_mummichog_files_by_contrast(files)

  expect_named(groups, "HL_vs_LL")
  expect_true(any(basename(groups[["HL_vs_LL"]]) ==
                    "ListOfEmpiricalCompounds.tsv"))
  expect_length(group_mummichog_files_by_contrast(character(0)), 0)
})


# ---------------------------------------------------------------------------
# agreement roll-up: four states with feature-level counts
# ---------------------------------------------------------------------------

test_that("the EC roll-up implements the four-state rule", {
  roll <- .mmc_summarise_agreement

  expect_identical(roll(c("Match", "Match")), "Match")
  expect_identical(roll(c("Conflict", "Conflict")), "Conflict")
  expect_identical(roll(c("Match", "Conflict")), "Mixed")
  expect_identical(roll(c("Conflict", "Match")), "Mixed")
  # a not-assessed feature never overrides assessed evidence
  expect_identical(roll(c("Match", "Not assessed")), "Match")
  expect_identical(roll(c("Conflict", "Not assessed")), "Conflict")
  expect_identical(roll(c("Not assessed", "Not assessed")), "Not assessed")
  expect_identical(roll("Not assessed"), "Not assessed")
  # all three together is still Mixed
  expect_identical(roll(c("Match", "Conflict", "Not assessed")), "Mixed")
})

test_that("a Mixed EmpiricalCompound is reported as Mixed, with its counts", {
  # E196 is formed by two features; annotate one as the pathway-matching
  # candidate (AMP / C00020) and the other as a different metabolite, so the EC
  # rolls up to Mixed rather than hiding the disagreement behind Match.
  rd <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    feature_id = paste0("feat_", 1:5),
    Name       = c("AMP", "D-Glucose", NA, NA, NA),
    KEGG       = c("C00020", "C00031", NA, NA, NA),
    identification_level = c(1L, 1L, NA, NA, NA)
  )
  ev <- evidence_fixture(rd)$evidence
  row <- ev$ec_table[ev$ec_table$EmpiricalCompound == "E196" &
                       ev$ec_table$Pathway == "Fatty acid degradation", ]

  expect_identical(row$Agreement, "Mixed")
  expect_equal(row$n_match, 1)          # feat_1 -> AMP
  expect_equal(row$n_conflict, 1)       # feat_2 -> glucose, not AMP
  expect_equal(row$n_not_assessed, 0)
  # per-feature verdicts stay untouched and visible
  ft <- ev$feature_table[ev$feature_table$EmpiricalCompound == "E196" &
                           ev$feature_table$Pathway == "Fatty acid degradation", ]
  expect_setequal(ft$Agreement, c("Match", "Conflict"))
  # and the pathway summary counts it in the Mixed column, not Match
  fad <- ev$pathway_summary$Pathway == "Fatty acid degradation"
  expect_equal(ev$pathway_summary[["ECs Mixed"]][fad], 1)
  expect_equal(ev$pathway_summary[["ECs Match"]][fad], 0)
})

test_that("EC and feature agreement counts stay at their own grains", {
  ev <- evidence_fixture()$evidence

  # EC-level state counts sum to the EC count, per pathway
  s <- ev$pathway_summary
  ec_states <- s[["ECs Match"]] + s[["ECs Conflict"]] + s[["ECs Mixed"]] +
    s[["ECs Not assessed"]]
  expect_equal(ec_states, s[["Supporting ECs"]])

  # feature-level counts sum to the feature count, per pathway
  feat_states <- s[["features Match"]] + s[["features Conflict"]] +
    s[["features Not assessed"]]
  expect_equal(feat_states, s[["Supporting features"]])

  # and neither equals the overlap / detected pathway size by construction:
  # Fatty acid degradation has 2 supporting ECs but 3 supporting features
  fad <- s[s$Pathway == "Fatty acid degradation", ]
  expect_equal(fad[["Supporting ECs"]], 2)
  expect_equal(fad[["Supporting features"]], 3)

  # per-EC counts sum to that EC's feature count
  expect_equal(ev$ec_table$n_match + ev$ec_table$n_conflict +
                 ev$ec_table$n_not_assessed,
               ev$ec_table[["# Features"]])
})


# ---------------------------------------------------------------------------
# pathway-specific candidate intersection (the shape the review asked for)
# ---------------------------------------------------------------------------

test_that("two candidates of one EC in one pathway are kept, EC counted once", {
  # EC "Ea" carries candidates C1 and C2, both of which belong to P1.
  # EC "Eb" carries C3, which also belongs to P1.
  # => P1 has TWO supporting ECs (not three), and Ea lists both candidates.
  root   <- withr::local_tempdir()
  tables <- file.path(root, "mummichog_pinned", "C1", "v2", "1.run", "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "pathway\toverlap_size\tpathway_size\tp-value\toverlap_EmpiricalCompounds (id)\toverlap_features (id)\toverlap_features (name)",
    "P1\t2\t3\t0.01\tEa,Eb\tC9\tirrelevant-face-compound"
  ), file.path(tables, "mcg_pathwayanalysis_C1.tsv"))
  writeLines(c(
    "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names",
    "Ea\trow1\trow1_M+H[1+]\tC1;C2\tone$two",
    "Eb\trow2\trow2_M+H[1+]\tC3\tthree"
  ), file.path(tables, "ListOfEmpiricalCompounds.tsv"))
  writeLines(c(
    "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    "row1\tEa\trow1_M+H[1+]\tC1;C2\tone$two\trow1\t100\t1\t0.01\t2\tf1",
    "row2\tEb\trow2_M+H[1+]\tC3\tthree\trow2\t200\t2\t0.02\t1\tf2"
  ), file.path(tables, "userInput_to_EmpiricalCompounds.tsv"))

  model_path <- file.path(root, "m.json")
  jsonlite::write_json(list(
    metabolic_pathways = list(
      list(id = "p1", name = "P1", cpds = list("C1", "C2", "C3"))),
    dict_cpds_def = list(C1 = "one", C2 = "two", C3 = "three", C9 = "nine")
  ), model_path, auto_unbox = TRUE)

  files <- list_mummichog_files(root)
  ev <- build_mummichog_pathway_evidence(
    read_mummichog_pathways(files), files,
    read_mummichog_model_pathways(model_path),
    normalize_metab_annotation(NULL))

  ec <- ev$ec_table[ev$ec_table$Pathway == "P1", ]
  expect_equal(nrow(ec), 2L)                          # Ea and Eb, not 3 rows
  expect_setequal(ec$EmpiricalCompound, c("Ea", "Eb"))
  expect_equal(ev$pathway_summary[["Supporting ECs"]], 2)

  # Ea keeps BOTH of its pathway-matching candidates
  ea <- ec[ec$EmpiricalCompound == "Ea", ]
  expect_setequal(strsplit(ea[["Pathway-matching candidate(s)"]], "; ")[[1]],
                  c("C1", "C2"))
  expect_equal(ea[["# Features"]], 1)                 # one EC, one feature

  # mummichog's own overlap_features (id) column (C9, the arbitrary
  # face_compound) is irrelevant and must not appear anywhere
  expect_false(any(grepl("C9", unlist(ec), fixed = TRUE)))
  expect_false(any(grepl("C9", unlist(ev$feature_table), fixed = TRUE)))
})
