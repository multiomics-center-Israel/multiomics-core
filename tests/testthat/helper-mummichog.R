# tests/testthat/helper-mummichog.R
#
# Shared synthetic fixtures for the mummichog test files. Kept here (rather than
# duplicated per file) because the GSEA, parity and evidence suites all need the
# same faithful v2 result tree; testthat sources helper*.R before the tests.
#
# Everything is SYNTHETIC. Compound ids are real KEGG accessions only so the
# ID-comparison paths are exercised realistically; pathway membership is
# invented for the tests and is not a biological claim.

# A larger synthetic mummichog tree + model, big enough for a real fgsea run:
# 24 EmpiricalCompounds, each with one feature and one candidate compound, and
# three pathways of 8 compounds each. Pathway "PW_pos" holds the 8 highest
# statistics, "PW_neg" the 8 lowest, "PW_mid" the middle 8.
build_gsea_fixture <- function(root, contrast_dir = "HL_vs_LL") {
  n     <- 24L
  cpds  <- sprintf("C%05d", seq_len(n))
  eids  <- sprintf("E%03d", seq_len(n))
  feats <- sprintf("feat_%02d", seq_len(n))
  # statistics from +12 down to -12 (never 0), so the ranking is unambiguous
  stat  <- c(seq(12, 1), seq(-1, -12))
  mzs   <- 100 + seq_len(n)

  tables <- file.path(root, "mummichog_pinned", contrast_dir, "v2",
                      "1700000000.1.run", "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names",
    sprintf("%s\trow%d\trow%d_M+H[1+]\t%s\tcpd %s", eids, seq_len(n),
            seq_len(n), cpds, cpds)
  ), file.path(tables, "ListOfEmpiricalCompounds.tsv"))

  writeLines(c(
    "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    sprintf("row%d\t%s\trow%d_M+H[1+]\t%s\tcpd %s\trow%d\t%s\t1.0\t0.05\t%s\t%s",
            seq_len(n), eids, seq_len(n), cpds, cpds, seq_len(n),
            format(mzs), format(stat), feats)
  ), file.path(tables, "userInput_to_EmpiricalCompounds.tsv"))

  writeLines(c(
    "pathway\toverlap_size\tpathway_size\tp-value\toverlap_EmpiricalCompounds (id)\toverlap_features (id)\toverlap_features (name)",
    sprintf("PW_pos\t3\t8\t0.01\t%s\t\t", paste(eids[1:3], collapse = ",")),
    sprintf("PW_neg\t3\t8\t0.02\t%s\t\t", paste(eids[22:24], collapse = ","))
  ), file.path(tables, "mcg_pathwayanalysis_HL_vs_LL.tsv"))

  readr::write_tsv(
    data.frame(dir = contrast_dir, contrast = "HL_vs_LL",
               stringsAsFactors = FALSE),
    file.path(root, "mummichog_pinned", "contrasts.tsv"))

  # Model: three disjoint pathways over the same 24 compounds, in native
  # mummichog-2 shape (the compact one).
  model_path <- file.path(root, "model.json")
  jsonlite::write_json(list(
    metabolic_pathways = list(
      list(id = "p1", name = "PW_pos", cpds = as.list(cpds[1:8])),
      list(id = "p2", name = "PW_mid", cpds = as.list(cpds[9:16])),
      list(id = "p3", name = "PW_neg", cpds = as.list(cpds[17:24]))
    ),
    dict_cpds_def = stats::setNames(as.list(paste("cpd", cpds)), cpds)
  ), model_path, auto_unbox = TRUE)

  list(
    files    = list_mummichog_files(root),
    model    = read_mummichog_model_pathways(model_path),
    de_table = data.frame(feature_id = feats, logFC = stat / 4,
                          statistic = stat, P.Value = rep(0.05, n),
                          stringsAsFactors = FALSE),
    eids     = eids,
    stat     = stat
  )
}

