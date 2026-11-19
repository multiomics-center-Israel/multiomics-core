# R/domain/metabolomics/06l_mummichog_ma_ecs.R
#
# MetaboAnalyst-compatible mummichog, component 3: empirical compounds (ECs)
# and their pathway sets.
#
# An EC groups the matches that MetaboAnalyst takes to be one metabolite. The
# grouping works on the stacked match table of component 2, so in mixed mode a
# positive-ion and a negative-ion feature can land in ONE EC: matches are split
# by COMPOUND only — never by polarity — and a compound's matches stay together
# unless their retention times are far apart. The primary-ion rule then asks
# only that SOME member of the EC is a primary ion, so a negative [M-H]- can
# keep an EC whose positive member is a non-primary adduct.
#
# ===========================================================================
# BEHAVIOURAL REFERENCE — pinned upstream version
# ===========================================================================
#     repo:   xia-lab/MetaboAnalystR
#     commit: 398476ae2a0c996e390925fa62adc46b5ecde334
#     file:   R/peaks_to_function.R
#
#   mmc_ma_build_ecs()     .search.compoundLib()   L1806-L2075 (the
#                          retention-time + "v2" branch: ECs, then the
#                          expression and lookup dictionaries); helpers
#                          Convert2Dictionary() L4564, primary_ions L4417,
#                          melt.list() (R/depends_refs.R L22)
#   mmc_ma_pathway_ecs()   .setup.psea.library()   L1537-L1546 (minLib) and
#                          L1564-L1583 (pathway -> ECs)
#
# The steps, as upstream runs them:
#   1. Split the matches by compound. A compound with several matches is
#      ordered by RT and cut between neighbours whose RT differs by more than
#      rt_tol AND whose RT rank differs by more than rt_tol_rank
#      (= peak-list size * rt_frac). Each piece is an initial EC (EC1, EC2...).
#   2. Rows with the same m/z, ion form and RT (one feature read as one ion
#      form, matched to several isobaric compounds) are merged into one row
#      that carries all their ECs ("EC3;EC9") and compounds.
#   3. Primary-ion rule: an EC (as labelled after step 2) is kept when any of
#      its rows is a primary ion (M+H, M+Na, M-H2O+H, M-H, M-2H, M-H2O-H).
#   4. Per compound, ECs sharing a member are joined (connected components).
#   5. ECs are renamed EC0001, EC0002... in first-appearance order, and the
#      lookup dictionaries are built.
#
# Reproduced deliberately, including upstream quirks:
#   * groups, merges and dictionaries use split() / aggregate() / merge(),
#     which order by the CURRENT LOCALE's collation, so EC labels (not EC
#     membership) depend on the locale, exactly as upstream's do;
#   * step 5 renames with regular expressions built from EC labels, where "|"
#     acts as alternation — upstream's own matching, kept as is.
#
# Divergence, documented: upstream stops with an error when NO compound is
# split by retention time (its unlist(recursive = FALSE) then returns a
# character vector, which its melt() turns into a frame without the EC labels).
# We complete the grouping in that case with the same rules. No upstream result
# exists to compare it with.
#
# Not reproduced: upstream's console messages and its result-file writes.


# Upstream's primary ions (peaks_to_function.R L4417-L4418).
.MMC_MA_PRIMARY_IONS <- c("M+H[1+]", "M+Na[1+]", "M-H2O+H[1+]", "M-H[-]",
                          "M-2H[2-]", "M-H2O-H[-]", "M+H [1+]", "M+Na [1+]",
                          "M-H2O+H [1+]", "M-H [1-]", "M-2H [2-]",
                          "M-H2O-H [1-]")

#' Two-column table to a named list (MetaboAnalystR's Convert2Dictionary)
#'
#' Groups column 2 by column 1. Names are sorted as `split()` sorts them; every
#' id keeps all its values, in table order.
#'
#' @param data A data frame or matrix with at least two columns.
#' @return Named list: id -> vector of values.
#' @noRd
.mmc_ma_dict <- function(data) {
  all_ids <- data[, 1]
  dup <- duplicated(all_ids)
  if (sum(dup) > 0) {
    uniq_list <- split(data[!dup, 2], all_ids[!dup])
    uniq_id_list <- names(uniq_list)
    for (id in unique(all_ids[dup])) {
      uniq_list[[which(uniq_id_list == id)]] <- data[which(all_ids == id), 2]
    }
    return(uniq_list)
  }
  split(data[, 2], data[, 1])
}

#' Group matches into empirical compounds the way MetaboAnalystR does
#'
#' Port of the retention-time ("v2") branch of `.search.compoundLib()` after
#' matching: initial ECs by compound and RT, merging of rows shared by isobaric
#' compounds, the per-EC primary-ion rule, joining overlapping ECs per compound,
#' EC naming, and the dictionaries the pathway statistics read. Positive- and
#' negative-ion matches are grouped together; see the file header.
#'
#' @param matched Match table from `mmc_ma_match_peaks()`.
#' @param prep Peak list from `mmc_ma_prepare_peaks()` (RT tolerances, the
#'   statistic per m/z, the reference list size).
#' @param input_mzlist Significant m/z keys from `mmc_ma_input_mzlist()`.
#' @param force_primary_ion Apply the primary-ion rule (upstream default TRUE,
#'   `primary_ion = "yes"`).
#' @return A list:
#'   \describe{
#'     \item{matched_res}{One row per (feature, compound, ion form) kept, with
#'       `Query.Mass`, `Matched.Compound`, `Matched.Form`, `Retention.Time`,
#'       `Mass.Diff`, `Empirical.Compound` (upstream's `mum_res`).}
#'     \item{initial_ecs}{The initial ECs before merging and filtering, as
#'       upstream saves them (`ec`, `mz`, `form`, `rt`, `cpd`, `str_row_inx`).}
#'     \item{ec_exp_dict, ec_exp, cpd_exp_dict}{Statistic per EC (all
#'       members / their max) and per compound.}
#'     \item{cpd_ecpd_dict, ecpd_cpd_dict, mz2cpd_dict, mz2ec_dict,
#'       ec2mz_dict, cpd_form_dict}{Lookup dictionaries.}
#'     \item{cpd_ecpd_counts}{Per compound, its ECs after joining.}
#'     \item{input_ecpdlist}{ECs holding at least one significant feature.}
#'     \item{total_matched_ecpds}{All ECs: the statistical background.}
#'   }
#' @examples
#' \dontrun{
#' ecs <- mmc_ma_build_ecs(matched, prep, mmc_ma_input_mzlist(prep, 0.05))
#' }
mmc_ma_build_ecs <- function(matched, prep, input_mzlist,
                             force_primary_ion = TRUE) {
  matched_res <- matched
  rt_tol <- prep$rt_tol
  rt_tol_rank <- length(prep$ref_mzlist) * prep$rt_frac

  # 1. initial ECs: by compound, cut by RT ------------------------------------
  empirical_cpd_list <- split(matched_res[, c(1, 3, 5, 6)], matched_res[, 2])
  empirical_cpds2cpds <- lapply(seq_along(empirical_cpd_list), function(i) {
    ec_data <- empirical_cpd_list[[i]]
    cpd_name <- names(empirical_cpd_list)[i]
    mzs <- ec_data$Query.Mass
    ions <- ec_data$Matched.Form
    rts_char <- ec_data$Retention.Time
    rt_rank_char <- ec_data$RT.Rank
    if (length(mzs) == 1) {
      return(paste0(mzs, ";", ions, ";", rts_char, ";", cpd_name))
    }
    rts_num <- as.numeric(rts_char)
    rt_ranks_num <- as.numeric(rt_rank_char)
    ec_names <- paste0(mzs, ";", ions, ";", rts_char, ";", cpd_name)
    sort_idx <- order(rts_num)
    rts_sorted <- rts_num[sort_idx]
    ranks_sorted <- rt_ranks_num[sort_idx]
    names_sorted <- ec_names[sort_idx]
    split_inx <- c(0, cumsum((abs(diff(rts_sorted)) > rt_tol) &
                             (abs(diff(ranks_sorted)) > rt_tol_rank)))
    if (length(unique(split_inx)) > 1) {
      return(split(names_sorted, split_inx))
    }
    paste0(names_sorted, collapse = "__")
  })
  names(empirical_cpds2cpds) <- names(empirical_cpd_list)

  initial_ecs <- unlist(empirical_cpds2cpds, recursive = FALSE)
  # Upstream errors here when nothing was split (see the file header).
  if (!is.list(initial_ecs)) initial_ecs <- as.list(initial_ecs)
  names(initial_ecs) <- paste0("EC", seq_along(initial_ecs))

  # melt.list(): one row per string, labelled with its EC.
  try_value <- unlist(lapply(initial_ecs, as.character), use.names = FALSE)
  try_l1 <- rep(names(initial_ecs), lengths(initial_ecs))
  try2 <- strsplit(try_value, split = "__", fixed = TRUE)
  try2_df <- data.frame(value = unlist(try2),
                        L1 = rep(try_l1, sapply(try2, length)))
  info <- strsplit(as.character(try2_df[, 1]), split = ";")
  df_ecs <- data.frame(ec = as.character(try2_df[, 2]),
                       mz = sapply(info, `[[`, 1), form = sapply(info, `[[`, 2),
                       rt = sapply(info, `[[`, 3), cpd = sapply(info, `[[`, 4),
                       stringsAsFactors = FALSE)
  df_ecs$str_row_inx <- paste(df_ecs$mz, df_ecs$form, df_ecs$rt, sep = "___")

  # 2. merge rows shared by isobaric compounds ---------------------------------
  merged_ecs <- stats::aggregate(. ~ str_row_inx, df_ecs, paste, collapse = ";")
  cols_to_clean <- c("mz", "form", "rt")
  merged_ecs[cols_to_clean] <- lapply(merged_ecs[cols_to_clean], function(col) {
    vapply(strsplit(as.character(col), ";", fixed = TRUE),
           function(x) paste(unique(x), collapse = ";"), FUN.VALUE = character(1))
  })

  # 3. primary-ion rule, per EC label ------------------------------------------
  if (isTRUE(force_primary_ion)) {
    ecs <- unique(merged_ecs$ec)
    new_info <- lapply(ecs, function(x) {
      rows <- merged_ecs[which(merged_ecs$ec == x), ]
      if (length(intersect(rows$form, .MMC_MA_PRIMARY_IONS)) > 0) rows else NULL
    })
    kept <- do.call(args = new_info, what = rbind)
    if (is.null(kept)) {
      .mmc_stop("no empirical compound contains a primary ion (",
                paste(.MMC_MA_PRIMARY_IONS[1:6], collapse = ", "), "), so ",
                "none survives force_primary_ion. Upstream stops here too.")
    }
    final_ecs <- kept[, -1]
  } else {
    final_ecs <- merged_ecs[, -1]
  }
  colnames(final_ecs) <- c("Empirical.Compound", "Query.Mass", "Matched.Form",
                           "Retention.Time", "Matched.Compound")

  # long format: one row per (EC row, compound)
  cpd_split <- strsplit(as.character(final_ecs$Matched.Compound), ";", fixed = TRUE)
  reps <- pmax(lengths(cpd_split))
  df2 <- final_ecs[rep(seq_len(nrow(final_ecs)), reps), 1:4]
  df2$Matched.Compound <- unlist(mapply(function(x, y) c(x, rep(NA, y)),
                                        cpd_split, reps - lengths(cpd_split)))

  matched_res <- merge(matched_res, df2)
  matched_res <- matched_res[, -6]                     # drop RT.Rank
  matched_res[, 6] <- as.character(matched_res[, 6])

  # 4. per compound, join ECs that share a member ------------------------------
  my_final_cpds <- stats::aggregate(. ~ Matched.Compound, matched_res, paste,
                                    collapse = "___")
  my_final_cpds_list <- lapply(split(my_final_cpds$Empirical.Compound,
                                     my_final_cpds$Matched.Compound), unlist)
  cpd2ec1 <- lapply(seq_along(my_final_cpds_list), function(x) {
    ecs <- unlist(strsplit(my_final_cpds_list[[x]], "___", fixed = TRUE))
    if (length(ecs) > 1) {
      ecs_list <- as.list(strsplit(ecs, ";", fixed = TRUE))
      m <- sapply(ecs_list, function(a) sapply(ecs_list, function(b)
        length(intersect(a, b)) > 0))
      g <- igraph::groups(igraph::components(
        igraph::graph_from_adjacency_matrix(m)))
      ecs <- paste0(sapply(g, function(z) paste0(ecs[z], collapse = "|")),
                    collapse = "___")
    }
    ecs
  })
  names(cpd2ec1) <- names(my_final_cpds_list)

  update_ecs <- lapply(seq_along(cpd2ec1), function(z) {
    ecs_old <- unlist(strsplit(my_final_cpds_list[[z]], "___", fixed = TRUE))
    ecs_new <- unlist(strsplit(cpd2ec1[[z]], "___", fixed = TRUE))
    for (i in seq_along(ecs_new)) {
      pattern_vec <- unlist(strsplit(ecs_new[i], "\\|"))
      ecs_old[ecs_old %in% pattern_vec] <- paste0(unique(pattern_vec),
                                                  collapse = "|")
    }
    paste0(ecs_old, collapse = "___")
  })
  my_final_cpds$Empirical.Compound <- unlist(update_ecs)

  # data.table's by = Matched.Compound on one-row-per-compound input, in base R.
  expand <- function(r) {
    sp <- function(col) unlist(strsplit(as.character(my_final_cpds[[col]][r]),
                                        "___", fixed = TRUE))
    data.frame(Query.Mass = sp("Query.Mass"),
               Matched.Compound = my_final_cpds$Matched.Compound[r],
               Matched.Form = sp("Matched.Form"),
               Retention.Time = sp("Retention.Time"),
               Mass.Diff = sp("Mass.Diff"),
               Empirical.Compound = sp("Empirical.Compound"),
               stringsAsFactors = FALSE)
  }
  matched_res <- do.call(rbind, lapply(seq_len(nrow(my_final_cpds)), expand))
  rownames(matched_res) <- NULL

  # 5. EC names ----------------------------------------------------------------
  ec <- matched_res$Empirical.Compound
  ec_unique <- unique(matched_res$Empirical.Compound)
  for (i in seq_along(ec_unique)) {
    ec <- replace(ec, grep(paste0("\\b", ec_unique[i], "\\b"), ec, perl = TRUE),
                  paste0("EC000", i))
  }
  matched_res$Empirical.Compound <- gsub("\\|.*", "", ec)

  # dictionaries ---------------------------------------------------------------
  matched_ts <- prep$expr_dic[matched_res[, 1]]
  ec_exp_dict <- .mmc_ma_dict(data.frame(key = matched_res[, 6],
                                         value = as.numeric(matched_ts),
                                         stringsAsFactors = FALSE))
  mz2ec_dict <- .mmc_ma_dict(matched_res[, c(1, 6)])
  refmz_num <- suppressWarnings(as.numeric(names(mz2ec_dict)))
  input_num <- suppressWarnings(as.numeric(as.character(input_mzlist)))
  input_num <- unique(input_num[!is.na(input_num)])
  hits_index <- which(!is.na(refmz_num) & refmz_num %in% input_num)

  list(
    matched_res = matched_res,
    initial_ecs = df_ecs,
    ec_exp_dict = ec_exp_dict,
    ec_exp = unlist(lapply(ec_exp_dict, max)),
    cpd_exp_dict = .mmc_ma_dict(data.frame(key = matched_res[, 2],
                                           value = as.numeric(matched_ts))),
    cpd_ecpd_dict = .mmc_ma_dict(matched_res[, c(2, 6)]),
    ecpd_cpd_dict = .mmc_ma_dict(matched_res[, c(6, 2)]),
    mz2cpd_dict = .mmc_ma_dict(matched_res[, c(1, 2)]),
    mz2ec_dict = mz2ec_dict,
    ec2mz_dict = .mmc_ma_dict(matched_res[, c(6, 1)]),
    cpd_form_dict = .mmc_ma_dict(cbind(matched_res[, 2], matched_res[, 3])),
    cpd_ecpd_counts = cpd2ec1,
    input_ecpdlist = unique(unlist(mz2ec_dict[hits_index])),
    total_matched_ecpds = unique(as.vector(matched_res$Empirical.Compound))
  )
}

#' Pathway sets in empirical-compound space (MetaboAnalystR's library setup)
#'
#' Port of the pathway steps of MetaboAnalystR's `.setup.psea.library()`:
#' pathways with fewer than `min_lib` compounds are dropped, each remaining
#' pathway is mapped to the ECs of its compounds, and pathways with no EC are
#' dropped.
#'
#' @param lib Library from `mmc_ma_library_from_model()`.
#' @param ecs Output of `mmc_ma_build_ecs()`.
#' @param min_lib Smallest pathway kept, in library compounds (upstream
#'   `PerformPSEA(minLib = 3)`).
#' @return List with `cpds`, `name`, `emp_cpds` (upstream's three fields, same
#'   order) and `id` (ours: the model's pathway ids, which upstream drops here).
#' @examples
#' \dontrun{
#' pw <- mmc_ma_pathway_ecs(lib, ecs)
#' }
mmc_ma_pathway_ecs <- function(lib, ecs, min_lib = 3) {
  min_inx <- which(sapply(lib$pathways$cpds, length) >= min_lib)
  cpds <- lib$pathways$cpds[min_inx]
  ids  <- lib$pathways$id[min_inx]
  nms  <- lib$pathways$name[min_inx]
  d <- ecs$cpd_ecpd_dict
  emp <- lapply(cpds, function(x) {
    unique(unlist(d[stats::na.omit(match(x, names(d)))]))
  })
  null_inx <- sapply(emp, is.null)
  list(cpds = cpds[!null_inx], name = nms[!null_inx],
       emp_cpds = emp[!null_inx], id = ids[!null_inx])
}
