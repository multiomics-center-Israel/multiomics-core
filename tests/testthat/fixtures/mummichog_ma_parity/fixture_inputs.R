# tests/testthat/fixtures/mummichog_ma_parity/fixture_inputs.R
#
# The SYNTHETIC inputs shared by the MetaboAnalyst-compatible mummichog parity
# generators (components 2 onwards). Sourced by the generators only; each
# generator stores the inputs it used in its .rds, and the tests read them from
# there, so this file is never needed at test time.
#
# Nothing here is real data. Compound ids C91xxx are invented; C00002 (ATP) and
# C00008 (ADP) are real KEGG ids used only because they are on MetaboAnalyst's
# currency list, which the matcher must drop.
#
# Designed cases (feature ids name the case):
#   both_*    compound seen as [M+H]+ AND [M-H]- at nearly the same RT
#   rescue_*  positive [M+K]+ (not a primary ion) + negative [M-H]- at close RT
#   kplus_*   positive [M+K]+ only (no primary ion anywhere)
#   farrt_*   one compound as [M+H]+ at two distant retention times
#   samemz_*  a positive and a negative feature with the identical m/z value
#   ppm_*     one feature 8 ppm off [M+H]+ (matches at 10 ppm), one 15 ppm off
#   print_*   two distinct m/z doubles that print alike at 15 digits; upstream
#             rounds m/z to that text form first, so print_b is dropped as a
#             duplicate of print_a
#   tie_*     two features sharing a retention time (random tie-breaking)
#   range_*   m/z outside 50 < m/z < 2001
#   atp_*     ATP [M-H]- (currency)
#   noise_*   random m/z, mostly unmatched

ma_fixture_inputs <- function() {
  proton <- 1.00727646677
  # 30 invented compounds with masses spread over the 100-900 range.
  mw <- round(seq(121.0507, 887.3311, length.out = 30) + (1:30) * 0.00731, 6)
  compounds <- data.frame(
    id   = c(sprintf("C91%03d", 1:30), "C00002", "C00008"),
    name = c(sprintf("synthetic compound %02d", 1:30), "ATP", "ADP"),
    mw   = c(mw, 506.995745, 427.029414),
    stringsAsFactors = FALSE
  )
  pathways <- list(
    syn10001 = sprintf("C91%03d", 1:8),
    syn10002 = sprintf("C91%03d", c(5:14)),
    syn10003 = sprintf("C91%03d", c(15:22)),
    syn10004 = c(sprintf("C91%03d", c(20:30)), "C00002", "C00008"),
    syn10005 = sprintf("C91%03d", c(2, 9))          # 2 compounds: < minLib
  )
  pathway_names <- sprintf("Synthetic pathway %d", seq_along(pathways))

  m <- stats::setNames(compounds$mw, compounds$id)
  mh  <- function(id) m[[id]] + proton                 # [M+H]+
  mmh <- function(id) m[[id]] - proton                 # [M-H]-
  mk  <- function(id) m[[id]] + 37.9555 + proton       # [M+K]+

  f <- list()
  add <- function(id, mz, rt, pol) {
    f[[length(f) + 1L]] <<- data.frame(feature_id = id, mz = mz, rt = rt,
                                       polarity = pol, stringsAsFactors = FALSE)
  }
  for (k in 1:10) {
    cid <- sprintf("C91%03d", k)
    add(sprintf("both_%02d_pos", k), mh(cid),  1 + k, "positive")
    add(sprintf("both_%02d_neg", k), mmh(cid), 1 + k + 0.01, "negative")
  }
  for (k in 11:14) {
    cid <- sprintf("C91%03d", k)
    add(sprintf("rescue_%02d_pos", k), mk(cid),  12 + k, "positive")
    add(sprintf("rescue_%02d_neg", k), mmh(cid), 12 + k + 0.02, "negative")
  }
  for (k in 15:17) {
    add(sprintf("kplus_%02d", k), mk(sprintf("C91%03d", k)), 12 + k, "positive")
  }
  add("farrt_18_a", mh("C91018"), 2.5, "positive")
  add("farrt_18_b", mh("C91018") + 1e-6, 28.0, "positive")
  add("samemz_19_pos", 400.123456, 9.0, "positive")
  add("samemz_19_neg", 400.123456, 9.4, "negative")
  add("ppm_20_in",  mh("C91020") * (1 + 8e-6),  14.0, "positive")
  add("ppm_21_out", mh("C91021") * (1 + 15e-6), 14.5, "positive")
  add("print_a", 300.1234567890121, 6.0, "negative")
  add("print_b", 300.1234567890123, 6.5, "negative")
  add("tie_22", mh("C91022"), 17.25, "positive")
  add("tie_23", mmh("C91023"), 17.25, "negative")
  add("range_low", 45.0, 3.0, "positive")
  add("range_high", 2005.0, 3.0, "negative")
  add("atp_neg", mmh("C00002"), 7.0, "negative")

  designed <- do.call(rbind, f)
  # 90 noise features, seeded so the fixture is fixed.
  noise <- withr::with_seed(20261102, data.frame(
    feature_id = sprintf("noise_%03d", 1:90),
    mz = round(stats::runif(90, 100, 1000), 5),
    rt = round(stats::runif(90, 0.5, 30), 3),
    polarity = sample(c("positive", "negative"), 90, replace = TRUE),
    stringsAsFactors = FALSE))
  peaks <- rbind(designed, noise)

  # p-values: designed signal features significant, a few noise features too;
  # two exact p ties. Statistics signed, never zero.
  n <- nrow(peaks)
  sig_ids <- grepl("^(both|rescue|kplus|farrt|tie|ppm|atp)", peaks$feature_id)
  p <- withr::with_seed(7, ifelse(sig_ids, stats::runif(n, 1e-4, 0.04),
                                  stats::runif(n, 0.06, 1)))
  p[peaks$feature_id %in% c("noise_001", "noise_002", "noise_003")] <- 0.01
  p[peaks$feature_id == "samemz_19_pos"] <- 0.03
  p[peaks$feature_id == "samemz_19_neg"] <- 0.02    # neg survives the dedup
  p[peaks$feature_id %in% c("both_01_pos", "both_02_neg")] <- 0.005  # p tie
  peaks$p_value <- p
  peaks$statistic <- withr::with_seed(11, round(stats::rnorm(n, 0, 2), 4))
  peaks$statistic[peaks$statistic == 0] <- 0.1

  list(compounds = compounds, pathways = pathways,
       pathway_names = pathway_names, peaks = peaks,
       ppm = 10, cutoff = 0.05, rt_frac = 0.02)
}
