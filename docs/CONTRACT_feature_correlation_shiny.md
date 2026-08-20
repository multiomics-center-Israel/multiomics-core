# Feature correlation — Shiny handoff

**Audience:** whoever owns the Shiny app.
**Status:** engine implemented in `multiomics-core`; the GUI panel is not built.
**Payload impact:** **none.** No new keys, no contract change, nothing to
re-run. Everything below works against payloads the pipeline already writes.

---

## What this is

A researcher picks one metabolite in the GUI and wants to see which other
metabolites show the same expression pattern — features that rise and fall with
it across the samples.

The pipeline cannot precompute this: the feature is chosen interactively, so it
is unknown at `tar_make()` time. It also does not need to. One-vs-all is
`O(features × samples)` — **measured at ~0.26 s for 3,000 features × 24 samples**,
~0.36 s with 15% missing values. Compute it on demand in the reactive.

`multiomics-core` provides three functions. The GUI supplies a dropdown and a
table.

---

## The three functions

All live in `R/domain/metabolomics/09_feature_correlation.R` and
`R/domain/metabolomics/09b_feature_correlation_plots.R`.

```r
prepare_correlation_matrix(expr_mat, sample_meta, condition_col,
                           sample_col = "sample_id", qc_flag_column = NULL)

correlate_feature_vs_all(expr_mat, feature_id, min_n = 5, top_n = NULL)

plot_feature_correlation_profiles(expr_mat, meta, feature_id, partner_ids,
                                  group_col, sample_col = "sample_id",
                                  label_map = NULL)
```

### Worked example

```r
# Once per payload:
mat <- prepare_correlation_matrix(
  payload$expr_norm,
  payload$sample_meta,
  condition_col = payload$group
)

# Dropdown choices — see "Building the dropdown" below.
choices <- rownames(mat)

# Per selection:
res <- correlate_feature_vs_all(mat, feature_id = input$feature)

DT::datatable(head(res, 50))

plot_feature_correlation_profiles(
  mat, payload$sample_meta,
  feature_id  = input$feature,
  partner_ids = head(res$feature_id, 5),
  group_col   = payload$group
)
```

---

## ⚠️ Filter first — `expr_norm` still contains QC samples

`payload$expr_norm` is assigned straight from `pre$expr_work`, and **QC, blank
and pooled samples are still in it.** This pipeline filters at each point of
use, not upstream — `mod_metabolomics_clustering()` and the DE code both call
`filter_to_biological()` themselves.

Pooled QC samples sit at the average of everything, so leaving them in distorts
every coefficient in the table. **Always call `prepare_correlation_matrix()`
first.** It wraps the same `filter_to_biological()` the rest of the pipeline
uses, so the GUI does not need to know the QC/blank/pool naming conventions.

**One gap to be aware of:** the payload does not carry the project's
`qc.qc_flag_column` setting. Projects that mark technical samples through that
configured column — rather than through the condition value or the sample ID —
must pass `qc_flag_column = "<column>"` explicitly, or those samples will not be
recognised as technical. Surfacing that setting in the payload would be a
contract change and needs team sign-off; it was deliberately left out of scope.

---

## Result columns

| Column | Meaning |
|---|---|
| `feature_id` | The partner feature |
| `pearson_r`, `pearson_pvalue`, `pearson_padj` | **Primary** measure |
| `spearman_rho`, `spearman_pvalue`, `spearman_padj` | **Supporting** measure |
| `n_used` | Samples this pair actually shared |
| `direction` | `"positive"`, `"negative"`, `"none"` (exactly zero), or `NA` |

Rows are returned sorted by `abs(pearson_r)` descending, ties broken by
`feature_id`, with untested pairs last. The ordering is fully deterministic —
safe to display as-is.

### Why two coefficients

**Pearson is primary.** The matrix is already log2-transformed, and
`build_clustering_distance()` also uses Pearson — so this ranking and the
pipeline's own clustering agree on what "similar" means. Rank by
`abs(pearson_r)`.

**Spearman is a robustness check.** Read the pair together:

- **Both high** — a convincing hit.
- **Spearman high, Pearson clearly lower** — monotonic but non-linear, or a
  saturating response. Real, but not a straight line.
- **Pearson high, Spearman clearly lower** — usually one outlying sample is
  carrying the correlation. Worth a look at the scatter before believing it.

Consider showing both columns rather than hiding Spearman behind a toggle; the
disagreement is often the most informative thing in the row.

---

## Missing values and untested pairs

Metabolomics has **no imputation stage**, so real `NA`s reach the matrix. Each
pair is computed on its own mutually-observed samples, and `n_used` reports that
count — different rows legitimately rest on different numbers of samples.

Rows can come back with `NA` statistics for three reasons. They are **kept in
the table on purpose**, so a reader sees "not tested" rather than wondering
where a feature went:

1. **Fewer than `min_n` shared samples** (default 5). A correlation of 0.99 from
   4 points should never sit at the top of the table.
2. **Either feature is constant on the samples that pair shares** — which can
   happen even when both vary across the full matrix. The correlation is
   genuinely undefined there.
3. **A globally constant feature.**

`NA` in `pearson_padj` means *not tested*, never *not significant*. Consider
rendering those cells as "—" rather than blank. BH adjustment is computed across
testable pairs only, and the query feature is excluded before adjusting, so
q-values are not inflated by rows that were never tested.

### On the Spearman p-value

`spearman_pvalue` is an **asymptotic approximation**, not an exact permutation
test: the Pearson t-transform applied to the rank correlation, matching
`stats::cor.test(method = "spearman", exact = FALSE)`. It degrades at small n
and in the presence of ties. `pearson_pvalue` is the exact t-test under
bivariate normality. If the UI shows a tooltip on the Spearman column, say so
there.

---

## Building the dropdown

Use **`rownames(expr_norm)` as the source of truth** for which features are
selectable, and `feature_annot` only to look up display labels. Driving the list
from `feature_annot` instead risks silently hiding features that exist in the
matrix but are missing an annotation row.

For display names, follow the precedence the rest of the codebase uses
(`R/domain/multiomics/02_mae.R`):

```
gene_symbol → original_id → Name → HMDB_ID / HMDB → feature_id
```

Pass that mapping to the plot as `label_map` (a named character vector, IDs as
names). It may be partial — anything without an entry falls back to its raw ID.
An **unnamed** vector is ignored rather than trusted, since indexing by position
would mislabel every feature.

---

## Adding a "same cluster?" column

`payload$clust_partition` gives cluster membership per feature. Joining it onto
the result table gives independent corroboration: a partner with a high `r`
*that also landed in the same cluster* is a much stronger claim than either
signal alone.

**This column is three-state and must stay that way:**

| Value | Meaning |
|---|---|
| `TRUE` | Both features clustered, same cluster |
| `FALSE` | Both features clustered, different clusters |
| `NA` | One or both were **not part of clustering** |

Clustering runs on **DE features only**, so `NA` is common and expected. It
means *unknown* — never *not in the same cluster*. Collapsing `NA` to `FALSE`
would tell the user that a non-DE metabolite is dissimilar when nothing of the
sort was ever tested.

---

## Caveats worth putting in the UI

- **Small n makes `r` unstable.** With 8 samples, `|r| > 0.7` is around the
  p ≈ 0.05 mark. Show `n_used` next to the coefficient, not in a hidden column.
- **Normalization is compositional.** TSS and PQN can induce mild spurious
  *negative* correlation between features. Treat weak negative hits with more
  suspicion than weak positive ones.
- **A strong `r` may just mean "both respond to the treatment."** Correlations
  are computed across all biological samples, so two metabolites that both go up
  under treatment will correlate strongly even if biologically unrelated.

  This is **intended** — it is exactly the "same expression pattern" the feature
  was requested for. But it is co-*response*, not necessarily co-*regulation*,
  and the UI should say so plainly. Something like:

  > Correlations are computed across all samples, so features that respond
  > similarly to the treatment will rank highly. This shows shared behaviour,
  > not necessarily a direct biological link.

  If a user ever needs to distinguish co-regulation from shared treatment
  response, the follow-up is a within-condition correlation. Not implemented —
  ask before building it.

---

## Tests

`tests/testthat/test-metab-feature-correlation.R` — 100 assertions covering the
maths (agreement with `stats::cor.test` for both coefficients, including a
tied-value case), the pairwise-complete zero-variance guard in both directions,
BH scope, `min_n`, determinism, and the plot's edge cases. Run it before
changing anything here.
