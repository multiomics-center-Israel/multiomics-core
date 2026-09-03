# Feature correlation — Shiny handoff

**Audience:** whoever owns the Shiny app.
**Status:** engine implemented in `multiomics-core`; the GUI panel is not built.
**Payload impact:** **this feature requires payload schema 2.1.** It adds one
canonical key, `expr_norm_missing`. Runs exported at 2.0 still load, but their
correlation statistics are not trustworthy — see *Schema 2.1* below.

---

## What this is

A researcher picks one metabolite in the GUI and wants to see which other
metabolites show the same expression pattern — features that rise and fall with
it across the samples.

The pipeline cannot precompute this: the feature is chosen interactively, so it
is unknown at `tar_make()` time. It also does not need to. One-vs-all is
`O(features × samples)` — measured at ~0.26 s for 3,000 features × 24 samples,
~0.36 s with 15% missing. Compute it on demand in the reactive.

`multiomics-core` provides four functions. The GUI supplies a dropdown and a
table.

---

## The four functions

In `R/domain/metabolomics/09_feature_correlation.R` and
`09b_feature_correlation_plots.R`.

```r
restore_missing_values(expr_mat, missing_mask)

prepare_correlation_matrix(expr_mat, sample_meta, condition_col,
                           sample_col = NULL, qc_flag_column = NULL)

correlate_feature_vs_all(expr_mat, feature_id, min_n = 5, top_n = NULL)

plot_feature_correlation_profiles(expr_mat, meta, feature_id, partner_ids,
                                  group_col, sample_col = NULL,
                                  label_map = NULL)
```

### Worked example

```r
# Once per payload — restore, then filter. Order matters.
mat <- restore_missing_values(payload$expr_norm, payload$expr_norm_missing)
mat <- prepare_correlation_matrix(mat, payload$sample_meta,
                                  condition_col = payload$group)

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

Note `sample_col` is omitted in both calls. **Do not pass `"sample_id"`** — see
*Sample IDs* below.

---

## ⚠️ Step 1: restore missingness

`payload$expr_norm` is declared NA-free, and the metabolomics builder enforces
that by substituting values. Precisely what happens:

- **`NA` and `NaN` cells are replaced with the feature's row median.** They are
  then indistinguishable from real measurements.
- **`±Inf` cells are *not* replaced** — the fill only targets `NA`/`NaN`, so an
  `-Inf` (which `transform_metab()` can produce by log-transforming a
  non-positive value) remains in `expr_norm` as `-Inf`.

Both cases are recorded in `expr_norm_missing`, and
`restore_missing_values()` turns both back into `NA`. So after restoration the
matrix contains only finite values and `NA`, and every consumer — Pearson,
Spearman, and the profile plot — shares one definition of an observed value.

**Skipping this step does not error; it silently changes what the numbers mean.**
Correlations get computed on substituted values at the wrong effective n, which
is invalid as pairwise observed-data inference. It can attenuate `r` (filled
values pull toward the row centre) while simultaneously inflating the degrees of
freedom, so the error does not even bias p-values in a predictable direction.

---

## ⚠️ Step 2: filter QC samples

`expr_norm` still contains QC, blank and pooled samples. This pipeline filters
at each point of use, not upstream — `mod_metabolomics_clustering()` and the DE
code both call `filter_to_biological()` themselves. Pooled QC samples sit at the
average of everything, so leaving them in distorts every coefficient.

`prepare_correlation_matrix()` wraps the same `filter_to_biological()` the rest
of the pipeline uses, so the GUI does not need to know the QC/blank/pool naming
conventions.

**One gap:** the payload does not carry the project's `qc.qc_flag_column`
setting. Projects that mark technical samples through that configured column —
rather than through the condition value or the sample ID — must pass
`qc_flag_column = "<column>"` explicitly. Surfacing that setting in the payload
would be a further contract change and was left out of scope.

---

## Schema 2.1 and the `expr_norm_missing` key

`expr_norm_missing` is a logical matrix with the same dim and dimnames as
`expr_norm`, `TRUE` where the cell was not a usable observation before the fill.
**It has three meaningful states, and they are deliberately distinguishable:**

| Value | Meaning | What the GUI should do |
|---|---|---|
| matrix, any `TRUE` | missingness recorded and restorable | normal path; results are pairwise-complete |
| matrix, all `FALSE` | schema ≥ 2.1 verified the data complete | normal path; nothing to restore |
| `NULL` / absent | **pre-2.1 payload — provenance unknown** | warn the user; recommend re-export |

`NULL` means something different per omics: for **metabolomics** it means a
pre-2.1 payload whose missingness was erased; for **rnaseq and proteomics** it
means *not applicable*, since those builders perform no fill.

Metabolomics payloads at schema ≥ 2.1 are **required** to carry the key —
`assert_shiny_payload_contract()` fails them otherwise, so the version bump
enforces the guarantee rather than merely advertising it.

### Checking provenance

```r
payload$payload_version                 # "2.1" or later
m <- payload$expr_norm_missing
if (is.null(m)) {
  # legacy payload: correlations will run, but n_used is not pairwise-complete
} else {
  sprintf("%d of %d cells unobserved", sum(m), length(m))
}
```

**Do not infer provenance from `n_used`.** It is tempting to treat "every pair
reports the same `n_used`" as evidence the mask was missing — it is not. An
all-`FALSE` mask is a legitimate verified-complete payload, and constant
`n_used` is exactly the correct result for one. The mask and schema version are
the authoritative signal; `n_used` is a statistic, not a diagnostic.

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

Metabolomics has **no imputation stage**, so real missingness reaches the matrix
once restored. Each pair is computed on its own mutually-observed samples, and
`n_used` reports that count — different rows legitimately rest on different
numbers of samples.

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

## Sample IDs

**Both** `prepare_correlation_matrix()` and `plot_feature_correlation_profiles()`
default `sample_col` to `NULL`, meaning *take the IDs from
`rownames(sample_meta)`* — which is what the payload builder populates, from the
project's configured `effects$samples` column.

**Do not pass `"sample_id"`.** The shipped metabolomics and proteomics templates
both set `effects$samples: "SampleID"`, and the payload does not expose that
name. Omitting the argument is correct for any standard payload. Pass it
explicitly only when feeding a raw `pre$meta` whose rownames are not the IDs.

Resolved IDs are validated the same way whichever route they came from: no
missing or empty values, no duplicates, and full coverage of
`colnames(expr_mat)`. A mismatch errors with a message naming the fix rather
than mis-aligning samples silently.

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

## The profile plot

Overlays z-scored group-mean profiles: the chosen feature in black, its
correlates in colour. Z-scoring per feature is what lets a low-abundance and a
high-abundance metabolite be compared — the shape is the point, not the height.

**Gaps are real.** If a feature has no observed values in an entire group, that
profile point is drawn as a **break in the line**, never interpolated — a
fabricated midpoint would show a measurement nobody took. The function emits a
message naming how many feature × group cells are gaps, so a broken line is
explainable rather than mysterious.

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

- `tests/testthat/test-metab-feature-correlation.R` — the engine: agreement with
  `stats::cor.test` for both coefficients (including a tied-value case), the
  pairwise-complete zero-variance guard in both directions, non-finite handling
  across Pearson *and* Spearman, BH scope, `min_n`, determinism, sample-ID
  resolution, and the plot's edge cases including group gaps.
- `tests/testthat/test-metab-expr-norm-missing-mask.R` — the payload key: the
  three mask states, `-Inf`/`NaN` capture, contract validation, the metabolomics
  ≥ 2.1 requirement, and v2.0 backward compatibility.

Run both before changing anything here.
