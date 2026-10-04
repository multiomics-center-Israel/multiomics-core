# Contract: RNA-seq `final_results` columns and run artifacts

> Status: **IMPLEMENTED** — describes what the pipeline writes today
> Last updated: 2026-08-18
> Companion to `CONTRACT_metabolomics_de_and_final_results.md`, which covers a
> proposed renaming and is design-only. This document is descriptive, not a plan.

Written because a reader computed a fold change from the per-sample columns of
`Final_results_ALL_*.xlsx`, got `-0.35` in log2, and could not reconcile it with
the reported `linearFC` of `-1.61`. Nothing in the workbook explained the gap.
Both numbers were correct; they are different quantities.

---

## 1. Column blocks, in the order they appear

`Datasets/final_results.tsv` and the `Results` sheet of both Excel workbooks
carry the same blocks.

| Block | Columns | Scale | Written by |
|---|---|---|---|
| identifier | `Gene` | — | `build_final_results_rnaseq()` |
| measured | `<sample>` | raw counts, post-filter | `pre$expr_filt` |
| model input | `<sample>.norm` | DESeq2 normalised counts | `deseq2_normalized_counts()` |
| group summary | `Mean.<group>` | normalised counts | `compute_group_mean_columns()` |
| group summary | `CV.<group>` | percent, on linear CPM | `compute_group_cv_columns()` |
| statistics | `log2FC.<contrast>` | log2 | `build_rnaseq_summary_df()` |
| statistics | `log2FC_from_means.<contrast>` | log2 | `compute_naive_log2fc_columns()` |
| statistics | `linearFC.<contrast>` | signed linear | `build_rnaseq_summary_df()` |
| statistics | `pvalue.` / `padj.<contrast>` | — | DESeq2 |
| statistics | `upDown.` / `manual_cutoffs.<contrast>` | — | Excel layer |
| flag | `pass_any_contrast` | `1` or `NA` | `add_pass_any_contrast()` |

Proteomics uses the same layout with an `.imputs.` infix on the statistics and
the imputed log2 matrix in the `.norm` block. Metabolomics and lipidomics emit
neither the `.norm` block nor `log2FC`; the new blocks are opt-in parameters of
`build_final_results_generic()`.

---

## 2. The three fold-change columns

They answer different questions and are not interchangeable.

| Column | Definition | Exact? |
|---|---|---|
| `log2FC` | the DESeq2 negative-binomial GLM coefficient | the model's own estimate |
| `linearFC` | `2^log2FC` if `log2FC >= 0`, else `-1 / 2^log2FC` | exactly derived from `log2FC` |
| `log2FC_from_means` | `log2(Mean.<num> / Mean.<den>)` | exactly derived from the two `Mean` cells |

Three consequences worth knowing:

1. **`linearFC` is not a log.** `-1.61` means 1.61-fold lower, not `log2FC = -1.61`.
   It is a signed reciprocal below 1, so `abs()` is a magnitude and the sign is
   the direction.
2. **`log2FC` is stored unrounded on purpose.** Rounding it made
   `signif(2^log2FC, 3)` disagree with the reported `linearFC` for about 1% of
   features. See `test-fc-provenance.R`, section P4.
3. **The raw `<sample>` columns will not reproduce `log2FC`.** They are not
   corrected for library size, so any difference in sequencing depth between the
   groups shifts the ratio. On the reference run the two groups differed by
   about 25% in depth, worth 0.41 in log2.

`log2FC` and `log2FC_from_means` normally agree closely: median difference
0.0008 on the reference run, 95th percentile 0.016. They diverge for low-count
features, because one is a GLM coefficient and the other a ratio of arithmetic
means.

---

## 3. Unbounded fold changes

A feature with no detected counts in one group has a `Mean` of exactly zero, so
its ratio has no finite value. The pipeline still reports a number: on the
reference run one gene shows `linearFC = 2.97e+07`. That is not a magnitude, it
records that nothing was detected on the other side.

20 of the 61 differentially expressed genes on that run were of this kind. When
quoting fold changes, split them out. `log2FC_from_means` is `NA` for these
rather than `-Inf`.

---

## 4. Detecting a collapsed fold change

`check_log2fc_shrinkage()` runs at export and writes
`Datasets/log2fc_shrinkage_check.tsv` on every run. It compares the model
estimate against the model-free one, over features where
`|log2FC_from_means| >= 0.5`.

| median `\|log2FC\| / \|log2FC_from_means\|` | reading |
|---|---|
| about 1.00 | no shrinkage applied |
| about 0.85 | ordinary `betaPrior = TRUE` behaviour |
| below 0.50 | flagged `shrunk` |
| near 0 | flagged `collapsed` |

The 0.5 floor is load-bearing. For a null feature both estimates are near zero
and their ratio is noise: gene `evm.TU.ptg000675l_np1212.2` on the reference run
shows a ratio of 0.22 arising from an absolute difference of 0.004 at
`padj = 0.9996`. Judge agreement by the difference, never by the ratio alone.

The check also counts significant features whose `log2FC` is flat. Those are the
points that draw a vertical stripe at x = 0 in a volcano, which otherwise reads
as "nothing changed" rather than "the estimates were flattened".

---

## 5. Other artifacts a run writes

| File | Contents |
|---|---|
| `results_fact_sheet.tsv` | every headline number with the file it can be checked against |
| `Datasets/log2fc_shrinkage_check.tsv` | the shrinkage diagnostic, per contrast |
| `Datasets/de_summary_counts.tsv` | up / down / total per contrast |
| `de_summary_counts.tsv` (mode root) | the same, **different header** — see below |
| `Final_results_ALL_*.xlsx` | all features, plus `Cutoffs` and `How to read` sheets |
| `Final_results_DE_*.xlsx` | passing features, ordered by the clustering |

### Known inconsistency

`de_summary_counts.tsv` is written by two code paths under the same name with
different headers:

| Location | Header | All-contrasts row |
|---|---|---|
| mode root | `contrast`, `up`, `down`, `total` | `any` |
| `Datasets/` | `Name`, `up`, `down`, `any` | `pass_any` |

Any consumer must accept both. `.read_de_summary_counts()` in
`R/domain/rnaseq/12_fact_sheet.R` normalises them. Unifying the two writers is
worth doing and has not been done.

---

## 6. How to check a row by hand

```r
fr <- read.delim("outputs/<run>/rna/Datasets/final_results.tsv", check.names = FALSE)

# the model-free estimate is reproducible from the two Mean cells
stopifnot(all.equal(fr$log2FC_from_means.<contrast>,
                    log2(fr$Mean.<num> / fr$Mean.<den>)))

# and linearFC is reproducible from log2FC, exactly
l <- fr$log2FC.<contrast>
stopifnot(identical(fr$linearFC.<contrast>,
                    signif(ifelse(l >= 0, 2^l, -1 * (2^-l)), 3)))
```

The `How to read` sheet in each workbook states the same rules for readers who
never open R.
