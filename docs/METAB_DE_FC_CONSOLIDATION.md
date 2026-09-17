# Three fixes to the metabolomics / lipidomics fold change

One branch, three commits, in dependency order. They are one change: commit 1
reads the label commit 3 corrects, so shipping 1 without 3 multiplies every
metabolomics log2FC on a decadic-transform lane by `log2(10)` = 3.32.

| # | Commit | What it fixes | What moves numerically |
|---|--------|---------------|------------------------|
| 1 | `fix(metab,lipid): the fold change comes from the matrix that was tested` | DE ran the test on the normalised matrix but took the fold change from the **raw filtered** matrix (`mat_for_fc = pre$expr_filt`), so the reported direction could contradict its own p-value. Both domains now take the difference on `mat_for_test`, and `fc_to_log2_units()` puts it in log2 units. | `logFC` and `AveExpr` on **every** metabolomics and lipidomics lane. Sign flips wherever the sample normalisation mattered. |
| 2 | `fix(metab): DE tests the matrix that carries the sample normalisation` | On a variance-scaled metabolomics lane, DE reached for `pre$expr_log` — which in metabolomics is `transform_metab()` of the filtered matrix, the **transform alone**. PQN / TSS / the median shift run on parallel targets, so those p-values were computed on data no sample normalisation had ever touched. `mod_met_corrected()` now also emits `mat_pre_scale`, carried as `expr_pre_scale`. | p-values (and everything downstream) on metabolomics lanes with `scaling` in `auto` / `pareto` / `range`. Lanes with `scaling: none` are untouched — that branch already used `expr_work`. |
| 3 | `fix(metab): the recorded transform names the one that was APPLIED` | `info$normalization$transform` reported the **configured** transform. TSS / PQN / both EigenMS targets normalise the linear matrix and then log2-transform it themselves, so a lane configured "PQN + glog10" was a log2 matrix wearing a glog10 label. Commit 1's unit conversion believes that label. Each normalisation target now declares the transform it applied, and `mod_met_corrected()` reads it off the chosen target. | `logFC` and `AveExpr` on metabolomics lanes that are TSS/PQN/EigenMS **and** configured `log10`/`glog10`: they were 3.32x too large. `transform_configured` keeps the configured value. |

## Rebased onto `main` (2026-09-17)

The branch was written on `19d985c` and replayed onto `main` 67 commits later.
Two files conflicted, both resolved to keep upstream's behaviour *and* the fix:

- **`R/domain/metabolomics/03_differential.R`** — upstream had already moved
  limma to `mat_for_fc = NULL` and had added, around `mat_raw`, a re-alignment
  to the biological samples the fit keeps (`0f6f4fe`) and a scale caveat for
  `chosen_norm = "none"`. Commit 1 deletes `mat_raw` outright, so the
  re-alignment goes with it — nothing but `mat_for_test` is read any more. The
  caveat is **kept and narrowed**: taking the fold change from the tested matrix
  in its declared unit answers it for every transform except `transform: "none"`,
  where the recorded value means "no transform applied" while
  `fc_to_log2_units()` reads it as "the values are linear". It now fires only on
  `chosen_norm = "none"` **and** `transform = "none"` — and for limma too, whose
  coefficient goes through the same conversion.
- **`R/modules/metabolomics/00_mod_preprocessing.R`** — upstream added two more
  normalisations to the `chosen_norm` switch: `none` (`36c779a`) and
  `bio_factor` (`3e79cb9`). The switch now selects the whole target, not its
  matrix, for all seven branches. `bio_factor` normalises on the linear matrix
  and log2s it itself, like tss/pqn, so it declares `transform = "log2"`;
  `none` reads `mod_met_log()`, the one target whose applied transform *is* the
  configured one, so it declares the configured value and means it. Missing
  either would have mislabelled that lane and multiplied its fold changes by
  `log2(10)` on a decadic config — the exact defect commit 3 fixes.

`tests/testthat/test-metab-transform-recorded-as-applied.R` now enumerates the
`chosen_norm` branches straight out of the parsed body of `mod_met_corrected()`
and fails if one of them is not covered here, so the *next* normalisation cannot
be added silently either.

Commit 2 also had to compose with `0f6f4fe`, which excludes QC/blank samples
from the fit. `expr_pre_scale` goes through the same post-normalisation pool
filter and the same drift correction as `expr_work`, so it stays column-aligned
with `pre$meta`: the in-fit exclusion lands once, on columns whose condition is
known, and nothing is filtered twice. Section 3 of
`tests/testthat/test-metab-de-tests-normalised-matrix.R` pins both halves.

## How to check it — a slope, not a correlation

Commit 3's defect scales every `logFC` by one constant, so `cor(reported, truth)`
stays exactly 1.0 and a correlation gate sees nothing. The check that sees it is

```
median( logFC / (rowMeans(tested[, numerator]) - rowMeans(tested[, denominator])) )
```

which reads **3.3219** against the defect and **1** when it is fixed.
Measured read-only on Yossi Tam A03 (PQN + glog10 + auto, 1887 features,
6 vs 6), re-measured on the rebased branch: slope `1.000000000000`,
`r = 1.000000000000` — and `r = 1.000000000000` against the defect too, which is
the whole point. `max|logFC - implied| = 0`. `mod_met_corrected()` reproduces
that lane's shipped `expr_tested.tsv` and `expr_normalized.tsv` at
`max|diff| = 0`, records `transform = log2` / `transform_configured = glog10`,
and its P-values match the shipped table exactly; the shipped `logFC` column
differs by up to `7.1e-3` only because the export path round-trips it through
`signif(linearFC, 3)` (`max|shipped - round_trip(ours)| = 0`).
`tests/testthat/test-metab-transform-recorded-as-applied.R` asserts the slope
and states the correlation, so nobody reaches for the gate that missed it.

## Blast radius

Nothing recomputes on its own. These are MOcore compute-time fixes: an already
imported project keeps the numbers it shipped with until it is deliberately
re-run and re-imported.

Seven imported projects carry a metabolomics or lipidomics layer.

| Project | Layers | Hit by | Status |
|---|---|---|---|
| `yossi-tam-masld` | metab + lipid | 1, 2, 3 | **Already correct.** Imported 2026-09-17 from the Yossi Tam clone with all three applied. Verified read-only: recorded `transform = log2`, `transform_configured = glog10`, slope 1.0. Nothing to do. |
| `elah-pick` | metab + lipid | 1 | Needs a re-run. Measured in commit 1: lipidomics 72 of 575 sign flips, flag 166 -> 110; metabolomics 9 of 662, flag 164 -> 159. (Not re-measured here — that data is Ifat's.) |
| `elad-chiel-spalangia` | metab | 1 | Needs a re-run. 441 features, `with_bacteria_vs_No_bacteria`: 77-101 sign flips, median abs delta 0.21-0.31 log2, flag 22 -> 14 / 22 -> 17 / 26 -> 21 across the A01-A03 variants; the EigenMS-forced A04/A05 variants set no FC cutoff, so their counts do not move (19 and 1 sign flips). |
| `liza-barki-cox2-exosomes` | metab | 1 | Needs a re-run. 310 features, `NT_vs_FL`: 85-108 sign flips and a median abs delta of 0.77-1.11 log2 — the largest directional change of any lane measured — while the flag count barely moves (14 -> 14, 19 -> 22, 20 -> 20). Directions on this project are the thing to re-check, not counts. |
| `lital-kalich-asd` | metab | 1 | Needs a re-run. 365 features, `ASD_vs_CONTROL`: 29 sign flips, median abs delta 0.047; no FC cutoff is configured so the significant count is unchanged (141 / 129 / 129 across m1-m3). |
| `itzik-kehat-mouse` | metab | 1 | Needs a re-run. 312 features, `F-mut_vs_Ctr` + `skNAC-KO_vs_Ctr`. The source MOcore run is not on this box, so the delta is not measured — it is in the same class as every other metabolomics lane. |
| `ascelegans` (Amir Sapir) | metab + lipid | 1 | Needs a re-run. Metabolomics 463 features, lipidomics 1794, `ppm1.56 vs ppm0`. The lipidomics lane is pareto-scaled, so its tested matrix is not among the shipped TSVs and the delta is not measured here. |

Only **one** lane in the whole fleet is hit by commits 2 and 3: Yossi Tam A03,
the only metabolomics run with variance scaling and the only one with a
non-log2 transform. Of 58 metabolomics/lipidomics MOcore runs found on disk,
58 are hit by commit 1, 1 by commit 2, 1 by commit 3.

Method: for every metabolomics lane with `scaling: none`, the shipped
`Datasets/expr_normalized.tsv` **is** the matrix the run tested, so the corrected
`logFC` is exactly its group-mean difference and the delta is computed from
shipped files with nothing re-run. Scaled lanes cannot be measured that way —
`expr_normalized.tsv` is the scaled matrix — which is why the two lipidomics
projects are unquantified above.
