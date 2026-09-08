# Mummichog GSEA parity reference

`reference.rds` / `reference.tsv` hold the output of the **pinned MetaboAnalystR
Peaks-to-Pathways GSEA engine**, run on the fixture in `generate_reference.R`.
`test-mummichog-gsea-parity.R` asserts our implementation reproduces it.

| | |
|---|---|
| Repo | `xia-lab/MetaboAnalystR` |
| Commit | `398476ae2a0c996e390925fa62adc46b5ecde334` |
| Files | `R/util_fgsea.R`, `R/peaks_to_function.R` |
| Upstream entry | `PerformPSEA` → `.init.RT.Permutations` → `.compute.mummichog.RT.fgsea` → `fgsea2` → `my.fgsea` → `.run_fgsea_inner` |
| Config path | retention time present + version `"v2"` (what our pinned mummichog stage produces) |
| `permNum` | 100 (`PerformPSEA()` default) |
| `gseaParam` / `minSize` / `maxSize` | 1 / 1 / `Inf` (`.run_fgsea_inner()` defaults) |
| Seed | `set.seed(123)` before `sample.int(10^9, n_batches)` |
| fgsea at generation | 1.39.4 → the `> 1.24.0` branch (`stats` re-sorted decreasing after `abs()`) |
| fgsea in `renv.lock` | 1.36.2 → the same `> 1.24.0` branch |

Both reference files are **byte-identical between the pinned commit and upstream's
current default-branch head**, so this pin is the live code path, not a stale one.

## What the generator actually runs

`.run_fgsea_inner()` is `source()`d from the pinned checkout and called directly —
it is the entire GSEA engine (ES, permutations, NES, p, BH, leading edge). Only
its two `ov_qs_*` calls and `AddErrMsg` are stubbed; those are web-session I/O,
not maths. The ranked-input construction (pinned `peaks_to_function.R` lines
3470–3483) is copied verbatim into the generator **and cross-checked against our
own `mmc_gsea_ranked_inputs()`**, so a drift in either side fails generation.

## Fixture design

16 EmpiricalCompounds with signed scores from `+4.0` to `-4.0`: positives and
negatives, an exact tie (`E02`/`E03` both `3.2`), and a `+/-` pair of identical
magnitude (`E01` `+4.0`, `E16` `-4.0`) to exercise the ordering. Seven pathway
sets: top-of-ranking, bottom-of-ranking, middle, two overlapping sets, one of
size 1 (kept because `minSize = 1`), and one with no detected EC (dropped).

## Two upstream behaviours the fixture deliberately pins

1. **Tied EC scores ⇒ `ES = 0`.** Pathway members become tie-group *indices* and
   are never de-duplicated, so a pathway holding two ECs with the same score
   passes a non-strictly-increasing `selectedStats`, which
   `fgsea::calcGseaStat()` rejects (`all(head(S,-1) < tail(S,-1))`). Upstream's
   `tryCatch` swallows that and substitutes `list(res = 0, leadingEdge = c())`.
   `PW_top` and `PW_overlap_b` are those pathways, and the reference really does
   carry `ES = 0`, `NES = 0`, empty leading edge and one warning each
   (`upstream_warnings`). We reproduce it and flag it in the results table so a
   defaulted `0` is never read as a measured score.
2. **Leading edges are often empty.** `abs()` is applied *before* the sort, so
   pathway indices taken from the signed ordering address a differently-ordered
   vector; the leading edge is then intersected with the pathway's own ECs and
   frequently comes back empty. Reproduced, not "fixed".

## Regenerating

Only when deliberately moving to a new upstream version — never in CI:

```sh
git clone https://github.com/xia-lab/MetaboAnalystR /tmp/mar
git -C /tmp/mar checkout 398476ae2a0c996e390925fa62adc46b5ecde334
Rscript tests/testthat/fixtures/mummichog_gsea_parity/generate_reference.R /tmp/mar
```

Record the new commit here and in `R/domain/metabolomics/06g_mummichog_gsea.R`.
