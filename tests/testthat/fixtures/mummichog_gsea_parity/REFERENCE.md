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
| fgsea at generation | **1.36.2** → the `> 1.24.0` branch (`stats` re-sorted decreasing after `abs()`) |
| fgsea in `renv.lock` | 1.36.2 — same 1.36 series |

Both reference files are **byte-identical between the pinned commit and upstream's
current default-branch head**, so this pin is the live code path, not a stale one.

## Why the fgsea version is pinned too

The engine calls fgsea **internals** (`calcGseaStatCumulativeBatch`), so "same
`> 1.24.0` branch" is not sufficient evidence of numerical parity — the exact
build matters. Three states, enforced rather than described:

| installed fgsea vs `renv.lock` | generator | parity test |
|---|---|---|
| exact match | generates, `fgsea_exact_match = TRUE` | compares |
| same `x.y`, different patch | **refuses** unless `MMC_PARITY_ALLOW_FGSEA_PATCH_MISMATCH=1`, and stores `fgsea_mismatch_note` | compares, and asserts the note exists |
| different `x.y` series | **refuses** | **skips** — a cross-series match is not parity |

A further test asserts that whenever the suite runs somewhere the **locked build
is installed**, the fixture must have been generated on it. The committed fixture
now satisfies that requirement.

### Current locked-version state

`fgsea_exact_match` is **TRUE**: the committed fixture was regenerated under the
exact `renv.lock` build, **fgsea 1.36.2**.

The regenerated fixture produced a reference table that is identical to the
previous 1.36.0 fixture (`all.equal(old$reference, new$reference) == TRUE`).
The MetaboAnalystR commit, permutation count, GSEA parameters and seed are unchanged;
only fixture provenance changed (`fgsea_version`, `fgsea_exact_match`,
`fgsea_mismatch_note`, and generation timestamp).

The earlier 1.36.0 work remains useful as historical validation context.
In the original authoring environment, the exact 1.36.2 build could not be obtained:

| route | outcome |
|---|---|
| `bioconductor.org` (release + Archive) | denied |
| `git.bioconductor.org` | denied |
| `cloud.r-project.org`, `packagemanager.posit.co` | denied |
| `bioc.r-universe.dev`, `api.anaconda.org` | denied |
| `depot.galaxyproject.org` (bioconda's source mirror) | denied |
| `alserglab/fgsea` `RELEASE_3_22` branch | not mirrored (branches stop at `RELEASE_3_21`) |
| `Bioconductor/fgsea`, `Bioconductor-mirror/fgsea`, `ctlab/fgsea` on GitHub | do not exist / no such branch |

Historical validation performed before the exact 1.36.2 build became available:

* the two primitives this engine calls were code-identical from 1.36.0 through
  1.39.4: `calcGseaStat` body md5 `9cac7ddd36b45885` at both revisions
  (pure R, no `.Call`), and `src/fastGSEA.cpp` with comments stripped md5
  `8060947310be0650` at both. The one post-branch-point commit touching
  `fastGSEA.cpp` (`7207d6b`) adds three comment lines and nothing else;
* regenerating the whole fixture under 1.39.4 and under 1.36.0 gave **bitwise
  identical** ES, NES, p, padj, `nMoreExtreme`, sizes, pathway order and leading
  edges.

The exact 1.36.2 regeneration now supersedes that historical bounding evidence:
the committed reference table is identical to the previous fixture while
`fgsea_exact_match` is `TRUE`.

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
   (`upstream_warnings`). We reproduce it and flag it in the results table with a
   **generic** `ES defaulted by reference implementation` column plus an
   `ES fallback reason` classified from the actual condition — so a defaulted `0`
   is never read as a measured score, and a non-tie failure is never mislabelled
   as a tie.
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

Install the exact `renv.lock` fgsea build first — the generator aborts otherwise
(and only accepts a patch-level difference with
`MMC_PARITY_ALLOW_FGSEA_PATCH_MISMATCH=1`, which records the deviation in the
fixture). Record any new MetaboAnalystR commit here and in
`R/domain/metabolomics/06g_mummichog_gsea.R`.
