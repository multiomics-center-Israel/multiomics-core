# Non-model organism multi-ORA — end-to-end GMT workflow

**Project:** `multiomics-core` multi-omics enrichment.
**Companion repo:** [`multiomics-annotation-prep`](https://github.com/multiomics-center-Israel/multiomics-annotation-prep)
**Updated:** 2026-07-26

For an organism **without a KEGG code / AnnotationDbi `OrgDb`** (e.g. a Trinity
transcriptome of *Coelastrella sp.*), `run_multi_ora()` cannot map features to
ENTREZ/KEGG. Instead of returning nothing it falls back to `run_multi_ora_gmt()`
(added in PR #119), which runs over-representation analysis against **custom GMT
gene sets** supplied per omic. This document walks the full path: from the
organism-specific KAAS/Trinotate inputs, through building the GMTs in the
companion repo, to wiring them into a `multiomics-core` config.

The GMTs are produced by `multiomics-annotation-prep`; `multiomics-core` only
consumes them.

---

## Overview

```
KAAS query.ko.txt ─┐
                   ├─▶ multiomics-annotation-prep ─▶ KEGG_pathway.gmt + GO_*.gmt
Trinotate GO ──────┘        (run where rest.kegg.jp is reachable)
                                          │
                                merge_gmt() + write_gmt()   (multiomics-core, R)
                                          │
                                one GMT per omic ─▶ modes.<omic>.pathway.gmt_file
                                          │
                                run_multi_ora() ─▶ run_multi_ora_gmt()
                                          │
                                multi_ora_results.csv + plots
```

---

## Stage 0 — organism-specific inputs (you provide these)

**`query.ko.txt`** — KAAS output, tab-separated `feature <TAB> K-number`:

```
TRINITY_DN101478_c0_g1_i1	K00844
TRINITY_DN10259_c5_g1_i1	K01810
```

**`trinotate_go.txt`** — Trinotate-style GO table (optional, for the GO sets):

```
TRINITY_DN101478_c0_g1	GO:0006096^biological_process^glycolytic process`GO:0005524^...
```

> KAAS entries are per-**isoform** (`_i1`). `prepare_kegg_nonmodel.py` strips the
> `_i\d+` suffix by default, so GMT members come out at the **gene** level
> (`TRINITY_DN101478_c0_g1`). This is the ID your DE `feature_id` must match
> (see Stage 3).

---

## Stage 1 — build the GMTs (`multiomics-annotation-prep`)

These commands **download live from KEGG (`rest.kegg.jp`) and GO**, so run them
where that host is reachable — your local machine or the repo's
`publish-organism-artifacts.yml` GitHub Action, **not** a restricted sandbox
(which returns a proxy `403`).

```bash
python -m venv .venv && source .venv/bin/activate
pip install requests pyyaml

# KEGG: KAAS KO list -> pathway gene sets
python scripts/run_kegg_nonmodel.py --kaas query.ko.txt --out results --cache data

# GO: Trinotate table -> GO gene sets (names + hierarchy expanded)
python scripts/run_go.py --go-table trinotate_go.txt --out results --cache data
```

Output in `results/`: `KEGG_pathway.gmt`, `GO_BP.gmt`, `GO_MF.gmt`, `GO_CC.gmt`
— standard GMT (`TERM_ID <TAB> description <TAB> gene1 <TAB> gene2 …`), directly
readable by `multiomics-core`'s `read_gmt()`.

---

## Stage 2 — merge into one GMT per omic (`multiomics-core`, R)

The GMT fallback reads a **single** `gmt_file` per omic, so combine the KEGG and
GO sets with `merge_gmt()`:

```r
# devtools::load_all(".")   # load R/core helpers
gmts <- file.path("results",
                  c("KEGG_pathway.gmt", "GO_BP.gmt", "GO_MF.gmt", "GO_CC.gmt"))
merged <- merge_gmt(gmts)                     # returns pathways with a descriptions attr
write_gmt(merged, "data/rna/pathways_coel.gmt")   # descriptions taken from the attr
```

Repeat (or reuse the same file) for each gene-based omic — see Stage 3.

---

## Stage 3 — align gene IDs (the critical step)

GMT members are gene-level Trinity IDs (`TRINITY_..._g1`). The fallback matches
DE `feature_id` values **directly** against them, so the ID schemes must agree:

- A leading **`Gene:`** prefix on RNA `feature_id`s is handled automatically —
  `run_multi_ora_gmt()` strips `^Gene:` before matching, mirroring the
  single-omics RNA pathway module.
- If **proteomics** is quantified against the same ORFs/genes (same Trinity
  IDs), point `modes.proteomics.pathway.gmt_file` at the **same** merged GMT.
  Otherwise build a proteomics-specific GMT keyed on the **same pathway IDs**
  (so per-omic hits pool under shared terms).

---

## Stage 4 — wire into the config (`multiomics_config.yaml`)

```yaml
global:
  organism: "Coelastrella sp."      # non-model -> get_kegg_organism / get_organism_db
                                    # return NULL -> triggers the GMT fallback

modes:
  rna:
    pathway:
      enabled: true
      gmt_file: "rna/pathways_coel.gmt"          # relative to paths.raw (data/), or absolute
  proteomics:
    pathway:
      enabled: true
      gmt_file: "proteomics/pathways_coel.gmt"    # or the same file if IDs are shared
```

`gmt_file` is resolved with `resolve_input_path()`: an **absolute** path is used
as-is, a **relative** one resolves under `paths.raw` (`data/`).

---

## Stage 5 — run

`run_multi_ora()` detects the missing KEGG/OrgDb organism, dispatches to
`run_multi_ora_gmt()`, runs pooled + per-omic ORA via
`clusterProfiler::enricher()`, and writes `multi_ora_results.csv` plus the
dot / support / pooled-barplot figures to the multi-ORA output directory.

---

## Constraints to know up front

1. **Two gene-based GMT omics are required.** A cross-omics ORA needs at least
   two omics with a usable GMT + significant features; with only one (e.g.
   RNA-only, or RNA + metabolomics where metabolomics is not covered by this
   gene-based fallback) `run_multi_ora_gmt()` declines and returns `NULL`. For a
   single omic, use the standard single-omics pathway path instead.
2. **Local path, not a URL.** The fallback checks `file.exists(gmt_file)`, so a
   GMT published as a GitHub Release must be **downloaded locally** and the
   config pointed at that path (unlike the mummichog model, which is fetched
   from a pinned URL).
3. **Metabolomics is a separate path.** The compound-centric mummichog model and
   the `*.compound_pathway.gmt` compound sets from the companion repo feed the
   metabolomics enrichment / mummichog steps — not this gene-based multi-ORA.

---

## References

- `R/domain/multiomics/07b_multigsea_plots.R` — `run_multi_ora()`,
  `run_multi_ora_gmt()`, `gmt_to_term2gene()`, `run_multi_ora_enricher()`.
- `R/core/13_gmt_utils.R` — `merge_gmt()`, `write_gmt()`, `filter_gmt_by_size()`.
- `R/core/09_enrichment.R` — `read_gmt()`.
- `multiomics-annotation-prep` — `scripts/run_kegg_nonmodel.py`,
  `scripts/run_go.py`, and `MODEL_CONTRACT.md` for the metabolomics artifacts.
