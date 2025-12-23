# CONTRIBUTING — multiomics-core

This document defines **how to contribute safely** to `multiomics-core` without degrading the architecture or breaking reproducibility.

For design philosophy and architecture, see 📘 `docs/developer_guide.md`

------------------------------------------------------------------------

## Scope

This guide applies to:

-   bug fixes
-   new pipeline steps (QC, preprocessing, DE, clustering)
-   refactors
-   config schema extensions

It does **not** cover how to *run* analyses (see onboarding).

------------------------------------------------------------------------

## Golden Rules (Non-Negotiable)

1.  **No analysis decisions in code**

    -   All tunable parameters must come from YAML config.
    -   Code implements *mechanics*, not *choices*.

2.  **No hidden I/O**

    -   Logic functions do not read/write files.
    -   File writing is explicit and tracked via `{targets}`.

3.  **Fail fast**

    -   Validate inputs and outputs immediately.
    -   Silent recycling or partial mismatches are bugs.

4.  **Do not break contracts**

    -   Existing function outputs must remain valid.
    -   Additive changes only, unless explicitly discussed.

------------------------------------------------------------------------

## Typical Contribution Types

### Bug fix

-   Add or strengthen validation
-   Fix logic without changing public API
-   Add a regression test *if possible*

### New functionality

-   New pure function in `R/`
-   Optional new target in `_targets.R`
-   Optional new config fields **with defaults**
-   Documentation update (required)

### Refactor

-   No behavior change without justification
-   Preserve function signatures or provide compatibility layer
-   Prefer small, reviewable diffs

------------------------------------------------------------------------

## Pull Request Checklist

Before opening a PR, confirm:

-   [ ] Code follows **single responsibility**
-   [ ] No file I/O in logic functions
-   [ ] New logic is callable interactively *and* via `{targets}`
-   [ ] All critical objects are validated
-   [ ] Naming conventions are respected
-   [ ] No breaking changes to existing outputs
-   [ ] New config fields have defaults
-   [ ] `tar_make()` runs end-to-end
-   [ ] Documentation updated if behavior changed

PRs that violate these rules will be requested for revision.

------------------------------------------------------------------------

## Versioning Expectations

-   **Patch (`x.y.z`)** — bug fixes, internal refactors
-   **Minor (`x.y`)** — new features, backward compatible
-   **Major (`x`)** — breaking API or output changes (rare)

If you are unsure, assume **minor** and ask.

------------------------------------------------------------------------

## What *Not* to Do

-   ❌ Copy legacy scripts verbatim
-   ❌ Add “helper” functions that do multiple things
-   ❌ Introduce mode-specific hacks into core utilities
-   ❌ Add randomness without config-controlled seeds
-   ❌ Add targets just to “make it run”

------------------------------------------------------------------------

## When in Doubt

Open a draft PR or ask:

-   “Where should this live?”
-   “Is this a function or a target?”
-   “Does this break a contract?”

Early questions save refactors later.
