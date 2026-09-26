---
name: run-claims-audit
description: Audit claims about what a pipeline run actually did, where configured intent and executed behaviour can diverge. Use when writing or reviewing report or Methods prose that describes a run; when one config or statistical decision is read by more than one consumer; when a review finding suggests a described behaviour was never executed; or when scoping remediation after a statistical bug. Not for ordinary config-driven code whose behaviour is not being described or claimed elsewhere.
---

# Run-claims audit

A claim about a run is only as good as the thing you checked it against. The
failure mode is not a wrong number — it is a sentence, a branch, or a test that
reads as verified because it names a config key, while the run did something
else.

`CLAUDE.md` carries the standing rules: configuration is intent not evidence,
one semantic decision gets one resolver, reports describe rather than decide,
and wording may not exceed what the implementation supports. This skill is how
to apply them.

## 1. Intent vs behaviour

For every sentence or branch keyed on configuration:

- Find the consumer that acts on the value. Copy its resolution exactly — key,
  default, fallback. A default in a template is not the default in the code,
  and the two disagree more often than you would expect.
- Check whether the step can decline to run despite the flag being set.
- Record which parts you verified against source and which you assumed. An
  assumption stated is fine; an assumption presented as a check is not.

## 2. What counts as evidence

| Signal | Proves |
|---|---|
| Config flag set | Intent only |
| Runtime result or returned state | Behaviour |
| Produced artefact | Behaviour, but only under both conditions below |

An artefact is evidence only when:

1. **Its ownership and provenance tie it to the current run.** A file in the
   expected location is not automatically this run's output. Establish that the
   run being described is the run that wrote it.
2. **Stale output cannot survive a skipped rerun.** Confirm the producing step
   removes the artefact on every path that can skip rewriting it. Where it does
   not, a previous run's output will read as proof that a computation ran — and
   will often be displayed as though it were current.

If either condition is missing, the artefact proves only that the file exists —
not that the current run produced it or that it represents current behaviour.

## 3. Shared contracts

Before changing how one consumer reads a setting, enumerate **every** consumer
of that decision — statistical behaviour, exports, and reader-facing text
alike.

Fixing the site where a bug surfaced while leaving other consumers reading raw
keys does not close the problem; it converts one silent error into a different
silent error, usually a harder one to notice.

Treat these as four separate concerns sharing one source of truth beneath them:

- compatibility aliases;
- validation and deprecation policy;
- runtime resolution (pure — no warnings, so repeated calls stay quiet);
- reader-facing reporting.

They drift when any consumer resolves the decision for itself.

## 4. Config allowlists

Derive a permitted-key list from a **runtime-consumer audit**, never from a
template alone. Templates routinely omit supported keys, so a template-derived
allowlist rejects working configurations.

Use the template as an *additional* regression fixture: every key in it must
validate. It is a lower bound on the schema, not the schema.

A config key that changes scientific behaviour must be validated. An
unrecognised key in such a block should fail, or have an explicitly supported
compatibility path — never be silently ignored, because the silent value is
frequently not even the documented default.

## 5. Matrix provenance

Keep analysis matrices semantically distinct and name which one you mean:
measured values, the matrix a model was actually fitted on, and any matrix
imputed for quality control, visualisation or an analysis stage are different
objects with different provenance.

Never let one stand in for another — in code, in a plot, in an export column,
or in a sentence. Where a figure or table mixes them, say so explicitly rather
than letting the reader assume a single origin.

## 6. Two gates, separately

Source correctness and render/output QA fail differently, and neither
substitutes for the other. Code can be correct and render wrongly; output can
look right and be produced by the wrong branch.

State which gate you performed. If you could not perform the other, say so
plainly rather than wording the report so that both appear covered.

## 7. Assertions that can actually fail

Before trusting a regression test, ask what incorrect output it would reject.

- **Negative substring traps.** A forbidden string must not be a substring of
  the permitted one, nor the reverse. `"adjusted"` is contained in
  `"unadjusted"`, so the obvious negative assertion either fails on correct
  text or passes on incorrect text depending on direction. Check both.
- **Prefix traps.** A pinned prefix of a longer correct sentence also matches
  the shorter wrong one, so it pins nothing.
- **Producibility.** When pinning generated text, confirm the literal is
  actually emittable by the generator — with interpolation applied — before
  relying on it.
- **Source-scanning tests.** A test that greps source for a call site can be
  satisfied by a *comment* containing the same text. Check what the pattern
  matches, not what you meant it to match.

After changing generated wording, re-check the assertions that already existed,
not only the ones you added.

## 8. After a statistical bug, separate three questions

They differ in cost and in necessity, and conflating them either over- or
under-scopes remediation:

1. **Recomputing the decision** from retained intermediates.
2. **Re-fitting the statistical model** — often unnecessary when sufficient
   per-stage outputs were retained and only a selection rule was wrong.
3. **Regenerating downstream deliverables** — tables, figures, enrichment
   inputs, exports, rendered reports. These may need regenerating even when (2)
   is avoided, because the decision they were built from changed.

Establish the direction of the error where you can: an error that can only
widen or only narrow a result set is far cheaper to bound than one that can do
either.
