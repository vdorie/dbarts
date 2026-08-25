# Pre-RC surface freeze: decisions

Status: DECIDED 2026-08-25 (VD, adopting the orchestrator's
recommendations as written; no work started). Work items carry TODO
entries named below. Evidence: the three pre-RC lens reports tracked at
docs/plans/review-2026-08-24/memos/prerc-lens1-surface.md,
prerc-lens2-backlog.md, prerc-lens3-external.md (findings cited by
their labels there).

Why now: 1.0-0 locks the R API and the flat C API. Each decision below
is a change that would be BREAKING after release and is cheap before
it. The utility rule applies: surface-bearing items with nameable value
land before the RC; additive items do not need to.

## D1. predict() argument shape (lens 1 B1, B2)

Decision: one signature order across every predict method,
`(object, newdata, type, ...)`, with the offset argument named
identically everywhere (`offset` for the training-shape offset,
`offset.test` retired as a spelling); `group.by` on rbart stays named,
never positional third. Rejected: refusing unknown names only (that
half landed as a defect fix; it does not remove the positional split).
Cost ~35 R + Rd + tests. TODO: predict-signature-unification.

## D2. keepTrees default (lens 2 P2)

Decision: keep `keepTrees = FALSE`; make predict()'s refusal on a fit
without trees name the one-argument cure (`keepTrees = TRUE`) in the
message. Rejected: flipping the default (every fit carries its trees).
Cost ~10 R. TODO: predict-refusal-names-cure.

## D3. Stub API-hash check (lens 1 A1)

Decision: the stub path checks what the header promises -
`apiMajorVersion() == MAJOR && apiMinorVersion() >= MINOR` - and the
hash equality check becomes the documented opt-in for lockstep
consumers. Rejected: leaving lockstep as the default (an additive
MINOR release would hard-error every stub consumer). Cost ~15 header
lines + the consumer-compile gate. TODO: stub-version-check.

## D4. Header type and naming drift (lens 1 A2, A6, A8, A9, A10)

Decision, all five before the freeze, one hash re-bake, one lockstep
rebuild of stan4bart and treatSens: `const int*` -> `const int32_t*`
where the same struct mixes them; `printEvery` -> `size_t`; the `get`
prefix dropped from `getLatents`/`getTrees` so the five readers share
one form (or added to all five - pick the form that keeps
`dbarts_sampler_getTrees`'s R twin readable); `printTrees` gains
`useLiveTrees` on `getTrees`' contract; a family accessor
`dbarts_sampler_family` returning the enum that replaces the
stringly-typed family. Rejected: deferring any (each is a MAJOR bump
after release). Cost ~25 header/bridge lines + tests. TODO:
dbarts-h-freeze-fixes.

## D5. Deprecation shims (lens 1 B3)

Decision: delete `predict.rbart`'s `value=` and `"post-mean"` shims;
nothing is released to deprecate from. Cost ~10 R. TODO:
deprecation-shim-removal.

## D6. Unadjudicated compositions (lens 2 P5)

Decision: grouped + `variance=` and heteroscedastic + `group.by` are
refused as validation errors by name, formals stay, with a door memo
naming what an adjudication would need. Rejected: shipping acceptance
of a composition nobody has checked. Cost ~20 R + memo. TODO:
composition-refusals.

## D7. BCF equivalence baseline format (lens 2 P7)

Decision: exempt the snapshot channels (mu, tau, glue, varcount,
forestFits, accepted, installed) under a cross-host flag so the
draws-axis channels alone gate statistically off-host; re-record
bcf-equivalence at the RC tip from the recording host with summaries.
Rejected: converting to draws-axis recordings (loses the bitwise
within-host channels). TODO: bcf-baseline-cross-host (absorbs the
existing cross-host half of equivalence-harness-statistical-mode).

## D8. NA in newdata on a training-complete column (lens 2 N1)

Decision: refuse by name at predict (the training column carried no
NA, so no missingness route was learned); `?dbarts` states the rule.
Rejected: documenting the silent left route. Cost ~20 R + ~10 C++.
TODO: predict-na-refusal.

## D9. Smaller surface items (lens 1 B5, B6, C1)

Decision: one `forest` convention across the four accessors (default
`NULL` = all, stacked, as `getForestAmplitudes`); `fitted()`'s `type`
vocabulary equals `predict()`/`extract()`'s on the same class; the 15
S3 methods with `\alias` and no `\usage` get usage entries (made moot
for the spellings by D1). Cost ~40 R/Rd. TODO: surface-smalls.

## Sequencing

D4 and D3 first (one hash re-bake; consumer rebuilds once), then D1,
D9, D5, D2, D6, D8, D7 at the RC tip. Each slice: Sonnet implementer
with the full gate battery; D4 an Opus design read on the enum first.
Post-1.0 by rule (additive): flat-API readers for heteroscedastic,
ordinal cutpoints, groupEffects count, setSigma getter; data-handle
serialization; pdbart on the new fit classes; variable-selection
inference and random-effects breadth (docs/plans/roadmap-survey.md).
