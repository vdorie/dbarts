# consumer-spec-surface

agent: Sonnet (R surface)
rng: neutral
budget: ~470 lines moved, ~120 new

## Goal

`dbartsSpec()` is exported: it resolves the `(control, model, data)` triple and
family token that describe a sampler, without constructing one. A
`LinkingTo: dbarts` package can then build a specification through supported
surface only - `dbartsData()` for the response half, `dbartsSpec()` for the rest -
and hand it to `dbarts_sampler_create`. Nothing about `dbarts()`'s behavior
changes.

## Context

Design: docs/design/consumer-spec-surface.md (the problem, the one-internal
factoring, the family-token trap, rejected alternatives).

`dbarts()` was the only place the resolution existed, and `bart2()`/`rbart_vi()`
reach it by `redirectCall`. The moved block ran from the classification-family
resolution through the variance-forest attribute, ending at the
`new("dbartsSampler", ...)` construction.

The motivating consumer is stan4bart, which hand-built a `dbartsModel` and called
`dbarts:::parsePriors` through a `quoteInNamespace` shim; see
[[dbarts-ecosystem-1.0-finalization]] gotcha 3 for what that cost at the 1.0 port.

## Constraints

- Draw-neutral. The move is contiguous and verbatim; the only edits inside it are
  parameter renames (`formals(dbarts)` -> `callFormals`).
- `dbarts()`'s own behavior is frozen, including the unconditional `n.cuts` and
  `sigma` overwrite it does before resolution (those lines stay in `dbarts()`;
  `dbartsSpec()` makes both conditional, which is a new-surface choice, not a
  change to the old one).
- Out of scope: multi-forest specifications (BCF, multinomial), the hazard
  family's person-period expansion, and any `dbarts_results` change.

## Steps

1. Move the resolution block from `dbarts()` into `resolveSamplerSpec()`
   (R/spec.R), taking the caller's `match.call()` and formals so the NSE prior
   parse still resolves in the caller's environment. Rewire `dbarts()` to call it.
2. Add `dbartsSpec()` over the same internal: validate `data`/`control`, resolve
   the aft `survival` status, honor an explicit `sigma`/`seed`, delegate.
3. Export it; document it (man/dbartsSpec.Rd, `_pkgdown.yml`, NEWS.Rd).
4. Cover it (inst/tinytest/test-spec.R), including the bitwise agreement with
   `dbarts()` and the features the surface newly reaches.

## Verification

- `R CMD INSTALL .`; `tinytest::test_package("dbarts")` - full suite green.
- `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-7903855.rds`
  - 27/27 identical draws (same RNG stream). This is the gate that proves the
  move was verbatim.
- `codetools::findGlobals(dbarts:::resolveSamplerSpec, merge = FALSE)$variables`
  names nothing but `parsePriors` and `pi`. A moved block silently captures its
  old frame otherwise; `dispersion` did exactly that on the first pass and only
  the nbinom equivalence scenario caught it.
- `air format --check .`; `R CMD check --as-cran`.

## Landing

LANDED 2026-07-25. Full suite 3468 (30 new), equivalence 27/27 bitwise,
findGlobals clean, check 0 errors / 0 warnings / 1 standard maintainer NOTE.

Follow-up 2026-07-27 (b96d3bb): `control` and `data` dropped from the prior
evaluation environment - they shadowed a caller's own variables inside prior
arguments and no default referenced them (design doc section 7). Regression
test in test-model-priors.R; suite 3469, equivalence 27/27 identical.

Follow-on, NOT part of this item: the per-observation working-weights/variance
field on `dbarts_results` (unblocks robust-t and heteroscedastic BART for an
embedded-Gibbs consumer) and a flat-API entry point for the BCF two-forest
sampler. Both are additive-MINOR by the ABI's own rules and can land in any
lockstep release; see stan4bart's TODO for the consumer side.
