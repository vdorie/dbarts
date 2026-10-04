# input-guards: refuse infinite case weights and a cut count below one; keep pdbart's burn-in sigma

Status: LANDED 2026-10-04 (8090138b to 7d383099).

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. Each change refuses an input that never produced a valid fit, or fills a result component that
was NULL; no accepted call draws differently.
window: pre-release.
budget: ~250 lines (R ~40, C++ ~30, tinytest ~120, tests/cpp ~40, records ~20). Plans have run 1.5-2x low.

## Goal

A case weight of `Inf` is refused by name on every entry that takes case weights, as the 1.0-0 NEWS already
states ("negative or non-finite weights are refused everywhere"). A per-column cut count below one is refused
wherever the cap of 65533 is checked. A pdbart or pd2bart result carries `first.sigma` after a burn-in.

## Context

- Infinite weights today: [`validateXYWeights`](../../R/data.R) checks type and length only;
  [`refuseNonCountWeights`](../../R/spec.R) tests `w != round(w)`, which `Inf` passes. On a logistic fit the
  weight reaches the Polya-Gamma count and one sweep does not finish (observed: a 20-row, 10-sweep fit still
  running after 20 seconds; the interrupt poll runs between sweeps). On a gaussian fit the call fails with
  "unable to obtain a starting estimate of sigma", which does not name the weights. The sampler's
  [`dbartsSampler$setWeights`](../../man/dbartsSampler-class.Rd) and the bridge's own weight checks must
  refuse it as well.
- Cut count below one: the dbartsControl validity check refuses `n.cuts <= 0`, but the dbartsData `n.cuts`
  slot is checked for length only, and [`dbartsSpec`](../../R/spec.R) keeps a data object's slot. The bridge
  and [`maxNumCutsRepresentable`](../../src/bartcore/data.hpp)'s checks bound it above only. With zero in
  quantile mode the grid divides by zero (SIGFPE on x86-64); in uniform mode the sampler builds and then
  refuses its own stored state. 0.9-34 shared the division.
- pdbart: TODO item pdbart-burn-in-sigma. The burn-in run's sigma is stored as `first.sigma` and the result
  builder in [partialDependence.R](../../R/partialDependence.R) reads it under another name, so the
  component is NULL. 0.9-34 has the same line.

## Constraints

- No accepted input changes its draws; the seeded snapshot files, the equivalence baselines and the exact
  gates are untouched.
- Messages name the argument and the rule, in the package's refusal style.
- Out of scope: an upper bound on finite count weights or on the negative-binomial shape (a decision for the
  maintainer); the quantile-refresh cut thinning (state-frame-prior); guessNumCores returning NA.
- No NEWS entry: the weights rule is already stated there, the cut-count case is reachable only by editing a
  slot, and nothing in the package reads pdbart's `first.sigma`.

## Steps

1. Refuse non-finite case weights in the R validators every weighted entry reaches (`bart`, `dbarts`,
   `xbart`, `bartBT`, `rbart_vi`, `dbartsData`, `$setWeights`) and in the bridge, beside the existing negative
   check. Tinytest: `Inf` and `-Inf` on gaussian, logistic and probit, each refused by name before any sweep.
2. Refuse a cut count below one in the dbartsData validity check, the bridge's cap checks and
   ColumnStore::build. Tinytest: the slot route above refused; tests/cpp: build refuses zero.
3. Read the burn-in sigma under the name it is stored by; drop the TODO item. Tinytest: `pdbart` with
   `nskip > 0` on a continuous response returns a non-NULL `first.sigma` of the burn-in length.

## Verification

- `R CMD INSTALL` into the slice's own library, single job; the new and touched tinytest files pass;
  `tests/cpp` builds and `./test_bartcore` passes.
- `lintr::lint_package()`, `air format --check .`, `tools/check-doc-freshness.R` clean.

## Landing note

Landed 2026-10-04 as 8090138b (infinite case weights), fcf79de9 (cut count below one), 4a89236e (pdbart
burn-in sigma), 9a74bbf4 and 2605d509 (review fixes: a data object's weights slot reaching `dbartsSpec` or
`dbarts`, and posterior predictive weights on gaussian, student and aft fits), edbb0595 (comments), with
records 893ec0f8, a023ba46 and 7d383099. Review SOUND WITH CORRECTIONS twice, corrections taken. Full tinytest
12394 results, 0 failures; tests/cpp, lintr, air and the tools/ checks clean; R CMD check --as-cran one NOTE
(incoming feasibility). Out of scope and open: a huge finite logistic count weight or negative-binomial shape
still makes a sweep that cannot be interrupted.
