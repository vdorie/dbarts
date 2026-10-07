# state-family-record-removal: a state stops recording its family and is not refused by family

Status: LANDED 2026-10-07 (14cfd187 to 7359f4bc; dec-B283, revising dec-A174).

agent: opus implementer, one (bridge, tests, records).
rng: NEUTRAL, bit for bit: no sweep changes. A stored state loses one attribute and `setState`, `copy` and a
reload lose one refusal. Proved by a seeded round trip on the tip's build and the slice's (Verification).
window: pre-release, before the merge to main. One bridge file, in the state's writer and reader.
budget: ~300 changed lines, mostly removals and tests.

## Goal

`setState` does what it can with no guarantee (dec-B283). A state carries no `family` attribute, and
`setState`, `copy` and a reload install a state whose blocks fit the sampler whatever family stored it. They
still refuse what is plainly broken: scales or Polya-Gamma variates that are not positive and finite.

## Context

[cross-family-state-install.md](cross-family-state-install.md) landed a family record on every state with a
refusal by name, and the floor on precisions, [`ResponseModel::canHoldLatents`](../../src/bartcore/model.hpp).
Only the bridge held the record; the engine and tests/cpp hold only the floor. The record was never released.

## Constraints

- The floor and its tests stay as they are. No new mechanism: no redraw on install, no warning.
- A state that carries a `family` attribute still installs: the attribute is not read.
- No engine, R or NEWS change; [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move.

## Steps

1. Bridge: [`storeState`](../../src/R_interface_bartcore.cpp) stops writing the attribute and
   [`setState`](../../src/R_interface_bartcore.cpp) stops reading it: the file is what it was before the record.
2. tinytest. [test-state-family.R](../../inst/tinytest/test-state-family.R) drops the record, the refusals
   by name and the malformed records; keeps the installs within one family, the floor by hand, by chain,
   under a mask and through `copy` and a reload, and the warm start; and gains: a state names no family; a
   `family` attribute of any value is not read; each of the 15 pairs that broke a sampler is refused and
   installs once its latents are made positive. The pin
   ["a gaussian state leaves a probit sampler's sigma, pinned at 1, where it is"](../../inst/tinytest/test-state-not-model.R)
   installs its state as stored again.
3. Records: the sampler's help, [state-not-model.md](../design/state-not-model.md),
   [public-surface.md](../design/public-surface.md) and the plan above say what stands.

## Verification

- Round trip on both builds, two chains, for a gaussian, probit, logistic, Student-t and negative-binomial
  sampler: draws after a store, a `setState` into a second sampler and a `copy` are identical across builds.
- The 400 ordered pairs of twenty kinds of sampler through `setState`, `copy` and a reload, each in a process
  of its own under a time limit: 325 of the 380 across kinds refused, 55 run, every one returns, fits finite.
- tests/cpp and the full tinytest suite at home; lintr, air, the rc-codoc, win-drift and doc-freshness checks
  and the mutation battery's anchors; `R CMD check --as-cran` on a clean export; stan4bart's suite.

## Landing note

Landed 2026-10-07 as 14cfd187 to 7359f4bc, 247 lines added and 432 removed. The bridge file is byte for
byte what it was before the record was added; the engine and tests/cpp are untouched, every line the
earlier slice added there serving the floor. The independent review, told to refute, found nothing
blocking: over 400 ordered pairs of twenty kinds of sampler, by `setState`, `copy` and a reload, none
hung and none gave fits that were not finite; a state written with the record installs with it
ignored; and two mutations of the floor each failed the tests. What the review found and this slice
leaves: a value edited by hand to the edge of what a double holds, a Polya-Gamma variate of 1e-320 or
the largest finite double as a scale, passes the floor and breaks a sweep, within one family as
across two (root TODO, precision-floor-edge-values). Gates on a clean copy of the rebased branch:
install, tests/cpp 351, the tinytest suite at home 16370 results and none failed, lintr, air,
rc-codoc, win-drift, anchors, build and check with the Date note alone. stan4bart's suite against the
build: 582 results, none failed.
