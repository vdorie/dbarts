# restore-empty-leaf: a state install merges leaves no row reaches, as 0.9-34 did

Status: IMPLEMENTED 2026-10-04 (276eddce, 84babf15), on wt/restore-empty-leaf, not landed.

agent: opus implementer, one; opus reviewer.
rng: SHIFTING. Draws change only for a same-grid warm start whose donor leaves a leaf with no rows, which
today installs the empty leaf and now merges it; the posterior does not change. Every `setState`, `copy` and
reload accepted today has no empty leaf, so nothing merges there and its draws are bitwise unchanged; installs
refused today now proceed.
window: pre-release. A regression against 0.9-34.
budget: ~500 lines (C++ ~150, tests/cpp ~150, tinytest ~150, records ~50). Plans have run 1.5-2x low.

## Goal

Installing a saved state whose live trees route no row of the current data to some bottom node succeeds: the
empty nodes are merged into their parents, the handling a forced `setPredictor` already applies. A sampler can
always restore its own state and `copy()` itself after its predictors changed.

## Context

- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) refuses any live tree, mean or variance, with an
  unoccupied bottom node. After `setPredictor(x, forceUpdate = TRUE)` on an equal cut grid, a state stored
  before the change is refused by `setState` and by `copy()` ("state is not consistent with this sampler"):
  0 of 10 restores succeed on the bartcore tip, 10 of 10 on 0.9-34, whose restore ended in
  `forceUpdateTrees` and collapsed empty nodes.
- The collapse a forced predictor update runs, including the merge of constant, vector and variance leaves,
  gp per-observation fits, and [`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp) for a
  monotone forest, is the handling to reuse.
- Install paths: `setState`, `copy`, reload from a stored state, `installTrees` and `warm.start`.
- The warm start does not refuse: a same-grid donor whose mean trees leave a leaf with no rows of the
  destination's data is installed with the empty leaf live
  ([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) only repartitions; the same-grid occupancy check
  in [`Sampler::installForests`](../../src/bartcore/sampler.hpp) covered the variance forest alone). The empty
  leaves persist through later sweeps, and the warm-started sampler then cannot restore or `copy()` itself.

## Constraints

- No accepted `setState`, `copy` or reload changes: the seeded snapshot files, equivalence baselines and exact
  gates are untouched. A same-grid warm start that installs an empty leaf today merges it instead.
- No new user signal (warning or return value) in this item: 0.9-34 merged silently, and how a caller learns
  that a restore was not exact is a separate maintainer decision.
- A state that is malformed for other reasons (tree shape, sizes, cut indices out of range) is still refused.
- No NEWS: 0.9-34 restored the same way, and `warm.start` and `installTrees` are new in 1.0-0 (not in
  0.9-34).
- Out of scope: states drawn on a different cut grid or standardization (state-frame-prior).

## Steps

1. Let the install accept an unoccupied bottom node and collapse after repartition on every install path,
   reusing the forced-update merge and the monotone reseed. tests/cpp: collapse on constant, linear, gp,
   variance and monotone forests, and on a multi-forest sampler.
2. tinytest: a sampler's own stale state after a forced `setPredictor` restores through `setState`, `copy()`
   and a save-reload, then runs; a multinomial and a BCF sampler do the same; a restore accepted before the
   change is bitwise identical in its next draws. A same-grid warm start leaving an empty leaf yields a
   sampler with no empty leaf that restores and copies itself.
3. Mutation check: reinstating the refusal makes the new tests fail.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes; the
  four seeded snapshot files pass unchanged.
- The equivalence compare in statistical (z) mode against the current baseline: no scenario warm-starts, so
  every scenario is expected to report identical draws.
- `lintr::lint_package()`, `air format --check .`, and the three tools/ checks clean.
