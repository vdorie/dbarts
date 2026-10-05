# restore-status: setState returns whether the restore was exact

Status: PLANNED.

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. Only a return value is added; no install changes.
window: pre-release (dec-B234).
budget: ~350 lines (C++ ~120, bridge and R ~40, tests/cpp ~80, tinytest ~80, manual ~30). Plans have run
1.5-2x low.

## Goal

`setState` invisibly returns `TRUE` when nothing had to be changed to install the state, and `FALSE`
otherwise. Code that stores the state, proposes a change and restores on rejection can test the value: an
inexact rejection there is quietly wrong. No warning is raised. `copy()` and a reload, which have no return
value to carry it, stay silent.

## Context

- A restore can leave a chain that is not the stored one when the sampler's data changed after the state was
  stored: leaves no row reaches are merged, a missing-value direction a column no longer routes is dropped,
  and a state stored under another response transform is converted into the sampler's.
- `setState` returns `NULL` invisibly today. No test asserts that value and no R code reads it.
- The value exists for code that restores in order to reject a proposal, and such code needs "exact", not
  "to rounding": a conversion therefore gives `FALSE` as a merge does. That is the orchestrator's refinement
  of the ruling, told to the maintainer.

## What gives FALSE

For any live tree of any chain, mean forests and the variance forest alike:
- a bottom node no row reaches was merged into its parent (a monotone tree reseeded after it included);
- a rule's missing-value direction was dropped because its column holds no missing value now;
- any chain's leaf values were converted between response units.

Two more causes arrive with later work and report through the same flag: a split moved onto another grid, and
slopes converted between covariate standardizations. The flag is built so that a later install path sets it
without touching the bridge or R again.

## What does not

`TRUE` means "the chain as stored", not "bit for bit as if never stored". These leave it `TRUE`:
- latents re-derived because the case weights differ, or redrawn because an aft status differs;
- a value the sampler holds fixed, or draws but the state lacks, left as the sampler's;
- a generator of another kind not installed;
- a pooled categorical rule keeping a missing-value bit its column no longer routes;
- anything about kept draws, which are copied and never rebuilt.

## Constraints

- Only the live build counts. The checks that build scratch trees before a live tree is touched must not set
  the flag, and a refused install returns nothing (it stops).
- The flag covers every chain: one inexact chain gives `FALSE`.
- No install changes: draws after any `setState`, `copy` or reload are bit for bit what they are now.
- `installTrees` and a warm start keep returning what they return: a warm start is a seed, not a restore.
- The flat C header has no state entry and does not change.

## Steps

1. Engine: one out-flag from the sampler's state install
   ([`Sampler::setState`](../../src/bartcore/sampler.hpp)), set by the merge branches of the live rebuild
   and its variance twin ([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp)), by the units pass
   ([`Chain::convertStateUnits`](../../src/bartcore/chain.hpp)) when the units differ, and by
   [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) when it clears a direction on a live build. The facade
   virtual changes, so the build is `--preclean`. tests/cpp: the flag set by a merge, a dropped direction, a
   variance tree's merge and a units conversion; clear for a sampler's own fresh state on each leaf model, on
   a multi-chain sampler, and across a weights mismatch.
2. Bridge and R: the bridge's state install returns the logical; `setState` returns it invisibly; `copy` and
   the re-creation after a reload discard it.
3. tinytest: a value beside each existing silent restore (test-state-empty-leaf-merge.R,
   test-heteroscedastic-mutation.R and test-state-missing-direction.R among them): `FALSE` where it merges,
   drops or converts units, `TRUE` otherwise; `withVisible` shows it invisible; a restore that puts the
   predictor back first with `forceUpdate = TRUE` gives `TRUE`, and the other order `FALSE` when it merges.
4. Mutation check: a flag that is always clear fails the new tests.
5. Manual, `setState`: the value, what gives `FALSE`, and that `TRUE` does not promise the latents or the
   generator.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes, clean
  under ASan and UBSan.
- The four seeded snapshot files pass unchanged on a reference build.
- The equivalence compares in statistical (z) mode against the current baselines: every scenario identical.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean.
