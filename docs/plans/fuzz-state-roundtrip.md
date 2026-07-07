# fuzz-state-roundtrip

agent: opus
rng: neutral (both fixes affect only degenerate/mutated-state paths;
     verify no live scenario shifts)
budget: ~120 lines

## Goal

getState always round-trips: any state a live sampler can legitimately
reach restores through setState. The two edges the mutation fuzzer
found (mutation-fuzzing.md landing note, 44cee98) are fixed at their
sources.

## Context

- Degenerate ordinal grid: a forced setPredictor/updatePredictor with
  updateCutPoints over a near-constant column runs fillCutsUniformly
  over a zero range, installing a non-strictly-ascending grid;
  cutsWouldRemainValid never rejects in uniform mode. setState then
  refuses the state its own getState produced. Fix direction: a
  zero-range column gets the single-cut (or empty) grid the fresh-
  build path would produce; decide whether forceUpdate should refuse
  instead when the re-cut cannot produce a valid grid.
- Categorical gauge: category-changing per-observation/setData
  mutations can leave a live rule's mask non-canonical; flatten
  preserves it and buildFromFlat's canonical-gauge check refuses. Fix
  direction: canonicalize the mask at the mutation site (the gauge is
  an engine invariant, not a serialization detail) - flatten-side
  normalization would mask the live-state violation instead.
- Repro: tests/cpp fuzzer seeds 1 (probit/logistic) and 2
  (categorical); traces in the fuzzer output.

## Constraints

- setState's reject-with-nothing-changed contract stays.
- Out of scope: relaxing buildFromFlat validation.

## Steps

1. Component tests reproducing both edges directly (not via fuzzer).
2. Fix at the sources per the directions above.
3. Fuzzer seeds 1-100 clean; the two direct tests pass.

## Verification

- tests/cpp full suite + fuzzer 100 seeds; full tinytest;
  equivalence exact vs current baseline (paths are degenerate-only).
