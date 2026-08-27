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

## Landing note (2026-07-07, a833c09)

Both edges fixed at their sources. Bug A: refreshCutsForColumn refuses
a degenerate uniform re-cut (valuesAreDegenerate) and keeps the old
ascending grid - forced updates route new values through the retained
grid and collapse what empties (the plan's single-cut alternative
perturbed a live path, missing seed 21, and was dropped). Bug B:
dropStaleMissingDirections clears a rule's missing-direction bit once
its column stops holding missing values, called from the three
post-success refresh paths; the non-forced rollback now also restores
hasMissing so rejection stays gauge-consistent. Two direct component
tests reproduce the edges without the fuzzer and fail with the fixes
disabled. Gates: fuzzer clean over 100 seeds, full C++ suite,
tinytest 2468 ok, equivalence EXACT 18/18 (reviewer re-verified).
While probing a stricter round-trip assertion the agent surfaced a
THIRD pre-existing edge, out of scope here: setCutPoints shrinking a
grid can orphan a split whose out-of-range index makes flatten read a
stale value. Triaged to cutpoints-shrink-orphan.
