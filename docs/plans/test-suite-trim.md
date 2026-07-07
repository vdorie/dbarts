# test-suite-trim

agent: sonnet
rng: neutral (test-only; no engine change)
budget: ~150 lines of test edits

## Goal

The tinytest suite spends its time where it buys coverage. The top
four files carry 52% of a 33.8s suite (profile 2026-07-06:
updatePredictorPerObservationJointly 6.9s, multithreaded 3.9s,
simd 3.5s, setPredictorPerObservation 3.4s); the target is the same
asserted contracts at roughly 20s total.

## Context

- test-sampler-updatePredictorPerObservationJointly.R: a 30-trial
  property loop (line 218) plus ten short run() calls over small
  fits; the trials drive a statistical rollback check.
- test-simd.R: every kernel crossed with every runtime instruction
  set; sizes are the multiplier.
- test-multithreaded.R: real thread pools; setup cost dominates.
- The wider time sinks are NOT tinytest: tests/cpp build (~3 min),
  equivalence record (~12 min), bench-sampler (~10 min) are the long
  gates and are inherent.

## Constraints

- An assertion's false-positive and false-negative rates must not
  silently worsen: reducing trials on a calibrated statistical bound
  requires recomputing the bound, or leaving that check alone.
- Deterministic contracts (rollback exactness, well-formedness,
  reshape behavior) keep full coverage - trim repetition, not cases.
- simd: keep every instruction set (the dispatch is the thing under
  test); trim vector sizes to one short and one cutoff-straddling
  length per kernel.
- No engine or R source changes; inst/tinytest only.

## Steps

1. Re-profile per file; confirm the top four.
2. Trim iteration/replication fat per the constraints; recompute any
   statistical bound whose trial count changes (document the
   calibration in the test comment).
3. Full suite passes; total wall time reported before/after in the
   landing note.

## Verification

- tinytest::test_package("dbarts") all ok; per-file profile shows
  the top-four share materially reduced with case coverage intact.
