# csc-code-validation

agent: opus
rng: neutral
budget: ~250 lines (two bridge validators + tinytest; no engine touch)

## Goal

Every categorical test-ingestion entrance bounds category codes against
the column's K (store numCuts), CSC-backed training columns included.
The typed-ingestion probe (typed-ingestion.md Survey) verified three
silent-acceptance consequences on the declared-K path: deterministic
mis-bin (inline, K <= 63), shift UB at code >= 64 (maskTestBit,
tree.hpp), and a heap-buffer-overflow read (pooled, K > 63; ASAN,
standalone). All become the same refusal dense columns already get.

## Context

- validateCategoricalPredictors (src/R_interface_bartcore.cpp,
  creation-time x.test): skips any categorical column whose TRAINING
  source is CSC; the skip is keyed on the training column, so it also
  drops validation of a DENSE x.test matrix against that column.
- parseTestContainer (same file): bounds only cscReferenceCodes into
  [0, K); the stored values themselves are never bounded.
- validateColumnValues (same file; the mutation/test/predict entries):
  already bounds against store numCuts - the model to unify on. For a
  CSC column numCuts IS the declared K, so in-range
  declared-but-unobserved codes (probe G2) stay accepted.
- Existing out-of-range coverage is dense-only
  (inst/tinytest/test-bartcore.R, "existing category codes" block).

## Decision

Bridge-only, or also an engine-side refusal inside
quantizeCscColumnInto (cold ingestion path)? Recommend bridge-only:
validation lives at the boundary by design (the rc.a pattern), the
engine assumes valid codes everywhere else, and the typed-ingestion
view will centralize entrance validation anyway - a second engine
check would be a third place to keep in sync. Evidence that would
change it: the step 1 enumeration finding an entrance that cannot see
numCuts, or a non-bridge host (tests/cpp shim, future pybind) feeding
codes the bridge never sees.

## Constraints

- Refuse, never clamp: clamping reproduces the mis-bin class.
- The bound is numCuts (declared K for CSC, inferred max+1 for dense);
  do not tighten to observed support - probe G2/I1-in-range acceptance
  is correct behavior and must survive.
- Error text matches the existing dense refusals.
- rng neutral: previously-valid inputs take identical paths; only
  silently-accepted invalid inputs turn into errors.
- Out of scope: the dense declared-K asymmetry (typed-ingestion slice
  1), factors="indicators" collapse, engine defense in depth (per the
  Decision).

## Steps

1. Enumerate every entrance that installs categorical test/predictor
   codes (callers of parseTestContainer, validateCategoricalPredictors,
   validateColumnValues; the predict/setTestPredictor/setTestData/
   setData/setPredictor entries) and record which validator covers
   each; the fix must leave none uncovered.
2. validateCategoricalPredictors: remove the CSC-training skip; bound
   against the column's numCuts.
3. parseTestContainer: bound every stored categorical value into
   [0, K) beside the existing reference-code check; NA handling
   mirrors the dense validator.
4. tinytest: out-of-range refused at creation (dense x.test vs CSC
   train), via container setTestPredictor (inline and pooled K > 63),
   and via a foreign container with larger declared K; in-range
   declared-but-unobserved stays accepted (G2, I1 in-range).
5. Re-run the probe's silent-acceptance scripts (scratchpad
   typed-ingestion-probe/ 03, 04, 06, I1/K1 probes) and confirm each
   invalid case now refuses and each valid case still fits.

## Verification

- R CMD INSTALL .; full tinytest (3484 + new); air format --check .
  and local R CMD check (tinytest files are R).
- Equivalence trio bitwise via the dedicated harnesses (neutral
  commit).
- tests/cpp make && ./test_bartcore unchanged-green (no engine touch).
- CI including sanitizers to green before landing more on top.
