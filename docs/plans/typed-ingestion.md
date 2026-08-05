# typed-ingestion

agent: opus
rng: per slice (0 and 2 neutral; 1 shifting for gap data only)
budget: staged; slice 0 spun out to csc-code-validation.md

## Goal

One borrowed typed PredictorSource view (per column: dense-double |
dense-integer | CSC; ordinal | categorical with declared K) taken by
creation, test ingestion, and mutation, collapsing SamplerOptions' 8
ingestion fields, the 7-arg setTestData, and the dense-only mutation
entries into one spelling. Declared K becomes universal ahead of the
view (the TODO's pre-registered rule, triggered by the probe below).

## Survey (probe record, 2026-08-05)

- REACHABLE. The default surface keeps factors unexpanded
  (dbartsData factors="categorical", drop.unused.levels FALSE), so a
  declared-but-unobserved top level is ordinary. Every dense entrance
  then REFUSES test/mutation rows carrying it - the bound is inferred
  max+1 (buildCutsForColumn, data.hpp) while the R layer correctly
  maps against declared levels (mapFactorColumnsToTrainingLevels,
  R/utility.R). A refusal defect, not the premised mis-bin. Only the
  top level matters; gaps below are covered.
- MIRROR HOLE, severe: the CSC declared-K path never bounds stored
  codes. validateCategoricalPredictors (R_interface_bartcore.cpp)
  skips any column whose TRAINING source is CSC - dense x.test against
  a CSC column included - and parseTestContainer bounds only the
  reference code. Verified: inline K <= 63 silent mis-bin (descends
  LEFT, bit-identical to a real category); code >= 64 shift UB
  (maskTestBit, tree.hpp); pooled K > 63 heap-buffer-overflow read
  (ASAN, standalone repro of the exact store/tree state). Spun out as
  the URGENT slice 0: csc-code-validation.md.
- Exactly one declared-K channel exists (ColumnSource::cscCategoryCount
  via SamplerOptions.cscCategoryCounts; provenance sparseFactor@levels,
  R/mixedMatrix.R) - the right shape, CSC-only today. Dense declared
  levels die at the boundary: makeCategoricalModelMatrix collects the
  tables (attr "factor.levels", R/utility.R) but parseData never reads
  them; dbartsData has no predictor-levels slot.
- Plumb-through touch list: R makeCategoricalModelMatrix ->
  assembleDenseColumnMatrix (R/mixedMatrix.R) -> dbartsData initialize
  (R/data.R); bridge parseData -> optionsFromParsed ->
  validateCategoricalPredictors; engine SamplerOptions dense analogue
  of cscCategoryCounts -> ColumnStore::build -> buildCutsForColumn
  branching on declared-vs-inferred. dbarts.h has zero categorical
  surface: no ABI cost.
- Recorded, out of scope: factors="indicators" silently collapses an
  unobserved top level onto the reference (pre-existing classic
  behavior); engine ingestion itself is unchecked (codeFor is a bare
  cast; categoricalValueIsValid is bridge-only).

## Decision

Slice order, per the pre-registered rule (reachable -> declare K ahead
of the full view): 0. csc-code-validation.md (bounding holes; memory
safety; own plan, first). 1. dense declared-K plumb-through along the
touch list above, making G2/G4 symmetric (dense factors with an
unobserved top level accepted like CSC). rng: datasets whose factors
observe every level keep an identical cut grid (K == max+1), bitwise
neutral - verify the equivalence-trio data qualifies; gap datasets
gain the declared grid, class shifting, and are the point. 2. the full
PredictorSource view plus the sparse-extensions mutation-shape lift
(CSC-backed mutation = ownership policy for the borrowed triple + a
slice resize). Slices 1 and 2 get step lists here before
implementation; VD forks on this ordering.

## Constraints

- Typing stops at the boundary: the hot layer stays quantized codes +
  double cuts; response/offset/weights stay double*.
- Validation keyed on the view, never on the storage kind (the CSC
  skip is the standing counterexample).
- Slice 1 must not change the grid when declared K == observed max+1.

## Verification

Per slice; slice 0's lives in csc-code-validation.md. Slice 1 adds:
equivalence trio bitwise (after verifying baseline data observes all
levels), a gap-data tinytest (train with unobserved top level, test
row carrying it accepted, fits finite), R CMD check (R surface moves).
