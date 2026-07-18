# data-store-consolidation

agent: Opus
rng: neutral (pure refactor; every gate bitwise)
budget: staged; each stage its own commit, ~300 lines max apiece

## Goal

Consolidate the C++ data layer's accreted structure, per the
2026-07-17 four-reviewer data review: the train/test twinning
collapses into one per-side block over a shared cut grid, per-column
storage kind becomes explicit, and the duplicated scaffolding around
builds and transactions is written once.

## Findings this addresses (all verified in the review)

1. Test side is a copy-adapted twin, not a mirror: 6+ train/test pairs
   in data.hpp (quantizeCscColumn vs quantizeTestCscColumn 537-566 vs
   598-621; rank-bitmap builds 768-786 vs 933-951; columnIsCscBacked
   242-245 vs 570-573; sparse accessors 1217-1227 vs 1229-1241), plus
   the bridge's duplicated train/test container assembly
   (R_interface_bartcore.cpp: source-map loops 585-595 vs 417-425,
   CSC resolution 655-673 vs 453-466 with a byte-identical predicate,
   factor assembly in three copies). Carve: a CodeBlock (codes,
   offsets, sparse slots, slices, reference codes for one row set)
   instantiated train + test over one shared CutGrid (types, numCuts,
   cutPoints, maxNumCuts); one bridge container-parsing helper.
2. Per-column storage kind is smeared across five parallel arrays plus
   three bools, discriminated by empty-vector sentinels
   (columnIsCscBacked infers from builtFromCsc + emptiness). Carve:
   one per-column source descriptor (dense-owned | dense-borrowed |
   csc-rank | csc-densified, + slice + refCode).
3. Build-reset boilerplate restated per builder (~20 field resets x3
   at data.hpp:646-674, 702-739, 974-1021; test resets x4): 
   resetTrainStorage()/resetTestStorage() helpers.
4. Sampler transaction scaffolding duplicated (sampler.hpp
   setPredictor vs updatePredictor: precheck loop, snapshot set, and
   the requantize-test epilogue x4): one transaction helper over a
   column list; fold the store-owned rollback record here (the review
   found rollback invariants maintained by convention across layers).
5. isView conflates provenance with capability (view == no-raw ==
   refuse-mutation, data.hpp:129-133); split before any aliased
   sharing (data-ownership open plan 4) lands.
6. DataHandle gathers ALL raw columns unconditionally
   (R_interface_bartcore.cpp:1891-1896): ~400MB at n=1e6 p=50 that a
   constant-leaf consumer never reads. Gate the gather on the
   prospective leaf actually reading raw.
7. Linear-leaf u_ is column-major against per-member row access
   (model.hpp:389,408): row-major fits the gather; linear leaf only.
8. Honesty: quantizeColumn vs setColumnJournaled share their body -
   an observer-parameterized quantize core serves both (deferred from
   the remediation arc); the design docs' u8 memory claim is 2x
   optimistic while codes stay u16 (hot-layer-u8 remains open).

## Constraints

- Every gate bitwise at every stage (equivalence NOW 23/23 on
  31dc05a, bcf, multinomial, suite); no dbarts.h change.
- Stage order REVISED (2026-07-18, maintenance block): smallest-first
  4 -> 3 -> 2 -> 1. The transaction helper (4) and build-reset
  helpers (3) are low-risk and independently valuable; the storage
  descriptor (2) then makes the de-twinning (1) simpler, not the
  reverse. 5-8 ride their stages opportunistically.
- Stages 4 and 3 ran under VD's 2026-07-18 maintenance directive.
  Stages 2 and 1 APPROVED by VD after the post-maintenance
  discussion: both, serialized, 2 then 1; planner-first, one
  implementer at a time, every gate bitwise per stage.

## Landings

Stage 4 6841e6b (2026-07-18). runPredictorTransaction owns precheck,
gatheredRaw snapshot, forceUpdate collapse, validate-or-rollback, and
the requantize-test epilogue (was 4x); setPredictor/updatePredictor
are two-line shells over WholeMatrixUpdate (move-swap) / SubsetUpdate
(journaled, owns the ColumnCodeRollback records). Oracle mutation
tests + fuzzer UNMODIFIED. Net +24 (accepted: the strategy structs
cost more than two call sites save; the shape amortizes if stage 1
adds a third transaction form). All six gates bitwise.

Stage 3 ad78c06 (2026-07-18). resetTrainStorage (14 fields) /
resetTestStorage (8 fields) in data.hpp; reset sites: build,
buildMixed, buildFromParent / buildMixed tail, clearTest. Kept
divergences, verified field-for-field: buildMixed installs its CSC
kind after the reset; build repopulates gathered; buildFromParent
copies grid fields from the parent; active test BUILDERS stay inline
(they size, not clear, and buildTest deliberately preserves
testOffset). Net -6, no deviations, all six gates bitwise.

Stage 2 4f03854 (2026-07-17). ColumnSourceKind (denseOwned |
denseBorrowed | cscRank | cscDensified) + ColumnSource (kind, rankSlot,
denseRaw, slice, cscCategoryCount, refCode) in data.hpp; per-side
sources/testSources vectors, always sized numPredictors on a built
store, replace sparseSlot, cscSlices, mixedRawColumns,
cscCategoryCounts, cscReferenceCodes and their four test twins;
testBuiltFromCsc deleted (its only reader subsumed by kind). The
planner corrected two stale review premises: the test side already
had a full CSC/rank twin, and builtFromCsc/hasSparse have external
readers (bridge refuseViewSampler, tests/cpp) so both stay as
store-level flags; columnIsCscBacked and the sparse accessors keep
their signatures, so no file outside data.hpp changed. Finding 5
(isView) deferred: a whole-store capability axis, sequenced against
the data-ownership plan, orthogonal to per-column kind. Deviation
accepted: cscCategoryCount/refCode assigned only in the CSC branch
(dense-backed categoricals take the max-value scan and never read
them). Descriptor designed per-row-set so stage 1's CodeBlock absorbs
it verbatim. Net +16 (224 lines changed), under budget. All six gates
bitwise; orchestrator re-ran component binary, suite (3098/0),
equivalence (23/23 identical draws), bcf, and multinomial compares
independently. Note: git dates stages 4/3 2026-07-17; the labels
above were written a day ahead.

Stage 1 (train/test de-twinning) remains open - approved, next in the
serialized order.
