# data-store-residuals

agent: opus
rng: neutral - refactor/layout/capability passes; draws unchanged; every
  gate bitwise.
budget: ~250-350 lines across ~8 files (data.hpp, model.hpp,
  R_interface_bartcore.cpp, the R handle surface,
  tests/cpp/test_data.cpp, docs/design/data-store.md).

## Goal

The four findings deferred out of the 2026-07-17 data review are
discharged: the journaled column set shares the quantize core, isView
carries provenance only with capability made explicit, DataHandle
gathers raw only for declared leaf-covariate columns, and the linear
leaf's u_ matches its row access. data-store.md reflects each pass.

## Context

docs/design/data-store.md is the spec. Anchors re-verified 2026-07-18:

- Quantize twin: quantizeColumn (data.hpp:644-655) rides
  quantizeDenseInto; setColumnJournaled (data.hpp:1134-1162)
  re-implements the same per-cell loop inline with journal bookkeeping.
  Sole caller sampler.hpp:1096 (SubsetUpdate), reached once per external
  iteration via updatePredictor/updatePredictorPerObservationJointly.
- isView: field data.hpp:214; set at :720, :768, :1001 (buildFromParent);
  provenance reader suppliedStandardization (data.hpp:342, plus the
  linear/GP reinitialize in model.hpp); capability readers
  refuseViewSampler (R_interface_bartcore.cpp:1479-1487, ORs
  builtFromCsc, 5 mutation entry points) and refuseViewSamplerOnly
  (:1491-1496; setCutPoints/setState/installForests);
  tests/cpp/test_data.cpp:72. Doc language: data-store.md:84-85,
  :181-189.
- Handle gather-all: R_interface_bartcore.cpp:2026-2034 in
  bartcore_createDataHandle (:1997). The sampler ctor
  (sampler.hpp:100-127) and per-fold views (createFromHandle :2145-2169
  -> buildFromParent data.hpp:995) already gather conditionally on
  options.leafCovariateColumns; the handle is the last unconditional
  site. buildFromParent already refuses designating ungathered columns.
- u_: column-major (model.hpp:218 resize; layout comment :538);
  per-member row gather at :384 and :401 strides numObservations_.
  Linear leaf only; GP leaf out of scope.

## Decision (handle gather)

Handles serve views whose specs are unknown at creation, so the gather
cannot be inferred: grow an optional leaf-covariate-columns argument on
handle creation (bridge + R surface), default empty -> no raw gathered.
Misuse is already refused downstream (ungathered -> rawColumn null ->
designation refused); internal callers (xbart folds) pass the union
their spec needs. Alternative - keep gather-all behind an opt-out -
rejected: it preserves the ~400MB waste (n=1e6, p=50) in the common
constant-leaf case.

## Constraints

- Bitwise at every pass: no draw-path arithmetic or RNG-order change
  (the u_ flip reorders storage, not the row fill arithmetic).
- The quantize core's observer is a template parameter, inlined; no
  std::function or virtual on the cell loop.
- data-store.md updated with each pass it describes (:84-85 and
  :181-189 for isView; the handle/view sections; the :358-361 residuals
  list shrinks per pass).
- No dbarts.h change. Pass order as in Steps (model.hpp last).

## Steps

1. Observer-parameterized quantize core: extract the cell loop;
   quantizeColumn = core + no-op observer; setColumnJournaled = core +
   journal observer.
2. isView split: isView becomes pure provenance; capability predicates
   match the two refusal granularities (no raw at all vs retains a
   re-quantize source); rewire both helpers, their call sites,
   test_data.cpp:72, and the doc language.
3. Conditional handle gather per the Decision block; a tinytest covers
   default-empty refusal and declared-column success.
4. u_ row-major flip: storage and both gathers (model.hpp).

## Verification

- Per pass: R CMD INSTALL --preclean .; tests/cpp from make clean; full
  tinytest suite (3267 expected); pause for orchestrator review before
  the next pass.
- Landing (orchestrator): three equivalence compares bitwise
  (equivalence-f494156, bcf-99205ee, multinomial-8c2b5fc); air format
  --check . for the R touches; no bench (no tree-loop hot path).
