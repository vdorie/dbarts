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

- Quantize twin: quantizeColumn ([src/bartcore/data.hpp:644-655](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L644-L655)) rides
  quantizeDenseInto; setColumnJournaled ([src/bartcore/data.hpp:1134-1162](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L1134-L1162))
  re-implements the same per-cell loop inline with journal bookkeeping.
  Sole caller [src/bartcore/sampler.hpp:1096](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/sampler.hpp#L1096) (SubsetUpdate), reached once per external
  iteration via updatePredictor/updatePredictorPerObservationJointly.
- isView: field [src/bartcore/data.hpp:214](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L214); set at [src/bartcore/data.hpp:720](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L720), [src/bartcore/data.hpp:768](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L768), [src/bartcore/data.hpp:1001](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L1001) (buildFromParent);
  provenance reader suppliedStandardization ([src/bartcore/data.hpp:342](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L342), plus the
  linear/GP reinitialize in model.hpp); capability readers
  refuseViewSampler ([src/R_interface_bartcore.cpp:1479-1487](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/R_interface_bartcore.cpp#L1479-L1487), ORs
  builtFromCsc, 5 mutation entry points) and refuseViewSamplerOnly
  ([src/R_interface_bartcore.cpp:1491-1496](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/R_interface_bartcore.cpp#L1491-L1496); setCutPoints/setState/installForests);
  [tests/cpp/test_data.cpp:72](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/tests/cpp/test_data.cpp#L72). Doc language: data-store.md "CodeBlock (train
  and test)" and "View semantics (buildFromParent)".
- Handle gather-all: [src/R_interface_bartcore.cpp:2026-2034](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/R_interface_bartcore.cpp#L2026-L2034) in
  bartcore_createDataHandle ([src/R_interface_bartcore.cpp:1997](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/R_interface_bartcore.cpp#L1997)). The sampler ctor
  ([src/bartcore/sampler.hpp:100-127](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/sampler.hpp#L100-L127)) and per-fold views (createFromHandle [src/R_interface_bartcore.cpp:2145-2169](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/R_interface_bartcore.cpp#L2145-L2169)
  -> buildFromParent [src/bartcore/data.hpp:995](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/data.hpp#L995)) already gather conditionally on
  options.leafCovariateColumns; the handle is the last unconditional
  site. buildFromParent already refuses designating ungathered columns.
- u_: column-major ([src/bartcore/model.hpp:218](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L218) resize; layout comment [src/bartcore/model.hpp:538](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L538));
  per-member row gather at [src/bartcore/model.hpp:384](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L384) and [src/bartcore/model.hpp:401](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L401) strides numObservations_.
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
- data-store.md updated with each pass it describes ([src/bartcore/model.hpp:84-85](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L84-L85) and
  [src/bartcore/model.hpp:181-189](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L181-L189) for isView; the handle/view sections; the [src/bartcore/model.hpp:358-361](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/src/bartcore/model.hpp#L358-L361) residuals
  list shrinks per pass).
- No dbarts.h change. Pass order as in Steps (model.hpp last).

## Steps

1. Observer-parameterized quantize core: extract the cell loop;
   quantizeColumn = core + no-op observer; setColumnJournaled = core +
   journal observer.
2. isView split: isView becomes pure provenance; capability predicates
   match the two refusal granularities (no raw at all vs retains a
   re-quantize source); rewire both helpers, their call sites,
   [tests/cpp/test_data.cpp:72](https://github.com/vdorie/dbarts/blob/4a52176065a83f73248b4334217683a8e9c10491/tests/cpp/test_data.cpp#L72), and the doc language.
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

## Landing

All four passes landed 2026-07-18, each a byte-identical refactor, in
order (serialized, one implementer, orchestrator review between passes):

1. afadbda - observer-parameterized quantize core: quantizeDenseObserved
   carries the per-cell loop with a compile-time observer; quantizeColumn
   is the no-op instantiation, setColumnJournaled supplies the journal
   observer (quarter-column fallback and missing tracking preserved).
2. bd7688e - isView split: isView is now pure provenance (gating only the
   inherited standardization constants); refusal reads two capability
   predicates, hasRequantizeSource and acceptsNewRawPredictors, and the
   bridge helpers renamed to refusePredictorMutation /
   refuseRequantizeWithoutSource. Every refusal fires where it did.
3. 26d1726 - conditional handle gather: bartcore_createDataHandle takes
   the 1-based leaf-covariate columns (default empty), owning raw for
   those alone instead of every column (~400MB of unread raw at n=1e6
   p=50 for a constant-leaf consumer); xbart declares its node prior's
   columns; an undeclared designation is refused with a message naming
   the handle declaration. New internal formal leafCovariateColumns on
   the unexported bartcoreDataHandle (no Rd). Bridge arity 2 -> 3; flag
   rchk on the next scheduled run. No dbarts.h change.
4. 170be0b - linear-leaf u_ row-major: matches accumulateNodeStatistics'
   per-observation gather; GP-leaf u_ left column-major (out of scope).

Landing gates (orchestrator, final combined state): install clean;
tests/cpp all pass from make clean; equivalence 26/26 bitwise
(equivalence-f494156) plus BCF 6-channel and multinomial 5-channel
bitwise; tinytest 3269/0; air format --check clean. Neutral throughout -
no baseline re-record. The data-store-residuals TODO entry is retired;
all four 2026-07-17 data-review residual findings are discharged.
