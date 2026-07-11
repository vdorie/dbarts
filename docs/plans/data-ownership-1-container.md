# data-ownership-1-container

agent: opus
rng: neutral - codes are integer comparison operands; quantization
  order, cut construction, and every draw byte-for-byte unchanged.
  u8/u16 is storage width only. Gate: equivalence 22/22 IDENTICAL at
  EVERY step; tests/cpp bit-exact; full tinytest.
window: plan 1 of docs/design/data-ownership.md (FROZEN 2026-07-11);
  plans 2-5 build on the container shape.
budget: ~1200-1800 lines across ~12 files, five staged commits; the
  u8 SIMD kernel family (step 4) is separable and AUTHORIZED TO SLIP
  to a plan-1b if it overruns ~1.5x budget (VD).
bench: VD granted the quiet machine. A fresh speed baseline
  bench-sampler-32fc7c8.csv is recorded before implementation starts
  (compare vs 31a4c01 first); step 1 gets a confirmatory compare and
  step 4 a full compare against it.

## Goal

Replace ColumnStore's borrowed `const double* x` + uniform u16 codes
with an owned, per-column-width, quantized container. The engine
borrows REAL(x) (or CSC slices) only DURING build/re-quantize calls
and retains no predictor raw (leaf-covariate columns keep an owned
gather, as today). Dense columns store codes at u8 (<= 254 cuts;
the default n.cuts = 100 path) or u16. Off-hot-path raw consumers
(setCutPoints, setData, setState cross-grid, getTrees saved-tree
replay) receive raw at call time. Bit-identical fits.

## Decisions folded in (VD, 2026-07-11 evening)

- ABI freeze LIFTED: stan4bart is the only dbarts.h consumer and we
  own it - "do whatever is clean and update the packages in
  lockstep." So getTrees GAINS an explicit training-replay data
  parameter (NULL = the retained creation spec) and setState may
  take raw for the cross-grid case, instead of hidden internal
  threading. Record the exact stan4bart lockstep delta in the
  landing note (the c-api-growth precedent). PROT_DATA stays
  regardless: it is the creation contract and the flat-C GC anchor
  (a C consumer holds only the external pointer). This also
  dissolves the draft's Q5 limitation - a flat-C caller passes
  current data explicitly.
- Mutation raw-source model (Q2): the sampler INSTALLS updated
  column vectors into its stored data object by reference and lets
  R handle GC - O(spine), not O(n x p) - once storage is
  column-oriented (plan 2 holds the frame). A perf-sensitive caller
  may update its own vector in place and re-install the same vector
  (the supported pattern; document it). PLAN-1 INTERIM while data@x
  is still a matrix: the R helper copy-modifies data@x for
  column-subset/per-observation updates - a temporary CoW cost,
  correctness unaffected. CONSEQUENCE TO FLAG at plan 3: the
  engine-owned mutable-raw flag's CoW rationale dissolves under
  reference-install; plan 3's convergence decides whether the
  mutable flag survives at all.
- State format (Q3): pre-release, compatibility breakable; keep the
  cutPoints-only format anyway on simplicity (codes never
  serialized; width recomputed from numCuts). No version bump.
- u8 SIMD x86 (Q1 addendum): this machine is arm64 and no x86 box
  exists yet (VD will set one up). Step 4 validates scalar + NEON
  here; x86 epi8 kernel sources may be written but the simd.c
  dispatch must NOT enable unvalidated u8 kernels - x86 entries
  stay routed to scalar until the component gate (u8 kernel ==
  scalar reference, bitwise) runs on VD's x86 machine; file the
  enable-x86-u8 follow-up on TODO at landing. (The existing u16
  x86 kernels are CRAN-proven and untouched.)

## Container design

Keep the ColumnStore name and struct-of-arrays layout (every
data.column(j)/codeAt/types[j]/numCuts[j] call site keeps its
shape). Changes:

1. Delete resident `const double* x` / `x_test`. build/buildMixed/
   buildTest take raw as a call argument, quantize, gather flagged
   columns, retain nothing. denseSourceColumn becomes build-local;
   mixedRawColumns/cscSlices borrows narrow to the call.
2. Per-column width: codeWidth[j] (1 or 2) computed from
   numCuts[j]/type (ordinal u8 iff numCuts <= 254; categorical u8
   iff inline-mask levels fit; sparse/CSC/densified stay u16 this
   pass); owned bytes in one codeBytes pool with codeByteOffset;
   typed columnU8/columnU16 accessors + width-generic codeAt for
   cold paths. naCode 0xFF / 0xFFFF. Width never serialized.
3. Leaf-covariate owned gather at the top-level store (reuse the
   view machinery's gatheredRaw* buffers; the top-level store is "a
   view of itself" for its leaf columns); rawColumn(j) serves owned
   memory after the borrow releases; regathered on setData.
4. CSC/sparse columns keep current machinery (plan 5); coexist.
5. Cut/level tables already owned; unchanged.

## Width dispatch

scan.hpp and the scalar partition helpers (partitionIndicesMIA/
ByMask/ByWideMask) become width-templated (one source each; each
instantiation a plain typed gather, no per-element branch).
tree.hpp partitionChildren's ordinal-dense arm gains the width
axis: u16 -> existing misc kernels; u8 -> step 4's epi8 family
(scalar until then). misc partition_body.c grows a u8 variant per
ISA; simd.c dispatch gains the (op, width) dimension
(kernel-vocabulary.md's planned table). Phase-1 measured no u8
SIMD throughput win, so step 4 exists to prevent a scalar-u8
regression on the default path, not to gain speed; the memory win
(u8 codes, ~8x under the double borrow) lands at step 3.

## Raw source per consumer

- setData: re-quantize the transient borrow of the new spec's @x;
  re-pin PROT_DATA to the new spec.
- setPredictor full-matrix: engine re-quantizes the transient
  borrow; rollback restores snapshotted codes (the data_.x = oldX
  write-back deletes); R side already re-assigns data@x.
- updatePredictor column-subset / per-observation partial: codes
  update from supplied values; write-through
  (sampler.hpp:858/861/1059, data.hpp:784/801 const_cast) DIES;
  the R helper keeps data@x current per the Q2 interim; the
  installed-mask contract unchanged.
- setCutPoints: R method passes current data@x; bridge re-quantizes
  the transient borrow (setCutPointsForColumn gains a raw param).
- setState: same-spec restores skip re-quantization (per-column
  cutPoints equality guard - codes already correct, no raw);
  cross-grid takes raw (R passes data@x; the dbarts.h signature
  change per the lifted freeze, or the retained spec - implementer
  picks the cleaner and records the lockstep delta).
- getTrees saved-tree replay: gains the explicit replay-data
  parameter (NULL = retained creation spec); the resident store.x
  read deletes. FlatNode stores resolved double cuts, so replay
  stays grid-independent; behavior unchanged.

## Steps

1. Container ownership at uniform u16: delete resident x/x_test,
   thread raw as call arguments, top-level leaf gather,
   buildFromParent adapted. tests/cpp: new-vs-old codes bitwise,
   view codes bitwise, leaf gather bitwise. Equivalence IDENTICAL;
   confirmatory bench-sampler compare vs bench-sampler-32fc7c8.csv.
2. Bridge + R rewire: PROT_PREDICTORS/TEST delete; getTrees replay
   parameter (+ dbarts.h signature update + stan4bart lockstep
   note); setState skip-guard + cross-grid raw; setCutPoints raw
   param; write-through killed; R-side data@x maintenance. Full
   tinytest (mutation, save/load, getTrees files named in the
   session draft).
3. Per-column width + typed accessors + codeBytes pool; width-
   templated scan/partition helpers; u8 columns on the SCALAR
   partition first. tests/cpp routing-equality: u8 column vs same
   data forced u16, bitwise identical partitions AND draws.
   Equivalence IDENTICAL.
4. u8 SIMD kernels (scalar + NEON validated here; x86 sources
   dispatch-gated until validated on VD's x86 machine) + (op,width)
   dispatch. Component: u8 kernel == scalar bitwise (per available
   ISA). bench-sampler compare vs bench-sampler-32fc7c8.csv: no u16
   regression; record the u8 result. SEPARABLE (slip = plan-1b).
5. Landing note in the design note; NEWS/Rd for the write-through
   death + the supported re-install pattern; file the x86-u8
   follow-up.

## Verification

- tests/cpp component gates above + existing suite (delete stale
  binaries; preclean installs - facade/struct changes).
- Full tinytest per step from a preclean install.
- Equivalence 22/22 IDENTICAL at every step (the frozen-default
  gate; any deviation = defect, stop).
- bench-sampler: fresh baseline recorded pre-implementation
  (32fc7c8); compares at steps 1 and 4 (quiet machine, granted).
- rchk when next scheduled (step 2 touches PROT machinery).

## Landing notes (steps 1-2)

Container: ColumnStore dropped the resident x/x_test borrows. build/
buildMixed/buildFromCsc take the raw as a call argument and quantize into
owned codes; a per-store gather designation (leaf covariates for a sampler,
every column for a data handle) has build own raw copies of those columns so
rawColumn serves them after the borrow releases. Test-side raw is owned in
full (ownedTestValues) so cut changes re-quantize the test codes without a
retained borrow; rawTestColumn serves it. Mixed/CSC builds keep their
mixedRawColumns/cscSlices (reachable through PROT_DATA); a mixed column's
rawColumn still returns its retained dense slice, so no gather is needed
there. Views set isView, which now discriminates the mutation/re-quantize
refusals in place of the deleted x==NULL test.

setState cross-grid decision (VD's "pick the cleaner"): the retained-spec
route, NOT a dbarts.h signature change. The same-spec continuation contract
(the only supported restore) matches the live grid column for column, so a
per-column skip guard (state.cutPoints[j] == live) re-quantizes nothing and
needs no raw. Only an off-contract cross-grid restore needs raw; the engine's
setState gained an internal currentPredictors argument the callers fill - the
R method from sampler$data@x, the flat-C path from the retained creation
spec's @x - so dbarts_sampler_setState keeps its 2-argument signature. A
flat-C caller that mutated predictors then restores a cross-grid state sees
the pre-mutation spec (documented plan-1 limitation, closed by plan 3's
extract verb).

stan4bart lockstep delta (do whatever is clean, update in lockstep):
- dbarts_sampler_getTrees gained a trailing `const double* trainingData`
  argument (column-major numObservations x numPredictors, or NULL to replay
  the creation spec). Every stan4bart call site adds one trailing argument:
  the current predictor matrix pointer, or NULL. No other behavior changes;
  live-tree getTrees ignores it.
- No other dbarts.h signature changed. dbarts_sampler_setState,
  setPredictor, updatePredictor, setTestPredictors are all unchanged.
- The .Call bridge entry points DID grow arguments (getTrees +trainingData,
  setState/setCutPoints +currentPredictors), but those are internal R-.Call
  symbols, not the LinkingTo C ABI, and their R callers are updated here.

Write-through death (step 2): the const_cast writers (data.hpp setColumns/
setCell), the sampler rollback write-back and the data_.x = oldX line, and
commitObservation's write-through are gone. Predictor mutation updates only
the owned codes (and a leaf column's gathered raw); rollback restores the
snapshotted codes and gathered raw. The R layer maintains sampler$data@x
itself for column-subset and per-observation updates (copy-modify), the
plan-1 interim while data@x is a matrix - a temporary copy-on-write cost on
the per-iteration updateX/IRT path, correctness unaffected (plan 3 removes
it). CONSEQUENCE, flagged: a shallow copy no longer shares a predictor
mutation with the original (write-through backed that sharing); the R-side
copy-modify diverges them per R copy-on-write. test-sampler-predictors.R's
shallow/deep-copy assertion was updated to the new model.

PROT slots: PROT_PREDICTORS and PROT_TEST_PREDICTORS are deleted (the engine
borrows no predictor/test matrix past a call). PROT_DATA stays as the
creation contract and the flat-C GC anchor; PROT_RESPONSE/OFFSET/WEIGHTS/
TEST_OFFSET stay (those borrows persist).
