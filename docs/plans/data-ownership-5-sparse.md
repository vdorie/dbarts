# data-ownership-5-sparse

agent: opus (the partition seam, buildMixed, the bridge parse path, and the
  sparseFactor R surface are one contract; one owner keeps the CSC-categorical
  code assignment coherent across all three).
rng: neutral - a sparse-categorical column carries the SAME level-order codes
  and category count as a dense factor of the same values, and the categorical
  split proposal draws over the rule-derived reachable mask, never over
  observed codes (collectAvailableVariables tree.hpp:481-490,
  categoryDirectionsForPattern model.hpp:1507), so storage cannot move a draw.
  Every existing path is guarded: the partitionChildren categorical branch
  gains a columnIsSparse sub-dispatch, buildMixed's CSC branch a categorical
  variant, both inert for today's columns. Gate: equivalence 22/22 IDENTICAL
  vs equivalence-ac6ec2c.rds at EVERY commit + a NEW sparse-cat-vs-dense-cat
  bitwise component test (the ordinal landing's device); full tinytest 2786
  from a preclean install, new tests ADD, no snapshot regen.
window: plan 5 of docs/design/data-ownership.md (FROZEN 2026-07-11); builds on
  plan-1 container, plan-2 ingestion (6c90507), plan-3 mutation (a141338),
  plan-4 views (49bb1d4). Closes the container's fourth kind (categorical x
  CSC, data-ownership.md:73-87). No plan follows it.
budget: ~700-1100 lines across ~14 files (src/bartcore/tree.hpp, data.hpp,
  src/R_interface_bartcore.cpp, R/mixedMatrix.R, R/data.R, R/utility.R,
  R/A_class.R already carries the class, tests/cpp/test_model.cpp, tinytest
  test-sparse-factor.R, man/*.Rd, the design landing note). Header edits ->
  --preclean; delete stale tests/cpp binaries.
bench: memory-first (the kernel's point). One container-bytes comparison
  (deterministic, any box) + one bench-sampler compare vs
  bench-sampler-4008675.csv on the quiet arm64 box (dense paths must not
  regress; a sparse-cat arm at parity). No x86 SIMD gate: the membership
  partition is engine-side scalar (dense categorical membership already is);
  the dense SIMD kernels are untouched, covered by the standing cross-ISA gate.

## Goal

A high-cardinality unordered factor stored sparsely (CSC over level codes, the
reference level implicit) becomes a first-class predictor: sparseFactor() is
ACCEPTED where it is refused today (sparseColumnSlices, mixedMatrix.R:36), so a
frame mixing dense numerics, dense factors, sparse ordinal columns, and sparse
factors ingests in one call. The win is memory on one-hot / high-cardinality
factor designs (the IRT / bairrtt item-design shape): a dense categorical
column costs 2 bytes/row (n*p u16 codes), a rank-bitmap sparse one ~0.19 + 2f
bytes/row (docs/design/sparse-columns.md results), a 7-10x column shrink at
f <= 0.05 while sampling stays at dense parity. User-visible contract:
sparseFactor accepted as a data-frame column via the x/y interface; sparse
categorical in x.test densifies and predicts; the raw-x MUTATION surface stays
refused for sparse sources (sparse designs fix at creation, the standing rule),
as does rbart_vi's R-loop path; linear/gp leaves stay refused on sparse-backed
columns (plan-4 precedent, R_interface_bartcore.cpp:463).

## Binding contracts inherited (plans 1-4, do not reopen)

- data@x is the COLLECTED RAW SOURCE. A sparseFactor rides data@x by R
  reference as the GC anchor (isSparseDataFrameColumn recognizes it,
  mixedMatrix.R:26); it is never a view of engine state. Per-draw mutation of a
  sparse column stays REFUSED (installPredictorColumns comment, mixedMatrix.R:
  278-280): sharing/swapping whole data replaces the source, in-place mutation
  does not.
- The engine keeps NO predictor raw except the leaf-covariate gather
  (data.hpp:193). Do not add engine-owned mutable raw; the mutable-raw flag was
  killed, never built.
- Views DENSIFY: buildFromParent (data.hpp:742) gathers through the
  storage-aware codeAt (data.hpp:828,837,924), so a sparse-categorical column
  folds into a dense categorical child for free - no view-side sparse code.
- extract(sampler, "predictors") materializes the NUMERIC code matrix over
  data@x (factors as level-order codes); a sparse-categorical column densifies
  to its codes there (implicit rows = the reference code). No engine reroute.
- dbarts.h is UNCHANGED: the CSC-categorical path is internal (engine .hpp +
  SamplerOptions, freely extensible + R surface). No stan4bart lockstep (the
  plan-1/2/4 precedent). PROT_DATA stays the creation contract and GC anchor.
- S-CAT (data-ownership-2-ingestion.md:63-70): plan 2 landed sparseFactor as
  R surface only and refused it before the engine. This plan lifts exactly that
  refusal and supplies the storage + kernel it was deferred against.

## Context (what already holds, read in code)

- The container crosses {ordinal, categorical} x {dense, CSC}; 3 of 4 cells
  ship, sparse categorical is the gap (data-ownership.md:73-87). ColumnType is
  the kind (data.hpp:83); a rank column lives in a SparseColumnData slot
  (data.hpp:95-108: bits, wordRanks, nzCodes, zeroCode, at(i)), keyed by
  sparseSlot (data.hpp:142); columnIsSparse/sparseColumn (data.hpp:917-922)
  select it; codeAt (data.hpp:924) is storage-aware already.
- buildMixed dispatches per column (data.hpp:640-690): dense-backed columns
  slice denseValues, CSC columns build a rank bitmap or densify at
  sparseDensityThreshold=0.2 (data.hpp:89,653). buildFromCsc is buildMixed with
  every column CSC (data.hpp:702). The scatter/quantize is storage-aware
  (quantizeCscColumn, data.hpp:498).
- The categorical partition is ENGINE-SIDE and scalar, NOT a misc.a kernel:
  partitionIndicesByMask (inline <=63 levels, tree.hpp:533) and
  partitionIndicesByWideMask (pooled >63, tree.hpp:554), both reading a dense
  const xint_t* column. partitionChildren dispatches
  categorical -> sparse-ordinal -> dense-ordinal (tree.hpp:637,650,665); the
  categorical branch reads data.column(variable) (dense), with NO sparse arm.
  The sparse-ordinal siblings that DO read through SparseColumnData::at already
  exist: misc_partitionIndicesSparse (tree.hpp:659) and the engine-side
  partitionIndicesSparseMIA (tree.hpp:575).
- Categorical splits are proposed via the MH move path over the rule-derived
  reachable mask (model.hpp:1507,1586-1611), never scanned (grow.hpp:55-59,100)
  and never read from observed codes; availability is mask-derived
  (tree.hpp:481-490). So codes are consumed ONLY at partition and (indirectly,
  via partitioned index ranges) at suffstat - the basis of the bitwise gate.
- The refusal choke points, all lifted here: sparseColumnSlices
  (mixedMatrix.R:36, the single point every sparse column passes), the
  formula-path pre-scans (data.R:410-414,680-681), and the bridge's
  all-ordinal-sparse guard (R_interface_bartcore.cpp:459-464,884).
- sparseFactor exists (R/sparseFactor.R constructor; A_class.R:513-554 class +
  validity): slots i (0-based ascending rows), values (1-based level codes),
  levels, reference (implicit level), length. It forbids NA (sparseFactor.R:47),
  so a sparse-categorical column is never MIA - no missing-direction sibling.

## THE GATE: bitwise sparse-vs-dense (PROMISABLE)

For any dataset expressible both ways, a sparse-categorical-ingested sampler
draws BITWISE-IDENTICALLY to the same data ingested dense. Why it holds:

- IDENTICAL CODES. A dense factor codes level c as c-1 (level order). A
  sparseFactor over the same data with the same `levels` also codes each stored
  entry in level order; the implicit rows carry the reference level, whose code
  is its OWN level-order index. So the CSC column's zeroCode is the reference
  level's level-order code (generally NOT numeric 0), and at(i) == the dense
  code[i] for every row. The reference choice is a STORAGE decision (which
  level is implicit), orthogonal to the codes.
- IDENTICAL CATEGORY COUNT. numCuts[j] = number of levels = K, the dense factor's
  count; the reachable-mask arithmetic (tree.hpp:481-490) is thus identical.
- IDENTICAL RNG STREAM. The proposal draws over the reachable mask, not codes;
  suffstats accumulate over partitioned index RANGES (y, not codes). Only the
  partition reads codes, and identical codes give an identical partition.
- The densified tier IS the gate device: a sparse-categorical column above the
  density threshold scatters to dense categorical codes (reference code in the
  implicit rows) and runs the existing dense categorical kernel bitwise. The
  rank tier reproduces the same partition through at(i). This is the ordinal
  landing's testSparseColumnStore / testSparseEndToEnd device (test_model.cpp)
  extended to the categorical kind.

Refinement forced on the design: data-ownership.md:82-85 reads "the implicit
zero as the reference level," which suggests code 0; the bitwise gate requires
the reference's ACTUAL level-order code as zeroCode. The store carries it per
CSC-categorical column (bridge metadata); zeroCode is set from it, not from
codeFor(j, 0.0).

## Scope decisions (recorded, decided from the design docs)

1. RESIDENT DENSE BLOCK revisit: IN SCOPE. Plan 2 landing (commit 2, c82f954)
   deferred it here explicitly ("plan 5 revisits when the CSC-categorical
   kernel lands"): the sparse-bearing mixed flavor keeps a resident cbound
   dense block because buildMixed borrows lifetime slices of it. With the mixed
   build already open for CSC-categorical, step 4 unifies the flavor onto the
   plan-2 dense per-column list + transient block, gathering leaf-covariate raw
   owned (as build already does) so no lifetime borrow remains. Byte-identical
   codes; memory hygiene only. (Open Q2 offers to spin it out if it bloats.)
2. TEST-SIDE SPARSE: resident sparse test STORAGE stays DEFERRED, consumer-
   gated per sparse-extensions.md; sparse-categorical x.test DENSIFIES (the
   container's as.matrix -> dense categorical test codes over training levels),
   which reuses the whole existing dense test-code path and predicts correctly.
   data-ownership.md:180 item 5 says "test-side sparse (with sparse-extensions)":
   read as "sparse categorical must WORK on the test side," satisfied by
   densification; resident sparse testCodes waits for a named memory-bound test
   workload (matching sparse-columns.md "x.test stays a dense matrix" for the
   ordinal kind). See Open Q4.
3. POOLED (>63-level) sparse categorical: IN SCOPE (Open Q1). High cardinality
   is the stated motivation; the sparse wide-mask sibling is a mechanical
   mirror of partitionIndicesByWideMask over at(i).

## Relationship to sparse-extensions.md

This plan neither supersedes nor is superseded by sparse-extensions.md. That
doc is the standing consumer-gated meta-plan for the FOUR deferred sparse
extensions (in-place nonzero mutation, resident sparse x.test, a streaming
range kernel, dense-column mutation in mixed stores); the sparse-CATEGORICAL
kernel is not among them - it is data-ownership item 5, a new container cell,
and this is its plan. It CONSUMES sparse-extensions' dense-equality
component-test device (dgCMatrix-vs-dense bitwise identity) and COORDINATES on
test-side sparse: resident sparse x.test stays deferred THERE (decision 2),
untouched by this plan. The other three extensions stay closed and
consumer-gated.

## Constraints

- Draw-neutral on every existing path: the categorical branch's sparse arm and
  buildMixed's categorical-CSC arm are reached only for a CSC-categorical
  column, which no current input produces. Equivalence 22/22 is the gate; any
  deviation on a dense/ordinal-sparse path = defect, stop.
- Single-writer, creation-fixed: sparse sources refuse the raw-x mutation +
  re-quantize surface exactly as today; do not weaken it, do not add
  engine-owned mutable raw.
- Fast over safe in C/C++: the sparse membership partition is a two-pointer
  swap over at(i), no allocation, mirroring the dense sibling.
- OUT of scope: resident sparse x.test (deferred, decision 2); in-place sparse
  mutation (sparse-extensions, no consumer); rbart_vi sparse (its R loop
  predicts over data@x); linear/gp leaves on sparse-backed columns (refused,
  plan-4 precedent). dbarts.h unchanged.

## Steps

1. Sparse-categorical membership partition + code assignment (engine). Add
   partitionIndicesSparseByMask (inline) and ...ByWideMask (pooled) in tree.hpp
   next to the dense siblings (tree.hpp:533,554) and the sparse MIA
   (tree.hpp:575), reading through SparseColumnData::at. partitionChildren's
   categorical branch (tree.hpp:637) gains a columnIsSparse sub-dispatch. Teach
   buildCutsForColumn (data.hpp:392: today derefs a null column for CSC
   categorical) to take K from metadata for a CSC-categorical column, and
   quantizeCscColumn (data.hpp:498-517) to seed zeroCode from the reference code
   (not codeFor(j,0.0)) and scatter/densify the implicit rows to it. Store the
   per-column reference code (a small vector, or reuse the built zeroCode set
   from SamplerOptions). Gate: NEW tests/cpp - a rank AND a densified
   sparse-categorical column bin bitwise vs a dense factor of the same values
   (codes, partition membership, an end-to-end fit draw); equivalence 22/22;
   full tinytest 2786 no regen. Abort: any dense/ordinal-sparse draw moves - the
   sub-dispatch leaked.
2. Ingest sparseFactor (R + bridge). sparseColumnSlices (mixedMatrix.R:33) stops
   refusing sparseFactor and emits its CSC slices (i, values-1 as 0-based
   codes), flagging the column categorical and carrying its reference code + K;
   extend the mixed container with a per-sparse-column reference vector; drop the
   x/y-path S-CAT refusals (data.R, utility.R:229,284). parseData
   (R_interface_bartcore.cpp:459-464) allows categorical + CSC when the
   reference metadata is present and passes reference code + K through
   SamplerOptions. Keep the formula-path S4-column limitation (plan-2: bare S4
   terms die in model.frame; sparseFactor enters via x/y like sparseVector).
   Gate: NEW/extended tinytest (ingest, fit, save/load via getPointer
   re-creation, refusal messages narrowed); equivalence 22/22; full tinytest.
   Abort: a sparseFactor reaches the engine on a path lacking reference metadata.
3. Sparse-categorical x.test densification (R). A sparseFactor test column
   coerces to dense categorical test codes over the TRAINING level table
   (validateXTest / the model-matrix replay), so predict runs the existing dense
   test path; no resident sparse testCodes (decision 2). Gate: tinytest -
   predictions on a sparse-cat test frame match predictions on the densified
   frame bitwise; full tinytest. Abort: test codes disagree with training
   levels.
4. Unify the mixed flavor onto a transient dense block (engine + R,
   decision 1). Convert the sparse-bearing mixed flavor to the plan-2 dense
   per-column list; parseData assembles the transient block; buildMixed gathers
   leaf-covariate raw owned instead of borrowing resident slices. Gate: existing
   testMixedColumnStore/EndToEnd/LinearLeaves bitwise unchanged; equivalence
   22/22; full tinytest. Abort: any code differs from the resident-block
   reference - the transient assembly is wrong. (Spin to a follow-up per Q2.)
5. Docs + landing. sparseFactor.Rd (the memory-win reference guidance),
   dbartsData.Rd (sparse-categorical accepted), kernel-vocabulary.md addition 2
   (record the engine-side realization; see drift below), the design landing
   note (steps 1-2 style: the level-order-code refinement, the bitwise gate,
   decisions 1-3). Gate: R CMD check man; full tinytest.

## Verification

- Full tinytest 2786 per commit from a preclean install (--preclean after the
  header/bridge edits; delete stale tests/cpp binaries).
- Equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at EVERY commit (the
  frozen-default gate; any deviation is a defect, stop). New sparse-cat tests
  ADD to the suite; no existing snapshot regenerates.
- tests/cpp: rank-tier and densified-tier sparse-categorical bitwise vs the
  dense factor (codes, membership, fit draw); the pooled variant if Q1 = yes.
  If a new test draws from the shared rngState, snapshot/restore it around its
  own draws (the plan-4 step-3 finding, test_data.cpp precedent).
- Memory demonstration: a container-bytes comparison, dense vs sparse-cat, on a
  high-cardinality design (a factor with K ~ 100-500, f ~ 0.05, n ~ 1e5, the
  bairrtt item shape) showing the 7-10x column shrink. Deterministic, any box.
- One bench-sampler compare vs bench-sampler-4008675.csv on the quiet arm64 box
  (dense + the sparse-cat arm; no dense-path regression, sparse-cat at parity).
  Sub-ms setPredictor metrics carry 6-8% noise; re-run a single flag.
- rchk on the next scheduled run (steps 2, 4 touch the bridge).
- dbarts.h unchanged; no stan4bart lockstep (bridge growth is internal .Call).

## Open questions for VD

- Q1 (pooled >63-level sparse categorical in v1). RECOMMEND include it: the
  sparse wide-mask sibling is a ~15-line mirror of partitionIndicesByWideMask
  over at(i), and high cardinality is the whole motivation - refusing >63-level
  factors would gut the feature. What would change it: shipping inline-only
  (<=63) first is faster and a >63 sparseFactor is refusable with a clear
  message, but it excludes exactly the designs the kernel targets.
- Q2 (resident-dense-block revisit here or as a follow-up). RECOMMEND land it as
  step 4: plan 2 assigned it here and the mixed build is already open, so
  unifying the flavor completes the columnar-source memory story in one place.
  What would change it: spinning it to a tiny standalone follow-up ships the
  kernel (steps 1-3) sooner, at the cost of leaving mixed frames materializing a
  resident dense block until then.
- Q3 (reference-level default). RECOMMEND keep the constructor's levels[1]
  default and DOCUMENT that the reference should be the MOST COMMON level to
  realize the memory win (storage holds the non-reference rows). What would
  change it: auto-picking the modal level would maximize sparsity but needs a
  construction-time scan of x and would surprise users expecting baseline-
  contrast semantics; codes and draws are correct under any reference.
- Q4 (resident sparse x.test). RECOMMEND defer, per sparse-extensions.md and the
  ordinal precedent: sparse-categorical x.test densifies and predicts, and
  resident sparse testCodes waits for a named memory-bound test workload. What
  would change it: a huge sparse test set where the densified test matrix is the
  memory ceiling would justify mirroring the training-side kernel onto testCodes.

## Design-vs-code drift to record in the landing note

- kernel-vocabulary.md planned addition 2 proposes categorical membership as
  misc.a kernels (misc_partitionRangeCat / misc_partitionIndicesCat). Reality:
  dense categorical membership landed ENGINE-SIDE (tree.hpp partitionIndicesBy
  Mask/ByWideMask), and this plan's sparse-categorical membership follows suit -
  no misc.a entry lands. Update addition 2 to record the engine-side
  realization.
- data-ownership.md:82-85 "the implicit zero as the reference level" reads as
  code 0; the bitwise gate forces zeroCode = the reference level's actual
  level-order code. The landing note records the refinement.
- data-ownership.md:180 item 5 bundles "test-side sparse (with
  sparse-extensions)"; this plan satisfies it by densification and leaves
  resident sparse x.test deferred in sparse-extensions.md (decision 2) - not a
  contradiction, but the wording could be read as requiring resident sparse
  testCodes, which it does not.
