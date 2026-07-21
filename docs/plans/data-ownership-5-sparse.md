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
  CSC, data-ownership.md). No plan follows it.
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
- S-CAT (data-ownership-2-ingestion.md): plan 2 landed sparseFactor as
  R surface only and refused it before the engine. This plan lifts exactly that
  refusal and supplies the storage + kernel it was deferred against.

## Context (what already holds, read in code)

- The container crosses {ordinal, categorical} x {dense, CSC}; 3 of 4 cells
  ship, sparse categorical is the gap (data-ownership.md). ColumnType is
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

Refinement forced on the design: data-ownership.md reads "the implicit
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
   data-ownership.md "Implementation record" item 5 says "test-side sparse (with sparse-extensions)":
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
- data-ownership.md "the implicit zero as the reference level" reads as
  code 0; the bitwise gate forces zeroCode = the reference level's actual
  level-order code. The landing note records the refinement.
- data-ownership.md "Implementation record" item 5 bundles "test-side sparse (with
  sparse-extensions)"; this plan satisfies it by densification and leaves
  resident sparse x.test deferred in sparse-extensions.md (decision 2) - not a
  contradiction, but the wording could be read as requiring resident sparse
  testCodes, which it does not.

## Landing notes

Commit 1 = 1d02907. Engine kernel. partitionIndicesSparseByMask (inline,
<= 63 levels) and partitionIndicesSparseByWideMask (pooled, > 63 levels)
added in tree.hpp next to the dense siblings and the sparse MIA partition,
reading through SparseColumnData::at; partitionChildren's categorical
branch gains a columnIsSparse sub-dispatch ahead of the dense arm.
buildCutsForColumn takes K from the bridge-supplied metadata for a
CSC-categorical column instead of dereferencing a null column;
quantizeCscColumn seeds zeroCode from the per-column REFERENCE CODE (see
the refinement below), not codeFor(j, 0.0), and scatters/densifies the
implicit rows to it. SamplerOptions.cscCategoryCounts/.cscReferenceCodes
carry the per-column K and reference code down from the bridge. New
tests/cpp: a rank-tier and a densified-tier sparse-categorical column
store test plus an end-to-end fit, both bitwise vs. a dense factor of the
same values. Gate: equivalence 22/22 identical; full tinytest 2786, no
regen.

Commit 2 = 54460f8. Ingestion. sparseColumnSlices (mixedMatrix.R) stops
refusing sparseFactor and emits its CSC slices - i, and values - 1 as
0-based level-order codes - flagging the column categorical and carrying
its reference code (0-based, level order) and K; assembleMixedMatrix
threads these into sparseReference/sparseCategoryCount alongside the map.
parseData (R_interface_bartcore.cpp) resolves this per-column metadata and
threads it through SamplerOptions; validateCategoricalPredictors skips
CSC-backed categorical columns (the engine, not the R-side scan, now owns
their validity); the x/y-path S-CAT refusals (data.R, utility.R) drop,
the formula-path pre-scan does not - a bare S4 column still cannot survive
model.frame, so sparseFactor stays x/y-only by construction, not by an
arbitrary restriction. New inst/tinytest/test-sparse-factor.R carries the
R-level BITWISE gate: sparse-vs-dense identical() on sigma and train fits
across three tiers - rank (dominant reference, stored fraction under the
0.2 densification threshold), densified (rare reference, over threshold),
and pooled (K = 80, > 63 levels, exercising the wide-mask kernel). Gate:
equivalence 22/22; full tinytest 2798, no regen.

Commit 3 = 34e87a6. Test-side densification + a materializer fix needed to
support it. as.matrix.dbartsMixedMatrix's sparse branch previously
converted a sparse column's dgCMatrix slice straight to dense
(as.matrix(x$sparse[...])), which fills every implicit row with numeric
zero - correct for a sparse ORDINAL column, wrong for a sparse
CATEGORICAL one whose implicit rows carry the reference level's code
(generally nonzero, and possibly coinciding with an explicit entry's own
code, which can legitimately be 0). Fixed to walk the CSC structure
directly (p/i/x on x$sparse) per column: fill every row with the
reference code FIRST, then scatter only the explicitly-stored entries
over their known row positions - a value-based zero-check would have
misclassified a legitimate explicit 0. This also fixes
extract(sampler, "predictors") for the same reason, since it shares the
materializer. A sparseFactor x.test column densifies to dense categorical
test codes over the TRAINING level table via densifySparseFactorColumns,
running ahead of the model-frame replay so prediction takes the existing
dense test-code path unchanged; predictions on a sparse-cat test frame
verified bitwise against its densified twin. The formula-path refusal
message was reworded to "must be supplied through the x/y interface"
(previously implied but not stated). Gate: equivalence 22/22; full
tinytest 2803, no regen.

Commit 4 = 5461f41. Mixed-flavor unification (decision 1 / Q2), WITH A
DEVIATION worth recording precisely. The plan's step-4 text called for the
engine to gather leaf-covariate raw owned instead of borrowing resident
dense slices; en route it surfaced that tests/cpp's testMixedColumnStore
pins rawColumn(0) == the source block by POINTER IDENTITY - a borrow
contract the test enforces deliberately, not an oversight - so the engine
side was left BYTE-FOR-BYTE UNCHANGED rather than reworked to own the
gather. Instead, ownership of the assembled dense block moved to the
BRIDGE: DataHandle and BartcoreHolder each gained an ownedMixedDense
member (mirroring ownedResponse/ownedTreatment), and parseData's
transient denseAssembly is moved into it after construction - a
std::vector move preserves the buffer address, so the store's borrowed
dense slices stay valid for the sampler's/handle's lifetime without
reworking the engine's borrow contract. On the R side,
assembleMixedMatrix now stores dense columns as a per-column list for
BOTH flavors (dense-only and mixed), so the two flavors are told apart by
whether $sparse is populated, not by the type of $dense; parseData
assembles the transient denseAssembly for the mixed flavor exactly as
plan 2 already did for the dense-only one. No R-side resident dense block
remains in either flavor. Gate: existing testMixedColumnStore /
EndToEnd / LinearLeaves unchanged and passing; equivalence 22/22; full
tinytest 2803; tests/cpp component tests unchanged.

Commit 5 = this commit (docs). man/sparseFactor.Rd: sparseFactor is now
accepted as a predictor through the x/y interface (a frame may mix dense
and sparse, ordinal and categorical, columns freely); the formula
interface still refuses it (a bare S4 column cannot survive
model.frame); the reference argument's doc records Q3's resolution - the
most common level as reference maximizes the memory win, correctness
holds under any choice. man/dbartsData.Rd: the data-frame ingestion table
updated to reflect sparseFactor acceptance and x/y-only entry; a note that
a sparseFactor test column densifies to codes over the training level
table before prediction. docs/design/kernel-vocabulary.md addition 2:
updated to record that categorical membership partitions, dense AND
sparse, landed engine-side in tree.hpp rather than as misc.a kernels (the
plan's drift item 1) - the planned misc_partitionRangeCat/
misc_partitionIndicesCat signatures never landed.
docs/design/data-ownership.md: item 5 of the implementation split marked
LANDED through 5461f41, including the zeroCode refinement (drift item 2)
below. Gate: tools::checkRd clean on both Rd files; full tinytest 2803, no
regen.

Refinement recorded (drift item 2): data-ownership.md's column-kinds sketch
read "the implicit zero as the reference level," suggesting numeric code 0.
The bitwise gate forced a correction: a dense factor codes its reference
level by LEVEL ORDER like every other level, so the CSC column's zeroCode
must be the reference level's ACTUAL level-order code (generally not 0) for
at(i) to equal the dense code at every row. Storage carries this as the
bridge-supplied reference code (SamplerOptions.cscReferenceCodes), and
quantizeCscColumn seeds zeroCode from it rather than from codeFor(j, 0.0).
data-ownership.md's sketch is corrected accordingly, not merely annotated.

Open-question resolutions, all landed as recommended:

- Q1 (pooled > 63-level sparse categorical): INCLUDED. The wide-mask
  sibling (partitionIndicesSparseByWideMask, commit 1) and a K = 80 pooled
  tier in test-sparse-factor.R (commit 2) cover it; no >63-level refusal
  was added.
- Q2 (resident-dense-block revisit here vs. a follow-up): LANDED HERE, as
  commit 4, with the pointer-identity deviation recorded above.
- Q3 (reference-level default): KEPT levels[1] as the constructor default;
  documented (sparseFactor.Rd) that the most common level should be chosen
  as the reference to realize the memory win, correctness holding under
  any choice.
- Q4 (resident sparse x.test): DEFERRED, per sparse-extensions.md and the
  ordinal precedent - sparse-categorical x.test densifies and predicts
  (commit 3). This densification interim is now explicitly SUPERSEDED by
  the queued backlog item docs/plans/test-data-parity.md (VD directive,
  2026-07-13: training and test sides converge to one data model), not
  just deferred a second time - a future test-data-parity plan replaces
  the densify-at-ingestion step with resident sparse testCodes end to end,
  training and test sides sharing one data model.

Memory demonstration (verification item, computed + measured, arm64,
deterministic seed 20260714): n = 1e5, one sparseFactor with K = 200
levels and a dominant (95%) reference level (f = 0.0500 non-reference,
5005 stored entries), plus two dense numeric columns - the bairrtt
item-design shape - ingested two ways: as a data frame carrying the
sparseFactor, and as the same values with a dense factor column. The
sparse ingestion path was confirmed taken (the dbartsMixedMatrix
container's $sparse slot populated for the sparseFactor frame, NULL for
the dense frame) and both fits ran to completion with finite draws;
sigma draws were additionally checked bitwise-identical between the two
ingestions (the gate this plan built, spot-checked here on a larger K
than the tinytest tiers exercise). Container-level object.size (R
containers, which retain the whole source frame as data@x's GC anchor,
not just engine storage) came to 2,479,352 bytes sparse vs. 2,830,744
bytes dense, a ~12.4% reduction - the weaker of the two signals, since
most of both totals is the retained frame rather than engine state. The
arithmetic column cost from docs/design/sparse-columns.md's formula is
the meaningful number: dense costs 2 bytes/row (200,000 bytes for this
column), rank-bitmap sparse costs 0.1875 + 2f bytes/row = 0.2876 bytes/row
at f = 0.05 (28,760 bytes) - a 6.95x shrink on the categorical column
alone, consistent with the plan's stated 7-10x range at f <= 0.05 (this
run sits just under the range's low end because f = 0.05 is exactly the
range's upper, least-sparse boundary).

Bench (quiet arm64 box): bench-sampler compare vs
bench-sampler-4008675.csv reported "OK: no metric regressed more than
5%" - zero flags, every metric at or below baseline, so the dense paths
carry no cost from the sparse arms' guards. Sampling parity on the
memory-demonstration shape (n = 1e5, K = 200, f = 0.05, single chain,
m = 75): sparse ingestion 13.9 ms/iteration vs dense 14.4 (ratio 0.964,
parity within noise). All Verification items are confirmed; the plan is
CLOSED, and with it the data-ownership program (plans 1-5) is COMPLETE.
