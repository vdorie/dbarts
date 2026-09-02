# test-data-parity

agent: opus (the typed test store, the descent-through-codeAt change, the
  bridge test-container parse, and the R test-ingestion surface are one
  contract; one owner keeps test codes bitwise-equal to the densified path
  at every seam).
rng: neutral, EVERY step. Test data never enters a draw: tree proposals and
  suffstats read training data only (chain.hpp draw loop), test rows are
  routed down already-grown trees and their fits RECORDED after the draws
  ([[chain.hpp:2396-2410@4a521760]]). Storage and codes may move; no draw can. Gate:
  equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds (guards training
  draws) + the test-fit BITWISE oracle below (guards test values) + full
  tinytest 2803 from a preclean install, new tests ADD, no snapshot regen.
window: supersedes data-ownership-5-sparse.md decision 2 / Q4 (test-side
  densification interim) and data-ownership.md "Implementation record" item 5's "test-side sparse
  by densification"; builds on the plans-1..5 container (5461f41). Closes the
  test half of the owned data model. No plan follows it.
budget: ~800-1300 lines across ~15 files (src/bartcore/data.hpp, chain.hpp,
  tree.hpp, sampler.hpp, facade.hpp, src/R_interface_bartcore.cpp,
  R/data.R, R/mixedMatrix.R, R/utility.R, R/A_class.R, R/bartcore.R,
  tests/cpp/test_data.cpp, inst/tinytest test-sparse-factor.R + a test-side
  file, man/*.Rd, the design landing note). Header/facade edits -> --preclean;
  delete stale tests/cpp binaries.

## Goal

After it lands, test data is the same owned, typed, quantized, per-column
container the training side became: a frame / columnar source, per column
{dense, sparse} x {ordinal, categorical}, resident end to end. A
sparse-categorical or sparse-ordinal x.test is stored sparsely (rank bitmap),
not densified at ingestion; a mixed test frame ingests in one call; the engine
descends test rows through storage-aware codeAt, so the row-major dense
testCodes buffer and the dense ownedTestValues copy both retire. Training and
test sides share ONE data model.

## Binding contracts inherited (plans 1-5 + dbarts.h freeze; do not reopen)

- dbarts.h is FROZEN and its test entries are DENSE double* by signature:
  dbarts_sampler_setTestPredictors ([[dbarts.h:185@4a521760]]), _predict (197),
  _numTestObservations (257), the results `test` slab (87). The dense test
  path MUST keep accepting a dense matrix and produce byte-identical codes /
  fits; the sparse/frame test path is a NEW INTERNAL .Call + facade entry,
  exactly as the training dbartsMixedMatrix container is internal while
  dbarts.h create stays dense. stan4bart (dense test only) is unaffected.
- The engine keeps test raw only where a leaf model reads it. Linear/gp leaves
  gather standardized test covariates (uTest_, [[model.hpp:266@4a521760]], [[model.hpp:669@4a521760]]) via
  rawTestColumn ([[data.hpp:247@4a521760]]); a sparse-BACKED test column REFUSES the leaf
  designation, the plan-4/5 precedent (R_interface_bartcore.cpp leaf-covariate
  refusal on sparse-backed columns).
- Views DENSIFY. buildFromParent gathers a fold's test rows through the
  storage-aware codeAt ([[data.hpp:859-864@4a521760]]) into dense codes; a fold has no
  external x.test (createFromParent, [[R_interface_bartcore.cpp:1532-1690@4a521760]], takes
  testRows indexing the parent). So the fold test surface needs the
  descent/storage refactor but NO new ingestion.
- Codes are identical across storage. A dense factor and a sparseFactor of the
  same values code each row in level order; the CSC column's zeroCode is the
  reference level's ACTUAL level-order code (quantizeCscColumn, [[data.hpp:514@4a521760]],
  the plan-5 refinement), so codeAt(j,i) == the densified code at every row.
  This is why a resident sparse test store fits BITWISE-identically to today's
  densified test store - the densified path IS the oracle.

## Context (seams, read in code)

- Test storage today ([[data.hpp:184-193@4a521760]]): numTestObservations, ownedTestValues
  (owned col-major raw copy, all columns), testCodes (owned ROW-major, all
  columns), testOffset (borrowed). No per-column type, no sparse tier - the
  classic dense pair the training container already replaced.
- buildTest ([[data.hpp:750@4a521760]]) copies the dense raw and quantizes every column into
  row-major testCodes; quantizeTestColumn (541) re-quantizes one column on a
  cut change (setCutPointsForColumn, 481); testRow(i) (940) hands descent a
  materialized row.
- Descent reads a materialized row: findBottomNodeForRow(data, xt) indexes
  xt[rule.variableIndex] ([[tree.hpp:846@4a521760]]); its six call sites pass data_.testRow(i)
  ([[chain.hpp:1107@4a521760]], [[chain.hpp:1120@4a521760]], [[chain.hpp:1137@4a521760]], [[chain.hpp:2144@4a521760]], [[chain.hpp:2176@4a521760]], [[chain.hpp:2216@4a521760]]). A sparse test column has no
  contiguous row, so descent must move to codeAt(var,i) - the change that lets
  test hold typed columns.
- Training container the test side mirrors: ColumnStore per-column dense codes,
  SparseColumnData rank bitmap ([[data.hpp:95-108@4a521760]]), cscSlices for re-quantize
  (144), codeAt storage-aware (951), quantizeCscColumn (510). buildMixed /
  buildFromCsc build it against COMPUTED cuts; test must build against the
  training cut grid (types, numCuts, cutPoints already shared).
- Bridge: parseData x.test is a dense real matrix (R_interface_bartcore.cpp:
  533-547); the training dbartsMixedMatrix branch (361-473) assembles a
  transient dense block into DataHandle::ownedMixedDense (1243-1247), gathers
  the dgCMatrix + per-sparse-column sparseReference/sparseCategoryCount
  (461-471). Categorical test codes are validated against training counts
  (969-988). bartcore_setTestPredictor (2094) / _setTestOffset (2123) are
  whole-object and dense; bartcore_predict (raw-double FlatNode) holds nothing.
- R: @x.test is matrixOrNULL ([[A_class.R:402@4a521760]]), validity ncol==ncol(x)
  (469-487). validateXTest ([[data.R:38-126@4a521760]]) densifies sparse-factor test columns
  (densifySparseFactorColumns, [[utility.R:420@4a521760]] - the interim this plan retires)
  then builds a dense code matrix over the TRAINING levels
  (mapFactorColumnsToTrainingLevels [[utility.R:437@4a521760]], makeCategoricalModelMatrix
  274). The training container is assembled by assembleMixedMatrix
  ([[mixedMatrix.R:115@4a521760]]: dense list + dgCMatrix + map + sparse metadata);
  as.matrix (296) reconstructs, filling implicit rows with the reference code.

## Constraints

- Draw-neutral on every existing path (see rng). The dense test path is the
  guard: any dense-x.test fit that moves means the layout refactor leaked -
  stop.
- Dense x.test remains ACCEPTED and byte-identical (matrix input, dbarts.h
  callers, xbart folds). Sparse/frame test is additive and internal.
- Single-writer, creation-fixed for sparse sources: a sparse test container is
  set whole-object (setTestPredictors replaces the test store); no per-cell
  sparse test mutation (the training-side rule).
- Fast over safe in C/C++: descent reads codeAt only for the columns on its
  path (cheaper than materializing a full row); the sparse test membership
  read reuses SparseColumnData::at.
- OUT of scope: resident sparse storage for the stateless predict() FlatNode
  path (frozen dense, holds nothing - densify at the R boundary, Open Q3);
  weights.test (a vector, no per-column kind); rbart_vi's R-loop predict (over
  data@x). dbarts.h unchanged.

## Scope decisions (recorded)

1. ENGINE REPRESENTATION: the test data is a SECOND ColumnStore-shaped store
   that SHARES the training cut grid (types, numCuts, cutPoints copied at
   build) and OWNS its raw. Rationale: it reuses buildMixed / buildFromCsc /
   codeAt / quantizeCscColumn wholesale - literally one data model - and owning
   raw (dense full + sparse nnz values+rows) keeps buildTest's "copies,
   nothing pinned" contract, so no ownedTestMixedDense borrow and no new PROT
   slot are needed even for a mutated sparse test set. Test owns where training
   borrows because test sets are smaller and re-quantize from the owned copy.
   (Open Q1 offers the parallel-fields alternative.)
2. buildTest DENSE stays (dbarts.h). A new internal build variant takes the
   test container's dense block + CSC slices + per-sparse-column
   reference/K, builds the typed test store against the training cuts. Sparse
   test columns take rank-bitmap or densify at sparseDensityThreshold, the
   training tier rule.
3. predict() DENSIFIES at the R boundary (Open Q3): frame/sparse test input is
   coded against training levels then materialized to a numeric matrix for the
   frozen raw-double FlatNode call; no resident sparse predict store.

## Steps

1. Typed per-column test store + storage-aware descent (engine; dense input,
   BYTE-IDENTICAL refactor). Replace ownedTestValues + row-major testCodes with
   the per-column test store (scope decision 1): dense codes per column, a
   SparseColumnData slot per sparse column, storage-aware testCodeAt(var,i).
   Change findBottomNodeForRow to descend via testCodeAt ([[tree.hpp:846@4a521760]]) and its
   six chain.hpp sites; buildFromParent writes the typed store's dense codes
   (folds unchanged, still densified). buildTest still takes a dense matrix and
   fills dense per-column codes - all sparseSlot=-1, so codes and fits are
   byte-identical. rawTestColumn serves dense-backed columns as today.
   Files: data.hpp, tree.hpp, chain.hpp, sampler.hpp, model.hpp (rawTestColumn
   unchanged interface). Tests: tests/cpp - testCodeAt equals the old
   row-major code at every (var,i) for a dense build; view test codes
   unchanged; buildTest/rawTestColumn contracts updated ([[test_data.cpp:258-267@4a521760]],
   testColumnStoreView:88-93). Gate class: RNG-neutral - equivalence 22/22 +
   tinytest 2803 no regen + tests/cpp. Size: L. Abort: any dense-x.test fit
   moves.
2. Sparse test build path (engine; INERT until step 3). Add the internal build
   variant (scope decision 2): dense block + CSC test slices + per-column
   reference/K -> typed test store against the training cuts; quantizeTestColumn
   becomes storage-aware (dense from owned raw, sparse from the owned nnz copy,
   mirroring quantizeCscColumn:510). Sparse-backed test columns refuse the leaf
   designation. facade.hpp gains the internal setTestData(dense,csc,meta)
   overload (dbarts.h dense entry untouched). Files: data.hpp, facade.hpp,
   sampler.hpp. Tests: NEW tests/cpp - a rank-tier and a densified-tier sparse
   test column bin BITWISE vs a dense test matrix of the same values (codes,
   descent leaf index, one recorded test fit), the plan-5 device on the test
   side. Gate class: RNG-neutral - equivalence 22/22 + tinytest 2803 (unchanged,
   no R producer yet) + NEW tests/cpp. Size: M. (Foldable into step 1 if the
   combined diff stays reviewable.)
3. Ingest a test container (R + bridge). @x.test accepts a dbartsMixedMatrix
   (A_class.R slot union + validity by ncol); validateXTest builds a test
   container coding factors/sparseFactors against the TRAINING level table
   (reuse assembleMixedMatrix / the makeCategoricalModelMatrix path), RETIRING
   densifySparseFactorColumns; dense-matrix and dense-frame x.test stay
   byte-identical. parseData gains an x.test container branch mirroring the
   training branch (transient dense block owned by the store per decision 1,
   CSC slices, reference/K through the internal setTestData); category
   validation reads the container. Files: R/A_class.R, R/data.R, R/utility.R,
   R/mixedMatrix.R, R_interface_bartcore.cpp. Tests: the existing plan-5
   test-side gate ([[test-sparse-factor.R:282-322@4a521760]], sparse-cat x.test ==
   densified) now runs RESIDENT-sparse and must stay identical(); NEW tinytest
   for a mixed/ordinal-sparse x.test ingest + fit; refusal messages narrowed.
   Gate class: RNG-neutral - equivalence 22/22 + tinytest (grows) + the
   test-fit oracle. Size: L. Abort: a resident-sparse test fit differs from its
   densified twin.
4. setTestPredictor / setTestOffset container mutation (R + bridge).
   bartcore_setTestPredictor (2094) gains a container branch: validate +
   rebuild the typed test store via setTestData; the R method
   ([[bartcore.R:326@4a521760]], [[bartcore.R:602@4a521760]]) routes a frame/sparse test set through the step-3
   builder and installs @x.test as the container. The EXISTING per-column
   update surface (setTestPredictor(x.test, column), index or name, dense
   values, [[bartcore.R:326-359@4a521760]]) stays unchanged on a dense-matrix x.test:
   copy-modify R-side, whole-object engine call, today's semantics exactly
   (Q2 resolution). A container-backed x.test takes whole-object replacement
   only - the training-side creation-fixed rule for sparse sources; no
   per-cell sparse test mutation. Offset length checks unchanged. Files:
   R/bartcore.R, R_interface_bartcore.cpp. Tests: NEW/extended tinytest -
   setTestPredictor with a sparse/frame container fits identically to the
   dense equivalent; dense matrix still accepted, including a per-column
   update; removal (NULL) clears. Gate class: RNG-neutral - equivalence
   22/22 + tinytest. Size: M.
5. predict() frame/sparse parity (R; densify-at-boundary, decision 3). Route
   predict's x.test through the container-aware validateXTest (coding against
   training levels), then materialize the numeric matrix the frozen raw-double
   engine predict consumes; no resident predict store. Files: R/bartcore.R,
   R/data.R (shared validateXTest). Tests: tinytest - predict on a sparse/frame
   test set == predict on the dense equivalent, BITWISE. Gate class:
   RNG-neutral - tinytest. Size: S.
6. Docs + landing. sparseFactor.Rd / dbartsData.Rd (sparse/frame x.test now
   resident, densification retired); the design landing note (this plan
   supersedes plan-5 decision 2 / Q4 and data-ownership.md "Implementation record" item 5's
   densification reading); mark the "sparse x.test" bullet REMOVED from
   docs/plans/sparse-extensions.md scope (the other three extensions stay - do
   not edit that file's steps beyond the scope note). Files: man/*.Rd,
   docs/design/data-ownership.md, docs/plans/sparse-extensions.md (scope line
   only). Gate class: R CMD check man + tinytest. Size: S.

## Verification

- Full tinytest 2803 per commit from a preclean install (--preclean after the
  header/facade edits; delete stale tests/cpp binaries).
- Equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at EVERY commit -
  guards that training draws never move (the equivalence datasets exercise no
  test path, so identity here is necessary, not sufficient).
- THE TEST-FIT ORACLE (sufficiency): for any dataset expressible both ways, the
  recorded test fits of a resident-sparse / frame x.test are BITWISE-identical
  to the dense-matrix / densified x.test - [[test-sparse-factor.R:282-322@4a521760]] extended
  to resident storage, plus the tests/cpp store-level check (step 2). If a new
  cpp test draws from the shared rngState, snapshot/restore around its draws
  (the plan-4 finding).
- Memory: a container-bytes comparison, dense vs resident-sparse x.test, on a
  high-cardinality test frame (bairrtt shape) - the win the densification
  interim gave up. Deterministic, any box.
- rchk on the next scheduled run (steps 3, 4 touch the bridge).
- dbarts.h unchanged; no stan4bart lockstep (the sparse test path is internal
  .Call + facade). One confirmatory bench-sampler compare vs the current
  baseline (dense test path must not regress; a resident-sparse arm at parity).

## Open questions for VD

RESOLVED (VD, 2026-07-14): Q1, Q3, Q4 as recommended. Q1 addendum, decided
at step-1 implementation per the recorded caveat: parallel per-column test
fields inside ColumnStore (testCodes/testCodeOffsets/testSparseSlot/
testSparseColumns), not a second store - sharing the cut grid by identity
keeps the setCutPointsForColumn -> quantizeTestColumn contract with no grid
copy to sync, and a second store would have carried the cut-building and
mutation surface the test side never uses (the caveat's trigger). Byte
identity was unaffected either way. Q2 amended: support
whatever was supported previously - the per-column dense update surface
(setTestPredictor(x.test, column), index or name) stays exactly as today;
container-backed x.test is whole-object only (new surface, nothing prior
to preserve). Step 4 records the amendment.

- Q1 (engine representation: second ColumnStore vs parallel test fields).
  RECOMMEND the second store sharing the training cut grid: it reuses
  buildMixed/buildFromCsc/codeAt/quantizeCscColumn verbatim, which is the
  directive's "one data model" in the most literal form, and confines the test
  half to one build variant. What would change it: if sharing the store struct
  drags in cut-building / mutation surface the test side never uses, parallel
  per-column test fields (testTypes, testSparseColumns, testCodeAt) inside the
  existing ColumnStore are lighter, at the cost of duplicating the quantize
  logic. Either is byte-identical; this is a code-shape call, decidable at
  implementation from the actual diff.
- Q2 (setTestPredictor granularity: whole-object vs per-column-kind mutation).
  RECOMMEND whole-object: sparse sources are creation-fixed everywhere else in
  the container (single-writer rule), a test set is set once, and whole-object
  reuses the step-3 builder with no new invariant. What would change it: a named
  caller that swaps one test column per outer iteration (an IRT test-scoring
  loop) would justify a per-column setTestColumn, but none exists; add it when
  one does, as a sparse-extensions-style consumer-gated follow-up.
- Q3 (predict() resident sparse vs densify-at-boundary). RECOMMEND densify at
  the R boundary: dbarts_sampler_predict is frozen dense, predict is stateless
  and holds nothing, and the FlatNode path reads resolved double cut values, not
  codes - a container-aware predict would duplicate the descent for no resident
  memory saved. What would change it: a huge sparse test set predicted where the
  transient dense matrix is the memory ceiling would motivate an internal
  container-aware predict entry (not touching dbarts.h).
- Q4 (test levels absent from training). RECOMMEND keep today's dense semantics:
  test factor / sparseFactor levels are coded against the TRAINING level table,
  a level absent from training maps to the missing-category code exactly as the
  dense path does (mapFactorColumnsToTrainingLevels), and the engine's
  NA-direction handling routes it. This preserves the sparse-vs-dense oracle; a
  sparseFactor forbids stored NA, but an absent-level -> NA happens at the R
  coding step, before the engine, identical to dense.

## Relationship to sparse-extensions.md and drift to record

- sparse-extensions.md lists "sparse x.test" among four consumer-gated deferred
  extensions (context bullet, referencing sparse-columns.md). This plan
  DELIVERS it; step 6 marks that bullet removed from that doc's scope, leaving
  the other three (in-place nonzero mutation, streaming range kernel,
  dense-column mixed mutation) closed and consumer-gated. Do not edit
  sparse-extensions.md beyond the scope note.
- data-ownership.md "Implementation record" (item 5) reads "test-side sparse ... satisfied by
  densification, not resident sparse testCodes (the test-data-parity backlog
  item supersedes that interim)"; the landing note records that this plan IS
  that item and that resident sparse testCodes now ship, superseding plan-5
  decision 2 / Q4.

## Landing notes

Commit 1 = 14bef56. Store test codes per column and descend storage-aware
(engine, dense-input, byte-identical refactor). ownedTestValues + row-major
testCodes retired for parallel per-column test fields inside ColumnStore
(testCodes/testCodeOffsets/testSparseSlot/testSparseColumns) sharing the
training cut grid by identity - the Q1 addendum, decided from the actual
diff shape rather than a second store (recorded inline above, Open
questions for VD). findBottomNodeForRow and its six chain.hpp call sites
descend via testCodeAt instead of a materialized row; buildFromParent's
fold test surface is unchanged, still densified (no external x.test on a
fold). Gate: equivalence 22/22 identical + full tinytest 2803 (no regen)
+ tests/cpp (testCodeAt equals the old row-major code at every (var,i)).

Commit 2 = 95bdf00. Build sparse test columns against the training cuts
(engine, inert until step 3). buildTestMixed added: a dense block + CSC
test slices + per-column reference/K build the typed test store against
the training cuts (shared by identity), rank-bitmap or densified per
column at sparseDensityThreshold (0.2), mirroring quantizeCscColumn.
sampler.hpp gains the internal setTestData overload, refusing (store
untouched) a leaf-covariate column that would be CSC-backed. DEVIATION
recorded: setState's cut-restore snapshot (sampler.hpp) saved
oldSparseColumns for rollback but had no analogous save for the new
testSparseColumns field; oldTestSparseColumns was added alongside it so a
rejected cut change rolls the test-side rank columns back too - a gap the
new sparse test build surfaced (inert on every existing dense-test path,
since testSparseColumns stays empty there). Gate: equivalence 22/22 + full
tinytest 2803 (unchanged, no R producer yet) + NEW tests/cpp (a rank-tier
and a densified-tier sparse test column bin bitwise vs. a dense test
matrix of the same values).

Commit 3 = a4ec3cd. Ingest mixed test containers and retire test-side
densification (R + bridge). dbartsData's x.test slot widened from
matrixOrNULL to ANY with a validity guard (is.matrix(x.test) ||
inherits(x.test, "dbartsMixedMatrix")) INSTEAD OF A CLASS UNION - Matrix
classes cannot appear in a slot union without Matrix loaded at package
build time, the same constraint @x already works around. validateXTest's
sparse branch recodes sparseFactor/sparse-ordinal test columns against
training levels and keeps them in a container
(remapSparseFactorToTrainingLevels replaces densifySparseFactorColumns,
which is deleted); parseData gains a dbartsMixedMatrix x.test branch
mirroring the training branch, reusing the sparseReference/
sparseCategoryCount plumbing, and installs through the new setTestData.
BEHAVIOR NARROWING recorded: a sparse-bearing test frame against a
formula-trained design with no unexpanded factor levels (factorLevels
NULL - a fully numeric or already-dummy-expanded design) now REFUSES
cleanly ("sparse test predictor columns require a categorical training
design; supply 'x' through the x/y interface") instead of silently
densifying then replaying through model.frame, which the old code's
termLabels branch allowed unconditionally. LATENT NULL-DEREF FIX recorded:
validateCategoricalPredictors' test-code bound-check loop called
rawTrainingColumn(data, j) unconditionally for every categorical column
and indexed data.x_test directly; a CSC-backed training-mixed categorical
column paired with any dense x.test (already reachable pre-plan, since
dbarts.h's test surface was dense-only but could pair with an R-side
mixed x/y training container) would have dereferenced a null pointer -
never exercised by an existing test, so latent rather than observed. Fixed
by skipping CSC-backed training columns explicitly and by routing test
reads through the new rawParsedTestColumn helper, which returns null
(skip) for a CSC-backed test column instead of indexing data.x_test. Gate:
equivalence 22/22 + tinytest (grows) + the test-fit oracle
([[test-sparse-factor.R:282-322@4a521760]] extended to resident storage, plus new
mixed-container and leaf-covariate-refusal cases).

Commit 4 = 7d781a8. Take mixed containers through the test mutation
surface (R + bridge). setTestPredictor/setTestPredictorAndOffset gain a
dbartsMixedMatrix branch. DEVIATION recorded: a4ec3cd had inlined the
container-parse logic once, in parseData; commit 4 factors it into a
SHARED parseTestContainer helper (plus an installTestContainer wrapper
around setTestData) used by all three call sites - parseData,
setTestPredictor, and setTestPredictorAndOffset - rather than duplicating
the parse a second and third time. Each mutation validates and rebuilds
the typed test store before any store mutation, so a refused container
(leaf-covariate-on-sparse) leaves the engine-side store untouched. R-side
(bartcoreSamplerSetTestPredictor): a container-backed x.test refuses a
per-column update ("cannot update a single column of a sparse test
matrix; replace the whole test matrix instead"); a whole-object
replacement installs the new @x.test optimistically, then a tryCatch
ROLLS BACK @x.test/@offset.test to their prior values if the bridge call
errors, keeping the R5 object and the engine's prior store consistent on
refusal - the same discipline setPredictor already used. bartcore_setData's
existing container refusal ("requires a dense predictor matrix; sparse
predictors fix the design at creation") gained a parallel test-side
refusal KEPT from a4ec3cd ("requires a dense test matrix; a sparse test
set fixes the design at creation") - the symmetric fixed-at-creation
contract already applied training-side. Gate: equivalence 22/22 +
tinytest (grows): setTestPredictor/AndOffset with a sparse/frame container
fits identically to the dense equivalent; per-column update (by index and
by name) and NULL removal on a dense x.test unchanged; leaf-covariate
refusal on mutation is inert (the prior store, and the sampler, remain
usable afterward).

Commit 5 = 22d7116. Take frames and sparse test sets through predict and
tree replay (R; densify-at-boundary, decision 3). predict and
getTrees(newdata = ) route their argument through the same
container-aware validateXTest used at creation (coding against training
levels), then materialize it to a numeric matrix before the frozen
raw-double engine calls - no resident sparse predict store, per Q3's
resolution. This is the same boundary fix in both entry points (getTrees'
newdata shares predict's validated/materialized path). Gate: tinytest
(grows): predict/getTrees on a sparse/frame test set bitwise-identical to
the dense equivalent, including a column-reorder case and an
absent-training-level refusal identical to the dense path's.

Commit 6 = this commit (docs). man/sparseFactor.Rd, man/dbartsData.Rd, and
man/dbartsSampler-class.Rd updated to record resident test-side storage
through creation and setTestPredictor, and the predict/getTrees(newdata =)
densify-at-boundary; man/dbarts.Rd's formula-argument paragraph corrected
too ("test matrices remain dense" was stale, found beyond the three named
files). docs/design/data-ownership.md item 5 records that this plan
supersedes the densification interim. docs/plans/sparse-extensions.md's
sparse-x.test context bullet marked delivered (scope line only, the other
three extensions untouched). Gate: R CMD INSTALL (Rd files parse clean);
full tinytest 2832 PASS, 0 FAIL (docs-only, no growth from commit 5);
tools::checkRd clean (no NOTEs/warnings) on every man/ file.

Gates held at every step: equivalence 22/22 identical throughout (training
draws never moved - the guard this plan's rng section required at every
commit, run by each step's implementer and re-run independently at
review for the engine/bridge commits); full tinytest grew 2803 ->
2832 across commits 3-5 with no snapshot regeneration; the tests/cpp
component suite gained the step-1 testCodeAt-equals-old-row-major-code
check and the step-2 rank-tier/densified-tier sparse test-store bins,
both passing at HEAD.

Memory measurement (verification item, run for this commit). A throwaway
script (not in the repo) built n.train = 2000, n.test = 1e5, one dense
numeric column, and a K = 200-level sparseFactor test column at ~5%
non-reference (nonzero) density (the bairrtt-shape high-cardinality test
frame), once as a dense x.test matrix and once as the sparse-bearing
data-frame form. Confirmed the dense build stays a plain matrix and the
sparse build installs a resident, sparse-backed dbartsMixedMatrix (not
densified), with bitwise-identical test fits between the two arms.
Measured object.size(data@x.test): 1,600,568 bytes (dense matrix) vs.
876,048 bytes (resident container) - a 1.83x ratio, modest because both
representations carry the shared dense numeric column (800,000 bytes)
undiminished; the categorical column is what the sparse path actually
shrinks. Isolating it analytically (dense test codes are u16, 2 bytes/row;
a rank-bitmap test column costs 0.1875 bytes/row fixed plus 2 bytes per
stored code, measured f = 0.0496, well under the 0.2 sparseDensityThreshold
that would fall back to a densified tier) gives 200,000 bytes (dense) vs.
28,672 bytes (rank) for the categorical column alone - a 6.98x code-shrink,
the win the densification interim gave up. Engine-side resident bytes were
not cheaply queryable from R (no exposed accessor), so this is the R-side
container ratio plus the per-column code-shrink argument, as anticipated
by the plan's memory verification text.

Confirmatory bench (quiet window, arm64). bench-sampler compare vs
benchmarks/baselines/bench-sampler-4008675.csv: "OK: no metric regressed
more than 5%", zero flags - the dense test path did not move. The
resident-sparse test arm (n.train = 1e4, K = 200 sparseFactor at f = 0.05
plus a dense column, 75 trees, test fits recorded every iteration,
interleaved repeats): 1.04x dense at n.test = 1e4 (5.26 vs 5.07 ms/iter)
and 1.11x at n.test = 1e5 (40.2 vs 36.3 ms/iter), the rank-bitmap
membership lookup on test-row descent vs a direct dense code read. Not
the strict parity the verification text hoped for, but the cost scales
with the test-row share of the workload, is opt-in (a dense x.test keeps
dense storage and today's speed exactly), and buys the 6.98x per-column
code shrink above; a speed-sensitive caller can densify sparse test input
themselves. The plan is CLOSED.
