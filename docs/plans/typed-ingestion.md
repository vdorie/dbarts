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

## Slice 1 steps (dense declared-K plumb-through)

Baseline check (2026-08-05): every dense factor in the equivalence
scenarios is factor(sample(...)) - levels computed from the observed
data, declared K == max+1 by construction - and the BCF/multinomial
harnesses carry no factors, so the trio must stay bitwise. Gap
datasets (declared > observed) shift, and are the point.

Critique record (2026-08-05, refuting): NO BLOCKER, six amendments
folded in below. Its counterfactual (buildCutsForColumn probe-patched
to take a declared K): neutrality bitwise by SHA-1 at equal K; gap
declared-4/observed-3 accepted end-to-end at every entrance; declared
70/observed 60 genuinely enters the pooled tier with store and trees
agreeing. Empty categories are ALREADY engine-native - the CSC path
fits them today, empty-leaf-veto.md documents the state, both scoring
paths veto, and every xbart fold view holds them (buildFromParent
copies numCuts verbatim). K reaches every consumer solely through
numCuts; categorical grids never rebuild (setData/refreshCuts/
ScopedCutGrid/setState all skip, setCutPoints refuses).

1. Engine: generalize the declared-count channel - SamplerOptions
   carries per-column category counts for ANY categorical column (the
   CSC-only cscCategoryCounts becomes the general spelling; the
   storage-keyed SCAN selection in the slice-0 validators stays -
   only the bound generalizes; the channel is creation-time only, the
   ctor nulls it). buildCutsForColumn takes max(declared, inferred) -
   never a declared K below observed, which would reopen the slice-0
   class via hand-flipped varTypes.
2. Bridge: parseData reads the declared counts off attr factor.levels
   lengths, gated on varTypes[j] == CATEGORICAL (ordered factors
   carry a level table but are ordinal); absent = infer (0 is the
   existing absent spelling). optionsFromParsed publishes;
   trainingCategoryBound prefers the declared count for dense
   columns, and the dense TRAINING bound tightens to it (today it
   bounds only by maxCategories).
3. R: NO new dbartsData slot and no initialize change - attr(x,
   "factor.levels") already carries the complete per-column level
   table in creation order on BOTH container flavors
   (assembleDenseColumnMatrix AND assembleMixedMatrix), preserved by
   subsetting and already consumed by decodeCategoricalSplits. The
   bridge reads it directly.
4. Docs in the SAME commit: dbartsData.Rd/dbarts.Rd note declared
   factor levels are honored end-to-end; NEWS bullet covers the
   acceptance change AND the restore break: a state saved pre-slice-1
   for a factor whose declared K crosses the 63/64 tier (>= 64
   declared, fewer observed) refuses to restore ("missing required
   block 'tree.masks'"); same-tier states restore silently and
   correctly.
5. tinytest: gap factor accepted with finite fits at creation,
   setTestData, predict, and setPredictor - INCLUDING a dense factor
   inside a MIXED container (the assembleMixedMatrix flavor, which a
   dense-only test would miss); the G2/G4 asymmetry test flips to
   symmetric; re-point the planted-code gate (test-sparse-factor.R
   dense-count matching) at a genuine gap; the tier-crossing restore
   refusal and a same-tier restore. NO in-file bitwise snapshot (the
   trio is the neutrality gate; in-file snapshots are the fragile
   spelling per CLAUDE.md).
6. RNG-shift inventory before landing: run benchmarks/R/
   categorical-exact.R (exact-posterior categorical gate - expected
   to shift only for gap designs, verify) and enumerate any RNG-locked
   tinytest whose data carries a gap factor (expected none; verify).

## Slice 1 verification

R CMD INSTALL --preclean (engine headers move); tests/cpp clean-build
(cut-grid cases); full tinytest; equivalence trio bitwise;
categorical-exact.R; air format --check .; full local R CMD check (R
and Rd move); CI incl. sanitizers green.

## Slice 1 landing

LANDED 0d3914c 2026-08-06, all CI green incl. sanitizers. 421+/78-:
SamplerOptions.cscCategoryCounts generalized to categoryCounts
(0 = infer, creation-time only); ColumnSource.declaredCategoryCount;
buildCutsForColumn takes max(declared, inferred) on BOTH storage
kinds (a slice-0-validated container has inferred <= declared, so the
CSC extension is bitwise-neutral - trio confirmed); parseData reads
attr(x, "factor.levels") lengths gated on CATEGORICAL varTypes with
CSC-backed columns skipped (their container stays authoritative),
factorLevelsExpr PROTECTed outright (rchk-friendly, gctorture-checked);
trainingCategoryBound and the dense training bound tighten to the
declared count. No R/ change (the attribute reaches .Call intact on
all three container flavors). Docs: Rd notes + NEWS bullet covering
the acceptance change and the 63/64 tier-crossing restore refusal.
Tests: test-data-categorical-declared.R (22 assertions, suite 3531);
the planted-code gate re-pointed at a genuine gap; both files verified
to FAIL against the pre-slice-1 build. RNG-shift inventory: none -
categorical-exact.R byte-identical (hand-typed matrix carries no
factor.levels), no RNG-locked test carries a gap factor (verified by
seed replay). Deviations accepted at review: declared counts wider
than maxCategories are ignored (inferred governs; R caps at 65535
anyway); the tier-crossing restore test constructs the pre-slice-1
state shape rather than replaying a literal old state.

Slice 2 remains, SUB-SLICED 2026-08-06 (orchestrator, under the
grant): 2a the mutation-shape lift + the alignment residual - R5-side
maintenance of mixed containers under mutation (installPredictorColumns
maintains only pure dgCMatrix today; R/bartcore.R refuses mixed
wholesale), the engine pattern rule keyed on refCode for categorical
CSC columns (mutateCscColumnFromDense's {i : value != 0} is the
ordinal rule; categorical is {i : code != refCode}), the slice-0
refusals lifted once both hold, and foreign-container level-order
alignment via the training factor.levels tables (the R-side
mapFactorColumnsToTrainingLevels treatment, today applied only to
data.frame test inputs at R/data.R). 2b the PredictorSource view
consolidation proper (the 9 SamplerOptions ingestion fields, the
7-arg setTestData, the mutation entries -> one spelling). Each
sub-slice: step list here + refuting critique before implementation,
2a first (user-facing; 2b is hygiene).

## Slice 2a steps (mutation-shape lift + alignment)

Design decisions carried in: the R5 value spelling for mutating a
sparse CATEGORICAL column is the dense one - a whole-column vector of
codes in the column's FIXED level table (declared K and reference are
creation-pinned, the setResponse scale-pinning analog); factor-typed
inputs re-coded by the wrapper wait for the view (2b). Neutrality:
every currently-accepted path stays bitwise (trio); currently-refused
paths become capability.

Critique record (2026-08-06, refuting): BLOCKER x2 + 7 weaknesses,
folded in below. Verified sound: the refCode-keyed engine rule on
both storage tiers (16/16 probe assertions incl. NaN; identity
bitwise inert; dense-equivalent exact since codeFor(refCode) ==
refCode); rollback (snapshotCscColumn restores byte-for-byte);
the transaction is all-or-nothing (runPredictorTransaction prechecks
then applies-or-restores all - no partial state); per-observation and
setData refusals survive the lift independently (bridge
columnIsCscBacked guards); a LEVEL-PERMUTED container needs NO
pattern rebuild - re-coding is a bijection on levels, so
non-reference cells stay non-reference (the orchestrator's rebuild
premise refuted; C1/C4 do not regress). The blockers: Matrix [<-
runs drop0 MATRIX-WIDE (writing one column corrupts every explicit
zero in the container mirror, untouched categorical columns
included), and a remapped container also needs its declared K lifted
or resolveCscCategoricalReferences falsely refuses.

1. Engine: mutateCscColumnFromDense keys the nonzero pattern on the
   column's kind - ordinal {i : value != 0} unchanged, categorical
   {i : code != refCode} with NA stored (NaN != refCode). tests/cpp
   asserts on CODES, not the stored pattern (a hand-built
   non-canonical sparseFactor makes an identity re-install a
   pattern change while staying code-inert): identity bitwise
   no-op; genuine change matches the dense-equivalent build.
2. R5: installPredictorColumns dispatches PER COLUMN on
   sign(x$map[columns[k]]) (a mixed updatePredictor can name both
   kinds in one call); the sparse-backed branch is DIRECT SLOT
   SURGERY on x$sparse@i/@p/@x - splice the column's entries out and
   the new {i : code != sparseReference[k]} entries in, shifting @p -
   NEVER Matrix [<- (its drop0 is matrix-wide and strips the
   explicit zeros a non-first-level reference requires; disqualified
   for the ordinal branch too). sparseReference/sparseCategoryCount
   stay creation-pinned. Compute the new slots BEFORE the .Call so
   the post-acceptance install cannot throw and leave data@x stale.
   The wholesale mixed refusal (R/bartcore.R ~57) lifts for
   column-granular setPredictor/updatePredictor; the WHOLE-MATRIX
   refusal (~127) STAYS - its install path would replace the
   container with a bare matrix (recorded door, 2b); its test in
   test-data-mixed.R changes text only. The per-observation refusal
   NARROWS to sparse-backed columns (x$map[column] < 0): the
   engine deliberately keeps dense-backed per-observation mutation
   open (the IRT latent case, sparse-columns.md Extension (i)).
3. Bridge/C API: refuseCscCategoricalMutation (slice 0) deletes from
   bartcore_setPredictor/updatePredictor and both dbarts_sampler_*
   entries once 1-2 hold; its refusal tests flip to acceptance cases
   with fits equal to the dense-equivalent design.
4. Alignment (the slice-0 residual): one container branch at the
   validateXTest funnel (R/data.R) covers every entrance (creation,
   setTestPredictor, predict, getTrees) by construction. Re-code a
   foreign container's factor columns - sparse AND dense-backed, the
   latter stored as real factors in $dense - against the training
   factor.levels tables via the remapSparseFactorToTrainingLevels
   treatment (rebuild with levels = trainingLevels), which lifts
   codes, reference, AND sparseCategoryCount together; genuinely-new
   levels refuse with the existing message; the level-permuted case
   then fits IDENTICALLY to its aligned equivalent (critique p4
   verified end-to-end). Sparse @x remaps by slot surgery, not [<-.
5. tinytest: column-granular mutation over mixed designs x
   {dense-backed, sparse ordinal, sparse categorical} x
   {setPredictor(column), multi-column updatePredictor incl. mixed
   kinds in one call}; a sparseFactor with reference != levels[1]
   (the only shape with explicit zeros - no existing test has one);
   mutate a sparse ordinal column beside a sparse categorical one
   and assert the categorical data@x slice bit-unchanged; a
   save/re-create round trip (fresh sampler from mutated data@x
   matches the mutated sampler - the only real mirror gate);
   rollback leaves data@x untouched; alignment acceptance + refusal;
   the flat C API mirror.
6. Docs in the SAME commit: NEWS bullet (mixed-container column
   mutation now supported; the stale "deferred at the R level"
   bullet corrected); man/dbartsSampler-class.Rd mutation prose;
   fix man/dbarts.Rd's stale claim that the dgCMatrix mutation
   surface is fixed at creation (extension (i) falsified it).
7. Out of scope, unchanged: whole-matrix setPredictor on mixed
   designs (door above), sparse x.test as bare dgCMatrix, the
   streaming range kernel, per-column u8 widths, the 2b view.

## Slice 2a verification

Slice-1 battery minus categorical-exact (no grid change): R CMD
INSTALL --preclean; tests/cpp clean-build + ASAN; full tinytest; trio
bitwise; air format --check .; full local R CMD check (R, Rd, and
tests move); CI incl. sanitizers green.

## Slice 2a landing

LANDED 95db5d7 2026-08-06, all CI green incl. sanitizers. 793+/162-:
mutateCscColumnFromDense keys the pattern on column kind (categorical
implicit = refCode); the slice-0 refuseCscCategoricalMutation guards
deleted (bridge + C API); installPredictorColumns dispatches per
column with the sparse branch as direct @i/@p/@x slot surgery
(replaceSparseColumn; Matrix [<- disqualified, its drop0 is
matrix-wide), computed before the .Call; the wholesale mixed refusal
lifted for column-granular paths, whole-matrix stays refused (2b
door), per-observation narrowed to sparse-backed columns (the
dense-backed IRT case opens); foreign containers align at the
validateXTest funnel (alignContainerFactorLevels lifts codes,
reference, and declared K together; new levels refuse). Suite 3531 ->
3581 (0 deleted; 3 refusals flipped to acceptance with
dense-equivalence bitwise gates, 11 container code-bound refusals now
refuse at the alignment funnel with the levels message - an
out-of-range code in a container IS an unseen level - with two new
container cases keeping the slice-0 code-bound lock-in alive, 1
whole-matrix text change). tests/cpp 181 incl. sparse categorical
mutation on both tiers asserted on CODES. Trio bitwise; R CMD check
OK; probe battery: the B1 cross-column corruption case clean through
the R5 path, p4b raw foreign container accepts and fits identical to
its data.frame equivalent, p5b dense-backed container factor
re-coded to truth. Deviation accepted at review: the critique's
p4(iii)/(iv) hand-remapped objects are self-inconsistent (codes
remapped, level table not) and correctly refuse; the specified case
is the raw container, which accepts. Save/load reproduces to
tolerance, pre-existing (verified against the pre-change build).

## Slice 2b census and defect record (2026-08-06)

Census (design agent, orchestrator-verified anchors): SamplerOptions
nulls 9 borrowed ingestion fields post-construction - 8 typed
(columnTypes, cscColumnPointers/RowIndices/Values, mixedDenseValues,
columnSources, categoryCounts, cscReferenceCodes) plus
maxNumCutsPerVariable (grid, stays) - and the x ctor argument they
modulate; the nulling repeats across 4 Sampler ctors and 3 Chain
copies (chain.hpp:459/644/682). 7 ColumnStore builders; 8 facade
virtuals in scope; bridge parseData is 3 x flavors x 2 x.test
flavors plus createDataHandle's own 3-way dispatch; 7 of 35 flat C
API entries touch ingestion/test/mutation, all SEXP or dense double*.

REACHABLE DEFECT, newly R-facing via 2a's lift, probe-verified and
independently re-run by the orchestrator: a mixed container's
dense-backed columns are denseBorrowed, and ColumnSource::denseRaw is
written ONLY by buildMixed (data.hpp:884) / buildTestMixed (:1030) -
no mutation path repoints or writes through it (setPredictors,
setColumns, setColumnJournaled, setCell write codes only). Every raw
reader sees creation-time values after a mutation: (1)
setCutPointsForColumn -> rawColumnForRequantize silently REVERTS the
column's codes (R passes currentPredictors = NULL for a
sparse-bearing container, rawPredictorMatrix in R/bartcore.R) -
probe: re-installing the store's own grid is an exact no-op on a
plain matrix (ssr 1993.6985504995 both sides) but moves the mixed
container 2121.78 -> 404.34, the pre-mutation fit; (2) linear/GP
regatherTrainingCovariates/reinitialize (model.hpp) regress on the
old covariate for the rest of the run - the IRT case; R/model.R
permits mixed containers and refuses only sparse-backed covariates.
Reachable through setPredictor(column=), updatePredictor, and
updatePredictorPerObservation on dense-backed columns; a THIRD
entrance (critique): setState replays cut points through the same
NULL-currentPredictors funnel (R/dbarts.R:1241/1399 -> bridge
:4146). NOT reachable before 2a. Pure dgCMatrix (ownedCsc*) and
dense columnar containers unaffected.

Decision (orchestrator, under the grant; slice-0 precedent): fix
FIRST in its own commit (FIX-0 below) - the fix IS 2b's ownership
rule for the mapped-dense case, and a correctness fix must not ride
a large hygiene diff. Folding it into 2b step 1 recorded and
declined.

Critique record (2026-08-06, refuting, with a counterfactual FIX-0
build): BLOCKER x4 + 7 weaknesses, folded in below. Verified sound:
the defect record, both legs re-derived independently (max|mutated -
never-mutated| = 0 on the mixed container against a moving dense
control; the linear-leaf leg isolated by zeroing the covariate's
split probability) plus the setState third entrance above; NO
restore leg needed (saved state carries no raw; re-create rebuilds
from data@x); no test-side staleness (no in-place test writer
exists - container column updates refuse, matrix updates reinstall
wholesale); the buildTest/buildTestMixed collapse bitwise (codes,
offsets, owned buffers, NaN, declared-K, both quantile modes);
buildFromCsc deletion bitwise engine-side; the xIsSparse early exit
dead (parseData:816 refuses categorical + dgCMatrix on every
entrance); the door already accepted below the R layer (whole-matrix
on mixed, bare dgCMatrix, and dense-backed-factor containers return
TRUE, re-create OK, out-of-range codes still refused; one-pass
splice is speed only); the options/ParsedData collapse neutral
(chain.hpp copies null only already-null fields; nothing reads the
8 fields post-construction); dbarts.h names no internal type.
Counterfactual: full tinytest 3581/0 under the probe FIX-0;
tests/cpp green except the address-identity assert
(test_model.cpp:3040, an expected test change, folded in). The
blockers: the dense-backed rollback carries no raw
(ColumnCodeRollback/WholeMatrixUpdate - a rolled-back transaction
left the owned block holding the REJECTED values, and the next
setCutPoints installed them); setCell and setColumnJournaled were
unnamed as writers; the mutation-refusal predicate admitted a
mapped dense view (every kernel indexes values + k * n); the step-6
gate named the function step 1 deletes.

## Slice 2b-pre steps (FIX-0: the store owns its mixed dense block)

1. Engine: buildMixed copies its dense block into a new
   ColumnStore::ownedDenseValues and points denseRaw there - the
   buildTestMixed treatment (ownedTestValues), which is already
   correct. Delete BartcoreHolder::ownedMixedDense and
   DataHandle::ownedMixedDense (no readers; unshipped types).
   Memory shape stated honestly: steady-state zero net, but
   creation PEAKS at 2x the dense block (the copy overlaps the
   parse-local assembly that today is moved). ALL FOUR writers
   write through the owned block beside the quantize:
   setPredictors, setColumns, setColumnJournaled (whose stale
   "Dense stores only" comment goes), and setCell (the
   per-observation IRT path, sampler.hpp:1352 - one cell, not a
   column memcpy). Transactional rollback: the dense-backed legs
   (ColumnCodeRollback/WholeMatrixUpdate) carry NO raw, so
   runPredictorTransaction snapshots and restores the owned raw
   PER TOUCHED COLUMN beside oldGatheredRaw (a whole-block
   snapshot would make the IRT per-sweep transactional update
   O(n * p_dense)); restore by memcpy, NEVER move-assign, which
   relocates the buffer and dangles every denseRaw. Add
   ColumnStore(const ColumnStore&) = delete: denseRaw is now
   self-referential and a copy would silently alias (production
   only moves, which keep the heap buffer stable).
2. tests/cpp: mutate a dense-backed column of a mixed store, then
   setCutPointsForColumn with the identical grid - codes unchanged;
   a linear-leaf covariate over a mutated dense-backed mixed column
   matches the dense-equivalent build after regather; a ROLLED-BACK
   transactional update leaves the owned raw byte-unchanged (the
   B1 regression). Update the address-identity assert
   (test_model.cpp:3040-3043, rawColumn(0) == denseSource.data())
   to value equality - the expected test change under ownership.
3. tinytest: the probe regression (mutate with updateCutPoints =
   TRUE, hand the store back its own grid, fits unchanged), a
   linear-leaf mixed-covariate mutation case, a refused/rolled-back
   updatePredictor followed by setCutPoints (fits unchanged), and
   the setState entrance (mutate, save, setState - fits unchanged);
   the first two verified to FAIL against the pre-fix build.
4. Docs in the SAME commit: NEWS bullet (mutating a dense-backed
   column of a mixed container no longer leaves cut-point
   requantization, setState cut-point replay, and linear/GP leaf
   covariates reading creation-time values).

## Slice 2b-pre verification

R CMD INSTALL --preclean; tests/cpp clean build + local ASAN (new
reachable engine code); full tinytest; trio bitwise (the owned copy
is value-identical; mutation is outside the trio); air format
--check . on any R touch; full local R CMD check; CI incl.
sanitizers green BEFORE 2b lands on top.

## Slice 2b-pre landing (FIX-0)

LANDED 8872dd2 2026-08-06, all CI green incl. sanitizers. 502+/48-:
buildMixed copies its dense block into ColumnStore::ownedDenseValues
and points denseRaw (now double*) there;
BartcoreHolder::ownedMixedDense and DataHandle::ownedMixedDense
deleted; all four writers (setPredictors, setColumns,
setColumnJournaled, setCell) write the new raw through the owned
block; runPredictorTransaction snapshots the touched columns' raw
per column (OwnedDenseRollback) beside oldGatheredRaw and restores
by memcpy on the reject leg; ColumnStore move-only (a copy would
alias the self-referential denseRaw). tests/cpp 182
(+testMixedDenseOwnership: grid re-install code-inert, linear-leaf
mixed/dense bitwise after covariate mutation, rollback restores raw
byte-for-byte on both whole-matrix and single-column granularity,
refused values never reach codes; the address-identity assert
flipped to value equality) - the new cases verified to FAIL on the
pre-fix engine and again with restoreOwnedDenseColumns removed.
tinytest 3590 (+9, test-data-mixed-mutation.R: cut-point
regression, rollback, setState replay, linear-leaf regather, with
dense controls; 4 verified to FAIL on the pre-fix build). Trio
bitwise 27/27 + every BCF and multinomial channel; R CMD check
--as-cran OK. Deviations accepted at review: denseRaw const
double* -> double* (the store owns both sides; removes writer
const_casts); two adjacent stale comments corrected (the mutation
section's dense-only claim, buildMixed's all-borrowed claim).

## Slice 2b steps (PredictorSource view consolidation)

Design decisions carried in (forks decided under the grant,
alternatives recorded): ONE borrowed POD, PredictorSource (data.hpp,
beside ColumnSource), trivially destructible with a static_assert -
the SamplerShape precedent (a bridge entry point Rf_errors past
destructors, so no owning field). Fields: numRows, numColumns;
denseValues + the CSC triple + columnSources; the typing channel
columnTypes, categoryCounts, referenceCodes. Per-column rule, the
existing convention plus one case: columnSources == nullptr means
column j is dense column j of the call's transient block; a
nonnegative entry names a dense column of a host-pinned denseValues;
a negative one the complement of a CSC column. That one convention
absorbs all train and test builders and retires buildFromCsc (the
bridge fills a ~j map for a bare dgCMatrix). Ownership stated once
on the struct: borrowed for the consuming call; a mapped source's
CSC slices are retained by the train store (2a's ownedCsc*), its
dense block is copied (FIX-0), so no host pins anything for a mixed
build; a null map retains nothing (the store re-quantizes from the
caller's later matrices, byte-for-byte today's dense build); the
test store copies everything as now. The grid spec (maxNumCuts,
maxNumCutsPerVariable, useQuantiles) stays OUT of the view - prior,
not data; a buildFromParent view carries its parent's grid. The x
ctor argument STAYS as the transient block a null map reads (the
fold-into-view alternative touches ~84 hand-edited construction
sites for no capability; declined). C API verdict: dbarts.h
UNCHANGED - all 7 touched entries carry SEXPs or dense double*, so
C_interface.cpp builds a stack dense view per call; no ABI
checklist, no stan4bart lockstep, no hash re-bake; typed C API
entries (a CSC-triple spelling, dbarts_sampler_setTestData-shaped)
recorded as a door - cost is a re-bake plus consumer lockstep with
no consumer asking. Neutrality is the gate: every currently-accepted
path stays bitwise (trio); the one capability opened (the door,
step 5) cannot perturb an accepted path.

1. Engine, view + builders: define PredictorSource with a
   sourceOf(j) accessor (null map = identity) so no caller
   dereferences a null map. ColumnStore::build takes (const
   PredictorSource&, grid spec, gatherColumns), keeping its two fill
   loops selected by columnSources == nullptr; buildMixed folds in;
   buildFromCsc is deleted, and BOTH its remaining call sites move
   to the view (createDataHandle, R_interface_bartcore.cpp:2549,
   and Sampler::setData, sampler.hpp:947). buildTest and
   buildTestMixed collapse to
   one buildTest(const PredictorSource&) - by inspection an identity
   map reproduces ownedTestValues, codeOffsets, empty
   sparseColumns/ownedTestCsc*, and the same quantize addresses, so
   the dense test path is bit-identical (gated in tests/cpp, step
   6).
2. Engine, SamplerOptions: the 8 typed ingestion fields collapse to
   one PredictorSource predictors; maxNumCutsPerVariable stays as
   the grid override. The Sampler ctors' 3-way dispatch becomes one
   build call; the nulling lines become options_.predictors = {}
   (sampler.hpp, and the Chain copies at chain.hpp:459/644/682).
3. Engine, test + mutation entrances: SamplerBase keeps ONE virtual
   bool setTestData(const PredictorSource&); setTestPredictors
   becomes a non-virtual inline building the dense view - no call
   site moves, the virtual surface shrinks by one.
   setPredictor/updatePredictor virtuals take the view too, dense
   convenience spellings kept as non-virtual inlines, and one shared
   predicate ACCEPTS only a null/identity map with NO CSC triple -
   refusing "CSC-valued" alone is a hole: a mapped dense view would
   be silently mis-consumed, every mutation kernel indexes
   values + k * n (sampler.hpp:1157/1209/1228/1265). Facade
   overrides of the view virtuals add using SamplerBase::setPredictor
   (etc.) so the dense inlines stay visible past name hiding. The
   alternative - mutation stays const double*, 2b scopes to
   creation + test - costs a second facade-virtual change later, and
   every such change forces --preclean rebuilds downstream;
   declined. Sampler::setTestData's leaf-covariate refusal reads
   sourceOf(), never the bare map.
4. Bridge: ParsedData carries one PredictorSource (plus one for the
   test container) replacing x/xIsSparse/xIsMixed/the CSC slots/
   mixedDenseValues/columnSources/cscReferenceCodes/categoryCounts;
   parseData fills the ~j map for a bare dgCMatrix so its three x
   flavors converge (xIsSparse/xIsMixed survive as parse-time
   booleans for refusal texts). optionsFromParsed publishes one
   field; parseTestContainer fills a view, installTestContainer
   passes it; createDataHandle's 3-way dispatch collapses.
   Validation stays keyed on the view and gets MORE uniform:
   validateCategoricalPredictors loses its storage-keyed
   if (data.xIsSparse) return early exit (a pure-CSC design is
   all-ordinal, the loop was already empty) and reads
   sourceOf(j) < 0; trainingCategoryBound, declaredCategoryCount,
   readDeclaredCategoryCounts, rawTrainingColumn,
   rawParsedTestColumn key on the same accessor.
   validateTestContainerAgainstStore and validateColumnValues are
   store-keyed already and do not move. Two spellings are
   load-bearing and stay keyed on the parse-time booleans, NEVER
   the map: parseData:816's categorical + dgCMatrix refusal and
   bartcore_setData:3284 (slice-0 memory safety). facade.hpp:556's
   pure-CSC predicate is RE-DERIVED, not renamed - a bare dgCMatrix
   now carries a ~j map, so the old spelling goes false (the
   per-column check at :565 happens to preserve the outcome, but
   the predicate must not silently flip).
5. R, the whole-matrix-on-mixed door: OPEN it. Refused today
   (R/bartcore.R ~127) only because the whole-matrix branch installs
   data@x <- x, replacing the container with a bare matrix; the
   engine and flat C API accept the mutation since extension (i).
   CONTAINER designs only: compute newX <-
   installPredictorColumns(data@x, NULL, seq_len(p), x) BEFORE the
   .Call (the 2a discipline) and assign on acceptance; the engine
   borrows the argument, so no pre-install/revert dance. The DENSE
   whole-matrix branch keeps its existing spelling untouched -
   installPredictorColumns drops dimnames where x[, columns] <-
   values preserves them, and by-name setPredictor/setCutPoints
   depend on that. Add a plural replaceSparseColumns splicing every
   named column into @i/@p/@x in ONE pass and route the
   multi-column case through it (the per-column loop would make
   whole-matrix O(p * nnz); one-pass is speed only, verified). Note
   in the docs step: replacing a sparse column densifies that
   column's storage permanently (nnz grows to n). DEPENDS ON FIX-0.
   The keep-shut alternative leaves a documented C API capability
   unreachable from R with no engine reason; the open cost is the
   caller materializes an n x p dense matrix, which column-granular
   mutation avoids; opened.
6. tests/cpp: collapsed builders bit-identical to the ones replaced -
   dense build via null-map view vs. pre-change build (codes, cuts,
   numCuts), dense TEST build via null-map view, and an
   all-negative-map view vs. the DENSE-EQUIVALENT build (the
   test_model.cpp:2189 fixture shape; buildFromCsc is deleted in
   step 1 and cannot be the oracle) - on both storage tiers, with a
   categorical column carrying declared K.
7. tinytest: the door's acceptance (whole-matrix setPredictor over a
   mixed container and a bare dgCMatrix, with 2a's save/re-create
   mirror gate), rollback leaving data@x untouched, the refusals in
   test-data-sparse.R:133-137 and test-data-mixed.R:171-179 flipped
   to acceptance. NO new in-file bitwise snapshots - the trio is the
   neutrality gate.
8. Docs in the SAME commit: sparse-columns.md extension (i) loses
   the whole-matrix clause (SUPERSEDED note gains the door's
   landing); man/dbarts.Rd drops "fixed at creation" for
   whole-matrix replacement; man/dbartsSampler-class.Rd mutation
   prose, including the sparse-column densification note (step 5);
   NEWS bullet for the door.
9. Out of scope, recorded: sparse-VALUED mutation (the view is its
   spelling; kernels not written); typed flat C API entries (the
   hash re-bake door); setData on a CSC/mixed store and container
   x.test through setData; per-observation mutation of CSC-backed
   columns; the grid spec in the view; per-column u8 widths; the
   streaming range kernel.

## Slice 2b verification

R CMD INSTALL --preclean (facade virtuals move - stale objects
bus-error); tests/cpp make clean build + local ASAN leg; full
tinytest; equivalence trio bitwise via the dedicated harnesses
(MANIFEST-current baselines) - the gate for the whole slice; air
format --check .; full local R CMD check (R, Rd, tests move); CI
incl. sanitizers green. No categorical-exact.R (no grid change).
