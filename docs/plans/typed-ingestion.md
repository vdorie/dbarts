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
