# data-ownership-2-ingestion

agent: opus (commit 2 touches the bridge parse path; commits 1, 3, 4 are
  R-surface, Sonnet-appropriate, but one owner keeps @x's contract coherent).
rng: neutral - frame ingestion reproduces today's categorical codes and
  quantization order byte-for-byte (equivalence gate). @x's storage vehicle
  changes; draws do not. Gate: equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds at EVERY commit; full tinytest from a preclean
  install.
window: plan 2 of docs/design/data-ownership.md (FROZEN 2026-07-11);
  builds on the plan-1 container (2e2b1c9). Plan 3 rewires mutation; plan 5
  ships the sparse-categorical kernel + test-side sparse.
budget: ~700-1100 lines across ~13 files. No engine (.hpp) change; the
  bridge grows one internal parseData branch; dbarts.h untouched.
bench: ingestion is construction-time; no partition/draw path is touched.
  One confirmatory bench-sampler compare at the end (build + sample time, no
  regression). Sub-ms setPredictor metrics carry 6-8% run-to-run noise -
  re-run a single-run flag before belief.

## Goal

data.frame/formula input ingests without ever retaining a cbound double
matrix: @x holds the ingested columnar source, and a transient contiguous
block is assembled only inside the create/setData call, byte-identical to
today's makeCategoricalModelMatrix output. sparseFactor() is exported as the
CSC-categorical wrapper (recognized, refused until plan 5). extract(sampler,
"predictors") is the public materialization verb and the getTrees-replay
source. Every internal data@x consumer is routed or retained by verdict.

## Baseline (what already holds at plan 1)

- build/buildMixed/buildFromCsc take raw as CALL arguments, quantize into
  owned codes, retain nothing unflagged ([[data.hpp:550@7842e2f2]], [[data.hpp:601@7842e2f2]], [[data.hpp:702@7842e2f2]]).
- @x is a matrix, dgCMatrix, or dbartsMixedMatrix; @varTypes carries
  ORDINAL/CATEGORICAL ([[A_class.R:394@7842e2f2]]; data.hpp ColumnType). Dense
  categorical already rides buildMixed's per-column dispatch.
- factors = "categorical" is ALREADY the default in dbartsData / dbarts /
  bart2 / rbart ([[data.R:273@7842e2f2]], [[dbarts.R:188@7842e2f2]], [[bart.R:377@7842e2f2]], [[rbart.R:47@7842e2f2]]). The
  dummy -> single-categorical model change landed with the categorical
  engine, NOT here; plan 2 preserves it bitwise. "indicators" remains the
  opt-out. No new posterior decision.
- getTrees replay, setState cross-grid, decodeCategoricalSplits already read
  data@x as their raw/metadata source ([[dbarts.R:1085@7842e2f2]], [[dbarts.R:1104@7842e2f2]], [[dbarts.R:1233@7842e2f2]], [[dbarts.R:1238@7842e2f2]]).
- PROT_DATA is the GC anchor and the call-time raw source (design, resolved).

## Ingestion type mapping (end to end)

    frame column            -> store path            -> engine kind    status
    numeric / logical       -> dense (build/Mixed)   -> ordinal        ships
    ordered factor          -> dense, codes as u_-1  -> ordinal        ships
    unordered factor        -> dense, codes as u_-1  -> categorical    ships
    matrix col (poly())     -> dense, spliced        -> ordinal        ships
    I() sparseVector/dgC    -> CSC (buildMixed)      -> ordinal        ships
    sparseFactor()          -> CSC categorical       -> (plan 5)       REFUSE

Factor/ordered columns supply INTEGER factor codes; the bridge coerces them
to the transient double code-vector (as.integer - 1), the exact values
makeCategoricalModelMatrix writes today. Numeric columns feed their REAL
pointer. buildMixed already accepts dense categorical; only the SOURCE of the
dense block moves from a retained R cbind to a transient bridge assembly.

## Named decision S-CAT (sparse categorical scope)

sparseFactor() lands as R SURFACE ONLY: the class, its export, Rd, print, and
dbartsData recognition. CSC-categorical CODE STORAGE and the membership
kernel are plan 5 (the container has no CSC-categorical path; parseData's
all-ordinal-sparse guard, [[R_interface_bartcore.cpp:414@7842e2f2]], [[R_interface_bartcore.cpp:419@7842e2f2]], already refuses
it). A sampler/data built on a sparseFactor column REFUSES at construction
with an explicit "sparse categorical predictors are not yet supported"
message. I() sparse ORDINAL is unaffected and keeps working.

## Formula-path rework

makeModelMatrix ([[data.R:286@7842e2f2]]) stops producing a retained cbound matrix.
makeCategoricalModelMatrix returns the columnar container (below) instead of
do.call(cbind, columns); makeModelMatrixFromDataFrame ("indicators", makeind,
the dummy path) is UNCHANGED - it still returns a plain matrix. The
model.frame -> model.matrix double detour is gone: term.labels/factor.levels/
varTypes ride on the container as today. Backward compatibility: matrix and
"indicators" input keep a plain-matrix @x; only dense-frame-with-factors @x
changes class (matrix -> container). Downstream matrix-interface users are
untouched.

## The columnar source container

Reuse dbartsMixedMatrix (mixedMatrix.R): it already carries dim/dimnames/
[/as.matrix + the term.labels/varTypes/factor.levels/drop attributes and
satisfies the inherits(@x, "dbartsMixedMatrix") guards ([[model.R:167@7842e2f2]], [[model.R:210@7842e2f2]]).
Extend it to keep DENSE columns as a per-column list (not cbound) with a
categorical-capable map; as.matrix materializes on demand (the getTrees /
partialDependence path). Pure-dense frames now yield a container, not a
matrix. Bridge: parseData gains a container branch that assembles the
transient contiguous dense block (coercing factor codes) + gathers CSC slots,
then dispatches to build/buildMixed exactly as today - byte-identical codes.

## extract(sampler, "predictors")

New extract.dbartsSampler(object, type = "predictors"); S3method(extract,
dbartsSampler) (extract is already the S3 generic, [[generics.R:3@7842e2f2]] - no new
generic). Returns the NUMERIC code matrix (factors as integer codes),
matching the historical @x contract and the getTrees-replay source; pure R
over @x (as.matrix of the container, or the matrix itself), no engine call in
plan 2. Plan 3 reroutes it to owned raw for mutable columns.

## data@x consumer census + drop verdict

Census: ~40 data@x reads across R (bartcore.R, dbarts.R, partialDependence.R,
plotTree.R, bart.R, rbart.R, model.R, utility.R, updatePredictorPerObserv...).
Kinds: metadata (ncol/nrow/colnames/factor.levels - served by the container's
dim/dimnames/attrs), value reads (partialDependence quantiles, plotTree,
getTrees replay, setState raw, estimateSigmaFromLinearModel), and the
plan-1 mutation interim write-back ([[bartcore.R:94@7842e2f2]], [[bartcore.R:171@7842e2f2]]; updatePredictor...:87).

VERDICT: @x is RETAINED, its role narrowed to "the ingested source". Dropping
it does NOT close cleanly - it is the PROT_DATA GC anchor, the call-time raw
source for re-cut/setState/replay, and the materialization root for
partialDependence/plotTree; the engine deliberately keeps no raw to
reconstruct from. Value consumers route through extract()/as.matrix (commit
1); metadata consumers read the container directly. The mutation interim
write-back stays until plan 3. Design mandate ("see if you can drop it")
answered: keep, as the design's fallback ("keep a snapshot slot only if a
consumer genuinely needs a materialized matrix") applies.

## sparseFactor()

S4 class: slots i (0-based rows), values (integer level codes), levels
(character), reference (implicit-zero level), length. Constructor
sparseFactor(x, levels, reference) validates reference %in% levels, codes in
1..K, i in range. print shows levels/reference/nnz. Rd (man/sparseFactor.Rd),
NAMESPACE export. isSparseDataFrameColumn ([[mixedMatrix.R:10@7842e2f2]]) recognizes it;
dbartsData refuses per S-CAT.

## Commits

1. extract(sampler, "predictors") + consumer routing (R). Add the method +
   S3method; route getTrees replay / decodeCategoricalSplits / setState raw /
   partialDependence / plotTree value-reads through extract or as.matrix so
   a later @x-class change is invisible to them. Gate: full tinytest;
   equivalence 22/22. Abort: any draw or snapshot moves - revert, the routing
   is not neutral.
2. Columnar container + frame-direct ingestion (R + bridge). Extend
   dbartsMixedMatrix to un-cbound dense; makeCategoricalModelMatrix /
   data.R / formula path store the container in @x; parseData container
   branch assembles the transient dense block + CSC dispatch. Gate: tests/cpp
   codes bitwise vs the cbind path; structural @x-class assertions in
   test-data-mixed / test-data-categorical / test-makeModelMatrix updated to
   the container contract (NOT RNG snapshots); equivalence 22/22; full
   tinytest. Abort: codes differ from the cbind reference - stop, the
   transient assembly is wrong.
3. sparseFactor() wrapper + refuse (R). Class, constructor, validation,
   print, Rd, export; dbartsData recognition + construction refusal (S-CAT).
   Gate: new tinytest (construct/validate/print/refuse); equivalence 22/22;
   full tinytest. Abort: a sparseFactor column reaches the engine - the
   refusal leaks.
4. Docs + landing. dbartsData.Rd (frame contract, factors), extract Rd
   ("predictors"), the design landing note (steps 1-2 style) recording the
   @x-source verdict and the S-CAT boundary. Gate: R CMD check man; full
   tinytest. Abort: none (docs).

## Verification

- Full tinytest per commit from a preclean install (R CMD INSTALL . ;
  --preclean after commit 2's header/bridge edits).
- Equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at every commit (the
  frozen-default gate; any deviation = defect, stop).
- tests/cpp bitwise codes: container-fed transient block == the cbind
  reference (commit 2); delete stale test binaries after the bridge edit.
- One bench-sampler compare at the end vs the current recorded baseline
  (construction + sample time; no regression expected, ingestion is not on a
  hot path). Sub-ms setPredictor noise caveat above.
- rchk on the next scheduled run (commit 2 touches the bridge).
- dbarts.h unchanged; no stan4bart lockstep delta (bridge growth is internal
  .Call only, the plan-1 precedent).

## Open questions for VD

- extract(sampler, "predictors") return VALUE for factor columns: the plan
  returns the NUMERIC code matrix (replay-canonical, historical @x contract).
  If the user-facing verb should instead return the original frame (factors
  as labels), the plan adds a form/type distinction. Default: proceed with
  codes. (Factor semantics and the @x-source verdict are resolved above, not
  reopened.)

## Landing notes

Commit 1 = 6dda9a7. extract(sampler, "predictors") + consumer routing.
rawPredictorMatrix added as a bridge-signal helper: it returns the dense
matrix, or the NULL a sparse source always signals - the bridges already
dispatch dense/sparse on that NULL, so routing through it costs no new
branch. decodeCategoricalSplits was reviewed and ruled a METADATA consumer
(it depends on the factor.levels attribute, not predictor values), so it
was left unrouted. Neutral: equivalence 22/22 identical, tinytest 2728 no
regen.

Commit 2 = c82f954. dbartsMixedMatrix's dense flavor stores columns as a
per-column list (factors keep their integer codes) instead of a cbound
matrix; the bridge assembles the transient contiguous block inside
parseData, owned by the parse result. DEVIATION adopted: the sparse-bearing
mixed flavor keeps its resident cbound dense block - buildMixed retains
borrowed per-column slices of it, and a transient block would dangle once
parseData returns; plan 5 revisits when the CSC-categorical kernel lands.
DEVIATION adopted: predictor-mutation write-backs got a container-aware
interim, assignIntoPredictorSource (frame-built samplers now carry
containers in data@x, not matrices); plan 3's mutation rewire replaces it.
One structural test change: test-dart-mixed-columns assumed anyNA on a list
recurses into elements, which it does not - the read now materializes
through the extraction verb first. Gates: tests/cpp bitwise
transient-assembly check vs the cbind reference; equivalence 22/22;
tinytest 2728 no regen; a frame fit matches a hand-coded-matrix fit draw
for draw.

Commit 3 = 20c90d3. sparseFactor class/constructor/show/length, Rd, export.
Refusal lives at sparseColumnSlices, the single choke point every sparse
column passes through, plus a formula-path pre-scan - model.frame would
otherwise die on a bare S4 column with an opaque type error before the real
refusal message fires. The length method is REQUIRED, not cosmetic:
without it NROW cannot size a sparseFactor column, so the class could not
be a data.frame column at all. Gates: equivalence 22/22; tinytest 2763 (35
new) no regen; Rd tooling clean.

Commit 4 = this commit (docs). extract.dbartsSampler.Rd and dbartsData.Rd
finished (frame contract, factor mapping, the no-retained-matrix note);
dbarts.Rd picked up the same touch (it repeats parts of the contract).
Gates: R CMD INSTALL clean; checkRd/undoc/codoc clean; full tinytest 2763,
no regen; bench-sampler compare vs bench-sampler-32fc7c8.csv shows no new
regression.

Resolutions: the open extract-return-form question above is RESOLVED to
the numeric code matrix - factors as integer codes, matching the historical
@x contract and getTrees' replay source. The @x verdict (KEEP, role
narrowed to "the ingested source") and the S-CAT boundary (sparse
categorical storage deferred to plan 5) both stood as planned; neither was
reopened.
