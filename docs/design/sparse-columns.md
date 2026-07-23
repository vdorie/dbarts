# Sparse columns

Prototype study settling the representation question for sparse predictor
columns (core-generalization.md data model; kernel-vocabulary.md planned
addition 4; public-surface.md section 7), including whether the engine
needs order-preserving partitions. Prototyped 2026-07-04 in
benchmarks/kernels/sparse.c; the landing plan and landing notes follow the
prototype sections.

## Problem

The hot per-node operation partitions an arbitrary, scrambled index set by
a column's quantized codes: today a random read into a dense u16 array per
member observation. Dense storage costs 2 bytes per row per column, which
is the pain point for wide mostly-zero designs (one-hot expansions, text
features, genomics); the per-observation work is already index-set-driven,
so sparsity's win is memory first and speed only incidentally.

## Candidates

- **dense** (baseline): n u16 codes, misc_partitionIndices, NEON.
- **rank**: bitmap of nonzero rows (1 bit/row), cumulative popcount per
  64-row word (0.5 bit/row), packed nonzero codes (2 bytes/nonzero).
  code(i) = bit clear ? zero code : nzCodes[rank(i)], rank in O(1) from
  the word count plus a masked popcount. Keeps the unstable two-pointer
  partition: nothing about the engine's index handling changes.
- **bsearch**: CSC (sorted nonzero rows + codes), code(i) by binary
  search. The naive layout.
- **merge**: CSC walked in tandem with the index segment, requiring the
  segment SORTED - the "order-preserving partition" hypothesis: if every
  partition in the engine were stable, segments would stay ascending and
  sparse columns could stream. Prototyped as a stable scratch partition.

## Results

Apple M-series (arm64), n = 262144, nonzero codes uniform 1..250, cut at
the nonzero median, zeros left; ns per segment element, scrambled segments
(sorted for merge). f = nonzero fraction, k = segment size.

| f    | k      | dense | rank  | bsearch | merge |
|------|--------|-------|-------|---------|-------|
| 0.50 | 262144 | 2.03  | 8.48  | 50.4    | 3.97  |
| 0.50 | 16384  | 0.89  | 3.11  | 49.8    | 10.8  |
| 0.50 | 1024   | 0.71  | 1.14  | 46.5    | 53.0  |
| 0.10 | 262144 | 1.31  | 1.87  | 27.4    | 1.53  |
| 0.10 | 16384  | 0.76  | 0.95  | 27.4    | 6.98  |
| 0.10 | 1024   | 0.46  | 0.61  | 25.6    | 13.2  |
| 0.05 | 262144 | 0.98  | 1.17  | 23.9    | 1.24  |
| 0.05 | 16384  | 0.76  | 0.92  | 23.8    | 4.91  |
| 0.05 | 1024   | 0.39  | 0.65  | 22.3    | 6.58  |
| 0.01 | 262144 | 0.76  | 0.72  | 22.7    | 1.00  |
| 0.01 | 16384  | 0.78  | 0.75  | 22.4    | 1.42  |
| 0.01 | 1024   | 0.38  | 0.63  | 20.9    | 2.65  |

Memory per row per column: dense 2 bytes; rank 0.1875 + 2f bytes (0.21 at
f = 0.01, 0.29 at f = 0.05, roughly 7-10x smaller where sparsity is
real); CSC 6f bytes.

## Conclusions

- **Representation: rank-bitmap.** In the sparsity regime that motivates
  the feature (f <= 0.05) it runs at 0.95-1.7x the dense kernel's time -
  at parity on large segments, modestly behind on small cache-resident
  ones - while storing the column in a tenth the memory, and it is a
  drop-in per-column layout: same unstable partitions, same index
  handling, O(1) random code access for the tree-descent and
  per-observation paths. The prototype kernel is scalar against dense's
  NEON, so there is headroom, but none is needed to justify the layout.
- **Order-preserving partitions are NOT required** - the question the
  prototype existed to settle. The merge layout is competitive only on
  root-sized segments and collapses on small ones (the CSC scan pays
  O(touched nonzeros) regardless of k), and adopting it would force every
  partition engine-wide to become stable: scratch traffic on the dense
  fast path, and a changed within-leaf accumulation order that breaks
  every bitwise baseline. Rejected.
- **Binary search is out** (25-60x): no niche once rank exists.

## Integration sketch (future work, phase-4 data model)

- ColumnStore (or its BartData successor) gains a per-column storage
  kind; the partition dispatch in Tree::partitionChildren and the code
  accessors (findBottomNodeForRow/ForObservation, quantize paths) branch
  per column exactly as ordinal/categorical/pooled already do. Kernel
  vocabulary gains misc_partitionIndicesSparse over the rank layout with
  a scalar reference first (kernel-vocabulary.md conventions).
- Cuts are computed over the full value distribution, zeros included;
  the zero code is codeFor(0.0) at build. Missing values compose as
  explicit stored entries carrying the reserved NA code, so MIA needs
  nothing new.
- Ingestion: accept Matrix::dgCMatrix in dbartsData, choosing rank
  storage per column below a density threshold (say f < 0.2, where the
  table crosses); dense numeric input never changes representation.
- Mutation: values of stored nonzeros can change in place (the codes
  array re-quantizes); a nonzero-pattern change rebuilds the column's
  bitmap and rank index in O(n / 64 + nnz), the same cost class as a
  dense column's re-quantize. setData rebuilds wholesale. Views
  (buildFromParent) can gather into either layout; simplest is to
  densify views until sparse proves itself there.

## Landing plan (2026-07-04)

The prototype's integration sketch made concrete. Scope: rank-bitmap
storage for sparse ordinal columns, entered through Matrix::dgCMatrix
ingestion; codes stay u16 throughout. Per-column code widths (u8 hot
layers) are NOT in this landing: the machinery the two share is the
per-column storage dispatch, which this delivers, while width variants
multiply the SIMD kernel surface and stand alone as future work.

### Data model

ColumnStore keeps one contiguous codes vector but gains per-column
offsets into it (codeOffsets[j]; j * n for dense-matrix builds, packed
among densified columns for CSC builds, unused for rank columns), so the
dense layout, arithmetic, and kernels are unchanged where they run today.
Rank columns live in SparseColumnData slots: bitmap words (ceil(n/64)),
per-word cumulative popcounts (u32; dgCMatrix row indices are int, so n
fits), packed nonzero codes in ascending-row order, the zero code, plus
borrowed CSC value/row slices for re-quantization. The rank build
scatters through the bitmap so it never assumes sorted or unique row
indices; a total-popcount != nnz check rejects duplicates.

A dgCMatrix build (buildFromCsc) stores columns with nnz/n <= 0.2 (the
prototype table's crossing region, internal and tunable) as rank
bitmaps and densifies the rest into packed dense codes - densified
columns run the existing SIMD kernels bitwise-identically to a dense
build of the same values, which is also the equivalence test. Raw
doubles are never materialized: the store's x stays null and rawColumn
serves nothing, so linear-leaf designations are refused up front.

Cuts are computed over the full logical distribution: uniform grids
min/max over nonzero values folding in 0.0 when implicit zeros exist
(nnz < n), quantile grids over the sorted unique nonzeros plus that
zero. Both reproduce the dense builders' arithmetic exactly, so a CSC
build and a dense build of the same matrix carry identical cut points
and identical codes (component-tested). Missing values arrive as stored
NaN entries (the Matrix convention), quantize to the reserved codes,
and set hasMissing: MIA composes with no new machinery.

### Kernels and engine dispatch

misc_partitionIndicesSparse joins the vocabulary as a dispatched pointer
with the scalar reference installed at every ISA level (kernel-
vocabulary.md planned addition 4): the two-pointer partition over
(bits, wordRanks, nzCodes, zeroCode), O(1) rank access per element. The
engine uses it for every sparse ordinal partition including the root -
misc_partitionRange assumes identity index content, which holds only on
the dense path - so a streaming range variant remains headroom. Columns
with missing values take an engine-side scalar MIA sibling
(partitionIndicesSparseMIA), mirroring partitionIndicesMIA over the rank
accessor. Tree descent (findBottomNodeForObservation) and state-restore
partition checks (setPartitionsFromOrderedIndices) go through a
storage-aware codeAt(j, i). Test codes stay dense row-major and x.test
stays a dense matrix: quantization only reads cut points.

Views densify (the sketch's recommendation): buildFromParent gathers
parent codes through codeAt into a fully dense child store, so the
data-handle and xbart paths compose without sparse-specific code; folds
of a wide sparse design pay a transient dense-codes copy.

### Mutation surface and state

The raw-x mutation surface - setPredictor, updatePredictor, the
per-observation sessions, setData - is refused for CSC-built stores in
v1 ("sparse predictors fix the design at creation; make a new sampler
instead", the grouped-random-effects precedent). The existing view
guard already fires on x == null; the bridge distinguishes the two
cases for the message and, unlike views, ALLOWS setState and
setCutPoints on CSC stores: quantizeColumn is storage-aware, and
getPointer() re-creation from the stored dbartsData (which holds the
dgCMatrix) rebuilds the sampler, so save/load works exactly as it does
for grouped fits. In-place nonzero-value mutation and pattern rebuilds
(the sketch's O(n/64 + nnz) path) wait for a consumer.

### R surface

dbartsData accepts a dgCMatrix as x (the x slot widens to ANY with the
class checked in validity; x.test stays dense). All columns are ordinal
- factors do not ride numeric sparse matrices - and the sigma estimate
falls back to sd(y). dbarts, bart2, rbart_vi, and xbart inherit
acceptance through dbartsData; node.prior = linear() with sparse x is
refused when columns are designated. The bridge's parseData reads the
i/p/x slots directly (borrowed for the sampler's lifetime, like dense x
today) and hands them to the engine through SamplerOptions. Matrix goes
to Suggests.

### Stages and gates

1. Kernel: misc_partitionIndicesSparse + simd.c registration +
   kernel-vocabulary.md; correctness vs the dense kernel in tests/cpp.
2. Engine: ColumnStore storage dispatch, buildFromCsc, cut/quantize
   variants, tree.hpp dispatch, view densification; component tests
   (CSC-vs-dense cuts/codes equality, partition membership equality,
   densified-tier bitwise fit equality, rank-tier end-to-end recovery,
   MIA, state round trip).
3. R and bridge: parseData CSC branch, refusal wiring, dbartsData
   acceptance, tinytest test-data-sparse.R, docs and NEWS.

Gates per the standing checklist; equivalence must report identical
draws (dense paths are layout-refactored only), speed must hold at
baseline on every metric.

## Landing notes (2026-07-04)

Landed per the plan; deltas and specifics:

- misc_partitionIndicesSparse lives in src/misc/partition.c as a plain
  exported function (no dispatch pointer until a SIMD variant exists),
  mirroring the scalar two-pointer body exactly. MIA columns take the
  engine-side partitionIndicesSparseMIA in tree.hpp.
- ColumnStore: codeOffsets indexes the single packed codes vector
  (j * n for dense builds, packed among densified columns for CSC
  builds); rank columns live in SparseColumnData slots keyed by
  sparseSlot, with borrowed CscColumnSlice values/rows per column for
  re-quantization. codeAt(j, i) is the storage-aware accessor used by
  tree descent, restore validation, and view gathering. The bitmap
  build scatters, so sorted rows are not assumed by the engine - but
  the bridge validates rows strictly increasing, unique, and in range
  ("malformed sparse predictor matrix") because the engine trusts its
  input beyond that.
- Cut construction: fillCutsUniformlyCsc folds an implicit zero into
  the min/max seed; quantileGridForCscColumn collects stored non-NA
  values plus one zero when implicit zeros exist. Both reproduce the
  dense builders bitwise (component-tested under both cut modes), so
  an all-densified dgCMatrix fit draws bitwise-identically to the
  dense-matrix fit of the same values - that is the equivalence-style
  component test (testSparseEndToEnd part 1).
- setState snapshots sparseColumns alongside codes for its rollback;
  save/load works through getPointer() re-creation from the stored
  dbartsData holding the dgCMatrix (the grouped/n.cuts precedent),
  with the bridge's view refusal split so setState/setCutPoints stay
  open on CSC-built samplers (refuseViewSamplerOnly) while the raw-x
  surface refuses with "sparse predictors fix the design at creation".
- R surface: dbartsData's x slot widened to ANY with validity checking
  matrix-or-dgCMatrix (a slot union cannot name a Suggests class at
  load); the sigma estimate falls back to sd(y - offset); bart pins
  missing = "error" so NA-bearing sparse fits go through bart2/dbarts.
  Matrix sits in Suggests; the bridge reads the i/p/x slots via
  Rf_getAttrib without touching the Matrix namespace.
- xbart/data-handle: bartcore_createDataHandle builds the parent store
  from CSC and fold views densify through codeAt - no other change.
- Component tests: testSparseKernel (dense-reference partition at
  boundary cuts), testSparseColumnStore (bitwise cuts/codes vs the
  dense builder, tier split, view densification), testSparseEndToEnd
  (densified-tier bitwise match; rank-tier signal recovery with MIA,
  sigma 0.295 vs truth 0.3), testSparseStateRoundTrip (bitwise
  continuation, 2 chains). tinytest test-data-sparse.R (21 asserts)
  covers ingestion, fits, MIA, refusals, save/load, xbart.
- Gates: tests/cpp 37 tests green; tinytest 2285 across 63 files, 0
  failures; equivalence identical draws on all nine scenarios (dense
  paths are layout-refactored only); speed at baseline; R CMD check
  --as-cran Status OK.

Still open, by design: sparse x.test, rbart_vi and linear-leaf
support, per-column u8 code widths, a streaming range kernel for
root-sized segments, and any public exposure of the density threshold.
(In-place nonzero-value mutation and pattern rebuilds LANDED as extension
(i), 2026-07-22 - see the section below.)

## Planned: mixed dense and sparse predictors (Vincent, 2026-07-04)

Requested surface: a data frame holding ordinary numeric and factor
columns alongside sparse columns - Matrix::sparseVector entries and
matrix-valued columns including dgCMatrix (both enter data frames
I()-wrapped, matrix columns natively) - through the existing
data-frame interfaces, so one input mixes both worlds. The v1 landing
above is deliberately all-or-nothing per input object; the per-column
machinery it introduced is the enabler here, and the restrictions v1
imposes globally are per-SOURCE facts, so a mixed store relaxes them
where the source allows:

- Engine: ColumnStore::buildMixed(denseX, cscTriple, columnSources,
  ...) assigns each column either a dense column-major slice
  (quantized exactly as build() does; categorical allowed; rawColumn
  serves it) or a CSC slice (exactly as buildFromCsc: threshold
  chooses rank bitmap vs densified codes; ordinal only). Nothing
  downstream changes: codeAt, partitionChildren, the cut builders, and
  quantizeColumn already dispatch per column. Dense-backed columns in
  a mixed store keep categorical splits (factors compose with sparse
  features, which v1 cannot do at all) and are designatable as linear
  leaf covariates (the facade check moves from store-wide to
  per-designated-column via rawColumn(j) != null, the
  createSamplerOverStore precedent).
- Ingestion: the model-matrix builders (makeCategoricalModelMatrix,
  makeModelMatrixFromDataFrame) gain a sparse pass-through - a
  sparseVector column contributes one ordinal column, a dgCMatrix
  column contributes its columns, dense columns behave exactly as
  today including the factors= treatment. Primary surface is the
  data-frame/xy path; the formula path routes through the same
  builders but model.frame's handling of S4 columns needs validation
  during implementation.
- Container: dbartsData@x (already ANY) holds an internal mixed
  container - the dense columns as one numeric matrix, the sparse
  columns cbound into one dgCMatrix, plus a source map and dimnames,
  with ncol/colnames/nrow served so the R-side consumers (varcount
  naming, plotTree, partialDependence) keep working. A plain R object,
  so getPointer() re-creation and save/load work unchanged.
- Bridge: parseData recognizes the container and passes the dense
  pointer, the CSC slots, and the per-column source map through
  SamplerOptions (internal, freely extensible); validation as v1 plus
  map consistency.
- Test data: stays dense in v1 of mixed too - a test data frame
  expands through the same builders to a dense matrix over all
  columns, which validateXTest already produces.
- Mutation: store-wide refusal stays while any CSC-backed column
  exists (the transactional paths are column-granular, so a later
  relaxation to dense-backed columns is mechanical, but it needs a
  consumer). rbart_vi keeps its refusal (the R loop predicts over
  data@x) until the in-core path takes sparse.

Order of work when picked up: buildMixed + component test (mixed store
bitwise vs a two-store reference), the container + ingestion
pass-through, bridge plumbing, then the per-column linear-leaf and
categorical relaxations with tests. No kernel work is involved.

## Mixed input landing notes (2026-07-04)

Implemented per the plan above; deltas and specifics:

- Engine: ColumnStore::buildMixed(denseValues, cscTriple, columnSources,
  ...) with columnSources per column: nonnegative names a dense-slice
  column, negative the complement (~) of a CSC column. buildFromCsc is
  now a thin wrapper (an all-complement map), so the pure-CSC path runs
  the same code. Per-column raw access went through a new
  denseSourceColumn(j) (the x matrix on dense builds, mixedRawColumns[j]
  on mixed builds, null otherwise) and the quantize/cut dispatch through
  columnIsCscBacked(j); every reader of x + j*n in data.hpp rewired.
  rawColumn(j) serves dense-backed mixed columns, which gives linear
  leaves and view raw-gather for free; buildFromParent now gathers
  through parent.rawColumn per column, skipping columns the parent
  cannot serve (the facade then refuses those designations), and
  inherits parent standardization constants when the parent is itself a
  view. SamplerOptions grew mixedDenseValues/columnSources; the facade
  refuses linear leaves per designated CSC-backed column instead of
  store-wide.
- I() does not survive on S4 objects and data.frame(...) rejects
  sparseVector/dgCMatrix arguments (NROW sees 1); columns must be
  assigned into an existing frame (df$s <- sv). model.frame() rejects S4
  columns outright, so - like the v1 dgCMatrix - mixed input enters
  through the x/y interface only. Documented in dbarts.Rd/bart.Rd.
- Container: class "dbartsMixedMatrix" (R/mixedMatrix.R), a plain list
  of dense matrix (NULL when the frame is all-sparse), one cbound
  dgCMatrix, a 1-based signed map (+k dense, -k sparse; R's -k IS the
  engine's ~(k-1), so the bridge maps almost verbatim), and the column
  names. S3 dim/dimnames/[/as.matrix registered; [ keeps the container
  (and builder attributes) on row subsetting, densifies on column
  selection with matrix drop semantics. as.matrix makes validateXTest
  work unchanged (its !is.matrix coercion densifies containers).
- Builders: makeCategoricalModelMatrix passes sparse columns through as
  ordinal (sparseVector one column, dgCMatrix its columns).
  makeModelMatrixFromDataFrame expands dense input columns ONE AT A
  TIME through the existing C builder - it treats columns independently,
  so blocks and drop entries match a whole-frame call - and splices
  sparse columns in place with all-FALSE drop entries, so the recorded
  drop pattern replays over a fully dense test frame.
- dbartsData: the xy data-frame branch skips complete-case dropping for
  containers (the sparse-branch precedent); the NA validation gains a
  container branch checking both parts. dbartsData validity accepts the
  container; estimateSigmaFromLinearModel's !is.matrix fallback covers
  it.
- Bridge: parseData gained the container branch (map validated per
  entry against the parts; the dgCMatrix structure validation is now
  the shared parseCscMatrix); categorical columns must be dense-backed
  ("sparse predictor columns must be ordinal");
  validateCategoricalPredictors reads through rawTrainingColumn;
  createDataHandle builds mixed parents; setData refuses mixed
  replacement data. refuseViewSampler needed no change (mixed stores
  are builtFromCsc, so the sparse message applies).
- resolveLeafCovariates allows containers but refuses sparse-backed
  designations ("leaf covariates must be dense columns ..."); rbart_vi
  and the setPredictor guard refuse via the existing !is.matrix checks.
- Tests: component testMixedColumnStore/EndToEnd/LinearLeaves/
  StateRoundTrip (41 total now) - the store and both samplers (constant
  and linear leaf) verify BITWISE against a dense build of the same
  values when all CSC columns densify; tinytest test-data-mixed.R (37
  results). Suite 2322 across 64 files.
- Observed while landing, pre-existing and out of scope: predict() on a
  bart2 fit after saveRDS/readRDS fails ("bartcore function called on
  NULL external pointer") for DENSE inputs too - bart2 does not store
  sampler state. Sampler-level save/load (storeState + saveRDS) works
  and is what the tests cover.

## Extension (i): sparse-column in-place mutation (landed)

The deferred in-place mutation from the integration sketch, made concrete.
Scope: setPredictor and updatePredictor (the whole-matrix and subset-column
surfaces, transactional and forced) now mutate CSC-backed columns, lifting the
store-wide refusal the v1 landing imposed. The motivating consumer is IRT-style
between-sweep mutation on a sparse design (a dense latent covariate mutated
while sparse item indicators sit in the store, and mutation of the sparse
columns themselves).

- The mutation surface hands the engine a DENSE column even for sparse storage
  (the bridge/flat C API re-quantize from a dense matrix). The nonzero pattern
  is {i : value != 0}, so a stored NaN (missing) stays stored and the minimal
  pattern a dense equivalent carries is produced - codes stay bitwise identical
  to a dense build of the same values. When the pattern is unchanged the
  nonzero values re-quantize IN PLACE (rank tier: nzCodes + zeroCode over the
  fixed bitmap; densified tier: the codes segment); when it changes the rank
  bitmap and index REBUILD, O(n/64 + nnz). ColumnStore::mutateCscColumnFromDense
  dispatches; setPredictors/setColumns route CSC-backed columns to it.
- The storage tier (rank vs densified) is FIXED at build and never flips on
  mutation, keeping the codeOffsets/codes layout and rank-slot assignment
  stable (a rank column whose nnz grows past the threshold stays rank -
  correct, just less memory-efficient; codes match dense regardless of tier).
- Owned re-quantize sources: a CSC column borrows R's dgCMatrix i/x slices at
  build; the first mutation copies the new nonzeros into engine-owned
  ownedCscRows/Values[j] and repoints the column's slice, since the borrow no
  longer reflects the live values (setCutPoints and state restore re-quantize
  from the slice). cscColumnOwned[j] tracks the borrowed-vs-owned split.
- Transactional roll-back mirrors the dense paths column-granularly: the subset
  path (updatePredictor) snapshots a CscColumnRollback per touched CSC column
  (source descriptor, rank storage or codes segment, owned buffers, ownership
  flag) alongside the dense journal for dense columns; the whole-matrix path
  (setPredictor) bulk-snapshots sparseColumns/sources/owned buffers when
  builtFromCsc, beside the moved codes. On reject both restore byte-for-byte
  and repoint owned slices (a moved-then-restored buffer need not sit at the
  pre-change address; a borrowed slice keeps its valid R pointer).
- Out of scope, refused by name: per-observation replacement of a CSC-backed
  column (updatePredictorPerObservation writes cells one at a time, which rank
  storage cannot take without an O(nnz) shift per cell - replace the whole
  column with updatePredictor; dense-backed columns of a mixed store, the IRT
  latent case, stay open); whole-matrix and per-observation replacement of a
  sparse design at the R5 level (data@x maintenance); mixed-container mutation
  at the R5 level (the container's sparse block is not yet rebuilt R-side - the
  engine and flat C API accept it); setData whole-data replacement of a CSC
  store; and CSC-backed categorical mutation (R never builds one). Save/load
  after a mutation still re-creates from the stored dgCMatrix, so data@x must
  reflect the mutation R-side (the dgCMatrix path maintains it via
  installPredictorColumns; the container path is deferred).

Gates: tests/cpp testSparseMutation (store codes bitwise-match the dense
builder after pattern-preserving and pattern-changing mutation across both
tiers; densified-tier draws bitwise-match dense after mutation; a rejected
transactional update restores codes/cuts/rank-storage/tree-fits
byte-for-byte). tinytest test-data-sparse.R covers the dgCMatrix column
mutation and the remaining refusals; test-data-mixed.R the mixed-container
refusal. Full suite green, all pre-existing snapshots intact (rng-neutral).

## Status

LANDED 2026-07-04 (plan and landing notes above); mixed dense/sparse
input LANDED 2026-07-04 (plan and landing notes above); extension (i)
sparse-column in-place mutation LANDED 2026-07-22 (section above).
