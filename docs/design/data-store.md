# Data store

The standing reference for the C++ predictor data layer: `ColumnStore` in
`src/bartcore/data.hpp` and the sampler transaction that mutates it
(`src/bartcore/sampler.hpp`). Required reading before any data-adjacent
work. This note documents the cross-cutting invariants and contracts the
per-function Doxygen comments cannot carry; where it cites `file:line`,
the anchor is HEAD at the time of writing and should be treated as a
starting point, not a guarantee.

The layer's one governing idea: the engine owns quantized codes, not raw
predictors. It borrows raw doubles only for the duration of a build or a
re-quantize call, then releases them. Trees consume codes; raw values
survive only where a downstream consumer genuinely needs them (re-cutting
an unchanged column, leaf models that read covariates). The whole design
- the cut grid, the code blocks, the source descriptors, the mutation
transaction, the R-protection story - follows from that.

## Cut grid (store level, singular)

`ColumnStore` IS the cut grid. The per-column quantization parameters live
directly on the store as parallel vectors, one entry per predictor:

- `types` - `ColumnKind::numeric`, `categorical`, or `orderedFactor`
  (`data.hpp:715`). Splitting keys on the DERIVED predicate
  `splitsBySubset(j)`, not on the kind: only grid construction, ingestion
  validation and reporting read the kind itself.
- `numCuts[j]` - the threshold count: `n.cuts` (or the quantile-induced
  count) for a numeric column, K - 1 for an ordered factor, and 0 for a
  categorical column (which has no cut grid at all).
  `refreshCutsForColumn` (the mutation re-cut) keeps a numeric column's
  count fixed and leaves a factor column of either kind alone, and only
  `setCutPointsForColumn` (the setCutPoints surface, refused on a factor
  column) changes it.
- `categoryCounts[j]` - the fixed level count K of a factor column of
  either kind, 0 for a numeric one. Fixed at build: every mask tier,
  reserved missing code and category histogram width derives from it.
- `cutPoints[j]` - the ascending thresholds (empty for categoricals). An
  ordered factor's are the K - 1 midpoints between consecutive declared
  level codes, so its code is its own level index and every adjacent
  level pair is separable.
- `maxNumCuts[j]` - the cap on quantile-induced counts, itself capped at
  `maxNumCutsRepresentable` so the reserved missing code `naCode` (0xFFFF)
  never collides with a real code. An ordered factor RAISES it to K - 1
  where the level table asks for more, keeping
  `numCuts[j] <= maxNumCuts[j]` while the grid follows the levels.
- `useQuantiles` - grid mode (store-wide).
- `hasMissing[j]`, `hasPooledCategorical` - derived flags gating the
  NA-aware partition kernel and the pooled-mask machinery.

There is deliberately no nested `CutGrid` type. The grid is physically
singular on the store and field-read at 150+ call sites across scan, tree,
model, and facade; nesting it would buy zero deduplication and force a
rename everywhere. The store is the shared grid, and the train and test
code blocks quantize against it by identity - a view or a test block bins
identically to its parent by construction because it copies these fields,
not because any code re-derives them.

`codeFor(j, value)` (`data.hpp:910`) is the one quantizer: ordinal values
map by `lower_bound` over `cutPoints[j]` (a value above every cut takes
code `numCuts[j]`, always right of any split); categorical values are
their own integer code; missing values take `naCode` (ordinal) or
`missingCategoryCode(categoryCounts[j])` (categorical, position 63 inline
or K pooled).

## CodeBlock (train and test)

A `CodeBlock` (`data.hpp:630`) holds one row set's codes over the store's
grid. `ColumnStore` instantiates two: `train` and `test`. Each owns:

- `codes` - the packed dense codes, a single contiguous `xint_t` vector.
- `codeOffsets[j]` - column j's start in `codes` (`j * numRows` for a
  dense-matrix build, packed among densified columns for a CSC build,
  unused for a rank-stored column).
- `sparseColumns` - the rank-bitmap storage of the sparse columns.
- `sources[j]` - the per-column source descriptor (next section).

The block's pure accessors read whichever storage a column takes:
`column(j)` for dense codes, `codeAt(j, i)` / (store-forwarded)
`testCodeAt(j, i)` for storage-aware single-code access that routes a tree
descent without materializing a row. The grid-dependent operations that
FILL a block - `quantizeDenseInto`, `quantizeCscColumnInto`,
`buildRankStorageInto` - are `ColumnStore` methods parameterized by the
block and its row count, so train and test run one body each.

What stays store-level, not block-level, and why:

- `hasMissing[j]` is train-only. The flag gates prior draws and the
  partition kernel; the test side gates no draws, so `quantizeTestColumn`
  passes a null `hasMissingOut` and tracks nothing.
- `gathered*` (raw copies, means, sds) is store-level: the leaf-covariate
  raw spans both sides, each with its own column list, not a per-block
  code concern.
- `ownedTestValues`, `ownedTestCscValues`, `ownedTestCscRows`,
  `testOffset` are store-level: the test store owns every raw it keeps
  (below), and the offset is a fit concern, not a code concern.
- `isView` (provenance), `builtFromCsc`, `hasSparse` are store-level flags
  with external readers; the bridge's refusal helpers read the capability
  predicates `hasRequantizeSource` and `acceptsNewRawPredictors` derived
  from them, not the flags directly.

## ColumnSource and its five kinds

Per-column storage is one explicit descriptor, `ColumnSource`
(`data.hpp:273`), carried in `CodeBlock::sources` and sized to
`numPredictors` on any side that has rows (train always after a build;
test whenever `numTestObservations > 0`; both empty on a reset test
store). `ColumnSourceKind` (`data.hpp:246`) discriminates five kinds; each
reads only the descriptor fields it owns:

The discriminator names WHICH POOL the column's raw lives in, not who owns
the codes: every kind's codes are the store's. "This column has no double"
is a kind, not a null a reader must know to expect.

- `denseCallSupplied` - dense codes in `codes[]`; the raw arrives with the
  call (train: the call-time `x`; test: nothing - the test store owns every
  raw it keeps). The build-reset baseline: every column is
  `denseCallSupplied` until a builder overwrites it.
- `denseResident` - dense codes in `codes[]`; re-quantizes from the
  `residentRaw` pointer the descriptor holds (a mixed build's store-owned
  dense slice - the REAL-VALUED dense columns are copied into
  `ownedDenseValues` at build, not borrowed - or a test build's owned slice
  into `ownedTestValues`). `residentRaw` is the only kind that serves
  `rawColumn`/`rawTestColumn` a raw pointer directly, and the only one a
  mutation's write-through and rollback touch.
- `denseCodesOnly` - dense codes in `codes[]` and NO raw anywhere. What a
  FACTOR column of either kind takes wherever its storage is dense: its
  cells are level codes, which the codes already carry, and its grid
  follows the level table rather than its values, so it never re-quantizes.
  A designated leaf covariate among them is served by the gather below.
- `cscRank` - rank-bitmap storage in `sparseColumns[rankSlot]`;
  re-quantizes from the retained `slice`. Below the
  `sparseDensityThreshold` (0.2) nonzero fraction.
- `cscDensified` - dense codes in `codes[]`; re-quantizes from the
  retained `slice`. Above the threshold. Runs the existing SIMD kernels
  bitwise-identically to a dense build of the same values.

Of the two remaining descriptor fields, only `refCode` is CSC-categorical
only. `declaredCategoryCount` (the fixed level count K a host declares -
unrecoverable from observed codes alone when a level has zero training
rows) is train-side only - the test side reuses the store's fixed
`categoryCounts[j]` for K - but rides every storage kind, not just CSC: a
dense-backed factor inside a mixed container declares a level table the
same way a CSC-backed one does. `refCode` is read on BOTH sides: its
test-side value comes from the test view's `referenceCodes`
(`PredictorSource`, `data.hpp`), and `quantizeCscColumnInto` reads it
block-parametrically for the implicit rows (`data.hpp:1296`). It holds the
reference level's level-order code - not the sparse storage's structural
zero; it may be any valid code, including 0. A dense factor codes its
reference level by level order too, which the bitwise-vs-dense gate
forces. A CSC-backed categorical column's implicit rows take `refCode` as
their code; an ordinal one's take the quantized zero.

Re-quantization resolves the raw source per kind through three accessors
whose fallback orders are the contract:

- `rawColumnForRequantize(j, x)` (`data.hpp:909`): factor -> null (it never
  re-quantizes); CSC-backed -> null (the slice serves it); `denseResident`
  -> `residentRaw`; else `x + j*n` (or null if `x` is null). This is the
  mutation/setCutPoints path.
- `rawColumn(j)` (`data.hpp:923`, owned training raw for leaf models):
  gathered slot -> `gatheredRawValues`; `denseResident` -> `residentRaw`;
  else null.
- `rawTestColumn(j)` (`data.hpp:938`): if `test.sources` is populated,
  `denseResident` -> `residentRaw` and CSC-backed -> null (sparse storage
  serves no dense test covariate); else the gathered slot ->
  `gatheredRawTestValues`; else null.

A null from the latter two is not a case a reader handles: the factory that
admits a leaf covariate refuses a designation the store cannot serve raw
for, so a correctly wired leaf model never meets one.

## The borrowed view's two value channels

A host hands the store its values through one borrowed view,
`PredictorSource`, whose dense storage has two channels: `denseValues`,
column-major doubles, and `denseCodes`, column-major int32 level codes with
`naDenseCode` (the minimum int32) for missing. `denseChannels` says which
channel each dense-backed column sits in and where within it - `k >= 0` is
double column `k`, `k < 0` is code column `~k` - so the two index
independently and each is packed over the columns it serves. A null
`denseChannels` is the single-block layout every entrance took before the
code channel existed, which is what a plain matrix (and a data object saved
before it) still uses.

The channel is a STORAGE fact, not a typing one: it says how a host holds a
column, while `ColumnKind` says how the store reads it. A factor column
whose host holds it as integers - which is what an R factor is - crosses as
integers, so no cell is widened at the bridge and narrowed again at the
quantize. Both `build` arms and `buildTest` consume either channel, reading
the view a column at a time: a real-valued column is copied (or widened
once, since its grid is over real values) into the block the side keeps, and
a factor column quantizes straight off its codes and keeps none. A ROWLESS
view - what dropping the test data hands `buildTest` - carries no pointer in
either dense channel, and a view holding no CSC nonzeros at all carries none
for those either, so every one of those copies is skipped rather than run at
zero length: a zero-byte copy from a null source is undefined all the same.
Only the MUTATION entrances still need the dense values as one block,
which their kernels index column-major, and they lay their own out.

The flat replay reads either channel too.
`PredictorSourceColumnReader` - the reader `predict`, the saved-tree replay
and the test-side refusals all take a borrowed view through - has three arms:
the double channel, the code channel, and a CSC rank lookup. The descent asks
`codes()` once per NODE and routes a coded column off its int32 cells against
the largest code the cut admits, so a factor test set costs no per-row widen.

Two conventions are NOT the same and must not be confused. The channel's
missing marker is `naDenseCode`; the code a missing value takes IN THE
STORE is `missingCategoryCode(K)` for a subset-splitting column and
`naCode` otherwise. The forward map is total; the inverse needs the
column's K, which is why the store's codes never cross the view's code
channel: the channel's missing marker is not theirs, and a reader has no
K to tell a reserved code from a level.

## Codes vs raw: the gathered mechanism

Quantized codes are what trees consume. Raw doubles exist for exactly two
reasons: re-quantizing an unchanged column (bin counts do not locate new
quantile values), and leaf models that read covariates
(linear/GP leaves). The engine keeps no resident predictor matrix, so a
consumer that needs raw after the build borrow releases must have asked
for an owned copy.

`gatheredRawColumns` names those columns. `build`/`setData` gather a
sampler's designated leaf covariates - or, for a data handle, the
leaf-covariate columns declared at its creation (empty for a
constant-leaf consumer) - into `gatheredRawValues`
(column-major, `numObservations x q`), refreshed in the same pass as each
column quantizes (`quantizeColumn`, `data.hpp:1416`). `rawColumn` then
serves owned memory for the store's lifetime, borrow long since released.
Few columns are gathered, so the slot lookup is a linear scan
(`gatheredSlotForColumn`).

**What a build gathers is what it cannot otherwise serve.** A dense build
retains no raw at all, so it copies every designated column. A MAPPED build
already serves its real-valued dense-backed columns from `ownedDenseValues`
and takes only the FACTOR ones among the designation - never a CSC-backed
column, which serves no dense raw at all and which the quantize would leave
a slot of zeros; the view built from that store then finds `rawColumn` null
there and its designation refused. That rule is what keeps an ordered factor
admissible as a leaf covariate on a mixed store: the block has no slice for
it, so the gather does.

The test-side twin `gatheredRawTestValues` carries its OWN column list,
`gatheredRawTestColumns`, because the two sides gather different subsets: a
top-level test store owns every real-valued dense test column already and
gathers only the designated FACTOR columns its block cannot serve, while a
view gathers each side of the same designation. `buildTest` fills it; a
reset test store drops it with the rest of the test store.

## View semantics (buildFromParent)

A view (`buildFromParent`, `data.hpp:1823`) is a row- and column-subset of
a built parent store, used by xbart folds and the data-handle path. It:

- copies the parent's grid fields (`types`, `cutPoints`, `numCuts`,
  `categoryCounts`, `maxNumCuts`) for the spanned columns, so it bins
  identically by construction;
- DENSIFIES: gathered train and test codes are fully dense whatever the
  parent's per-column storage, so the sparse-specific paths never run in a
  fold (every `sources` entry stays `denseCallSupplied`, `sparseColumns`
  empty);
- gathers raw only for designated leaf-covariate columns the parent can
  serve (`parent.rawColumn(...) != null`), on BOTH sides and into both
  column lists, inheriting the parent's standardization constants
  (`suppliedStandardization`) so a full-rows view standardizes
  bit-identically to a sampler over the raw data. Columns the parent cannot
  serve are left ungathered; the view's `rawColumn` returns null there and
  the facade refuses the designation. A factor column is one the parent can
  serve exactly when the parent gathered it, which is why a data handle
  declares its leaf covariates whatever its container's shape.

The view is self-contained: nothing references the parent after
`buildFromParent` returns.

`isView` is pure provenance - this store was built from a parent - and
gates only the parent-derived standardization constants
(`suppliedStandardization`). Refusal reads two capability predicates on
the store, one per granularity: `acceptsNewRawPredictors` (whole-data
replacement; false for views, which hold no raw, and for CSC-built
stores, which fix the design at creation - the bridge's
`refusePredictorMutation`) and `hasRequantizeSource` (column mutation,
cut installation and state restore; false only for views, which retain
no re-quantize source at all - the bridge's `refuseMutationOnView`).

## Predictor mutation transaction

Mutation is transactional and lives in the sampler, not the store: the
store exposes primitive re-quantize/rollback methods, and
`runPredictorTransaction` (`sampler.hpp:1850`) sequences them. Two
strategies parameterize it:

- `WholeMatrixUpdate` (`sampler.hpp:1730`), driving `setPredictor`: moves
  the whole live `train.codes` aside into `oldCodes`, rebuilds into fresh
  storage, snapshots `hasMissing` and (if cuts refresh) `cutPoints`. A
  reject swaps the codes back by move; an accept drops them. No
  whole-matrix copy survives.
- `SubsetUpdate` (`sampler.hpp:1786`), driving `updatePredictor`:
  journals each touched column cell-by-cell via `setColumnJournaled`,
  recording each changed cell's old code into a `ColumnCodeRollback`
  (`data.hpp:2137`). Past a quarter of the column changed, the journal
  falls back to a whole pre-change column copy and stops journaling. Per
  column it also snapshots `hasMissing[j]` and (if cuts refresh)
  `cutPoints[j]`.

The transaction sequence (`sampler.hpp:1850`):

1. Precheck every column with `cutsWouldRemainValid` - quantile
   feasibility for numeric columns when cuts refresh, level-code validity
   always for a factor column of either kind. A failure returns
   `invalidCutPoints`, nothing mutated.
2. `forceUpdate`: apply, `forceRefreshTrees` (every forest, collapsing
   emptied leaves), requantize the test columns, accept.
3. Otherwise snapshot-apply, then `revalidateAllChains`. The raw a reject
   must put back is snapshotted HERE, in the transaction: the gathered
   leaf-covariate copies (`oldGatheredRaw`) and, on a mixed store, the owned
   dense block of the touched columns (`oldOwnedDense`), both of which the
   store's re-quantize writes in place. On success requantize test and
   accept; on failure restore both, `strategy.restore`, repartition every
   chain, and return `rolledBack`.

Snapshot ownership: the strategy owns the codes, missing flags, cut grids
and CSC storage it moves; the transaction owns the two raw snapshots. This
is the whole rollback record - there is no store-resident undo log.

Refusal contract: `refuseMutationOnView` guards these paths -
`setPredictor`, `updatePredictor`, the per-observation entries, cut
installation, state restore, forest installation - and refuses only a
data-handle view, which retains no re-quantize source. CSC and mixed
stores keep all of them, re-quantizing from their retained slices; a
per-observation entry additionally refuses a CSC-BACKED target column
by name (a dense column of a mixed store stays open).
`refusePredictorMutation` guards whole-data `setData` alone, which
refuses views and CSC/mixed both ("sparse predictors fix the design at
creation"). Multi-forest samplers carry no predictor guard at any entry.

## Ownership and borrow rules

Borrowed (valid only for the call, or anchored by R protection for the
store's lifetime):

- the call-time training `x` (dense builds): borrowed for the build or
  re-quantize call only; the store retains nothing unflagged.
- a mixed build's CSC slices (`slice`): borrowed for the store's lifetime,
  pointing into R-owned memory or a bridge-owned assembly.

Owned by the store (survive the borrow):

- all `codes`, `sparseColumns`, cut grids, `gatheredRaw*`;
- a mixed build's dense block (`ownedDenseValues`, which
  `denseResident.residentRaw` points into): copied at build, not borrowed,
  and holding the REAL-VALUED dense-backed columns alone, packed per
  predictor over the columns it serves - a map naming one dense column
  twice therefore gives it a slot per predictor;
- the whole test store's raw: `ownedTestValues` (the real-valued dense
  columns, packed the same way) and
  `ownedTestCscValues`/`ownedTestCscRows` (a mixed/CSC test build packs
  every CSC-backed test column's nonzeros so each slice points into
  storage that never reallocates). The engine borrows no test matrix.

Owned by nobody, because it does not exist: a FACTOR column's doubles, on
either side. What was retained per factor column - 8 bytes a cell in each
owned block - is gone, and a designated leaf covariate among them costs the
gather instead.

R-protection lifetime anchoring (`R_interface_bartcore.cpp`): the external
pointer carries a fixed set of protection slots (`PROT_DATA`,
`PROT_RESPONSE`, ...; `line 41`). Predictors are NOT among them: the
engine owns its codes and the R data object (`PROT_DATA` at creation, plus
the live `sampler$data` the R methods hold) is the sole predictor GC
anchor - there is no `PROT_PREDICTORS`. The mixed container needs no
lifetime special-casing here: the store COPIES what it keeps of its
transiently assembled dense values (`ParsedData.denseAssembly`) into
`ownedDenseValues` during `build`, so the assembly need only survive the
call that builds the store, which it already does as an entrance-scoped
local - no holder/handle field extends its lifetime. Every container
assembles TWO transients, one per value channel (above), and a factor
column costs 4 bytes a cell in the code one rather than 8 in the double
one: the dense columnar and mixed creation funnels, the three test-store
funnels, and the READ-ONLY test funnel (`parseTestSource`, shared by
`predict`, `predictPerForest` and the saved-tree replay) alike, since every
consumer reads the view a column at a time. The mutation entrances are the
exception and assemble one block of doubles, which their kernels index
column-major. The CSC slots borrow R container memory instead, valid while
`dataExpr` stays protected.

This borrow-and-anchor arrangement is an UNDOCUMENTED CONTRACT for any
non-R host. The `python-bindings` TODO entry records it: host-side
ingestion validation (CSC structure, sparse metadata assembly, and the
level-table bounds the store's own entrance checks do not cover) lives in
the R bridge and would need reimplementing, and
mixed/CSC stores hold lifetime borrows anchored by R protection that
another host must reproduce. Do not assume the store keeps a mixed/CSC
predictor alive on its own.

## Role contracts

### Leaf model implementer

You MAY read raw covariates through the store's raw accessors and nothing
else. `rawColumn(j)` gives owned training raw for a gathered column;
`rawTestColumn(j)` gives the test twin; `suppliedStandardization(j, ...)`
gives parent-derived constants on a view (compute your own from
`standardizationMomentsForColumn` when it returns false). This is exactly
what `LinearGaussianLeaf` (`model.hpp:220`, `242`, `262`) and
`GPGaussianLeaf` (`model.hpp:620`) do.

Guarantees: for a column your model designated as a covariate, `rawColumn`
and `rawTestColumn` are both non-null on a top-level store (the build
gathers, or the store's own block serves, whichever the column's values
call for) and non-null on a view built WITH that designation; they are null
when the parent could not serve raw, and the facade refuses the designation
upstream in that case, so a null never reaches a correctly-wired leaf. Missing covariate
values enter at the standardized mean (zero); rules on the same column
still route the missingness. You MUST NOT touch codes, `sources`, the cut
grid, or any block internal - the store re-quantizes and rolls back
underneath you, and the transaction re-gathers your covariates
(`regatherTrainingCovariates` on in-place change, `reinitialize` on
whole-data replacement) at the right moments.

### Response family implementer

You touch NONE of the data layer. A response family owns the response
channel - `workingResponse()`, `workingWeights()`, `latents()`,
`workingWeightsVaryPerSweep()` - which the `Chain` scan consumes
(representative sites: `chain.hpp:585`, `588`, `726`, `781`); the
predictor store is a separate object (`data_`) the family never reads. Ordinal (cumulative-probit,
`docs/design/ordinal.md`) and robust Student-t errors
(`docs/design/robust-errors.md`) are the precedent trail: both added a
family entirely through the response channel and its state, with zero data
layer change. What you consume instead of predictors: the working
response/weights the family maintains, the latents it draws, and the
`SamplerStateData` block it serializes.

### Host/bridge implementer

You own container assembly and the ingestion checks the engine cannot
make, but the store now guards its own entrances: `ColumnStore::build`
refuses a factor column whose codes are fractional or out of range, or
whose level count passes the code type's capacity (65535 categorical,
65534 ordered - one lower, since an ordered factor spends a code on the
upper bin of its K - 1 cut grid), and `buildTest` refuses a test code
outside the training level table, leaving the test store untouched. Both
answer `false`, the factories answer null, and the host raises. `setData`
guards its own entrance too, on both layers: the store's checks that every
factor cell of a whole-data replacement is a level code of the table its
pinned count fixes, before it writes a code, and the sampler's checks the
test matrix against that same table before anything moves - so a test
matrix it cannot ingest refuses the whole call rather than leaving new
training values beside a silently dropped test set. The split: the R
bridge validates CSC structure (rows strictly increasing, unique, in
range - "malformed sparse predictor matrix") and categorical placement - a fully-sparse dgCMatrix x
refuses categoricals outright ("sparse predictor matrices must be
entirely ordinal", `R_interface_bartcore.cpp:1139`), while a mixed
container ADMITS CSC-backed categoricals given the container's reference
metadata (`cscReferenceMeta`/`cscCategoryCountMeta`, resolved per
predictor through `resolveCscCategoricalReferences`). The bridge also
assembles the mixed container's source map; the engine's builders assume
all of that holds. Use the composed helpers - `parseCscMatrix`, `parseMixedContainerBlock`
(`mapColumnSources` + `codeDenseColumn`/`codeDenseCodeColumn`),
`resolveCscCategoricalReferences` - which serve both the training and test
container parses from one body.

State blocks are read BY NAME and defaulted when absent
(`R_interface_bartcore.cpp:6927`): a newer release's additive blocks load
into an older reader, and only an encoding below
`minReadableStateFormatVersion` is refused. Preserve that when you add a
block - name it, default it, do not reorder. Predictors are not serialized
in state; `getPointer()` re-creation rebuilds the store from the stored
data object (the grouped/CSC precedent), which is why CSC samplers
save/load without resident raw.

## Where deeper detail lives

- `docs/design/sparse-columns.md` - the rank-bitmap representation study,
  the density threshold, the CSC ingestion and mixed-input landing notes.
- `docs/design/data-ownership.md` - the owned-container program:
  call-time raw supply, the five implementation plans (container,
  ingestion, mutation, views/sharing, sparse categorical), and the
  test-data-parity superseding note. The "why" behind every ownership
  decision above.
- `docs/plans/archive/data-store-consolidation.md` - the four-stage consolidation
  (transaction helper, reset helpers, `ColumnSource` descriptor,
  `CodeBlock` de-twinning) whose landings carry the invariant decisions,
  including why there is no nested `CutGrid` type.
- `docs/design/pooled-masks.md`, `docs/design/mia-missingness.md` -
  categorical mask tiers and the missing-value routing the codes encode.

Per-function truth remains the Doxygen comments in `data.hpp` and
`sampler.hpp`; this note does not restate them.
