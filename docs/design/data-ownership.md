# Owned predictor storage

Status: FROZEN 2026-07-11, approved for pre-release implementation
(VD). Direction and panel resolutions accepted 2026-07-06
(never-retain default, mutable-column CoW resolution, shared-handle
views); the 2026-07-11 convergence (recorded in "Resolved
considerations" below) sharpened never-retain into call-time raw
supply, which dissolved the re-cuttable flag, the epoch machinery,
and every lifetime pin. Supersedes the rejected copy-raw plan;
docs/plans/data-ownership.md points here.

## Motivation

The engine borrows REAL(x) for the sampler's lifetime and services
mutation through PROT_* pinning plus const_cast write-through into R
memory. The replacement: an owned, typed, quantized container
(XGBoost-DMatrix-style; mixed categorical/ordinal, dense/sparse
columns), ingested directly from a data.frame so the loosely packed
R-side double matrix need never exist. At n = 1e6, p = 50 with default
n.cuts = 100: ~50 MB of owned u8 codes versus a 400 MB double borrow.

## Panel findings (three independent reviews, verified in code)

Unanimous default: the engine owns quantized data only; no implicit
raw retention. The verified residual consumers of stored raw values:

- Re-cutting an UNCHANGED column: setCutPoints re-quantization
  (data.hpp:450) and quantile refresh (data.hpp:411). Not recoverable
  from codes (bin counts do not locate new quantile values).
- Linear/gp leaves: gather owned standardized copies at (re)initialize
  only (model.hpp:210,231,498,523); no per-draw raw reads. Leaves can
  own a raw gather too (q <= 8 columns), removing store dependence.
- getTrees saved-tree replay reads store.x (R_interface_bartcore.cpp:
  3069); routable from codes while the cut grid is unchanged since
  save; needs a route-or-refuse decision otherwise.
- dbarts.h exposes no raw-x getter; the C ABI is unaffected.

Proposed resolutions: re-cuttability is an explicit creation-time
column flag (default off) so refusals reflect declared state - the
sparse-column precedent (facade.hpp:363) - not policy surprise.
Split panel call (2/3): keep a READ-ONLY borrow of REAL(x) as the
pure-continuous matrix fast path; write-through dies regardless.

## Open considerations (VD, 2026-07-06)

1. updateX and copy-on-write. Engine-side ownership originally gave
   live updates a home that dodges R's copy-on-write; a naive "write
   updated columns back into data@x R-side" (the panel's stale-data
   fix) would reintroduce CoW - a full-matrix duplication per Gibbs
   iteration. Proposed resolution: columns declared MUTABLE at
   creation get engine-owned raw storage; the updateX family lands
   there in place, CoW-free, O(changed cells), with codes updated in
   the same pass. data@x is then a creation-time snapshot BY
   DEFINITION (no write-back, no staleness ambiguity); an extraction
   verb returns owned current raw for mutable columns and the
   snapshot for immutable ones - correct in both cases. Undeclared
   columns refuse mutation at the call with the declared-state error.
2. Sharing across multiple BART models (BCF's prognostic/treatment
   forests; sum-of-BART families with per-model column subsets).
   Proposed resolution: the standalone data handle already planned in
   core-generalization.md:168-172 plus COLUMN-SUBSET views. Kernels
   consume one column at a time (codes.data() + j*n), so a view is a
   column-index list over the shared store - no contiguous block per
   model is needed. Cut tables and codes are shared when grids match
   (the default); a per-model grid override allocates only the
   diverging column. Mutation under sharing follows the single-writer
   rule; a shared mutable column updated once is seen by every model,
   which collapses bairrtt's two-copy setPredictorJointly workaround
   into one update. Aggressive quantization also caps the worst case:
   a fallback private copy of codes is ~8x smaller than the double
   matrix it replaces.

## Column kinds (VD, 2026-07-06)

The container crosses kind {ordinal, categorical} with storage
{dense, CSC} as orthogonal per-column properties. Today 3 of 4 cells
ship (data.hpp:83 ColumnType; buildMixed per-column dispatch,
sparse-columns.md mixed landing notes): dense ordinal (numerics and
ordered factors), dense categorical (unordered factors, membership
splits, <= 65535 levels), sparse ordinal (CSC, rank bitmap or
densified codes). Sparse categorical is the gap; its design is the
ordinal sparse kernel's sibling - CSC over level codes with the
implicit zero as the reference level, membership masks deciding the
implicit rows by whether they contain that level, occupancy of the
reference level counted as n - nnz. Today's store-wide restrictions
(mutation refused while any CSC column exists; test data densified)
are per-source facts the owned container relaxes per column.

## Container design (proposal, 2026-07-07)

One owned BartData replaces ColumnStore's borrow. Per column: a kind
(ordinal or categorical), a storage class (dense codes or CSC), a
code width chosen by cardinality (u8 for <= 255 cuts - the default
n.cuts = 100 fits - u16 above; hot-layer-u8's per-column widths land
here, not as a separate retrofit), the cut table or level table, an
NA policy, and two creation-time flags: mutable (engine keeps an
owned raw column; the updateX family lands there in place,
copy-on-write-free, O(changed cells)) and leaf-covariate (linear/gp
leaves gather their own raw + working buffers at designation, as
today). No flag, no raw: the default column owns codes only.
Undeclared capabilities refuse at the call with the declared-state
error, the sparse precedent.

CALL-TIME RAW SUPPLY (VD, 2026-07-11) replaces both the re-cuttable
flag and the epoch machinery. Every off-hot-path consumer of raw
values takes them as call arguments assembled by the R layer from
the data object it holds: setCutPoints re-quantization, setData
rebuild, and getTrees saved-tree replay (the bridge's existing
newdata replay path, now the only replay path; FlatNode stores
resolved double cut values, so replay is grid-independent and
getTrees behavior is UNCHANGED in every corner, including after a
re-cut and after mutation - the replay matrix is sourced through the
extraction verb, so mutable columns arrive current and immutable
ones as the snapshot, matching today's write-through-kept-current
store.x). No cut-grid versioning exists: one current cut table per
column, the state format's existing snapshot at save, nothing else.

## Ingestion

dbartsData ingests a data.frame directly: numeric -> dense ordinal,
unordered factor -> dense categorical, ordered factor -> ordinal
codes, I()-wrapped sparseVector / dgCMatrix columns -> CSC ordinal,
and sparseFactor(), a small exported wrapper class carrying levels +
a reference level -> CSC categorical (wrapper decided VD 2026-07-06;
exported with the Matrix-convention name sparseFactor, VD
2026-07-11). The R-side double matrix never materializes for
data-frame input; formula input defers factor expansion to ingestion
(factors = "categorical" already avoids dummy columns; the
model.matrix double detour goes away). Matrix input keeps a
READ-ONLY borrow of REAL(x) during CONSTRUCTION ONLY (VD,
2026-07-11): the borrow serves quantization and the gathering of
flagged columns, then releases - no lifetime pin anywhere; the R
data object holding the matrix/frame is what keeps it GC-alive, and
the C++ side retains nothing unflagged. The data object HOLDS the
ingested frame (or matrix): it is the call-time raw source for
re-cut/setData/replay and the GC anchor. Plan 2 INVESTIGATES
DROPPING @x outright (VD: "see if you can drop it"): route every
internal data@x consumer through the extraction verb or codes, drop
the slot if that closes cleanly, keep a snapshot slot only if some
consumer genuinely needs a materialized matrix. The extraction verb
is extract(sampler, "predictors") (no new generic): owned current
values for mutable columns, the snapshot/frame-derived values
otherwise - correct in both cases, and the canonical source for
getTrees replay.

## Sharing

The standalone data handle (core-generalization.md:168-172) owns the
container once; samplers and forests attach through column-subset
views (a column-index list - kernels consume one column at a time, so
views need no contiguity). Cut tables and codes are shared when grids
match; a per-model grid override allocates only the diverging column.
Mutation follows the single-writer rule and one update is visible to
every attached model, collapsing bairrtt's setPredictorJointly
two-copy workaround; BCF's prognostic/treatment forests are the
first multi-model consumer (forest-split-bcf).

## Compatibility

- The PROT_* protection-slot machinery is deleted wholesale
  (construction-only borrow needs no lifetime pin; the R data object
  is the GC anchor); the const_cast writers and rollback write-back
  (sampler.hpp:684,869) are deleted with the write-through.
- dbarts.h: no signature changes; the creation-contract comment
  updates. The C API keeps matrix semantics (borrow fast path).
- State format: the container serializes per-column metadata + codes;
  one version bump when the first implementation lands.
- mutation-journal's build-new-and-swap and cell-granular journal
  operate on the owned codes + flagged raw columns; sparse-extensions'
  per-column relaxations (mutation with CSC present, sparse x.test)
  become flag questions instead of store-wide ones.

## Implementation split (plans to be written when scheduled)

1. container + dense ingestion + borrow fast path (engine core);
2. data.frame-direct ingestion + wrapper class (R surface);
3. flags + mutation surface rewire + data@x snapshot semantics;
4. views/sharing + standalone handle (blocks forest-split-bcf's
   multi-forest data story);
5. sparse categorical kernel + test-side sparse (with
   sparse-extensions).

## Resolved considerations (VD, 2026-07-11 convergence)

- Call-time raw supply replaces the re-cuttable flag: the R layer
  holds the frame/matrix and assembles raw values as arguments for
  re-cut, setData, and getTrees replay. Engine-side raw narrows to
  mutable-flagged columns only.
- Epoch machinery DROPPED: its sole consumer was the route-or-refuse
  decision for replay without resident raw, which call-time supply
  makes moot. getTrees behavior is unchanged in every corner.
- C++ protection is construction-only; the R data object holding the
  frame/matrix is the GC anchor. All PROT_* slots delete.
- data@x: plan 2 attempts to drop the slot (route internal consumers
  through the extraction verb or codes); the data object holds the
  ingested frame.
- Wrapper class: EXPORTED as sparseFactor() (Matrix-package naming
  convention, recognizable to sparse-data users; carries levels +
  reference level).
- Extraction verb: extract(sampler, "predictors"), no new generic.

## Resolved considerations, second round (VD, 2026-07-11 evening)

- dbarts.h freeze LIFTED for this program: stan4bart is the only ABI
  consumer and we own it - do whatever is clean and update in
  lockstep. getTrees gains an explicit training-replay data
  parameter (NULL = the retained creation spec); setState may take
  raw for cross-grid restores. PROT_DATA stays as the creation
  contract and the flat-C GC anchor.
- Mutation raw-source model: the sampler installs updated column
  VECTORS into its stored data object by reference and lets R
  handle GC - O(spine) once storage is column-oriented (plan 2's
  frame), no O(n x p) copy-on-write. A perf-sensitive caller may
  update its own vector in place and re-install the same vector
  (the supported pattern). Plan-1 interim on matrix-held data@x:
  R-side copy-modify, temporary. OPEN AT PLAN 3: reference-install
  removes the CoW rationale for the engine-owned mutable-raw flag;
  plan 3's convergence decides whether the mutable flag survives.
- State format free to change pre-release (VD); plan 1 keeps the
  cutPoints-only encoding anyway on simplicity.
- u8 SIMD on x86: no x86 machine yet (VD will set one up); x86 epi8
  kernels stay dispatch-gated to scalar until the bitwise component
  gate runs on that machine. NEON + scalar validate now.
- Speed reference: bench-sampler-32fc7c8.csv recorded on the granted
  quiet machine before plan-1 implementation; the pre-record compare
  vs 31a4c01 found NO real regression across the season (one noisy
  flag on the sub-ms setPredictor metrics did not reproduce;
  operational note - those metrics show ~6-8% run-to-run noise, so a
  single-run 5% flag on them warrants a re-run before belief).
