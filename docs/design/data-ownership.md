# Owned predictor storage

Status: COMPLETE. All five implementation plans landed (2026-07-11 through
2026-07-14); docs/plans/data-ownership.md's landing note records "the
data-ownership program (plans 1-5) is COMPLETE." The design below converged
and FROZE 2026-07-11 - that freeze marks the design settling, not
unfinished work; the program itself finished landing three days later.
Supersedes the original copy-raw plan (rejected, VD, 2026-07-06 - see
"Considered and rejected"); docs/plans/data-ownership.md tracks the
plan-by-plan history and points here for the design record.

Summary: the engine now owns predictor data as a typed, quantized,
XGBoost-DMatrix-style container (mixed categorical/ordinal, dense/sparse
columns) ingested directly from a data.frame, replacing the old borrow of
REAL(x) for the sampler's lifetime plus PROT_* pinning and const_cast
write-through. By default a column owns only its quantized codes; the R
layer supplies raw values as call-time arguments (from the data.frame or
matrix it already holds) wherever the engine still needs them - re-cutting,
setData, and getTrees replay - so no engine-side raw retention, versioning,
or lifetime pin is needed for those paths. Two more elaborate mechanisms for
the same problem were designed and then designed OUT once call-time supply
proved sufficient: a re-cuttable creation-time flag plus cut-grid "epoch"
versioning (replaced by call-time supply directly), and an engine-owned raw
copy for "mutable" columns to dodge R's copy-on-write cost (replaced by
reference-install, which turned out to make the CoW hedge unnecessary).
Both are kept below under "Considered and rejected," with the reasoning,
since the record matters even though neither mechanism shipped.

## Motivation

The engine borrows REAL(x) for the sampler's lifetime and services mutation
through PROT_* pinning plus const_cast write-through into R memory. The
replacement: an owned, typed, quantized container (XGBoost-DMatrix-style;
mixed categorical/ordinal, dense/sparse columns), ingested directly from a
data.frame so the loosely packed R-side double matrix need never exist. At
n = 1e6, p = 50 with default n.cuts = 100: ~50 MB of owned u8 codes versus a
400 MB double borrow.

Three independent panel reviews (verified in code) converged on one
unanimous default: the engine owns quantized data only, with no implicit
raw retention. That default only works because every residual consumer of
stored raw values was tracked down and given an explicit, narrow home:

- Re-cutting an UNCHANGED column: setCutPoints re-quantization (data.hpp:928)
  and quantile refresh (data.hpp:897). Not recoverable from codes (bin
  counts do not locate new quantile values).
- Linear/gp leaves: gather owned standardized copies at (re)initialize only
  (model.hpp:998,1013,1362,1381); no per-draw raw reads. Leaves can own a raw
  gather too (q <= 8 columns), removing store dependence.
- getTrees saved-tree replay (bartcore_getTrees,
  R_interface_bartcore.cpp:5956) read store.x; it was routable from codes
  while the cut grid was unchanged since save. As shipped it needs neither:
  the replay reads the training predictors the R method supplies (data@x)
  and the engine keeps no matrix (R_interface_bartcore.cpp:6024-6026).
- dbarts.h exposes no raw-x getter; the C ABI is unaffected.

The panel's proposed resolutions, from that same 2026-07-06 synthesis:
re-cuttability as an explicit creation-time column flag (default off), so
refusals reflect a column's declared state rather than a policy surprise -
the sparse-column precedent (facade.hpp:810-814); and, on a split panel call
(2/3), keeping a READ-ONLY borrow of REAL(x) as the pure-continuous matrix
fast path, with write-through dying regardless of that split. The read-only
borrow survived into the shipped design, narrowed to construction only (see
Ingestion, below); the re-cuttable flag did not - see "Considered and
rejected."

## The adopted design

### Container and columns

One owned BartData replaces ColumnStore's borrow (container design,
2026-07-07). The container crosses kind {ordinal, categorical} with storage
{dense, CSC} as orthogonal per-column properties; per column it also
carries a code width chosen by cardinality (u8 for <= 255 cuts - the
default n.cuts = 100 fits - u16 above; hot-layer-u8's per-column widths
land here, not as a separate retrofit), a cut table or level table, and an
NA policy.

Coverage of the four kind x storage cells is now complete. Three shipped
before this design (data.hpp:92 ColumnType; buildMixed per-column dispatch,
sparse-columns.md mixed landing notes): dense ordinal (numerics and ordered
factors), dense categorical (unordered factors, membership splits, <= 65535
levels), and sparse ordinal (CSC, rank bitmap or densified codes). Sparse
categorical was the gap this design closed (plan 5, landed - see
"Implementation record"): CSC over level codes, with the reference level
implicit and membership masks deciding the implicit rows by whether they
contain that level, occupancy of the reference level counted as n - nnz.
The design sketch called the implicit value "numeric 0"; the shipped kernel
refined that to the reference level's ACTUAL level-order code (never
literally 0), because the bitwise-vs-dense equality gate forced it - a
dense factor codes its reference level by level order too, and the two
representations must match bitwise. The store-wide restrictions that
predated this design (mutation refused while any CSC column exists; test
data densified) are now per-source facts the owned container relaxes per
column instead.

Exactly one creation-time column flag survives into the shipped design:
leaf-covariate, which lets linear/gp leaves gather their own raw values and
working buffers at designation, as they always have. (The container-design
proposal paired this with a second flag, "mutable," for engine-owned raw
storage on updatable columns; that flag was designed, then killed before
being built - see "Considered and rejected." Mutation instead works by
reference-install, below.) With no flag set, a column owns codes only;
undeclared capabilities refuse at the call with a declared-state error -
the same precedent already used for sparse columns (facade.hpp:810-814), so
refusals reflect a column's declared state rather than being a policy
surprise.

### Ingestion

dbartsData ingests a data.frame directly: numeric -> dense ordinal,
unordered factor -> dense categorical, ordered factor -> ordinal codes,
I()-wrapped sparseVector / dgCMatrix columns -> CSC ordinal, and
sparseFactor() - a small wrapper class carrying levels plus a reference
level (proposed VD, 2026-07-06) - -> CSC categorical. It ships EXPORTED,
named sparseFactor per the Matrix package's naming convention so it reads
familiarly to sparse-data users (VD, 2026-07-11). The R-side double matrix
never materializes for data.frame input; formula input defers factor
expansion to ingestion (factors = "categorical" already avoids dummy
columns, so the model.matrix double detour goes away).

Matrix input keeps a READ-ONLY borrow of REAL(x), but only during
CONSTRUCTION (VD, 2026-07-11): the borrow serves quantization and the
gathering of flagged columns, then releases - no lifetime pin anywhere. The
R data object holding the matrix or frame is what keeps it GC-alive; the
C++ side retains nothing that was not explicitly flagged. That data object
also HOLDS the ingested frame (or matrix), which makes it both the GC
anchor and the call-time raw source the rest of this design draws on.

@x was a candidate for outright removal: plan 2 was charged with "see if
you can drop it" (VD), routing every internal @x consumer through the
extraction verb or the codes and dropping the slot if that closed cleanly.
It did not close cleanly enough to drop - plan 2 landed keeping @x,
narrowed to the ingested source (data-ownership-2-ingestion.md). Plan 3
then settled its final semantics: @x reconciles to the collected raw
source - the per-column values R mapped at creation - with an accepted
mutation swapping the affected column(s) in by reference at the mutating
call; it is never a representation of the engine's quantized state and is
never written during sampling (data-ownership-3-mutation.md). The
extraction verb, extract(sampler, "predictors") (no new generic), stays
the on-demand, container-authoritative materializer for both
reference-installed and untouched columns, and is the canonical source for
getTrees replay.

### Call-time raw supply

Three consumers above need raw values the container does not keep by
default: re-cutting an unchanged column, setData, and getTrees replay.
Call-time raw supply (VD, 2026-07-11) is the single mechanism that answers
all three: every off-hot-path consumer of raw values takes them as call
ARGUMENTS, assembled by the R layer from the data object it already holds
(setCutPoints re-quantization, setData rebuild, and getTrees saved-tree
replay - the bridge's existing newdata replay path, now the only replay
path). No cut-grid versioning exists anywhere in the design: one current
cut table per column, plus the state format's existing snapshot at save,
and nothing else.

FlatNode stores resolved double cut values, so replay is grid-independent
and getTrees behavior is UNCHANGED in every corner, including after a
re-cut and after mutation: the replay matrix is sourced through the
extraction verb, so any reference-installed column arrives current and
untouched columns arrive as the creation-time snapshot - matching the
behavior of today's write-through-kept-current store.x, just without the
write-through.

### Mutation: reference-install

Live updates (the updateX family - dbartsSampler's whole reason for
existing) needed a home that dodges R's copy-on-write. Engine-side
ownership originally gave live updates such a home; a naive "write updated
columns back into data@x R-side" (the panel's stale-data fix) would
reintroduce CoW - a full-matrix duplication per Gibbs iteration. The
shipped mechanism is reference-install: the sampler installs the updated
column VECTOR into its stored data object BY REFERENCE and lets R's
garbage collector take it from there - O(spine) once storage is
column-oriented (plan 2's frame), not O(n x p) copy-on-write. Concretely,
an R-surface helper, installPredictorColumns(x, rows, columns, values),
installs a full column or a per-observation slice by reference and returns
the container; codes update in the same pass. A perf-sensitive caller may
update its own vector in place and re-install the same vector - the
supported pattern, with no aliasing hazard, since R's own copy-on-write
makes any OTHER live reference to that vector copy rather than silently
desync.

Sharing (below) follows a single-writer rule; undeclared columns refuse
mutation at the call with the same declared-state error used elsewhere.

### Sharing across models

BCF's prognostic/treatment forests, and sum-of-BART families with
per-model column subsets, need one container shared by several samplers.
The shipped design: a standalone data handle (core-generalization.md:
181-185) owns the container once; samplers and forests attach through
COLUMN-SUBSET views. Kernels already consume one column at a time
(codes.data() + j*n), so a view is just a column-index list - no
contiguous block per model is required. Cut tables and codes are shared
when grids match (the default); a per-model grid override allocates only
the diverging column. Mutation under sharing follows the single-writer
rule; a shared column updated once is visible to every attached model,
which collapses bairrtt's two-copy setPredictorJointly workaround into a
single update. Aggressive quantization caps the worst case too: a fallback
private copy of codes is ~8x smaller than the double matrix it replaces.

This landed as plan 4 (data-ownership-4-views.md), with two items left
open rather than resolved: shared MUTABLE codes were deferred (sharing
settled to READ-ONLY single-writer for now, Q3), and the standalone handle
stays unserialized and unexported (Q4). Commit-level detail is in
"Implementation record" below. BCF's prognostic/treatment forests are the
first multi-model consumer (forest-split-bcf).

## Compatibility

The PROT_* protection-slot machinery is deleted wholesale: a
construction-only borrow needs no lifetime pin, and the R data object is
the GC anchor instead. The const_cast writers and the rollback write-back
(sampler.hpp) are deleted along with the write-through they
served.

dbarts.h's freeze was LIFTED for this program specifically (VD): stan4bart
is the only ABI consumer and dbarts owns it, so the two update in
lockstep rather than the C API being held rigid. Concretely: getTrees
gains an explicit training-replay data parameter (NULL means the retained
creation spec), and setState may take raw values for cross-grid restores.
PROT_DATA itself stays, as the creation contract and the flat-C GC anchor.

The state format changes too: the container serializes per-column
metadata plus codes, with one version bump when the first implementation
landed. (The format was free to change pre-release regardless; plan 1 kept
the cutPoints-only encoding anyway, on simplicity.)

mutation-journal's build-new-and-swap and cell-granular journal operate on
the owned codes plus whatever columns reference-install has touched (not a
flag-declared subset - the mutable flag never shipped; see "Considered and
rejected"). sparse-extensions' per-column relaxations (mutation with a CSC
column present, sparse x.test) likewise become per-column facts instead of
store-wide ones.

## Considered and rejected

Three mechanisms were designed for this problem and did not ship. The
record is kept because each one closes a real design question, even
though none of them survived contact with a simpler alternative.

- The original copy-raw plan (rejected, VD, 2026-07-06). Before the
  owned-quantized direction was approved, a simpler plan proposed just
  copying the raw R matrix rather than owning quantized codes. Rejected in
  favor of the quantized container (this design), which is the whole
  reason the ~50 MB vs ~400 MB comparison in Motivation matters: a raw
  copy would not have bought the memory win.
- The re-cuttable flag and cut-grid epoch machinery. The first-pass panel
  synthesis (2026-07-06) proposed re-cuttability as an explicit
  creation-time column flag (default off), on the sparse-column
  declared-state precedent, plus - to let getTrees replay without
  resident raw values decide whether it could route through codes (grid
  unchanged since save) or had to refuse (grid changed) - a cut-grid
  "epoch" version machinery. Both were REPLACED, not merely deferred, once
  call-time raw supply (above) converged 2026-07-11: if the R layer always
  supplies raw values as call arguments from the data object it holds,
  there is no "replay without resident raw" case left to route or refuse,
  so the epoch machinery has no consumer, and no flag is needed to declare
  which columns tolerate re-cutting. getTrees behavior is unchanged in
  every corner as a result - the strongest evidence the replacement is
  complete, not a narrower stand-in.
- The engine-owned mutable-raw flag. The open-considerations proposal (VD,
  2026-07-06), formalized in the container design (2026-07-07) as one of
  two creation-time flags alongside leaf-covariate: columns declared
  MUTABLE would get engine-owned raw storage, with the updateX family
  landing there in place, CoW-free, O(changed cells), and codes updated in
  the same pass. Under that proposal, data@x would have been a
  creation-time snapshot BY DEFINITION (no write-back, no staleness
  ambiguity), with the extraction verb returning owned current raw for
  mutable columns and the snapshot for immutable ones. This stayed OPEN
  through plan 3's convergence (2026-07-13), which KILLED it - never
  built. The flag was a hedge against R-side copy-on-write cost, and
  reference-install (above) already cuts that cost to O(spine) with no
  second engine-owned copy; the hedge turned out to be unneeded once
  reference-install existed, so the container itself - not a flagged
  subset of it - is the current raw store. The final model differs from
  the proposal accordingly: @x reconciles to the collected raw source for
  ANY reference-installed column, not just ones that would have been
  flagged mutable, and extract(sampler, "predictors") stays the on-demand,
  container-authoritative materializer, with no engine-side raw to
  reroute to.

## Standing numerical guardrail

(review-6): the constant leaf's centered-SS and the BCF 1e-9 floor are
numerically safe ONLY because responses arrive pre-standardized to O(1);
no ownership path may feed raw-scale residuals to the leaves.

## Implementation record

Five plans, landed sequentially:

1. container + dense ingestion + borrow fast path (engine core).
2. data.frame-direct ingestion + wrapper class (R surface). LANDED through
   20c90d3 + docs (data-ownership-2-ingestion.md): @x kept, narrowed to the
   ingested source; sparse-categorical storage (S-CAT) deferred to plan 5.
3. flags + mutation surface rewire + data@x snapshot semantics. LANDED
   through 0ff4842 + docs (data-ownership-3-mutation.md): the engine-owned
   mutable-raw flag killed (never built); reference-install replaces the
   copy-modify interim; data@x reconciled to the collected raw source,
   never the engine's quantized state; extract stays
   container-authoritative.
4. views/sharing + standalone handle (blocks forest-split-bcf's
   multi-forest data story). LANDED through 49bb1d4 + docs
   (data-ownership-4-views.md): per-forest column-availability mask
   (f7763ba), BCF treatment-forest moderator list (4e1fb5b), handle view
   column axis (49bb1d4); sharing settles to READ-ONLY single-writer
   (shared mutable codes deferred, Q3); the standalone handle stays
   unserialized and unexported (Q4).
5. sparse categorical kernel + test-side sparse (with sparse-extensions).
   LANDED through 5461f41 + docs (data-ownership-5-sparse.md): sparseFactor
   accepted as a predictor via the x/y interface, CSC over level codes with
   the reference level implicit; refines the sketch above (Container and
   columns) from "the implicit zero as the reference level" to the
   reference level's ACTUAL level-order code, never numeric 0 - the
   bitwise-vs-dense gate forced it, since a dense factor codes its
   reference level by level order too. Test-side sparse was satisfied by
   densification over training levels at the time, not resident sparse
   testCodes; the mixed flavor's resident dense block (deferred here
   explicitly at plan 2) was also retired, ownership moving to the bridge
   holder/handle. That densification interim is now SUPERSEDED: LANDED
   through 14bef56..22d7116 (docs/plans/test-data-parity.md), the test
   store gained its own per-column typed fields sharing the training cut
   grid, a sparse/frame x.test stays rank-bitmap or densified per column
   (the training tier rule) resident end to end through creation and
   setTestPredictor, and only predict/getTrees(newdata) materialize a
   dense matrix at the R boundary - superseding plan-5 decision 2 / Q4.

Two operational notes from the close of the program (VD, 2026-07-11
evening):

- u8 SIMD on x86: no x86 machine existed yet at that point (VD was setting
  one up); x86 epi8 kernels stay dispatch-gated to scalar until the
  bitwise component gate runs on that machine. NEON and scalar validate
  now.
- Speed reference: bench-sampler-32fc7c8.csv was recorded on the granted
  quiet machine before plan-1 implementation began; the pre-record compare
  against 31a4c01 found NO real regression across the season (one noisy
  flag on the sub-ms setPredictor metrics did not reproduce - those
  metrics show ~6-8% run-to-run noise, so a single-run 5% flag on them
  warrants a re-run before belief).
