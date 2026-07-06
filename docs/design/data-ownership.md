# Owned predictor storage

Status: draft under discussion (2026-07-06); direction approved and
panel resolutions accepted by VD (never-retain default, declared
column flags, read-only borrow fast path, mutable-column CoW
resolution, shared-handle views), not yet approved for implementation.
Supersedes the rejected copy-raw plan; docs/plans/data-ownership.md
points here.

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

## Still open

- data.frame ingestion surface: what dbartsData grows, and whether
  formula processing can defer factor expansion to ingestion. Sparse
  categorical columns declare themselves via a small wrapper class
  carrying levels (decided, VD 2026-07-06); the wrapper's exact
  surface is design detail.
- getTrees saved-replay: route from codes vs refuse after grid change.
- Relationship to hot-layer-u8 (same container decides code widths),
  sparse-extensions, mutation-journal (journal lives on owned
  storage), forest-split-bcf (first multi-model consumer of views).
