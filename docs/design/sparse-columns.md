# Sparse columns

Prototype study settling the representation question for sparse predictor
columns (core-generalization.md data model; kernel-vocabulary.md planned
addition 4; public-surface.md section 7), including whether the engine
needs order-preserving partitions. Prototyped 2026-07-04 in
benchmarks/kernels/sparse.c; nothing lands in the engine from this pass.

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

## Status

Prototype only: the benchmark and this writeup. The engine is untouched;
landing sparse columns waits on the per-column-width data model work it
shares machinery with (u8 hot layers, per-column code widths).
