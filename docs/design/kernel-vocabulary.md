# Kernel vocabulary

The contract between the generic BART core and the compiled kernel library
(`misc.a`). Everything the sampler does per observation must be expressible
as calls into this vocabulary; the generic layer only decides which kernel to
run on which column over which index set, a decision made once per node
operation. This document is normative for the generalized core: adding or
changing a kernel requires updating it and running benchmarks/kernels. For
where the engine calls into this vocabulary, see docs/architecture.md.

## Dispatch mechanism

- SIMD variants live in separate translation units compiled with per-ISA
  flags (src/misc/Makefile); the dispatched entry points are function
  pointers (e.g. `misc_partitionRange`) installed by `misc_simd_init()` at
  library load (src/misc/simd.c, called from `R_init_dbarts`).
- `misc_simd_setSIMDInstructionSet(i)` reinstalls pointers for a given level
  (not thread safe; benchmarking/testing use only).
  `misc_simd_getMaxSIMDInstructionSet()` reports the host's ceiling.
- Every kernel handles unaligned input; there is no aligned-only variant and
  no alignment contract for callers to satisfy.

## Types

- `misc_xint_t` (misc/types.h): predictor code type, configure-selected via
  `--with-xint-size`, default `uint16_t`. The generalized data model moves
  this to a per-column property (u8/u16); kernels gain width-suffixed
  variants selected through the same tables.
- `misc_size_t` (misc/stddef.h): `size_t`; observation indices.

## Current vocabulary

### Partition (misc/partition.h, src/misc/partition_body.c)

    size_t misc_partitionRange  (const misc_xint_t* x, misc_xint_t cut,
                                 size_t* indices, size_t length);
    size_t misc_partitionIndices(const misc_xint_t* x, misc_xint_t cut,
                                 size_t* indices, size_t length);

Contract:

- Left = `code <= cut`, right = `code > cut` (matches `Rule::goesRight`).
- Returns the number of observations on the left; on exit
  `indices[0 .. left)` are left observations, the rest right.
- `Range`: `x` points at a contiguous segment of `length` codes; `indices`
  is output only, filled with a permutation of `0 .. length-1` (root-node
  case, or any contiguously stored node).
- `Indices`: `indices` is an existing arbitrary index set into the full
  column `x`; permuted in place.
- **Unstable**: two-pointer swap; relative order is not preserved. The one
  candidate consumer of a stable variant (sparse columns) was prototyped
  and rejected in favor of a rank-bitmap layout that keeps this contract
  (docs/design/sparse-columns.md).
- Codes are assumed valid (no NA sentinel yet; see planned additions).
- ISA variants: C, SSE2, SSE4.1, AVX2, NEON.

Sparse sibling (landed 2026-07-04, docs/design/sparse-columns.md):

    size_t misc_partitionIndicesSparse(const uint64_t* bits,
                                       const uint32_t* wordRanks,
                                       const misc_xint_t* nzCodes,
                                       misc_xint_t zeroCode, misc_xint_t cut,
                                       size_t* indices, size_t length);

Same contract as misc_partitionIndices over the rank-bitmap column layout:
code(i) is zeroCode when bit i is clear, else nzCodes[rank(i)] with rank(i)
= wordRanks[i / 64] + popcount of the word's lower bits. Scalar only and a
plain function (no dispatch pointer until a SIMD variant justifies one);
missing-value columns use an engine-side MIA sibling instead, matching the
dense split.

### Sufficient statistics / moments (misc/stats.h, src/misc/moments.c)

Family: mean, variance (given or computing mean), sum of squared residuals;
each crossed with {full-vector, indexed} x {unweighted, weighted}:

    misc_compute[Indexed][Weighted]Mean
    misc_compute[Indexed][Weighted]VarianceForKnownMean
    misc_computeVariance / misc_computeIndexedVariance
    misc_compute[Weighted]SumOfSquaredResiduals

Weighted mean returns w'x / w'1 and reports w'1 (the effective count) through
an out-parameter; variance is ssr / (n - 1). These compute exactly the
constant-leaf sufficient statistic (sum w, sum wz, sum wz^2 in normalized
form); the generalized `LeafModel::accumulate` for the constant leaf is a
thin wrapper over them.

Threaded variants: `misc_mt_*` (flat thread manager) and `misc_htm_*`
(hierarchical thread manager, taskId-scoped so per-chain tasks can subdivide;
misc/thread.h). Single-thread, mt, and htm variants must agree numerically
apart from reduction order.

IMPORTANT: the serial fast path is the `misc_compute*Fast` family: the
unrolled accumulators without the thread-manager indirection. bartcore calls
them directly; the htm entry points delegate to them when run without a
manager (or single-threaded). The plain `misc_compute*` moment functions
dispatch to slower online algorithms at different cutoffs and remain only
for standalone callers.

### Vector operations (misc/linearAlgebra.h)

Backfitting state maintenance (`totalFits` updates), all over `double`:

- axpy family: `misc_subtractVectors`, plus the in-place
  `misc_addVectorsInPlace[WithMultiplier]` and
  `misc_subtractVectorsInPlace` (dispatched pointers).
- scalar ops: `misc_addScalarToVectorInPlace`, `misc_setVectorToConstant`,
  `misc_setIndexedVectorToConstant`, scalar multiply.
- elementwise: `misc_hadamardMultiplyVectors[InPlace]`.
- reductions: `misc_sumVectorElements`, `misc_sumIndexedVectorElements`.
- layout: `misc_multiplyMatrixIntoVector`.

### Support (not hot-path)

Adaptive radix tree (misc/adaptiveRadixTree.h), binary IO, string utilities,
thread managers themselves. Not part of the per-observation contract.

## Planned additions (generalized core)

Each lands with a scalar reference implementation first (same table slot),
SIMD specializations only when profiling justifies them.

1. **Width variants**: partition and (where profitable) scan kernels for u8
   codes alongside u16; table index becomes (rule kind, width).
2. **Categorical membership partition**: LANDED, but not as a misc.a kernel -
   dense categorical membership landed engine-side as
   partitionIndicesByMask (inline, <= 63 levels) and
   partitionIndicesByWideMask (pooled, > 63 levels) in src/bartcore/tree.hpp,
   reading a dense `const xint_t*` column. The sparse-categorical sibling
   (docs/plans/data-ownership-5-sparse.md) mirrors both as
   partitionIndicesSparseByMask / partitionIndicesSparseByWideMask, reading
   through SparseColumnData::at instead. Both pairs live next to the sparse
   MIA partition (partitionIndicesSparseMIA, tree.hpp) and are dispatched
   from partitionChildren's categorical branch (dense vs. columnIsSparse),
   scalar throughout - no table-lookup/shuffle SIMD variant was built.
3. **NA-aware variants**: a reserved per-column NA code plus a
   goes-left/right flag folded into the rule encoding; kernels take the
   encoded rule rather than a bare cut once missingness lands (phase 4).
4. **Sparse partition**: LANDED 2026-07-04 as misc_partitionIndicesSparse
   (see the partition section above); a streaming range variant for
   root-sized segments remains headroom.
5. **Cut-scan histogram**: LANDED 2026-07-10, but not as a misc.a kernel -
   the cut-scan is header-only in bartcore (src/bartcore/scan.hpp
   `scanOrdinalCuts` + src/bartcore/grow.hpp `growTreeFromRoot`), leaf-model
   templated rather than a fixed-suffstat entry point in this vocabulary.
   See docs/design/grow-from-root.md.
6. **Fused suffstat updates** as profiling demands (e.g. weighted count+sum
   in one pass for Polya-Gamma iterations).

## Invariants

- No kernel allocates, throws, or calls back into generic code.
- Index arrays are `size_t`; kernels never reorder anything but their own
  `indices` argument.
- Reduction order may differ between ISA variants and thread counts; exact
  cross-variant equality is not part of the contract (only tolerance-level
  agreement), which is consistent with the no-bit-parity decision.
- Instruction-set switching is a test/bench facility only; production
  installs once at load.
- misc.a is R-free: output goes through the `misc_printf`/`misc_flushOutput`
  hooks (misc/io.h), which default to stderr; the package points them at
  `Rprintf`/`R_FlushConsole` on load. Standalone consumers need no stubs.
