# Kernel vocabulary

The contract between the generic BART core and the compiled kernel library
(`misc.a`). Everything the sampler does per observation must be expressible
as calls into this vocabulary; the generic layer only decides which kernel to
run on which column over which index set, a decision made once per node
operation. This document is normative for the generalized core: adding or
changing a kernel requires updating it and running benchmarks/kernels.

## Dispatch mechanism

- SIMD variants live in separate translation units compiled with per-ISA
  flags (src/misc/Makefile); the dispatched entry points are function
  pointers (e.g. `misc_partitionRange`) installed by `misc_simd_init()` at
  library load (src/misc/simd.c, called from `R_init_dbarts`).
- `misc_simd_setSIMDInstructionSet(i)` reinstalls pointers for a given level
  (not thread safe; benchmarking/testing use only).
  `misc_simd_getMaxSIMDInstructionSet()` reports the host's ceiling.
- `misc_simd_alignment` is the required alignment for `*Aligned*` variants;
  aligned allocations come from misc/memalign.h.

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
- **Unstable**: two-pointer swap; relative order is not preserved. Anything
  that needs order (candidate: sparse columns) requires a new stable variant.
- Codes are assumed valid (no NA sentinel yet; see planned additions).
- ISA variants: C, SSE2, SSE4.1, AVX2, NEON.

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

IMPORTANT: the plain `misc_compute*` moment functions use online (per-element
division) algorithms and are several times slower than the unrolled
accumulators behind the `misc_htm_*` entry points, which accept a null
manager and taskId 0 for serial use. The htm entry points with a null
manager ARE the fast path; bartcore calls them that way. Exposing the fast
serial variants directly (without the manager indirection) is a cleanup
candidate.

### Vector operations (misc/linearAlgebra.h)

Backfitting state maintenance (`totalFits` updates), all over `double`:

- axpy family: `misc_addVectors[WithMultiplier]`, `misc_subtractVectors`,
  in-place and `*Aligned*` in-place variants (dispatched pointers).
- scalar ops: `misc_addScalarToVectorInPlace`, `misc_setVectorToConstant`,
  `misc_setIndexedVectorToConstant`, scalar multiply.
- elementwise: `misc_hadamardMultiplyVectors[InPlace]`.
- reductions: `misc_sumVectorElements`, `misc_sumIndexedVectorElements`.
- layout: `misc_transposeMatrix` (dispatched), `misc_multiplyMatrixIntoVector`.

### Support (not hot-path)

Adaptive radix tree (misc/adaptiveRadixTree.h), binary IO, string utilities,
thread managers themselves. Not part of the per-observation contract.

## Planned additions (generalized core)

Each lands with a scalar reference implementation first (same table slot),
SIMD specializations only when profiling justifies them.

1. **Width variants**: partition and (where profitable) scan kernels for u8
   codes alongside u16; table index becomes (rule kind, width).
2. **Categorical membership partition**:

       size_t misc_partitionRangeCat  (const misc_xint_t* x,
                                       const uint64_t* directions,  // bitset over level codes
                                       size_t* indices, size_t length);
       size_t misc_partitionIndicesCat(...same...);

   Left = bit `code` of `directions` unset, right = set. Scalar first;
   candidate SIMD via table lookup / shuffle for <= 16-level factors.
3. **NA-aware variants**: a reserved per-column NA code plus a
   goes-left/right flag folded into the rule encoding; kernels take the
   encoded rule rather than a bare cut once missingness lands (phase 4).
4. **Sparse partition**: partition an index set by a CSC column (sorted
   nonzero row ids + codes for nonzeros, implicit zero code). Semantics,
   including whether an order-preserving partition is required generally,
   to be settled by prototype before the interface freezes.
5. **Cut-scan histogram** (enables grow-from-root samplers and exhaustive
   change-rule proposals; XGBoost's histogram trick):

       void misc_scanColumn(const misc_xint_t* x, const size_t* indices,
                            size_t length, const double* z, const double* w,
                            double* restrict out);  // out[3 * numCodes]: n, sum, ssq per code

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
- misc.a is not fully R-free: the hierarchical thread manager prints through
  `Rprintf`/`R_FlushConsole`. Standalone consumers (benchmarks/kernels) stub
  these; an injectable output hook is a candidate cleanup for the
  generalized core.
