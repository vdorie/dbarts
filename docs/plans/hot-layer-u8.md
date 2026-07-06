# hot-layer-u8

agent: opus (phase 2); sonnet can run phase 1
rng: neutral (codes are comparison operands; identical comparisons,
     identical draws - no snapshot churn)
budget: phase 1 ~150 lines; phase 2 ~600 lines

## Goal

Columns whose cut count fits eight bits store u8 codes: half the hot-
layer memory and double the SIMD lanes in partition, targeted at the
DRAM-bound n >= 1e5 regime (confirmed common, VD 2026-07-06). The
design doc's own per-column-width hot layer
(docs/design/core-generalization.md:155-156), unblocked.

## Context

- Global width today: src/bartcore/data.hpp:19 (uint16_t, "matches the
  classic engine"), one flat codes vector (data.hpp:139); kernels typed
  on the single misc_xint_t (src/include/misc/partition.h:14-15).
- Default n.cuts = 100 (R/dbarts.R:17); u8 with reserved NA 0xFF caps
  at 254 cuts. Pooled categorical columns carry codes to 65535 - u16
  stays for them, hence per-column.
- kernel-vocabulary.md already plans the (op, width) table; SSE2 has
  unsigned u8 min/max natively (the u8 compare is easier than u16's).
- benchmarks/kernels/bench.c benches the partition kernels.

## Constraints

- Phase gate: phase 2 proceeds only if phase 1 shows a partition win
  or the memory argument alone convinces VD (record either way).
- testCodes stay u16 uniformly (test routing is cold; see
  test-fit-parallel). Sparse nzCodes stay u16 in this pass.
- Zero-regression bench bar for the u16 path (the dispatch must not
  slow the existing width).
- Out of scope: per-column widths for sparse columns
  (sparse-extensions), moment kernels (they consume doubles).

## Steps

Phase 1 (measure):
1. u8 variants of misc_partitionRange/Indices (scalar + one SIMD arm)
   in benchmarks/kernels only; bench u8 vs u16 at n in {1e4, 1e5, 1e6},
   segment sizes down to leaf scale.
2. Report table; VD go/no-go recorded in this file.

Phase 2 (land):
3. Width-suffixed kernels in misc.a per ISA; dispatch table gains the
   width dimension per kernel-vocabulary.md.
4. ColumnStore: per-column width tag + byte-strided storage; accessors
   return typed pointers; tree.hpp partition dispatch switches once per
   node op on (type, width). NA code 0xFF for u8 columns; cut caps
   already leave the reserved code free (data.hpp:24 pattern).
5. Mutation paths (snapshots, setCell, sessions) go width-generic.

## Verification

- Component tests: routing equality u8-column vs the same data forced
  u16 (bitwise identical partitions and draws).
- Full tinytest unchanged (neutrality is the gate); equivalence exact.
- bench-sampler compare at n in {1e4, 1e5}: no regression; record the
  n = 1e5 improvement in the phase-1 table.
