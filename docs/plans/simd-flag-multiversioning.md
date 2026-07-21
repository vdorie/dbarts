# simd-flag-multiversioning

status: EVALUATION COMPLETE (2026-07-16). Verdict: NO-GO as an end-to-end
sampler lever. Flag-built wide-ISA variants are a real, cheap mechanism -
and are ALREADY the mechanism the shipped elementwise linalg uses - but the
one function family VD named (unrolled Mean/Variance) has zero sampler
callers, the hot draw-path reductions are bitwise-locked to scalar, and the
elementwise linalg is already 256-bit at -O2. No flag variant recovers a
material end-to-end win. bench-sampler was deliberately NOT run (no material
candidate; the box stays free for real gates).

VD's ask: evaluate SIMD compilation flags on high-impact functions that
do not currently have wide-ISA variants, as a complement to the
intrinsics-based multi-versioning (same function under different symbols
per instruction set, dispatched at load by src/misc/simd.c). Named
example: see whether the unrolled means and variances perform better
with an AVX2-compiled version on the x86 bench box.

Priors: docs/plans/simd-survey.md (hotness map; draw-path reductions stay
scalar for bitwise repro; Mean/Variance are NOT draw-path), docs/plans/
x86-simd.md and x86-simd-plan.md (the x86 backlog: three "stubbed"
intrinsics bodies, the setIndexedVectorToConstant gap, the SSR leaf).

## Host / method

dbarts-bench: x86-64 Ubuntu 24.04, 16 core, AVX2+FMA, NO AVX512.
gcc 13.3.0 (Ubuntu). clang is NOT installed on the box - the clang half of
the mechanism-portability question (below) is reasoned, not measured.
R 4.3.3. Package builds at R's default -O2 (baseline ISA = x86-64 == SSE2,
no -mavx); per-file wide-ISA flags reach only src/{misc,external,rc} via
misc/Makefile (SSE2_FLAG/AVX_FLAG/AVX2_FLAG). Note: the worktree tip already
carries the AVX2-detection fix (simd.c:107 __get_cpuid_count(leaf,0,...)),
so AVX2 dispatch is live on this box.

Microbenches: faithful copies of the scalar unrolled moment bodies
(moments.c) and the SSE2 intrinsics bodies (moments_sse2.c), plus hand AVX2
and elementwise linalg, in ~/simd-flag-bench on the box. clock_gettime
MONOTONIC, best (min) and median over repeats with warmup, taskset -c 3,
ns/element. The addv elementwise numbers were re-measured with warmup +
interleave after an uncontrolled first-call turbo-ramp artifact was caught
(min-of-7 without warmup inflated the first kernel of the process ~1.7x);
moment numbers were confirmed clean by the same harness. No repo modified.

## THE CENTRAL RESULT: reductions vs elementwise autovectorize differently

At the package's -O2, gcc's "very-cheap" vectorizer cost model governs.

- ELEMENTWISE (linalg add/sub/axpy): a MANUALLY UNROLLED body (the shipped
  unroll-by-8) autovectorizes to full width from the per-file flag alone,
  no fast-math, bitwise-identical across widths. -O2 -mavx already emits
  256-bit `vaddpd ymm`. A PLAIN (non-unrolled) loop does NOT vectorize at
  -O2 even with -mavx2 (cost model); it needs -O3.
- REDUCTION (mean/variance/suffstat/SSR): the per-file flag alone does
  NOTHING. `sum += x[i]` stays scalar `vaddsd` at -O2 -mavx2, at
  -O2 -mavx2 -mfma, AND at -O2 -mavx2 -mfma -ffast-math. Two things must
  hold together to vectorize a reduction: (1) reassociation permission
  (-ffast-math / -funsafe-math / -fassociative-math), because lane-splitting
  a running sum changes rounding; (2) a cost-model override
  (-fvect-cost-model=dynamic, or -O3), because the -O2 very-cheap model
  rejects the horizontal-reduction epilogue. Confirmed by codegen:

    flags on the unrolled mean/var body      256-bit packed adds emitted?
    -O2 -mavx2                                no  (scalar)
    -O2 -mavx2 -mfma -ffast-math              no  (scalar; cost model)
    -O2 -mavx2 -mfma -ffast-math -fvect-cost-model=dynamic   YES
    -O3 -mavx2 -mfma -ffast-math              YES
    -O3 -mavx2 -mfma          (no fast-math)  no  (needs reassociation)

## Measurement: unrolled Mean / Variance (contiguous), ns/element (min)

"avx2" = -O2 -mavx2 (flag only). "avx2+fm" = -O2 -mavx2 -mfma -ffast-math
(still -O2 cost model). "flag-AVX2" = the working flag build,
-O2 -mavx2 -mfma -ffast-math -fvect-cost-model=dynamic (genuine 256-bit
autovec, per-file mechanism). "sse2" = the shipped hand SSE2 intrinsics.
"hand-AVX2" = 4x256-bit hand intrinsics (ceiling).

MEAN
  n       base   avx2   avx2+fm  flag-AVX2  sse2   hand-AVX2
  1e2     0.144  0.145  0.217    0.073      0.092  0.072
  1e3     0.150  0.140  0.262    0.042      0.072  0.046
  1e4     0.141  0.143  0.277    0.059      0.061  0.060
  1e5     0.145  0.149  0.283    0.080      0.088  0.079
  1e6     0.149  0.149  0.281    0.091      0.106  0.092
  1e7     0.470  0.464  0.463    0.389      0.443  0.392   (DRAM-bound)

VARIANCE (known mean)
  n       base   avx2   avx2+fm  flag-AVX2  sse2   hand-AVX2
  1e2     0.292  0.288  0.329    0.107      0.209  0.114
  1e3     0.246  0.260  0.368    0.065      0.133  0.080
  1e4     0.247  0.257  0.373    0.060      0.136  0.075
  1e5     0.251  0.252  0.371    0.086      0.144  0.089
  1e6     0.261  0.257  0.375    0.092      0.152  0.101
  1e7     0.517  0.510  0.519    0.374      0.428  0.373   (DRAM-bound)

Speedup of flag-AVX2 vs scalar base / vs shipped SSE2:
  mean  small(1e3) 3.6x / 1.7x   med(1e5) 1.8x / 1.1x   large(1e7) 1.2x / 1.1x
  var   small(1e3) 3.8x / 2.0x   med(1e5) 2.9x / 1.7x   large(1e7) 1.4x / 1.1x

Reads: flag-AVX2 == hand-AVX2 (autovec reaches the hand ceiling once the two
gates are set). It beats the shipped SSE2 by ~1.1-2x, biggest cache-resident,
shrinking to ~1.1x once DRAM-bound. "avx2" alone is a no-op; "avx2+fm" at -O2
is a REGRESSION (fast-math perturbs the scalar chain, still no vectors).
Weighted mean (MAC reduction) behaves identically (flag-AVX2 ~3x over scalar
cache-resident, ~1.1x DRAM). Indexed (gather) mean: flag/sse2/hand all tie
scalar (gather-bound, ~0.26 -> 5.5 ns/elem as n outgrows cache); hardware
AVX2 gather (_mm256_i64gather_pd) is 1.15-2.2x SLOWER than scalar loads
everywhere - do not use it. (Matches x86-simd-plan.md.)

## Measurement: elementwise linalg (addv y[i]+=x[i]), ns/element (min)

  n        base(-O2)  -O2 -mavx2  -O3 -mavx2  hand-NT-store
  20000    0.304      0.280       0.118       0.600
  4e6      1.030      1.032       0.956       0.828

  - Plain loop: -O2 -mavx2 == base (NOT vectorized at -O2); -O3 gives 2.6x
    cache-resident, ~1.08x DRAM.
  - The SHIPPED bodies are unroll-by-8, so they DO vectorize at -O2:
    -O2 -mavx -> 256-bit `vaddpd ymm` already; -O2 -mavx2 is bit-identical
    codegen (AVX2 adds nothing for double elementwise - it is an integer/
    gather extension). So the "stubbed" _sse2/_avx linalg arms are NOT
    silently scalar; the manual unroll + per-file -m flag already delivers
    the wide vectors.
  - The ONE thing the flag/autovec path cannot emit is the non-temporal
    (_mm256_stream_pd) store the deleted intrinsics used. NT stores win only
    DRAM-bound and only for store-heavy ops: setVectorToConstant 1.6x at
    n=1e7 (0.47 vs 0.74), addv ~1.24x (0.83 vs 1.03); they LOSE cache-
    resident. That is a hand-intrinsics call, not a flag one.

## Mechanism comparison

(1) Per-file -m flag via misc/Makefile (STATUS QUO for the wide bodies).
    Perf: full width for elementwise (with manual unroll) and for reductions
    (with -ffast-math + cost override). Maintenance: one Make line per TU;
    the kernel stays plain C. Portability: CRAN-clean (already how partition/
    linalg/moments SSE2 ship); one flag set per TU; runtime dispatch keeps
    working; Rtools gcc on Windows honors the same flags. Draw-safety:
    elementwise SAFE (bit-identical across widths); reductions need fast-math
    -> ISA-dependent rounding, so ONLY for waived non-draw kernels.

(2) __attribute__((target("avx2,fma"))) in a single TU at -O2.
    Perf: for a reduction this needs optimize("O3","fast-math") ON THE
    FUNCTION as well - target/fast-math alone stays scalar at the TU's -O2
    (measured: packed=0). With optimize("O3","fast-math") it does vectorize.
    Maintenance: no build change, but two pitfalls bit in testing: (a) the
    `optimize` attribute is a gcc-ism clang IGNORES (clang has no equivalent;
    it needs #pragma clang loop + -ffast-math-scope), so the reduction recipe
    is gcc-only - a portability hole for a CRAN package built by both; (b)
    the CROSS-TARGET INLINING BARRIER: a shared `inline` helper (default
    target) will NOT inline into an avx2-target caller. Measured a
    target-avx2 mean that calls a shared add() helper at 1.16-1.20 ns/elem
    vs 0.49-0.74 for the plain version - 1.6-2.4x SLOWER, from the
    per-element call the barrier forces. gcc also documents `optimize` as not
    for production. Verdict: strictly worse than (1) for this codebase.

(3) Hand intrinsics (STATUS QUO for the vectorized reduction bodies).
    Perf: deterministic full width at plain -O2, no semantics change; only
    mechanism that can emit NT stores and is the only lane-parallel gather
    option (though gather loses here). Maintenance: highest - per-ISA source,
    the moments_sse2.c bulk. Portability: maximal (plain -msse2/-mavx flags,
    every compiler). Draw-safety: an intrinsics reduction is as ISA-dependent
    as an autovec one unless written to a fixed cross-ISA lane layout.

Bottom line: (1) is the right mechanism and is already in use; (2) buys
nothing over (1) here and adds a gcc-only + inlining-barrier tax; (3) stays
justified only where NT stores or a fixed-lane draw kernel are actually
wanted.

## Per-function go/no-go

- Unrolled Mean/Variance/WeightedMean/WeightedVariance (the named target):
  NO-GO. `git grep` finds ZERO callers in src/R/inst outside moments.c /
  stats.h - the sampler's live reduction is the SEPARATE scalar undispatched
  SufficientStatisticsFast family (tree.hpp:505-518), not these. A flag-AVX2
  variant is a measured 1.1-2x over SSE2 on code the sampler never runs =
  ~0% end-to-end. Do not add it. If a standalone/utility consumer ever
  profiles hot, the drop-in is a moments_avx2.c compiled
  `$(AVX2_FLAG) -ffast-math -fvect-cost-model=dynamic` (or -O3), matching
  hand-AVX2 for free - record the recipe, do not build speculatively.

- SufficientStatisticsFast (draw-path leaf reduction) and SSR (sigma draw):
  NO-GO by invariant. Flag-building them requires -ffast-math reassociation,
  which makes MCMC draws ISA-dependent and breaks machine-independent bitwise
  reproducibility. Off-limits to flags. (Any SIMD here is the separate
  fixed-lane-layout + re-record decision in simd-survey.md #3, unchanged.)

- Elementwise linalg (add/sub/axpy/setVectorToConstant/transpose): NO flag
  WIN available. Already 256-bit AVX at -O2 via unroll + -mavx; AVX2 == AVX
  for doubles. The only unrealized win is NT stores at DRAM-bound n
  (~1.2-1.6x on store-bound ops), which is HAND intrinsics, not a flag, and
  is exactly the x86-simd.md "restore-or-delete" item - keep it there, decide
  on its own measurement. Nothing for this doc to add.

- setIndexedVectorToConstant (fit scatter): NO-GO for flags. Scatter has no
  pre-AVX512 SIMD form; memory-bound. Flags cannot help (as x86-simd-plan.md
  R3 found). Its dispatch-table gap is a correctness/completeness item, not
  a perf one.

## Constraint (restated, load-bearing)

Draw-path reductions STAY SCALAR AND UNDISPATCHED. The MCMC draw is bitwise
machine-independent because every dispatched draw-path kernel is elementwise
or a permutation (width-invariant), and the hot reductions
(SufficientStatisticsFast, SSR) are scalar. A flag-built AVX2 reduction only
vectorizes with -ffast-math reassociation, whose lane-order summation differs
per width -> ISA-dependent draws. The waiver in kernel-vocabulary.md
("reduction order may differ between ISA variants ... exact cross-variant
equality is not part of the contract") is CONFIRMED still present in the tree
and applies ONLY to the standalone Mean/Variance moment family, which is off
the draw path. No recommendation in this doc puts a fast-math reduction on
the draw path.

## Recommended landing shape

None is a go, so nothing to land. The single durable artifact is the recipe,
not code: reductions need per-file `-ffast-math -fvect-cost-model=dynamic`
(or -O3) on top of the ISA flag, and only for waived off-draw-path kernels;
elementwise ops need only the ISA flag plus a manual unroll and are already
done. If the x86-simd.md backlog is picked up, the sole flag-relevant note is
that restoring the three linalg arms needs NO intrinsics for the compute
width (the unroll+flag already vectorizes 256-bit) - restore intrinsics ONLY
if the DRAM-bound NT-store win is judged worth it, and delete the
commented-out bodies either way.
