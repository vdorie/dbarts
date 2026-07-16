# simd-flag-multiversioning

status: investigation in flight (queued 2026-07-15 by VD; the evaluation
lands in this file)

VD's ask: evaluate SIMD compilation flags on high-impact functions that
do not currently have wide-ISA variants, as a complement to the
intrinsics-based multi-versioning (same function under different symbols
per instruction set, dispatched at load by src/misc/simd.c). Named
example: see whether the unrolled means and variances perform better
with an AVX2-compiled version on the x86 bench box.

Priors: docs/plans/simd-survey.md (the hotness map and the
draw-path-reductions-stay-scalar bitwise invariant; the dispatched
Mean/Variance kernels are NOT draw-path), docs/plans/x86-simd.md (the
existing x86 backlog: three stubbed intrinsics bodies, the
setIndexedVectorToConstant dispatch gap, the SSR SIMD leaf).
Measurement host: dbarts-bench (x86-64, AVX2+FMA, 16 cores).
