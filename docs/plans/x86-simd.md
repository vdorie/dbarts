# x86-simd

agent: sonnet (bench runs are maintainer-assisted on x86 hardware)
rng: neutral (kernel outputs identical to scalar by contract)
budget: ~250 lines

## Goal

The x86 dispatch arms carry real payloads or stop existing: the three
scalar-bodied linalg kernels are restored (if measurement justifies)
or their x86 dispatch arms deleted; the two dispatch-table gaps close;
Windows-arm64 stops silently running scalar kernels.

## Context

- misc_addVectorsInPlace, misc_subtractVectorsInPlace,
  misc_setVectorToConstant: SSE2/AVX bodies are scalar loops with
  intrinsics commented out (src/misc/linearAlgebra.c ~:173-215 and
  siblings) - predates the rewrite (verified on main). The originals
  used non-temporal stores (_mm_stream_pd), which autovectorization
  will not emit; in the DRAM-bound n >= 1e5 regime that difference is
  plausibly real, so "scalar autovectorizes fine" is unproven there.
- misc_setIndexedVectorToConstant never entered the dispatch table;
  the SSR kernel has no SIMD leaf unlike its mean/variance neighbors.
- Windows-arm64: config.h.win never defines COMPILER_SUPPORTS_NEON and
  Makevars.win's ARCH check has no aarch64 branch.
- Speed baselines are ARM-recorded; there is no current x86
  measurement of anything.

## Constraints

- Measure before choosing: benchmarks/kernels on one x86 machine
  (maintainer-run), scalar vs restored-intrinsic bodies, sizes spanning
  cache-resident to DRAM-bound.
- Exact-equality contract: kernel variants must produce identical
  results to the scalar reference on the linalg ops (no reassociation);
  the moment kernels' cross-variant tolerance rule
  (kernel-vocabulary.md) does not apply here.
- Out of scope: new kernels beyond the two named gaps; the u8 widths
  (hot-layer-u8).

## Steps

1. Restore the commented intrinsic bodies on a branch; bench both on
   x86 (add the three ops to benchmarks/kernels if absent); record the
   table here.
2. Keep whichever wins per op/size; delete dead arms outright where
   scalar wins (dispatch falls through to _c) - no commented-out code
   remains either way.
3. Wire misc_setIndexedVectorToConstant into the dispatch table (scalar
   reference at minimum); add an SSR SIMD leaf only if the same bench
   shows headroom (it may not - record).
4. Windows-arm64: define COMPILER_SUPPORTS_NEON in the .win config for
   aarch64 and add the Makevars.win ARCH branch compiling the _neon
   TUs; verify on win-builder or GitHub's windows-11-arm runner.

## Verification

- tests/cpp component tests (kernel-equality tests cover the variants
  the runtime dispatches to).
- Full tinytest on x86 and arm; equivalence exact on arm (recording
  machine).
- The bench table recorded in this file; benchmarks/README notes the
  x86 data point.
