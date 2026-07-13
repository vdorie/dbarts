# within-chain-threading

agent: opus
rng: shifting if the reduction blocking changes accumulation order
     (fixed-block design below keeps thread-count invariance)
budget: design note + ~500 lines if measurements say go

## Goal

Single-chain large-n workloads (n >= 1e5, confirmed common - VD
2026-07-06) get within-chain parallelism over the O(n) passes, without
giving up the engine's thread-count-invariant results.

## Context

- Prior assessment (pre-rewrite TODO): backfitting serializes trees;
  parallelism lives inside each O(n) pass at ~3 barrier syncs x
  numTrees per sweep (~300 at defaults) against 2.1 ms/iteration at
  n = 1e4 - overhead on the order of the win; headroom at n >= 1e5
  (DRAM-bound, 31 ms/iteration). The dormant htm machinery
  (src/misc/hierarchicalThreadManager.c; zero live callers) is the
  nominated substrate, though C++20 std::barrier + std::thread may be
  simpler - substrate choice is step 1's output.
- Interaction: hot-layer-u8 attacks the same DRAM-bound regime by
  halving code bytes; measure threading on top of (or against) the u8
  result so wins are not double-counted.
- Determinism: fixed-size block reduction (constant block size, partial
  results combined in block-index order) makes sums independent of
  which worker computes which block - thread-count invariance holds by
  construction, extending the existing bitwise invariance component
  test. Versus today's serial order it is one-time shifting.
- The embedded single-chain case (run(0, 1) per outer Gibbs sweep) is
  the flagship consumer; the worker pool must persist across run()
  calls to amortize thread startup there.
- Competing mechanism: blocked-jacobi-trees.md - exact parallel tree
  updates via noise-splitting augmentation, m/b barriers per sweep
  instead of ~3m. Step 2's measurements run head-to-head against its
  experiment on the same hardware; one mechanism wins and both plans
  record the outcome.

## Constraints

- Workers never touch the R API (existing rule; the ProgressSink
  pattern is the precedent).
- Component-test bar: results bitwise-identical across n.threads in
  {1, 2, 8} for the parallel paths.
- Go/no-go after step 2's measurements; "still not worth it" is a
  valid outcome recorded here.
- Out of scope: cross-tree parallelism (backfitting is sequential);
  test fits (test-fit-parallel).

## Steps

1. Design note: which passes parallelize (partition, suffstat
   accumulation, fit updates, latent refresh - the last is
   embarrassingly parallel and RNG-consuming, needing per-block
   substreams: scope it explicitly or defer it), substrate choice,
   fixed-block scheme and block size.
2. Prototype the suffstat + fit-update passes only; measure at n in
   {1e4, 1e5, 1e6} x threads {2, 4, 8} on top of current u16 and (if
   landed) u8 hot layers.
3. If go: land behind the existing n.threads control (single chain
   uses the worker budget chains would have); snapshot regeneration
   for the one-time accumulation-order shift; re-record equivalence.

## Verification

- Thread-count bitwise invariance component test extended to the new
  paths.
- Full tinytest after regeneration; equivalence z-mode vs prior
  baseline.
- bench-sampler at n = 1e5: recorded speedup; no regression at n = 1e4
  with threading active but below its size cutoff.

## Landing notes

- Step 1 landed cb284eb: design note docs/design/within-chain-threading.md.
  Substrate DECIDED: a new persistent std::thread + std::barrier pool.
  Neither dormant manager is revived and misc_mt is not extended - both
  are condvar fork-join, measured 5-10x costlier per sync than
  std::barrier (arm64 microbench, relative-only), which ~3m syncs/sweep
  cannot absorb; testFitPool_ stays misc_mt (cold, once-per-run).
  Fixed-block reduction: B = 1024 doubles, scalar block interiors,
  serial path shares the scheme - bitwise across n.threads by
  construction; one-time shifting re-record deferred to step 3. Latent
  refresh DEFERRED (needs per-block RNG substreams, a separate RNG
  architecture change; under 1% of parallelizable work). Predicted
  ceiling ~1.5-1.7x at n = 1e5, ~1.8-2.1x at n = 1e6 - modest; the
  memory wall that killed block fusion caps this too, but the downside
  is bounded at neutral-below-cutoff and threading reaches the
  latency-bound gather that bandwidth amortization could not. Go/no-go
  for step 2: >= 1.4x net at n = 1e5 x 4 threads plus the hard bitwise
  invariance gate; head-to-head vs blocked-jacobi scored on ESS/sec.
  Step 2 awaits direction given the revised (lower) ceiling.
