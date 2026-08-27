# rbart-loop-fix

agent: opus
rng: neutral vs the current loop (same operations, same order; a
     bitwise comparison test is the gate)
budget: ~200 lines

## Goal

rbart_vi's custom-prior path runs one sampler$run() per chain with a
per-sweep R hook instead of n.samples run(0L, 1L) round trips,
reclaiming the measured 15-76% loop overhead (rbart-loop-profile's
table). Draws are unchanged: the hook performs the identical
operations in the identical order.

## Context

- The loop: R/rbart.R:490-529 (one run(0L, 1L) per kept sample at
  :505, ranef/tau draws between, offset reinstalled through the
  mutation surface); the warm-up sweeps at :775 and :809 stay as-is.
- The hook: SweepCallback (src/bartcore/chain.hpp, landed e7fc71a)
  fires on the calling thread before every sweep, true stops. The C
  API's dbarts_sampler_setCallback is for C consumers; rbart needs a
  bridge-internal trampoline that evaluates an R closure.
- Bridge RNG sync: bartcore_run wraps GetRNGstate/PutRNGstate
  (src/R_interface_bartcore.cpp:975-982). An R closure drawing from
  R's stream mid-run must see a consistent .Random.seed, so the
  trampoline entry point must PutRNGstate before each callback and
  GetRNGstate after (or skip the outer pair entirely - the engine
  never touches R's stream during sampling, af04d0c).

## Constraints

- An R error inside the closure must not longjmp through Chain::run's
  C++ frames: evaluate under R_tryEval (or R_UnwindProtect), convert
  failure to a cooperative stop, and re-raise after run returns.
- The callback path must reproduce the loop bitwise on a fixed seed:
  R draws and sampler draws interleave in the same order. A tinytest
  asserts a custom-prior fit equals the loop-driven implementation
  (keep the old loop reachable internally until the test pins it,
  then delete it - no dead code lands).
- The per-sample bookkeeping the loop does R-side (assignInPlace
  slices at :529 and kin) moves into the closure unchanged.
- Out of scope: exposing R-level callbacks publicly; the grouped
  (non-custom) rbart_vi path, which already runs in-core.

## Steps

1. Bridge: C_dbarts_bartcore_runWithCallback(ptr, burnIn, samples,
   closure, rho) registering a trampoline SweepCallback; RNG sync and
   error containment per the constraints.
2. R: rewrite the custom-prior sample loop as a closure over `state`,
   handing it to the new entry point; delete the run(0L, 1L) loop
   once the equality test passes.
3. tinytest: bitwise equality old-vs-new on a fixed seed (both paths
   driven explicitly); the existing test-rbart-custom-prior.R
   regression stays green unchanged (it pins the loop's semantics).
4. Re-run the rbart-loop-profile cells; record the reclaimed overhead
   in this file.

## Verification

- Component tests; full tinytest (test-rbart-custom-prior.R
  untouched and green); equivalence exact (non-rbart paths).
- The step-3 equality test is the rng-neutrality gate; the step-4
  numbers are the payoff record.

## Result (2026-07-07)

Landed: bridge entry dbarts_bartcore_runWithCallback, the closure-driven
rbart_vi_run (loop deleted), and test-rbart-loop-callback.R, which pins
the deleted loop's exact draws (verified bitwise loop-vs-callback across
thin 1/2/3, gaussian/binary, test fits, user callback, and dart before
deletion). Full tinytest 2425/0; equivalence 18/18 identical.

Payoff: none measurable. Same-build A/B (rbart-loop-profile method: 5
interleaved reps, 200 Gibbs iters/run, medians), overhead vs the in-core
floor (A - core)/A for the deleted loop then the callback path:

| n | groups | loop | callback | in-core (s) |
|---|---|---|---|---|
| 1e3 | 10 | 63.7% | 61.9% | 0.037 |
| 1e3 | 100 | 77.4% | 75.6% | 0.038 |
| 1e4 | 10 | 28.4% | 30.4% | 0.336 |
| 1e4 | 100 | 42.7% | 42.7% | 0.343 |
| 1e5 | 10 | 12.3% | 11.0% | 3.537 |
| 1e5 | 100 | 26.7% | 27.3% | 3.486 |

Loop and callback sit within measurement noise (the profile itself saw
up to 1.38x within-arm swing). The overhead over the in-core floor is
the R-side Gibbs work the closure still runs unchanged - the O(n)
group-mean (sapply/mean), the tau slice sampler's per-draw optim, and
the residual/offset arithmetic - not the per-sweep run(0, 1) round trip
this item removed (the profile's own Rprof attribution already put the
round trip at a small slice of the 57% non-engine time). Removing the
per-sweep tryCatch changed nothing, confirming the bridge cost was never
the bottleneck. Reclaiming the rest needs those blocks in C, the
profile's scoped-out follow-up, which this callback entry now hosts.
