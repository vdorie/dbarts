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
