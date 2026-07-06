# rbart-loop-profile

agent: sonnet
rng: neutral (measurement; any fix is a follow-up plan)
budget: measurement + report; no product code

## Goal

A measured answer to whether the rbart_vi custom-prior R loop's
per-sample bridge round trip matters, so the item either closes or
gets a targeted fix plan.

## Context

- The loop pays setOffset + run(0L, 1L) per posterior sample with ~10
  fresh SEXPs per bartcore_run call (R/rbart.R); built-in tau priors
  bypass via the in-core decorator, so only custom-prior users pay.
- capi-callbacks' conditioning hook is the natural fix shape if the
  overhead is real - note the dependency, do not build it here.

## Constraints

- Profile realistic sizes: n in {1e3, 1e4, 1e5}, groups in {10, 100},
  against the in-core path as the floor.
- Out of scope: any change to rbart.R beyond instrumentation that is
  removed afterward.

## Steps

1. Rprof + a wall-clock A/B (custom prior expressing the built-in
   model vs the built-in path) at the sizes above.
2. Report: fraction of wall time in bridge overhead vs tree sweeps.
   Threshold for action: > 15% overhead at any measured size.
3. Record the numbers and the go/no-go in this file; if go, write the
   follow-up plan (likely riding capi-callbacks).

## Verification

- The recorded table; no product diff.
