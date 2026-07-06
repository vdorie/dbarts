# gp-followups

agent: opus
rng: posterior-changing (new sampling machinery)
budget: decision memo first; implementation planned separately if approved

## Goal

The two deferred GP-leaf extensions have a recorded go/no-go each:
sampled lengthscales and low-rank kernels. Neither proceeds without a
consumer or a decision memo.

## Context

- Sampled lengthscales: docs/design/gp-leaves.md stage-4 addendum
  already retracted the "mechanical follow-up" framing - needs a theta
  prior choice, per-sample theta channels in the state and flatten
  formats (which flat-format-v2 must not preclude), and costs multiples
  of the draw phase per sweep (kernel cache invalidation per theta
  draw).
- Low-rank kernels: consumer-gated in the same doc; no consumer named.

## Constraints

- Blocked on: flat-format-v2 (channel shape), a named consumer or VD
  wanting it for its own sake.
- Out of scope until the memo is approved: any code.

## Steps

1. When picked up: write the decision memo (theta prior candidates,
   cache-cost measurement on the stage-3 benchmark, channel design
   sketch) as a gp-leaves.md addendum.
2. If approved, write a fresh implementation plan; this file then
   points to it.

## Verification

- Memo reviewed by VD; TODO updated with the outcome.
