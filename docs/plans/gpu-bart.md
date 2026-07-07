# gpu-bart

agent: opus
rng: n/a (survey memo; any mechanism it endorses gets its own plan)
budget: memo only

## Goal

A considered survey of innovative ways to use GPUs for BART - VD wants
the thinking to happen at some point (2026-07-07), with no clock and no
commitment to any mechanism in advance.

## Context

- Two candidate mechanisms already live in the backlog, each with its
  own item: grow-from-root (XBART-style root-down sampling on the
  cut-scan kernel; batched scans are GPU-shaped) and blocked-jacobi-
  trees (the b = m noise-splitting form updates every tree
  concurrently - the one classic-BART formulation with enough
  concurrent work per step to saturate a device). The memo should
  weigh them but is NOT limited to them.
- Directions the memo should at least consider: prediction/test-fit
  offload (embarrassingly parallel, low risk, modest ceiling);
  many-chain parallelism (device-per-chain or chains-as-batch);
  particle/SMC tree samplers; full-posterior reformulations that trade
  per-step exactness for device-scale throughput (recorded as
  posterior-changing if pursued); mixed CPU/GPU splits where the
  device owns the O(n) kernels (sufficient statistics, partition) and
  the host owns tree logic.
- Constraints any proposal inherits: R package distribution (CRAN has
  no CUDA toolchain - a GPU path is an optional backend or a separate
  package), the engine's thread-count-invariant reduction discipline,
  and the exact-equivalence testing culture (a GPU backend needs its
  own gate class).

## Steps

1. Memo (docs/design/gpu-bart.md): the design space above, an honest
   assessment of each direction's ceiling and cost, and a
   recommendation for which one (if any) earns a prototype.
2. VD reviews; a "nothing earns it yet" outcome is a valid resolution
   and closes the item until the landscape changes.

## Verification

- The memo exists, cites the relevant literature, and ends in a
  recommendation VD has acted on (approved a prototype or recorded a
  documented no).
