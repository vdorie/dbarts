# grow-from-root

agent: opus
rng: posterior-changing (a different sampler, same model)
budget: memo + kernel prototype first; implementation planned on
        approval

Note (2026-07-08): VD greenlit steps 1-2. Two corrections since this
plan was written: MoveStrategy is design vocabulary with no code
symbol (docs/architecture.md) - the memo must name the actual
insertion seam (the move dispatch in Chain::metropolisJumpForTree and
the CGMTreePrior/DartPrior structs are what exist); and the cut-scan
kernel is the shared primitive of parallel-bart-frontier.md 3.1
(informed proposals), so the prototype's cost table should report the
numbers that falsifier needs too (cost multiplier vs the current
single-cut move).

## Goal

A recorded, measured decision on adding XBART-style root-down tree
sampling (He & Hahn) as a second MoveStrategy, targeting the n >= 1e5
regime where per-iteration backfitting is DRAM-bound. The cut-scan
histogram kernel it needs is already on the kernel vocabulary's
planned list.

## Context

- Provisioned: "Grow-from-root / particle Gibbs, warm starts:
  cut-scan kernel; tree state import/export"
  (docs/design/core-generalization.md, extensions table); MoveStrategy
  was designed as the seam (implementation 1 is the conjugate MH port).
- The cut-scan kernel (one pass over a node's rows accumulating
  per-cut sufficient statistics for all cuts of all variables) is the
  same histogram primitive GBM implementations build - sequential
  memory streams, SIMD-friendly, and the one BART operation that maps
  to GPU-style parallelism if that is ever wanted (see the GPU note in
  the TODO discussion; not in scope).
- Interaction: grow-from-root changes sampler behavior per iteration
  (it is closer to a stochastic tree re-fit than an MH walk); it
  complements, not replaces, backfitting - surface as an opt-in
  MoveStrategy or a warm-start phase (XBART's own usage: XBART init,
  BART refinement; connects to warm-starts).

## Constraints

- Memo must quantify before any engine work: cut-scan kernel
  prototype in benchmarks/kernels, cost per node-expansion vs the
  current partition + per-cut scoring at n in {1e4, 1e5, 1e6}.
- Statistical validity story required in the memo: exact XBART is not
  a posterior sampler in the MH sense; decide whether dbarts ships it
  as (a) warm-start only, (b) a documented approximate sampler, or
  (c) a Gibbs-valid variant - the exact-posterior gates arbitrate (b)
  and (c).
- Out of scope until approval: everything else.

## Steps

1. Kernel prototype + bench table (as above); record here.
2. Memo (docs/design/grow-from-root.md): validity story, surface
   choice, MoveStrategy fit, warm-start interaction; VD go/no-go.
3. Fresh implementation plan on approval.

## Verification

- The bench table and memo exist; TODO updated with the outcome.
