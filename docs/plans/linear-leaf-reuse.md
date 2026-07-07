# linear-leaf-reuse

agent: opus
rng: neutral if landed (same math, cached)
budget: measurement first; ~100 lines if it proceeds

## Goal

A measured decision on the two linear-leaf micro-optimizations: reusing
the score-phase crossproducts in the draw phase, and reusing chol(V)
within a sweep. Neither lands without numbers showing they matter.

## Context

- Crossproduct reuse needs a persistent per-node stats cache the design
  deliberately avoids ("move rollbacks need no cache invalidation",
  src/bartcore/model.hpp:162-165) - the old "free win" note was wrong.
- chol(V) reuse was measured second-order after the kernel cache
  (docs/design/gp-leaves.md kernel-cache landing, 5.85s -> 4.58s).

## Constraints

- Measure on a linear-leaves workload at n in {1e4, 1e5} with realistic
  column counts; only proceed if the combined win clears ~5% end to
  end.
- Cache correctness burden: any cache must be provably coherent under
  move rollback and the mutation surface, or invalidated wholesale on
  mutation (cheap, coarse, acceptable).
- Out of scope: constant-leaf paths; gp kernel cache (already landed).

## Steps

1. Profile: fraction of linear-leaf sweep time in score-phase
   crossproducts and in chol(V), via benchmarks/ harness additions.
2. Record numbers here; VD go/no-go.
3. If go: cache keyed by (node, columns-epoch), invalidated on any
   mutation entry point; component test that mutation invalidates.

## Measurement (steps 1-2, benchmarks/kernels/linear_leaf.cpp)

Real sweeps, numTrees 50, single chain/thread, leaf basis capped at the
engine's 8 covariates (so "10/50 columns" is infeasible; measured q in
{4, 8}). Denominator is uninstrumented sweep wall time; the leaf-method
time is clocked directly in a shadowing leaf (no cost model), then
attributed by standalone chol and a with/without-U'WU crossproduct pass.
R end-to-end sweep times agree within 3%.

    config                crossproduct   chol(V)    cacheable U'WU -> ceiling
    n=1e4 q=4 linear         91.5%       0.028%       37.0%  ->  33.9%
    n=1e4 q=8 linear         93.0%       0.045%       49.3%  ->  45.9%
    n=1e5 q=4 linear         88.7%       0.003%       39.4%  ->  35.0%
    n=1e5 q=8 linear         91.6%       0.005%       50.1%  ->  45.9%
    n=1e4 q=8 step (deep)    88.2%       0.054%       48.7%  ->  43.0%
    n=1e5 q=8 step (deep)    90.6%       0.008%       49.0%  ->  44.4%

Crossproduct accumulation is 88-93% of every sweep; chol(V) is under
0.06%. BART trees stay shallow (1-3 leaves/tree), so the draw pass
accumulates U'WU + U'Wz + z'Wz over nearly all n each sweep. Only U'WU is
reusable (U'Wz and z'Wz need the observation pass every sweep, since the
residual moves), and the strided u gather it shares is ~half the cost, so
a persistent stats cache's realistic ceiling is 34-46%, not 90%.

Recommendation: chol(V) reuse is NO-GO (< 0.06%, cannot clear 5%).
Crossproduct (U'WU) reuse is GO: 34-46% ceiling, far above 5%, but bears
the plan's cache-coherence burden under move rollback and mutation. VD
decides step 3.

Resolution (VD, 2026-07-07): GO on the U'WU cache; step 3 is
unblocked. chol(V) reuse is closed NO-GO. The step-1 harness
(benchmarks/kernels/linear_leaf.cpp) is the win-verification bench;
bench-sampler is not implicated (its scenarios are constant-leaf).

## Verification

- Component tests; full tinytest; equivalence exact.
- bench on the linear-leaves workload showing the measured win.
