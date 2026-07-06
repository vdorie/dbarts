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

## Verification

- Component tests; full tinytest; equivalence exact.
- bench on the linear-leaves workload showing the measured win.
