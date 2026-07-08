# parallel-falsifiers

agent: opus
rng: neutral (nothing lands in the engine; instrumentation stays in
     the worktree, only drivers/records land)
budget: ~150 lines landed (benchmarks drivers + design-note verdicts);
        instrumentation diffs in-worktree are unbounded and discarded

## Goal

Measured verdicts on the first three falsifiers of
docs/design/parallel-bart-frontier.md section 5, each with its kill
condition applied and recorded in that note: E2 (field-fraction
profile), E1 (atom census), and stale-residual agreement logging.
These decide construction 3.4 (block-fused atom sub-sweeps, the
flagship CPU candidate) and 3.2a (batched stale scoring with delayed
acceptance) before either earns a prototype. "Dead, with the number
that killed it" is a complete outcome.

## Context

- The claims under test are quantified in the frontier note: 3.4
  models a ~93/7 split between residual-field maintenance (leaf-draw
  fit writes, residual rolls, rejected-move bookkeeping - the O(n)
  streams) and cut/move scan work, and needs real-forest atom counts
  well below n at b in {4, 8}; 3.2a needs stale-vs-true acceptance
  agreement high but survivor rate meaningfully below 1.
- The sweep to instrument is Chain::run ->
  Chain::metropolisJumpForTree (src/bartcore/chain.hpp); move scoring
  is logLikelihoodForBranch and friends (src/bartcore/moves.hpp);
  residual maintenance is the totalFits/residual roll in the sweep
  loop.
- Atoms: for a block of b consecutive trees, the distinct tuples of
  per-tree leaf assignments over observations. Leaf routing for saved
  trees exists as flat replay (tests/cpp uses
  countFlatObservationsBelow; src/bartcore/tree.hpp
  findBottomNodeForObservation routes live trees). Real fitted forests
  (defaults and a deep-tree config, post burn-in), NOT simulated ones.
- Driver pattern for standalone measurement code: benchmarks/kernels.

## Constraints

- The shipped engine does not change: instrumentation lives behind the
  worktree only; what lands is (a) any standalone driver worth keeping
  under benchmarks/, (b) the verdicts written into
  parallel-bart-frontier.md section 5, (c) this plan's Status.
- Measurements at n in {1e4, 1e5} (the confirmed-common regime is
  n >= 1e5, single chain), default m = 75 trees plus one deep-tree
  config; fixed seeds; record the configs with the numbers.
- Kill conditions verbatim from the note: E2 kills 3.4 if field
  maintenance is not the dominant share (the model needs ~90%; falling
  under ~70% collapses the 7-10x DRAM claim); E1 kills 3.4 if atoms
  approach n at b in {4, 8} on real forests; the logging kills 3.2a if
  stale-vs-true decisions disagree often OR stage-1 survivors are
  near-universal (either way the 2.5-6x serial-path claim fails).

## Steps

1. E2: count per-sweep element traffic (reads + writes) attributable
   to field maintenance vs scan work, via counters compiled into the
   sweep in the worktree; report the split at both n, both configs.
2. E1: from post-burn-in kept forests, compute distinct-atom counts
   for consecutive blocks at b in {4, 8, 16}, several sweeps apart;
   report atoms vs n and per-atom mean occupancy.
3. Stale-residual logging: per proposed move, score against
   start-of-sweep frozen residuals AND true residuals; report
   agreement rate and the stage-1 survivor fraction by move type.
4. Write the three verdicts into parallel-bart-frontier.md section 5
   (kill or proceed per construction), update the TODO's gpu-bart
   entry, and record here.

## Verification

- Three numbers with configs, seeds, and kill conditions applied;
  frontier note updated; the engine tree diff-clean apart from
  benchmarks/ and docs/.
