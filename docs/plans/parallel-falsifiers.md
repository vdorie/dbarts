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

## Status (DONE 2026-07-08)

All three measured; both flagship candidates SURVIVE (nothing killed).
Driver: benchmarks/R/parallel-falsifiers.R (e1|e2|e3|all). It drives an
INSTRUMENTED engine build - counters compiled into src/bartcore
(falsifier.hpp + chain/tree hooks, read via BC_FALSIFIER_MODE, CSV out to
BC_FALSIFIER_OUT) - that is reverted before commit and never lands, so the
script is inert against the shipped engine and only reproduces numbers on
a checkout that re-adds the instrumentation. Common config: constant-leaf
bart, m = 75, Friedman-1 (p = 10, data seed 11), sampler seed 99, single
thread/chain, n in {1e4, 1e5}; two forests - default (k 2, base 0.95,
power 2) and deep (power 1). Raw CSVs are regenerable (not committed);
`Rscript benchmarks/R/parallel-falsifiers.R all <outdir>` reprints them.

E2 field-fraction profile (kills 3.4 if field < ~70%). Config: ndpost 200,
nskip 200 (bart splits keepTrees=FALSE into two run() calls of 200 sweeps
each; both rows measured). Per-sweep DRAM-byte model: FIELD = residual
roll (32 B/obs/tree) + whole-tree suffstat recompute + fit scatter +
totalFits rebuild + change/swap index-segment snapshot/restore; SCAN =
partition (x + index reads, index-swap writes as an upper band) +
affected-subtree suffstat recompute. Result: field share 82-88% across all
n/config (scan 12-18%), stable, +2pt as n rises 1e4->1e5. Dominant, clear
of the 70% floor, below the ~93% model; removable field ~85% => DRAM drop
~5.5-6.6x (scan-upper band). VERDICT: SURVIVES (does not kill 3.4).

E1 atom census (kills 3.4 if atoms approach n at b in {4,8}). Config:
keeptrees, ndpost 10, nskip 200, keepevery 20 (10 post-burn kept forests),
block starts {0, m-16}, b in {4,8,16}; leaf tuples via
findBottomNodeForObservation on the live post-draw forest, distinct count
per block. Mean atoms (atoms/n): b=4 43-113 (0.0007-0.005); b=8 312-1596
(0.011-0.031); b=16 2196-20755 (0.126-0.229). Atoms 2-4 orders below n at
b in {4,8}; atoms/n FALLS with n (sub-linear growth); saturation only past
b=16. VERDICT: SURVIVES decisively. 3.4 needs BOTH E1 and E2 => SURVIVES.

E3 stale-residual agreement (kills 3.2a if disagree often OR survivors
near 1). Config: ndpost 100, nskip 100. Per move, decision under frozen
start-of-sweep residual (y - totalFits_start + oldFit_t) vs true rolled
residual, same proposal randomness (rng serialized-state save/restore
around a full tree snapshot; the true decision is the stock chain).
Overall agreement 0.966-0.978 (birth 0.92-0.98 low, change 0.98-0.996
high); stage-1 survivor rate 0.068-0.125 (~= true accept rate), below the
modeled 20-40%. High agreement AND far-from-universal survivors => neither
kill trips. VERDICT: SURVIVES; batched-scoring delayed acceptance earns a
prototype.

Instrumentation reverted (git checkout src/bartcore + rm falsifier.hpp)
before commit; pristine tree rebuilt and gated (R CMD INSTALL --preclean;
tests/cpp make && ./test_bartcore). Only benchmarks/ and docs/ land.
