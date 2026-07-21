# grow-from-root-warm-start

agent: opus
rng: neutral - existing paths draw unchanged; the producer consumes RNG
  only when invoked, yielding a different valid STARTING state (the
  warm-starts.md precedent). Default fits stay bit-identical.
budget: ~700-900 lines across ~10-12 files, staged into four
  commit-sized steps. Full review draft with seam anchors:
  the session scratch gfr-warm-start-plan-draft.md (this file is the
  trimmed shipped form).
window: the cut-scan is the shared primitive of
  the informed-proposal kernel and its prototype (parallel-bart-frontier.md
  3.1 / next-action 4); it lands ONCE here (step 1)
  as a self-contained header that frontier work includes unchanged.

## Goal

Build an initial forest by XBART-style root-down stochastic tree
construction (He, Yalov and Hahn 2019) as a WARM START: promote the
cut-scan kernel from benchmarks/kernels/grow_from_root.c into the
engine as an occupancy-aware, leaf-model-templated header; add
growTreeFromRoot recursing on the scan's per-cut integrated
likelihoods; run k grow sweeps in place, then the exact MH sampler
owns the forest. The posterior is UNCHANGED once sampling begins
(memo role (a), docs/design/grow-from-root.md); warm-start only,
never a standalone posterior sampler (memo NO-GO stands).

## Constraints

- Default path FROZEN: no grow request => zero new RNG draws and the
  chain.hpp run sweep untouched. The grow phase is a SEPARATE sweep
  function (duplicated body, not a policy branch in the hot loop);
  that keeps bench-sampler unnecessary - argue from the diff.
- The promoted scan is OCCUPANCY-AWARE: any cut with an empty side
  gets exactly zero selection weight (the MH -1e7 empty-leaf veto in
  logLikelihoodForBranch never runs on this path - TODO's review-6
  note). This contract gets its own component gate.
- The scan omits sumWZ2 (dead weight for the constant leaf: additive
  over any partition, cancels in every within-node comparison).
  Templated on the leaf suffstat so linear/GP adopt it later;
  constant-leaf specialization only in v1. Other leaves fall back to
  prior-grown init.
- Grown trees satisfy the SAME structural invariants MH enforces
  (non-empty leaves, rules in gauge, availability/depth vetoes); a
  grown forest is a legal chain state.
- Split weights use MH's own CGM prior factors (growthProbability,
  P(var), P(cut)) times exp of the leaf's integrated log-likelihood;
  no-split weight (1 - growthProbability) * L(node). Categorical
  handling (prior-mixed vs ordinal-only scan) is the implementer's
  call, documented and gated.

## Steps

1. scan.hpp: leaf-model-templated cut-scan (per-code histogram of
   (count, sumW, sumWZ), prefix-scan to per-cut collapsed left/right
   suffstats, per-cut integrated log-likelihood). Occupancy-aware.
   Component tests: scan == brute-force per-cut recompute (bitwise);
   empty-side cut weight is the never-selected sentinel. RNG-free.
2. grow.hpp growTreeFromRoot (draw {no-split, (var,cut)} from the
   assembled weights on the chain's own rng_, missing-direction coin
   as drawRuleForVariable, birth, recurse) + Chain::growForestFromRoot
   (numSweeps): a separate sweep loop mirroring the run sweep
   (residual roll, grow in place of metropolisJumpForTree,
   sampleParametersAndSetFits, sigma/k/DART updates). Component
   tests: seeded determinism with documented draw count; grown
   forest well-formed and legal; grow-then-MH consistency.
3. Sampler::growFromRoot(numSweeps) across chains on the thread pool
   + facade virtual + bridge entry + R5 method growFromRoot(n.sweeps,
   updateState = FALSE) mirroring sampleTreesFromPrior. tinytest:
   grow/continue/well-formed; cross-sampler donor$growFromRoot ->
   target$installTrees(donor) round trip (no new install code).
4. bart2 surface at the bart.R:453 init fork + Rd + end-to-end
   tests (grow-initialized fit converges and beats prior-init on
   early-iteration train RMSE; seeded reproducibility) + landing
   note in docs/design/grow-from-root.md. DECIDED (VD 2026-07-10,
   approving the plan's recommendations): bart2 gains the count
   argument n.grow.sweeps = 0L (0 = today's prior init; k > 0 runs
   k grow-from-root sweeps then samples as usual) - no strategy
   string, no hidden internal count; the R5 method is
   growFromRoot(n.sweeps = 2L, updateState = FALSE); an explicit
   warm.start together with n.grow.sweeps > 0 is an ERROR
   (ambiguous initialization, the survival family-conflict
   precedent); this item owns the one-time scan landing that
   the informed-kernel prototype (parallel-bart-frontier next-action 4) later
   includes.

## Verification

- tests/cpp component gates above; full tinytest; benchmarks/R/
  equivalence.R vs equivalence-de67cbb.rds must stay IDENTICAL
  (confirmation the default path is untouched - no re-record).
- No exact-posterior gates (posterior unchanged once MH starts);
  bench-sampler not required while the grow sweep stays a separate
  duplicated function off the hot path.

## Status

- 2026-07-10: LANDED as c8c7764 (squash of wt/gfr-warm-start; the four
  staged commits followed the plan's steps exactly). All constraints
  held: the scan (src/bartcore/scan.hpp) is occupancy-aware with its
  own poison-proven component gate, omits sumWZ2 by evaluating the
  marginal with a zero sumWeightedResponseSq (the constant identical
  across every candidate of a node, so the normalized draw is exact),
  and counts members rather than weights; growTreeFromRoot
  (src/bartcore/grow.hpp) uses MH's own CGM prior factors with the
  documented exact draw count (one discrete draw per positive-growth
  node plus one missing coin per split on a missing column);
  categoricals are ordinal-only in v1 (never split by the grower, gated);
  Chain::growForestFromRoot duplicates the run sweep body so the
  default path is untouched by construction. Surface as decided:
  bart2 n.grow.sweeps = 0L (conflict with warm.start errors), R5
  growFromRoot(n.sweeps = 2L, updateState = FALSE) refusing
  linear/gp, bridge backstop. Warm-start claim measured: grow-init
  beats prior-init early-iteration train RMSE at ratios 0.67-0.88
  across seeds (asserted < 0.9x). Landing note in
  docs/design/grow-from-root.md section 7.
- Gates: tests/cpp all pass (scan == brute force bitwise; occupancy
  sentinel poison-proven; seeded determinism; grown-forest legality;
  grow-then-MH consistency); tinytest 2675/0 (2661 + 14); equivalence
  21/21 identical vs equivalence-de67cbb.rds (the frozen-default
  gate, run at steps 3 and 4 and re-run by the reviewer on the landed
  tree); air + lintr clean; pkgdown clean; no dbarts.h ABI change.
- Remaining scope elsewhere: frontier item 4 (informed proposals)
  consumes scan.hpp unchanged; linear/GP scan bins are a later
  adoption behind the same interface.
