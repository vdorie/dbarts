# monotone-prior-draw

agent: opus
rng: shifting on monotone samplers' prior-draw path only; every other
     leaf model bitwise-unchanged (the new branch is if constexpr gated)
budget: ~60-120 engine lines + tests

Found by the 2026-08-04 SBC-extension blind critique while auditing
prior-draw entry points; fix item, no TODO entry (same-session).

## Goal

`samplePriorPredictive` on a monotone fit draws node parameters from the
constrained prior and leaves the sampler in a monotone-feasible state.
Today `Chain::sampleNodeParametersFromPrior` (chain.hpp) installs
unconstrained draws: its scalar branch calls `forest.leaf.drawFromPrior`
per leaf, and `MonotoneConstantGaussianLeaf::drawFromPrior` (model.hpp)
delegates to the plain constant leaf. The drawn `paramByNode` becomes
live fits, so the prior predictive is drawn from the wrong prior AND any
subsequent `run()` starts from a monotone-infeasible state.

## Context

- The constrained machinery exists and is used elsewhere:
  `monotoneNeighborBounds`, `drawTruncatedNormal`, `drawOneLeaf`
  (model.hpp); other chain.hpp sites branch
  `if constexpr (TreeDrawLeafModel<L>)` where the coupled draw matters.
- Exactness subtlety: the monotone leaf prior is iid normals truncated
  to the feasibility cone. A single sequential sweep of the leaf full
  conditionals is NOT an exact joint draw (the conditionals of the
  truncated joint are not the unconstrained-neighbor truncations).
  docs/design/monotone.md section 6 records the constraint machinery.
- `samplePriorPredictive` is new in 1.0-0 (unreleased), so no NEWS
  bullet; fix silently.

## Decision

Exact per-tree rejection: draw the tree's leaf vector iid from the
unconstrained prior, accept iff monotone-feasible (check via the same
neighbor-bounds walk), retry to a cap; on cap exhaustion raise the
engine's usual error rather than install an infeasible or law-bending
state (sorting draws into feasibility changes the law; do not).
Justification: prior trees under the CGM(0.95, 2) depth prior are
shallow, so per-tree acceptance is high; the cap only guards
pathological warm-started structures. Alternative rejected: refusing
`samplePriorPredictive` for monotone fits - drops a shipped verb over a
fixable defect. Evidence that would change this: measured acceptance
low enough at realistic depths that rejection stalls; then revisit with
an exact sequential construction.

## Constraints

- Non-monotone leaf models bitwise-unchanged: guard the new path with
  `if constexpr`, do not touch the existing scalar branch otherwise.
- No public-surface change; no dbarts.h change.
- Out of scope: `sampleTreesFromPrior` (structure draws carry no leaf
  values), the SBC harness, heteroscedastic variance-forest prior draws.

## Steps

1. In `sampleNodeParametersFromPrior`'s scalar branch add the
   `TreeDrawLeafModel<L>` path: per tree, rejection-draw the leaf vector
   from the unconstrained prior until feasible under the tree's monotone
   bounds; cap retries (report the chosen cap); error on exhaustion per
   the engine's existing error convention.
2. tests/cpp: a monotone sampler's post-draw state is feasible (walk
   `monotoneNeighborBounds` over every leaf of every tree); a
   non-monotone sampler's draw is bitwise-identical to a recorded
   pre-change value.
3. tinytest: `samplePriorPredictive` on a monotone `dbarts` fit
   succeeds, subsequent `$run()` proceeds, and the drawn prior means are
   monotone across the constrained column (statistical check, loose).
4. Measure and report per-tree acceptance at the default depth prior
   (temporary instrumentation, removed before landing).

## Verification

- `R CMD INSTALL --preclean .`; `cd tests/cpp && make clean && make &&
  ./test_bartcore`; full `tinytest::test_package("dbarts")`.
- `Rscript benchmarks/R/equivalence.R compare
  benchmarks/baselines/equivalence-7903855.rds` - 27/27 identical (the
  default path must be untouched).
- Local ASAN tests/cpp build (new engine code); CI sanitizer watched to
  green after landing.

## Landing

LANDED 173a710 (2026-08-04). `MonotoneConstantGaussianLeaf::
drawFromPriorForTree` (rejection against `monotoneTreeIsFeasible`,
cap 1e6) + a `TreeDrawLeafModel` branch in
`sampleNodeParametersFromPrior`; law-level test compares the
constrained 1-D/2-D marginals to independent quadrature. Gates re-run
independently at landing: tests/cpp all pass, tinytest 3474/0,
equivalence 27/27 identical draws plus BCF and multinomial trios
bitwise, air clean; implementer ASAN/UBSAN model suite clean. Measured
per-tree acceptance: 63%/try with one of five axes constrained, ~1/L!
over L leaves with all axes constrained (mean 2.5 leaves on prior
trees), P(cap exhausted) ~2e-11 at 8 leaves. Follow-up in the landing
batch: the pre-existing `testMonotoneInteractionCoexistence` monotone
assertion was vacuous (`buildFromFlat` leaves partitions stale);
`repartitionSubtree` now precedes the feasibility walk.
