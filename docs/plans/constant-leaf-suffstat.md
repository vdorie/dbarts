# constant-leaf-suffstat

agent: opus
rng: shifting, hot-path
budget: ~250 lines

## Goal

ConstantGaussianLeaf consumes the design's own sufficient statistic
(sum w, sum wz, sum wz^2) accumulated in one order-insensitive pass,
instead of the classic (mean, effective count, variance) triple that
requires a two-pass, index-order-dependent reduction. This is the
prerequisite that lets state-continuation drop the per-tree index
buffers.

## Context

- Current consumer: ConstantGaussianLeaf::logIntegratedLikelihood
  (src/bartcore/model.hpp:102-145) fed by computeNodeStatistics /
  computeVariance (src/bartcore/tree.hpp:471-513) via
  misc_htm_computeMean / computeVarianceForKnownMean.
- Design target: docs/design/core-generalization.md:120 "(sum w,
  sum wz, sum wz^2), delegating to existing moment kernels".
  LinearGaussianLeaf already uses the crossproduct moment form
  (model.hpp:348-377).
- The two forms are algebraically equivalent for the marginal
  likelihood and posterior draw; only rounding differs.

## Constraints

- Suffstat accumulation is THE hot loop: bench-sampler compare is
  mandatory, zero-regression bar.
- Kernel vocabulary growth follows kernel-vocabulary.md's contract
  (scalar reference first; SIMD only if profiling justifies). A fused
  (sum w, sum wz, sum wz^2) kernel may reuse the existing moment
  kernels' per-ISA structure.
- Out of scope: changing the marginal-likelihood math itself; state
  format changes (state-continuation owns those).

## Steps

1. Derive logIntegratedLikelihood and drawFromPosterior in the
   (Sw, Swz, Swz2) parameterization; exact algebra check against the
   current form in a component test (same inputs, equal to tolerance,
   plus a hand-computed reference).
2. Add the fused accumulation kernel (scalar reference; weighted and
   unweighted, contiguous and indexed variants) to misc.a per the
   vocabulary contract; wire dispatch table entries.
3. Replace computeNodeStatistics/computeVariance's mean/variance passes
   with the single pass; delete the now-unused node.average plumbing if
   nothing else reads it (printTrees reads occupancy - check).
4. Regenerate RNG-locked snapshots (replay whole files); re-record
   equivalence baseline; statistical mode against the previous one.

## Verification

- tests/cpp: new exact-algebra component test; full suite.
- Full tinytest after snapshot regeneration.
- Equivalence z-mode vs prior baseline passes.
- bench-sampler compare on a quiet machine: no metric > 5% slower.
