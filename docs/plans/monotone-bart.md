# monotone-bart

agent: opus
rng: posterior-changing (constrained prior)
budget: design note + ~500 lines

## Goal

Users can declare monotone-increasing/decreasing relationships per
predictor (mBART, Chipman et al. 2022), enforced through the
tree-level leaf-parameter draw hook the design provisioned.

## Context

- Provision: "Leaf draws are at tree granularity (vector of leaf
  parameters given the tree), defaulting to independent per-leaf
  draws; this is the hook for monotone/shape-constrained BART"
  (docs/design/core-generalization.md:128-131).
- mBART's mechanics: the constraint couples leaves that are neighbors
  along a constrained variable's axis; the tree-level draw becomes a
  truncated multivariate draw (sequential conditional truncated
  normals in practice), and the branch marginal likelihood changes
  (integration over the constrained region) - the design note must
  decide between exact constrained marginals and mBART's approach.
- Surface sketch: monotone = c(x1 = "+", x3 = "-") on dbarts()/bart2,
  resolved into per-variable flags on the model spec.

## Constraints

- Constant leaves only in v1 (linear/gp leaf constraints are a
  different problem); categorical predictors refuse the constraint.
- Exact-posterior gate: a single-tree monotone problem is enumerable
  with truncated-normal leaf marginals via quadrature - build it, per
  the established gate pattern.
- Out of scope: convexity/general shape constraints; multivariate
  monotonicity beyond per-variable.

## Steps

1. Design note (docs/design/monotone.md): constraint representation,
   the tree-level draw algorithm, marginal-likelihood treatment,
   interaction with birth/death (a birth can violate feasibility -
   how mBART handles the proposal), prior calibration under
   truncation.
2. Engine: leaf-draw hook implementation + constrained branch scoring
   per the note; model spec plumbing; R surface argument.
3. Gates: exact-posterior on the toy problem; component tests
   (feasibility invariants under every move); recovery test
   (monotone truth recovered, non-monotone truth shows the expected
   bias).

## Verification

- Exact-posterior gate to MC error; component tests; full tinytest.
- bench-sampler: unconstrained paths unaffected (the hook must not
  slow the default).
