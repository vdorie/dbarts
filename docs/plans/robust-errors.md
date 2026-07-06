# robust-errors

agent: opus
rng: posterior-changing (new residual model)
budget: design note + ~350 lines

## Goal

Continuous responses can carry Student-t residuals: the classic
scale-mixture augmentation (per-observation latent precisions) riding
the workingWeights hook built for logistic, giving outlier-robust fits
- a direct mitigation for range scaling's outlier sensitivity.

## Context

- Mechanism: eps_i ~ t_nu is eps_i | lambda_i ~ N(0, sigma^2 /
  lambda_i), lambda_i ~ Gamma(nu/2, nu/2); the per-iteration lambda
  draw is conjugate given residuals, and lambda enters the tree stage
  exactly as logistic's omega does (ResponseModel::workingWeights,
  docs/design/core-generalization.md phase 3).
- nu: fixed user value first (nu = 4 default is conventional); sampled
  nu on a grid (the DART alpha-grid pattern,
  src/bartcore/model.hpp DartPrior) is the natural follow-up the note
  scopes.
- Composition: user weights multiply the mixture precisions (w_i *
  lambda_i) - the weighted kernels take the product; sigma's prior
  anchoring is unchanged (lambda has mean 1).

## Constraints

- Exact-posterior gate: single-tree problem with t likelihood via
  quadrature over the leaf and lambda (or fat-grid integration) per
  the established pattern.
- Surface: resid.dist or a family value - follow family-on-model's
  outcome; probit/logistic unaffected.
- setResponse/setData semantics for lambda (re-draw, like probit
  latents) defined in the note.
- Out of scope: sampled nu (follow-up), asymmetric/quantile errors
  (its own item if ever wanted).

## Steps

1. Design note (docs/design/robust-errors.md): parameterization, nu
   default and its eventual prior, calibration interaction with k and
   sigquant (t variance is nu/(nu-2), not 1).
2. TResponse (or a decorator over GaussianResponse - the grouped
   decorator is the pattern); lambda draws; surface plumbing + Rd.
3. Gates: exact-posterior, component tests (lambda conjugacy on fixed
   inputs), recovery under contaminated synthetic data, mutation
   smoke.

## Verification

- Exact-posterior gate to MC error; component tests; full tinytest;
  bench on gaussian unchanged.
