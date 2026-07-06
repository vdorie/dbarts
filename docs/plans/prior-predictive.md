# prior-predictive

agent: sonnet
rng: neutral (new draws from existing prior-sampling entry points)
budget: ~150 lines R + Rd

## Goal

Prior calibration is a documented workflow: an R-level verb draws
forests and leaf parameters from the prior and returns prior
predictive f(x) (and y) draws over the training or supplied grid.

## Context

- The machinery is ported and tested: sampleTreesFromPrior +
  sampleNodeParametersFromPrior exist on the sampler and in dbarts.h
  (docs/design/core-generalization.md, prior sampling section); what
  is missing is a user-facing composition that loops draw -> predict
  and packages the result.
- Use case: choosing k/power/base/node.scale by looking at prior
  f(x) spread on the actual data scale - pairs with prior-constants'
  documentation.

## Constraints

- No engine change; R composition over existing R5 methods only. If a
  loop over draws proves slow at defaults, note it and stop - speeding
  it up is a follow-up, not scope creep.
- Respects family: binary families return latent-scale f and
  probability-scale transforms via the existing
  probabilityFromLatents.
- Out of scope: prior predictive checks/plots (users own those; give
  them draws).

## Steps

1. dbarts_prior_predict(sampler-or-spec, ndraws, newdata = NULL) (name
   per VD's taste at review): builds/uses a sampler, loops prior draws,
   returns draws x n (x chains if applicable) with sigma draws from
   the residual prior for y-scale simulation.
2. Rd with a k-calibration example; one paragraph cross-referenced
   from vignette-refresh.
3. tinytest: moments of prior f(x) match the node prior's analytic
   spread (the existing prior-sampling component test logic, from R).

## Verification

- Full tinytest; R CMD check clean; equivalence untouched.
