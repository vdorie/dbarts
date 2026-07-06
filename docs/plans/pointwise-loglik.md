# pointwise-loglik

agent: opus
rng: neutral (pure output; computed only when requested)
budget: ~250 lines engine + R

## Goal

Fits can hand loo/WAIC machinery what they need: per-observation
log-likelihood draws as an opt-in Results channel (n x samples x
chains), with an R accessor in the shape the loo package documents.

## Context

- Results fields are caller-owned with null-means-skip
  (src/bartcore/chain.hpp; dbarts.h dbarts_results) - the pattern
  extends directly; storeSample is the fill site.
- Per-family pointwise terms: gaussian
  dnorm(y_i | f_i + offset, sigma, log) with weights entering as
  precision; probit pbinom via log Phi((2y-1)(f+offset)); logistic
  log plogis likewise - all computable from quantities live at
  storeSample time. Grouped fits include the ranef in the linear
  predictor (matches how train fits record).
- loo integration: Suggests-level; the accessor returns the
  (samples*chains) x n matrix loo::loo expects plus r_eff helpers via
  the chain dimension.

## Constraints

- Zero cost when not requested (null skip verified in bench).
- Binary families report the y-scale likelihood, not the latent-scale
  density (document; it is what LOO model comparison wants).
- C API: new dbarts_results field rides capi-callbacks' version bump
  if both land; otherwise its own additive bump.
- Out of scope: computing loo inside dbarts; test-set log-lik
  (predict-side follow-up if wanted).

## Steps

1. Engine: Results::logLikelihood + per-family fill in storeSample.
2. Bridge + R: run results gain "loglik" when requested (control flag
   or run argument); packaged fits (bart2 keepevery-consistent) store
   it opt-in; extractor logLik-style method returning loo-shaped
   matrix.
3. Rd + one worked loo example (Suggests, conditional).
4. tinytest: gaussian values equal dnorm recomputation from stored
   sigma/yhat draws; probit equals pnorm recomputation; shapes.

## Verification

- Full tinytest incl. the recomputation equalities; component tests.
- Equivalence exact (neutral); bench with the channel off unchanged.
