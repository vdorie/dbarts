# convergence-diagnostics

agent: sonnet
rng: neutral (accessors over existing draws)
budget: ~250 lines R + Rd

## Goal

Multi-chain fits answer "did my chains converge" without the user
hand-rolling array surgery: draws come out in posterior-package
shape, and summary() reports R-hat/ESS for the scalar parameters.

## Context

- Chain-dimensioned draws already exist on fits (sigma, k, varcount,
  varprobs, yhat.*) but as raw arrays with dbarts' native ordering;
  no as_draws support, no diagnostics anywhere.
- The posterior package is the ecosystem standard (rstan, brms, cmdstanr
  interoperate through it); it must be Suggests, not Depends.

## Constraints

- No hard dependency: as_draws methods registered conditionally;
  summary() computes rank-normalized split-R-hat/ESS via posterior when
  available, and degrades to a plain quantile summary with a note when
  not.
- Scalars first: sigma, k, tau (grouped), leaf-count/tree-depth
  summaries if cheaply derivable from stored trees; f(x) draws are
  large - provide the accessor, not automatic diagnostics over n.
- Out of scope: within-fit warnings ("your chains have not converged"
  hard errors); plotting.

## Steps

1. as_draws_array/as_draws_df methods (or a draws() extractor) for
   bart/bart2/rbart fits mapping the native (par x sample x chain)
   arrays to posterior's (iteration, chain, variable) convention -
   note this is the same transpose yhat-shape already pays; do it
   lazily.
2. summary() method gains an R-hat/ESS block for the scalar
   parameters; document the thresholds it flags (rhat > 1.01 noted,
   not enforced).
3. Rd + a vignette-refresh cross-reference (one paragraph there).
4. tinytest: shape round-trips; diagnostics agree with posterior's
   reference values on a fixed fit; graceful degrade without posterior
   (simulate via with_mock or an environment check).

## Verification

- Full tinytest with and without posterior installed.
- R CMD check clean with posterior in Suggests.
