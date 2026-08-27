# sigma-df-zero-weights

agent: opus
rng: posterior-changing (zero-weight fits only; all-positive-weight bit-identical)
budget: small (engine df term + count tracking + two tests)

## Goal

The gaussian `sigma` posterior's degrees of freedom count only positive-weight
rows. Zero-weight observations are documented as ignored - they drop from the
weighted residual sum of squares and every leaf conditional - so they no longer
inflate the posterior df. Fits with no zero weights (including the unweighted
default) draw bit-for-bit as before.

## Derivation

The `sigma^2` Gibbs step is the conjugate update for a scaled-inverse-chi-squared
prior `nu_0, lambda_0`: `sigma^2 | . ~ (nu_0 lambda_0 + S) / chisq(nu_0 + N)`,
where `S = sum_i w_i (y_i - f_i)^2`. The likelihood that produces this update is
`prod_i N(y_i | f_i, sigma^2 / w_i)`, so a row enters the posterior df exactly
when it enters `S`. A `w_i = 0` row is not a zero-precision observation of a
coherent Gaussian (that is an improper, information-free density, not the
"ignored" the docs promise); under the documented semantics it is simply absent.
The correct posterior df is therefore `nu_0 + N_pos`, `N_pos` the count of rows
with `w_i > 0`. Counting all `N` deflates `sigma`: the draw shrinks by roughly
`sqrt((nu_0 + N_pos) / (nu_0 + N))` in the conditional and compounds in
stationarity as trees absorb the residual. Unweighted and all-positive fits have
`N_pos == N`, so nothing changes there. The scale term `nu_0 lambda_0` and `S`
were already correct.

## Fix

- `ChiSquaredScalePrior::drawSigmaSqFromPosterior` takes `numPositiveWeights`
  and uses it (not `numObservations`) for the posterior df; the SSR still scans
  all rows, where zero weights zero their own contribution (src/bartcore/model.hpp).
- `GaussianResponse` tracks `numPositiveWeights_` wherever weights are set:
  construction, `setData`, `setWeights` (helper `countPositiveWeights`, returning
  `n` for null weights). It is the sole caller of the prior draw; the grouped
  decorator and BCF both reach `sigma` through it, so they inherit the count.

## Context

- src/bartcore/model.hpp: ChiSquaredScalePrior (df term); GaussianResponse
  (weight-setting paths + drawSigma).
- src/bartcore/chain.hpp: response_->drawSigma on combinedFits (single + BCF).
- R/A_class.R: the weights validator warning ("weights of 0 will be ignored").
- correctness-audit.md Block 4 ADJUDICATED FINDING (verified z = -264).

## Test

- tests/cpp testSigmaPosteriorDf: n = 20, n_pos = 10, prior df = 3; 4e5 draws
  pin the sample mean and variance of `sigma^2` to the scaled-inverse-chi-squared
  moments at posterior df 13, separating n_pos from n (df 23) decisively.
  Zero-weight rows carry garbage responses to confirm they leave the SSR too.
- inst/tinytest/test-zero-weights.R: paired duplicate design - fit A on n rows
  all w = 1, fit B on the same rows plus exact w = 0 duplicates; asserts the
  sigma-mean ratio B/A in (0.7, 1.4) (unfixed code gives ~0.13).

## Verification / Status

Gate 1 R CMD INSTALL --preclean .: PASS.
Gate 2 tests/cpp (make && ./test_bartcore): PASS ("all tests passed";
"ok: sigma posterior positive-weight df").
Gate 3 tinytest::test_package: PASS (2465 tests, 0 failures; reproducibility
and other positive-weight snapshots unchanged).
Gate 4 equivalence compare equivalence-31a4c01.rds: PASS as designed. 17/18
scenarios identical draws; the `zeroweights` scenario shifts legitimately
(sigma.mean z = -61.74, max |z| = 61.74 - sigma no longer deflated).
Baseline NOT re-recorded.
Gate 5 fixed-build fit-B/fit-A sigma ratio: 1.0000 (was ~0.13 unfixed).
