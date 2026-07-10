# pointwise-loglik

agent: fable
rng: neutral - pure post-processing of stored draws; no sampling-path or
  RNG-consuming code changes, so equivalence is byte-identical for every
  scenario and even ppd draws are untouched.
budget: ~75 lines in generics.R, ~200 lines of tests, two Rd items.

## Goal

Per-observation, per-draw log-likelihood for fitted models (issues #25/#37):
loo::loo and loo::waic consume an S x n pointwise log-likelihood matrix, and
dbarts had no way to produce one without hand-rolling the density math and
the draw-pairing conventions.

## Scope decision

Implemented R-SIDE as a new extract type, computed from the stored draws.
The engine-side version (a per-family Results::logLikelihood channel filled
at storeSample time, as this file previously planned) is deferred to the
decided c-api-growth item - the C results struct cannot gain fields until
then. The R computation fully serves the demand; test-set log-likelihood
remains a predict-side follow-up if ever wanted.

## Design

extract gains type = "loglik" for bart/bart2/rbart_vi fits, training sample
only (there is no stored test response; sample = "test" errors). It is not
added to predict (no y at test points) or fitted (no per-draw meaning);
those reject "loglik" with their standard vocabulary error. Draw-dimension
conventions follow the other types: S x n combined, chains split when
combineChains = FALSE, shaped by the existing machinery. Binary families
report the y-scale likelihood, not the latent-scale density - what LOO
model comparison wants.

Per draw s, observation i, deriving each family's density from what the fit
stores (stored train fits include any offset; extract "ev" is the full
per-draw location/probability):

- gaussian: dnorm(y_i, ev_si, sigma_s / sqrt(w_i), log = TRUE); weights are
  precision per the documented y | x ~ N(f(x), sigma^2 / w).
- probit: dbinom(y_i, 1, p_si, log = TRUE) with p from extract "ev"
  (includes binaryOffset). Probit fits never store weights (all-1 dropped).
- logistic (bart2 family = "logistic"): integer frequency weights mean w_i
  iid trials with the same recorded outcome, so w_i * dbinom(y_i, 1, p_si,
  log = TRUE).
- rbart_vi: extract "ev" is yhat.train + the drawn group intercepts, so the
  log-likelihood conditions on all drawn location parameters and sigma_s.

## Implementation

R/generics.R:

- pointwiseLogLikelihood(object, ev): shared helper. ev enters with chains
  split ((n.chains x) n.samples x n.obs) - obtained by a recursive
  extract(type = "ev", combineChains = FALSE) - because in that layout
  as.vector(ev) enumerates draws chain-fastest, which is the order
  as.vector(object$sigma) yields in BOTH of sigma's storage layouts (the
  n.chains x n.samples matrix of combineChains = FALSE fits and the
  chain-interleaved combined vector of combineChains = TRUE fits, verified
  empirically), so sigma pairs by plain recycling. The caller reshapes the
  result to the requested chain convention with combineOrUncombineChains.
  Errors informatively if the fit lacks $y.
- extract.bart / extract.rbart: "loglik" added to the type vocabulary, a
  train-only guard after sample validation, and a branch that computes ev
  recursively, applies the helper, and reshapes. All availability errors
  (keepTrainingFits etc.) come from the recursive ev extraction.

Observed while pairing sigma (out of scope, left untouched): a fit stored
with combineChains = TRUE has chain-INTERLEAVED sigma but chain-BLOCKED
yhat rows, so sampleFromPPD's positional gaussian sigma pairing is
misaligned for multi-chain combined fits - statistically invisible since
sigma draws are near-exchangeable, but worth its own item.

man/bart.Rd and man/rbart.Rd: "loglik" in the extract usage vocabulary and
the type item - definition per family, the weighted semantics, and one
sentence noting the combined S x n matrix feeds WAIC/PSIS-LOO packages
(loo) directly.

## Test

inst/tinytest/test-pointwise-loglik.R, bit-agreement (expect_identical)
against hand-computed densities on the stored draws:

1. gaussian unweighted, single chain: entries equal dnorm(y_i, ev_si,
   sigma_s, log = TRUE); dims match "ev".
2. gaussian weighted: a w = 4 column uses sigma / 2, a w = 1 column plain
   sigma.
3. probit: dbinom(y, 1, ev, log = TRUE) columnwise.
4. weighted logistic: a w = 5 column equals 5 * the unweighted form on its
   ev; a w = 1 column the unweighted form.
5. rbart_vi: the location is extract "bart" + extract "ranef" at the
   observation's group, checked entrywise.
6. multi-chain: dims match "ev" for combineChains TRUE and FALSE; each
   chain pairs with its own sigma row (fit stored uncombined, so sigma's
   C x S layout is unambiguous); combined rows are chain blocks; a
   same-seed fit stored combined yields identical loglik (layout
   invariance, covering the interleaved-sigma path).
7. errors: predict and fitted reject type = "loglik" with the standard
   vocabulary error; extract(sample = "test") errors for both classes.

## Verification

Gates, worktree at .claude/worktrees/pointwise-loglik, private library
/Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-land:

1. R CMD INSTALL -l <lib> (no --preclean).
2. tinytest::test_package("dbarts") - baseline 2557 passed / 0 failed plus
   the new checks.
3. Equivalence vs benchmarks/baselines/equivalence-de67cbb.rds - all 21
   scenarios "identical draws (same RNG stream)".
4. air format --check and lintr on touched R files - 0 lints.
5. tools::parse_Rd on man/bart.Rd and man/rbart.Rd - parse, ASCII-clean.
6. loo::waic smoke call if loo is installed (skipped otherwise).

## Status

- 2026-07-10: implemented on wt/pointwise-loglik; all gates green
  (tinytest 2581 passed / 0 failed = 2557 baseline + 24 new checks;
  equivalence 21/21 "identical draws (same RNG stream)", 21 compared /
  0 skipped; air format --check and lintr 0 lints on both touched R
  files; bart.Rd/rbart.Rd parse, all touched files ASCII-clean;
  loo::waic consumed the combined 400 x 100 loglik matrix of a 2-chain
  gaussian smoke fit directly).

- 2026-07-10: LANDED as 8688e12 (squash of wt/pointwise-loglik).
  Reviewer independently verified the layout claim underpinning the
  sigma pairing (combined sigma IS chain-interleaved while combined
  yhat is chain-blocked; probed empirically on same-seed 2-chain
  fits) and re-ran gates on the landed shared-library build:
  tinytest 2581/0, equivalence 21/21 identical. The pre-existing
  sampleFromPPD gaussian sigma misalignment this work surfaced is
  filed as TODO ppd-sigma-pairing.
