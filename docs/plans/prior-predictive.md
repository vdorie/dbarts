# prior-predictive

agent: sonnet
rng: neutral - no engine or sampling-path change; the new verb only calls the
  existing sampleTreesFromPrior/sampleNodeParametersFromPrior/predict entry
  points on a private sampler copy, plus R-side rnorm/rbinom/rchisq for the
  observation layer. Equivalence is byte-identical.
budget: ~90 line R function, ~100 line Rd, ~125 line test, one cross-ref
  line (task budget: 60-100/~80/~120; landed slightly over the Rd/test
  guidance, within the stated 1.5x tolerance).

## Goal

An R-level verb for drawing from the BART prior (issue #31, the most
requested item): given a sampler, repeatedly draw trees + node parameters
from the prior and evaluate the forest at x.test - optionally adding
response noise - for prior calibration before ever touching y (e.g. the
prior on a treatment effect f(x1) - f(x0)). The building blocks
(sampleTreesFromPrior, sampleNodeParametersFromPrior, predict, copy)
already exist on dbartsSampler; this composes them into one function and
never touches the caller's sampler.

## Context

sampleTreesFromPrior/sampleNodeParametersFromPrior exist on the sampler and
are exercised directly by inst/tinytest/test-sampler-prior.R and by
benchmarks/R/sbc.R's simulation-based calibration harness; what was missing
was a user-facing composition that loops draw -> predict, adds the
observation layer, and packages the result in the package's own return-
shape conventions. Decided API (per task brief, not re-litigated here):
`samplePriorPredictive(sampler, x.test = NULL, n.samples = 200L, type =
c("ev", "ppd"), offset.test = NULL, n.threads = ...)`. No naming deviation
from the requested shape.

## Design

Added to R/dbarts.R (a standalone function, not an RC method) right before
the dbartsSampler class body, alongside the other sampler-adjacent helper
warmStartState. Exported via NAMESPACE; `qchisq`/`rchisq` added to the
stats importFrom line (previously unused by the package).

Isolation: operates on a freshly constructed private sampler,
`dbartsSampler$new(newControl, sampler$model, sampler$data)`. It gets its
own external pointer and engine state; it borrows the caller's data object,
which is safe because the verb never mutates predictors (the only in-place
mutation path, setPredictor, is never called). Confirmed empirically:
predict() and getSigmas() on the original sampler are bit-identical before
and after a samplePriorPredictive call (test (b) below), including for a
keepTrees + already-run sampler.

DEFECT CAUGHT IN REVIEW (fixed before landing): the first cut used
`sampler$copy(shallow = FALSE)`. copy() installs the donor's cached state
via setState - INCLUDING the engine RNG state - and the caller's sampler
never advances, so every call replayed the same frozen stream: two
consecutive no-seed calls returned bit-identical matrices (a user drawing
200 samples then 200 more would silently get the same 200 forests twice).
The fix drops copy() for the fresh construction above, which skips
setState; no donor state is needed because the prior draws overwrite all
tree state, sigest lives on data@sigma, and predict() reads the live trees
just drawn. Fresh creation seeds the chain RNGs per
src/R_interface_bartcore.cpp createChainRngs (L955-989): with
control.haveRngSeed a dedicated Mersenne-twister seed generator hands each
chain its seed (pinned, reproducible without R's stream); without one, the
no-seed branch (L973-982) draws each chain seed from R's stream via
GetRNGstate() + `unif_rand() * 4294967295.0` + PutRNGstate(), so a
set.seed() beforehand governs the engine draw and successive calls are
independent by default.

keepTrees handling: if the private sampler inherited control@keepTrees =
TRUE, predict() would serve saved posterior samples instead of the live
trees this function just drew (bartcore_predict: capacity > 0 reads saved
trees, capacity == 0 reads live trees). The function forces keepTrees off
on the control before constructing the private sampler, which also skips
the saved-tree storage allocation.

Multi-chain samplers: sampleTreesFromPrior loops over every chain
internally (sampler.hpp), so a multi-chain sampler still advances every
chain's engine RNG each iteration. Prior draws are chain-free (no
mixing/convergence concept applies to independent prior draws), so the
function keeps only the first chain's predict() column each iteration and
discards the rest - the natural single-stream path, at the cost of wasted
compute on the unused chains. Documented in the Rd's Details.

x.test default: NULL evaluates at `draw$data@x` (the sampler's own training
predictors, already validated/typed); an explicit x.test is passed straight
through predict(), which does its own validateXTest.

type = "ev": the forest sum on the response scale for a gaussian sampler
(predict()'s own scale, per R_interface_bartcore.cpp's bartcore_predict
comment), or `probabilityFromLatents(result, list(family = ...))` (the
same helper extract.bart/predict.bart use) for a binary one - reusing
`draw$model@family` (already resolved to "probit"/"logistic" by dbarts(),
never "auto" for a real binary sampler).

type = "ppd" sigma draw (gaussian only; binary reuses the ev probabilities
through rbinom): dispatches on `draw$model@resid.prior`'s class, the only
two concrete dbartsResidPrior subclasses in the codebase.
  - dbartsChiSqPrior: the reported-scale scaled-inverse-chi-squared verified
    against benchmarks/R/sbc.R's sbcSigmaDraw (mirrored inline, not
    sourced - the package must not depend on benchmarks/): sigma^2 ~ df *
    sigest^2 * rawScale / chisq(df), rawScale = qchisq(1 - quantile, df) /
    df, sigest = draw$data@sigma, df/quantile = resid.prior@df/@quantile.
  - dbartsFixedPrior: rather than re-deriving fixedSigmaSq's scale (it is
    stored as an internal-scale sigma^2 in ParsedModel, per
    R_interface_bartcore.cpp - not obviously the reported scale), the
    function reads `draw$getSigmas()[1L]` once and repeats it for every
    sample. getSigmas() ("current residual error term on original,
    standard deviation scale") already does the internal-to-reported
    conversion (chain.hpp: `sigma_ * response_->sigmaScale()`) regardless
    of which prior set sigma_, and a fixed prior never updates sigma_ after
    construction, so this is exact and sidesteps the scale question
    entirely. Verified empirically: fixed(0.25) sampler reports
    getSigmas() == 0.5 == sqrt(0.25), and ppd - ev variance over 320 draws
    landed at 0.243 (target 0.25).
  - any other dbartsResidPrior subclass: stop() naming the class, rather
    than guess at an unverified scale convention.

Return shape: n.samples rows x observations columns, no chain dimension -
verified against extract.bart's own combined-chains convention
(convertSamplesFromDbartsToBart + combineChains: n.samples x n.obs for a
single chain, or n.samples x n.chains x n.obs uncombined - prior draws have
no second axis to combine, so the function always returns the combined
2-D shape directly via `do.call(rbind, results)`).

Seeding semantics (documented in the Rd Details, tested in (c)): with
rngSeed NA on the control, successive calls are independent by default and
set.seed() reproduces the whole draw (engine forests via createChainRngs'
no-seed branch, observation layer via rnorm/rbinom/rchisq); with a
control-fixed rngSeed, the engine stream is pinned and every call replays
the identical forest sequence, only the R-side ppd layer still following
set.seed().

Cross-call pairing note (Rd Details + example): two separate
samplePriorPredictive calls draw two independent sets of forests, so
f(x1) minus f(x0) computed from two separate calls is meaningless unless
rngSeed is fixed on the control. Stacking both settings as rows of one
x.test pairs the forests within a single call instead (used in the Rd
example) - simpler and does not depend on the rngSeed subtlety.

## Test

inst/tinytest/test-prior-predictive.R, 17 checks: (a) shapes for ev/ppd x
gaussian/binary, both default (x.test = NULL) and an explicit x.test,
binary ev draws in [0, 1] and ppd draws in {0, 1}; (b) the caller's
sampler's predict()/getSigmas() bit-identical before and after; (c)
seeding semantics - two no-seed calls NOT identical (the frozen-replay
regression check), set.seed() before each of two calls -> identical (both
families; real only because the engine seed comes off R's stream), and a
control-fixed rngSeed -> two ev calls identical with no set.seed; (d) a
scale-robust statistical check: aggregate ppd variance >= ev variance on
the PINNED sampler, whose replayed stream pairs the forests across the two
calls, making the inequality structural (ppd = the same forests + noise)
rather than a comparison of independent forest sets - the unpinned form
flaked exactly that way when first run; (e) not-a-sampler and bad-type
validation errors. Manually
verified beyond the automated file: multi-chain (n.chains = 3) still
returns n.samples x n.obs with no chain axis; a fixed resid.prior sampler's
ppd variance lands near its fixed sigma^2; the Rd example runs clean under
pdf(NULL).

## Verification

Gates (worktree, private library
/Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-land):
1. R CMD INSTALL -l <lib> . (R-only change, no --preclean) - DONE, clean.
2. tinytest::test_package("dbarts") - 2547 passed / 0 failed (2530 baseline
   + 17 new checks).
3. Equivalence vs benchmarks/baselines/equivalence-de67cbb.rds - all 21
   scenarios "identical draws (same RNG stream)", 21 compared / 0 skipped.
4. air format --check on R/dbarts.R, NAMESPACE,
   inst/tinytest/test-prior-predictive.R - clean.
5. lintr::lint() on R/dbarts.R and inst/tinytest/test-prior-predictive.R -
   0 lints each.
6. man/samplePriorPredictive.Rd and the dbartsSampler-class.Rd seealso edit:
   ASCII-clean (grep for non-ASCII bytes, empty), tools::parse_Rd +
   Rd2txt render cleanly, example code run manually - executes without
   error or warning.

## Status

- 2026-07-10: landed on wt/prior-predictive. samplePriorPredictive added to
  R/dbarts.R, exported; man/samplePriorPredictive.Rd added; a seealso line
  added to dbartsSampler-class.Rd; inst/tinytest/test-prior-predictive.R
  added. All gates green (verbatim above).
- 2026-07-10: review caught the frozen-state replay defect (copy() +
  setState reinstalls the engine RNG, so no-seed calls repeated
  identically); fixed by fresh construction without setState, seeding
  semantics re-documented in the Rd, regression check added to the tests
  (17 checks now), and all gates re-run green. Amended into the single
  landing commit.
