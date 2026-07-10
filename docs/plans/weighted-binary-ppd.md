# weighted-binary-ppd

agent: sonnet
rng: neutral for the sampler - no sampling-path changes, so equivalence is
  byte-identical for every scenario. The ppd post-processing draws themselves
  do change, but only for weighted-binary fits (a new RNG call shape,
  rbinom(, size, prob) instead of rbinom(, 1, prob) followed by a multiply);
  unweighted binary, gaussian, and weighted-gaussian ppd draws are untouched.
budget: ~15 lines net in generics.R, ~90 lines of tests, one doc sentence.

## Goal

Fix the last live defect of issue #79: dbarts's binary families already
handle weights correctly at the sampler level (probit refuses weights;
logistic implements exact integer frequency weights, verified equivalent to
physical row replication), but the posterior-predictive draw path drew a
degenerate two-point distribution for weighted binary rows instead of the
coherent binomial posterior predictive.

## Context

sampleFromPPD (R/generics.R) is the single function behind every ppd path:
predict.bart, extract.bart, predict.rbart, and extract.rbart all call it
(fitted.bart and fitted.rbart delegate to extract, so they inherit the same
fix without their own call site). For a weighted binary row, the old code
drew `rbinom(length(ev), 1L, ev)` - one bernoulli trial per posterior draw -
and then multiplied by the row's weight, landing on {0, w} only. A weight-w
row is w iid bernoulli trials at the fit's per-trial probability (issue #79's
sampler-level semantics, already exact and unchanged here); the coherent ppd
draw is the number of successes, binomial(w, p), support 0..w. The old
draw's mean (w * p) was coincidentally correct - it's what made the defect
easy to miss - but its intervals and quantiles were wrong. rbart_vi has no
logistic family (binary rbart_vi fits are always probit, which refuses
weights outright), so the weighted-binary branch is only reachable through
bart / bart2 with family = "logistic" and integer weights.

## Fix

R/generics.R, sampleFromPPD, the weighted responseIsBinary branches only
(unweighted branches and the weighted-gaussian branch, which already
implements sigma / sqrt(w) precision semantics correctly, are untouched):

- 2-d ev (samples x obs, obs in columns): replace
  `rbinom(length(ev), 1L, ev)` + `t(weights * t(result))` with a single
  `rbinom(length(ev), size, ev)` where
  `size <- rep(weights, each = nrow(ev))` - the same per-column recycling
  the old code's transpose dance achieved, now supplied as the binomial size
  instead of a post-hoc multiplier.
- 3-d ev (obs in the last dimension): same idea, with
  `size <- rep(weights, each = prod(dim(ev)[1:2]))`, replacing the
  aperm-multiply-aperm recycling.
- rbinom's return type (integer) is unchanged from the pre-existing
  unweighted branch, so no new storage-type inconsistency is introduced.
- object$seed save/restore at the top and bottom of sampleFromPPD is
  untouched.

man/bart.Rd's weights item gains one sentence documenting the ppd semantics
for a weighted logistic fit (binomial(w, p)); rbart.Rd needs no change since
rbart_vi's binary fits are always probit and probit refuses weights, so the
weighted-binary ppd draw is unreachable there.

## Test

inst/tinytest/test-weighted-binary-ppd.R, five bart2(family = "logistic")
fits:

1. Main fit, weights = rep_len(c(1, 3, 5), n), moderate slope (p away from
   0/1), 500 samples: every column's draws stay in [0, w]; every w = 5
   column hits at least one intermediate value (1..4) across the 500 draws;
   colMeans(ppd) ~= w * colMeans(ev) at a loose tolerance (this held even
   under the old code - it only guards the fix against a prob/size mixup);
   a set.seed/extract pair repeated verbatim agrees exactly.
2. A dedicated small fit with weights = c(5, 1, 1, ..., 1): only column 1
   may exceed 1; this specifically guards the rep(weights, each = nrow(ev))
   recycling direction (a transposed or mis-recycled size vector would leak
   the weight-5 support into the wrong column).
3. An unweighted fit: ppd draws stay in {0, 1}, confirming the untouched
   branch is unaffected.
4. A two-chain fit with extract(..., combineChains = FALSE): the 3-d path
   (obs in the last dimension) is reachable through this public surface: dim
   length 3, correct obs-dimension extent, and 0..w support checked once for
   shape/support rather than statistical power.

## Verification

Gates, worktree at .claude/worktrees/weighted-binary-ppd, private library
/Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-land:

1. R CMD INSTALL -l <lib> (R-only, no --preclean) - DONE (dbarts).
2. tinytest::test_package("dbarts") - 2557 passed / 0 failed (2547 baseline
   + 10 new checks).
3. Equivalence vs benchmarks/baselines/equivalence-de67cbb.rds - all 21
   scenarios "identical draws (same RNG stream)", 21 compared / 0 skipped.
4. air format --check and lintr on both touched R files (R/generics.R,
   inst/tinytest/test-weighted-binary-ppd.R) - clean, 0 lints.
5. tools::parse_Rd on man/bart.Rd - parses; file remains ASCII-only.

## Status

- 2026-07-10: fix + tests + doc landed on wt/weighted-binary-ppd; all gates
  green (2557/0 tinytest, 21/21 identical-draws equivalence, 0 lints).
