# extract scalar types

agent: Sonnet (R and man only; no engine code)
rng: neutral
window: after the fitted-argument-order slice lands (see Constraints)
budget: about +260/-230 over R/ (2 files), man/ (3 edited, 1 deleted),
inst/tinytest (6 files), inst/NEWS.Rd, one vignette, _pkgdown.yml, NAMESPACE,
four docs

Written as if it lives at docs/plans/extract-scalar-types.md; cites resolve
from there.

## Goal

extract returns a fit's scalar parameter draws in the shape stan4bart's and
bartCause's extract methods use. The exported draws() generic and its five
methods are gone. summary keeps its own computation internally, stops listing
the dead "tau" among its default variables, and labels a multinomial fit's
pooled probabilities prob[level].

## Context

- Ruling: dec-A68 in docs/decisions.md (TODO "extract scalar types"); dec-B99
  and dec-A24 for why posterior went. dec-A17 (a fit carries only what it
  produced) sets the absent-channel rule below.
- Today's type lists: [`extract.bart`](../../R/generics.R) ev, ppd, bart,
  loglik, trees, forest; [`extract.bartOrdinal`](../../R/generics.R) and
  [`extract.bartNegbin`](../../R/generics.R) ev, ppd, bart, loglik;
  [`extract.bartMultinomial`](../../R/generics.R) ev, ppd, bart, forest,
  loglik; [`extract.bartHurdle`](../../R/generics.R) ev, ppd, prob, bart,
  loglik. All take sample = c("train", "test") and combineChains = TRUE.
- Stored layouts, probed on the installed build. Scalars (sigma, k,
  dispersion): with combineChains = TRUE, a chain-major vector of length
  chains*samples; with FALSE, a chains x samples matrix; at one chain, always a
  vector. varcount: (chains*samples) x p, or chains x samples x p, with the
  predictor names on the p margin. A multi-forest fit's varcount has an extra
  trailing margin named forest1..forestK, and a multinomial fit's has one named
  by the levels. thresholds: (chains*samples) x (K-1), or chains x samples x
  (K-1), with no dimnames; column 1 is pinned at 0.
- Reshape helpers: [`combineChains`](../../R/bart.R),
  [`uncombineChains`](../../R/bart.R) and
  [`combineOrUncombineChains`](../../R/generics.R). The last does nothing to a
  vector and collapses only 3-D arrays, so it cannot serve the scalars or a
  4-D varcount. [`reshapeChainedChannel`](../../R/bart.R) takes a `trailing`
  margin count and covers varcount and thresholds.
- Summary internals [`presentDrawsVars`](../../R/diagnostics.R) and
  [`resolveDrawsVars`](../../R/diagnostics.R) are cited by
  docs/design/aft-variance-forest.md; keep their names.

### Shape contract (settled)

1. Orientation (agent-made, adjudicated by the orchestrator; VD approved
   chains x samples in dec-A68): vector combined, chains x samples
   uncombined, matching bartCause's BCF and dbarts' own chain-major layout.
   stan4bart differs: its sigma and k are samples x chains with dimnames
   (iterations, chain), its varcount p x samples x chains, and bartCause's
   extract.bartcFit transposes stan4bart's sigma to chains x samples.
2. Scalars (sigma, k, dispersion). combineChains = TRUE gives a vector of
   length chains*samples in chain-major order, the same row order as
   yhat.train. FALSE gives a chains x samples matrix with NULL dimnames. At one
   chain, uncombined follows extract(type = "ev", combineChains = FALSE)
   today, which returns samples x n with no chain margin (probed: 20 x 120
   at one chain, 2 x 20 x 120 at two), so a scalar is a vector and varcount
   samples x p. The result does not depend on the combineChains the fit was
   made with.
3. varcount is per predictor, so it follows yhat.train's convention, which
   bartCause's bartBCF varcount already uses: (chains*samples) x p combined,
   chains x samples x p uncombined, with the predictor names kept on the p
   margin. A multi-forest or multinomial fit keeps its trailing forest or
   level margin. Reshape through `reshapeChainedChannel` with trailing = 1, or
   2 when that fourth margin exists. Any other extension would give varcount a
   layout unlike the fit's other per-column channels.
4. thresholds are shaped like varcount (trailing = 1): (chains*samples) x (K-1)
   or chains x samples x (K-1). The pinned-at-0 first column is kept, as
   stored, and no dimnames are added.
5. A fixed k is an error, not a constant vector: a fixed-k fit carries no k
   draws (dec-A17), and stan4bart refuses the same way ("model was not fit
   with end-node sensitivity as a modeled parameter"). Wording: "cannot extract
   'k': this fit's k was fixed, not sampled". Suggest a chi hyperprior only
   where one is allowed: [`resolveNodeHyperprior`](../../R/model.R) refuses
   chi under a monotone constraint, and multi-forest fits pin k, so no
   advice is given on those fits.
6. A fit with no sigma is an error. Binary fits (probit, logistic, hazard, and
   the hurdle's occupancy part) carry none; bartCause's wording is "binary
   response model does not have a residual standard deviation parameter
   (sigma)". Wording: "cannot extract 'sigma': a <family> fit has no residual
   scale parameter". A heteroscedastic fit
   ([`fitIsHeteroscedastic`](../../R/generics.R)) stores a constant sigma
   ([`resolveDrawsVars`](../../R/diagnostics.R)'s reason), so it errors
   instead: "... a heteroscedastic fit has no scalar residual scale; its
   per-observation scale draws are the fit's 's.train'".
7. varcount is always present (keepFits = FALSE keeps it; every class
   carries it). At most a plain, untested `is.null` guard; no dedicated error.
8. sample. The scalar types, varcount and thresholds are not per-observation.
   A `sample` the caller supplies (checked with `missing(sample)`) is refused
   by name, as [`refuseTreesArguments`](../../R/generics.R) does for trees.
   The sister packages ignore it silently; refusing inert arguments by name
   is this package's rule (docs/plans/surface-refusals.md). `contribution`
   stays refused. `forest` stays refused by
   [`refuseForestSelectionOutsideForestArm`](../../R/generics.R), whose "every
   forest is already recombined into the location it reports" is false for
   these types: for sigma, k, dispersion and thresholds the text becomes
   "type = \"<type>\" is a model parameter, not a per-forest quantity", and
   for varcount "type = \"varcount\" keeps every forest on its trailing
   margin; subset that margin". A varcount forest selector is a follow-up.
9. Meaning, stated in man: on a weighted fit sigma is the scale at weight 1
   (row i's is sigma / sqrt(w_i)); on a student() fit it is the t scale, not
   the standard deviation.

### Types per class

- bart (bartBT fits share the class): + "sigma", "k", "varcount".
- bartOrdinal: + "thresholds", "varcount".
- bartNegbin: + "dispersion", "varcount".
- bartMultinomial: + "varcount". It has no scalar parameter (dec-A68), but it
  carries varcount, and "on bart fits" reads as every bart() fit. VD may veto
  this.
- bartHurdle, the "matching pieces". The components are occupancy, a probit
  fit with k sampled under its default chi(1.5, 2), and positive, a gaussian
  fit on log y with sigma present and k sampled only on request. Each has its
  own varcount.
  "sigma" is positive$sigma, the only one. "k" and "varcount" are a list
  named occupancy and positive, each element shaped per the contract
  (bartBCF's varcount, a list by forest, is the precedent); a fixed-k
  component is left out of "k", which errors when both are fixed (dec-A17).
  n.chains from [`hurdleNChains`](../../R/generics.R). See Decision.
- ordinal and negbin fits drop a modeled k when they are packaged
  ([`packageOrdinalResults`](../../R/bart.R) and
  [`packageNegbinResults`](../../R/bart.R) never read samples$k; probed with
  k = chi(1, Inf), the request is accepted and silently discarded). So "k"
  is not offered on either class. File a TODO line; it is outside this item.

## Decision (needs VD)

A hurdle fit's "k" and "varcount" are a list keyed occupancy and positive
(recommended; bartCause's bartBCF varcount is the precedent). The alternative
is prefixed types, "occupancy.k", "positive.k", "occupancy.varcount" and
"positive.varcount", matching summary's labels but adding four type names. If
VD rules for prefixes, only Step 2's hurdle arm changes.

## Constraints

- RNG neutral: extract reads stored channels and draws nothing. New fits in
  the test files go at the END of each file, or reuse fits already made, so
  the draws that later snapshot pins see do not move.
- Overlap with the concurrent fitted slice (dec-A15). Both edit
  R/generics.R: that slice edits the fitted methods and fittedForeignReasons,
  this one the extract methods, which do not conflict. In man/bart.Rd and
  man/bartBT.Rd, though, each \method{extract} usage block sits directly above
  its \method{fitted} block, and both slices rewrite the shared `sample` item.
  Branch after that slice lands, or expect to resolve conflicts in those two
  Rd files. Do not touch fittedForeignReasons or extractForeignReasons.
- draws() never reached main (`git show origin/main:NAMESPACE` has only
  extract), so there is no tombstone. There is also no NEWS entry for the
  removal (news-scope-main). The 1.0-0 NEWS items that describe draws()
  describe unshipped surface, so rewrite them in place.
- Out of scope: warmup draws (first.sigma, first.k; stan4bart's
  include_warmup), resid.df, a forest selector on varcount, k on ordinal or
  negbin, and the test-scaffolding consolidation (dec-A52; the
  dbarts:::bartDrawsArray calls in test-convergence-diagnostics.R stay).

## Steps

1. Add a scalar reshape helper to R/generics.R next to
   `combineOrUncombineChains`, with this contract: if the input is a vector,
   return it when combine is TRUE and `uncombineChains(x, n.chains)` otherwise;
   if it is a matrix, return `as.vector(t(x))` when combine is TRUE and the
   matrix otherwise. At n.chains <= 1 it returns its input unchanged.
2. Extend the five type vectors and add one early branch per class. Put it
   after validateType and the `...` refusal and before any sample or test
   channel check, so a fit with keepTrainingFits = FALSE still serves sigma.
   The branch applies contract items 5-8 and returns. It tests
   `missing(sample)` and refuses a passed sample by name before validateSample
   runs, so the resolved default never reaches it.
3. R/diagnostics.R: delete the five draws methods (and the generic at the
   top of R/generics.R); keep `bartDrawsArray`, `hurdleDrawsArray`,
   `toDrawsArray` and summary's `vars`. Drop "tau" from the four summary
   default vectors (bart, bartOrdinal, bartNegbin, bartHurdle) and "tau"/"first.tau" from
   [`scalarFields`](../../R/diagnostics.R): no shipped fit carries tau, the
   removed grouped model's group spread, so summary's table is unchanged and
   only the empty-fit message loses the name. Rewrite comments naming draws().
4. Fix [`ordinalThresholdsArray`](../../R/diagnostics.R). Its 3-D branch
   assumes (K-1) x samples x chains, but the stored layout is
   chains x samples x (K-1): probed at two chains and K = 5, summary on the
   combineChains = FALSE fit reports two "thresholds" that are chain slices.
   Permute c(2, 1, 3) and fix the comment. This also repairs
   [`plot.bartOrdinal`](../../R/plot.R), which reads the same array. The
   ordinal family is unreleased, so no NEWS entry.
5. multinomial. [`multinomialDrawsArray`](../../R/diagnostics.R) labels
   become `paste0("prob[", levels, "]")`, and
   [`summary.bartMultinomial`](../../R/diagnostics.R) sets `vars = "prob"`.
   Reword the "mean-probability" in multinomialSummaryVarsReason and the
   comments to match.
6. NAMESPACE: remove export(draws) and the five S3method(draws, ...) lines.
   _pkgdown.yml: remove the `- draws` entry under Diagnostics.
7. man: delete man/draws.Rd. [`summary.bart`](../../man/summary.bart.Rd)
   absorbs the vars rule it now defers to draws.Rd (field naming,
   heteroscedastic mean.s, per-class vocabularies), drops tau, writes
   prob[level] and points chain-separated draws at extract; its usage drops
   tau. man/bartBT.Rd and man/bart.Rd: the new types in the usage blocks and
   `type` items, with the sample refusal and contract item 9. man/bart.Rd
   carries 5 draws links in 4 Value paragraphs; rewrite all four:
   multinomial (its meanProb[level] and "sigma/k/tau summary" mentions),
   ordinal and negbin (their draws defaults quoting tau), hurdle (two links).
   Its summary usage and `vars` item drop tau.
8. Tests (item 7 needs none), draws() -> extract (callers under Verification):
   convergence-diagnostics checks extract(fit, "sigma", combineChains =
   FALSE) against t(bartDrawsArray(fit, "sigma")[, , 1]); hurdle checks
   "sigma" and names(extract(fit, "k")) == "occupancy"; nbinom "dispersion";
   ordinal "thresholds" with column 1 == 0; multinomial's draws() block becomes
   varcount dimnames, and its meanProb pins
   (["sort(sm$stats$variable)"](../../inst/tinytest/test-multinomial-generics.R)
   and the printed grepl) become prob. Appended to
   test-convergence-diagnostics.R: same-seed fits made with combineChains TRUE
   and FALSE return identical extract output for each new type under both
   settings (this catches Step 4's bug class); the one-chain vector case;
   varcount dimnames; each error in contract items 5, 6 and 8; and summary equal
   across the two ordinal fits.
9. inst/NEWS.Rd, 1.0-0 (draws never reached main: remove, never add a
   removal entry): rewrite the draws() sentences in the combineChains item
   and the summary.bart item, and delete the third passage ("draws covers
   all four families too") in the loglik item. The summary.bart item gains:
   extract offers "sigma", "k", "varcount", "thresholds", "dispersion" and
   the hurdle's, as a vector or chains x samples. Drop tau from the
   summary-methods item.
10. vignettes/gibbs_sampler_mixture_model.Rmd: the sentence that names
    draws() instead names extract(fit, "sigma", combineChains = FALSE).
11. Docs: `retired:` on the draws cite in
    docs/design/retire-grouped-random-effects.md and in
    docs/plans/interfaces-and-dependencies.md's S1 landing note, whose
    man/draws.Rd link becomes a history cite at
    62311f95c7a73cf2484e1721d91aae2100ddcb7b; amend chg-U82 and chg-U48 in
    docs/plans/bartcore-landing/changes.md; TODO closes the entry and gains
    the ordinal/negbin k gap: a requested k = chi(1, Inf) is accepted and
    sampled, then silently discarded from the fit.
12. Commits: Step 4 with its regression test lands first, alone. Steps 3 and
    5-10 then land as one commit (method deletion with NAMESPACE, the prob
    label with its test pins); Steps 1-2 go in that commit or before it.

## Verification

- draws() callers before the change, all gone after (cites marked retired
  so this file stays fresh once they are):
  - R: retired: [`draws.bart`](../../R/diagnostics.R) and four siblings, the
    generic in R/generics.R, and 4 comments.
  - tinytest, 5 calls:
    retired: ["d <- draws(fit, "](../../inst/tinytest/test-convergence-diagnostics.R),
    retired: ["adNames <- dimnames(draws(fit))"](../../inst/tinytest/test-hurdle.R),
    retired: ["d <- draws(fitCombined)"](../../inst/tinytest/test-multinomial-generics.R),
    retired: ["d <- draws(fit)"](../../inst/tinytest/test-nbinom.R),
    retired: ["adNames <- dimnames(draws(fit))"](../../inst/tinytest/test-ordinal.R).
  - man: the draws.Rd example, links in summary.bart.Rd (3) and bart.Rd
    (5, in 4 paragraphs); one vignette sentence; 3 NEWS passages;
    _pkgdown.yml; NAMESPACE (6).
  - None in inst/common, benchmarks/R, README.md, or the sister packages
    (stan4bart bartcore, bartCause and treatSens dbarts-1.0 import only
    extract). After: `git grep -n "draws(" -- R inst man vignettes` finds none.
- Probe on the slice's library: for a two-chain gaussian fit,
  identical(extract(f, "sigma"), as.vector(t(extract(f, "sigma",
  combineChains = FALSE)))). The same for varcount via
  reshapeChainedChannel, and for the combineChains = FALSE fit twin.
- Gates (neutral): tests/cpp and the full tinytest suite;
  `lintr::lint_package()` (NAMESPACE edited); `air format --check .`;
  tools/check-rc-codoc.R, check-win-drift.R and check-doc-freshness.R, each
  on its own exit status; `pkgdown::check_pkgdown(".")` (a topic removed); the
  NEWS parse gate; `R CMD check --as-cran` from a clean staged tarball (catches
  a stale draws alias). Not the exact gates: what a fit carries is unchanged.
- Downstream: bartCause's dbarts-1.0 and stan4bart's bartcore suites do not
  call draws(). Run bartCause's extract tests once against the slice's
  library.
