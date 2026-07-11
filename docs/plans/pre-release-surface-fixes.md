# pre-release-surface-fixes

agent: opus
rng: neutral for the MCMC anchor (equivalence must stay 21/21
  identical). The binary PPD fix changes seeded sampleFromPPD draws
  in one chain layout (post-hoc draws, outside the anchor) - affected
  test snapshots regenerate by whole-file replay.
budget: ~150-250 lines R + Rd + tests; dbarts.h comment-only edits.

## Context

A three-lens pre-release surface audit (2026-07-10, reports in the
session records; findings F1-F6/D2-D3) of the season's additions
found one shipping defect and a batch of freeze-regret paper-cuts.
VD approved the recommended fixes 2026-07-10.

## Fixes

1. aft loglik (THE defect, all three lenses; loo demo elpd -1389
   wrong vs -8.8 right): extract(type = "loglik") on a family "aft"
   fit scores censored rows with the gaussian density; the correct
   contribution is the log survival upper tail, as the engine's
   AFTResponse::computeLogLikelihood already does. Fix for real:
   store the 0/1 status vector on aft fit objects (bart2/bart return
   value; NULL for other families), make pointwiseLogLikelihood
   (R/generics.R) dispatch on object$family - gaussian, binary, aft
   (events dnorm, censored pnorm upper tail, both on the log-time
   scale with the paired sigma draw), and ERROR on any unrecognized
   family (future NB/ordinal cannot silently repeat this). Verify
   against the engine channel on a censored fit (agreement to
   tolerance) and against hand-computed tails.
2. Binary multi-chain PPD layout invariance (ecosystem lens): extend
   the ppd-sigma-pairing treatment to the binary families so
   combined and split-chain sampleFromPPD draws agree under a fixed
   seed, matching the gaussian guarantee landed at 03a5b85. Test
   combined == combineChains(split) for probit and weighted
   logistic.
3. growFromRoot updateState default FALSE -> NA (align with
   sampleTreesFromPrior / sampleNodeParametersFromPrior; NA =
   respect control@updateState). Rd usage line updated; bart2 passes
   FALSE explicitly so its path is unchanged.
4. state/setState trap (downstream lens): sampler$state <- x
   silently no-ops. Document setState in dbartsSampler-class.Rd
   (alias + usage + a sentence that assigning the state field does
   not restore). INVESTIGATE making direct assignment error via an
   active binding; take it ONLY if the internal storeState write
   path stays untouched and clean - otherwise Rd-only.
5. Doc batch: dbarts.h logLikelihood comment names aft (event
   density / censored survival tail; grouped composes) and
   dbarts_sampler_setState notes the encoding floor + by-name
   refusal (COMMENT-ONLY header edits, no ABI); Rd lines for aft
   ppd/prior-predictive draws being log-scale; samplePriorPredictive
   weight-blindness vs sampleFromPPD's binomial(w, p) (one sentence
   each); bart.Rd type entry cross-references the binomial(w, p)
   count scale; growFromRoot Rd rewords the linear/gp sentence so
   the hard error does not read as a silent fallback; extract
   loglik Rd documents the aperm(x, c(2, 1, 3)) recipe for
   loo::relative_eff (house chains-first convention kept).

## Verification

- Full tinytest (baseline 2675 + additions); equivalence vs
  equivalence-de67cbb.rds 21/21 IDENTICAL; air + lintr on touched R;
  pkgdown clean; tests/cpp untouched (rebuild-and-run once since the
  header comments force a recompile check).
- Poison-proof fix 1 once: feed the old dnorm-everywhere matrix to
  the new test and watch it fail.

## Status

- 2026-07-10: approved; implementation dispatched.
- 2026-07-10: implemented on wt/surface-fixes (6 commits, 5883722..8cdeee6).
  All five fixes landed. Fix 4 shipped Rd-only: an active binding on the
  state field cannot coexist with the lazy-state delayedAssign and would
  block the internal selfEnv$state <- write path identically to an external
  one (empirically confirmed), so blocking direct assignment would require a
  shadow-field refactor of the state machinery - rejected per the plan's
  zero-contortion criterion. Gates green: install --preclean OK; tests/cpp
  rebuilt clean, all pass; tinytest 2694 passed / 0 failed (2675 baseline +
  19); equivalence 21/21 identical draws (same RNG stream), 21 compared / 0
  skipped; air format --check + lintr 0 lints on all touched R; pkgdown
  check clean; dbarts.h diff comment-only (verified). Fix-1 poison: the new
  aft censored-loglik test run against the reverted dnorm-everywhere
  pointwiseLogLikelihood fails exactly the 4 censoring/family-guard checks
  (30/34 pass - event rows and gaussian/probit/logistic untouched), restored
  after.
