# bcf-testfits-guard

agent: opus
rng: neutral (no draw touched; only the previously-garbage BCF testFits
     recording changes, plus new bridge Rf_errors that fire for BCF only.
     equivalence exact 18/18, bcf-exact still matches)
budget: engine ~20 lines, bridge guard ~12 lines + 4 call sites, docs,
        tests/cpp BCF case + one tinytest pair

## Goal

Under BCF (two forests) the recorded test-fit and out-of-sample
prediction surface is single-forest-only, so it silently misreports: it
records the bare prognostic (mu) forest with the combined-response
centering, with no a, no b_z tau, and no test treatment vector to blend
with. Stop emitting that garbage - NaN the recorded channel and refuse
the bridge entries that would reach it - and document the mu-forest-only
diagnostics.

## Context

- Defect: correctness-audit block-6/7 finding (latent). In
  src/bartcore/chain.hpp storeSample, the trainingFits branch carries the
  BCF blend `scale * (a mu + b_z tau) + shift`, but the testFits branch
  unconditionally recorded `scale * forests_[0].totalTestFits + shift`
  (+ testOffset) - the bare mu forest.
- A correct test blend is ill-defined: the API carries no test treatment
  vector. BCF prediction recombines per forest through getForestFits
  (facade forestTotalFits) + the BCF glue (bcfGlue), docs/design/bcf.md.
- Reachable: bartcore_setTestPredictor (and setTestOffset,
  setTestPredictorAndOffset) had no forest-count guard, so a BCF sampler
  pointer handed test predictors silently recorded garbage; a dbarts.h
  consumer would too. The facade predict() path
  (predictFromSavedSample / predictFromCurrentTrees) also sums only
  forests_[0], so bartcore_predict on a BCF pointer had the same defect.
- Creation-path confirmation: createBCFHolder
  (src/R_interface_bartcore.cpp ~1185) wires NO test data - it never calls
  setTestPredictors, unlike the non-BCF createHolder (~1171). Any
  creation-time x.test in the parsed data is parsed but ignored. So the
  only way to get test data onto a BCF sampler was the setters, now
  guarded.

## Fix

1. Engine (chain.hpp storeSample testFits branch): when bcf_ != nullptr,
   fill the sample's testFits slots with quiet_NaN() and skip the offset
   add; the non-BCF path is unchanged (the else branch is the original
   code verbatim). Added <limits>.
2. Bridge (R_interface_bartcore.cpp): new helper refuseBCFTestSurface
   (anonymous namespace, beside refuseViewSampler), Rf_error when
   numForests() >= 2 - the same predicate bartcore_setTreatment uses -
   with a message pointing at getForestFits + the BCF glue. Called from
   bartcore_setTestPredictor, bartcore_setTestOffset,
   bartcore_setTestPredictorAndOffset (each after its null-clear early
   return, so clearing test data stays a harmless no-op), and
   bartcore_predict (the out-of-sample analog, same single-forest defect).
3. Docs: the Results struct doc (chain.hpp) now states that under BCF
   testFits is NaN and k / variableCounts / splitProbabilities report the
   prognostic (mu) forest only; the treatment forest is reached through
   the per-forest channels.

Deviation from the literal spec: item 1 asked to check the facade
test-fits accessor and "treat it consistently". The equivalent path is
the facade predict() (forests_[0]-only). Rather than NaN it engine-side
(it is a general single-forest primitive, and a forest-selective BCF
predict could reuse it later), I guarded it at the bridge
(bartcore_predict), consistent with the setTestPredictor guard and the
"defense at both levels" framing. forestTotalFits / bcfGlue are the
documented per-forest channels and were left untouched. The R BCF surface
never calls bartcore_predict on a BCF sampler (it uses bartcoreForestFits
+ bartcoreBCFGlue), so nothing R-side regresses.

## Verification

- tests/cpp testBCFTwoForest (test_sampler.cpp) extended: sets test
  predictors before the run, asserts every recorded testFit is NaN
  ("BCF test fits come back NaN"), and asserts the final recorded
  trainingFits equals the fitScale * (a mu + b_z tau) + shift blend by
  recovering the shift from one row and confirming the linear structure
  over the rest ("BCF training fits are the a*mu + b_z*tau blend") - which
  would fail for the bare mu forest.
- tinytest test-bcf.R: two expect_error calls assert bartcoreSetTestPredictor
  and bartcorePredict both error (message "BCF") on a BCF sampler.

## Status (2026-07-08, LANDED)

Gates (from the worktree):
- R CMD INSTALL --preclean . : clean (DONE).
- tests/cpp: rm binaries, make clean && make && ./test_bartcore ->
  all tests passed, including "ok: BCF two-forest sampler" and
  "ok: BCF fixed glue".
- tinytest::test_package("dbarts") -> 2498 pass / 0 fail (2496 baseline
  + the 2 new BCF error expectations).
- equivalence compare equivalence-cf99a00.rds -> 18 compared / 0 skipped,
  every arm "identical draws (same RNG stream)"; "posteriors statistically
  indistinguishable". No divergence.
- bcf-exact.R -> "OK: BCF sampler matches the exact posterior" (consumes
  the per-forest channels, unaffected).
