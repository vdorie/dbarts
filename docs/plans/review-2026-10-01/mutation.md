# Review 3 - mutation leg (test strength)

Tree: wt/r3-mut at 01dee4b4; battery extended and committed as 4a0f92e7
(benchmarks/R/mutation-battery.R, entries m28-m89). Not rebased onto
8e7a3d19; the one site riding on code it rewrote (m40) was dropped.

Method. 62 one-token defects planted in code changed since 7ad0bbea, chosen by
risk. Each was applied to a `git archive HEAD` copy, built --preclean, checked
fresh, and run against its designated killers (owning tinytest files, tests/cpp,
cheap quick exact gates). Every survivor was then replanted against the FULL
tinytest suite (202 files; pristine build 0 failing), tests/cpp, and ten cheap
quick gates (hazard-exact, hazard-reduction, hurdle-exact, hurdle-reduction,
negbin-exact, ordinal-exact, t-exact, aft-exact, multinomial-exact,
monotone-reference). verify-anchors clean before and after.

Not covered: the 20-minute monotone-exact-enumeration gate, equivalence.R, the
reproducibility snapshot files (they exit off the reference build), the flat C
API consumer build.

Totals. 62 run. 42 killed by designated killers; 4 more (m50, m51, m81, m88)
killed only by the replant (full suite or a gate) - their killers were added to
the battery and they re-ran KILLED. 16 survived everything: 2 equivalent (m36,
m55), 1 dropped (m40), 13 confirmed gaps, now SURVIVE_DOCUMENTED.

Observation. tests/cpp carries most of the monotone engine: m29, m32-m35, m44
die only there; test-monotone.R and the quick monotone gates miss all of them.

## Per site

Killers: tt = tinytest file, cpp = tests/cpp, gate = quick exact gate.

| id | site | mutation | killed by |
|---|---|---|---|
| m28 | MonotoneLeafGeometry::relation | `hi + 1 == lo` -> `hi == lo` | tt monotone, cpp, gate monotone-reference |
| m29 | MonotoneLeafGeometry::share | `<=` -> `<` on interval overlap | cpp |
| m30 | MonotoneLeafGeometry::build | missing-side flag `!=` -> `==` | tt monotone, cpp |
| m31 | MonotoneLeafGeometry::share | factor level overlap `!= 0` -> `== 0` | tt monotone, cpp |
| m32 | monotoneLogExtensions | forward count `+=` -> `=` | cpp |
| m33 | monotoneLogPositionLaw | drop backward exponent | cpp |
| m34 | monotoneLogAdjacentShare | binomial `i+j-2` -> `i+j-1` | cpp |
| m35 | monotoneLogPairRatio | m drops second.size | cpp |
| m36 | monotoneMovePair | flip `< 0` -> `> 0` | EQUIVALENT (see below) |
| m37 | decideNormalizedMove | Z_T ratio sign swapped | cpp, gates monotone-reference + successive-conditional |
| m38 | monotoneDrawPriorLeaves | isolated leaf draws constrainedSd | SURVIVED full suite |
| m39 | MonotoneConstantGaussianLeaf::priorSd | c-inflation on free leaves | cpp, gate monotone-reference |
| m40 | redrawAfterBirth | floor `max(aU, aL)` -> `aU` | DROPPED (rides on 8e7a3d19's rewritten inversion) |
| m41 | redrawAfterBirth | lower child not capped by sibling | tt monotone, cpp |
| m42 | jointPriorAccepts | cone test `>` -> `<` | SURVIVED full suite |
| m43 | Chain ctor | joint/leaf switch `== 1` -> `== 0` | tt monotone, cpp, gate monotone-reference |
| m44 | Chain::sampleTreesFromPrior | joint accept under leaf | cpp |
| m45 | Chain::setForestMapSd | leaf scale not re-derived | tt multiforest-leaf-prior-writer |
| m46 | writeForestSpreads (R) | forest index 1-based | tt multiforest-leaf-prior-writer |
| m47 | writeForestSpreads (R) | mirror channel `> 0` -> `< 0` | tt multiforest-leaf-prior-writer |
| m48 | Chain::setForestFixedK | `!=` -> `==`, k never written | tt multiforest-leaf-prior-writer, multinomial-r5-surface |
| m49 | resolveHazardGrid | `left.open = FALSE` | tt hazard-factors, gate hazard-exact |
| m50 | expandDiscreteTimeHazard | offset by period not subject | full suite: tt family-mutation-parity (now designated) |
| m51 | hazardSurvivalProbabilities | horizon `<=` -> `<` | gate hazard-exact only (now designated); full tinytest clean |
| m52 | hazardSurvivalProbabilities | period column `each` -> `times` | tt hazard, hazard-factors |
| m53 | combineHurdleChannel | ev drops the 0.5 | tt hurdle, gate hurdle-exact |
| m54 | hurdleLogLik | Jacobian `-` -> `+` | tt hurdle |
| m55 | AFTResponse::setSurvivalStatus | drop `logT_[i] = observed` | EQUIVALENT (see below) |
| m56 | AFTResponse::redrawCensored | variance used as sd | tt aft-heteroscedastic, gate aft-exact |
| m57 | multinomial PG draw | `c < trials` -> `<=` | gate multinomial-exact only |
| m58 | composeEffectiveRows | zero-trial rows never masked | tt multinomial-zero-trials |
| m59 | resolveMultinomialCounts (R) | numeric codes not shifted to 1-based | SURVIVED full suite |
| m60 | multinomialLogLik (R) | drop multinomial coefficient denominator | tt multinomial-generics, multinomial-zero-trials |
| m61 | NBDispersionPrior::computeKernel | histogram ignores mask | tt active-rows-pins, cpp |
| m62 | NBResponse::setActiveRows | kernel not rebuilt over subsample | tt active-rows-pins, cpp |
| m63 | NBResponse::computeLogLikelihood | `logOnePlusExp(-eta)` -> `(eta)` | SURVIVED full suite |
| m64 | negbinLogLik (R) | `each` -> `times` | tt nbinom |
| m65 | TResponse::refreshLatents | lambda ignores weight | cpp only (t-exact quick clean) |
| m66 | TResponse::refreshLatents | nu stats count masked rows | cpp only |
| m67 | TResponse::computeLogLikelihood | `/ sqrt(w)` -> `/ w` | SURVIVED full suite |
| m68 | ordinalThresholdLogAcceptance | top gap `<` -> `<=` | tt ordinal, gate ordinal-exact, cpp |
| m69 | ordinalThresholdLogAcceptance | masked rows in acceptance | tt active-rows-pins, cpp |
| m70 | ordinalLogLik (R) | `each` -> `times` | tt ordinal |
| m71 | runPredictorTransaction | gathered raw not restored on reject | SURVIVED full suite |
| m72 | runPredictorTransaction | no repartition on reject | cpp only |
| m73 | runPredictorTransaction | factor precheck only with cut refresh | cpp only |
| m74 | SubsetUpdate::restore | hasMissing not restored | SURVIVED full suite |
| m75 | ColumnStore::categoricalValueIsValid | `<` -> `<=` | tt data-categorical-declared, cpp |
| m76 | inferredCategoryCountCsc | drop `splitsBySubset(j) &&` | SURVIVED full suite |
| m77 | SparseColumnData::at | rank counts own bit | tt data-sparse, sparse-factor, cpp |
| m78 | rowsWithMissingPredictors (R) | dgC row index not +1 | SURVIVED full suite |
| m79 | rowsWithMissingPredictors (R) | mixed: dense NA rows dropped | SURVIVED full suite |
| m80 | mapFactorColumnsToTrainingLevels (R) | drop `levels = factorLevels[[j]]` | SURVIVED full suite |
| m81 | remapSparseFactorToTrainingLevels (R) | `<` -> `<=` | full suite: tt sparse-factor-na (now designated) |
| m82 | validateResponseSupport | ordinal `< 1` -> `< 0` | tt ordinal |
| m83 | validateResponseSupport | nbinom integrality dropped | tt nbinom |
| m84 | refuseNonBinaryMask | `!= 1` -> `> 1` | SURVIVED full suite |
| m85 | validateCategoryOffset | dim `[1]` -> `[0]` | tt multinomial-category-offset, multinomial-test-offset |
| m86 | predictFromSource | draw axis = capacity | tt tree-store-order |
| m87 | posteriorInterval (R) | upper `/2` dropped | SURVIVED full suite |
| m88 | reshapeScalarChannel (R) | `t(x)` dropped | full suite: tt convergence-diagnostics (now designated) |
| m89 | resolveOrdinalResponse (R) | `+ 1L` dropped (control) | tt ordinal |

Equivalent. m36: a split on a constrained axis always leaves its two children
in one component (adjacent on the axis, identical elsewhere), where theta =
e(C0)/e(C*) is symmetric in the pair; on a free axis direction is 0 and flip is
false either way. m55: setSurvivalStatus is reached only from
bartcore_setResponse, which memcpys the new y over logT_ immediately after, so
the restore is a dead store on every reachable path.

## Confirmed gaps and the test that would kill each

Probes on the pristine build confirm m59, m78, m80, m87 change observable
output (m59: y = c(0,1,2,2) codes rows 1-4 as 0-trial, cat 0, cat 1, cat 1
instead of cat 0, 1, 2, 2; m78: NA at row 3 flags row 2; m80: test factor
c("b","c") codes 1,2 instead of 2,3; m87: ci.upper 1.74 instead of 2.01 on
N(0,1) draws).

- m38 isolated monotone leaf drawn at the c-inflated sd. tests/cpp/test_monotone.cpp:
  draw prior leaves many times on a tree with a free split (one isolated leaf);
  assert the isolated leaf's sample sd is scale/k (not cInflation*scale/k).
- m42 joint prior cone test inverted. tests/cpp/test_monotone.cpp: on a stump
  split on a monotone axis, jointPriorAccepts over many draws accepts with rate
  1/2 (two related iid leaves); the mutant accepts 1/2 too only by symmetry, so
  assert on a 3-leaf chain (rate 1/6) or assert every accepted draw is in the
  cone via monotoneTreeIsFeasible.
- m59 numeric multinomial codes. inst/tinytest/test-multinomial-surface.R: fit
  with y = c(0, 1, 2, ...) numeric and the same as factor(y); assert identical
  fits (same seed) or that the count matrix has one 1 per row.
- m63 engine nbinom loglik. tests/cpp/test_model.cpp: NBResponse
  computeLogLikelihood at a few (y, eta, r) equals lgamma(y+r) - lgamma(r) -
  lgamma(y+1) + y log p + r log(1-p), p = plogis(eta); or the C API consumer
  (inst/tinytest/capi) compares results.logLikelihood to dnbinom.
- m67 engine t loglik with weights. tests/cpp/test_model.cpp: TResponse with
  weights w != 1; assert out[i] = log dt((y-mu)/s, nu) - log s with
  s = sigma*range/sqrt(w).
- m71 rejected predictor update on a leaf-covariate sampler. inst/tinytest/test-linear-leaves.R:
  $setPredictor that empties a leaf (rolledBack) on a linear-leaf sampler, then
  run; assert fits equal a twin that never attempted the change.
- m74 subset rollback missingness. tests/cpp/test_sampler.cpp (or test_data.cpp):
  updatePredictor on one column introducing NA where the change is rejected;
  assert data.hasMissing[j] equals its pre-call value.
- m76 sparse ordered-factor level count. inst/tinytest/test-sparse-factor.R:
  an ordered sparseFactor whose reference level code is above every stored
  code; assert getCutPoints reports the stored levels' count (not one more).
- m78 sparse NA row index. inst/tinytest/test-sparse-factor-na.R (or
  test-na-action.R): dgCMatrix with an NA at row 3, na.action = na.omit; assert
  the fit's na.action names row 3.
- m79 mixed container NA rows. inst/tinytest/test-na-action.R: mixed container
  with a dense NA at row i and a sparse block without NA; na.omit drops row i.
- m80 test factor with a level subset. inst/tinytest/test-data-categorical.R:
  predict on newdata whose factor holds only some training levels (its own
  levels() a subset); assert predictions equal those with the full level set.
- m84 fractional active-row mask. inst/tinytest/test-active-rows-pins.R:
  expect_error($setActiveRows(c(0.5, rep(1, n - 1))), "exactly 0 or 1").
- m87 interval upper bound. inst/tinytest/test-generics-intervals.R: assert
  ci.upper equals quantile(draws, (1 + ci.level) / 2) for one observation.

Weaker spots worth noting (killed, but narrowly): m51 dies only in the
hazard-exact gate (no tinytest pins S(t) at a grid time); m57 only in
multinomial-exact; m65, m66, m72, m73, m44 only in tests/cpp; t-exact quick
missed both t lambda/nu mutations.
