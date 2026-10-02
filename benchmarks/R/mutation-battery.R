#!/usr/bin/env Rscript

# Repeatable package-scale mutation harness (docs/plans/release-candidate-
# review.md wave 3). Seeded from the July poison sweep (docs/plans/
# gate-blindspot-audit.md "## Status", POISON SWEEP table) and its follow-up
# gates (docs/plans/gate-hardening-1.0.md), extended to R/ and to three
# SURVIVE_DOCUMENTED entries from the P7 branch-reach feed (an untracked
# session file). Every site
# below was RE-DERIVED against the current tree by symbol, not copied from
# the July file:line anchors, which are months stale.
#
# NEVER touches the working tree: each mutation is applied to a fresh
# `git archive HEAD` copy, built with a --preclean install into one lib
# shared across the run, verified fresh (tools/check-build-freshness.R, the
# wave-0 stale-install guard), run against its designated killer(s), then
# discarded. KILL_EXPECTED entries must be caught by at least one killer;
# SURVIVE_DOCUMENTED entries must be caught by none (they document a real,
# currently-untested gap rather than assert one that doesn't exist).
#
# Modes:
#   list                    print the inventory
#   verify-anchors          resolve every anchor against the CURRENT tree,
#                           no build - the cheap drift guard
#   run all|id[,id...]      apply, build, and gate the selected entries
#     --keep-going          keep processing entries after one comes back
#                           wrong or errors, instead of stopping there
#
# Scratch build dirs default to tempdir(); override with the
# MUTATION_BATTERY_SCRATCH_DIR env var. Equivalence killers scope to one
# scenario via EQUIVALENCE_SCENARIOS so a mutation run stays minutes, not the
# full battery of scenarios equivalence.R carries.

`%||%` <- function(a, b) if (is.null(a)) b else a

scriptDir <- dirname(sub(
  "--file=",
  "",
  grep("--file=", commandArgs(), value = TRUE)
))
repoRoot <- normalizePath(file.path(scriptDir, "..", ".."))
# kEquiv's compare stays statistical (no --bitwise): the battery installs the
# shipped build, which owes the reference-build baselines only the
# statistical match, and its equivalence killers catch a posterior shift.
equivBaseline <- "benchmarks/baselines/equivalence-54d8054b.rds"

## ---- mutation-list constructors -------------------------------------------

# `text` serves both as the drift-detection anchor and the exact span
# replaced by `mutant`; keeping them the same string is what makes
# verify-anchors a faithful preflight for run's own substitution.
mk <- function(id, file, text, mutant, class, killers, note) {
  list(
    id = id,
    file = file,
    anchor = text,
    original = text,
    mutant = mutant,
    class = class,
    killers = killers,
    note = note
  )
}
kScript <- function(path, args = character(0), env = NULL) {
  list(list(argv = c("Rscript", path, args), cwd = ".", env = env))
}
kEquiv <- function(scenario) {
  list(list(
    argv = c("Rscript", "benchmarks/R/equivalence.R", "compare", equivBaseline),
    cwd = ".",
    env = paste0("EQUIVALENCE_SCENARIOS=", scenario)
  ))
}
kCpp <- function() {
  list(list(
    argv = c("sh", "-c", "make -j4 && ./test_bartcore"),
    cwd = "tests/cpp"
  ))
}
kTinytest <- function(testFile) {
  expr <- sprintf(
    'suppressPackageStartupMessages(library(dbarts));res<-tinytest::run_test_file("%s",verbose=0);q(status=if(length(res)>0 && all(as.logical(res))) 0L else 1L)',
    testFile
  )
  list(list(argv = c("Rscript", "-e", expr), cwd = "."))
}

## ---- the mutation list -----------------------------------------------------
## KILL_EXPECTED entries m01-m16 are the 16 poisons of the July sweep, all
## still live (none retired: every site below was re-derived by symbol and
## still encodes the same semantic breakage, though several moved file or
## function under refactors - e.g. the BCF glue left chain.hpp for the new
## combiner.hpp). m17-m20 extend the battery to R/. m21-m23 are the
## SURVIVE_DOCUMENTED trio. m24-m25 are the perturb kernel's two poisons and
## m26-m27 the rule_gibbs kernel's.

mutations <- list(
  mk(
    "m01",
    "src/bartcore/model.hpp",
    "return base / std::pow(1.0 + static_cast<double>(tree.depthOf(nodeIndex)), power);",
    "return base / std::pow(1.0 + static_cast<double>(tree.depthOf(nodeIndex) + 1), power);",
    "KILL_EXPECTED",
    kScript("benchmarks/R/bd-balance.R"),
    "poison 1: CGM growth probability off-by-one on depth (was model.hpp:1398)"
  ),

  mk(
    "m02",
    "src/bartcore/moves.hpp",
    paste0(
      "    Node oldNode = tree.at(nodeToChange);\n",
      "    tree.orphanChildren(nodeToChange);\n\n",
      "    BranchScore newScore =\n",
      "      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);\n",
      "    double transitionProbabilityOfBirthStepReverse =\n",
      "      probabilityOfBirthStep(ctx, tree, true);\n",
      "    double reverseTransitionProbabilityOfSelectingNodeForBirth =\n",
      "      probabilityOfSelectingNodeForBirth(ctx, tree);"
    ),
    paste0(
      "    Node oldNode = tree.at(nodeToChange);\n",
      "    double reverseTransitionProbabilityOfSelectingNodeForBirth =\n",
      "      probabilityOfSelectingNodeForBirth(ctx, tree);\n",
      "    tree.orphanChildren(nodeToChange);\n\n",
      "    BranchScore newScore =\n",
      "      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);\n",
      "    double transitionProbabilityOfBirthStepReverse =\n",
      "      probabilityOfBirthStep(ctx, tree, true);"
    ),
    "KILL_EXPECTED",
    kScript("benchmarks/R/bd-balance.R"),
    "poison 2: death's reverse node-selection count moved onto the pre-death tree (was moves.hpp ~245); gate-hardening-1.0 sub-item 5 added bd-balance.R for exactly this"
  ),

  mk(
    "m03",
    "src/bartcore/moves.hpp",
    paste0(
      "  } else if (!newIsCategorical) {\n",
      "    logProposalCorrection =\n",
      "      std::log(static_cast<double>(forwardValid)) -\n",
      "      std::log(static_cast<double>(forwardInterval));\n",
      "  } else if (!oldIsCategorical) {"
    ),
    paste0(
      "  } else if (!newIsCategorical) {\n",
      "    logProposalCorrection = 0.0;\n",
      "  } else if (!oldIsCategorical) {"
    ),
    "KILL_EXPECTED",
    kScript("benchmarks/R/change-balance.R"),
    paste0(
      "poison 3: change move's forward (new-side) proposal correction ",
      "dropped (was moves.hpp ~459); MIXED arm only visits this single-",
      "sided branch, not the both-ordinal one MAIN/CONTROL use - full mode ",
      "(quick's MCse leaves |z| 3.6, under threshold)"
    )
  ),

  mk(
    "m04",
    "src/bartcore/moves.hpp",
    paste0(
      "  if (!newIsCategorical && !oldIsCategorical) {\n",
      "    logProposalCorrection =\n",
      "      std::log(static_cast<double>(reverseInterval)) -\n",
      "      std::log(static_cast<double>(forwardInterval)) +\n",
      "      std::log(static_cast<double>(forwardValid)) -\n",
      "      std::log(static_cast<double>(reverseValid));"
    ),
    paste0(
      "  if (!newIsCategorical && !oldIsCategorical) {\n",
      "    logProposalCorrection =\n",
      "      std::log(static_cast<double>(forwardValid)) -\n",
      "      std::log(static_cast<double>(forwardInterval));"
    ),
    "KILL_EXPECTED",
    kScript("benchmarks/R/change-balance.R", "quick"),
    "poison 4: change move's reverse (old-side) proposal correction dropped (was moves.hpp ~474)"
  ),

  mk(
    "m05",
    "src/bartcore/moves.hpp",
    paste0(
      "  bool swapIsSensible = ruleIsValid(ctx, tree, parent, childRule.variableIndex);\n",
      "  if (childRule.variableIndex != parentRule.variableIndex && swapIsSensible)\n",
      "    swapIsSensible = ruleIsValid(ctx, tree, parent, parentRule.variableIndex);\n",
      "  // interaction is a WHOLE-subtree, all-variables property the per-variable\n",
      "  // ruleIsValid checks above cannot see (the swap sibling-strand break): a\n",
      "  // swap that lifts x2 above x3 co-occurs a forbidden pair with neither\n",
      "  // swapped variable equal to x3. Score it the -1.0 no-op (pi(T') = 0).\n",
      "  if (swapIsSensible) swapIsSensible = tree.interactionSubtreeIsValid(parent);"
    ),
    "  bool swapIsSensible = true;",
    "KILL_EXPECTED",
    kCpp(),
    "poison 5: swap move's descendant-validity walk skipped outright (was moves.hpp ~634; July's cpp arm caught this as a crash)"
  ),

  mk(
    "m06",
    "src/bartcore/model.hpp",
    "    double posteriorDegreesOfFreedom =\n      degreesOfFreedom + static_cast<double>(numPositiveWeights);",
    "    double posteriorDegreesOfFreedom =\n      degreesOfFreedom + static_cast<double>(numObservations);",
    "KILL_EXPECTED",
    kEquiv("zeroweights"),
    "poison 6: sigma posterior df counts zero-weight rows too (was model.hpp:1756)"
  ),

  mk(
    "m07",
    "src/bartcore/model.hpp",
    paste0(
      "    double sumOfSquaredResiduals = weights == nullptr\n",
      "      ? misc_computeSumOfSquaredResiduals(y, numObservations, totalFits)\n",
      "      : misc_computeWeightedSumOfSquaredResiduals(y, numObservations, weights,\n",
      "                                                  totalFits);"
    ),
    "    double sumOfSquaredResiduals =\n      misc_computeSumOfSquaredResiduals(y, numObservations, totalFits);",
    "KILL_EXPECTED",
    kEquiv("weighted"),
    "poison 7: sigma posterior SSR drops the per-row weight (was model.hpp:1751)"
  ),

  mk(
    "m08",
    "src/bartcore/model.hpp",
    "    double shape = 0.5 * (totalNumLeaves + degreesOfFreedom);",
    "    double shape = 0.5 * (totalNumLeaves + 2.0 * degreesOfFreedom - 1.0);",
    "KILL_EXPECTED",
    kEquiv("chik2"),
    "poison 8: chi-k hyperprior shape mislabeled (was model.hpp:1729); gate-hardening-1.0 added the chik2 scenario and its disjoint-seed-range channel for exactly this"
  ),

  mk(
    "m09",
    "src/bartcore/model.hpp",
    "    double rate = 0.5 * sumSquaredParams / (leafScale * leafScale);",
    "    double rate = 0.5 * sumSquaredParams;",
    "KILL_EXPECTED",
    kEquiv("chik"),
    "poison 9: chi-k hyperprior rate drops the /leafScale^2 term (was model.hpp:1731)"
  ),

  mk(
    "m10",
    "src/bartcore/model.hpp",
    "      double draw = ext_rng_simulateGamma(\n        rng, alpha / p + static_cast<double>(splitCounts[j]), 1.0);",
    "      double draw = ext_rng_simulateGamma(\n        rng, alpha / p + static_cast<double>(splitCounts[j]) + 1.0, 1.0);",
    "KILL_EXPECTED",
    kEquiv("dart"),
    "poison 10: DART's per-variable Dirichlet shape adds a spurious +1 to the split count (was model.hpp:1680)"
  ),

  mk(
    "m11",
    "src/bartcore/model.hpp",
    "        omega += ext_rng_simulatePolyaGamma(rng, psi);\n      omega_[i] = omega;\n      double weight = weights_ != nullptr ? weights_[i] : 1.0;",
    "        omega += ext_rng_simulatePolyaGamma(rng, psi);\n      omega_[i] = omega * omega;\n      double weight = weights_ != nullptr ? weights_[i] : 1.0;",
    "KILL_EXPECTED",
    kEquiv("logistic"),
    "poison 11: logistic response reports omega^2 as its working weight, though the working response itself still divides by the true omega (was model.hpp:2180)"
  ),

  # m12 retired, id left unassigned rather than reused: it targeted the
  # grouped-intercept precision accumulation in GroupedResponse, which is
  # deleted. No surviving site carries the same defect.

  mk(
    "m13",
    "src/bartcore/combiner.hpp",
    "    double priorPrecision = 1.0 / glue_.prior[f].variance;",
    "    double priorPrecision = 0.0;",
    "KILL_EXPECTED",
    kScript("benchmarks/R/bcf-exact-weak.R"),
    "poison 13: the amplitude full conditional drops its prior precision (was chain.hpp:2049, then combiner.hpp's two-scalar a-glue draw; retargeted to drawForestAmplitude's prior seed when that specialized path was deleted, which is the same defect at every forest rather than at a alone); gate-hardening-1.0 sub-item 1 added bcf-exact-weak.R for exactly this"
  ),

  mk(
    "m14",
    "src/bartcore/combiner.hpp",
    "          crossproduct[j * q + k] += wi * row[j] * row[k] * invSigmaSq;",
    "          crossproduct[j * q + k] += wi * row[j] * invSigmaSq;",
    "KILL_EXPECTED",
    kScript("benchmarks/R/bcf-exact.R", "quick"),
    "poison 14: the amplitude precision accumulates w*x instead of w*x^2 (was chain.hpp:2022, then combiner.hpp's two-scalar a-glue draw; retargeted to drawForestAmplitude's crossproduct when that specialized path was deleted)"
  ),

  mk(
    "m15",
    "src/bartcore/model.hpp",
    "                             projection);\n\n    double ridge = (k / scale) * (k / scale) * residualVariance;",
    "                             projection);\n\n    double ridge = (k / scale) * (k / scale);",
    "KILL_EXPECTED",
    kScript("benchmarks/R/linear-exact.R"),
    "poison 15: linear leaf's branch marginal (score side only) drops sigma^2 from the ridge (was model.hpp:304); gate-hardening-1.0 sub-item 4 added linear-exact.R for exactly this"
  ),

  mk(
    "m16",
    "src/bartcore/model.hpp",
    "      double w = weights == nullptr ? 1.0 : weights[i];\n      double noise = residualVariance / w;",
    "      double w = weights == nullptr ? 1.0 : weights[i];\n      double noise = residualVariance;",
    "KILL_EXPECTED",
    kEquiv("wtgp"),
    "poison 16: GP leaf's weighted score-path nugget drops /w_i (was model.hpp:721, the zero-weight fallback path stays clean); gate-hardening-1.0 sub-item 2 added the wtgp scenario for exactly this"
  ),

  mk(
    "m17",
    "R/validateComposition.R",
    "  agrees <- is.null(reference) ||\n    (length(value) == length(reference) &&\n      identical(names(value), names(reference)))",
    "  agrees <- is.null(reference) ||\n    length(value) == length(reference)",
    "KILL_EXPECTED",
    kTinytest("inst/tinytest/test-validate-composition.R"),
    "R/: compositionFunctionals stops checking that a renamed functional is still the ranked one, only its length"
  ),

  mk(
    "m18",
    "R/validateComposition.R",
    "  whole <- is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)\n  if (!whole || x < minimum) {",
    "  whole <- is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)\n  if (!whole) {",
    "KILL_EXPECTED",
    kTinytest("inst/tinytest/test-validate-composition.R"),
    "R/: compositionCount stops enforcing its per-argument minimum (n.replications >= 2, n.thin >= 1, n.burn >= 0)"
  ),

  # m19 retired, id left unassigned rather than reused: it targeted
  # rejectUnknownDotsArgs, which refused a typo'd or retired bart2/
  # rbart_vi argument name by checking it against the family's known
  # formals. Both that function and bart2/rbart_vi's own '...' formal are
  # gone; an unrecognized name now hits R's own base "unused argument"
  # error, which is interpreter behavior, not source this harness can
  # plant a mutation in. No other site is checked by the same gate, so
  # there is nothing live to repoint m19 at.

  mk(
    "m20",
    "R/utility.R",
    "checkMissingPolicy <- function(data, hasMissing, what) {\n  if (data@missing == \"error\" && hasMissing) {",
    "checkMissingPolicy <- function(data, hasMissing, what) {\n  if (data@missing != \"error\" && hasMissing) {",
    "KILL_EXPECTED",
    kTinytest("inst/tinytest/test-data-missing.R"),
    paste0(
      "R/: checkMissingPolicy's policy comparison inverted, so a ",
      "missing = \"error\" sampler no longer refuses new missing values ",
      "(dropping hasMissing instead is undetectable: the one existing test's",
      " scenario always has hasMissing TRUE)"
    )
  ),

  ## SURVIVE_DOCUMENTED: the setState / readWarmStartState SEXP-parsing
  ## validation clusters are the two largest UNTESTED-PATH clusters in the
  ## whole scoped codebase per p7-branch-reach.md (128 and 86 never-hit lines
  ## and branch arms respectively) - "a long sequence of per-block malformed-
  ## input checks... tests/cpp exercises the well-formed round trip and a
  ## handful of hand-picked malformations", none of which are these three.
  ## Their closest existing gates are tests/cpp (which never touches the R
  ## SEXP parser at all - it drives bartcore::SamplerBase::setState directly
  ## in C++) and the state-round-trip tinytest files (which never construct a
  ## state THIS malformed). Expected verdict: SURVIVED by every killer named -
  ## this executably documents the gap rather than asserting one that isn't
  ## there, and becomes a P5b/future test target.

  mk(
    "m21",
    "src/R_interface_bartcore.cpp",
    paste0(
      "    if (!Rf_isInteger(sampleNumExpr) || Rf_xlength(sampleNumExpr) != 1 ||\n",
      "        INTEGER(sampleNumExpr)[0] < 0)\n",
      "      errorMessage = \"malformed sample number in bartcore state\";"
    ),
    paste0(
      "    if (!Rf_isInteger(sampleNumExpr) || Rf_xlength(sampleNumExpr) != 1)\n",
      "      errorMessage = \"malformed sample number in bartcore state\";"
    ),
    "SURVIVE_DOCUMENTED",
    c(kCpp(), kTinytest("inst/tinytest/test-sampler-state-format.R")),
    "P7 setState cluster (R_interface_bartcore.cpp:6398): a negative currentSampleNum is no longer refused and wraps to a huge size_t on restore"
  ),

  mk(
    "m22",
    "src/R_interface_bartcore.cpp",
    paste0(
      "      if (!Rf_isReal(dartProbabilitiesExpr) || !Rf_isReal(dartAlphaExpr) ||\n",
      "          Rf_xlength(dartAlphaExpr) != 1 || !Rf_isInteger(dartSkippedExpr) ||\n",
      "          Rf_xlength(dartSkippedExpr) != 1 || INTEGER(dartSkippedExpr)[0] < 0) {\n",
      "        errorMessage = \"malformed dart state in bartcore state\";"
    ),
    paste0(
      "      if (!Rf_isReal(dartProbabilitiesExpr) || !Rf_isReal(dartAlphaExpr) ||\n",
      "          Rf_xlength(dartAlphaExpr) != 1 || !Rf_isInteger(dartSkippedExpr) ||\n",
      "          Rf_xlength(dartSkippedExpr) != 1) {\n",
      "        errorMessage = \"malformed dart state in bartcore state\";"
    ),
    "SURVIVE_DOCUMENTED",
    c(kCpp(), kTinytest("inst/tinytest/test-sampler-state-format.R")),
    "P7 setState cluster (R_interface_bartcore.cpp:6398): a negative dart.updates.skipped is no longer refused"
  ),

  mk(
    "m23",
    "src/R_interface_bartcore.cpp",
    paste0(
      "      SEXP kExpr = rc_getListElement(forestExpr, \"k\");\n",
      "      if (!Rf_isReal(kExpr) || Rf_xlength(kExpr) != 1) {\n",
      "        errorMessage = \"malformed parameters in warm-start donor\";\n",
      "        break;\n",
      "      }\n",
      "      fs.k = REAL(kExpr)[0];"
    ),
    paste0(
      "      SEXP kExpr = rc_getListElement(forestExpr, \"k\");\n",
      "      if (!Rf_isReal(kExpr)) {\n",
      "        errorMessage = \"malformed parameters in warm-start donor\";\n",
      "        break;\n",
      "      }\n",
      "      fs.k = REAL(kExpr)[0];"
    ),
    "SURVIVE_DOCUMENTED",
    c(kCpp(), kTinytest("inst/tinytest/test-warm-start.R")),
    "P7 readWarmStartState cluster (R_interface_bartcore.cpp:6709): a wrong-length k in a warm-start donor forest is silently truncated to its first element instead of refused"
  ),

  ## m24-m25 are the perturb kernel's two poisons: the whole Hastings term of a
  ## cut displacement is the window ratio, which fires only at the ends of the
  ## descendant-valid interval, so both defects are invisible in the middle of
  ## it and only the cut-law gate reads them.

  mk(
    "m24",
    "src/bartcore/moves.hpp",
    paste0(
      "  int32_t reverseCount = std::min(upper, target + perturbWidth) -\n",
      "                         std::max(lower, target - perturbWidth);\n",
      "  double logProposalCorrection =\n",
      "    std::log(static_cast<double>(forwardCount)) -\n",
      "    std::log(static_cast<double>(reverseCount));"
    ),
    "  double logProposalCorrection = 0.0;",
    "KILL_EXPECTED",
    kScript("benchmarks/R/perturb-balance.R"),
    "poison 24: perturb move's window proposal correction dropped; the uncorrected chain is reversible for pi(c)|W(c)| and the boundary cuts starve"
  ),

  mk(
    "m25",
    "src/bartcore/moves.hpp",
    "  int32_t forwardLow = std::max(lower, current - perturbWidth);",
    "  int32_t forwardLow = current;",
    "KILL_EXPECTED",
    kScript("benchmarks/R/perturb-balance.R"),
    "poison 25: perturb move's window made one-sided (+w only), so every proposal is c -> c+1 and the cut is absorbed at the top of the interval"
  ),

  mk(
    "m26",
    "src/bartcore/moves.hpp",
    paste0(
      "        double below =\n",
      "          std::log(1.0 -\n",
      "                   ctx.treePrior.growthProbability(tree, data, leftChild)) +\n",
      "          std::log(1.0 - ctx.treePrior.growthProbability(tree, data,\n",
      "                                                         leftChild + 1));"
    ),
    "        double below = 0.0;",
    "KILL_EXPECTED",
    kScript("benchmarks/R/rule-gibbs-balance.R"),
    "poison 26: rule_gibbs neighbourhood weights lose the two log(1 - growth(child)) terms, so a candidate that strands a grandchild is no longer favoured"
  ),

  mk(
    "m27",
    "src/bartcore/moves.hpp",
    paste0(
      "    double logRulePrior = -std::log(static_cast<double>(high - low + 1)) -\n",
      "                          (doubled ? std::log(2.0) : 0.0);"
    ),
    "    double logRulePrior = (doubled ? -std::log(2.0) : 0.0);",
    "KILL_EXPECTED",
    kScript("benchmarks/R/rule-gibbs-balance.R"),
    "poison 27: rule_gibbs neighbourhood weights lose the 1/|SI_v| rule factor, the low-cardinality bias change-balance.R's own gate repaired"
  )
)

## ---- third review: code changed since 7ad0bbea ----------------------------
## m28-m89 are one-token defects in what landed after the second review: the
## monotone engine, the multi-forest leaf-prior writer, the survival, hurdle,
## zero-trial, nbinom, t and ordinal families, predictor-mutation rollback,
## data ingestion, the bridge's validation, and the generics. Killers are the
## tinytest files, tests/cpp and quick exact gates that own each site.
## Of the thirteen that first survived the full tinytest suite, tests/cpp and
## the cheap quick gates, twelve (m38, m42, m59, m63, m67, m71, m74, m76, m78,
## m79, m84, m87) each have a killing test now and are KILL_EXPECTED; r3s stays
## for the next gap found. m80 (mapFactorColumnsToTrainingLevels coding a test
## factor against its own levels) is equivalent: every validateXTest path then
## builds a container that alignContainerFactorLevels re-codes by label, and
## the unseen-level refusal it would skip raises there with the same message.
## Ids left
## unassigned: m36 (monotoneMovePair flip) and m55 (aft setSurvivalStatus
## restoring logT_) are equivalent on every reachable path - a split on a
## constrained axis keeps its children in one component, where the ratio is
## symmetric, and setResponse overwrites logT_ right after the status
## install; m40 (redrawAfterBirth's pair floor) rides on drawPairUpper's
## inversion, which 8e7a3d19 rewrote.

kTests <- function(...) {
  do.call(
    c,
    lapply(c(...), function(f) kTinytest(file.path("inst/tinytest", f)))
  )
}
kGate <- function(name) kScript(file.path("benchmarks/R", name), "quick")
r3 <- function(id, file, text, mutant, killers, note) {
  mk(id, file, text, mutant, "KILL_EXPECTED", killers, note)
}
r3s <- function(id, file, text, mutant, killers, note) {
  mk(id, file, text, mutant, "SURVIVE_DOCUMENTED", killers, note)
}
kMonotone <- c(
  kTests("test-monotone.R"),
  kCpp(),
  kGate("monotone-reference.R")
)
mh <- "src/bartcore/model.hpp"

mutations <- c(
  mutations,
  list(
    r3(
      "m28",
      mh,
      "bool rightAbove = hi[il] + 1 == lo[ir];",
      "bool rightAbove = hi[il] == lo[ir];",
      kMonotone,
      "monotone geometry: the adjacency test along an ordered axis loses its +1"
    ),
    r3(
      "m29",
      mh,
      "return std::max(lo[il], lo[ir]) <= std::min(hi[il], hi[ir]) ||",
      "return std::max(lo[il], lo[ir]) < std::min(hi[il], hi[ir]) ||",
      kMonotone,
      "monotone geometry: two leaves sharing one cut cell no longer share the axis"
    ),
    r3(
      "m30",
      mh,
      "if (tree.ruleMissingGoesRight(data, rule) != isRight) leafMissing[a] = 0;",
      "if (tree.ruleMissingGoesRight(data, rule) == isRight) leafMissing[a] = 0;",
      c(kMonotone, kTests("test-data-missing.R")),
      "monotone geometry on a missing axis: the reaches-missing flag clears on the wrong side"
    ),
    r3(
      "m31",
      mh,
      "if ((wl[w] & wr[w]) != 0) return true;",
      "if ((wl[w] & wr[w]) == 0) return true;",
      c(kMonotone, kTests("test-data-categorical.R")),
      "monotone geometry on a factor axis: level-set overlap inverted"
    ),
    r3(
      "m32",
      mh,
      "next.forward[monotoneLayerInsert(next, s.key.data(), words)] += count;",
      "next.forward[monotoneLayerInsert(next, s.key.data(), words)] = count;",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone order count: the down-set DP overwrites rather than sums its forward counts"
    ),
    r3(
      "m33",
      mh,
      "up.backwardExponent) * std::numbers::ln2 -",
      "0) * std::numbers::ln2 -",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone position law drops the backward counts' binary exponent (matters once a layer rescales)"
    ),
    r3(
      "m34",
      mh,
      "return firstLaw[i - 1] + secondLaw[j - 1] + logChoose(i + j - 2, i - 1) +",
      "return firstLaw[i - 1] + secondLaw[j - 1] + logChoose(i + j - 1, i - 1) +",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone two-component ratio: interleaving binomial off by one"
    ),
    r3(
      "m35",
      mh,
      "result = std::log(static_cast<double>(first.size + second.size)) +",
      "result = std::log(static_cast<double>(first.size)) +",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone two-component ratio: m loses the second component's size"
    ),
    r3(
      "m37",
      "src/bartcore/moves.hpp",
      "*ratio *= std::exp(birth ? logRatio : -logRatio);",
      "*ratio *= std::exp(birth ? -logRatio : logRatio);",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone leaf prior: Z_T ratio applied with the wrong sign on both moves"
    ),
    r3(
      "m38",
      mh,
      "drawn[next++] = freeSd * ext_rng_simulateStandardNormal(rng);",
      "drawn[next++] = constrainedSd * ext_rng_simulateStandardNormal(rng);",
      kMonotone,
      "monotone prior leaf draw: an isolated leaf takes the c-inflated sd; killed by test_monotone.cpp's exact prior draw, whose leaf variances now match rejection sampling"
    ),
    r3(
      "m39",
      mh,
      "return (constrained ? cInflation : 1.0) * scale / k;",
      "return (constrained ? 1.0 : cInflation) * scale / k;",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone prior sd: c-inflation applied to the free leaves instead of the constrained"
    ),
    r3(
      "m41",
      mh,
      "mu[lower] = drawTruncatedNormal(rng, mL, sL, aL, std::min(bL, muUpper));",
      "mu[lower] = drawTruncatedNormal(rng, mL, sL, aL, bL);",
      kMonotone,
      "monotone birth redraw: the lower child is not capped by its drawn sibling"
    ),
    r3(
      "m42",
      mh,
      "if (lower > upper) return false;",
      "if (lower < upper) return false;",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone joint prior: the structure draw's cone test inverted; equivalent in the tree law (negated iid draws of one sd map the cone onto its reverse), killed by test_monotone.cpp's accept-is-cone-test check against the point oracle's pairs"
    ),
    r3(
      "m43",
      "src/bartcore/chain.hpp",
      "forest.leaf.prior = options.monotonePrior == 1 ? MonotonePrior::joint",
      "forest.leaf.prior = options.monotonePrior == 0 ? MonotonePrior::joint",
      c(kMonotone, kGate("monotone-successive-conditional.R")),
      "monotone prior switch: joint and leaf exchanged at creation"
    ),
    r3(
      "m44",
      "src/bartcore/chain.hpp",
      "if (forest.leaf.prior != MonotonePrior::joint ||",
      "if (forest.leaf.prior == MonotonePrior::joint ||",
      c(kMonotone, kTests("test-calibration-prior-draws.R")),
      "monotone prior tree draw: the joint prior's accept step runs under leaf instead"
    ),
    r3(
      "m45",
      "src/bartcore/chain.hpp",
      paste0(
        "    nodeScaleFactors_[f] = sd;\n",
        "    forests_[f].leaf.scale = mapLeafScale(f);\n"
      ),
      "    nodeScaleFactors_[f] = sd;\n",
      kTests(
        "test-multiforest-leaf-prior-writer.R",
        "test-calibration-midchain.R",
        "test-forest-basis-r5.R"
      ),
      "leaf-prior writer: forest(sd = ) on a fixed-variance map forest records the factor but never re-derives the leaf scale"
    ),
    r3(
      "m46",
      "R/dbarts.R",
      ".Call(C_dbarts_bartcore_setForestSd, ptr, index - 1L, sds[[index]])",
      ".Call(C_dbarts_bartcore_setForestSd, ptr, index, sds[[index]])",
      kTests(
        "test-multiforest-leaf-prior-writer.R",
        "test-calibration-midchain.R",
        "test-forest-basis-r5.R"
      ),
      "leaf-prior writer: forest index passed 1-based to the engine"
    ),
    r3(
      "m47",
      "R/dbarts.R",
      "params[[if (params[[7L]] > 0) 7L else 4L]] <- sds[[index]]",
      "params[[if (params[[7L]] < 0) 7L else 4L]] <- sds[[index]]",
      kTests(
        "test-multiforest-leaf-prior-writer.R",
        "test-calibration-midchain.R",
        "test-forest-basis-r5.R"
      ),
      "leaf-prior writer: the control mirror writes the wrong channel, so a re-creation restates the creation sd"
    ),
    r3(
      "m48",
      "src/bartcore/chain.hpp",
      "if (k != forests_[f].k) forests_[f].k = k;",
      "if (k == forests_[f].k) forests_[f].k = k;",
      kTests(
        "test-multiforest-leaf-prior-writer.R",
        "test-multinomial-r5-surface.R"
      ),
      "leaf-prior writer: a multinomial normal(k = ) restatement never writes k"
    ),
    r3(
      "m49",
      "R/dbarts.R",
      "findInterval(time, periods, left.open = TRUE) + 1L,",
      "findInterval(time, periods, left.open = FALSE) + 1L,",
      c(
        kTests("test-hazard.R", "test-hazard-factors.R"),
        kGate("hazard-exact.R"),
        kGate("hazard-reduction.R")
      ),
      "hazard grid: a time on a grid point lands in the next period"
    ),
    r3(
      "m50",
      "R/dbarts.R",
      "result$offset <- offset[subjectOf]",
      "result$offset <- offset[periodOf]",
      c(
        kTests(
          "test-hazard.R",
          "test-hazard-factors.R",
          "test-family-offset.R",
          "test-family-mutation-parity.R"
        ),
        kGate("hazard-reduction.R")
      ),
      "hazard expansion: offsets replicated by period instead of subject"
    ),
    r3(
      "m51",
      "R/bart.R",
      "m <- sum(periods <= times[j])",
      "m <- sum(periods < times[j])",
      c(
        kTests(
          "test-hazard-grid-horizon.R",
          "test-hazard.R",
          "test-hazard-factors.R",
          "test-predict-na-action.R"
        ),
        kGate("hazard-exact.R")
      ),
      "hazard survivalProbabilities: a horizon on a grid point excludes its own period; test-hazard-grid-horizon.R pins S(t) on and off the grid against the (1 - hazard) product"
    ),
    r3(
      "m52",
      "R/bart.R",
      "bigX[[\"period\"]] <- rep(seq_len(K), each = n)",
      "bigX[[\"period\"]] <- rep(seq_len(K), times = n)",
      kTests(
        "test-hazard.R",
        "test-hazard-factors.R",
        "test-predict-na-action.R"
      ),
      "hazard survivalProbabilities on a data.frame newdata: period column misaligned with the period-major rows"
    ),
    r3(
      "m53",
      "R/generics.R",
      "ev = piVec * exp(fVec + 0.5 * sigmaVec^2),",
      "ev = piVec * exp(fVec + sigmaVec^2),",
      c(
        kTests("test-hurdle.R", "test-hurdle-surface.R"),
        kGate("hurdle-exact.R")
      ),
      "hurdle ev: lognormal mean loses the 1/2 on sigma^2"
    ),
    r3(
      "m54",
      "R/generics.R",
      "    ) -\n    log(yRep[positive])",
      "    ) +\n    log(yRep[positive])",
      kTests(
        "test-hurdle.R",
        "test-hurdle-surface.R",
        "test-pointwise-loglik.R"
      ),
      "hurdle loglik: lognormal Jacobian added instead of subtracted"
    ),
    r3(
      "m56",
      mh,
      "variance_ != nullptr ? std::sqrt(variance_[i]) * varianceScale : sd;",
      "variance_ != nullptr ? variance_[i] * varianceScale : sd;",
      c(kTests("test-aft-heteroscedastic.R"), kGate("aft-exact.R")),
      "heteroscedastic aft: censored redraw uses the variance as the sd"
    ),
    r3(
      "m57",
      "src/bartcore/combiner.hpp",
      "for (int c = 1; c < trials_[i]; ++c)",
      "for (int c = 1; c <= trials_[i]; ++c)",
      c(
        kTests(
          "test-multinomial-counts-mutation.R",
          "test-multinomial-surface.R"
        ),
        kGate("multinomial-exact.R")
      ),
      "multinomial counts: one extra Polya-Gamma draw per row"
    ),
    r3(
      "m58",
      "src/bartcore/combiner.hpp",
      "    while (i < n && trials_[i] != 0) ++i;\n    if (i == n) return;",
      "    while (i < n && trials_[i] != 0) ++i;\n    if (i != n) return;",
      kTests(
        "test-multinomial-zero-trials.R",
        "test-multinomial-counts-mutation.R"
      ),
      "multinomial zero-trial rows never composed into the effective mask"
    ),
    r3(
      "m59",
      "R/data.R",
      "codes <- if (is.factor(y)) as.integer(y) else as.integer(y) + 1L",
      "codes <- if (is.factor(y)) as.integer(y) else as.integer(y)",
      kTests(
        "test-multinomial-numeric-codes.R",
        "test-multinomial-surface.R",
        "test-multinomial-generics.R"
      ),
      "multinomial numeric category codes not shifted to 1-based"
    ),
    r3(
      "m60",
      "R/generics.R",
      "logCoef <- lgamma(n + 1) - rowSums(lgamma(counts + 1))",
      "logCoef <- lgamma(n + 1)",
      kTests(
        "test-multinomial-generics.R",
        "test-multinomial-zero-trials.R",
        "test-pointwise-loglik.R"
      ),
      "multinomial loglik: count-matrix multinomial coefficient loses its denominator"
    ),
    r3(
      "m61",
      mh,
      paste0(
        "      if (active != nullptr && active[i] == 0.0) continue;\n",
        "      histogram["
      ),
      "      histogram[",
      c(kTests("test-nbinom.R", "test-active-rows-pins.R"), kCpp()),
      "nbinom shape kernel: count histogram ignores the active-row mask"
    ),
    r3(
      "m62",
      mh,
      "rPrior_.computeKernel(y_, numObservations_, activePointer());",
      "rPrior_.computeKernel(y_, numObservations_);",
      c(kTests("test-nbinom.R", "test-active-rows-pins.R"), kCpp()),
      "nbinom setActiveRows: kernel not rebuilt over the subsample"
    ),
    r3(
      "m63",
      mh,
      "y * logOnePlusExp(-psi) - r_ * logOnePlusExp(psi);",
      "y * logOnePlusExp(psi) - r_ * logOnePlusExp(psi);",
      c(kTests("test-nbinom.R", "test-pointwise-loglik.R"), kCpp()),
      "nbinom pointwise loglik: log p sign flipped (re-pointed to the log-mean form, psi = log mu - log r); killed by test_model.cpp's nb log-mean anchor, which checks it against dnbinom at mu"
    ),
    r3(
      "m64",
      "R/generics.R",
      "    rep(y, each = n.draws),\n    size = disp,",
      "    rep(y, times = n.draws),\n    size = disp,",
      kTests("test-nbinom.R", "test-pointwise-loglik.R"),
      "nbinom R loglik: response replicated in the wrong layout"
    ),
    r3(
      "m65",
      mh,
      "ext_rng_simulateGamma(rng, shape, 2.0 / (nu_ + w * r * r / sigmaSq));",
      "ext_rng_simulateGamma(rng, shape, 2.0 / (nu_ + r * r / sigmaSq));",
      c(kTests("test-bart-weights-parity.R"), kGate("t-exact.R"), kCpp()),
      "student t: lambda draw ignores the user weight"
    ),
    r3(
      "m66",
      mh,
      "if (w * a > 0.0) {",
      "if (w > 0.0) {",
      c(kTests("test-active-rows-pins.R"), kGate("t-exact.R"), kCpp()),
      "student t: nu statistics count masked rows"
    ),
    r3(
      "m67",
      mh,
      "? sigmaOriginal / std::sqrt(userWeights_[i])",
      "? sigmaOriginal / userWeights_[i]",
      c(
        kTests("test-pointwise-loglik.R", "test-bart-weights-parity.R"),
        kCpp()
      ),
      "student t loglik: weighted scale divides by w rather than sqrt(w); killed by test_model.cpp's weighted t log-likelihood against the closed-form dt"
    ),
    r3(
      "m68",
      mh,
      "if (s < numCategories_ - 1)  // finite upper gap only below the top cutpoint",
      "if (s <= numCategories_ - 1)  // finite upper gap only below the top cutpoint",
      c(kTests("test-ordinal.R"), kGate("ordinal-exact.R"), kCpp()),
      "ordinal cutpoint MH: top cutpoint scored against an infinite upper gap"
    ),
    r3(
      "m69",
      mh,
      "if (!isActive(i)) continue;  // the target is the subsample's likelihood",
      "// the target is the subsample's likelihood",
      c(kTests("test-ordinal.R", "test-active-rows-pins.R"), kCpp()),
      "ordinal cutpoint MH: masked rows enter the acceptance likelihood"
    ),
    r3(
      "m70",
      "R/generics.R",
      "  idx <- rep(seq_len(nObs), each = n.draws)\n  result <- log(",
      "  idx <- rep(seq_len(nObs), times = n.draws)\n  result <- log(",
      kTests("test-ordinal.R", "test-pointwise-loglik.R"),
      "ordinal R loglik: observation index replicated in the wrong layout"
    ),
    r3(
      "m71",
      "src/bartcore/sampler.hpp",
      "      data_.gatheredRawValues = std::move(oldGatheredRaw);\n",
      "",
      c(
        kTests(
          "test-linear-leaves.R",
          "test-gp-leaves.R",
          "test-composition-sequences.R"
        ),
        kCpp()
      ),
      "predictor rollback: leaf-covariate raw copies not restored; killed by test-linear-leaves.R's refused replacement against a twin that never attempted it"
    ),
    r3(
      "m72",
      "src/bartcore/sampler.hpp",
      "      for (auto& chain : chains_) chain->repartitionTrees();\n",
      "",
      c(
        kTests("test-sampler-predictors.R", "test-data-mixed-mutation.R"),
        kCpp()
      ),
      "predictor rollback: trees not repartitioned after the restore"
    ),
    r3(
      "m73",
      "src/bartcore/sampler.hpp",
      "if ((updateCutPoints || data_.isFactor(j)) &&",
      "if (updateCutPoints &&",
      c(
        kTests("test-data-categorical.R", "test-data-categorical-declared.R"),
        kCpp()
      ),
      "predictor mutation: factor level-code precheck skipped without a cut refresh"
    ),
    r3(
      "m74",
      "src/bartcore/sampler.hpp",
      "        data.hasMissing[j] = oldHasMissing[k];\n",
      "",
      c(kTests("test-data-missing.R", "test-sampler-predictors.R"), kCpp()),
      "subset predictor rollback: missingness flags not restored; killed by test_sampler.cpp's subset rollback missingness"
    ),
    r3(
      "m75",
      "src/bartcore/data.hpp",
      "value < static_cast<double>(categoryCounts[variable]) &&",
      "value <= static_cast<double>(categoryCounts[variable]) &&",
      c(
        kTests("test-data-categorical.R", "test-data-categorical-declared.R"),
        kCpp()
      ),
      "factor ingestion: a code one past the level table is accepted"
    ),
    r3(
      "m76",
      "src/bartcore/data.hpp",
      "splitsBySubset(j) && source.slice.numNonzero < numObservations",
      "source.slice.numNonzero < numObservations",
      c(kTests("test-sparse-factor.R", "test-data-sparse.R"), kCpp()),
      "sparse ingestion: an ordered factor's reference folds into its level count; killed by test_data.cpp's undeclared CSC ordered factor (unreachable from R, which declares K or passes reference 0)"
    ),
    r3(
      "m77",
      "src/bartcore/data.hpp",
      paste0(
        "    return nzCodes[wordRanks[i >> 6] +\n",
        "                   static_cast<size_t>(std::popcount(word & (bit - 1u)))];"
      ),
      paste0(
        "    return nzCodes[wordRanks[i >> 6] +\n",
        "                   static_cast<size_t>(std::popcount(word & ((bit << 1) - 1u)))];"
      ),
      c(kTests("test-data-sparse.R", "test-sparse-factor.R"), kCpp()),
      "sparse rank storage: a row's nonzero rank counts its own bit"
    ),
    r3(
      "m78",
      "R/data.R",
      "rows[x@i[missingEntries] + 1L] <- TRUE",
      "rows[x@i[missingEntries]] <- TRUE",
      kTests(
        "test-na-action-sparse.R",
        "test-sparse-factor-na.R",
        "test-predict-na-action.R",
        "test-na-action.R"
      ),
      "sparse missing rows: 0-based row index used as 1-based"
    ),
    r3(
      "m79",
      "R/data.R",
      "rows <- rows | rowsWithMissingPredictors(x$sparse)",
      "rows <- rowsWithMissingPredictors(x$sparse)",
      kTests(
        "test-na-action-sparse.R",
        "test-sparse-factor-na.R",
        "test-predict-na-action.R",
        "test-na-action.R"
      ),
      "mixed container missing rows: dense-column NA rows dropped when a sparse block is present"
    ),
    r3(
      "m81",
      "R/utility.R",
      "referenceTaken <- length(column@i) < column@length",
      "referenceTaken <- length(column@i) <= column@length",
      kTests(
        "test-sparse-factor.R",
        "test-sparse-factor-frames.R",
        "test-predict-sparse.R",
        "test-sparse-factor-na.R"
      ),
      "sparseFactor test remap: a fully stored column still demands its reference be a training level"
    ),
    r3(
      "m82",
      "src/R_interface_bartcore.cpp",
      "y[i] != std::floor(y[i]) || y[i] < 1.0 ||",
      "y[i] != std::floor(y[i]) || y[i] < 0.0 ||",
      kTests("test-ordinal.R", "test-sampler-bridge-errors.R"),
      "bridge: ordinal response category 0 accepted"
    ),
    r3(
      "m83",
      "src/R_interface_bartcore.cpp",
      "if (!std::isfinite(y[i]) || y[i] < 0.0 || y[i] != std::floor(y[i]))",
      "if (!std::isfinite(y[i]) || y[i] < 0.0)",
      kTests("test-nbinom.R", "test-sampler-bridge-errors.R"),
      "bridge: nbinom fractional counts accepted"
    ),
    r3(
      "m84",
      "src/R_interface_bartcore.cpp",
      "if (active[i] != 0.0 && active[i] != 1.0)",
      "if (active[i] != 0.0 && active[i] > 1.0)",
      kTests("test-active-rows-pins.R", "test-sampler-bridge-errors.R"),
      "bridge: fractional active-row mask accepted (the R method refuses first, so the killer calls the bridge directly)"
    ),
    r3(
      "m85",
      "src/R_interface_bartcore.cpp",
      "      static_cast<size_t>(INTEGER(dimsExpr)[1]) != K)\n    Rf_error(\"%s: requires a real matrix",
      "      static_cast<size_t>(INTEGER(dimsExpr)[0]) != K)\n    Rf_error(\"%s: requires a real matrix",
      kTests(
        "test-multinomial-category-offset.R",
        "test-multinomial-test-offset.R"
      ),
      "bridge: category offset column count checked against the row dim"
    ),
    r3(
      "m86",
      "src/R_interface_bartcore.cpp",
      paste0(
        "  // capacity: a run short of capacity reports the draws it made\n",
        "  size_t numSamples = capacity > 0 ? shape.numSavedDraws : 1;"
      ),
      paste0(
        "  // capacity: a run short of capacity reports the draws it made\n",
        "  size_t numSamples = capacity > 0 ? capacity : 1;"
      ),
      kTests(
        "test-predict-forest.R",
        "test-tree-store-order.R",
        "test-bartcore-keepfits.R"
      ),
      "bridge predict: draw axis sized to capacity, not recorded draws"
    ),
    r3(
      "m87",
      "R/generics.R",
      "probs <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)",
      "probs <- c((1 - ci.level) / 2, 1 - (1 - ci.level))",
      kTests("test-generics-intervals.R", "test-hurdle.R", "test-ordinal.R"),
      "posteriorInterval: upper quantile loses the /2"
    ),
    r3(
      "m88",
      "R/generics.R",
      "if (combine) as.vector(t(x)) else x",
      "if (combine) as.vector(x) else x",
      kTests(
        "test-one-chain-dimension.R",
        "test-nbinom.R",
        "test-utility-chains.R",
        "test-convergence-diagnostics.R"
      ),
      "reshapeScalarChannel: a split scalar channel combined sample-major"
    ),
    r3(
      "m89",
      "R/data.R",
      "    codes <- as.integer(data@y) + 1L",
      "    codes <- as.integer(data@y)",
      kTests("test-ordinal.R"),
      "ordinal R ingestion: 0-based factor codes handed to the engine (control: the bridge refuses)"
    )
  )
)
names(mutations) <- vapply(mutations, `[[`, character(1), "id")

## ---- shared plumbing --------------------------------------------------

countOccurrences <- function(text, needle) {
  m <- gregexpr(needle, text, fixed = TRUE)[[1]]
  if (identical(as.vector(m), -1L)) 0L else length(m)
}
readWhole <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

verifyOne <- function(mut, root) {
  path <- file.path(root, mut$file)
  if (!file.exists(path)) {
    return(list(ok = FALSE, msg = "file missing"))
  }
  n <- countOccurrences(readWhole(path), mut$anchor)
  if (n == 1L) {
    list(ok = TRUE, msg = "OK")
  } else {
    list(ok = FALSE, msg = sprintf("%d hits (expected 1)", n))
  }
}

applyMutation <- function(mut, root) {
  path <- file.path(root, mut$file)
  content <- readWhole(path)
  n <- countOccurrences(content, mut$original)
  if (n != 1L) {
    stop(sprintf(
      "[%s] anchor %s in %s: expected 1 hit, found %d - drifted, refusing to mutate blind",
      mut$id,
      if (n == 0L) "absent" else "ambiguous",
      mut$file,
      n
    ))
  }
  writeLines(sub(mut$original, mut$mutant, content, fixed = TRUE), path)
}

archiveHeadInto <- function(dest) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  tarPath <- tempfile(fileext = ".tar")
  on.exit(unlink(tarPath))
  if (
    system2("git", c("-C", repoRoot, "archive", "HEAD", "-o", tarPath)) != 0L
  ) {
    stop("git archive HEAD failed")
  }
  if (system2("tar", c("-x", "-f", tarPath, "-C", dest)) != 0L) {
    stop("tar extraction failed")
  }
}

# Runs one non-killer build step (install, freshness guard); ok is TRUE iff
# it exited 0.
runStep <- function(cmd, args) {
  out <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  list(ok = (attr(out, "status") %||% 0L) == 0L, output = out)
}
installBuild <- function(root, lib) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  runStep("R", c("CMD", "INSTALL", "--preclean", "-l", lib, root))
}
checkFreshness <- function(root, lib) {
  runStep(
    "Rscript",
    c(file.path(root, "tools", "check-build-freshness.R"), lib, root)
  )
}

# One killer's verdict. `pass` is the KILLER SCRIPT's own exit status
# (0 = it saw nothing wrong, nonzero = it caught something) - the caller
# maps that onto CAUGHT/CLEAN against the entry's class. A shared, heavily-
# loaded dev machine occasionally fails to even fork/exec the subprocess
# ("error in running command", no exit status at all, not the killer script
# saying anything) - that is an infrastructure hiccup, not a verdict, so it
# gets a few retries before being reported as a genuine harness error.
runKillerOne <- function(killer, root, lib, retries = 3L) {
  oldwd <- setwd(file.path(root, killer$cwd %||% "."))
  on.exit(setwd(oldwd))
  env <- c(paste0("R_LIBS=", lib), killer$env)
  # system2() does not itself quote `args` before handing the assembled line
  # to the shell env forces it through - an unquoted arg containing shell
  # metacharacters (every -e expression here) corrupts the command line, so
  # every arg is quoted explicitly (the command itself must stay bare: a
  # quoted command is looked up literally, quote marks and all, and fails).
  run <- function() {
    tryCatch(
      list(
        errored = FALSE,
        out = suppressWarnings(system2(
          killer$argv[1],
          shQuote(killer$argv[-1]),
          env = env,
          stdout = TRUE,
          stderr = TRUE
        ))
      ),
      error = function(e) list(errored = TRUE, out = conditionMessage(e))
    )
  }
  result <- run()
  for (attempt in seq_len(retries - 1L)) {
    if (!result$errored) {
      break
    }
    Sys.sleep(5)
    result <- run()
  }
  out <- if (result$errored) {
    structure(paste("harness error:", result$out), status = 1L)
  } else {
    result$out
  }
  status <- attr(out, "status") %||% 0L
  list(pass = status == 0L, argv = killer$argv, output = out)
}

runEntry <- function(mut, root, lib) {
  killerResults <- lapply(mut$killers, runKillerOne, root = root, lib = lib)
  caught <- !vapply(killerResults, `[[`, logical(1), "pass")
  anyCaught <- any(caught)
  good <- if (mut$class == "KILL_EXPECTED") anyCaught else !anyCaught
  list(killerResults = killerResults, anyCaught = anyCaught, good = good)
}

## ---- modes -------------------------------------------------------------

modeList <- function(muts) {
  for (mut in muts) {
    cat(sprintf(
      "%-4s %-18s %-14s %s\n",
      mut$id,
      mut$class,
      basename(mut$file),
      mut$note
    ))
  }
  invisible(0L)
}

modeVerifyAnchors <- function(muts) {
  ok <- TRUE
  for (mut in muts) {
    v <- verifyOne(mut, repoRoot)
    cat(sprintf(
      "%-4s %-4s %s (%s)\n",
      mut$id,
      if (v$ok) "OK" else "FAIL",
      mut$file,
      v$msg
    ))
    ok <- ok && v$ok
  }
  if (ok) {
    cat("\nverify-anchors: all anchors resolve cleanly\n")
  } else {
    cat("\nverify-anchors: DRIFT DETECTED (see FAIL lines above)\n")
  }
  if (ok) 0L else 1L
}

modeRun <- function(muts, keepGoing) {
  scratchRoot <- Sys.getenv("MUTATION_BATTERY_SCRATCH_DIR", tempdir())
  lib <- file.path(scratchRoot, "mutation-battery-lib")
  overallOk <- TRUE
  for (mut in muts) {
    cat(sprintf("\n==== %s [%s] %s ====\n", mut$id, mut$class, mut$note))
    root <- tempfile(
      pattern = paste0("mutation-battery-", mut$id, "-"),
      tmpdir = scratchRoot
    )
    entryOk <- tryCatch(
      {
        archiveHeadInto(root)
        applyMutation(mut, root)
        cat(sprintf("[%s] mutated %s, installing...\n", mut$id, mut$file))
        inst <- installBuild(root, lib)
        if (!inst$ok) {
          cat(paste(inst$output, collapse = "\n"), "\n")
          stop(sprintf("[%s] R CMD INSTALL --preclean failed", mut$id))
        }
        fresh <- checkFreshness(root, lib)
        if (!fresh$ok) {
          cat(paste(fresh$output, collapse = "\n"), "\n")
          stop(sprintf(
            "[%s] build-freshness guard failed - stale install",
            mut$id
          ))
        }
        result <- runEntry(mut, root, lib)
        for (kr in result$killerResults) {
          how <- if (kr$pass) "clean (exit 0)" else "CAUGHT (nonzero exit)"
          cat(sprintf(
            "[%s] killer %s: %s\n",
            mut$id,
            paste(kr$argv, collapse = " "),
            how
          ))
          cat(paste(" ", tail(kr$output, 15L)), sep = "\n")
        }
        verdict <- if (result$anyCaught) "KILLED" else "SURVIVED"
        tag <- if (result$good) {
          "(as expected)"
        } else {
          "<- UNEXPECTED, this is news"
        }
        cat(sprintf("[%s] %s: %s %s\n", mut$id, mut$class, verdict, tag))
        result$good
      },
      error = function(e) {
        cat(sprintf("[%s] ERROR: %s\n", mut$id, conditionMessage(e)))
        FALSE
      }
    )
    unlink(root, recursive = TRUE, force = TRUE)
    overallOk <- overallOk && entryOk
    if (!entryOk && !keepGoing) {
      cat(sprintf("\nstopping at %s (no --keep-going)\n", mut$id))
      break
    }
  }
  cat(sprintf("\nmutation-battery: %s\n", if (overallOk) "PASS" else "FAIL"))
  if (overallOk) 0L else 1L
}

## ---- entry point ---------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
keepGoing <- "--keep-going" %in% args
args <- setdiff(args, "--keep-going")
mode <- if (length(args) >= 1L) args[[1L]] else ""

status <- if (mode == "list") {
  modeList(mutations)
} else if (mode == "verify-anchors") {
  modeVerifyAnchors(mutations)
} else if (mode == "run") {
  selector <- if (length(args) >= 2L) args[[2L]] else "all"
  ids <- if (identical(selector, "all")) {
    names(mutations)
  } else {
    strsplit(selector, ",", fixed = TRUE)[[1L]]
  }
  unknown <- setdiff(ids, names(mutations))
  if (length(unknown) > 0L) {
    stop("unknown mutation id(s): ", paste(unknown, collapse = ", "))
  }
  modeRun(mutations[ids], keepGoing)
} else {
  cat(
    "usage: mutation-battery.R list | verify-anchors | run all|id[,id...] [--keep-going]\n"
  )
  2L
}
quit(status = as.integer(status))
