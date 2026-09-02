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
equivBaseline <- "benchmarks/baselines/equivalence-02d41365.rds"

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
## SURVIVE_DOCUMENTED trio.

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

  mk(
    "m12",
    "src/bartcore/model.hpp",
    "    weightScratch[j] += w;",
    "    weightScratch[j] += 1.0;",
    "KILL_EXPECTED",
    kEquiv("grouped"),
    "poison 12: grouped-intercept precision counts members, not their weight (was model.hpp:2334); gate-hardening-1.0 added the grouped scenario for exactly this"
  ),

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
