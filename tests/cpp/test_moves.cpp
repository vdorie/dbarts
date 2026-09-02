#include "common.hpp"

static void testLogisticMutation(ext_rng* rng) {
  const size_t n = 200, n2 = 260;
  std::vector<double> x, f;
  makeMutationData(x, f, n);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = runif01() < 1.0 / (1.0 + std::exp(-f[i])) ? 1.0 : 0.0;

  SamplerOptions options;
  options.numTrees = 25;
  options.nodeScale = 3.0;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                         &rng);
  Results empty;
  sampler.run(50, 0, empty);

  // transactional predictor swap under per-iteration weights
  std::vector<double> xCopy(x);
  check(sampler.setPredictor(xCopy.data(), false, false) ==
          PredictorUpdateResult::accepted,
        "logistic identity setPredictor accepted");

  // whole-data replacement with a resize
  std::vector<double> x2, f2;
  makeMutationData(x2, f2, n2);
  std::vector<double> y2(n2);
  for (size_t i = 0; i < n2; ++i)
    y2[i] = runif01() < 1.0 / (1.0 + std::exp(-f2[i])) ? 1.0 : 0.0;
  check(sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr,
                        0),
        "logistic setData ingests the replacement");
  check(sampler.numObservations() == n2, "logistic setData resizes");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "logistic setData leaves no empty leaves");

  std::vector<double> trainingFits(n2 * 5);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(0, 5, results);
  bool finite = true;
  for (double v : trainingFits) finite &= std::isfinite(v);
  check(finite, "logistic sampler runs after mutation");

  printf("ok: logistic mutation\n");
}

static void testDartUpdate(ext_rng* rng) {
  const size_t p = 4;
  DartPrior dart;
  dart.alpha = 2.0;
  dart.updateAlpha = false;
  dart.initialize(p);

  uint32_t counts[p] = {10, 5, 0, 1};
  double totalCount = 16.0;

  // s | counts ~ Dirichlet(alpha/p + m); independent draws given fixed counts
  const int numDraws = 50000;
  double sums[p] = {0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < numDraws; ++i) {
    dart.update(rng, counts);
    for (size_t j = 0; j < p; ++j) sums[j] += dart.probabilities[j];
  }
  for (size_t j = 0; j < p; ++j) {
    double expected = (dart.alpha / (double) p + (double) counts[j]) /
                      (dart.alpha + totalCount);
    checkNear(sums[j] / numDraws, expected, 0.01, "dart dirichlet mean");
  }

  // with alpha updates on, concentrated counts should drive alpha down
  DartPrior adaptive;
  adaptive.updateAlpha = true;
  adaptive.initialize(20);
  uint32_t sparseCounts[20] = {40, 35, 0};  // rest zero-initialized
  double alphaSum = 0.0;
  for (int i = 0; i < 2000; ++i) {
    adaptive.update(rng, sparseCounts);
    alphaSum += adaptive.alpha;
  }
  check(alphaSum / 2000.0 < 5.0, "dart alpha adapts downward under sparsity");

  printf("ok: dart updates\n");
}

static void testDartSparsityRecovery(ext_rng* rng) {
  const size_t n = 250, p = 25;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i] + 3.0 * x[i + n] + 0.5 * normal;  // vars 0 and 1 only
  }

  auto signalMass = [&](bool useDart) {
    SamplerOptions options;
    options.numTrees = 50;
    options.useDart = useDart;
    options.dart.updateDelay = 100;  // half of burn-in, BART-package style
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
                           1.0, 3.0, 0.37804942330213542, options, &rng);
    const size_t numSamples = 300;
    std::vector<uint32_t> varcount(p * numSamples);
    std::vector<double> splitProbs(p * numSamples, -1.0);
    Results results;
    results.variableCounts = varcount.data();
    results.splitProbabilities = splitProbs.data();
    sampler.run(200, numSamples, results);

    if (useDart) {
      // every recorded sample is a simplex over the predictors
      bool inRange = true, sumsToOne = true;
      for (size_t s = 0; s < numSamples; ++s) {
        double sum = 0.0;
        for (size_t j = 0; j < p; ++j) {
          double prob = splitProbs[j + s * p];
          inRange &= prob >= 0.0 && prob <= 1.0;
          sum += prob;
        }
        sumsToOne &= std::fabs(sum - 1.0) < 1.0e-10;
      }
      check(inRange, "dart varprobs in [0, 1]");
      check(sumsToOne, "dart varprobs sum to one");
    } else {
      // recording is dart-only; without it the buffer is left untouched
      bool untouched = true;
      for (double prob : splitProbs) untouched &= prob == -1.0;
      check(untouched, "varprobs untouched without dart");
    }

    double signal = 0.0, total = 0.0;
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t j = 0; j < p; ++j) {
        total += varcount[j + s * p];
        if (j < 2) signal += varcount[j + s * p];
      }
    return signal / total;
  };

  double uniformMass = signalMass(false);
  double dartMass = signalMass(true);
  check(dartMass > uniformMass + 0.15,
        "dart concentrates splits on signal variables");

  printf("ok: dart sparsity recovery (signal mass %.2f uniform -> %.2f dart)\n",
         uniformMass, dartMass);
}

static void testSetPredictorTransaction(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());

  // identity swap: new buffer, same values; must accept and preserve fits
  std::vector<double> xCopy(x);
  check(sampler.setPredictor(xCopy.data(), false, false) ==
          PredictorUpdateResult::accepted,
        "identity setPredictor accepted");
  check(storageDigest(sampler.data()) == codesBefore,
        "identity swap preserves codes");
  check(sampler.chain(0).treeFits() == treeFitsBefore, "identity swap preserves fits");

  // constant predictors empty one side of every split: must reject and
  // roll back completely
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "degenerate setPredictor rejected");
  check(storageDigest(sampler.data()) == codesBefore,
        "rollback restores codes");
  check(sampler.chain(0).treeFits() == treeFitsBefore, "rollback leaves fits untouched");
  for (size_t t = 0; t < 25; ++t)
    if (!sampler.chain(0).tree(t).bottomNodesAreOccupied()) {
      check(false, "rollback leaves a partition empty");
      break;
    }

  // rejected with updateCutPoints: cuts must be restored too
  std::vector<double> cuts0Before(sampler.data().cutPoints[0]);
  check(sampler.setPredictor(xConstant.data(), false, true) ==
          PredictorUpdateResult::rolledBack,
        "degenerate setPredictor with new cuts rejected");
  check(sampler.data().cutPoints[0] == cuts0Before, "rollback restores cuts");

  // the sampler remains usable
  std::vector<double> sigmaDraws(10);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 10, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after rollback");

  printf("ok: setPredictor transaction\n");
}

static void testSetPredictorForced(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  // force-install degenerate predictors: every split empties a side, so all
  // trees collapse to their roots with merged parameters
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), true, true) ==
          PredictorUpdateResult::accepted,
        "forced setPredictor accepted");

  bool allCollapsed = true, occupied = true;
  for (size_t t = 0; t < 25; ++t) {
    allCollapsed &= sampler.chain(0).tree(t).hasSingleNode();
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  }
  check(allCollapsed, "forced degenerate update collapses trees");
  check(occupied, "collapsed trees keep observations");

  // fits stay consistent: totalFits == sum of constant tree fits
  double fitSum = 0.0;
  for (size_t t = 0; t < 25; ++t) fitSum += sampler.chain(0).treeFits()[t * n];
  checkNear(sampler.chain(0).totalFits()[0], fitSum, 1e-10,
            "forced update keeps fit identity");

  std::vector<double> sigmaDraws(10);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 10, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after forced update");

  printf("ok: forced setPredictor\n");
}

// A forced re-cut over a near-constant column once built a non-ascending
// uniform grid getState serialized and setState then refused. The re-cut now
// refuses, keeping the old ascending grid, so the own-state round-trips.
static void testDegenerateReCutRoundTrips(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  std::vector<double> cutsBefore(sampler.data().cutPoints[0]);
  std::vector<double> xDegenerate(x);
  for (size_t i = 0; i < n; ++i) xDegenerate[i] = 0.5;  // column 0 constant
  check(sampler.setPredictor(xDegenerate.data(), true, true) ==
          PredictorUpdateResult::accepted,
        "forced re-cut over a constant column accepted");
  check(sampler.data().cutPoints[0] == cutsBefore,
        "degenerate re-cut keeps the old ascending grid");

  SamplerStateData st, st2;
  sampler.getState(st);
  // a same-spec restore skips re-quantization, so no raw is consulted
  check(sampler.setState(st, nullptr), "state over a degenerate column restores");
  sampler.getState(st2);
  check(statesAgree(st, st2), "degenerate-column state round trip agrees");
  printf("ok: degenerate re-cut round trip\n");
}

static int countCatMissingRight(const Tree& t, int32_t i, size_t col) {
  const Node& nd(t.at(i));
  if (nd.isBottom()) return 0;
  int here = (static_cast<size_t>(nd.rule.variableIndex) == col &&
              nd.rule.missingGoesRight()) ? 1 : 0;
  return here + countCatMissingRight(t, nd.leftChild, col) +
         countCatMissingRight(t, nd.leftChild + 1, col);
}

// A category-changing mutation once left a live rule routing a missing value
// right after its column stopped holding one, so getState emitted a mask
// outside the reachable gauge and setState refused. The mutation now drops the
// stale direction, so the own-state round-trips.
static void testCategoricalMissingRoundTrips(ext_rng* rng) {
  const size_t n = 300, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    bool missing = runif01() < 0.25;
    x[i] = missing ? std::nan("")
                   : static_cast<double>(static_cast<int>(runif01() * 3.0) % 3);
    x[i + n] = runif01();
    double cat = missing ? 1.5 : x[i];
    y[i] = 3.0 * (cat - 1.0) + 2.0 * (x[i + n] - 0.5) + 0.2 * (runif01() - 0.5);
  }
  x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;  // every category present

  std::vector<ColumnKind> types = {ColumnKind::categorical,
                                   ColumnKind::numeric};
  SamplerOptions options;
  options.numTrees = 50;
  options.predictors.columnTypes = types.data();
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(200, 0, empty);

  int missingRight = 0;
  for (size_t t = 0; t < options.numTrees; ++t)
    missingRight += countCatMissingRight(sampler.chain(0).tree(t), 0, 0);
  check(missingRight > 0, "burn-in routes a categorical missing value right");

  // drop every missing value: hasMissing flips false, leaving those rules
  // routing a category the column no longer reaches
  std::vector<double> noMissing(x.begin(), x.begin() + n);
  for (double& v : noMissing) if (std::isnan(v)) v = 1.0;
  size_t col0 = 0;
  check(sampler.updatePredictor(noMissing.data(), &col0, 1, true, false) ==
          PredictorUpdateResult::accepted,
        "forced drop of the categorical missing values accepted");

  SamplerStateData st, st2;
  sampler.getState(st);
  check(sampler.setState(st, nullptr),
        "categorical state restores after missing drop");
  sampler.getState(st2);
  check(statesAgree(st, st2), "categorical missing-drop round trip agrees");
  printf("ok: categorical missing round trip\n");
}

static void testUpdatePredictorColumns(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));

  // jittered column: usually accepted; either way the sampler must stay
  // consistent, and rejection must restore the column bitwise
  std::vector<double> column0Before(x.begin(), x.begin() + n);
  std::vector<double> jittered(column0Before);
  for (double& v : jittered) v += 0.001 * (runif01() - 0.5);
  size_t columnIndex = 0;
  bool accepted = sampler.updatePredictor(jittered.data(), &columnIndex, 1,
                                          false, false) ==
                  PredictorUpdateResult::accepted;
  // the store owns codes and no raw, so the invariant is on the codes: the
  // installed values quantized, or the snapshot restored bitwise on rollback
  bool columnMatches = true;
  for (size_t i = 0; i < n; ++i) {
    columnMatches &= sampler.data().codeAt(0, i) ==
      (accepted ? sampler.data().codeFor(0, jittered[i]) : codesBefore[i]);
  }
  check(columnMatches, "updatePredictor installs or restores the column");

  // degenerate single column: rejected, second column untouched
  std::vector<double> constantColumn(n, 0.25);
  check(sampler.updatePredictor(constantColumn.data(), &columnIndex, 1,
                                false, false) ==
          PredictorUpdateResult::rolledBack,
        "degenerate column update rejected");
  bool otherUntouched = true;
  for (size_t i = 0; i < n; ++i)
    otherUntouched &= sampler.data().codeAt(1, i) == codesBefore[i + n];
  check(otherUntouched, "rejected column update leaves other columns");

  printf("ok: updatePredictor columns\n");
}

static void testPerObservationUpdate(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  // identity: every observation trivially remains in its leaf
  std::vector<double> identity(x.begin(), x.begin() + n);
  std::unique_ptr<bool[]> installed(new bool[n]);
  check(sampler.updatePredictorPerObservation(identity.data(), 0,
                                              installed.get()),
        "identity per-observation update finalizes");
  bool allInstalled = true;
  for (size_t i = 0; i < n; ++i) allInstalled &= installed[i];
  check(allInstalled, "identity per-observation update installs all");

  // push everything far right: last occupants of left leaves must roll back
  std::vector<double> extreme(n, 10.0);
  check(sampler.updatePredictorPerObservation(extreme.data(), 0,
                                              installed.get()),
        "aggressive per-observation update finalizes");

  size_t numInstalled = 0;
  bool valuesConsistent = true;
  for (size_t i = 0; i < n; ++i) {
    // the store owns codes, not raw: an installed cell carries the new value's
    // code, a rolled-back one keeps the old value's code (per-cell)
    if (installed[i]) {
      ++numInstalled;
      valuesConsistent &=
        sampler.data().codeAt(0, i) == sampler.data().codeFor(0, 10.0);
    } else {
      valuesConsistent &=
        sampler.data().codeAt(0, i) == sampler.data().codeFor(0, identity[i]);
    }
  }
  check(numInstalled > 0 && numInstalled < n,
        "aggressive update installs some, rolls back last occupants");
  check(valuesConsistent, "per-observation rollback is per-cell");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "per-observation update keeps every leaf occupied");

  printf("ok: per-observation update (%zu/%zu installed)\n", numInstalled, n);
}

static void testJointPerObservationUpdate() {
  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 17) != 0 ||
      ext_rng_setSeed(rngB, 31) != 0) {
    check(false, "joint update rng creation");
    return;
  }

  const size_t n = 200;
  std::vector<double> xA, yA, xB, yB;
  makeMutationData(xA, yA, n);
  // second sampler shares column 0 but has its own response and second column
  xB = xA;
  yB.resize(n);
  for (size_t i = 0; i < n; ++i) {
    xB[i + n] = runif01();
    yB[i] = -3.0 * (xB[i] - 0.5) + xB[i + n] + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  SamplerFacade<ConstantGaussianLeaf> samplerA(
    xA.data(), yA.data(), n, size_t(2), nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
    0.37804942330213542, options, &rngA);
  SamplerFacade<ConstantGaussianLeaf> samplerB(
    xB.data(), yB.data(), n, size_t(2), nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
    0.37804942330213542, options, &rngB);
  Results empty;
  samplerA.run(100, 0, empty);
  samplerB.run(100, 0, empty);

  std::vector<double> extreme(n, 10.0);
  std::unique_ptr<bool[]> installed(new bool[n]);
  SamplerBase* samplers[] = {&samplerA, &samplerB};
  size_t columns[] = {0, 0};
  check(updatePredictorPerObservationJointly(samplers, 2, extreme.data(),
                                             columns, installed.get()),
        "joint per-observation update finalizes");

  // all-or-none: the shared column stays identical across samplers. Both were
  // built from the same column-0 values, so they share its cut grid and codes
  // agree; the store owns codes, not raw
  bool columnsAgree = true, valuesConsistent = true;
  size_t numInstalled = 0;
  for (size_t i = 0; i < n; ++i) {
    columnsAgree &=
      samplerA.impl().data().codeAt(0, i) == samplerB.impl().data().codeAt(0, i);
    if (installed[i]) {
      ++numInstalled;
      valuesConsistent &=
        samplerA.impl().data().codeAt(0, i) ==
        samplerA.impl().data().codeFor(0, 10.0);
    } else {
      valuesConsistent &=
        samplerA.impl().data().codeAt(0, i) ==
        samplerA.impl().data().codeFor(0, xB[i]);
    }
  }
  check(columnsAgree, "joint update keeps shared column identical");
  check(valuesConsistent, "joint update installs all-or-none per observation");
  check(numInstalled > 0 && numInstalled < n,
        "joint aggressive update is partial");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= samplerA.impl().chain(0).tree(t).bottomNodesAreOccupied() &&
                samplerB.impl().chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "joint update keeps every leaf occupied in both samplers");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);

  printf("ok: joint per-observation update (%zu/%zu installed)\n",
         numInstalled, n);
}

// Committing an NA into a previously NA-free ordinal column through a
// per-observation session must mark the column missing (setCell owns the flag),
// so the finalize repartition and a fresh descent route the naCode alike:
// a stale gauge lets dropStaleMissingDirections clear a bit the partition
// already routed by, splitting the two. A local generator and a restored
// rngState leave the shared stream untouched for the downstream snapshots.
static void testPerObservationMissingCommit(ext_rng* /*rng*/) {
  std::uint64_t savedRngState = rngState;
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 13579u);
  std::unique_ptr<ConstantLeafSampler> samplerPtr =
    makeBurnedInSampler(x, y, n, localRng);
  ConstantLeafSampler& sampler(*samplerPtr);
  check(!sampler.data().hasMissing[0], "column 0 starts NA-free");

  // drive column 0 fully missing; the guard installs the NAs that keep every
  // leaf occupied and rolls the rest back per cell
  std::vector<double> missing(n, std::nan(""));
  std::unique_ptr<bool[]> installed(new bool[n]);
  check(sampler.updatePredictorPerObservation(missing.data(), 0,
                                              installed.get()),
        "all-missing per-observation update finalizes");

  size_t numInstalled = 0;
  for (size_t i = 0; i < n; ++i) numInstalled += installed[i] ? 1 : 0;
  check(numInstalled > 0, "some NA cells install");
  check(sampler.data().hasMissing[0],
        "committing an NA marks the ordinal column missing");

  // a fresh descent must land every observation where the finalize repartition
  // counted it: a per-leaf descent tally equal to each reached node's stored
  // occupancy (dead arena slots are never descended into)
  bool consistent = true;
  const size_t numTrees = 25;
  for (size_t t = 0; t < numTrees && consistent; ++t) {
    const Tree& tree(sampler.chain(0).tree(t));
    std::vector<size_t> descentCount(tree.nodes.size(), 0);
    for (size_t i = 0; i < n; ++i)
      ++descentCount[static_cast<size_t>(
        tree.findBottomNodeForObservation(sampler.data(), i))];
    for (size_t i = 0; i < n; ++i) {
      int32_t leaf = tree.findBottomNodeForObservation(sampler.data(), i);
      consistent &= descentCount[static_cast<size_t>(leaf)] ==
        tree.at(leaf).numObservations();
    }
  }
  check(consistent, "descent tally matches the finalize partition after NA commit");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: per-observation missing commit (%zu/%zu installed)\n",
         numInstalled, n);
}

static void testQuantilePredictorUpdate(ext_rng* rng) {
  const size_t n = 200;
  // discrete predictors so quantile cut counts are small and controllable
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 10);
    x[i + n] = runif01();
    y[i] = 0.5 * x[i] + 2.0 * x[i + n] + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.useQuantiles = true;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(100, 0, empty);

  check(sampler.data().numCuts[0] == 9, "sampler builds quantile cuts");

  // a coarser column with updateCutPoints must be refused without mutating
  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));
  std::vector<double> coarse(n);
  for (size_t i = 0; i < n; ++i) coarse[i] = static_cast<double>(i % 4);
  size_t columnIndex = 0;
  check(sampler.updatePredictor(coarse.data(), &columnIndex, 1, false, true) ==
          PredictorUpdateResult::invalidCutPoints,
        "coarser quantile column update refused");
  check(storageDigest(sampler.data()) == codesBefore,
        "refused quantile update mutates nothing");

  // same column without cut refresh follows normal transaction semantics
  PredictorUpdateResult result =
    sampler.updatePredictor(coarse.data(), &columnIndex, 1, false, false);
  check(result == PredictorUpdateResult::accepted ||
          result == PredictorUpdateResult::rolledBack,
        "coarser column without cut refresh is transactional");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after quantile updates");

  printf("ok: quantile predictor updates\n");
}

static void testSetCutPoints(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

  // shrink column 0 to three cuts: codes re-quantize, out-of-range splits
  // collapse, and the fit identity holds
  double newCuts[] = {0.25, 0.5, 0.75};
  std::uint32_t numNewCuts = 3;
  const double* cutsByColumn[] = {newCuts};
  size_t columnIndex = 0;
  sampler.setCutPoints(cutsByColumn, &numNewCuts, &columnIndex, 1, x.data());

  check(sampler.data().numCuts[0] == 3, "setCutPoints installs the new count");
  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= sampler.data().codeAt(0, i) == sampler.data().codeFor(0, x[i]);
  check(codesMatch, "setCutPoints re-quantizes the column");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "setCutPoints leaves no empty leaves");

  bool fitIdentity = true;
  for (size_t i = 0; i < n && fitIdentity; i += 37) {
    double total = 0.0;
    for (size_t t = 0; t < 25; ++t) total += sampler.chain(0).treeFits()[t * n + i];
    fitIdentity = std::fabs(total - sampler.chain(0).totalFits()[i]) < 1e-10;
  }
  check(fitIdentity, "setCutPoints keeps the fit identity");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after setCutPoints");

  printf("ok: explicit setCutPoints\n");
}

// Stripping every column to zero cut points leaves the trees root-only with no
// available split variable anywhere. Birth/death must no-op each tree's move;
// unfixed, the forced birth draws a rule for invalidVariable and reads
// data.types out of bounds. A private rng and deterministic data leave both
// the shared ext_rng stream and the global runif01 state untouched, so the
// downstream suites' draws do not shift.
static void testDegenerateRootGuard() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 3000);

  const size_t n = 200;
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i) / static_cast<double>(n);
    x[i + n] = static_cast<double>((i * 7) % n) / static_cast<double>(n);
    y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n];
  }

  SamplerOptions options;
  options.numTrees = 25;
  ConstantLeafSampler sampler(x.data(), y.data(), n, size_t(2), nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);

  const double* noCuts[] = {nullptr, nullptr};
  std::uint32_t numCuts[] = {0, 0};
  size_t columns[] = {0, 1};
  sampler.setCutPoints(noCuts, numCuts, columns, 2, x.data());

  bool rootOnly = true;
  for (size_t t = 0; t < options.numTrees; ++t)
    rootOnly &= sampler.chain(0).tree(t).hasSingleNode();
  check(rootOnly, "zero cuts leave every tree at its root");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs with every column at zero cuts");

  bool stillRoots = true;
  for (size_t t = 0; t < options.numTrees; ++t)
    stillRoots &= sampler.chain(0).tree(t).hasSingleNode();
  check(stillRoots, "degenerate trees stay root-only across the run");

  ext_rng_destroy(rng);
  printf("ok: degenerate root guard\n");
}

static void testMultiChainMutation() {
  const size_t n = 200, numChains = 2;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  std::vector<ext_rng*> rngs(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 2000 + static_cast<uint_least32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = numChains;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, rngs.data());
  Results empty;
  sampler.run(100, 0, empty);

  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));
  std::vector<double> fitsBefore[numChains] = {sampler.chain(0).treeFits(),
                                               sampler.chain(1).treeFits()};

  // the transaction spans chains: degenerate predictors roll back everywhere
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "multi-chain degenerate setPredictor rejected");
  check(storageDigest(sampler.data()) == codesBefore,
        "multi-chain rollback restores codes");
  bool fitsUntouched = true, occupied = true;
  for (size_t c = 0; c < numChains; ++c) {
    fitsUntouched &= sampler.chain(c).treeFits() == fitsBefore[c];
    for (size_t t = 0; t < 25; ++t)
      occupied &= sampler.chain(c).tree(t).bottomNodesAreOccupied();
  }
  check(fitsUntouched, "multi-chain rollback leaves every chain's fits");
  check(occupied, "multi-chain rollback leaves valid partitions");

  // per-observation guard consults every chain's occupancy
  std::vector<double> extreme(n, 10.0);
  std::unique_ptr<bool[]> installed(new bool[n]);
  check(sampler.updatePredictorPerObservation(extreme.data(), 0,
                                              installed.get()),
        "multi-chain per-observation update finalizes");
  size_t numInstalled = 0;
  for (size_t i = 0; i < n; ++i) numInstalled += installed[i] ? 1 : 0;
  check(numInstalled > 0 && numInstalled < n,
        "multi-chain per-observation update is partial");
  occupied = true;
  for (size_t c = 0; c < numChains; ++c)
    for (size_t t = 0; t < 25; ++t)
      occupied &= sampler.chain(c).tree(t).bottomNodesAreOccupied();
  check(occupied, "multi-chain per-observation update keeps leaves occupied");

  std::vector<double> sigmaDraws(5 * numChains);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "multi-chain sampler runs after mutation");

  for (size_t c = numChains; c > 0; --c) ext_rng_destroy(rngs[c - 1]);

  printf("ok: multi-chain mutation (%zu/%zu installed)\n", numInstalled, n);
}

static void testCategoricalMutation(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x(n * 2), y(n);
  const double categoryMeans[] = {2.0, -1.0, 3.0, 0.0};
  for (size_t i = 0; i < n; ++i) {
    size_t category = i % 4;
    x[i] = static_cast<double>(category);
    x[i + n] = runif01();
    y[i] = categoryMeans[category] + 2.0 * x[i + n] +
           0.3 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  ColumnKind types[] = {ColumnKind::categorical, ColumnKind::numeric};
  options.predictors.columnTypes = types;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(100, 0, empty);

  // identity swap accepted; out-of-range category codes refused untouched
  std::vector<double> xCopy(x);
  check(sampler.setPredictor(xCopy.data(), false, false) ==
          PredictorUpdateResult::accepted,
        "categorical identity setPredictor accepted");
  std::vector<double> xBad(x);
  xBad[3] = 9.0;
  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));
  check(sampler.setPredictor(xBad.data(), false, false) ==
          PredictorUpdateResult::invalidCutPoints,
        "out-of-range category code refused");
  check(storageDigest(sampler.data()) == codesBefore,
        "refusal mutates nothing");

  // permuting the categories of some observations is transactional
  std::vector<double> permuted(xCopy.begin(), xCopy.begin() + n);
  for (size_t i = 0; i < n; i += 3)
    permuted[i] = static_cast<double>((static_cast<size_t>(permuted[i]) + 1) % 4);
  size_t columnIndex = 0;
  PredictorUpdateResult result =
    sampler.updatePredictor(permuted.data(), &columnIndex, 1, false, false);
  check(result == PredictorUpdateResult::accepted ||
          result == PredictorUpdateResult::rolledBack,
        "categorical column update is transactional");

  // per-observation updates route through the mask logic
  std::vector<double> newColumn(n);
  for (size_t i = 0; i < n; ++i)
    newColumn[i] = static_cast<double>((i + 1) % 4);
  std::unique_ptr<bool[]> installed(new bool[n]);
  check(sampler.updatePredictorPerObservation(newColumn.data(), 0,
                                              installed.get()),
        "categorical per-observation update finalizes");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "categorical mutation keeps leaves occupied");

  // whole-data replacement with a resize keeps category counts fixed
  const size_t n2 = 260;
  std::vector<double> x2(n2 * 2), y2(n2);
  for (size_t i = 0; i < n2; ++i) {
    size_t category = i % 4;
    x2[i] = static_cast<double>(category);
    x2[i + n2] = runif01();
    y2[i] = categoryMeans[category] + 2.0 * x2[i + n2] +
            0.3 * (runif01() - 0.5);
  }
  check(sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr,
                        0),
        "categorical setData ingests the replacement");
  check(sampler.numObservations() == n2 &&
          sampler.data().categoryCounts[0] == 4,
        "categorical setData resizes and keeps categories");
  occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "categorical setData leaves no empty leaves");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after categorical mutation");

  printf("ok: categorical mutation\n");
}

// An ordered factor's grid is its level table, so the mutation surface holds
// it to that table on BOTH cut-refresh settings: a value off the table is
// refused whether or not the caller asked for a refresh, and a refresh that
// is asked for leaves the midpoints in place rather than re-cutting over the
// replacement's own range.
static void testOrderedFactorMutation(ext_rng* rng) {
  const size_t n = 200;
  const std::uint32_t numLevels = 8;
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    std::uint32_t level = static_cast<std::uint32_t>(i) % numLevels;
    x[i] = static_cast<double>(level);
    x[i + n] = runif01();
    y[i] = 0.4 * static_cast<double>(level) + 2.0 * x[i + n] +
           0.3 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  ColumnKind types[] = {ColumnKind::orderedFactor, ColumnKind::numeric};
  options.predictors.columnTypes = types;
  options.maxNumCutsPerVariable = nullptr;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(100, 0, empty);
  std::vector<double> originalCuts(sampler.data().cutPoints[0]);
  check(sampler.data().numCuts[0] == numLevels - 1,
        "the ordered factor carries one cut per level boundary");

  // a value off the level table is refused on both settings, mutating nothing
  std::vector<double> xBad(x);
  xBad[3] = static_cast<double>(numLevels);
  std::vector<xint_t> codesBefore(storageDigest(sampler.data()));
  check(sampler.setPredictor(xBad.data(), false, false) ==
          PredictorUpdateResult::invalidCutPoints &&
        sampler.setPredictor(xBad.data(), false, true) ==
          PredictorUpdateResult::invalidCutPoints,
        "an off-table level code is refused with and without cut refresh");
  check(storageDigest(sampler.data()) == codesBefore,
        "the refusal mutates nothing");

  // a refresh over a replacement spanning half the table keeps the midpoints
  std::vector<double> narrowed(x);
  for (size_t i = 0; i < n; ++i)
    narrowed[i] = static_cast<double>(static_cast<std::uint32_t>(i) % 3u);
  PredictorUpdateResult result =
    sampler.setPredictor(narrowed.data(), true, true);
  check(result == PredictorUpdateResult::accepted,
        "a forced update over a subset of the levels is accepted");
  check(sampler.data().cutPoints[0] == originalCuts &&
          sampler.data().numCuts[0] == numLevels - 1,
        "the refresh keeps the level midpoints rather than re-cutting");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after ordered-factor mutation");

  printf("ok: ordered factor mutation\n");
}

// ---------------------------------------------------------------------------
// The empty-leaf veto counts POSITIVE-WEIGHT members, not members. A leaf all
// of whose rows carry weight zero enters no likelihood term, so a branch
// holding one scores the sentinel; a chain driven under such weights may never
// settle on one. With NO weight vector the veto stays the member count it has
// always been - the same decision and the same arithmetic
// (docs/design/empty-leaf-veto.md). A local generator leaves the shared stream
// untouched for the suites that follow.
// ---------------------------------------------------------------------------
static void testEmptyLeafVetoCountsWeight() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20260812u);
  const size_t n = 256, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 8) / 8.0;              // x0: codes 0..7
    x[i + n] = static_cast<double>((i / 8) % 8) / 8.0;    // x1: codes 0..7
    y[i] = 2.0 * x[i] - x[i + n] + 0.1 * static_cast<double>(i % 3);
  }
  ColumnStore store;
  built(store.build(x.data(), n, p, 7));

  // the lowest cut of x0 isolates the code-0 rows; weight zero exactly there
  std::vector<double> zeroed(n, 1.0), positive(n, 1.0);
  for (size_t i = 0; i < n; ++i)
    if (i % 8 == 0) zeroed[i] = 0.0;

  CGMTreePrior prior;
  prior.base = 0.95;
  prior.power = 2.0;
  ConstantGaussianLeaf leaf{0.5};
  MoveScratch scratch;
  const double sigma = 0.7, k = 2.0;

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  auto buildSplit = [&](const double* weights) {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), weights);
    Rule rule;
    rule.variableIndex = 0;
    rule.setSplitIndex(0);
    tree.birth(store, 0, rule, y.data(), weights);
  };

  buildSplit(zeroed.data());
  int32_t left = tree.at(0).leftChild;
  check(tree.at(left).numObservations() == n / 8,
        "veto fixture: the zero-weight leaf is occupied by count");
  check(tree.at(left).sumWeights == 0.0,
        "veto fixture: the zero-weight leaf holds no weight");

  MoveContext zeroCtx{store,      prior, 0.5, 0.1, 0.5,
                      zeroed.data(), k,   scratch};
  BranchScore zeroScore =
    logLikelihoodForBranch(zeroCtx, leaf, tree, 0, y.data(), sigma);
  check(zeroScore.rank == 1,
        "a leaf of only zero-weight rows vetoes its branch");
  check(std::isfinite(zeroScore.logLikelihood),
        "the vetoed branch's finite part skips the vetoed leaf");

  buildSplit(positive.data());
  MoveContext positiveCtx{store,           prior, 0.5, 0.1, 0.5,
                          positive.data(), k,     scratch};
  BranchScore positiveScore =
    logLikelihoodForBranch(positiveCtx, leaf, tree, 0, y.data(), sigma);
  check(positiveScore.rank == 0 && std::isfinite(positiveScore.logLikelihood),
        "the same leaf under positive weights scores finite");

  // the no-weights path: the same branch, scored with no weight vector, is
  // bitwise the sum of its leaves' marginals, and the veto there is still the
  // member count
  buildSplit(nullptr);
  MoveContext nullCtx{store, prior, 0.5, 0.1, 0.5, nullptr, k, scratch};
  std::vector<int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  double reference = 0.0;
  for (int32_t b : bottoms)
    reference +=
      leaf.logIntegratedLikelihoodForNode(tree, y.data(), nullptr, k,
                                          sigma * sigma, b);
  check(logLikelihoodForBranch(nullCtx, leaf, tree, 0, y.data(), sigma)
            .logLikelihood == reference,
        "the no-weights branch score is bitwise the leaf marginals");
  bool countLaw = true;
  for (size_t i = 0; i < tree.nodes.size(); ++i)
    countLaw &= tree.leafHasNoWeight(static_cast<int32_t>(i), nullptr) ==
                (tree.at(static_cast<int32_t>(i)).numObservations() == 0);
  check(countLaw, "with no weights, emptiness is exactly the member count");

  Node saved = tree.at(left);
  tree.at(left).end = tree.at(left).begin;  // strand the leaf outright
  check(logLikelihoodForBranch(nullCtx, leaf, tree, 0, y.data(), sigma).rank ==
          2,
        "a member-empty leaf still vetoes with no weights");
  tree.at(left) = saved;

  // the invariant a move chain must maintain under those weights, and the
  // non-vacuity measurement beside it: the same chain under the count law (no
  // weights installed) settles on leaves the weight law forbids
  auto driveChain = [&](const double* weights) {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), weights);
    MoveContext ctx{store, prior, 0.5, 0.1, 0.5, weights, k, scratch};
    size_t violations = 0, accepted = 0;
    for (int iter = 0; iter < 4000; ++iter) {
      bool stepTaken = false;
      StepType stepType;
      metropolisJumpForTree(ctx, leaf, rng, tree, y.data(), sigma, &stepTaken,
                            &stepType);
      accepted += stepTaken ? 1 : 0;
      bottoms.clear();
      tree.fillBottom(0, bottoms);
      for (int32_t b : bottoms)
        violations += tree.leafHasNoWeight(b, zeroed.data()) ? 1 : 0;
    }
    return std::pair<size_t, size_t>{violations, accepted};
  };

  auto weighted = driveChain(zeroed.data());
  check(weighted.second > 0, "the weighted move chain moves");
  check(weighted.first == 0,
        "no accepted move leaves a leaf of only zero-weight rows");
  auto counted = driveChain(nullptr);
  check(counted.first > 0,
        "non-vacuity: the count law does settle on such leaves");

  ext_rng_destroy(rng);
  printf("ok: empty-leaf veto counts weight (%zu count-law leaves of only "
         "zero-weight rows, %zu under the weight law)\n",
         counted.first, weighted.first);
}

// ---------------------------------------------------------------------------
// Weights do not ride the tree, so a vector installed on a GROWN tree can
// leave leaves the veto refuses - a CURRENT state outside the admissible set,
// which the veto owes a law. The rank supplies it: two vetoed branches compare
// by rank before likelihood, so no acceptance ratio is NaN, the structure
// keeps moving (under prior x transition, at constant likelihood), no move
// installs a MEMBER-empty leaf (the state law the rest of the engine enforces,
// which the second rank level keeps absolute), and a partially stranded tree
// is absorbed back into the admissible set. Total zeroing is the frozen-forest
// case outright; the partial arm is the one a masking host actually reaches.
// ---------------------------------------------------------------------------
static void testVetoRankUnfreezesStrandedTree() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20260817u);
  // continuous predictors on a fine grid, deliberately NOT the lattice the
  // veto-predicate fixture above uses: a cut interval is constrained by the
  // ancestors, never by occupancy, so a small node here has cuts on both sides
  // of its whole member set and MEMBER-empty proposals are common. That is
  // what makes the rank's second level testable at all.
  const size_t n = 256, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = ext_rng_simulateContinuousUniform(rng);
    x[i + n] = ext_rng_simulateContinuousUniform(rng);
    y[i] = 2.0 * x[i] - x[i + n] + 0.1 * ext_rng_simulateStandardNormal(rng);
  }
  ColumnStore store;
  built(store.build(x.data(), n, p, 40));

  CGMTreePrior prior;
  prior.base = 0.95;
  prior.power = 1.5;
  ConstantGaussianLeaf leaf{2.0};
  MoveScratch scratch;
  const double sigma = 0.25, k = 2.0;
  std::vector<index_t> indexBuffer(n);
  Tree tree;

  std::vector<double> ones(n, 1.0), zeros(n, 0.0);

  // every arm starts from the same tree, grown under positive weights
  auto growTree = [&](unsigned seed) {
    ext_rng_setSeed(rng, seed);
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), ones.data());
    MoveContext ctx{store, prior, 0.5, 0.1, 0.5, ones.data(), k, scratch};
    for (int iter = 0; iter < 3000; ++iter) {
      bool stepTaken = false;
      StepType stepType;
      metropolisJumpForTree(ctx, leaf, rng, tree, y.data(), sigma, &stepTaken,
                            &stepType);
    }
  };

  auto countVetoed = [&](const double* weights) {
    std::vector<int32_t> bottoms;
    tree.fillBottom(0, bottoms);
    size_t vetoed = 0;
    for (int32_t b : bottoms)
      vetoed += tree.leafVetoRank(b, weights) != 0 ? 1 : 0;
    return vetoed;
  };

  struct Driven {
    size_t nanAlpha = 0;
    size_t structureChanges = 0;
    size_t memberEmpty = 0;
    size_t vetoedAtEnd = 0;
    int absorbed = -1;  // first sweep at which no leaf is vetoed
  };
  auto drive = [&](const double* weights, int iterations) {
    MoveContext ctx{store, prior, 0.5, 0.1, 0.5, weights, k, scratch};
    Driven driven;
    std::uint64_t signature = treeStructureSignature(tree);
    std::vector<int32_t> bottoms;
    for (int iter = 0; iter < iterations; ++iter) {
      bool stepTaken = false;
      StepType stepType;
      double alpha = metropolisJumpForTree(ctx, leaf, rng, tree, y.data(),
                                           sigma, &stepTaken, &stepType);
      if (std::isnan(alpha)) ++driven.nanAlpha;
      std::uint64_t next = treeStructureSignature(tree);
      if (next != signature) {
        ++driven.structureChanges;
        signature = next;
      }
      bottoms.clear();
      tree.fillBottom(0, bottoms);
      size_t vetoed = 0;
      for (int32_t b : bottoms) {
        // membership is read straight off the node, not through the rank: the
        // state law must hold even if the rank is the thing that is wrong
        if (tree.at(b).numObservations() == 0) ++driven.memberEmpty;
        if (tree.leafVetoRank(b, weights) != 0) ++vetoed;
      }
      if (vetoed == 0 && driven.absorbed < 0) driven.absorbed = iter;
      driven.vetoedAtEnd = vetoed;
    }
    return driven;
  };

  // ---- every leaf stranded: the whole tree is out of the admissible set ----
  growTree(20260818u);
  check(!tree.hasSingleNode(), "the stranded-tree fixture grows past the root");
  std::vector<int32_t> grownBottoms;
  tree.fillBottom(0, grownBottoms);
  check(countVetoed(zeros.data()) == grownBottoms.size(),
        "non-vacuity: an all-zero vector strands every leaf of the grown tree");
  Driven stranded = drive(zeros.data(), 2000);
  check(stranded.nanAlpha == 0,
        "a wholly stranded tree reports no NaN acceptance probability");
  check(stranded.structureChanges > 0,
        "a wholly stranded tree keeps moving rather than freezing");
  check(stranded.memberEmpty == 0,
        "no move installs a member-empty leaf from a stranded state");

  // ---- partially stranded: the tree must find its way back to the set ----
  // strand every third leaf by construction, zeroing exactly its members: the
  // worst case of the install a masking host makes, and non-vacuous by
  // construction rather than by luck of the grown partition
  growTree(20260819u);
  std::vector<double> partial(n, 1.0);
  std::vector<int32_t> toStrand;
  tree.fillBottom(0, toStrand);
  for (size_t b = 0; b < toStrand.size(); b += 3) {
    const Node& node(tree.at(toStrand[b]));
    for (size_t j = node.begin; j < node.end; ++j)
      partial[tree.indices[j]] = 0.0;
  }
  size_t strandedLeaves = countVetoed(partial.data());
  check(strandedLeaves > 0,
        "non-vacuity: the half-space vector strands part of the grown tree");
  Driven partialDriven = drive(partial.data(), 2000);
  check(partialDriven.nanAlpha == 0,
        "a partially stranded tree reports no NaN acceptance probability");
  check(partialDriven.absorbed >= 0 && partialDriven.vetoedAtEnd == 0,
        "a partially stranded tree is absorbed back into the admissible set");
  check(partialDriven.memberEmpty == 0,
        "no move installs a member-empty leaf while stranded");

  ext_rng_destroy(rng);
  printf("ok: veto rank unfreezes a stranded tree (%zu leaves, %zu structure "
         "changes while wholly stranded; %zu leaves stranded by the partial "
         "vector, absorbed at sweep %d)\n",
         grownBottoms.size(), stranded.structureChanges, strandedLeaves,
         partialDriven.absorbed);
}

// ---------------------------------------------------------------------------
// EQUAL RANK 1 - both compared branches weight-vetoed - is the veto rank's one
// CHANGED comparison and the one that unfreezes a stranded forest. Its law
// (docs/design/empty-leaf-veto.md, "Is vetoed-vs-vetoed reachable?"): equal
// rank takes the acceptance arithmetic unchanged on the finite parts, and the
// finite part SKIPS the vetoed leaves rather than summing their marginals, so a
// wholly vetoed pair contributes exactly 1 and the tree mixes under prior x
// transition at constant likelihood.
//
// The fixture puts that law in closed form. One binary split column carries
// exactly one cut, so the only rule the prior can draw halves the members and
// neither child can split again; an all-zero weight vector then puts the root,
// the split branch and every leaf of both at rank 1, which the checks below
// read directly off the fixture. No comparison in the move arms can therefore
// be anything but equal rank 1 - the current state is occupied and weightless
// whatever the move does, and so is the one proposal available - and each
// prior x transition ratio is a hand-derived constant. A local generator
// leaves the shared stream untouched for the suites that follow.
// ---------------------------------------------------------------------------
static void testEqualRankOneComparison() {
  const size_t n = 64, p = 1;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 2);
    y[i] = 0.3 * x[i] + 0.01 * static_cast<double>(i) - 0.5;
  }
  ColumnStore store;
  // the store owns a raw copy of column 0, which the linear leaf below reads
  // as its covariate through rawColumn
  size_t leafGather[] = {0};
  built(store.build(x.data(), n, p, 1, false, nullptr, leafGather, 1));
  check(store.numCuts[0] == 1,
        "equal-rank fixture: the split column carries exactly one cut");

  std::vector<double> zeros(n, 0.0);
  const double sigma = 1.3, k = 2.0;
  ConstantGaussianLeaf constant{0.7};
  // the same branch under a leaf model whose zero-weight marginal is NOT zero:
  // over p parameters the linear leaf's is 0.5 p log(ridge) - p log(sqrt(ridge)),
  // which cancels in exact arithmetic and does not in doubles
  LinearGaussianLeaf linear;
  linear.scale = 0.7;
  size_t covariates[] = {0};
  linear.initialize(store, covariates, 1);

  CGMTreePrior growPrior;  // birth arm: prior ratio 0.25 / 0.75
  growPrior.base = 0.25;
  growPrior.power = 2.0;
  CGMTreePrior prunePrior;  // death arm: the same ratio, inverted
  prunePrior.base = 0.75;
  prunePrior.power = 2.0;

  MoveScratch scratch;
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  auto buildRoot = [&]() {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), zeros.data());
  };
  auto buildSplit = [&]() {
    buildRoot();
    Rule rule;
    rule.variableIndex = 0;
    rule.setSplitIndex(0);
    tree.birth(store, 0, rule, y.data(), zeros.data());
  };

  // ---- the resolution: equal rank hands the move its own operands ----
  double currentLogL, proposalLogL, rank0Current, rank0Proposal;
  resolveVetoRank({1, -2.5}, {1, 4.25}, &currentLogL, &proposalLogL);
  check(currentLogL == -2.5 && proposalLogL == 4.25,
        "equal rank 1 passes both operands through unchanged");
  resolveVetoRank({0, -2.5}, {0, 4.25}, &rank0Current, &rank0Proposal);
  check(currentLogL == rank0Current && proposalLogL == rank0Proposal,
        "equal rank 1 resolves bitwise as equal rank 0 does");
  resolveVetoRank({1, -2.5}, {0, 4.25}, &currentLogL, &proposalLogL);
  check(currentLogL == -HUGE_VAL && proposalLogL == 4.25,
        "a vetoed current against an admissible proposal loses outright");
  resolveVetoRank({0, -2.5}, {1, 4.25}, &currentLogL, &proposalLogL);
  check(currentLogL == -2.5 && proposalLogL == -HUGE_VAL,
        "a vetoed proposal against an admissible current loses outright");
  resolveVetoRank({1, -2.5}, {2, 4.25}, &currentLogL, &proposalLogL);
  check(currentLogL == -2.5 && proposalLogL == -HUGE_VAL,
        "a member-empty proposal loses to a weight-vetoed current");
  // the one -HUGE_VAL pair that is not the veto: a constrained leaf model's
  // feasibility sentinel on both sides rejects rather than reporting NaN
  resolveVetoRank({1, -HUGE_VAL}, {1, -HUGE_VAL}, &currentLogL, &proposalLogL);
  check(currentLogL == 0.0 && proposalLogL == -HUGE_VAL,
        "two sentinels at equal rank reject rather than differencing to NaN");

  // ---- the finite parts the equal-rank comparison differences ----
  buildSplit();
  int32_t leftChild = tree.at(0).leftChild;
  check(tree.at(leftChild).numObservations() == n / 2 &&
          tree.at(leftChild + 1).numObservations() == n / 2,
        "equal-rank fixture: the one rule halves the members");
  check(tree.leafVetoRank(leftChild, zeros.data()) == 1 &&
          tree.leafVetoRank(leftChild + 1, zeros.data()) == 1,
        "equal-rank fixture: both children are weight-vetoed, neither empty");

  MoveContext ctx{store, growPrior, 1.0, 0.0, 0.5, zeros.data(), k, scratch};
  BranchScore splitScore =
    logLikelihoodForBranch(ctx, constant, tree, 0, y.data(), sigma);
  check(splitScore.rank == 1 && splitScore.logLikelihood == 0.0,
        "a wholly vetoed branch scores rank 1 and exactly zero");
  double naive = linear.logIntegratedLikelihoodForNode(
                   tree, y.data(), zeros.data(), k, sigma * sigma, leftChild) +
                 linear.logIntegratedLikelihoodForNode(
                   tree, y.data(), zeros.data(), k, sigma * sigma,
                   leftChild + 1);
  check(naive != 0.0,
        "non-vacuity: the linear leaf's zero-weight marginals do not sum to 0");
  BranchScore linearScore =
    logLikelihoodForBranch(ctx, linear, tree, 0, y.data(), sigma);
  check(linearScore.rank == 1 && linearScore.logLikelihood == 0.0,
        "the skip holds for a leaf model whose vetoed marginal is not zero");

  buildRoot();
  BranchScore rootScore =
    logLikelihoodForBranch(ctx, constant, tree, 0, y.data(), sigma);
  check(rootScore.rank == 1 && rootScore.logLikelihood == 0.0,
        "the undivided branch is vetoed too, and scores exactly zero");
  resolveVetoRank(rootScore, splitScore, &currentLogL, &proposalLogL);
  check(std::exp(proposalLogL - currentLogL) == 1.0,
        "an equal-rank-1 pair contributes exactly 1 to the acceptance ratio");

  // ---- the moves, whose comparisons the fixture forces to equal rank 1 ----
  // From the root the birth step is forced (a single-node tree births with
  // probability 1) and the reverse death is forced back (the split tree has no
  // birthable node), so transitionRatio is exactly 1 and priorRatio is
  // base (1 - 0)(1 - 0) / (1 - base) - neither child can split again. At
  // base = 0.25 that is 0.25 / 0.75, and the equal-rank-1 likelihood is the
  // only remaining factor: prior x transition is the whole acceptance.
  const double expected = 0.25 / 0.75;
  struct MoveOutcome {
    double alpha;
    bool stepTaken;
    bool wasBirth;
  };
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  auto runBirth = [&](const auto& leafModel, const double* response,
                      double leafSigma, double leafK) {
    ext_rng_setSeed(rng, 20260819u);
    buildRoot();
    MoveContext armCtx{store,        growPrior, 1.0, 0.0, 0.5,
                       zeros.data(), leafK,     scratch};
    MoveOutcome out{0.0, false, false};
    out.alpha = birthOrDeathMove(armCtx, leafModel, rng, tree, response,
                                 leafSigma, &out.stepTaken, &out.wasBirth);
    return out;
  };

  MoveOutcome birth = runBirth(constant, y.data(), sigma, k);
  check(birth.wasBirth, "the equal-rank birth arm takes a birth step");
  check(birth.alpha == expected,
        "an equal-rank-1 birth accepts on prior x transition alone");
  std::vector<double> shifted(n);
  for (size_t i = 0; i < n; ++i) shifted[i] = 100.0 * y[i] + 7.0;
  check(runBirth(constant, shifted.data(), sigma, k).alpha == expected,
        "the equal-rank-1 acceptance does not read the response");
  check(runBirth(constant, y.data(), 4.0 * sigma, 0.5 * k).alpha == expected,
        "the equal-rank-1 acceptance does not read sigma or k");
  check(runBirth(linear, y.data(), sigma, k).alpha == expected,
        "the equal-rank-1 acceptance does not read the leaf model");

  // The death of a parent whose children are BOTH weight-vetoed inherits their
  // veto (its members are their union), so the collapse that repairs a stranded
  // pair is an equal-rank-1 comparison, not a rank improvement. Same closed
  // form with the prior inverted: 0.25 / 0.75 again at base = 0.75.
  auto runDeath = [&](const auto& leafModel) {
    ext_rng_setSeed(rng, 20260820u);
    buildSplit();
    MoveContext armCtx{store,        prunePrior, 1.0, 0.0, 0.5,
                       zeros.data(), k,          scratch};
    MoveOutcome out{0.0, false, false};
    out.alpha = birthOrDeathMove(armCtx, leafModel, rng, tree, y.data(), sigma,
                                 &out.stepTaken, &out.wasBirth);
    return out;
  };
  MoveOutcome death = runDeath(constant);
  check(!death.wasBirth, "the equal-rank death arm takes a death step");
  check(death.alpha == expected,
        "an equal-rank-1 death accepts on prior x transition alone");
  check(runDeath(linear).alpha == expected,
        "the equal-rank-1 death does not read the leaf model either");

  // the change move's own resolution: the single available rule is the current
  // one, so the subtree prior and the proposal correction both cancel and the
  // acceptance is the equal-rank-1 likelihood by itself
  ext_rng_setSeed(rng, 20260821u);
  buildSplit();
  MoveContext changeCtx{store,        growPrior, 0.0, 0.0, 0.0,
                        zeros.data(), k,         scratch};
  bool changeTaken = false;
  double changeAlpha = changeMove(changeCtx, constant, rng, tree, y.data(),
                                  sigma, &changeTaken);
  check(changeAlpha == 1.0,
        "an equal-rank-1 change move scores exactly 1, not 0 and not NaN");

  ext_rng_destroy(rng);
  printf("ok: equal-rank-1 veto comparison (birth and death alpha %.17g, "
         "change 1.0, linear-leaf vetoed marginals %.3e)\n",
         expected, naive);
}

static void testLinearLeafMutation(ext_rng* rng) {
  const size_t n = 200, p = 3;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    x[i + 2 * n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = (x[i] > 0.5 ? x[i + n] : 0.0) + 0.2 * normal;
  }

  SamplerOptions options;
  options.numTrees = 20;
  size_t covariates[] = {1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;

  std::unique_ptr<SamplerBase> sampler = createSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    1.0, 3.0, 0.37804942330213542, options, &rng);
  Results warmup;
  sampler->run(150, 0, warmup);

  // a constant split column empties a child of every tree that uses it:
  // the transaction rolls back and sampling continues
  std::vector<double> badX(x);
  for (size_t i = 0; i < n; ++i) badX[i] = 0.5;
  check(sampler->setPredictor(badX.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "a degenerate split column rolls back");

  // jittering the leaf covariate leaves every partition intact (codes are
  // driven by the unchanged split columns' cuts only when the jitter stays
  // within bins; force-update to be deterministic) and the fits pick up the
  // regathered covariates
  std::vector<double> newX(x);
  for (size_t i = 0; i < n; ++i) newX[i + n] = x[i + n] * 1.1;
  check(sampler->setPredictor(newX.data(), true, true) ==
          PredictorUpdateResult::accepted,
        "a forced covariate update installs");

  const size_t numSamples = 100;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler->run(50, numSamples, results);
  bool allFinite = true;
  for (double fit : trainingFits) allFinite = allFinite && std::isfinite(fit);
  check(allFinite, "fits stay finite through predictor mutation");

  // per-observation updates sweep the covariate column through the real
  // session machinery
  std::vector<double> newColumn(n);
  for (size_t i = 0; i < n; ++i) newColumn[i] = newX[i + n] + 0.01;
  std::vector<unsigned char> installedBytes(n);
  bool* installed = reinterpret_cast<bool*>(installedBytes.data());
  check(sampler->updatePredictorPerObservation(newColumn.data(), 1, installed),
        "per-observation update finalizes");

  // whole-data replacement re-standardizes and continues: the new data
  // doubles the varying slope, which the continued chain recovers
  const size_t n2 = 150;
  std::vector<double> x2(n2 * p), y2(n2), f2(n2);
  for (size_t i = 0; i < n2; ++i) {
    x2[i] = runif01();
    x2[i + n2] = 4.0 * runif01() - 2.0;
    x2[i + 2 * n2] = runif01();
    f2[i] = x2[i] > 0.5 ? 2.0 * x2[i + n2] : 0.0;
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y2[i] = f2[i] + 0.2 * normal;
  }
  check(sampler->setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr,
                         0, nullptr),
        "setData ingests the replacement");

  std::vector<double> fits2(n2 * numSamples);
  Results results2;
  results2.trainingFits = fits2.data();
  sampler->run(150, numSamples, results2);

  std::vector<double> posteriorMean(n2, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n2; ++i)
      posteriorMean[i] += fits2[i + s * n2] / (double) numSamples;
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, count = 0.0;
  for (size_t i = 0; i < n2; ++i) {
    if (x2[i] <= 0.5) continue;
    sx += x2[i + n2];
    sy += posteriorMean[i];
    sxx += x2[i + n2] * x2[i + n2];
    sxy += x2[i + n2] * posteriorMean[i];
    count += 1.0;
  }
  double slope = (sxy - sx * sy / count) / (sxx - sx * sx / count);
  check(slope > 1.5 && slope < 2.5, "setData recovers the doubled slope");

  // prediction from the live trees agrees with the recorded fits
  std::vector<double> livePredictions(n2);
  sampler->predict(x2.data(), n2, 1, livePredictions.data());
  const double* lastFits = fits2.data() + (numSamples - 1) * n2;
  bool liveMatches = true;
  for (size_t i = 0; i < n2 && liveMatches; ++i)
    liveMatches = std::fabs(livePredictions[i] - lastFits[i]) < 1e-8;
  check(liveMatches, "live-tree prediction matches the last recorded fits");

  printf("ok: linear leaf mutation (post-setData slope %.2f)\n", slope);
}

// leafOf must equal the obs-to-leaf map derived independently from each tree's
// fillBottom + index segments after EVERY sweep of a randomized run (accepted
// and rejected moves of all kinds land across the sweeps), and again after a
// wholesale prior-drawn structure reset. A local generator and a restored
// rngState leave the shared stream untouched, so the downstream snapshot tests
// are unperturbed.
static void testLeafOfConsistency(ext_rng* /*rng*/) {
  std::uint64_t savedRngState = rngState;
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 24680u);
  std::unique_ptr<ConstantLeafSampler> samplerPtr =
    makeBurnedInSampler(x, y, n, localRng);
  ConstantLeafSampler& sampler(*samplerPtr);

  const size_t numTrees = sampler.chain(0).numTrees();
  std::vector<int32_t> bottoms;
  std::vector<std::uint32_t> expected(n);
  bool allMatch = true;
  auto mapsMatch = [&]() {
    for (size_t t = 0; t < numTrees && allMatch; ++t) {
      const Tree& tree = sampler.chain(0).tree(t);
      bottoms.clear();
      tree.fillBottom(0, bottoms);
      for (int32_t b : bottoms) {
        const Node& node = tree.at(b);
        for (size_t m = node.begin; m < node.end; ++m)
          expected[tree.indices[m]] = static_cast<std::uint32_t>(b);
      }
      const std::uint32_t* actual = sampler.chain(0).leafOfForTesting(t);
      for (size_t i = 0; i < n; ++i)
        if (actual[i] != expected[i]) { allMatch = false; break; }
    }
  };

  Results empty;
  for (int sweep = 0; sweep < 25 && allMatch; ++sweep) {
    sampler.run(1, 0, empty);
    mapsMatch();
  }
  check(allMatch, "leafOf matches the derived map after every sweep");

  // wholesale reset: prior-drawn structures mark the map for rebuild, which
  // the next sweep must clear tree by tree
  size_t fusedBefore = sampler.chain(0).fusedSuffstatRunsForTesting();
  sampler.chain(0).sampleTreesFromPrior();
  sampler.run(1, 0, empty);
  size_t fusedAfterReset = sampler.chain(0).fusedSuffstatRunsForTesting();
  mapsMatch();
  check(allMatch, "leafOf matches after a prior-drawn structure reset");

  // the same state is the fused roll+suffstat's stale-map decline: over that
  // sweep the map still describes the PREVIOUS partition, so every tree falls
  // back to the stock pair, and the fusion resumes once the rebuild lands
  check(fusedAfterReset == fusedBefore,
        "a stale leafOf declines the fused suffstat for every tree");
  sampler.run(1, 0, empty);
  check(sampler.chain(0).fusedSuffstatRunsForTesting() - fusedAfterReset ==
          numTrees,
        "the fused suffstat resumes on the sweep after the leafOf rebuild");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: leafOf consistency\n");
}

// A wholesale prior structure draw returns the forest to the zero-fit state a
// freshly built chain carries: totalFits zero, so the cached-fits evaluator
// mu[leafOf] the next sweep's residual roll reads agrees with it (INV-1), and
// every map entry indexes the reset arena (INV-2, which a draw that SHRINKS a
// tree breaks unless the map moves with the values). A local generator and a
// restored rngState leave the shared stream untouched.
static void testPriorResetContract(ext_rng* /*rng*/) {
  std::uint64_t savedRngState = rngState;
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 13579u);
  std::unique_ptr<ConstantLeafSampler> samplerPtr =
    makeBurnedInSampler(x, y, n, localRng);
  ConstantLeafSampler& sampler(*samplerPtr);
  const size_t numTrees = sampler.chain(0).numTrees();

  // the burn-in leaves real fits behind, so the zeroing below is a change of
  // state rather than a no-op restatement of construction
  bool anyFit = false;
  for (double v : sampler.chain(0).totalFits())
    if (v != 0.0) { anyFit = true; break; }
  check(anyFit, "the burn-in left non-zero totalFits");

  sampler.chain(0).sampleTreesFromPrior();

  // every tree stays marked for rebuild: clearing the mark here would make the
  // fused roll+suffstat pass eligible one sweep early and move the draws of
  // every prior-initialized run
  bool allStale = true;
  for (size_t t = 0; t < numTrees; ++t)
    allStale = allStale && sampler.chain(0).leafOfStaleForTesting(t) != 0;
  check(allStale, "a prior reset leaves every tree marked for rebuild");

  bool totalZero = true;
  for (double v : sampler.chain(0).totalFits())
    totalZero = totalZero && v == 0.0;
  check(totalZero, "a prior reset zeroes totalFits");

  // one sweep on, totalFits sums the per-tree fits again: the roll retired
  // exactly the cached fits it read, so no displacement enters the chain. The
  // two sides are the same quantity by different summation orders (the roll's
  // running residual against a gather over mu[leafOf]), hence a last-ulp band
  // rather than equality.
  Results empty;
  sampler.run(1, 0, empty);
  std::vector<double> perTree = sampler.chain(0).treeFits();
  const std::vector<double>& total = sampler.chain(0).totalFits();
  bool invariantHolds = true;
  double worstDeviation = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double sum = 0.0;
    for (size_t t = 0; t < numTrees; ++t) sum += perTree[t * n + i];
    double deviation = std::fabs(sum - total[i]);
    if (deviation > worstDeviation) worstDeviation = deviation;
    invariantHolds = invariantHolds &&
      deviation <= 1.0e-12 * (1.0 + std::fabs(total[i]));
  }
  check(invariantHolds, "totalFits sums the per-tree fits after a reset sweep");

  ext_rng_destroy(localRng);

  // an arena SHRINK is the bounds half: growSubtreeFromPrior returns before its
  // Bernoulli when the growth probability is zero, so a flattened prior forces
  // bare roots and consumes no draws. Any observation whose pre-reset leaf id
  // survived would then index past the one-node arena. The check reads only
  // leafOf and nodes.size(), never mu, so it diagnoses the read without making
  // it.
  std::vector<double> xShrink, yShrink;
  makeMutationData(xShrink, yShrink, n);
  ext_rng* shrinkRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(shrinkRng, 86420u);
  std::unique_ptr<ConstantLeafSampler> shrunkPtr =
    makeBurnedInSampler(xShrink, yShrink, n, shrinkRng);
  ConstantLeafSampler& shrunk(*shrunkPtr);

  bool anyGrown = false;
  for (size_t t = 0; t < numTrees; ++t)
    anyGrown = anyGrown || shrunk.chain(0).tree(t).nodes.size() > 1;
  check(anyGrown, "the burn-in grew at least one tree");

  ModelParameters flat;
  flat.base = 0.0;
  shrunk.chain(0).setModel(flat);
  shrunk.chain(0).sampleTreesFromPrior();

  bool insideArena = true;
  for (size_t t = 0; t < numTrees; ++t) {
    const std::uint32_t* leafOf = shrunk.chain(0).leafOfForTesting(t);
    size_t arenaSize = shrunk.chain(0).tree(t).nodes.size();
    for (size_t i = 0; i < n; ++i)
      if (leafOf[i] >= arenaSize) { insideArena = false; break; }
  }
  check(insideArena, "the obs-to-leaf map stays inside the reset arena");

  ext_rng_destroy(shrinkRng);
  rngState = savedRngState;
  printf("ok: prior reset contract (INV-1 residual %.2e)\n", worstDeviation);
}

// The move-validity predicates, called directly. ruleIsValid and the two
// subtree predicates under it decide whether a swap or change proposal is
// SCORED or treated as a no-op, and they are reached only as a side condition
// of those moves: a wrong verdict moves a proposal between the two, which no
// structural check downstream can see. Both bounds of each are pinned here,
// against hand-built trees whose rules are set in place (the predicates read
// rules alone, so the stale member partition is not part of the question).
//
// Ordinal: a descendant splitting on the same variable must fall inside the
// interval its ancestors leave, INCLUSIVE at both ends - the cut one past
// either bound is the rule that would put two members of one leaf on opposite
// sides of an ancestor cut.
//
// Categorical: a rule's direction mask must be a nonempty STRICT subset of the
// categories reaching the node. Equality is the case a gauge that only tests
// containment admits, and it leaves the right child holding every reachable
// category and the left none.
static void testMoveValidityPredicates() {
  const size_t n = 200;
  std::vector<double> xOrdinal(n), xCategorical(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    xOrdinal[i] = static_cast<double>(i) / static_cast<double>(n - 1);
    xCategorical[i] = static_cast<double>(i % 4);
  }

  ColumnStore ordinalStore;
  built(ordinalStore.build(xOrdinal.data(), n, 1, 100));
  check(ordinalStore.numCuts[0] == 100, "the ordinal fixture has 100 cuts");

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(49);
  tree.birth(ordinalStore, 0, rule, y.data(), nullptr);
  int32_t left = tree.at(0).leftChild, right = left + 1;
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.setSplitIndex(30);
  tree.birth(ordinalStore, left, childRule, y.data(), nullptr);
  childRule.setSplitIndex(70);
  tree.birth(ordinalStore, right, childRule, y.data(), nullptr);

  // the root's rule leaves [0, 48] to the left child and [50, 99] to the right
  auto ordinalValid = [&]() {
    return ordinalRuleIsValid(tree, 0, 0, 0,
                              static_cast<int32_t>(ordinalStore.numCuts[0]) - 1);
  };
  tree.at(left).rule.setSplitIndex(48);
  check(ordinalValid(), "a descendant AT its ancestor's upper bound is valid");
  tree.at(left).rule.setSplitIndex(49);
  check(!ordinalValid(), "a descendant one past it is not");
  tree.at(left).rule.setSplitIndex(30);
  tree.at(right).rule.setSplitIndex(50);
  check(ordinalValid(), "a descendant AT its ancestor's lower bound is valid");
  tree.at(right).rule.setSplitIndex(49);
  check(!ordinalValid(), "a descendant one below it is not");
  tree.at(right).rule.setSplitIndex(70);

  // and the same through the dispatcher, which reads the interval itself
  MoveScratch scratch;
  CGMTreePrior prior;
  MoveContext ordinalCtx{ordinalStore, prior, 0.5, 0.1, 0.5, nullptr, 2.0,
                         scratch};
  check(ruleIsValid(ordinalCtx, tree, 0, 0),
        "ruleIsValid accepts the well-formed ordinal subtree");
  tree.at(left).rule.setSplitIndex(60);  // outside [0, 48]
  check(!ruleIsValid(ordinalCtx, tree, 0, 0),
        "ruleIsValid refuses a descendant outside its ancestor interval");

  ColumnKind type = ColumnKind::categorical;
  ColumnStore categoricalStore;
  PredictorSource source =
    densePredictorSource(xCategorical.data(), n, 1, &type);
  built(categoricalStore.build(source, nullptr, 100, false));
  check(categoricalStore.categoryCounts[0] == 4,
        "the categorical fixture has four levels");

  Tree categoricalTree;
  categoricalTree.initialize(indexBuffer.data(), n);
  categoricalTree.computeLeafStats(0, y.data(), nullptr);
  Rule categoricalRule;
  categoricalRule.variableIndex = 0;
  categoricalRule.setCategoryDirections(0x3u);
  categoricalTree.birth(categoricalStore, 0, categoricalRule, y.data(),
                        nullptr);
  std::uint64_t reachable =
    categoricalTree.reachableCategories(categoricalStore, 0, 0);
  check(reachable == 0xfu, "every category reaches the root");

  auto categoricalValid = [&](std::uint64_t directions) {
    categoricalTree.at(0).rule.setCategoryDirections(directions);
    return categoricalSubtreeIsValid(categoricalTree, 0, 0, reachable);
  };
  check(categoricalValid(0x3u),
        "a nonempty strict subset of the reachable categories is valid");
  check(!categoricalValid(0x0u), "an empty direction mask is not");
  check(!categoricalValid(reachable),
        "a mask equal to the reachable set is not: it empties the left child");
  check(!categoricalValid(0x13u),
        "a mask reaching outside the reachable set is not");

  categoricalValid(0x3u);
  MoveContext categoricalCtx{categoricalStore, prior, 0.5, 0.1, 0.5, nullptr,
                            2.0, scratch};
  check(ruleIsValid(categoricalCtx, categoricalTree, 0, 0),
        "ruleIsValid accepts the well-formed categorical subtree");
  categoricalValid(reachable);
  check(!ruleIsValid(categoricalCtx, categoricalTree, 0, 0),
        "ruleIsValid refuses a rule at the reachable set");

  printf("ok: move validity predicates\n");
}

void runMovesTests(ext_rng* rng) {
  testDartUpdate(rng);
  testLeafOfConsistency(rng);
  testPriorResetContract(rng);
  testDartSparsityRecovery(rng);
  testSetPredictorTransaction(rng);
  testSetPredictorForced(rng);
  testDegenerateReCutRoundTrips(rng);
  testCategoricalMissingRoundTrips(rng);
  testUpdatePredictorColumns(rng);
  testPerObservationUpdate(rng);
  testJointPerObservationUpdate();
  testPerObservationMissingCommit(rng);
  testQuantilePredictorUpdate(rng);
  testSetCutPoints(rng);
  testDegenerateRootGuard();
  testMultiChainMutation();
  testLogisticMutation(rng);
  testCategoricalMutation(rng);
  testOrderedFactorMutation(rng);
  testLinearLeafMutation(rng);
  testEmptyLeafVetoCountsWeight();
  testVetoRankUnfreezesStrandedTree();
  testEqualRankOneComparison();
  testMoveValidityPredicates();
}
