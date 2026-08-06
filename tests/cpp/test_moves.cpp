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
  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);
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

  std::vector<xint_t> codesBefore(sampler.data().train.codes);
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());

  // identity swap: new buffer, same values; must accept and preserve fits
  std::vector<double> xCopy(x);
  check(sampler.setPredictor(xCopy.data(), false, false) ==
          PredictorUpdateResult::accepted,
        "identity setPredictor accepted");
  check(sampler.data().train.codes == codesBefore, "identity swap preserves codes");
  check(sampler.chain(0).treeFits() == treeFitsBefore, "identity swap preserves fits");

  // constant predictors empty one side of every split: must reject and
  // roll back completely
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "degenerate setPredictor rejected");
  check(sampler.data().train.codes == codesBefore, "rollback restores codes");
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

  std::vector<ColumnType> types = {ColumnType::categorical,
                                   ColumnType::ordinal};
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

  std::vector<xint_t> codesBefore(sampler.data().train.codes);

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
    columnMatches &= sampler.data().train.codes[i] ==
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
    otherUntouched &= sampler.data().train.codes[i + n] == codesBefore[i + n];
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
        sampler.data().train.codes[i] == sampler.data().codeFor(0, 10.0);
    } else {
      valuesConsistent &=
        sampler.data().train.codes[i] == sampler.data().codeFor(0, identity[i]);
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
      samplerA.impl().data().train.codes[i] == samplerB.impl().data().train.codes[i];
    if (installed[i]) {
      ++numInstalled;
      valuesConsistent &=
        samplerA.impl().data().train.codes[i] == samplerA.impl().data().codeFor(0, 10.0);
    } else {
      valuesConsistent &=
        samplerA.impl().data().train.codes[i] == samplerA.impl().data().codeFor(0, xB[i]);
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
  std::vector<xint_t> codesBefore(sampler.data().train.codes);
  std::vector<double> coarse(n);
  for (size_t i = 0; i < n; ++i) coarse[i] = static_cast<double>(i % 4);
  size_t columnIndex = 0;
  check(sampler.updatePredictor(coarse.data(), &columnIndex, 1, false, true) ==
          PredictorUpdateResult::invalidCutPoints,
        "coarser quantile column update refused");
  check(sampler.data().train.codes == codesBefore,
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
    codesMatch &= sampler.data().train.codes[i] == sampler.data().codeFor(0, x[i]);
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

  std::vector<xint_t> codesBefore(sampler.data().train.codes);
  std::vector<double> fitsBefore[numChains] = {sampler.chain(0).treeFits(),
                                               sampler.chain(1).treeFits()};

  // the transaction spans chains: degenerate predictors roll back everywhere
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "multi-chain degenerate setPredictor rejected");
  check(sampler.data().train.codes == codesBefore,
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
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
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
  std::vector<xint_t> codesBefore(sampler.data().train.codes);
  check(sampler.setPredictor(xBad.data(), false, false) ==
          PredictorUpdateResult::invalidCutPoints,
        "out-of-range category code refused");
  check(sampler.data().train.codes == codesBefore, "refusal mutates nothing");

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
  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);
  check(sampler.numObservations() == n2 && sampler.data().numCuts[0] == 4,
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
  sampler->setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0,
                   nullptr);

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
  sampler->predict(x2.data(), n2, livePredictions.data());
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
  sampler.chain(0).sampleTreesFromPrior();
  sampler.run(1, 0, empty);
  mapsMatch();
  check(allMatch, "leafOf matches after a prior-drawn structure reset");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: leafOf consistency\n");
}

void runMovesTests(ext_rng* rng) {
  testDartUpdate(rng);
  testLeafOfConsistency(rng);
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
  testLinearLeafMutation(rng);
}
