// Component and smoke tests for bartcore against independently coded
// reference math. Exit code 0 on success.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <memory>
#include <vector>

#include <misc/simd.h>
#include <external/random.h>

#include <bartcore/bartcore.hpp>

using namespace bartcore;

// misc.a/external.a reach into the R runtime for output; stub it.
extern "C" {
void Rprintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}
void R_FlushConsole(void) { fflush(stdout); }
}

static int failures = 0;

static void check(bool condition, const char* what) {
  if (!condition) {
    ++failures;
    printf("FAIL: %s\n", what);
  }
}

static void checkNear(double actual, double expected, double tolerance,
                      const char* what) {
  if (!(std::fabs(actual - expected) <= tolerance)) {
    ++failures;
    printf("FAIL: %s (actual %.15g, expected %.15g)\n", what, actual, expected);
  }
}

static uint64_t rngState = 0x9E3779B97F4A7C15ull;
static double runif01() {
  rngState ^= rngState << 13;
  rngState ^= rngState >> 7;
  rngState ^= rngState << 17;
  return (double) (rngState >> 11) * 0x1.0p-53;
}

static void testColumnStoreCodes() {
  const size_t n = 500, p = 3;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  ColumnStore store;
  store.build(x.data(), n, p, 100);

  // reference: linear scan against uniformly spaced cuts
  for (size_t j = 0; j < p; ++j) {
    const double* column = x.data() + j * n;
    double xMin = column[0], xMax = column[0];
    for (size_t i = 1; i < n; ++i) {
      xMin = std::min(xMin, column[i]);
      xMax = std::max(xMax, column[i]);
    }
    double increment = (xMax - xMin) / 101.0;
    for (size_t i = 0; i < n; ++i) {
      uint32_t k = 0;
      while (k < 100 && column[i] > xMin + (double) (k + 1) * increment) ++k;
      if (store.codes[i + j * n] != (xint_t) k) {
        check(false, "column store code mismatch");
        return;
      }
    }
  }
  check(true, "");
  printf("ok: column store codes\n");
}

static void testColumnStoreView() {
  const size_t n = 120, p = 3;
  std::vector<double> x(n * p);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    x[i + 2 * n] = (double) (i % 4);  // categorical codes 0..3
  }
  // rows 0 and 1 carry column 0's extremes; the subset below excludes them,
  // so a store built over only the subset would bin that column differently
  x[0] = 0.0;
  x[1] = 1.0;

  std::vector<ColumnType> types = {ColumnType::ordinal, ColumnType::ordinal,
                                   ColumnType::categorical};
  ColumnStore parent;
  parent.build(x.data(), n, p, 25, false, types.data());

  std::vector<size_t> rows, testRows;
  for (size_t i = 2; i < n; i += 2) rows.push_back(i);
  testRows.push_back(0);  // the extremes land in the test rows
  testRows.push_back(1);
  for (size_t i = 3; i < n; i += 2) testRows.push_back(i);

  ColumnStore view;
  view.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                       testRows.size());

  check(view.x == nullptr && view.x_test == nullptr,
        "view holds no raw values");
  check(view.numObservations == rows.size() &&
        view.numTestObservations == testRows.size() &&
        view.numPredictors == p, "view dimensions");
  check(view.types == parent.types && view.numCuts == parent.numCuts &&
        view.cutPoints == parent.cutPoints &&
        view.maxNumCuts == parent.maxNumCuts,
        "view copies the parent cut grid");

  bool codesMatch = true;
  for (size_t j = 0; j < p && codesMatch; ++j)
    for (size_t i = 0; i < rows.size() && codesMatch; ++i)
      codesMatch =
        view.codes[i + j * rows.size()] == parent.codes[rows[i] + j * n];
  check(codesMatch, "view gathers subset codes");

  bool testCodesMatch = true;
  for (size_t i = 0; i < testRows.size() && testCodesMatch; ++i)
    for (size_t j = 0; j < p && testCodesMatch; ++j)
      testCodesMatch =
        view.testCodes[i * p + j] == parent.codes[testRows[i] + j * n];
  check(testCodesMatch, "view gathers test codes from parent rows");

  // demonstrate the property matters: a store built over the subset's raw
  // values bins column 0 onto a different grid than the view keeps
  std::vector<double> subsetX(rows.size() * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < rows.size(); ++i)
      subsetX[i + j * rows.size()] = x[rows[i] + j * n];
  ColumnStore rebuilt;
  rebuilt.build(subsetX.data(), rows.size(), p, 25, false, types.data());
  bool anyDiffer = false;
  for (size_t i = 0; i < rows.size() && !anyDiffer; ++i)
    anyDiffer = rebuilt.codes[i] != view.codes[i];
  check(anyDiffer, "subset-built store bins differently than the view");

  printf("ok: column store view\n");
}

static void testViewSamplerMatchesFull() {
  const size_t n = 300, p = 4;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = std::sin(3.0 * x[i]) + x[i + n] + 0.5 * runif01();

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 1234) != 0 ||
      ext_rng_setSeed(rngB, 1234) != 0) {
    check(false, "view sampler: rng creation");
    return;
  }

  SamplerOptions options;
  options.numTrees = 25;

  ClassicSampler full(x.data(), y.data(), n, p, nullptr, nullptr,
                      ResponseFamily::gaussian, 1.0, 3.0,
                      0.37804942330213542, options, &rngA);

  ColumnStore parent;
  parent.build(x.data(), n, p, options.maxNumCuts);
  std::vector<size_t> rows(n);
  for (size_t i = 0; i < n; ++i) rows[i] = i;
  ColumnStore store;
  store.buildFromParent(parent, rows.data(), n, nullptr, 0);
  ClassicSampler view(std::move(store), y.data(), nullptr, nullptr,
                      ResponseFamily::gaussian, 1.0, 3.0,
                      0.37804942330213542, options, &rngB);

  const size_t numBurnIn = 30, numSamples = 40;
  std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
  std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = fitsA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = fitsB.data();
  full.run(numBurnIn, numSamples, resultsA);
  view.run(numBurnIn, numSamples, resultsB);

  bool identical = true;
  for (size_t s = 0; s < numSamples && identical; ++s)
    identical = sigmaA[s] == sigmaB[s];
  for (size_t i = 0; i < n * numSamples && identical; ++i)
    identical = fitsA[i] == fitsB[i];
  check(identical, "full-rows view sampler bitwise-matches the classic path");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: view sampler matches full\n");
}

static void testIntegratedLikelihood() {
  ConstantGaussianLeaf leaf{0.5 / std::sqrt(200.0)};
  double k = 2.0, sigmaSq = 0.01;
  double average = 0.031, numEff = 47.0, variance = 0.0042;
  size_t numObs = 47;

  // independent transcription of the CGM marginal likelihood
  double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
  double posteriorPrecision = numEff / sigmaSq;
  double expected = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision))
    - 0.5 * (variance / sigmaSq) * (double) (numObs - 1)
    - 0.5 * ((priorPrecision * average) * (posteriorPrecision * average)) /
        (priorPrecision + posteriorPrecision);

  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, average, numEff, variance, numObs),
            expected, 1e-13, "integrated likelihood formula");
  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, average, numEff, variance, 0),
            0.0, 0.0, "empty leaf contributes zero");
  printf("ok: integrated likelihood\n");
}

static void testPosteriorDraw(ext_rng* rng) {
  ConstantGaussianLeaf leaf{0.5 / std::sqrt(50.0)};
  double k = 2.0, sigmaSq = 0.02, average = 0.12, numEff = 30.0;

  double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
  double posteriorPrecision = numEff / sigmaSq;
  double expectedMean = posteriorPrecision * average / (priorPrecision + posteriorPrecision);
  double expectedSd = 1.0 / std::sqrt(priorPrecision + posteriorPrecision);

  const int numDraws = 200000;
  double sum = 0.0, sumSq = 0.0;
  for (int i = 0; i < numDraws; ++i) {
    double draw = leaf.drawFromPosterior(rng, k, average, numEff, sigmaSq);
    sum += draw;
    sumSq += draw * draw;
  }
  double mean = sum / numDraws;
  double sd = std::sqrt(sumSq / numDraws - mean * mean);
  checkNear(mean, expectedMean, 5.0 * expectedSd / std::sqrt((double) numDraws),
            "posterior draw mean");
  checkNear(sd, expectedSd, 0.01 * expectedSd, "posterior draw sd");
  printf("ok: posterior draws\n");
}

// Build x with a known partition structure and verify splitting mechanics.
static void testTreeMechanics() {
  const size_t n = 100;
  std::vector<double> x(n);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = (double) i / (double) (n - 1);
    y[i] = x[i] < 0.5 ? 1.0 : 3.0;
  }

  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  checkNear(tree.at(0).average, 2.0, 1e-12, "root average");

  // split at the median cut
  Rule rule;
  rule.variableIndex = 0;
  rule.splitIndex = 49;
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  check(left.numObservations() + right.numObservations() == n,
        "birth partitions all observations");

  // verify against a direct count with the store's own codes
  size_t expectedLeft = 0;
  for (size_t i = 0; i < n; ++i)
    if (store.codes[i] <= 49) ++expectedLeft;
  check(left.numObservations() == expectedLeft, "partition count matches codes");

  double leftSum = 0.0, rightSum = 0.0;
  for (size_t i = 0; i < n; ++i)
    (store.codes[i] <= 49 ? leftSum : rightSum) += y[i];
  checkNear(left.average, leftSum / (double) left.numObservations(), 1e-12,
            "left child average");
  checkNear(right.average, rightSum / (double) right.numObservations(), 1e-12,
            "right child average");

  // split interval of the left child is bounded by the parent's rule
  int32_t lo, hi;
  tree.splitInterval(store, tree.at(0).leftChild, 0, &lo, &hi);
  check(lo == 0 && hi == 48, "left child split interval");
  tree.splitInterval(store, tree.at(0).leftChild + 1, 0, &lo, &hi);
  check(lo == 50 && hi == 99, "right child split interval");

  // orphanChildren merges with effective-observation weighting
  double expectedMerged =
    left.average * (left.numEffectiveObservations /
                    (left.numEffectiveObservations + right.numEffectiveObservations)) +
    right.average * (right.numEffectiveObservations /
                     (left.numEffectiveObservations + right.numEffectiveObservations));
  tree.orphanChildren(0);
  checkNear(tree.at(0).average, expectedMerged, 1e-12, "orphan merge average");
  check(tree.at(0).isBottom(), "orphan clears children");

  printf("ok: tree mechanics\n");
}

static void testTreePriorMath() {
  const size_t n = 64;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = (double) i / (double) (n - 1);

  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);

  CGMTreePrior prior;  // base 0.95, power 2

  checkNear(prior.growthProbability(tree, store, 0), 0.95, 1e-15,
            "root growth probability");

  Rule rule;
  rule.variableIndex = 0;
  rule.splitIndex = 49;
  tree.birth(store, 0, rule, y.data(), nullptr);
  int32_t left = tree.at(0).leftChild;

  checkNear(prior.growthProbability(tree, store, left), 0.95 / 4.0, 1e-15,
            "depth-1 growth probability");

  // single-split tree: log p = log(g0) + log(1/p) + log(1/numCuts)
  //                          + 2 log(1 - g1)
  double expected = std::log(0.95) + std::log(1.0) - std::log(100.0) +
                    2.0 * std::log(1.0 - 0.95 / 4.0);
  checkNear(prior.treeLogProbability(tree, store), expected, 1e-12,
            "tree log probability");

  printf("ok: tree prior math\n");
}

static void testEndToEndGaussian(ext_rng* rng) {
  const size_t n = 500, p = 5;
  std::vector<double> x(n * p), f(n), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    f[i] = std::sin(3.0 * x[i]) + 2.0 * x[i + n] * x[i + 2 * n];
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = f[i] + 0.3 * normal;
  }

  double yMean = 0.0, ySd = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= (double) n;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / (double) (n - 1));

  SamplerOptions options;
  options.numTrees = 75;
  ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
                         ySd, 3.0, 0.37804942330213542 /* qchisq(0.1, 3)/3 */,
                         options, &rng);

  const size_t numBurnIn = 200, numSamples = 300;
  std::vector<double> sigmaDraws(numSamples);
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.sigma = sigmaDraws.data();
  results.trainingFits = trainingFits.data();
  sampler.run(numBurnIn, numSamples, results);

  // posterior mean fit should track f substantially better than the mean
  std::vector<double> posteriorMean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i)
      posteriorMean[i] += trainingFits[i + s * n] / (double) numSamples;

  double sseFit = 0.0, sseMean = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sseFit += (posteriorMean[i] - f[i]) * (posteriorMean[i] - f[i]);
    sseMean += (yMean - f[i]) * (yMean - f[i]);
  }
  check(sseFit < 0.2 * sseMean, "end-to-end: fit explains most signal");

  double sigmaPosteriorMean = 0.0;
  for (double s : sigmaDraws) sigmaPosteriorMean += s / (double) numSamples;
  check(sigmaPosteriorMean > 0.2 && sigmaPosteriorMean < 0.45,
        "end-to-end: sigma near truth (0.3)");

  printf("ok: end-to-end gaussian (sigma posterior mean %.3f, sse ratio %.3f)\n",
         sigmaPosteriorMean, sseFit / sseMean);
}

static void testEndToEndProbit(ext_rng* rng) {
  const size_t n = 800, p = 3;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double eta = 2.0 * (x[i] - 0.5);
    double probability = 0.5 * std::erfc(-eta / std::sqrt(2.0));
    y[i] = runif01() < probability ? 1.0 : 0.0;
  }

  SamplerOptions options;
  options.numTrees = 50;
  options.nodeScale = 3.0;
  ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::probit,
                         1.0, 3.0, 1.0, options, &rng);

  const size_t numBurnIn = 150, numSamples = 200;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(numBurnIn, numSamples, results);

  // monotone signal in x1: mean latent fit should increase across quartiles
  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < numSamples; ++s) {
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < 0.25) { lowSum += trainingFits[i + s * n]; ++lowCount; }
      if (x[i] > 0.75) { highSum += trainingFits[i + s * n]; ++highCount; }
    }
  }
  check(highSum / (double) highCount > lowSum / (double) lowCount + 0.3,
        "end-to-end probit: monotone signal recovered");

  printf("ok: end-to-end probit\n");
}

static void testPolyaGamma(ext_rng* rng) {
  // moments of PG(1, c): mean tanh(c / 2) / (2 c), variance
  // (sinh(c) - c) sech^2(c / 2) / (4 c^3); the limits at 0 are 1/4 and 1/24.
  // psi values cover both proposal branches, the sign fold, and the tail.
  const double psis[] = {0.0, 1.5, -3.0, 8.0};
  const int numDraws = 50000;

  for (double psi : psis) {
    double c = std::fabs(psi);
    double expectedMean = c == 0.0 ? 0.25 : std::tanh(0.5 * c) / (2.0 * c);
    double expectedVariance = c == 0.0
      ? 1.0 / 24.0
      : (std::sinh(c) - c) / (4.0 * c * c * c * std::cosh(0.5 * c) *
                              std::cosh(0.5 * c));

    double sum = 0.0, sumSq = 0.0;
    for (int i = 0; i < numDraws; ++i) {
      double draw = ext_rng_simulatePolyaGamma(rng, psi);
      if (!(draw > 0.0)) {
        check(false, "polya-gamma draw not positive");
        return;
      }
      sum += draw;
      sumSq += draw * draw;
    }
    double mean = sum / numDraws;
    double variance = sumSq / numDraws - mean * mean;

    double meanSe = std::sqrt(expectedVariance / (double) numDraws);
    checkNear(mean, expectedMean, 5.0 * meanSe, "polya-gamma mean");
    checkNear(variance, expectedVariance, 0.05 * expectedVariance,
              "polya-gamma variance");
  }

  printf("ok: polya-gamma sampler\n");
}

static void testEndToEndLogistic(ext_rng* rng) {
  const size_t n = 800, p = 3;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double eta = 4.0 * (x[i] - 0.5);
    double probability = 1.0 / (1.0 + std::exp(-eta));
    y[i] = runif01() < probability ? 1.0 : 0.0;
  }

  SamplerOptions options;
  options.numTrees = 50;
  options.nodeScale = 3.0;
  ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                         ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                         &rng);

  const size_t numBurnIn = 150, numSamples = 200;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(numBurnIn, numSamples, results);

  // monotone signal in x1 on the log-odds scale
  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < numSamples; ++s) {
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < 0.25) { lowSum += trainingFits[i + s * n]; ++lowCount; }
      if (x[i] > 0.75) { highSum += trainingFits[i + s * n]; ++highCount; }
    }
  }
  double lowFit = lowSum / (double) lowCount;
  double highFit = highSum / (double) highCount;
  check(highFit > lowFit + 0.5, "end-to-end logistic: monotone signal recovered");

  // the fits should be calibrated to the quartiles' empirical log odds
  double lowRate = 0.0, highRate = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (x[i] < 0.25) lowRate += y[i] / ((double) lowCount / numSamples);
    if (x[i] > 0.75) highRate += y[i] / ((double) highCount / numSamples);
  }
  double lowProbability = 1.0 / (1.0 + std::exp(-lowFit));
  double highProbability = 1.0 / (1.0 + std::exp(-highFit));
  checkNear(lowProbability, lowRate, 0.1, "logistic low-quartile calibration");
  checkNear(highProbability, highRate, 0.1, "logistic high-quartile calibration");

  // omega latents are exposed and positive
  const double* omega = sampler.latents(0);
  bool omegaValid = omega != nullptr;
  for (size_t i = 0; i < n && omegaValid; ++i)
    omegaValid = omega[i] > 0.0 && std::isfinite(omega[i]);
  check(omegaValid, "logistic latents expose positive omega draws");

  printf("ok: end-to-end logistic\n");
}

static void makeMutationData(std::vector<double>& x, std::vector<double>& y,
                             size_t n);

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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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
    ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
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

static void testChiKHyperprior(ext_rng* rng) {
  ChiKHyperprior prior;
  prior.degreesOfFreedom = 1.25;

  double sumSquaredParams = 3.7, numLeaves = 320.0;
  double leafScale = 0.5 / std::sqrt(200.0);

  // k = sqrt(X), X ~ Gamma(shape, 1/rate): E[k] = rate^-1/2 G(a+1/2)/G(a)
  double shape = 0.5 * (numLeaves + 2.0 * prior.degreesOfFreedom - 1.0);
  double rate = 0.5 * sumSquaredParams / (leafScale * leafScale);
  auto expectedMean = [shape](double r) {
    return std::sqrt(1.0 / r) * std::exp(std::lgamma(shape + 0.5) - std::lgamma(shape));
  };
  auto drawMean = [&](int numDraws) {
    double sum = 0.0;
    for (int i = 0; i < numDraws; ++i)
      sum += prior.draw(rng, sumSquaredParams, numLeaves, leafScale);
    return sum / numDraws;
  };

  checkNear(drawMean(100000), expectedMean(rate), 0.005 * expectedMean(rate),
            "chi-k posterior draw mean, infinite prior scale");

  prior.scale = 5.0;  // finite prior scale enters the rate
  double rateFinite = rate + 0.5 / (prior.scale * prior.scale);
  checkNear(drawMean(100000), expectedMean(rateFinite),
            0.005 * expectedMean(rateFinite),
            "chi-k posterior draw mean, finite prior scale");

  printf("ok: chi-k hyperprior\n");
}

static void testColumnStoreMutation() {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();

  ColumnStore store;
  store.build(x.data(), n, p, 100);
  std::vector<xint_t> originalCodes(store.codes);
  std::vector<double> originalCuts0(store.cutPoints[0]);

  // column overwrite with cut refresh: column 0 codes untouched, column 1
  // re-quantized against cuts spanning the new range
  std::vector<double> newColumn(n);
  for (double& v : newColumn) v = 2.0 + 3.0 * runif01();
  size_t columnIndex = 1;
  store.setColumns(newColumn.data(), &columnIndex, 1, true);

  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i) {
    codesMatch &= store.codes[i] == originalCodes[i];
    codesMatch &= x[i + n] == newColumn[i];
    codesMatch &= store.codes[i + n] == store.codeFor(1, newColumn[i]);
  }
  check(codesMatch, "setColumns re-quantizes only the target column");
  check(store.cutPoints[1].front() > 2.0 && store.cutPoints[1].back() < 5.0,
        "setColumns recomputes cuts over the new range");
  check(store.cutPoints[0] == originalCuts0, "setColumns leaves other cuts");

  // single-cell overwrite against existing cuts
  xint_t before = store.codes[7];
  store.setCell(7, 0, x[8]);
  check(store.codes[7] == originalCodes[8] && x[7] == x[8],
        "setCell re-quantizes one cell");
  check(before == originalCodes[7], "");  // silence unused warning

  // whole-matrix pointer swap without cut refresh: quantized on old cuts
  std::vector<double> x2(n * p);
  for (double& v : x2) v = runif01();
  store.setPredictors(x2.data(), false);
  check(store.x == x2.data(), "setPredictors swaps the pointer");
  codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= store.codes[i] == store.codeFor(0, x2[i]);
  check(codesMatch, "setPredictors re-quantizes against existing cuts");

  printf("ok: column store mutation\n");
}

// A burned-in sampler for mutation tests: strong signal in both columns so
// trees certainly split.
static std::unique_ptr<ClassicSampler> makeBurnedInSampler(
  std::vector<double>& x, std::vector<double>& y, size_t n, ext_rng* rng) {
  SamplerOptions options;
  options.numTrees = 25;
  auto sampler = std::make_unique<ClassicSampler>(
    x.data(), y.data(), n, size_t(2), nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
    0.37804942330213542, options, &rng);
  Results empty;
  sampler->run(100, 0, empty);
  return sampler;
}

static void makeMutationData(std::vector<double>& x, std::vector<double>& y,
                             size_t n) {
  x.resize(n * 2);
  y.resize(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n] + 0.2 * normal;
  }
}

static void testSetPredictorTransaction(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  std::vector<xint_t> codesBefore(sampler.data().codes);
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());

  // identity swap: new buffer, same values; must accept and preserve fits
  std::vector<double> xCopy(x);
  check(sampler.setPredictor(xCopy.data(), false, false) ==
          PredictorUpdateResult::accepted,
        "identity setPredictor accepted");
  check(sampler.data().x == xCopy.data(), "accepted swap installs pointer");
  check(sampler.data().codes == codesBefore, "identity swap preserves codes");
  check(sampler.chain(0).treeFits() == treeFitsBefore, "identity swap preserves fits");

  // constant predictors empty one side of every split: must reject and
  // roll back completely
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "degenerate setPredictor rejected");
  check(sampler.data().x == xCopy.data(), "rejected swap keeps old pointer");
  check(sampler.data().codes == codesBefore, "rollback restores codes");
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
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

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

static void testUpdatePredictorColumns(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  std::vector<xint_t> codesBefore(sampler.data().codes);

  // jittered column: usually accepted; either way the sampler must stay
  // consistent, and rejection must restore the column bitwise
  std::vector<double> column0Before(x.begin(), x.begin() + n);
  std::vector<double> jittered(column0Before);
  for (double& v : jittered) v += 0.001 * (runif01() - 0.5);
  size_t columnIndex = 0;
  bool accepted = sampler.updatePredictor(jittered.data(), &columnIndex, 1,
                                          false, false) ==
                  PredictorUpdateResult::accepted;
  bool columnMatches = true;
  for (size_t i = 0; i < n; ++i) {
    double expected = accepted ? jittered[i] : column0Before[i];
    columnMatches &= sampler.data().x[i] == expected;
    columnMatches &= sampler.data().codes[i] ==
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
    otherUntouched &= sampler.data().codes[i + n] == codesBefore[i + n];
  check(otherUntouched, "rejected column update leaves other columns");

  printf("ok: updatePredictor columns\n");
}

static void testPerObservationUpdate(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

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
    if (installed[i]) {
      ++numInstalled;
      valuesConsistent &= sampler.data().x[i] == 10.0;
    } else {
      valuesConsistent &= sampler.data().x[i] == identity[i];
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

  // all-or-none: the shared column stays identical across samplers
  bool columnsAgree = true, valuesConsistent = true;
  size_t numInstalled = 0;
  for (size_t i = 0; i < n; ++i) {
    columnsAgree &= samplerA.impl().data().x[i] == samplerB.impl().data().x[i];
    if (installed[i]) {
      ++numInstalled;
      valuesConsistent &= samplerA.impl().data().x[i] == 10.0;
    } else {
      valuesConsistent &= samplerA.impl().data().x[i] == xB[i];
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

static void testQuantileCutPoints() {
  const size_t n = 200;
  // column 0: continuous with more uniques than cuts; column 1: 10 discrete
  // levels, so quantile mode induces exactly 9 midpoint cuts
  std::vector<double> x(n * 2);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = static_cast<double>(i % 10);
  }

  ColumnStore store;
  store.build(x.data(), n, 2, 20, true);

  check(store.numCuts[0] == 20, "quantile cuts capped at maxNumCuts");
  check(store.numCuts[1] == 9, "few uniques induce numUnique - 1 cuts");
  bool discreteCutsMatch = true;
  for (std::uint32_t k = 0; k < 9; ++k)
    discreteCutsMatch &= store.cutPoints[1][k] == static_cast<double>(k) + 0.5;
  check(discreteCutsMatch, "discrete quantile cuts are unique-value midpoints");
  bool discreteCodesMatch = true;
  for (size_t i = 0; i < n; ++i)
    discreteCodesMatch &= store.codes[i + n] == static_cast<xint_t>(i % 10);
  check(discreteCodesMatch, "discrete quantile codes are value ranks");

  // continuous column: reference the thinning directly
  std::vector<double> sorted(x.begin(), x.begin() + n);
  std::sort(sorted.begin(), sorted.end());
  size_t step = n / 20, offset = step / 2;
  bool continuousCutsMatch = true;
  for (std::uint32_t k = 0; k < 20; ++k) {
    size_t index = std::min(static_cast<size_t>(k) * step + offset, n - 2);
    continuousCutsMatch &=
      store.cutPoints[0][k] == 0.5 * (sorted[index] + sorted[index + 1]);
  }
  check(continuousCutsMatch, "continuous quantile cuts thin sorted uniques");

  // refresh feasibility: fewer uniques than existing cuts is invalid
  std::vector<double> coarse(n);
  for (size_t i = 0; i < n; ++i) coarse[i] = static_cast<double>(i % 4);
  check(!store.cutsWouldRemainValid(1, coarse.data()),
        "coarser column fails the quantile feasibility check");
  std::vector<double> finer(n);
  for (size_t i = 0; i < n; ++i) finer[i] = static_cast<double>(i % 25);
  check(store.cutsWouldRemainValid(1, finer.data()),
        "finer column passes the quantile feasibility check");

  printf("ok: quantile cut points\n");
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(100, 0, empty);

  check(sampler.data().numCuts[0] == 9, "sampler builds quantile cuts");

  // a coarser column with updateCutPoints must be refused without mutating
  std::vector<xint_t> codesBefore(sampler.data().codes);
  std::vector<double> coarse(n);
  for (size_t i = 0; i < n; ++i) coarse[i] = static_cast<double>(i % 4);
  size_t columnIndex = 0;
  check(sampler.updatePredictor(coarse.data(), &columnIndex, 1, false, true) ==
          PredictorUpdateResult::invalidCutPoints,
        "coarser quantile column update refused");
  check(sampler.data().codes == codesBefore,
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
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  // shrink column 0 to three cuts: codes re-quantize, out-of-range splits
  // collapse, and the fit identity holds
  double newCuts[] = {0.25, 0.5, 0.75};
  std::uint32_t numNewCuts = 3;
  const double* cutsByColumn[] = {newCuts};
  size_t columnIndex = 0;
  sampler.setCutPoints(cutsByColumn, &numNewCuts, &columnIndex, 1);

  check(sampler.data().numCuts[0] == 3, "setCutPoints installs the new count");
  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= sampler.data().codes[i] == sampler.data().codeFor(0, x[i]);
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

static void testMultiChain() {
  const size_t n = 300, numChains = 3, numSamples = 50;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  auto makeRngs = [](std::vector<ext_rng*>& rngs) {
    for (size_t c = 0; c < numChains; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], 1000 + static_cast<uint_least32_t>(c));
    }
  };

  auto runSampler = [&](size_t numThreads, std::vector<double>& sigma,
                        std::vector<double>& fits) {
    std::vector<ext_rng*> rngs(numChains);
    makeRngs(rngs);

    SamplerOptions options;
    options.numTrees = 25;
    options.numChains = numChains;
    options.numThreads = numThreads;
    ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                           1.0, 3.0, 0.37804942330213542, options,
                           rngs.data());

    sigma.assign(numSamples * numChains, 0.0);
    fits.assign(n * numSamples * numChains, 0.0);
    Results results;
    results.sigma = sigma.data();
    results.trainingFits = fits.data();
    sampler.run(100, numSamples, results);

    for (size_t c = numChains; c > 0; --c) ext_rng_destroy(rngs[c - 1]);
  };

  std::vector<double> sigmaSerial, fitsSerial, sigmaThreaded, fitsThreaded;
  runSampler(1, sigmaSerial, fitsSerial);
  runSampler(numChains, sigmaThreaded, fitsThreaded);

  // results depend only on the chain rngs, never on the thread count
  check(sigmaSerial == sigmaThreaded, "thread count does not change sigma");
  check(fitsSerial == fitsThreaded, "thread count does not change fits");

  bool sigmaSane = true;
  for (double s : sigmaSerial) sigmaSane &= std::isfinite(s) && s > 0.0;
  check(sigmaSane, "every chain draws finite positive sigma");

  // chains use distinct rngs, so their draws differ
  bool chainsIdentical = true;
  for (size_t s = 0; s < numSamples; ++s)
    chainsIdentical &= sigmaSerial[s] == sigmaSerial[s + numSamples];
  check(!chainsIdentical, "chains draw distinct samples");

  // every chain's posterior mean fit explains the signal
  for (size_t c = 0; c < numChains; ++c) {
    double sse = 0.0, sseMean = 0.0;
    double yMean = 0.0;
    for (size_t i = 0; i < n; ++i) yMean += y[i] / static_cast<double>(n);
    for (size_t i = 0; i < n; ++i) {
      double fitMean = 0.0;
      for (size_t s = 0; s < numSamples; ++s)
        fitMean += fitsSerial[i + (c * numSamples + s) * n] /
                   static_cast<double>(numSamples);
      sse += (fitMean - y[i]) * (fitMean - y[i]);
      sseMean += (yMean - y[i]) * (yMean - y[i]);
    }
    if (sse >= 0.5 * sseMean) {
      check(false, "multi-chain fit explains most signal");
      break;
    }
  }

  printf("ok: multi-chain run\n");
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, rngs.data());
  Results empty;
  sampler.run(100, 0, empty);

  std::vector<xint_t> codesBefore(sampler.data().codes);
  std::vector<double> fitsBefore[numChains] = {sampler.chain(0).treeFits(),
                                               sampler.chain(1).treeFits()};

  // the transaction spans chains: degenerate predictors roll back everywhere
  std::vector<double> xConstant(n * 2, 0.5);
  check(sampler.setPredictor(xConstant.data(), false, false) ==
          PredictorUpdateResult::rolledBack,
        "multi-chain degenerate setPredictor rejected");
  check(sampler.data().codes == codesBefore,
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

static void testMapOldCutPointsOntoNew() {
  // quantile grids over 1..8 give cuts {1.5, ..., 7.5}; hand-check the remap
  const size_t n = 8;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i + 1);

  ColumnStore store;
  store.build(x.data(), n, 1, 7, true);
  std::vector<std::vector<double>> oldCuts(store.cutPoints);

  // root splits at index 3 (cut 4.5), its left child at index 1 (cut 2.5)
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rootRule;  rootRule.variableIndex = 0;  rootRule.splitIndex = 3;
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;  leftRule.variableIndex = 0;  leftRule.splitIndex = 1;
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);
  int32_t leftChild = tree.at(0).leftChild;

  // new values 2..16 by twos: cuts {3, 5, ..., 15}; 4.5 is nearest 5 (index
  // 1), and 2.5 clamps into the left child's shrunken interval [0, 1)
  std::vector<double> x2(n);
  for (size_t i = 0; i < n; ++i) x2[i] = 2.0 * static_cast<double>(i + 1);
  store.setData(x2.data(), n);

  std::vector<double> params(tree.nodes.size(), 0.0);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex == 1, "root split remaps to nearest cut");
  check(tree.at(leftChild).rule.splitIndex == 0,
        "child split clamps into the ancestor-constrained interval");

  // shift the grid entirely above the old cuts: the root clamps to index 0,
  // leaving the left child no interval, so it collapses with plain-mean param
  std::vector<double> x3(n);
  for (size_t i = 0; i < n; ++i) x3[i] = 20.0 + 2.0 * static_cast<double>(i);
  oldCuts = store.cutPoints;
  int32_t oldRootIndex = tree.at(0).rule.splitIndex;
  check(oldRootIndex == 1, "");  // silence unused in release
  store.setData(x3.data(), n);

  params.assign(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(leftChild, bottoms);
  for (size_t k = 0; k < bottoms.size(); ++k)
    params[static_cast<size_t>(bottoms[k])] = static_cast<double>(k + 1);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex == 0, "root split clamps to the low end");
  check(tree.at(leftChild).isBottom(), "interval-starved subtree collapses");
  checkNear(params[static_cast<size_t>(leftChild)], 1.5, 1e-15,
            "collapsed subtree takes the plain mean of its leaf parameters");

  printf("ok: mapOldCutPointsOntoNew\n");
}

static void testSetData(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  double sigmaBefore = sampler.sigma(0);
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());
  std::vector<xint_t> codesBefore(sampler.data().codes);

  // identity replacement: same values in new buffers; the rebuilt cuts equal
  // the old ones, every split remaps onto itself, and fits are recovered
  // exactly
  std::vector<double> xCopy(x), yCopy(y);
  sampler.setData(xCopy.data(), yCopy.data(), n, nullptr, nullptr, nullptr, 0);
  check(sampler.data().x == xCopy.data(), "setData installs the new pointer");
  check(sampler.data().codes == codesBefore, "identity setData preserves codes");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "identity setData preserves fits");
  check(sampler.sigma(0) == sigmaBefore, "identity setData preserves sigma");

  // a linear response transform doubles the data range and so halves the
  // internal sigma; on the original scale it must not move
  std::vector<double> yScaled(n);
  for (size_t i = 0; i < n; ++i) yScaled[i] = 2.0 * y[i] + 3.0;
  sampler.setData(xCopy.data(), yScaled.data(), n, nullptr, nullptr, nullptr, 0);
  checkNear(sampler.sigma(0), sigmaBefore, 1e-12 * sigmaBefore,
            "setData preserves sigma on the original scale");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "setData leaves no empty leaves");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after setData");

  printf("ok: setData\n");
}

static void testSetDataResize(ext_rng* rng) {
  const size_t n = 200, n2 = 320, nTest = 50, nTest2 = 30;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();
  sampler.setTestPredictors(xTest.data(), nTest);

  // replace everything at once: more observations, a smaller test set, and a
  // shifted predictor range so cut points genuinely move
  std::vector<double> x2, y2;
  makeMutationData(x2, y2, n2);
  for (double& v : x2) v = 2.0 * v + 1.0;
  std::vector<double> xTest2(nTest2 * 2);
  for (double& v : xTest2) v = 2.0 * runif01() + 1.0;

  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, xTest2.data(),
                  nTest2);
  check(sampler.numObservations() == n2, "setData resizes observations");
  check(sampler.numTestObservations() == nTest2, "setData resizes test set");

  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "resized setData leaves no empty leaves");

  bool fitIdentity = true;
  for (size_t i = 0; i < n2 && fitIdentity; i += 41) {
    double total = 0.0;
    for (size_t t = 0; t < 25; ++t)
      total += sampler.chain(0).treeFits()[t * n2 + i];
    fitIdentity = std::fabs(total - sampler.chain(0).totalFits()[i]) < 1e-10;
  }
  check(fitIdentity, "resized setData keeps the fit identity");

  std::vector<double> sigmaDraws(5), testFits(nTest2 * 5);
  Results results;
  results.sigma = sigmaDraws.data();
  results.testFits = testFits.data();
  sampler.run(0, 5, results);
  bool finite = true;
  for (double s : sigmaDraws) finite &= std::isfinite(s) && s > 0.0;
  for (double f : testFits) finite &= std::isfinite(f);
  check(finite, "sampler runs after resized setData");

  // dropping the test set entirely
  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);
  check(sampler.numTestObservations() == 0, "setData clears the test set");

  printf("ok: setData resize\n");
}

static void testSetDataQuantileShrink(ext_rng* rng) {
  const size_t n = 200;
  // 10 discrete levels induce 9 quantile cuts
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 10);
    x[i + n] = runif01();
    y[i] = 0.5 * x[i] + 2.0 * x[i + n] + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.useQuantiles = true;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(100, 0, empty);
  check(sampler.data().numCuts[0] == 9, "quantile sampler starts with 9 cuts");

  // a coarser column is refused by updatePredictor but legal here: the cut
  // count shrinks and out-of-range splits remap or collapse
  std::vector<double> x2(x);
  for (size_t i = 0; i < n; ++i) x2[i] = static_cast<double>(i % 4);
  sampler.setData(x2.data(), y.data(), n, nullptr, nullptr, nullptr, 0);
  check(sampler.data().numCuts[0] == 3, "setData shrinks quantile cut counts");

  bool occupied = true, codesInRange = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  for (size_t i = 0; i < n; ++i)
    codesInRange &= sampler.data().codes[i] <= 3;
  check(occupied, "quantile shrink leaves no empty leaves");
  check(codesInRange, "quantile shrink re-quantizes codes");

  std::vector<double> sigmaDraws(5);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "sampler runs after quantile shrink");

  printf("ok: setData quantile shrink\n");
}

static void testSetDataProbit(ext_rng* rng) {
  const size_t n = 200, n2 = 150;
  std::vector<double> x, yContinuous;
  makeMutationData(x, yContinuous, n);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i) y[i] = yContinuous[i] > 0.0 ? 1.0 : 0.0;

  SamplerOptions options;
  options.numTrees = 25;
  options.nodeScale = 3.0;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::probit,
                         1.0, 3.0, 0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(50, 0, empty);

  std::vector<double> x2, y2Continuous;
  makeMutationData(x2, y2Continuous, n2);
  std::vector<double> y2(n2);
  for (size_t i = 0; i < n2; ++i) y2[i] = y2Continuous[i] > 0.0 ? 1.0 : 0.0;

  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);
  check(sampler.numObservations() == n2, "probit setData resizes");

  // latents cold-initialize to 2 y - 1, as the reference engine does
  bool latentsMatch = true;
  for (size_t i = 0; i < n2; ++i)
    latentsMatch &= sampler.latents(0)[i] == 2.0 * y2[i] - 1.0;
  check(latentsMatch, "probit setData cold-initializes latents");

  Results results;
  sampler.run(0, 5, results);
  bool occupied = true;
  for (size_t t = 0; t < 25; ++t)
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(occupied, "probit sampler runs after setData");

  printf("ok: setData probit\n");
}

static void testMultiChainSetData() {
  const size_t n = 200, n2 = 260, numChains = 2;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  std::vector<ext_rng*> rngs(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 3000 + static_cast<uint_least32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = numChains;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
                         1.0, 3.0, 0.37804942330213542, options, rngs.data());
  Results empty;
  sampler.run(100, 0, empty);

  std::vector<double> x2, y2;
  makeMutationData(x2, y2, n2);
  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);

  bool occupied = true;
  for (size_t c = 0; c < numChains; ++c)
    for (size_t t = 0; t < 25; ++t)
      occupied &= sampler.chain(c).tree(t).bottomNodesAreOccupied();
  check(occupied, "multi-chain setData leaves every chain valid");

  std::vector<double> sigmaDraws(5 * numChains);
  Results results;
  results.sigma = sigmaDraws.data();
  sampler.run(0, 5, results);
  bool sigmaFinite = true;
  for (double s : sigmaDraws) sigmaFinite &= std::isfinite(s) && s > 0.0;
  check(sigmaFinite, "multi-chain sampler runs after setData");

  for (size_t c = numChains; c > 0; --c) ext_rng_destroy(rngs[c - 1]);

  printf("ok: multi-chain setData\n");
}

static void testCategoricalMechanics() {
  const size_t n = 120;
  // column 0: categorical with 4 levels; column 1: continuous
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    x[i + n] = runif01();
    y[i] = static_cast<double>(i % 4);
  }

  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  check(store.numCuts[0] == 4, "categorical column counts its categories");
  check(store.cutPoints[0].empty(), "categorical column keeps no cut points");
  bool codesMatch = true;
  for (size_t i = 0; i < n; ++i)
    codesMatch &= store.codes[i] == static_cast<xint_t>(i % 4);
  check(codesMatch, "categorical codes are the values");
  check(store.categoricalValueIsValid(0, 3.0) &&
          !store.categoricalValueIsValid(0, 4.0) &&
          !store.categoricalValueIsValid(0, 1.5),
        "categorical value validity");

  // rule sending {1, 3} right partitions by bit test
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;
  rule.variableIndex = 0;
  rule.categoryDirections = (1u << 1) | (1u << 3);
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  check(left.numObservations() == n / 2 && right.numObservations() == n / 2,
        "categorical partition splits by mask");
  bool sidesMatch = true;
  for (size_t i = left.begin; i < left.end; ++i) {
    size_t category = tree.indices[i] % 4;
    sidesMatch &= category == 0 || category == 2;
  }
  check(sidesMatch, "left side holds only left-bound categories");
  checkNear(left.average, (0.0 + 2.0) / 2.0, 1e-12, "left average by category");

  // reachability filters through nested rules
  check(tree.reachableCategories(store, 0, 0) == 0xfu,
        "root reaches all categories");
  check(tree.reachableCategories(store, tree.at(0).leftChild, 0) ==
          ((1u << 0) | (1u << 2)),
        "left child reaches left-bound categories");
  check(tree.reachableCategories(store, tree.at(0).leftChild + 1, 0) ==
          ((1u << 1) | (1u << 3)),
        "right child reaches right-bound categories");
  check(tree.variableAvailable(store, tree.at(0).leftChild, 0),
        "two reachable categories keep the variable available");

  // a second split below exhausts one side
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.categoryDirections = 1u << 2;
  tree.birth(store, tree.at(0).leftChild, childRule, y.data(), nullptr);
  int32_t grandLeft = tree.at(tree.at(0).leftChild).leftChild;
  check(tree.reachableCategories(store, grandLeft, 0) == (1u << 0),
        "grandchild reaches a single category");
  check(!tree.variableAvailable(store, grandLeft, 0),
        "one reachable category exhausts the variable");

  // routing agrees with the partitions, including a code override
  bool routesMatch = true;
  for (size_t i = 0; i < n; i += 7) {
    int32_t leaf = tree.findBottomNodeForObservation(store, i);
    size_t inLeaf = 0;
    for (size_t k = tree.at(leaf).begin; k < tree.at(leaf).end; ++k)
      inLeaf += tree.indices[k] == i ? 1 : 0;
    routesMatch &= inLeaf == 1;
  }
  check(routesMatch, "categorical routing matches partitions");
  int32_t overridden = tree.findBottomNodeForObservation(
    store, 0, 0, static_cast<xint_t>(3));
  check(overridden == tree.at(0).leftChild + 1,
        "code override routes to the right side");

  printf("ok: categorical mechanics\n");
}

static void testCategoricalPriorMath(ext_rng* rng) {
  const size_t n = 120;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 5);

  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);

  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);

  CGMTreePrior prior;

  // R = 5 reachable: 2^5 - 2 = 30 valid assignments, drawn uniformly with
  // unreachable bits (none here) zero and neither side empty
  const int numDraws = 60000;
  std::vector<int> patternCounts(1u << 5, 0);
  for (int i = 0; i < numDraws; ++i) {
    Rule rule = prior.drawRuleForVariable(tree, store, rng, 0, 0);
    check(rule.categoryDirections > 0 && rule.categoryDirections < 31u + 1u &&
            rule.categoryDirections != 31u,
          "categorical draw leaves neither side empty");
    ++patternCounts[rule.categoryDirections];
  }
  bool uniform = patternCounts[0] == 0 && patternCounts[31] == 0;
  double expected = static_cast<double>(numDraws) / 30.0;
  for (std::uint32_t pattern = 1; pattern < 31; ++pattern)
    uniform = uniform &&
      std::fabs(patternCounts[pattern] - expected) < 5.0 * std::sqrt(expected);
  check(uniform, "categorical rule draw is uniform over valid assignments");

  Rule rootRule;
  rootRule.variableIndex = 0;
  rootRule.categoryDirections = (1u << 3) | (1u << 4);
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -std::log(30.0), 1e-13, "root categorical rule probability");

  // left child reaches {0, 1, 2}: R = 3 gives 2^3 - 2 = 6
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.categoryDirections = 1u << 1;
  tree.birth(store, tree.at(0).leftChild, childRule, y.data(), nullptr);
  checkNear(prior.ruleForVariableLogProbability(tree, store, tree.at(0).leftChild),
            -std::log(6.0), 1e-13, "child categorical rule probability");

  // full tree log probability: root splits (g0), children of root: left
  // splits (g1), right leaf (1 - g1'); right of root reaches {3,4} so g1'
  // uses availability; grandchildren reach single categories, forced leaves
  double g0 = prior.growthProbability(tree, store, 0);
  double g1 = prior.growthProbability(tree, store, tree.at(0).leftChild);
  double gRight = prior.growthProbability(tree, store, tree.at(0).leftChild + 1);
  int32_t leftChild = tree.at(0).leftChild;
  double gGrandLeft = prior.growthProbability(tree, store, tree.at(leftChild).leftChild);
  double gGrandRight = prior.growthProbability(tree, store, tree.at(leftChild).leftChild + 1);
  // the left child's rule sends category 1 right: its left grandchild
  // reaches {0, 2} and can grow, the right reaches only {1} and cannot
  check(gGrandLeft != 0.0, "two-category node can grow");
  check(gGrandRight == 0.0, "single-category node cannot grow");
  double expectedLogProbability =
    std::log(g0) - std::log(30.0) +
    std::log(g1) - std::log(6.0) +
    std::log(1.0 - gRight) +
    std::log(1.0 - gGrandLeft) + std::log(1.0 - gGrandRight);
  checkNear(prior.treeLogProbability(tree, store), expectedLogProbability,
            1e-12, "categorical tree log probability");

  printf("ok: categorical prior math\n");
}

static void testEndToEndCategorical(ext_rng* rng) {
  const size_t n = 400;
  // column 0: 4 categories with distinct means that are NOT ordered in
  // category-code order, so ordinal splits could not represent them
  // cheaply; column 1: continuous noise dimension
  std::vector<double> x(n * 2), y(n);
  const double categoryMeans[] = {2.0, -1.0, 3.0, 0.0};
  for (size_t i = 0; i < n; ++i) {
    size_t category = i % 4;
    x[i] = static_cast<double>(category);
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = categoryMeans[category] + 0.3 * normal;
  }

  SamplerOptions options;
  options.numTrees = 50;
  options.columnTypes = nullptr;  // set via ctor path below
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  options.columnTypes = types;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);

  const size_t numBurnIn = 150, numSamples = 200;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(numBurnIn, numSamples, results);

  bool meansRecovered = true;
  for (size_t category = 0; category < 4; ++category) {
    double sum = 0.0;
    size_t count = 0;
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t i = category; i < n; i += 4) {
        sum += trainingFits[i + s * n];
        ++count;
      }
    meansRecovered &=
      std::fabs(sum / static_cast<double>(count) - categoryMeans[category]) < 0.2;
  }
  check(meansRecovered, "categorical splits recover category means");

  printf("ok: end-to-end categorical\n");
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
  options.columnTypes = types;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  std::vector<xint_t> codesBefore(sampler.data().codes);
  check(sampler.setPredictor(xBad.data(), false, false) ==
          PredictorUpdateResult::invalidCutPoints,
        "out-of-range category code refused");
  check(sampler.data().codes == codesBefore, "refusal mutates nothing");

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

static void testWideCategorical(ext_rng* rng) {
  // masks are 64 bits wide, capped at 53 categories by the flattened
  // format's double encoding; exercise bit positions past 31 throughout
  const size_t K = 53;
  const size_t n = 8 * K;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    size_t category = i % K;
    x[i] = static_cast<double>(category);
    y[i] = (category >= 27 ? 2.0 : 0.0) + 0.3 * (runif01() - 0.5);
  }

  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);
  check(store.numCuts[0] == K, "wide categorical column counts its categories");
  check(store.categoricalValueIsValid(0, 52.0) &&
          !store.categoricalValueIsValid(0, 53.0),
        "wide categorical value validity");

  // partition and reachability with direction bits past 31
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;
  rule.variableIndex = 0;
  rule.categoryDirections = (1ull << 1) | (1ull << 40) | (1ull << 52);
  tree.birth(store, 0, rule, y.data(), nullptr);
  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  check(right.numObservations() == 8 * 3 &&
          left.numObservations() == n - 8 * 3,
        "wide mask partitions by high bits");
  bool sidesMatch = true;
  for (size_t i = right.begin; i < right.end; ++i) {
    size_t category = tree.indices[i] % K;
    sidesMatch &= category == 1 || category == 40 || category == 52;
  }
  check(sidesMatch, "right side holds only right-bound categories");
  check(tree.reachableCategories(store, 0, 0) == (1ull << K) - 1ull,
        "root reaches all 53 categories");
  check(tree.reachableCategories(store, tree.at(0).leftChild + 1, 0) ==
          rule.categoryDirections,
        "right child reaches the mask's categories");

  // equals must compare the full mask width (the swap move's same-rule
  // test); ordinal rules zero the high word at construction
  Rule wideA, wideB;
  wideA.variableIndex = wideB.variableIndex = 0;
  wideA.categoryDirections = 1ull | (1ull << 40);
  wideB.categoryDirections = 1ull;
  check(!wideA.equals(wideB), "equals sees high mask bits");
  Rule ordinalA, ordinalB;
  ordinalA.variableIndex = ordinalB.variableIndex = 1;
  ordinalA.splitIndex = 5;
  ordinalB.splitIndex = 5;
  check(ordinalA.equals(ordinalB), "fresh ordinal rules compare equal");

  // each pattern bit is marginally exactly 1/2 under the uniform prior on
  // nonempty assignments; the single range draw this replaced pinned low
  // pattern bits for wide masks
  Tree drawTree;
  drawTree.initialize(indices.data(), n);
  CGMTreePrior prior;
  const int numDraws = 20000;
  std::vector<int> bitCounts(K, 0);
  for (int d = 0; d < numDraws; ++d) {
    Rule drawn = prior.drawRuleForVariable(drawTree, store, rng, 0, 0);
    check(drawn.categoryDirections > 0 &&
            drawn.categoryDirections < (1ull << K) - 1ull,
          "wide draw leaves neither side empty");
    for (size_t bit = 0; bit < K; ++bit)
      bitCounts[bit] +=
        static_cast<int>((drawn.categoryDirections >> bit) & 1);
  }
  bool marginsMatch = true;
  double tolerance = 5.0 * std::sqrt(0.25 * numDraws);
  for (size_t bit = 0; bit < K; ++bit)
    marginsMatch &=
      std::fabs(bitCounts[bit] - 0.5 * numDraws) < tolerance;
  check(marginsMatch, "wide draw direction bits are marginally uniform");
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -std::log(std::pow(2.0, 53.0) - 2.0), 1e-13,
            "wide rule probability");

  // the widest valid mask round-trips through the double-valued flat format
  std::vector<FlatNode> flat;
  std::vector<double> paramByNode(tree.nodes.size(), 0.0);
  tree.flatten(store, paramByNode.data(), flat);
  check(flat[0].value == static_cast<double>(rule.categoryDirections),
        "wide mask flattens exactly");
  flat[0].value = static_cast<double>(((1ull << K) - 1ull) & ~1ull);
  check(flatTreeIsWellFormed(store, flat.data(), flat.size()),
        "all-but-one mask is well formed");
  Tree rebuilt;
  std::vector<size_t> rebuiltIndices(n);
  rebuilt.initialize(rebuiltIndices.data(), n);
  std::vector<double> rebuiltParams;
  check(rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams),
        "wide mask rebuilds from flat");
  check(rebuilt.at(0).rule.categoryDirections == (((1ull << K) - 1ull) & ~1ull),
        "wide mask round-trips exactly");
  flat[0].value = 9007199254740992.0;  // 2^53: past the exact-integer bound
  check(!flatTreeIsWellFormed(store, flat.data(), flat.size()),
        "mask past the double-exact bound is rejected");

  // end to end: category codes carry no order, so recovering the two
  // groups requires subset splits over high codes
  SamplerOptions options;
  options.numTrees = 50;
  options.columnTypes = types;
  ClassicSampler sampler(x.data(), y.data(), n, 1, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  const size_t numSamples = 100;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(100, numSamples, results);
  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i) {
      if (i % K >= 27) { highSum += trainingFits[i + s * n]; ++highCount; }
      else             { lowSum += trainingFits[i + s * n]; ++lowCount; }
    }
  check(std::fabs(highSum / static_cast<double>(highCount) - 2.0) < 0.2 &&
          std::fabs(lowSum / static_cast<double>(lowCount)) < 0.2,
        "splits over wide categories recover group means");

  printf("ok: wide categorical\n");
}

static void testFlattenRoundTrip() {
  // one categorical column (codes 0..3) and one ordinal, a hand-built tree
  // with one rule of each kind
  const size_t n = 12;
  std::vector<double> x(n * 2), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  for (size_t i = 0; i < n; ++i) x[i + n] = static_cast<double>(i) / (n - 1.0);
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};

  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rootRule;
  rootRule.variableIndex = 0;
  rootRule.categoryDirections = 0x6;  // categories 1, 2 right
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;
  leftRule.variableIndex = 1;
  leftRule.splitIndex = 4;
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);

  std::vector<double> params(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  for (size_t k = 0; k < bottoms.size(); ++k)
    params[static_cast<size_t>(bottoms[k])] = static_cast<double>(k + 1);

  std::vector<FlatNode> flat;
  std::vector<std::uint32_t> counts;
  tree.flatten(store, params.data(), flat, &counts);

  check(flat.size() == 5 && counts.size() == 5, "flatten emits pre-order");
  check(flat[0].variable == 0 && flat[0].value == 6.0,
        "flatten value-encodes a categorical mask");
  check(flat[1].variable == 1 && flat[1].value == store.cutPoints[1][4],
        "flatten value-encodes an ordinal cut");
  check(flat[2].variable == invalidVariable && flat[2].value == 1.0 &&
        flat[3].value == 2.0 && flat[4].value == 3.0,
        "flatten stores leaf parameters left-first");
  check(counts[0] == n &&
        counts[1] == static_cast<std::uint32_t>(
                       tree.at(tree.at(0).leftChild).numObservations()) &&
        counts[2] + counts[3] == counts[1],
        "flatten counts mirror the partitions");

  // replay against the raw predictors reproduces the live counts
  std::vector<std::uint32_t> replayed(flat.size());
  std::vector<size_t> replayIndices(n);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  countFlatObservationsBelow(flat.data(), store.types.data(), x.data(), n,
                             replayIndices.data(), 0, n, replayed.data());
  check(replayed == counts, "flat replay reproduces the live counts");

  // per-row prediction accumulates each row's leaf parameter
  std::vector<double> fits(n, 0.0);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  addFlatPredictionsBelow(flat.data(), store.types.data(), x.data(), n,
                          replayIndices.data(), 0, n, fits.data());
  bool fitsMatch = true;
  for (size_t i = 0; i < n; ++i) {
    int32_t leaf = tree.findBottomNodeForObservation(store, i);
    fitsMatch &= fits[i] == params[static_cast<size_t>(leaf)];
  }
  check(fitsMatch, "flat prediction routes rows to their leaves");

  // rebuild recovers the rules exactly
  std::vector<size_t> indices2(n);
  Tree tree2;
  tree2.initialize(indices2.data(), n);
  std::vector<double> params2;
  check(tree2.buildFromFlat(store, flat.data(), flat.size(), params2),
        "buildFromFlat accepts its own flatten");
  check(tree2.at(0).rule.variableIndex == 0 &&
        tree2.at(0).rule.categoryDirections == 0x6,
        "buildFromFlat recovers the categorical mask");
  check(tree2.at(tree2.at(0).leftChild).rule.splitIndex == 4,
        "buildFromFlat recovers the ordinal split index exactly");
  tree2.repartitionSubtree(store, 0);
  bool partitionsMatch = true;
  for (size_t i = 0; i < tree.nodes.size(); ++i)
    partitionsMatch &= tree2.at(static_cast<int32_t>(i)).numObservations() ==
                       tree.at(static_cast<int32_t>(i)).numObservations();
  check(partitionsMatch, "rebuilt tree partitions identically");
  bool paramsMatch = true;
  for (int32_t i : bottoms)
    paramsMatch &= params2[static_cast<size_t>(i)] ==
                   params[static_cast<size_t>(i)];
  check(paramsMatch, "buildFromFlat recovers leaf parameters");

  // malformed inputs: an ordinal value off the cut grid, a mask outside the
  // canonical gauge, and a truncated pre-order all refuse
  std::vector<FlatNode> bad(flat);
  bad[1].value += 1.0e-3;
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, bad.data(), bad.size(), params2),
        "buildFromFlat rejects a value off the cut grid");
  bad = flat;
  bad[1].variable = 0;
  bad[1].value = 2.0;  // category 1 goes right, but 1 is unreachable here
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, bad.data(), bad.size(), params2),
        "buildFromFlat rejects an out-of-gauge mask");
  tree2.initialize(indices2.data(), n);
  check(!tree2.buildFromFlat(store, flat.data(), flat.size() - 1, params2),
        "buildFromFlat rejects a truncated pre-order");
  check(flatTreeIsWellFormed(store, flat.data(), flat.size()) &&
        !flatTreeIsWellFormed(store, flat.data(), flat.size() - 1),
        "well-formedness check matches");

  printf("ok: flatten round trip\n");
}

static void testKeepTrees(ext_rng* rng) {
  const size_t n = 200, nTest = 20, numSamples = 4;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();

  SamplerOptions options;
  options.numTrees = 25;
  options.keepTrees = true;
  options.numSamplesToStore = numSamples;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  sampler.setTestPredictors(xTest.data(), nTest);

  Results empty;
  sampler.run(100, 0, empty);
  check(sampler.currentSampleNum() == 0, "burn-in does not advance the slot");

  std::vector<double> sigma(numSamples), testFits(nTest * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.testFits = testFits.data();
  sampler.run(0, numSamples, results);
  check(sampler.currentSampleNum() == 0, "a full run wraps the slot");

  // replaying the saved trees against the raw test rows must reproduce the
  // run's recorded test fits exactly: same parameters, same addition order
  std::vector<double> predicted(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, predicted.data());
  check(predicted == testFits, "saved-tree predictions equal the run's test fits");

  // a partial run overwrites the oldest slots and leaves the rest
  std::vector<double> sigma2(2), testFits2(nTest * 2);
  Results results2;
  results2.sigma = sigma2.data();
  results2.testFits = testFits2.data();
  sampler.run(0, 2, results2);
  check(sampler.currentSampleNum() == 2, "a partial run advances the slot");

  std::vector<double> predicted2(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, predicted2.data());
  bool overwritten = std::equal(testFits2.begin(), testFits2.end(),
                                predicted2.begin());
  bool preserved = std::equal(predicted.begin() + nTest * 2, predicted.end(),
                              predicted2.begin() + nTest * 2);
  check(overwritten, "new samples overwrite the oldest slots");
  check(preserved, "later slots survive a partial run");
  check(sampler.savedTreeCapacity() == numSamples,
        "keepTrees capacity comes from the options");

  printf("ok: keepTrees\n");
}

static void testPredictCurrentTrees(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  check(sampler.savedTreeCapacity() == 0, "keepTrees defaults off");

  // routing the training rows through the live trees agrees with the run's
  // recorded training fits (both are scale * sum-of-tree-fits + shift; only
  // the accumulation order differs)
  std::vector<double> trainingFits(n);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(0, 1, results);

  std::vector<double> predicted(n);
  sampler.predict(x.data(), n, predicted.data());
  bool match = true;
  for (size_t i = 0; i < n; ++i)
    match &= std::fabs(predicted[i] - trainingFits[i]) <= 1.0e-8;
  check(match, "live-tree predictions equal the recorded training fits");

  printf("ok: predict from current trees\n");
}

static void testStateRoundTripScaledOffset() {
  // setOffset(updateScale) moves the gaussian response transform after
  // creation; the state must carry it or a restored sampler mis-scales
  // every internal quantity (the classic engine forced hosts to export the
  // scale themselves)
  const size_t n = 200, numSamples = 5;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> offset(n);
  for (size_t i = 0; i < n; ++i)
    offset[i] = 2.0 * std::sin(0.1 * (double) i);  // widens y - offset

  SamplerOptions options;
  options.numTrees = 25;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 77) != 0 ||
      ext_rng_setSeed(rngB, 78) != 0) {
    check(false, "scaled state: rng creation");
    return;
  }

  ClassicSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngA);
  Results empty;
  original.run(20, 0, empty);
  original.setOffset(offset.data(), true);
  original.run(20, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(state.chains[0].fitMax > state.chains[0].fitMin,
        "gaussian state captures the response transform");

  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngB);
  // a host reinstalls the current offset but cannot reproduce the scale
  // trajectory that produced the state; the stored transform must win
  restored.setOffset(offset.data(), false);
  check(restored.setState(state), "scaled state restores");

  std::vector<double> sigmaA(numSamples), trainA(n * numSamples);
  Results resultsA;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  original.run(0, numSamples, resultsA);

  std::vector<double> sigmaB(numSamples), trainB(n * numSamples);
  Results resultsB;
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  restored.run(0, numSamples, resultsB);

  check(sigmaA == sigmaB, "scaled restore continues sigma bitwise");
  check(trainA == trainB, "scaled restore continues fits bitwise");

  std::vector<double> xTest(20 * 2);
  for (double& v : xTest) v = runif01();
  std::vector<double> predictionsA(20), predictionsB(20);
  original.predict(xTest.data(), 20, predictionsA.data());
  restored.predict(xTest.data(), 20, predictionsB.data());
  check(predictionsA == predictionsB,
        "scaled restore predicts on the original scale");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: state round-trip with a moved scale\n");
}

static void testStateRoundTrip() {
  // the strong gate: store the state, continue the original, and continue a
  // fresh sampler restored from the state; the draws must agree bitwise
  const size_t n = 200, numChains = 2, numSamples = 5;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  };

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = numChains;
  options.keepTrees = true;
  options.numSamplesToStore = 3;

  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 1001);
  ClassicSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(60, 2, empty);  // leave a nonzero slot position behind

  SamplerStateData state;
  original.getState(state);
  check(state.chains.size() == numChains && state.currentSampleNum == 2,
        "state captures every chain and the slot position");

  // different seeds on purpose: the serialized rng state must win
  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 9999);
  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state), "a stored state restores");

  std::vector<double> sigmaA(numSamples * numChains);
  std::vector<double> trainA(n * numSamples * numChains);
  Results resultsA;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  original.run(0, numSamples, resultsA);

  std::vector<double> sigmaB(numSamples * numChains);
  std::vector<double> trainB(n * numSamples * numChains);
  Results resultsB;
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  restored.run(0, numSamples, resultsB);

  check(sigmaA == sigmaB, "restored chains draw identical sigmas");
  check(trainA == trainB, "restored chains draw identical fits");

  // saved trees also round-tripped: predictions agree
  std::vector<double> xTest(20 * 2);
  for (double& v : xTest) v = runif01();
  size_t capacity = original.savedTreeCapacity();
  std::vector<double> predictA(20 * capacity * numChains);
  std::vector<double> predictB(20 * capacity * numChains);
  original.predict(xTest.data(), 20, predictA.data());
  restored.predict(xTest.data(), 20, predictB.data());
  check(predictA == predictB, "saved trees ride along with the state");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: state round trip\n");
}

static void testStateRoundTripLatents(ext_rng* rng) {
  // logistic Polya-Gamma latents and dart state must ride along too
  const size_t n = 150, numSamples = 4;
  std::vector<double> x(n * 2), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = (x[i] + 0.25 * (runif01() - 0.5) > 0.5) ? 1.0 : 0.0;

  SamplerOptions options;
  options.numTrees = 10;
  options.nodeScale = 3.0;
  options.useDart = true;
  options.dart.updateDelay = 10;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngA, 555);
  ClassicSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                          &rngA);
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(!state.chains[0].latents.empty(), "logistic state carries omega");
  check(!state.chains[0].dartProbabilities.empty(), "dart state is captured");

  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngB, 777);
  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options,
                          &rngB);
  check(restored.setState(state), "a binary+dart state restores");

  std::vector<double> trainA(n * numSamples), trainB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.trainingFits = trainA.data();
  resultsB.trainingFits = trainB.data();
  original.run(0, numSamples, resultsA);
  restored.run(0, numSamples, resultsB);
  check(trainA == trainB, "restored logistic + dart chains draw identically");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  (void) rng;
  printf("ok: state round trip with latents\n");
}

static void testStateValidation(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ClassicSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ClassicSampler& sampler(*samplerPtr);

  SamplerStateData state;
  sampler.getState(state);
  std::vector<double> treeFitsBefore(sampler.chain(0).treeFits());
  std::vector<std::vector<double>> cutsBefore(sampler.data().cutPoints);

  SamplerStateData bad(state);
  bad.chains.pop_back();
  check(!sampler.setState(bad), "setState rejects a chain-count mismatch");

  bad = state;
  bool corrupted = false;
  for (std::vector<FlatNode>& tree : bad.chains[0].trees) {
    if (tree.size() > 1) {
      tree[0].value += 1.0e-3;  // off the cut grid
      corrupted = true;
      break;
    }
  }
  check(corrupted && !sampler.setState(bad),
        "setState rejects a split value off the cut grid");
  check(sampler.data().cutPoints == cutsBefore,
        "a rejected state leaves the cuts untouched");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "a rejected state leaves the fits untouched");

  check(sampler.setState(state), "the original state still restores");
  check(sampler.chain(0).treeFits() == treeFitsBefore,
        "restoring the current state is an identity");

  printf("ok: state validation\n");
}

static void testSetWeightsAndTestOffset() {
  const size_t n = 200, nTest = 20;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();
  std::vector<double> w1(n), w2(n);
  for (double& v : w1) v = 0.5 + runif01();
  for (double& v : w2) v = 0.5 + runif01();
  std::vector<double> testOffset(nTest);
  for (size_t i = 0; i < nTest; ++i)
    testOffset[i] = static_cast<double>(i) - 10.0;

  SamplerOptions options;
  options.numTrees = 25;

  auto makeSeededRng = []() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 314);
    return rng;
  };

  // setWeights before any run is indistinguishable from creating with the
  // new weights: a bare pointer swap, nothing rescales
  ext_rng* rngA = makeSeededRng();
  ClassicSampler samplerA(x.data(), y.data(), n, 2, w2.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngA);
  ext_rng* rngB = makeSeededRng();
  ClassicSampler samplerB(x.data(), y.data(), n, 2, w1.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngB);
  samplerB.setWeights(w2.data());

  std::vector<double> sigmaA(5), sigmaB(5), trainA(n * 5), trainB(n * 5);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  samplerA.run(30, 5, resultsA);
  samplerB.run(30, 5, resultsB);
  check(sigmaA == sigmaB && trainA == trainB,
        "setWeights equals creation with the new weights");

  // a test offset shifts recorded test fits by exactly itself and leaves
  // everything else alone
  ext_rng* rngC = makeSeededRng();
  ClassicSampler samplerC(x.data(), y.data(), n, 2, w2.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngC);
  samplerC.setTestPredictors(xTest.data(), nTest);
  samplerC.setTestOffset(testOffset.data());

  ext_rng* rngD = makeSeededRng();
  ClassicSampler samplerD(x.data(), y.data(), n, 2, w2.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngD);
  samplerD.setTestPredictors(xTest.data(), nTest);

  std::vector<double> trainC(n * 5), trainD(n * 5);
  std::vector<double> testC(nTest * 5), testD(nTest * 5);
  Results resultsC, resultsD;
  resultsC.trainingFits = trainC.data();
  resultsC.testFits = testC.data();
  resultsD.trainingFits = trainD.data();
  resultsD.testFits = testD.data();
  samplerC.run(30, 5, resultsC);
  samplerD.run(30, 5, resultsD);

  check(trainC == trainD, "a test offset leaves training fits alone");
  bool shifted = true;
  for (size_t s = 0; s < 5; ++s)
    for (size_t i = 0; i < nTest; ++i)
      shifted &= testC[s * nTest + i] == testD[s * nTest + i] + testOffset[i];
  check(shifted, "a test offset shifts recorded test fits by itself");

  // clearing restores the unshifted fits, and the offset survives a
  // same-length predictor replacement
  samplerC.setTestOffset(nullptr);
  check(samplerC.data().testOffset == nullptr, "a null test offset clears");
  samplerC.setTestOffset(testOffset.data());
  samplerC.setTestPredictors(xTest.data(), nTest);
  check(samplerC.data().testOffset == testOffset.data(),
        "the test offset survives a predictor replacement");

  ext_rng_destroy(rngD);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: setWeights and test offsets\n");
}

static void testSetControlAndModel() {
  const size_t n = 200, nTest = 20;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = runif01();

  auto makeSeededRng = []() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 271828);
    return rng;
  };

  // setModel before any run is indistinguishable from creating with the new
  // model: fixed-k full swap
  double splitProbabilities[2] = {0.7, 0.3};
  SamplerOptions optionsA;
  optionsA.numTrees = 25;
  optionsA.k = 3.0;
  optionsA.nodeScale = 0.7;
  optionsA.base = 0.8;
  optionsA.power = 1.5;
  optionsA.birthOrDeathProbability = 0.6;
  optionsA.swapProbability = 0.1;
  optionsA.changeProbability = 0.3;
  optionsA.birthProbability = 0.4;
  optionsA.splitProbabilities = splitProbabilities;
  ext_rng* rngA = makeSeededRng();
  ClassicSampler samplerA(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 5.0, 0.9, optionsA,
                          &rngA);

  SamplerOptions optionsB;
  optionsB.numTrees = 25;
  ext_rng* rngB = makeSeededRng();
  ClassicSampler samplerB(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsB, &rngB);
  ModelParameters model;
  model.base = 0.8;
  model.power = 1.5;
  model.birthOrDeathProbability = 0.6;
  model.swapProbability = 0.1;
  model.changeProbability = 0.3;
  model.birthProbability = 0.4;
  model.nodeScale = 0.7;
  model.k = 3.0;
  model.sigmaEstimate = 1.0;
  model.sigmaDf = 5.0;
  model.sigmaRawScale = 0.9;
  model.splitProbabilities = splitProbabilities;
  samplerB.setModel(model);

  std::vector<double> sigmaA(5), sigmaB(5), trainA(n * 5), trainB(n * 5);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  samplerA.run(30, 5, resultsA);
  samplerB.run(30, 5, resultsB);
  check(sigmaA == sigmaB && trainA == trainB,
        "setModel equals creation with the new model");

  // the same when swapping in the chi hyperprior on k
  SamplerOptions optionsC;
  optionsC.numTrees = 25;
  optionsC.updateK = true;
  optionsC.kHyperprior.degreesOfFreedom = 2.0;
  optionsC.kHyperprior.scale = 5.0;
  ext_rng* rngC = makeSeededRng();
  ClassicSampler samplerC(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsC, &rngC);

  ext_rng* rngD = makeSeededRng();
  ClassicSampler samplerD(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsB, &rngD);
  ModelParameters modelK;
  modelK.sigmaEstimate = 1.0;
  modelK.sigmaDf = 3.0;
  modelK.sigmaRawScale = 0.37804942330213542;
  modelK.updateK = true;
  modelK.kHyperprior.degreesOfFreedom = 2.0;
  modelK.kHyperprior.scale = 5.0;
  samplerD.setModel(modelK);

  std::vector<double> kC(5), kD(5);
  Results resultsC, resultsD;
  resultsC.sigma = sigmaA.data();
  resultsC.k = kC.data();
  resultsD.sigma = sigmaB.data();
  resultsD.k = kD.data();
  samplerC.run(30, 5, resultsC);
  samplerD.run(30, 5, resultsD);
  check(sigmaA == sigmaB && kC == kD,
        "setModel swaps in the k hyperprior");

  // tree storage toggles on and off between runs, the classic setControl's
  // keepTrees flip: saved predictions reproduce the recorded test fits
  samplerB.setTestPredictors(xTest.data(), nTest);
  check(samplerB.savedTreeCapacity() == 0, "storage starts off");
  samplerB.setTreeStorage(true, 4);
  check(samplerB.savedTreeCapacity() == 4 && samplerB.currentSampleNum() == 0,
        "enabling storage allocates the slots");

  std::vector<double> testFits(nTest * 4);
  Results keepResults;
  keepResults.testFits = testFits.data();
  samplerB.run(0, 4, keepResults);
  std::vector<double> predicted(nTest * 4);
  samplerB.predict(xTest.data(), nTest, predicted.data());
  check(predicted == testFits,
        "post-toggle saved predictions equal the recorded test fits");

  samplerB.setTreeStorage(true, 4);
  check(samplerB.savedTreeCapacity() == 4,
        "an unchanged toggle preserves the storage");
  samplerB.setTreeStorage(false, 0);
  check(samplerB.savedTreeCapacity() == 0, "disabling storage frees it");
  std::vector<double> livePredictions(nTest);
  samplerB.predict(xTest.data(), nTest, livePredictions.data());
  bool liveFinite = true;
  for (double v : livePredictions) liveFinite &= std::isfinite(v);
  check(liveFinite, "live predictions work after disabling storage");

  // setNumThin matches creating with the thinning rate
  SamplerOptions optionsE;
  optionsE.numTrees = 25;
  optionsE.numThin = 3;
  ext_rng* rngE = makeSeededRng();
  ClassicSampler samplerE(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsE, &rngE);
  ext_rng* rngF = makeSeededRng();
  ClassicSampler samplerF(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsB, &rngF);
  samplerF.setNumThin(3);
  samplerF.setNumThreads(4);  // min(numThreads, numChains) still runs inline

  std::vector<double> sigmaE(3), sigmaF(3);
  Results resultsE, resultsF;
  resultsE.sigma = sigmaE.data();
  resultsF.sigma = sigmaF.data();
  samplerE.run(2, 3, resultsE);
  samplerF.run(2, 3, resultsF);
  check(sigmaE == sigmaF, "setNumThin equals creation with the rate");

  // sums of squared residuals descale exactly: the last recorded sample's
  // fits are the current state
  std::vector<double> offset(n);
  for (size_t i = 0; i < n; ++i) offset[i] = 0.5 + 0.01 * static_cast<double>(i);
  SamplerOptions optionsG;
  optionsG.numTrees = 25;
  ext_rng* rngG = makeSeededRng();
  ClassicSampler samplerG(x.data(), y.data(), n, 2, nullptr, offset.data(),
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsG, &rngG);
  std::vector<double> train(n);
  Results resultsG;
  resultsG.trainingFits = train.data();
  samplerG.run(10, 1, resultsG);

  double manual = 0.0;
  for (size_t i = 0; i < n; ++i) {
    // recorded training fits include the offset (original-scale convention)
    double residual = y[i] - train[i];
    manual += residual * residual;
  }
  double reported = samplerG.sumOfSquaredResiduals(0);
  check(std::fabs(reported - manual) <= 1.0e-8 * manual,
        "sums of squared residuals descale to the original scale");

  ext_rng_destroy(rngG);
  ext_rng_destroy(rngF);
  ext_rng_destroy(rngE);
  ext_rng_destroy(rngD);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: setControl and setModel\n");
}

static void testSampleFromPrior(ext_rng* rng) {
  const size_t n = 200, numTrees = 50, numReplications = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  SamplerOptions options;
  options.numTrees = numTrees;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  const ColumnStore& store(sampler.data());

  // structure draws: every tree well-formed with occupied leaves, and the
  // growth frequencies match the depth-decayed prior probabilities (root
  // children of a single ordinal split are never empty, so the root split
  // rate is exactly base; deeper nodes lose a little to empty-leaf collapse)
  size_t numRoots = 0, numRootSplits = 0;
  size_t numDepth1Nodes = 0, numDepth1Splits = 0;
  bool allWellFormed = true, allOccupied = true;

  // parameter draws: with fits rebuilt from the drawn parameters, the leaf
  // values recovered by flattening are the draws themselves
  double paramSum = 0.0, paramSumOfSquares = 0.0;
  size_t numParams = 0;

  std::vector<FlatNode> flat;
  std::vector<std::uint32_t> counts;
  for (size_t r = 0; r < numReplications; ++r) {
    sampler.sampleTreesFromPrior();
    sampler.sampleNodeParametersFromPrior();
    for (size_t t = 0; t < numTrees; ++t) {
      sampler.flattenTree(0, t, flat, counts);
      allWellFormed &= flatTreeIsWellFormed(store, flat.data(), flat.size());
      for (size_t l = 0; l < flat.size(); ++l) {
        if (flat[l].variable != invalidVariable) continue;
        allOccupied &= counts[l] > 0;
        paramSum += flat[l].value;
        paramSumOfSquares += flat[l].value * flat[l].value;
        ++numParams;
      }

      ++numRoots;
      if (flat[0].variable == invalidVariable) continue;
      ++numRootSplits;
      size_t numOnLeft = flatSubtreeLength(flat.data() + 1);
      numDepth1Nodes += 2;
      numDepth1Splits += (flat[1].variable != invalidVariable ? 1 : 0) +
                         (flat[1 + numOnLeft].variable != invalidVariable ? 1 : 0);
    }
  }

  check(allWellFormed, "prior-drawn trees are well-formed");
  check(allOccupied, "prior-drawn trees have occupied leaves");
  checkNear(static_cast<double>(numRootSplits) / static_cast<double>(numRoots),
            0.95, 0.015, "root splits at the prior rate");
  checkNear(
    static_cast<double>(numDepth1Splits) / static_cast<double>(numDepth1Nodes),
    0.95 / 4.0, 0.03, "depth-1 splits near the prior rate");

  double paramMean = paramSum / static_cast<double>(numParams);
  double paramSd = std::sqrt(
    paramSumOfSquares / static_cast<double>(numParams) - paramMean * paramMean);
  double priorSd = (0.5 / std::sqrt(static_cast<double>(numTrees))) / 2.0;
  checkNear(paramMean, 0.0, 0.001, "prior leaf parameters center at zero");
  checkNear(paramSd, priorSd, 0.002,
            "prior leaf parameters have the prior spread");

  // the sampler keeps running from a prior-drawn state
  std::vector<double> sigma(2), trainingFits(n * 2);
  Results results;
  results.sigma = sigma.data();
  results.trainingFits = trainingFits.data();
  sampler.run(0, 2, results);
  bool finite = true;
  for (double v : sigma) finite &= std::isfinite(v);
  for (double v : trainingFits) finite &= std::isfinite(v);
  check(finite, "a run continues from a prior-drawn state");

  printf("ok: sample from prior\n");
}

int main() {
  misc_simd_init();

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rng == NULL || ext_rng_setSeed(rng, 42) != 0) {
    printf("FAIL: could not create rng\n");
    return 1;
  }

  testColumnStoreCodes();
  testColumnStoreView();
  testViewSamplerMatchesFull();
  testIntegratedLikelihood();
  testPosteriorDraw(rng);
  testTreeMechanics();
  testTreePriorMath();
  testEndToEndGaussian(rng);
  testEndToEndProbit(rng);
  testDartUpdate(rng);
  testDartSparsityRecovery(rng);
  testChiKHyperprior(rng);
  testColumnStoreMutation();
  testSetPredictorTransaction(rng);
  testSetPredictorForced(rng);
  testUpdatePredictorColumns(rng);
  testPerObservationUpdate(rng);
  testJointPerObservationUpdate();
  testQuantileCutPoints();
  testQuantilePredictorUpdate(rng);
  testSetCutPoints(rng);
  testMultiChain();
  testMultiChainMutation();
  testMapOldCutPointsOntoNew();
  testSetData(rng);
  testSetDataResize(rng);
  testSetDataQuantileShrink(rng);
  testSetDataProbit(rng);
  testMultiChainSetData();
  testPolyaGamma(rng);
  testEndToEndLogistic(rng);
  testLogisticMutation(rng);
  testCategoricalMechanics();
  testCategoricalPriorMath(rng);
  testEndToEndCategorical(rng);
  testCategoricalMutation(rng);
  testWideCategorical(rng);
  testFlattenRoundTrip();
  testKeepTrees(rng);
  testPredictCurrentTrees(rng);
  testStateRoundTrip();
  testStateRoundTripScaledOffset();
  testStateRoundTripLatents(rng);
  testStateValidation(rng);
  testSetWeightsAndTestOffset();
  testSampleFromPrior(rng);
  testSetControlAndModel();

  ext_rng_destroy(rng);

  if (failures > 0) {
    printf("%d failure(s)\n", failures);
    return 1;
  }
  printf("all tests passed\n");
  return 0;
}
