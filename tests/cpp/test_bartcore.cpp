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

  Rule rule{0, 49};
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
  ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, false,
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
  ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, true,
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
    ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, false,
                           1.0, 3.0, 0.37804942330213542, options, &rng);
    const size_t numSamples = 300;
    std::vector<uint32_t> varcount(p * numSamples);
    Results results;
    results.variableCounts = varcount.data();
    sampler.run(200, numSamples, results);

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
    x.data(), y.data(), n, size_t(2), nullptr, nullptr, false, 1.0, 3.0,
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
    xA.data(), yA.data(), n, size_t(2), nullptr, nullptr, false, 1.0, 3.0,
    0.37804942330213542, options, &rngA);
  SamplerFacade<ConstantGaussianLeaf> samplerB(
    xB.data(), yB.data(), n, size_t(2), nullptr, nullptr, false, 1.0, 3.0,
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, false,
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
    ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, false,
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, false,
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, false,
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, true,
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
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, false,
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

int main() {
  misc_simd_init();

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rng == NULL || ext_rng_setSeed(rng, 42) != 0) {
    printf("FAIL: could not create rng\n");
    return 1;
  }

  testColumnStoreCodes();
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

  ext_rng_destroy(rng);

  if (failures > 0) {
    printf("%d failure(s)\n", failures);
    return 1;
  }
  printf("all tests passed\n");
  return 0;
}
