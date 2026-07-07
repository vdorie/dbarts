// Component and smoke tests for bartcore against independently coded
// reference math. Exit code 0 on success.

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <vector>

#include <misc/partition.h>
#include <misc/simd.h>
#include <external/random.h>

#include <bartcore/bartcore.hpp>

using namespace bartcore;

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
  // raw weighted suffstat over three responses (2 * 0.5 + 1 * -0.2 + 3 * 0.1)
  double sumW = 6.0, sumWZ = 1.1, sumWZSq = 0.61;

  // independent transcription of the CGM marginal in crossproduct form
  double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
  double posteriorPrecision = sumW / sigmaSq;
  double mean = sumWZ / sumW;
  double centered = sumWZSq - sumWZ * mean;
  double expected = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision))
    - 0.5 * centered / sigmaSq
    - 0.5 * ((priorPrecision * mean) * (posteriorPrecision * mean)) /
        (priorPrecision + posteriorPrecision);

  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ, sumWZSq),
            expected, 1e-13, "integrated likelihood formula");
  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, 0.0, 0.0, 0.0),
            0.0, 0.0, "empty leaf contributes zero");
  printf("ok: integrated likelihood\n");
}

static void testPosteriorDraw(ext_rng* rng) {
  ConstantGaussianLeaf leaf{0.5 / std::sqrt(50.0)};
  double k = 2.0, sigmaSq = 0.02, sumW = 30.0, sumWZ = 3.6;

  double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
  double posteriorPrecision = sumW / sigmaSq;
  double expectedMean =
    (sumWZ / sigmaSq) / (priorPrecision + posteriorPrecision);
  double expectedSd = 1.0 / std::sqrt(priorPrecision + posteriorPrecision);

  const int numDraws = 200000;
  double sum = 0.0, sumSq = 0.0;
  for (int i = 0; i < numDraws; ++i) {
    double draw = leaf.drawFromPosterior(rng, k, sumW, sumWZ, sigmaSq);
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

// The suffstat marginal must equal the classic (mean, effective count,
// centered variance) form the leaf used to consume, on the same raw data,
// weighted and unweighted, and against the node-context path.
static void testConstantLeafSuffstatEquivalence() {
  ConstantGaussianLeaf leaf{0.5 / std::sqrt(20.0)};
  double k = 1.7, sigmaSq = 0.05;
  const size_t n = 6;
  double z[n] = {0.35, -0.10, 0.22, 0.55, -0.42, 0.16};
  double w[n] = {1.0, 2.0, 0.5, 1.3, 0.8, 1.1};

  for (int weighted = 0; weighted < 2; ++weighted) {
    double sumW = 0.0, sumWZ = 0.0, sumWZSq = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double wi = weighted ? w[i] : 1.0;
      sumW += wi; sumWZ += wi * z[i]; sumWZSq += wi * z[i] * z[i];
    }
    double average = sumWZ / sumW;
    double variance = (sumWZSq - sumW * average * average) / (double) (n - 1);

    // the classic three-argument form, transcribed here
    double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
    double posteriorPrecision = sumW / sigmaSq;
    double classic =
      0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision))
      - 0.5 * (variance / sigmaSq) * (double) (n - 1)
      - 0.5 * ((priorPrecision * average) * (posteriorPrecision * average)) /
          (priorPrecision + posteriorPrecision);

    checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ, sumWZSq),
              classic, 1e-12,
              weighted ? "suffstat marginal equals classic, weighted"
                       : "suffstat marginal equals classic, unweighted");
  }

  // and through the node-context path, over a split's two children
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = (double) i;
  ColumnStore store;
  store.build(x.data(), n, 1, 100);
  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  std::vector<double> zv(z, z + n), wv(w, w + n);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(2);
  tree.birth(store, 0, rule, zv.data(), wv.data());
  int32_t left = tree.at(0).leftChild;

  auto reference = [&](int32_t nodeIndex) {
    const Node& node = tree.at(nodeIndex);
    double sumW = 0.0, sumWZ = 0.0, sumWZSq = 0.0;
    for (size_t m = node.begin; m < node.end; ++m) {
      size_t i = tree.indices[m];
      sumW += wv[i]; sumWZ += wv[i] * zv[i]; sumWZSq += wv[i] * zv[i] * zv[i];
    }
    return leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ, sumWZSq);
  };
  checkNear(leaf.logIntegratedLikelihoodForNode(tree, zv.data(), wv.data(), k,
                                                sigmaSq, left),
            reference(left), 1e-12, "node-context marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(tree, zv.data(), wv.data(), k,
                                                sigmaSq, left + 1),
            reference(left + 1), 1e-12, "node-context marginal, right child");
  printf("ok: constant leaf suffstat equivalence\n");
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
  checkNear(tree.at(0).sumWeightedResponse / tree.at(0).sumWeights, 2.0, 1e-12,
            "root average");

  // split at the median cut
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(49);
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
  checkNear(left.sumWeightedResponse, leftSum, 1e-12, "left child sum wz");
  checkNear(right.sumWeightedResponse, rightSum, 1e-12, "right child sum wz");

  // split interval of the left child is bounded by the parent's rule
  int32_t lo, hi;
  tree.splitInterval(store, tree.at(0).leftChild, 0, &lo, &hi);
  check(lo == 0 && hi == 48, "left child split interval");
  tree.splitInterval(store, tree.at(0).leftChild + 1, 0, &lo, &hi);
  check(lo == 50 && hi == 99, "right child split interval");

  // orphanChildren sums the children's sufficient statistics
  double mergedSumWZ = left.sumWeightedResponse + right.sumWeightedResponse;
  tree.orphanChildren(0);
  checkNear(tree.at(0).sumWeightedResponse, mergedSumWZ, 1e-12,
            "orphan merge sum wz");
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
  rule.setSplitIndex(49);
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

// The cooperative cancellation the R interrupt handler drives: run() polls the
// supplied predicate and returns true when it stops early, on both the
// single-chain (inline poll) and multi-chain (worker cancel flag) paths.
static void testRunCancellation(ext_rng* rng) {
  const size_t n = 300, p = 3;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) y[i] = x[i] + 0.1 * runif01();

  // a poll that never fires runs to completion and returns false, filling
  // the results as usual
  {
    SamplerOptions options;
    options.numTrees = 25;
    ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, options, &rng);
    const size_t numSamples = 20;
    std::vector<double> sigmaDraws(numSamples, 0.0);
    Results results;
    results.sigma = sigmaDraws.data();
    std::function<bool()> never = []() { return false; };
    bool cancelled = sampler.run(10, numSamples, results, never);
    check(!cancelled, "cancellation: non-firing poll completes");
    check(sigmaDraws[numSamples - 1] > 0.0,
          "cancellation: completed run filled results");
  }

  // a firing poll stops the single-chain run and returns true, after
  // consulting the predicate at least once
  {
    SamplerOptions options;
    options.numTrees = 25;
    ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, options, &rng);
    int calls = 0;
    std::function<bool()> always = [&calls]() { ++calls; return true; };
    std::vector<double> sigmaDraws(100000, 0.0);
    Results results;
    results.sigma = sigmaDraws.data();
    bool cancelled = sampler.run(0, 100000, results, always);
    check(cancelled, "cancellation: single-chain poll stops the run");
    check(calls >= 1, "cancellation: poll consulted");
  }

  // multi-chain: the poll fires, every worker stops, and run returns true
  {
    const size_t numChains = 4;
    std::vector<ext_rng*> rngs(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], 2000 + static_cast<uint_least32_t>(c));
    }
    SamplerOptions options;
    options.numTrees = 25;
    options.numChains = numChains;
    options.numThreads = numChains;
    ClassicSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, options, rngs.data());
    std::atomic<int> calls(0);
    std::function<bool()> always = [&calls]() {
      calls.fetch_add(1);
      return true;
    };
    std::vector<double> sigmaDraws(numChains * 100000, 0.0);
    Results results;
    results.sigma = sigmaDraws.data();
    bool cancelled = sampler.run(0, 100000, results, always);
    check(cancelled, "cancellation: multi-chain poll stops all workers");
    for (size_t c = numChains; c > 0; --c) ext_rng_destroy(rngs[c - 1]);
  }

  printf("ok: run cancellation\n");
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

static void testWeightedLogistic(ext_rng* rng) {
  const size_t n = 400, p = 3;
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

  // integer-count weights: a count-w observation's Polya-Gamma latent is
  // PG(w, psi), drawn as the sum of w PG(1, psi) variates. Correctness is
  // pinned by determinism here and against the replicated-rows fit in
  // test-weighted-logistic.R. (This does NOT assert weight-1 == unweighted:
  // that equality holds through the engine, but comparing two separately-
  // constructed samplers here is subject to the codebase's heap-layout /
  // SIMD-reduction-split sensitivity, so a bitwise cross-sampler check is
  // unreliable in this harness regardless of weights.)
  std::vector<double> w(n);
  for (size_t i = 0; i < n; ++i) w[i] = static_cast<double>(1 + (i % 3));

  // the weighted path is deterministic: identical weights and seed reproduce
  // the draw bit for bit
  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 908) != 0 ||
      ext_rng_setSeed(rngB, 908) != 0) {
    check(false, "weighted logistic: rng creation");
    return;
  }
  ClassicSampler a(x.data(), y.data(), n, p, nullptr, w.data(),
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngA);
  ClassicSampler b(x.data(), y.data(), n, p, nullptr, w.data(),
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngB);
  const size_t numBurnIn = 20, numSamples = 30;
  std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
  Results rA, rB;
  rA.trainingFits = fitsA.data();
  rB.trainingFits = fitsB.data();
  a.run(numBurnIn, numSamples, rA);
  b.run(numBurnIn, numSamples, rB);
  bool identical = true;
  for (size_t i = 0; i < n * numSamples && identical; ++i)
    identical = fitsA[i] == fitsB[i];
  check(identical, "weighted logistic is deterministic under a fixed seed");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);

  // the weighted fit recovers the monotone signal and its omega latents stay
  // positive and finite
  ClassicSampler wsampler(x.data(), y.data(), n, p, nullptr, w.data(),
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rng);
  const size_t wSamples = 200;
  std::vector<double> wfits(n * wSamples);
  Results wres;
  wres.trainingFits = wfits.data();
  wsampler.run(150, wSamples, wres);
  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < wSamples; ++s)
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < 0.25) { lowSum += wfits[i + s * n]; ++lowCount; }
      if (x[i] > 0.75) { highSum += wfits[i + s * n]; ++highCount; }
    }
  check(highSum / (double) highCount > lowSum / (double) lowCount + 0.5,
        "weighted logistic recovers monotone signal");
  const double* omega = wsampler.latents(0);
  bool omegaValid = omega != nullptr;
  for (size_t i = 0; i < n && omegaValid; ++i)
    omegaValid = omega[i] > 0.0 && std::isfinite(omega[i]);
  check(omegaValid, "weighted logistic omega positive and finite");

  printf("ok: weighted logistic\n");
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
  Rule rootRule;  rootRule.variableIndex = 0;  rootRule.setSplitIndex(3);
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;  leftRule.variableIndex = 0;  leftRule.setSplitIndex(1);
  tree.birth(store, tree.at(0).leftChild, leftRule, y.data(), nullptr);
  int32_t leftChild = tree.at(0).leftChild;

  // new values 2..16 by twos: cuts {3, 5, ..., 15}; 4.5 is nearest 5 (index
  // 1), and 2.5 clamps into the left child's shrunken interval [0, 1)
  std::vector<double> x2(n);
  for (size_t i = 0; i < n; ++i) x2[i] = 2.0 * static_cast<double>(i + 1);
  store.setData(x2.data(), n);

  std::vector<double> params(tree.nodes.size(), 0.0);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex() == 1, "root split remaps to nearest cut");
  check(tree.at(leftChild).rule.splitIndex() == 0,
        "child split clamps into the ancestor-constrained interval");

  // shift the grid entirely above the old cuts: the root clamps to index 0,
  // leaving the left child no interval, so it collapses with plain-mean param
  std::vector<double> x3(n);
  for (size_t i = 0; i < n; ++i) x3[i] = 20.0 + 2.0 * static_cast<double>(i);
  oldCuts = store.cutPoints;
  int32_t oldRootIndex = tree.at(0).rule.splitIndex();
  check(oldRootIndex == 1, "");  // silence unused in release
  store.setData(x3.data(), n);

  params.assign(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(leftChild, bottoms);
  for (size_t k = 0; k < bottoms.size(); ++k)
    params[static_cast<size_t>(bottoms[k])] = static_cast<double>(k + 1);
  tree.mapOldCutPointsOntoNew(store, oldCuts, params);
  check(tree.at(0).rule.splitIndex() == 0, "root split clamps to the low end");
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
  rule.setCategoryDirections((1u << 1) | (1u << 3));
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
  checkNear(left.sumWeightedResponse / left.sumWeights, (0.0 + 2.0) / 2.0, 1e-12,
            "left average by category");

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
  childRule.setCategoryDirections(1u << 2);
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
    check(rule.categoryDirections() > 0 && rule.categoryDirections() < 31u + 1u &&
            rule.categoryDirections() != 31u,
          "categorical draw leaves neither side empty");
    ++patternCounts[rule.categoryDirections()];
  }
  bool uniform = patternCounts[0] == 0 && patternCounts[31] == 0;
  double expected = static_cast<double>(numDraws) / 30.0;
  for (std::uint32_t pattern = 1; pattern < 31; ++pattern)
    uniform = uniform &&
      std::fabs(patternCounts[pattern] - expected) < 5.0 * std::sqrt(expected);
  check(uniform, "categorical rule draw is uniform over valid assignments");

  Rule rootRule;
  rootRule.variableIndex = 0;
  rootRule.setCategoryDirections((1u << 3) | (1u << 4));
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -std::log(30.0), 1e-13, "root categorical rule probability");

  // left child reaches {0, 1, 2}: R = 3 gives 2^3 - 2 = 6
  Rule childRule;
  childRule.variableIndex = 0;
  childRule.setCategoryDirections(1u << 1);
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
  rule.setCategoryDirections((1ull << 1) | (1ull << 40) | (1ull << 52));
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
          rule.categoryDirections(),
        "right child reaches the mask's categories");

  // equals must compare the full mask width (the swap move's same-rule
  // test); ordinal rules zero the high word at construction
  Rule wideA, wideB;
  wideA.variableIndex = wideB.variableIndex = 0;
  wideA.setCategoryDirections(1ull | (1ull << 40));
  wideB.setCategoryDirections(1ull);
  check(!wideA.equals(wideB), "equals sees high mask bits");
  Rule ordinalA, ordinalB;
  ordinalA.variableIndex = ordinalB.variableIndex = 1;
  ordinalA.setSplitIndex(5);
  ordinalB.setSplitIndex(5);
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
    check(drawn.categoryDirections() > 0 &&
            drawn.categoryDirections() < (1ull << K) - 1ull,
          "wide draw leaves neither side empty");
    for (size_t bit = 0; bit < K; ++bit)
      bitCounts[bit] +=
        static_cast<int>((drawn.categoryDirections() >> bit) & 1);
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
  check(flat[0].value == static_cast<double>(rule.categoryDirections()),
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
  check(rebuilt.at(0).rule.categoryDirections() == (((1ull << K) - 1ull) & ~1ull),
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

static void testPooledMaskMechanics(ext_rng* rng) {
  // columns past 63 categories store per-tree pooled mask words; exercise
  // the partition kernel, reachability, draws, the pool's compaction, and
  // the flattened side channel (docs/design/pooled-masks.md)
  const size_t K = 70;  // ceil(71 / 64) = 2 words
  const size_t n = 10 * K;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    size_t category = (i * 37) % K;
    x[i] = static_cast<double>(category);
    y[i] = (category % 3 == 0 ? 2.0 : 0.0) + 0.3 * (runif01() - 0.5);
  }

  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);
  check(store.numCuts[0] == K, "pooled column counts its categories");
  check(store.columnIsPooled(0) && store.columnHasWideMask(0) &&
          store.hasPooledCategorical && store.hasWideCategorical,
        "pooled tier predicates");
  check(maskWordsForCount(K) == 2, "70 categories need two words");
  check(missingCategoryCode(K) == 70 &&
          store.codeFor(0, std::nan("")) == 70,
        "a pooled column's missing code is K");

  // a hand-built rule whose direction bits straddle both words
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  size_t offset = tree.allocateMask(2);
  {
    std::uint64_t* words = tree.mutableMaskWordsFor(offset);
    maskSetBit(words, 1);
    maskSetBit(words, 40);
    maskSetBit(words, 64);
    maskSetBit(words, 69);
  }
  Rule rule;
  rule.variableIndex = 0;
  rule.setMaskOffset(offset);
  tree.birth(store, 0, rule, y.data(), nullptr);
  const Node& left(tree.at(tree.at(0).leftChild));
  const Node& right(tree.at(tree.at(0).leftChild + 1));
  size_t expectedRight = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t category = (i * 37) % K;
    if (category == 1 || category == 40 || category == 64 || category == 69)
      ++expectedRight;
  }
  check(right.numObservations() == expectedRight &&
          left.numObservations() == n - expectedRight,
        "pooled mask partitions across words");
  bool sidesMatch = true;
  for (size_t i = right.begin; i < right.end; ++i) {
    size_t category = (tree.indices[i] * 37) % K;
    sidesMatch &= category == 1 || category == 40 || category == 64 ||
                  category == 69;
  }
  check(sidesMatch, "right side holds only right-bound categories");

  std::uint64_t reachable[2];
  tree.reachableCategoriesWide(store, tree.at(0).leftChild + 1, 0, reachable);
  check(maskEquals(reachable, tree.maskWordsFor(tree.at(0).rule), 2),
        "right child reaches the mask's categories");
  tree.reachableCategoriesWide(store, tree.at(0).leftChild, 0, reachable);
  check(maskPopcount(reachable, 2) == K - 4 && !maskTestBit(reachable, 40) &&
          maskTestBit(reachable, 0) && maskTestBit(reachable, 68),
        "left child reaches the complement");

  // each pooled direction bit is marginally 1/2 under the uniform prior on
  // assignments that leave neither side empty
  Tree drawTree;
  drawTree.initialize(indices.data(), n);
  CGMTreePrior prior;
  const int numDraws = 20000;
  std::vector<int> bitCounts(K, 0);
  bool neverEmpty = true;
  for (int d = 0; d < numDraws; ++d) {
    size_t mark = drawTree.maskPoolMark();
    Rule drawn = prior.drawRuleForVariable(drawTree, store, rng, 0, 0);
    const std::uint64_t* words = drawTree.maskWordsFor(drawn);
    size_t numRight = maskPopcount(words, 2);
    neverEmpty &= numRight > 0 && numRight < K;
    for (size_t bit = 0; bit < K; ++bit)
      bitCounts[bit] += maskTestBit(words, static_cast<std::uint32_t>(bit))
        ? 1 : 0;
    drawTree.truncateMaskPool(mark);
  }
  check(neverEmpty, "pooled draws leave neither side empty");
  bool marginsMatch = true;
  double tolerance = 5.0 * std::sqrt(0.25 * numDraws);
  for (size_t bit = 0; bit < K; ++bit)
    marginsMatch &= std::fabs(bitCounts[bit] - 0.5 * numDraws) < tolerance;
  check(marginsMatch, "pooled draw direction bits are marginally uniform");
  checkNear(prior.ruleForVariableLogProbability(tree, store, 0),
            -(70.0 * std::log(2.0) + std::log1p(-std::pow(2.0, -69.0))),
            1e-13, "pooled rule probability uses the closed form");

  // compaction: strand garbage past the high-water mark, then verify the
  // live rule's words survive at a fresh offset and still partition
  std::vector<std::uint64_t> liveMask(tree.maskWordsFor(tree.at(0).rule),
                                      tree.maskWordsFor(tree.at(0).rule) + 2);
  for (int g = 0; g < 200; ++g) tree.allocateMask(2);
  tree.compactMaskPoolIfNeeded(store);
  check(tree.maskPool.size() == 2, "compaction keeps only live masks");
  check(maskEquals(tree.maskWordsFor(tree.at(0).rule), liveMask.data(), 2),
        "compaction preserves mask contents");
  tree.repartitionSubtree(store, 0);
  check(tree.at(tree.at(0).leftChild + 1).numObservations() == expectedRight,
        "compacted rules still partition");

  // the flattened side channel: offsets sequential in pre-order, category
  // bits only, gauge re-checked on rebuild
  std::vector<FlatNode> flat;
  std::vector<std::uint64_t> flatMasks;
  std::vector<double> paramByNode(tree.nodes.size(), 0.0);
  tree.flatten(store, paramByNode.data(), flat, nullptr, 1, nullptr,
               &flatMasks);
  check(flat.size() == 3 && flatMasks.size() == 2 && flat[0].value == 0.0,
        "pooled flatten emits the side channel");
  check(maskEquals(flatMasks.data(), liveMask.data(), 2),
        "flattened words match the live mask");
  check(flatTreeIsWellFormed(store, flat.data(), flat.size(),
                             flatMasks.data(), flatMasks.size()),
        "pooled flat tree is well formed");
  check(!flatTreeIsWellFormed(store, flat.data(), flat.size()),
        "a missing mask channel is rejected");

  std::vector<size_t> rebuiltIndices(n);
  std::vector<double> rebuiltParams;
  Tree rebuilt;
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                              1, nullptr, flatMasks.data(), flatMasks.size()),
        "pooled tree rebuilds from flat");
  check(maskEquals(rebuilt.maskWordsFor(rebuilt.at(0).rule),
                   liveMask.data(), 2),
        "pooled mask round-trips exactly");

  flat[0].value = 1.0;  // offsets must be the running cursor
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(!rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                               1, nullptr, flatMasks.data(),
                               flatMasks.size()),
        "a non-sequential mask offset is rejected");
  flat[0].value = 0.0;
  std::vector<std::uint64_t> badMasks(flatMasks);
  maskSetBit(badMasks.data(), 70);  // the missing position must stay clear
  rebuilt.initialize(rebuiltIndices.data(), n);
  check(!rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams,
                               1, nullptr, badMasks.data(), badMasks.size()),
        "mask bits past the categories are rejected");

  // missing values compose: the NA pseudo-category is bit K
  std::vector<double> xMissing(x);
  xMissing[0] = std::nan("");
  xMissing[7] = std::nan("");
  ColumnStore storeMissing;
  storeMissing.build(xMissing.data(), n, 1, 10, false, types);
  check(storeMissing.numCuts[0] == K && storeMissing.hasMissing[0] == 1 &&
          storeMissing.codes[0] == 70,
        "pooled NA takes code K");
  Tree missingTree;
  missingTree.initialize(indices.data(), n);
  missingTree.reachableCategoriesWide(storeMissing, 0, 0, reachable);
  check(maskPopcount(reachable, 2) == K + 1 && maskTestBit(reachable, 70),
        "the missing pseudo-category enters the reachable set");

  printf("ok: pooled mask mechanics\n");
}

static void testPooledMaskSampler(ext_rng* rng) {
  // end to end over one pooled column (K = 70) and one inline-band column
  // (K = 60, whose flattened masks also move to the side channel), with
  // keepTrees, saved-tree replay, live prediction, and a bitwise state
  // round trip mid-run
  const size_t K = 70, K2 = 60;
  const size_t n = 10 * K;
  std::vector<double> x(2 * n), y(n);
  for (size_t i = 0; i < n; ++i) {
    size_t category = (i * 37) % K;
    size_t category2 = (i * 13) % K2;
    x[i] = static_cast<double>(category);
    x[i + n] = static_cast<double>(category2);
    y[i] = (category % 3 == 0 ? 2.0 : 0.0) +
           (category2 >= 55 ? -1.5 : 0.0) + 0.3 * (runif01() - 0.5);
  }

  ColumnType types[] = {ColumnType::categorical, ColumnType::categorical};
  {
    ColumnStore probe;
    probe.build(x.data(), n, 2, 10, false, types);
    check(!probe.columnIsPooled(1) && probe.columnHasWideMask(1),
          "60 categories stay inline but flatten wide");
  }

  const size_t numSamples = 4, nTest = K;
  std::vector<double> xTest(2 * nTest);
  for (size_t i = 0; i < nTest; ++i) {
    xTest[i] = static_cast<double>(i % K);
    xTest[i + nTest] = static_cast<double>(i % K2);
  }

  SamplerOptions options;
  options.numTrees = 50;
  options.columnTypes = types;
  options.keepTrees = true;
  options.numSamplesToStore = numSamples;
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  sampler.setTestPredictors(xTest.data(), nTest);

  const size_t numRecorded = 100;
  std::vector<double> trainingFits(n * numRecorded);
  Results burnResults;
  burnResults.trainingFits = trainingFits.data();
  sampler.run(100, numRecorded, burnResults);

  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < numRecorded; ++s)
    for (size_t i = 0; i < n; ++i) {
      double centered = trainingFits[i + s * n] -
                        ((i * 13) % K2 >= 55 ? -1.5 : 0.0);
      if ((i * 37) % K % 3 == 0) { highSum += centered; ++highCount; }
      else                       { lowSum += centered; ++lowCount; }
    }
  check(std::fabs(highSum / static_cast<double>(highCount) - 2.0) < 0.2 &&
          std::fabs(lowSum / static_cast<double>(lowCount)) < 0.2,
        "pooled subset splits recover group means");

  // saved-tree replay equals the recorded test fits exactly
  std::vector<double> sigma(numSamples), testFits(nTest * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.testFits = testFits.data();
  sampler.run(0, numSamples, results);
  std::vector<double> predicted(nTest * numSamples);
  sampler.predict(xTest.data(), nTest, predicted.data());
  check(predicted == testFits,
        "pooled saved-tree predictions equal the run's test fits");

  // a stored state restores to a bitwise-identical continuation
  SamplerStateData state;
  sampler.getState(state);
  check(!state.chains[0].treeMasks.empty(),
        "wide-column states carry mask channels");
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 4242);
  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rng2);
  restored.setTestPredictors(xTest.data(), nTest);
  check(restored.setState(state), "a pooled state restores");

  std::vector<double> sigmaA(numSamples), trainA(n * numSamples);
  std::vector<double> sigmaB(numSamples), trainB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  sampler.run(0, numSamples, resultsA);
  restored.run(0, numSamples, resultsB);
  check(sigmaA == sigmaB && trainA == trainB,
        "restored pooled chains continue bitwise");

  // live-tree prediction agrees with the final tree fits
  std::vector<double> livePredict(n);
  {
    ClassicSampler live(x.data(), y.data(), n, 2, nullptr, nullptr,
                        ResponseFamily::gaussian, 1.0, 3.0,
                        0.37804942330213542, options, &rng2);
    Results liveResults;
    std::vector<double> liveSigma(1), liveTrain(n);
    liveResults.sigma = liveSigma.data();
    liveResults.trainingFits = liveTrain.data();
    live.setTreeStorage(false, 0);
    live.run(50, 1, liveResults);
    live.predict(x.data(), n, livePredict.data());
    bool livePredictionMatches = true;
    for (size_t i = 0; i < n; ++i)
      livePredictionMatches &=
        std::fabs(livePredict[i] - liveTrain[i]) < 1e-10;
    check(livePredictionMatches,
          "live pooled prediction matches the last recorded fits");
  }
  ext_rng_destroy(rng2);

  printf("ok: pooled mask sampler\n");
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
  rootRule.setCategoryDirections(0x6);  // categories 1, 2 right
  tree.birth(store, 0, rootRule, y.data(), nullptr);
  Rule leftRule;
  leftRule.variableIndex = 1;
  leftRule.setSplitIndex(4);
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
        tree2.at(0).rule.categoryDirections() == 0x6,
        "buildFromFlat recovers the categorical mask");
  check(tree2.at(tree2.at(0).leftChild).rule.splitIndex() == 4,
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


static void testMissingIngestion() {
  const size_t n = 100;
  double na = std::nan("");
  // column 0: ordinal with NAs; column 1: categorical with NAs; column 2:
  // complete ordinal
  std::vector<double> x(n * 3);
  for (size_t i = 0; i < n; ++i) {
    x[i] = i % 10 == 0 ? na : runif01();
    x[i + n] = i % 7 == 0 ? na : static_cast<double>(i % 4);
    x[i + 2 * n] = runif01();
  }

  ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical,
                        ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, 3, 10, false, types);

  check(store.hasMissing[0] == 1 && store.hasMissing[1] == 1 &&
          store.hasMissing[2] == 0,
        "hasMissing marks exactly the columns with NAs");
  check(store.numCuts[1] == 4,
        "categorical categories are counted over observed values");

  bool codesRight = true;
  for (size_t i = 0; i < n; ++i) {
    codesRight &= (store.codes[i] == naCode) == (i % 10 == 0);
    codesRight &=
      (store.codes[i + n] == static_cast<xint_t>(naCategory)) == (i % 7 == 0);
  }
  check(codesRight, "missing values take the reserved codes");

  double loObserved = 2.0, hiObserved = -1.0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 0) continue;
    if (x[i] < loObserved) loObserved = x[i];
    if (x[i] > hiObserved) hiObserved = x[i];
  }
  check(store.cutPoints[0].front() > loObserved &&
          store.cutPoints[0].back() < hiObserved,
        "uniform cuts span the observed range only");

  ColumnStore quantileStore;
  quantileStore.build(x.data(), n, 3, 10, true, types);
  check(quantileStore.numCuts[0] == 10 && quantileStore.hasMissing[0] == 1,
        "quantile grids skip NaN and keep the requested count");

  check(store.categoricalValueIsValid(1, na),
        "a missing categorical value is representable");
  printf("ok: missing ingestion\n");
}

static void testMissingMechanics() {
  const size_t n = 90;
  double na = std::nan("");
  std::vector<double> x(n * 2), y(n);
  size_t numMissingOrdinal = 0, numMissingCategorical = 0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 3 == 0) { x[i] = na; ++numMissingOrdinal; }
    else x[i] = runif01();
    if (i % 5 == 0) { x[i + n] = na; ++numMissingCategorical; }
    else x[i + n] = static_cast<double>(i % 4);
    y[i] = runif01();
  }
  ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 2, 10, false, types);

  // an ordinal rule routes the reserved code by its missing direction
  std::vector<size_t> indices(n);
  Tree tree;
  tree.initialize(indices.data(), n);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(4);
  rule.setMissingGoesRight(true);
  check(rule.splitIndex() == 4 && rule.missingGoesRight(),
        "rule accessors pack the split index and missing direction");
  tree.birth(store, 0, rule, y.data(), nullptr);

  const Node& right(tree.at(tree.at(0).leftChild + 1));
  size_t missingOnRight = 0;
  for (size_t i = right.begin; i < right.end; ++i)
    if (store.codes[tree.indices[i]] == naCode) ++missingOnRight;
  const Node& left(tree.at(tree.at(0).leftChild));
  size_t missingOnLeft = 0;
  for (size_t i = left.begin; i < left.end; ++i)
    if (store.codes[tree.indices[i]] == naCode) ++missingOnLeft;
  check(missingOnRight == numMissingOrdinal && missingOnLeft == 0,
        "missing ordinal rows all route by the rule's direction");

  std::vector<size_t> indices2(n);
  Tree treeLeft;
  treeLeft.initialize(indices2.data(), n);
  Rule ruleLeft;
  ruleLeft.variableIndex = 0;
  ruleLeft.setSplitIndex(4);
  treeLeft.birth(store, 0, ruleLeft, y.data(), nullptr);
  const Node& leftLeft(treeLeft.at(treeLeft.at(0).leftChild));
  size_t missingOnLeft2 = 0;
  for (size_t i = leftLeft.begin; i < leftLeft.end; ++i)
    if (store.codes[treeLeft.indices[i]] == naCode) ++missingOnLeft2;
  check(missingOnLeft2 == numMissingOrdinal,
        "the canonical zero direction sends missing rows left");

  // a categorical rule treats missing as one more category: the reachable
  // mask seeds it and a rule can isolate it
  std::vector<size_t> indices3(n);
  Tree catTree;
  catTree.initialize(indices3.data(), n);
  check(catTree.reachableCategories(store, 0, 1) ==
          (0xfull | Rule::missingDirectionBit),
        "the reachable mask includes the missing category");
  Rule catRule;
  catRule.variableIndex = 1;
  catRule.setCategoryDirections(Rule::missingDirectionBit);
  catTree.birth(store, 0, catRule, y.data(), nullptr);
  const Node& catRight(catTree.at(catTree.at(0).leftChild + 1));
  check(catRight.numObservations() == numMissingCategorical,
        "a mask can send exactly the missing rows down one side");
  check(catTree.reachableCategories(store, catTree.at(0).leftChild, 1) ==
          0xfull &&
        catTree.reachableCategories(store, catTree.at(0).leftChild + 1, 1) ==
          Rule::missingDirectionBit,
        "children inherit the filtered missing category");

  // the flat format moves the missing direction out of the value
  std::vector<double> params(tree.nodes.size(), 0.0);
  std::vector<FlatNode> flat;
  tree.flatten(store, params.data(), flat);
  check(flat[0].flags == flatMissingGoesRight &&
          flat[0].value == store.cutPoints[0][4],
        "an ordinal flat node carries its missing direction in flags");

  std::vector<double> catParams(catTree.nodes.size(), 0.0);
  std::vector<FlatNode> catFlat;
  catTree.flatten(store, catParams.data(), catFlat);
  check(catFlat[0].value == 0.0 && catFlat[0].flags == flatMissingGoesRight,
        "a missing-only mask flattens to value zero plus the flag");
  check(flatTreeIsWellFormed(store, catFlat.data(), catFlat.size()),
        "the missing-only mask is well formed");

  std::vector<size_t> indices4(n);
  Tree rebuilt;
  rebuilt.initialize(indices4.data(), n);
  std::vector<double> rebuiltParams;
  check(rebuilt.buildFromFlat(store, catFlat.data(), catFlat.size(),
                              rebuiltParams),
        "a flat tree with a missing direction rebuilds");
  check(rebuilt.at(0).rule.categoryDirections() == Rule::missingDirectionBit,
        "the rebuilt rule recovers the missing direction");

  // raw-value replay routes NaN by the flag
  std::vector<std::uint32_t> counts(flat.size());
  std::vector<size_t> replayIndices(n);
  for (size_t i = 0; i < n; ++i) replayIndices[i] = i;
  countFlatObservationsBelow(flat.data(), store.types.data(), x.data(), n,
                             replayIndices.data(), 0, n, counts.data());
  size_t expectedRight = 0;
  for (size_t i = 0; i < n; ++i)
    if (isNA(x[i]) || x[i] > store.cutPoints[0][4]) ++expectedRight;
  check(counts[2] == expectedRight,
        "replay against raw values routes NaN by the stored direction");

  FlatNode bad = catFlat[0];
  bad.flags = 0x2;
  std::vector<FlatNode> badTree(catFlat);
  badTree[0] = bad;
  check(!flatTreeIsWellFormed(store, badTree.data(), badTree.size()),
        "unknown flag bits are malformed");
  std::vector<FlatNode> badLeaf(catFlat);
  badLeaf[1].flags = flatMissingGoesRight;
  check(!flatTreeIsWellFormed(store, badLeaf.data(), badLeaf.size()),
        "a flagged leaf is malformed");
  printf("ok: missing mechanics\n");
}

static void testMissingEndToEnd() {
  // a fit on data with NAs learns a signal carried by missingness, routes
  // NaN test rows, and its state round-trips bitwise
  const size_t n = 300, numSamples = 4;
  double na = std::nan("");
  std::vector<double> x(n * 2), y(n);
  std::vector<bool> missing(n);
  for (size_t i = 0; i < n; ++i) {
    missing[i] = runif01() < 0.3;
    x[i] = missing[i] ? na : runif01();
    x[i + n] = runif01();
    y[i] = (missing[i] ? 2.0 : 0.0) + 0.25 * x[i + n] +
           0.1 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.numChains = 1;
  options.keepTrees = true;
  options.numSamplesToStore = numSamples;

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 5150);
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rng);
  std::vector<double> sigma(numSamples), train(n * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.trainingFits = train.data();
  sampler.run(150, numSamples, results);

  double missingMean = 0.0, observedMean = 0.0;
  size_t numMissing = 0;
  const double* lastFits = train.data() + n * (numSamples - 1);
  for (size_t i = 0; i < n; ++i) {
    if (missing[i]) { missingMean += lastFits[i]; ++numMissing; }
    else observedMean += lastFits[i];
  }
  missingMean /= static_cast<double>(numMissing);
  observedMean /= static_cast<double>(n - numMissing);
  check(missingMean - observedMean > 1.0,
        "the fit recovers a signal carried by missingness");

  // saved-tree prediction routes NaN rows to the missing side
  double xTest[] = {na, 0.5, 0.5, 0.5};  // column-major, two rows
  size_t capacity = sampler.savedTreeCapacity();
  std::vector<double> predictions(2 * capacity);
  sampler.predict(xTest, 2, predictions.data());
  double missingFit = 0.0, observedFit = 0.0;
  for (size_t k = 0; k < capacity; ++k) {
    missingFit += predictions[2 * k];
    observedFit += predictions[2 * k + 1];
  }
  check(missingFit / static_cast<double>(capacity) -
          observedFit / static_cast<double>(capacity) > 1.0,
        "prediction routes NaN rows by the learned directions");

  // bitwise state round trip carries the missing directions along
  SamplerStateData state;
  sampler.getState(state);
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 999);
  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rng2);
  check(restored.setState(state), "a state with missing directions restores");

  std::vector<double> sigmaA(numSamples), trainA(n * numSamples);
  Results resultsA;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  sampler.run(0, numSamples, resultsA);

  std::vector<double> sigmaB(numSamples), trainB(n * numSamples);
  Results resultsB;
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  restored.run(0, numSamples, resultsB);

  check(sigmaA == sigmaB && trainA == trainB,
        "restored chains continue bitwise with missing data");

  ext_rng_destroy(rng);
  ext_rng_destroy(rng2);
  printf("ok: missing end to end\n");
}

// Shared fixture for the linear-leaf component tests; the reference
// constants come from an independent R implementation
// (scratchpad linear_leaf_reference.R).
struct LinearLeafFixture {
  static constexpr size_t n = 8;
  std::vector<double> x, z, w;
  ColumnStore store;
  std::vector<size_t> indexBuffer;
  Tree tree;
  double k = 2.0, scale = 0.5 / std::sqrt(10.0), sigmaSq = 0.04;

  LinearLeafFixture() : indexBuffer(n) {
    double x1[] = {0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9};
    double u1[] = {1.3, -0.4, 0.7, 2.2, -1.5, 0.6, 0.1, -0.8};
    double u2[] = {0.2, 1.1, -0.6, 0.4, 0.9, -1.3, 0.5, -0.2};
    x.assign(3 * n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      x[i] = x1[i];
      x[i + n] = u1[i];
      x[i + 2 * n] = u2[i];
    }
    z = {0.35, -0.10, 0.22, 0.55, -0.42, 0.16, 0.03, -0.25};
    w = {1.0, 2.0, 0.5, 1.0, 1.5, 1.0, 0.8, 1.2};
    store.build(x.data(), n, 3, 100);
    tree.initialize(indexBuffer.data(), n);
  }
};

static void testLinearLeafMarginal() {
  LinearLeafFixture f;

  LinearGaussianLeaf leaf;
  leaf.scale = f.scale;
  size_t columns[] = {1};
  leaf.initialize(f.store, columns, 1);
  check(leaf.numParams() == 2, "linear leaf parameter count");

  // root marginals against the R reference
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, 0),
            -5.517188945419923, 1e-9, "linear marginal, weighted root");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), nullptr,
                                                f.k, f.sigmaSq, 0),
            -5.256174278535296, 1e-9, "linear marginal, unit-weight root");

  // q = 0 reduces exactly to the constant leaf's formula
  LinearGaussianLeaf interceptOnly;
  interceptOnly.scale = f.scale;
  interceptOnly.initialize(f.store, nullptr, 0);
  ConstantGaussianLeaf constant{f.scale};
  f.tree.setNodeAverages(f.z.data(), f.w.data());
  checkNear(interceptOnly.logIntegratedLikelihoodForNode(
              f.tree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0),
            constant.logIntegratedLikelihoodForNode(
              f.tree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0),
            1e-9, "q = 0 linear marginal equals the constant leaf's");

  // children of a split at x1 ~ 0.5 (cut 50 of the uniform grid)
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(50);
  f.tree.birth(f.store, 0, rule, f.z.data(), f.w.data());
  int32_t leftChild = f.tree.at(0).leftChild;
  check(f.tree.at(leftChild).numObservations() == 4 &&
        f.tree.at(leftChild + 1).numObservations() == 4,
        "fixture split partitions 4/4");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild),
            -3.760266419465231, 1e-9, "linear marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild + 1),
            -3.062089185586155, 1e-9, "linear marginal, right child");

  // two covariates exercise the 3x3 Cholesky
  LinearGaussianLeaf leaf2;
  leaf2.scale = f.scale;
  size_t columns2[] = {1, 2};
  leaf2.initialize(f.store, columns2, 2);
  Tree rootTree;
  std::vector<size_t> rootIndices(LinearLeafFixture::n);
  rootTree.initialize(rootIndices.data(), LinearLeafFixture::n);
  checkNear(leaf2.logIntegratedLikelihoodForNode(
              rootTree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0),
            -5.791887259328131, 1e-9, "linear marginal, two covariates");

  // standardization constants match R's mean/sd
  checkNear(leaf2.covariateMeans()[0], 0.275, 1e-12, "covariate mean");
  checkNear(leaf2.covariateSds()[0], 1.185326959112970, 1e-12, "covariate sd");

  printf("ok: linear leaf marginal\n");
}

static void testLinearLeafDraw(ext_rng* rng) {
  LinearLeafFixture f;
  LinearGaussianLeaf leaf;
  leaf.scale = f.scale;
  size_t columns[] = {1};
  leaf.initialize(f.store, columns, 1);

  // posterior moments for the weighted root, from the R reference
  const double expectedMean[2] = {0.023101719628228, 0.176899579133385};
  const double expectedCov[3] = {0.002628477539281, 0.000290149914052,
                                 0.002709159455639};

  const int numDraws = 200000;
  double sum[2] = {0.0, 0.0};
  double sumSq[3] = {0.0, 0.0, 0.0};
  double draw[2];
  for (int r = 0; r < numDraws; ++r) {
    leaf.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(), f.k,
                                  f.sigmaSq, 0, draw);
    sum[0] += draw[0];
    sum[1] += draw[1];
    sumSq[0] += draw[0] * draw[0];
    sumSq[1] += draw[0] * draw[1];
    sumSq[2] += draw[1] * draw[1];
  }
  double mean0 = sum[0] / numDraws, mean1 = sum[1] / numDraws;
  checkNear(mean0, expectedMean[0], 1e-3, "posterior draw mean, intercept");
  checkNear(mean1, expectedMean[1], 1e-3, "posterior draw mean, slope");
  checkNear(sumSq[0] / numDraws - mean0 * mean0, expectedCov[0], 1e-4,
            "posterior draw variance, intercept");
  checkNear(sumSq[1] / numDraws - mean0 * mean1, expectedCov[1], 1e-4,
            "posterior draw covariance");
  checkNear(sumSq[2] / numDraws - mean1 * mean1, expectedCov[2], 1e-4,
            "posterior draw variance, slope");

  // an empty leaf zeroes its block without consuming generator draws; no
  // valid rule can empty a child of this fixture, so fabricate the range
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(50);
  f.tree.birth(f.store, 0, rule, f.z.data(), f.w.data());
  int32_t rightChild = f.tree.at(0).leftChild + 1;
  f.tree.at(rightChild).begin = f.tree.at(rightChild).end;
  check(f.tree.at(rightChild).numObservations() == 0, "empty right child");
  draw[0] = draw[1] = 1.0;
  leaf.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(), f.k,
                                f.sigmaSq, rightChild, draw);
  check(draw[0] == 0.0 && draw[1] == 0.0, "empty leaf draws a zero block");

  printf("ok: linear leaf draw\n");
}

static void testLinearLeafEndToEnd(ext_rng* rng) {
  const size_t n = 400, p = 2;
  std::vector<double> x(n * p), f(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    f[i] = x[i] > 0.5 ? x[i + n] : 0.0;
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = f[i] + 0.2 * normal;
  }

  double yMean = 0.0, ySd = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= (double) n;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / (double) (n - 1));

  SamplerOptions options;
  options.numTrees = 30;
  size_t covariates[] = {1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;

  // invalid designations refuse cleanly at the factory
  {
    SamplerOptions bad = options;
    size_t outOfRange[] = {7};
    bad.leafCovariateColumns = outOfRange;
    check(createSampler(x.data(), y.data(), n, p, nullptr, nullptr,
                        ResponseFamily::gaussian, ySd, 3.0,
                        0.37804942330213542, bad, &rng) == nullptr,
          "out-of-range leaf covariate refused");
    ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical};
    bad = options;
    bad.columnTypes = types;
    check(createSampler(x.data(), y.data(), n, p, nullptr, nullptr,
                        ResponseFamily::gaussian, ySd, 3.0,
                        0.37804942330213542, bad, &rng) == nullptr,
          "categorical leaf covariate refused");
  }

  std::unique_ptr<SamplerBase> sampler = createSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);
  check(sampler != nullptr, "linear-leaf sampler creation");

  // the first five training rows double as test rows: their fits must agree
  const size_t numTest = 5;
  std::vector<double> xTest(numTest * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < numTest; ++i)
      xTest[i + j * numTest] = x[i + j * n];
  sampler->setTestPredictors(xTest.data(), numTest);

  const size_t numBurnIn = 200, numSamples = 300;
  std::vector<double> sigmaDraws(numSamples);
  std::vector<double> trainingFits(n * numSamples);
  std::vector<double> testFits(numTest * numSamples);
  Results results;
  results.sigma = sigmaDraws.data();
  results.trainingFits = trainingFits.data();
  results.testFits = testFits.data();
  sampler->run(numBurnIn, numSamples, results);

  std::vector<double> posteriorMean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i)
      posteriorMean[i] += trainingFits[i + s * n] / (double) numSamples;

  double sseFit = 0.0, sseMean = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sseFit += (posteriorMean[i] - f[i]) * (posteriorMean[i] - f[i]);
    sseMean += (yMean - f[i]) * (yMean - f[i]);
  }
  check(sseFit < 0.2 * sseMean, "linear end-to-end: fit explains most signal");

  double sigmaPosteriorMean = 0.0;
  for (double s : sigmaDraws) sigmaPosteriorMean += s / (double) numSamples;
  check(sigmaPosteriorMean > 0.15 && sigmaPosteriorMean < 0.3,
        "linear end-to-end: sigma near truth (0.2)");

  bool testFitsMatch = true;
  for (size_t s = 0; s < numSamples && testFitsMatch; ++s)
    for (size_t i = 0; i < numTest && testFitsMatch; ++i)
      testFitsMatch = std::fabs(testFits[i + s * numTest] -
                                trainingFits[i + s * n]) < 1e-8;
  check(testFitsMatch, "test fits of duplicated rows match training fits");

  // the slope structure: within x1 > 0.5, fits track x2 with slope ~1;
  // within x1 < 0.5, slope ~0 (simple regression of fits on x2)
  double slope[2];
  for (int side = 0; side < 2; ++side) {
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, count = 0.0;
    for (size_t i = 0; i < n; ++i) {
      if ((x[i] > 0.5) != (side == 1)) continue;
      sx += x[i + n];
      sy += posteriorMean[i];
      sxx += x[i + n] * x[i + n];
      sxy += x[i + n] * posteriorMean[i];
      count += 1.0;
    }
    slope[side] = (sxy - sx * sy / count) / (sxx - sx * sx / count);
  }
  check(std::fabs(slope[0]) < 0.25, "flat side has near-zero slope");
  check(slope[1] > 0.75 && slope[1] < 1.25, "steep side recovers unit slope");

  // the format-dependent surfaces work post-stage 3: an identical
  // predictor matrix revalidates, and a state round-trips into itself
  check(sampler->setPredictor(x.data(), false, true) ==
          PredictorUpdateResult::accepted,
        "reinstalling identical predictors is accepted");
  SamplerStateData state;
  sampler->getState(state);
  check(sampler->setState(state), "a state restores into its own sampler");
  check(sampler->savedTreeCapacity() == 0, "keepTrees stays off by default");

  // a short run with the k hyperprior exercises the all-coordinate
  // accumulation
  SamplerOptions kOptions = options;
  kOptions.numTrees = 20;
  kOptions.updateK = true;
  std::unique_ptr<SamplerBase> kSampler = createSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, kOptions, &rng);
  const size_t numKSamples = 100;
  std::vector<double> kDraws(numKSamples), kSigma(numKSamples);
  Results kResults;
  kResults.k = kDraws.data();
  kResults.sigma = kSigma.data();
  kSampler->run(100, numKSamples, kResults);
  bool kSane = true;
  for (double kDraw : kDraws)
    kSane = kSane && std::isfinite(kDraw) && kDraw > 0.05 && kDraw < 50.0;
  check(kSane, "sampled k stays finite and positive");

  printf("ok: linear leaf end to end (sse ratio %.3f, sigma %.3f, "
         "slopes %.2f/%.2f)\n",
         sseFit / sseMean, sigmaPosteriorMean, slope[0], slope[1]);
}

static void testLinearLeafFormats(ext_rng* rng) {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = (x[i] > 0.5 ? x[i + n] : 0.0) + 0.2 * normal;
  }

  SamplerOptions options;
  options.numTrees = 20;
  options.keepTrees = true;
  options.numSamplesToStore = 25;
  size_t covariates[] = {1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;

  std::unique_ptr<SamplerBase> sampler = createSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    1.0, 3.0, 0.37804942330213542, options, &rng);

  const size_t numTest = 6;
  std::vector<double> xTest(numTest * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < numTest; ++i)
      xTest[i + j * numTest] = x[i + j * n];
  sampler->setTestPredictors(xTest.data(), numTest);

  const size_t numBurnIn = 100, numSamples = 25;
  std::vector<double> trainingFits(n * numSamples), testFits(numTest * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  results.testFits = testFits.data();
  sampler->run(numBurnIn, numSamples, results);

  // saved-tree replay reproduces the test fits recorded during the run
  check(sampler->savedTreeCapacity() == numSamples, "saved-tree capacity");
  std::vector<double> predictions(numTest * numSamples);
  sampler->predict(xTest.data(), numTest, predictions.data());
  bool replayMatches = true;
  for (size_t s = 0; s < numSamples && replayMatches; ++s)
    for (size_t i = 0; i < numTest && replayMatches; ++i)
      replayMatches = std::fabs(predictions[i + s * numTest] -
                                testFits[i + s * numTest]) < 1e-10;
  check(replayMatches, "saved-tree replay matches recorded test fits");

  // flattened trees carry one slope per leaf, and the saved slopes agree
  std::vector<FlatNode> nodes;
  std::vector<uint32_t> counts;
  std::vector<double> slopes;
  sampler->flattenTree(0, 0, nodes, counts, &slopes);
  check(slopes.size() == (nodes.size() + 1) / 2,
        "flattened tree carries one slope per leaf");
  const std::vector<FlatNode>& saved(sampler->savedTree(0, numSamples - 1, 0));
  const std::vector<double>& savedSlopes(
    sampler->savedTreeSlopes(0, numSamples - 1, 0));
  check(savedSlopes.size() == (saved.size() + 1) / 2,
        "saved tree carries one slope per leaf");

  // bitwise state round trip through the flat + slopes format
  SamplerStateData state;
  sampler->getState(state);
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 4321);
  std::unique_ptr<SamplerBase> restored = createSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    1.0, 3.0, 0.37804942330213542, options, &rng2);
  check(restored->setState(state), "a linear-leaf state restores");
  restored->setTestPredictors(xTest.data(), numTest);

  std::vector<double> sigmaA(numSamples), trainA(n * numSamples);
  Results resultsA;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  sampler->run(0, numSamples, resultsA);

  std::vector<double> sigmaB(numSamples), trainB(n * numSamples);
  Results resultsB;
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  restored->run(0, numSamples, resultsB);

  check(sigmaA == sigmaB && trainA == trainB,
        "restored linear-leaf chains continue bitwise");

  // a state with mismatched slope shapes is refused
  SamplerStateData malformed = state;
  malformed.chains[0].treeParams[0].push_back(0.0);
  check(!sampler->setState(malformed),
        "a state with mismatched slopes is refused");

  ext_rng_destroy(rng2);
  printf("ok: linear leaf formats\n");
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

static void testLinearLeafViews() {
  const size_t n = 400, p = 3;
  std::vector<double> x(n * p), f(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = 2.0 * runif01() - 1.0;
    x[i + 2 * n] = runif01();
    f[i] = x[i] > 0.5 ? x[i + n] : 0.0;
    double u1 = runif01(), u2 = runif01();
    double normal =
      std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    y[i] = f[i] + 0.2 * normal;
  }

  SamplerOptions options;
  options.numTrees = 25;
  size_t covariates[] = {1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;

  ColumnStore parent;
  parent.build(x.data(), n, p, options.maxNumCuts);

  // view gather mechanics: subset raw values and parent-derived constants
  {
    std::vector<size_t> rows, testRows;
    for (size_t i = 0; i < n; ++i)
      (i % 4 == 0 ? testRows : rows).push_back(i);
    ColumnStore view;
    view.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                         testRows.size(), covariates, 1);

    const double* raw = view.rawColumn(1);
    const double* rawTest = view.rawTestColumn(1);
    bool valuesMatch = raw != nullptr && rawTest != nullptr;
    for (size_t i = 0; i < rows.size() && valuesMatch; ++i)
      valuesMatch = raw[i] == x[rows[i] + n];
    for (size_t i = 0; i < testRows.size() && valuesMatch; ++i)
      valuesMatch = rawTest[i] == x[testRows[i] + n];
    check(valuesMatch, "view gathers raw leaf covariates by row");
    check(view.rawColumn(0) == nullptr && view.rawTestColumn(0) == nullptr,
          "undesignated columns stay unserved on views");

    double mean, sd, parentMean, parentSd;
    standardizationMomentsForColumn(x.data() + n, n, &parentMean, &parentSd);
    check(view.suppliedStandardization(1, &mean, &sd) &&
            mean == parentMean && sd == parentSd,
          "view standardization constants come from the parent's full data");
    check(!parent.suppliedStandardization(1, &mean, &sd),
          "raw-data stores supply no gathered constants");
  }

  // a full-rows linear view runs bitwise identically to the raw-data path
  {
    ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 1234) != 0 ||
        ext_rng_setSeed(rngB, 1234) != 0) {
      check(false, "linear view: rng creation");
      return;
    }

    std::unique_ptr<SamplerBase> full = createSampler(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngA);

    std::vector<size_t> rows(n);
    for (size_t i = 0; i < n; ++i) rows[i] = i;
    ColumnStore store;
    store.buildFromParent(parent, rows.data(), n, nullptr, 0, covariates, 1);
    std::unique_ptr<SamplerBase> view = createSamplerOverStore(
      std::move(store), y.data(), nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, options, &rngB);
    check(view != nullptr, "linear view sampler creation");

    const size_t numBurnIn = 30, numSamples = 40;
    std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
    std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
    Results resultsA, resultsB;
    resultsA.sigma = sigmaA.data();
    resultsA.trainingFits = fitsA.data();
    resultsB.sigma = sigmaB.data();
    resultsB.trainingFits = fitsB.data();
    full->run(numBurnIn, numSamples, resultsA);
    view->run(numBurnIn, numSamples, resultsB);

    bool identical = true;
    for (size_t s = 0; s < numSamples && identical; ++s)
      identical = sigmaA[s] == sigmaB[s];
    for (size_t i = 0; i < n * numSamples && identical; ++i)
      identical = fitsA[i] == fitsB[i];
    check(identical,
          "full-rows linear view bitwise-matches the raw-data path");

    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
  }

  // a proper fold recovers the varying slope on its held-out rows, and raw
  // test predictors installed later supersede the gathered fold values
  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rng == NULL || ext_rng_setSeed(rng, 5678) != 0) {
      check(false, "linear fold view: rng creation");
      return;
    }

    std::vector<size_t> rows, testRows;
    for (size_t i = 0; i < n; ++i)
      (i % 4 == 0 ? testRows : rows).push_back(i);
    size_t numTest = testRows.size();
    std::vector<double> yTrain(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) yTrain[i] = y[rows[i]];

    ColumnStore store;
    store.buildFromParent(parent, rows.data(), rows.size(), testRows.data(),
                          numTest, covariates, 1);
    std::unique_ptr<SamplerBase> view = createSamplerOverStore(
      std::move(store), yTrain.data(), nullptr, nullptr,
      ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options, &rng);
    check(view != nullptr, "linear fold view creation");

    const size_t numBurnIn = 200, numSamples = 200;
    std::vector<double> testFits(numTest * numSamples);
    Results results;
    results.testFits = testFits.data();
    view->run(numBurnIn, numSamples, results);

    std::vector<double> posteriorMean(numTest, 0.0);
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t i = 0; i < numTest; ++i)
        posteriorMean[i] += testFits[i + s * numTest] / (double) numSamples;
    double sseFit = 0.0, sseMean = 0.0;
    for (size_t i = 0; i < numTest; ++i) {
      double truth = f[testRows[i]];
      sseFit += (posteriorMean[i] - truth) * (posteriorMean[i] - truth);
      sseMean += truth * truth;
    }
    check(sseFit < 0.35 * sseMean,
          "fold view recovers held-out signal through the linear leaf");

    // duplicated training rows through setTestPredictors: raw values now
    // come from x_test, and their fits must match the training fits
    const size_t numDup = 5;
    std::vector<double> xDup(numDup * p);
    for (size_t j = 0; j < p; ++j)
      for (size_t i = 0; i < numDup; ++i)
        xDup[i + j * numDup] = x[rows[i] + j * n];
    view->setTestPredictors(xDup.data(), numDup);

    std::vector<double> dupFits(numDup * numSamples);
    std::vector<double> trainingFits(rows.size() * numSamples);
    Results dupResults;
    dupResults.testFits = dupFits.data();
    dupResults.trainingFits = trainingFits.data();
    view->run(0, numSamples, dupResults);
    bool dupMatches = true;
    for (size_t s = 0; s < numSamples && dupMatches; ++s)
      for (size_t i = 0; i < numDup && dupMatches; ++i)
        dupMatches = std::fabs(dupFits[i + s * numDup] -
                               trainingFits[i + s * rows.size()]) < 1e-8;
    check(dupMatches,
          "raw test predictors on a view supersede the gathered fold");

    ext_rng_destroy(rng);
  }

  // invalid designations refuse cleanly at the store factory
  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rng == NULL || ext_rng_setSeed(rng, 91) != 0) {
      check(false, "linear view refusals: rng creation");
      return;
    }
    std::vector<size_t> rows(n);
    for (size_t i = 0; i < n; ++i) rows[i] = i;

    // built without the designation, the view holds no raw values to serve
    ColumnStore bare;
    bare.buildFromParent(parent, rows.data(), n, nullptr, 0);
    check(createSamplerOverStore(std::move(bare), y.data(), nullptr, nullptr,
                                 ResponseFamily::gaussian, 1.0, 3.0,
                                 0.37804942330213542, options, &rng) ==
            nullptr,
          "view without gathered covariates refused");

    ColumnStore store;
    store.buildFromParent(parent, rows.data(), n, nullptr, 0, covariates, 1);
    SamplerOptions bad = options;
    size_t outOfRange[] = {7};
    bad.leafCovariateColumns = outOfRange;
    check(createSamplerOverStore(std::move(store), y.data(), nullptr, nullptr,
                                 ResponseFamily::gaussian, 1.0, 3.0,
                                 0.37804942330213542, bad, &rng) == nullptr,
          "out-of-range designation refused over a store");
    ext_rng_destroy(rng);
  }

  printf("ok: linear leaf views\n");
}

static void testGroupedMath(ext_rng* rng) {
  // built-in tau priors and the tau posterior against R constants
  // (scratchpad grouped_reference.R): dcauchy/dgamma at scale 0.55, then
  // the R loop's posterior closure at J = 5, b.sq = 0.8
  checkNear(logTauPrior(TauPriorKind::cauchy, 0.55, 0.05),
            -0.55512338423029517, 1.0e-12, "cauchy tau prior at 0.05");
  checkNear(logTauPrior(TauPriorKind::cauchy, 0.55, 0.30),
            -0.80734814484534678, 1.0e-12, "cauchy tau prior at 0.30");
  checkNear(logTauPrior(TauPriorKind::cauchy, 0.55, 1.10),
            -2.1563307975278803, 1.0e-12, "cauchy tau prior at 1.10");
  checkNear(logTauPrior(TauPriorKind::gamma, 0.55, 0.05),
            -3.3745978698239458, 1.0e-12, "gamma tau prior at 0.05");
  checkNear(logTauPrior(TauPriorKind::gamma, 0.55, 0.30),
            -1.1415041205273175, 1.0e-12, "gamma tau prior at 0.30");
  checkNear(logTauPrior(TauPriorKind::gamma, 0.55, 1.10),
            -0.64712509887738068, 1.0e-12, "gamma tau prior at 1.10");
  checkNear(logTauPosterior(0.20, 5.0, 0.8, TauPriorKind::cauchy, 0.55),
            -2.6238937031546601, 1.0e-12, "tau posterior at 0.20");
  checkNear(logTauPosterior(0.45, 5.0, 0.8, TauPriorKind::cauchy, 0.55),
            0.9578598022153062, 1.0e-12, "tau posterior at 0.45");
  checkNear(logTauPosterior(0.90, 5.0, 0.8, TauPriorKind::cauchy, 0.55),
            -1.8162012038679745, 1.0e-12, "tau posterior at 0.90");
  check(logTauPosterior(0.0, 5.0, 0.8, TauPriorKind::cauchy, 0.55) ==
        -HUGE_VAL, "tau posterior refuses the boundary");

  // conjugate b update: empirical moments of drawGroupEffects against the
  // per-group posterior computed in R, weighted and unweighted
  const double z[] = {0.61, -0.32, 0.85, 0.10, -0.44, 0.92,
                      0.33, -0.15, 0.58, -0.71, 0.05, 0.40};
  const double w[] = {1.0, 0.5, 2.0, 1.5, 1.0, 0.25,
                      3.0, 1.0, 0.75, 1.25, 2.5, 1.0};
  const double fits[] = {0.20, -0.10, 0.40, 0.05, -0.30, 0.55,
                         0.25, -0.05, 0.35, -0.50, 0.10, 0.15};
  const std::uint32_t groups[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
  const double sigma = 0.7, tau = 0.4;
  const double weightedMean[] = {0.15813953488372096, 0.011127819548872186,
                                 0.0040875912408759162};
  const double weightedSd[] = {0.24652625377117468, 0.24279079146675359,
                               0.23922014416069304};
  const double unweightedMean[] = {0.097699115044247775,
                                   0.029734513274336287,
                                   0.031150442477876111};
  const double unweightedSd = 0.26340184314740722;

  const size_t numDraws = 100000;
  double effects[3], scratch[3];
  double sum[3], sumSq[3];
  for (int pass = 0; pass < 2; ++pass) {
    const double* weights = pass == 0 ? w : nullptr;
    for (size_t j = 0; j < 3; ++j) sum[j] = sumSq[j] = 0.0;
    for (size_t r = 0; r < numDraws; ++r) {
      drawGroupEffects(rng, z, weights, fits, groups, 12, 3, sigma, tau,
                       effects, scratch);
      for (size_t j = 0; j < 3; ++j) {
        sum[j] += effects[j];
        sumSq[j] += effects[j] * effects[j];
      }
    }
    for (size_t j = 0; j < 3; ++j) {
      double mean = sum[j] / static_cast<double>(numDraws);
      double sd = std::sqrt(sumSq[j] / static_cast<double>(numDraws) -
                            mean * mean);
      double expectedMean = pass == 0 ? weightedMean[j] : unweightedMean[j];
      double expectedSd = pass == 0 ? weightedSd[j] : unweightedSd;
      checkNear(mean, expectedMean, 4.0 * expectedSd / std::sqrt(100000.0),
                pass == 0 ? "weighted b posterior mean"
                          : "unweighted b posterior mean");
      checkNear(sd, expectedSd, 0.02 * expectedSd,
                pass == 0 ? "weighted b posterior sd"
                          : "unweighted b posterior sd");
    }
  }

  // the slice sampler recovers the tau posterior's moments (quadrature
  // reference from R) under both built-in priors
  const double sliceMean[] = {0.47739017758490643, 0.57907089567288195};
  const double sliceSd[] = {0.1836338997856671, 0.24745727859830713};
  for (int pass = 0; pass < 2; ++pass) {
    TauPriorKind kind = pass == 0 ? TauPriorKind::cauchy
                                  : TauPriorKind::gamma;
    auto logDensity = [kind](double t) {
      return logTauPosterior(t, 5.0, 0.8, kind, 0.55);
    };
    double x = 0.45, total = 0.0, totalSq = 0.0;
    for (size_t r = 0; r < numDraws; ++r) {
      x = sliceSampleOnce(rng, logDensity, x, 0.55, 0.0, HUGE_VAL);
      total += x;
      totalSq += x * x;
    }
    double mean = total / static_cast<double>(numDraws);
    double sd =
      std::sqrt(totalSq / static_cast<double>(numDraws) - mean * mean);
    checkNear(mean, sliceMean[pass], 0.01, "slice-sampled tau mean");
    checkNear(sd, sliceSd[pass], 0.01, "slice-sampled tau sd");
  }

  printf("ok: grouped math\n");
}

static void testGroupedEndToEnd(ext_rng* rng) {
  // gaussian grouped recovery: y = f(x) + b_g + eps; the sampler must
  // recover b (per-group posterior means), tau, sigma, and record fits
  // that exclude the random intercepts
  const size_t n = 800, numGroups = 20, numSamples = 300;
  std::vector<double> x(n * 2), y(n), f(n), b(numGroups);
  std::vector<std::uint32_t> groups(n);
  const double tauTrue = 0.8, sigmaTrue = 0.5;
  for (size_t j = 0; j < numGroups; ++j)
    b[j] = tauTrue * ext_rng_simulateStandardNormal(rng);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    groups[i] = static_cast<std::uint32_t>(i % numGroups);
    f[i] = x[i] < 0.5 ? -1.0 : (x[i + n] < 0.5 ? 0.5 : 2.0);
    y[i] = f[i] + b[groups[i]] +
           sigmaTrue * ext_rng_simulateStandardNormal(rng);
  }
  double yMean = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= static_cast<double>(n);
  double ySd = 0.0;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / static_cast<double>(n - 1));

  SamplerOptions options;
  options.numTrees = 50;
  options.groupIndices = groups.data();
  options.numGroups = numGroups;
  options.tauPriorKind = TauPriorKind::cauchy;
  options.tauPriorScale = ySd;
  options.tauSliceSteps = 5;

  ext_rng* chainRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                     NULL);
  ext_rng_setSeed(chainRng, 2026);
  ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                         ResponseFamily::gaussian, ySd, 3.0,
                         0.37804942330213542, options, &chainRng);

  std::vector<double> sigmaSamples(numSamples), tauSamples(numSamples);
  std::vector<double> effectSamples(numGroups * numSamples);
  std::vector<double> trainSamples(n * numSamples);
  Results results;
  results.sigma = sigmaSamples.data();
  results.tau = tauSamples.data();
  results.groupEffects = effectSamples.data();
  results.trainingFits = trainSamples.data();
  sampler.run(200, numSamples, results);

  double tauMean = 0.0, sigmaMean = 0.0;
  for (size_t s = 0; s < numSamples; ++s) {
    tauMean += tauSamples[s];
    sigmaMean += sigmaSamples[s];
  }
  tauMean /= static_cast<double>(numSamples);
  sigmaMean /= static_cast<double>(numSamples);
  checkNear(tauMean, tauTrue, 0.35, "grouped tau recovery");
  checkNear(sigmaMean, sigmaTrue, 0.1, "grouped sigma recovery");

  std::vector<double> bHat(numGroups, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t j = 0; j < numGroups; ++j)
      bHat[j] += effectSamples[j + s * numGroups];
  for (size_t j = 0; j < numGroups; ++j)
    bHat[j] /= static_cast<double>(numSamples);
  double bCor;
  {
    double mx = 0.0, my = 0.0;
    for (size_t j = 0; j < numGroups; ++j) { mx += b[j]; my += bHat[j]; }
    mx /= numGroups; my /= numGroups;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (size_t j = 0; j < numGroups; ++j) {
      sxy += (b[j] - mx) * (bHat[j] - my);
      sxx += (b[j] - mx) * (b[j] - mx);
      syy += (bHat[j] - my) * (bHat[j] - my);
    }
    bCor = sxy / std::sqrt(sxx * syy);
  }
  check(bCor > 0.9, "grouped b recovery");

  // recorded fits are f(x)-only: adding the ranef posterior means back must
  // shrink the residuals
  std::vector<double> fitMean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i) fitMean[i] += trainSamples[i + s * n];
  double sseWithout = 0.0, sseWith = 0.0;
  for (size_t i = 0; i < n; ++i) {
    fitMean[i] /= static_cast<double>(numSamples);
    double residual = y[i] - fitMean[i];
    sseWithout += residual * residual;
    residual -= bHat[groups[i]];
    sseWith += residual * residual;
  }
  check(sseWith < sseWithout, "recorded fits exclude the random intercepts");

  ext_rng_destroy(chainRng);
  printf("ok: grouped end to end (tau %.2f, sigma %.2f, b cor %.3f)\n",
         tauMean, sigmaMean, bCor);
}

static void testGroupedBinary(ext_rng* rng) {
  // probit and logistic compose with the decorator: the working weights
  // (unit / Polya-Gamma) enter the conjugate b update, rel.scale is the
  // binary 0.5, and the group intercepts are recovered on the latent scale
  const size_t n = 1200, numGroups = 8, numSamples = 250;
  std::vector<double> x(n * 2), b(numGroups);
  const double tauTrue = 0.9;
  for (size_t j = 0; j < numGroups; ++j)
    b[j] = tauTrue * ext_rng_simulateStandardNormal(rng);

  for (int pass = 0; pass < 2; ++pass) {
    ResponseFamily family = pass == 0 ? ResponseFamily::probit
                                      : ResponseFamily::logistic;
    std::vector<double> y(n);
    std::vector<std::uint32_t> groups(n);
    for (double& v : x) v = runif01();
    for (size_t i = 0; i < n; ++i) {
      groups[i] = static_cast<std::uint32_t>(i % numGroups);
      double eta = (x[i] < 0.5 ? -0.75 : 0.75) + b[groups[i]];
      if (pass == 0) {
        y[i] = eta + ext_rng_simulateStandardNormal(rng) > 0.0 ? 1.0 : 0.0;
      } else {
        y[i] = runif01() < 1.0 / (1.0 + std::exp(-eta)) ? 1.0 : 0.0;
      }
    }

    SamplerOptions options;
    options.numTrees = 25;
    // the R-side per-family node scales: probit 3, logistic pi sqrt(3)
    options.nodeScale = pass == 0 ? 3.0
                                  : 3.14159265358979324 * std::sqrt(3.0);
    options.groupIndices = groups.data();
    options.numGroups = numGroups;
    options.tauPriorKind = TauPriorKind::cauchy;
    options.tauPriorScale = 0.5;
    options.tauSliceSteps = 5;

    ext_rng* chainRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                       NULL);
    ext_rng_setSeed(chainRng, 3000 + pass);
    ClassicSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                           family, 1.0, 3.0, 1.0, options, &chainRng);

    std::vector<double> tauSamples(numSamples);
    std::vector<double> effectSamples(numGroups * numSamples);
    Results results;
    results.tau = tauSamples.data();
    results.groupEffects = effectSamples.data();
    sampler.run(150, numSamples, results);

    std::vector<double> bHat(numGroups, 0.0);
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t j = 0; j < numGroups; ++j)
        bHat[j] += effectSamples[j + s * numGroups];
    double mx = 0.0, my = 0.0;
    for (size_t j = 0; j < numGroups; ++j) {
      bHat[j] /= static_cast<double>(numSamples);
      mx += b[j]; my += bHat[j];
    }
    mx /= numGroups; my /= numGroups;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (size_t j = 0; j < numGroups; ++j) {
      sxy += (b[j] - mx) * (bHat[j] - my);
      sxx += (b[j] - mx) * (b[j] - mx);
      syy += (bHat[j] - my) * (bHat[j] - my);
    }
    double bCor = sxy / std::sqrt(sxx * syy);
    check(bCor > 0.85, pass == 0 ? "grouped probit b recovery"
                                 : "grouped logistic b recovery");
    for (double t : tauSamples)
      check(t > 0.0 && std::isfinite(t), "grouped binary tau draws valid");

    ext_rng_destroy(chainRng);
  }
  printf("ok: grouped binary families\n");
}

static void testGroupedStateRoundTrip() {
  // grouped state (b, tau) rides the chain state: a restored sampler
  // continues bitwise identically, and mismatched effect vectors are
  // refused
  const size_t n = 200, numGroups = 6, numChains = 2, numSamples = 5;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<std::uint32_t> groups(n);
  for (size_t i = 0; i < n; ++i)
    groups[i] = static_cast<std::uint32_t>(i % numGroups);

  SamplerOptions options;
  options.numTrees = 15;
  options.numChains = numChains;
  options.groupIndices = groups.data();
  options.numGroups = numGroups;
  options.tauPriorScale = 1.0;
  options.tauSliceSteps = 3;

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  };

  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 4242);
  ClassicSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(state.chains[0].groupEffects.size() == numGroups,
        "grouped state carries the effects");
  check(state.chains[0].groupTau > 0.0, "grouped state carries tau");

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 1111);
  ClassicSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state), "a grouped state restores");

  std::vector<double> sigmaA(numSamples * numChains);
  std::vector<double> tauA(numSamples * numChains);
  std::vector<double> trainA(n * numSamples * numChains);
  Results resultsA;
  resultsA.sigma = sigmaA.data();
  resultsA.tau = tauA.data();
  resultsA.trainingFits = trainA.data();
  original.run(0, numSamples, resultsA);

  std::vector<double> sigmaB(numSamples * numChains);
  std::vector<double> tauB(numSamples * numChains);
  std::vector<double> trainB(n * numSamples * numChains);
  Results resultsB;
  resultsB.sigma = sigmaB.data();
  resultsB.tau = tauB.data();
  resultsB.trainingFits = trainB.data();
  restored.run(0, numSamples, resultsB);

  check(sigmaA == sigmaB, "restored grouped chains draw identical sigmas");
  check(tauA == tauB, "restored grouped chains draw identical taus");
  check(trainA == trainB, "restored grouped chains draw identical fits");

  // an effects vector of the wrong length must be refused before anything
  // is overwritten
  SamplerStateData badState(state);
  badState.chains[0].groupEffects.resize(numGroups - 1);
  check(!restored.setState(badState),
        "a mismatched grouped state is refused");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: grouped state round trip\n");
}

// A logical matrix held both densely and as CSC arrays, for comparing the
// two build paths over identical values.
struct CscFixture {
  size_t n = 0, p = 0;
  std::vector<double> dense;   // column-major, zeros where nothing stored
  std::vector<int> pointers;   // p + 1
  std::vector<int> rows;
  std::vector<double> values;

  // fraction of rows stored per column; stored NaNs count as entries, the
  // Matrix convention for missing values
  void build(size_t n_, const std::vector<double>& nonzeroFractions,
             size_t numMissingPerColumn = 0) {
    n = n_;
    p = nonzeroFractions.size();
    dense.assign(n * p, 0.0);
    pointers.assign(p + 1, 0);
    rows.clear();
    values.clear();
    for (size_t j = 0; j < p; ++j) {
      size_t numMissing = 0;
      for (size_t i = 0; i < n; ++i) {
        if (runif01() >= nonzeroFractions[j]) continue;
        double value = 0.5 + runif01();
        if (numMissing < numMissingPerColumn) {
          value = std::nan("");
          ++numMissing;
        }
        dense[i + j * n] = value;
        rows.push_back(static_cast<int>(i));
        values.push_back(value);
      }
      pointers[j + 1] = static_cast<int>(rows.size());
    }
  }
};

static void testSparseKernel() {
  const size_t n = 1000;
  const xint_t zeroCode = 5;
  std::vector<xint_t> denseCodes(n, zeroCode);
  size_t numWords = (n + 63) / 64;
  std::vector<std::uint64_t> bits(numWords, 0);
  std::vector<std::uint32_t> wordRanks(numWords, 0);
  std::vector<xint_t> nzCodes;
  for (size_t i = 0; i < n; ++i) {
    if (runif01() >= 0.1) continue;
    xint_t code = static_cast<xint_t>(1.0 + runif01() * 249.0);
    denseCodes[i] = code;
    bits[i >> 6] |= std::uint64_t{1} << (i & 63u);
    nzCodes.push_back(code);
  }
  std::uint32_t runningRank = 0;
  for (size_t w = 0; w < numWords; ++w) {
    wordRanks[w] = runningRank;
    runningRank += static_cast<std::uint32_t>(std::popcount(bits[w]));
  }

  // a scrambled index segment, partitioned at several cuts against a dense
  // reference
  std::vector<size_t> segment(n);
  for (size_t i = 0; i < n; ++i) segment[i] = i;
  for (size_t i = n - 1; i > 0; --i) {
    size_t k = static_cast<size_t>(runif01() * static_cast<double>(i + 1));
    std::swap(segment[i], segment[k]);
  }

  const xint_t cuts[] = { 0, 4, 5, 100, 250 };
  bool allMatch = true;
  for (xint_t cut : cuts) {
    std::vector<size_t> indices(segment);
    size_t numOnLeft = misc_partitionIndicesSparse(
      bits.data(), wordRanks.data(), nzCodes.data(), zeroCode, cut,
      indices.data(), indices.size());

    size_t expected = 0;
    for (size_t i : segment) expected += denseCodes[i] <= cut ? 1 : 0;
    allMatch &= numOnLeft == expected;
    for (size_t k = 0; k < indices.size(); ++k)
      allMatch &= (denseCodes[indices[k]] <= cut) == (k < numOnLeft);
    std::vector<size_t> sorted(indices);
    std::sort(sorted.begin(), sorted.end());
    std::vector<size_t> sortedSegment(segment);
    std::sort(sortedSegment.begin(), sortedSegment.end());
    allMatch &= sorted == sortedSegment;
  }
  check(allMatch, "sparse kernel matches dense reference at every cut");
  check(misc_partitionIndicesSparse(bits.data(), wordRanks.data(),
                                    nzCodes.data(), zeroCode, 100, nullptr,
                                    0) == 0,
        "sparse kernel handles the empty segment");
  printf("ok: sparse partition kernel\n");
}

static void testSparseColumnStore() {
  const size_t n = 300;
  CscFixture fixture;
  // rank, densified, fully dense, rank-with-missing
  fixture.build(n, {0.05, 0.5, 1.0, 0.1}, 0);
  // mark ~5 of column 3's stored entries missing
  for (int k = fixture.pointers[3];
       k < fixture.pointers[3] + 5 && k < fixture.pointers[4]; ++k) {
    fixture.values[static_cast<size_t>(k)] = std::nan("");
    fixture.dense[static_cast<size_t>(fixture.rows[static_cast<size_t>(k)]) +
                  3 * n] = std::nan("");
  }

  for (bool useQuantiles : { false, true }) {
    ColumnStore fromCsc;
    fromCsc.buildFromCsc(fixture.pointers.data(), fixture.rows.data(),
                         fixture.values.data(), n, fixture.p, nullptr, 100,
                         useQuantiles);
    ColumnStore fromDense;
    fromDense.build(fixture.dense.data(), n, fixture.p, 100, useQuantiles);

    check(fromCsc.builtFromCsc && fromCsc.hasSparse,
          "CSC store flags itself");
    check(fromCsc.columnIsSparse(0) && !fromCsc.columnIsSparse(1) &&
          !fromCsc.columnIsSparse(2) && fromCsc.columnIsSparse(3),
          "density threshold splits the storage tiers");
    check(fromCsc.hasMissing[0] == 0 && fromCsc.hasMissing[1] == 0 &&
          fromCsc.hasMissing[2] == 0 && fromCsc.hasMissing[3] == 1,
          "stored NaN entries mark the column missing");

    bool cutsMatch = true;
    for (size_t j = 0; j < fixture.p; ++j) {
      cutsMatch &= fromCsc.numCuts[j] == fromDense.numCuts[j];
      cutsMatch &= fromCsc.cutPoints[j] == fromDense.cutPoints[j];
    }
    check(cutsMatch, useQuantiles
          ? "CSC quantile cuts bitwise-match the dense builder's"
          : "CSC uniform cuts bitwise-match the dense builder's");

    bool codesMatch = true;
    for (size_t j = 0; j < fixture.p; ++j)
      for (size_t i = 0; i < n; ++i)
        codesMatch &= fromCsc.codeAt(j, i) == fromDense.codes[i + j * n];
    check(codesMatch, "CSC codes match the dense builder's at every cell");

    check(fromCsc.sparseColumn(0).zeroCode == fromDense.codeFor(0, 0.0),
          "rank zero code is the quantized zero");

    // a view over the CSC parent densifies and gathers identical codes
    std::vector<size_t> viewRows;
    for (size_t i = 0; i < n; i += 3) viewRows.push_back(i);
    ColumnStore view;
    view.buildFromParent(fromCsc, viewRows.data(), viewRows.size(), nullptr,
                         0);
    check(!view.hasSparse && !view.builtFromCsc,
          "views over CSC parents densify");
    bool viewMatches = true;
    for (size_t j = 0; j < fixture.p; ++j)
      for (size_t i = 0; i < viewRows.size(); ++i)
        viewMatches &= view.codes[i + j * viewRows.size()] ==
          fromCsc.codeAt(j, viewRows[i]);
    check(viewMatches, "view gathers parent codes through the rank layout");
  }
  printf("ok: sparse column store\n");
}

static void testSparseEndToEnd() {
  // densified tier: a CSC build whose columns all densify runs the dense
  // kernels over identical codes, so draws are bitwise identical to the
  // dense-matrix path
  {
    const size_t n = 300, p = 3;
    CscFixture fixture;
    fixture.build(n, {1.0, 0.6, 0.9});
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i)
      y[i] = std::sin(3.0 * fixture.dense[i]) + fixture.dense[i + n] +
             0.5 * runif01();

    ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 1234) != 0 ||
        ext_rng_setSeed(rngB, 1234) != 0) {
      check(false, "sparse end-to-end: rng creation");
      return;
    }

    SamplerOptions options;
    options.numTrees = 25;
    ClassicSampler dense(fixture.dense.data(), y.data(), n, p, nullptr,
                         nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rngA);

    SamplerOptions cscOptions(options);
    cscOptions.cscColumnPointers = fixture.pointers.data();
    cscOptions.cscRowIndices = fixture.rows.data();
    cscOptions.cscValues = fixture.values.data();
    ClassicSampler sparse(nullptr, y.data(), n, p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, cscOptions, &rngB);
    check(!sparse.data().hasSparse && sparse.data().builtFromCsc,
          "all-densified CSC sampler holds no rank columns");

    const size_t numBurnIn = 30, numSamples = 40;
    std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
    std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
    Results resultsA, resultsB;
    resultsA.sigma = sigmaA.data();
    resultsA.trainingFits = fitsA.data();
    resultsB.sigma = sigmaB.data();
    resultsB.trainingFits = fitsB.data();
    dense.run(numBurnIn, numSamples, resultsA);
    sparse.run(numBurnIn, numSamples, resultsB);

    check(sigmaA == sigmaB && fitsA == fitsB,
          "all-densified CSC sampler bitwise-matches the dense path");
    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
  }

  // rank tier: signal riding sparse columns is recovered, missing entries
  // included (exercising the sparse MIA partition)
  {
    const size_t n = 600;
    CscFixture fixture;
    fixture.build(n, {0.1, 0.1, 0.1, 0.1});
    // NAs on a noise column so rules there route through the MIA sibling
    for (int k = fixture.pointers[2];
         k < fixture.pointers[2] + 8 && k < fixture.pointers[3]; ++k)
      fixture.values[static_cast<size_t>(k)] = std::nan("");

    std::vector<double> f(n), y(n);
    double yMean = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double x0 = fixture.dense[i];
      double x1 = fixture.dense[i + n];
      f[i] = 2.0 * x0 - 1.5 * x1;
      double u1 = runif01(), u2 = runif01();
      double normal =
        std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
      y[i] = f[i] + 0.3 * normal;
      yMean += y[i] / static_cast<double>(n);
    }
    double ySd = 0.0;
    for (size_t i = 0; i < n; ++i)
      ySd += (y[i] - yMean) * (y[i] - yMean);
    ySd = std::sqrt(ySd / static_cast<double>(n - 1));

    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rng == NULL || ext_rng_setSeed(rng, 2026) != 0) {
      check(false, "sparse end-to-end: rng creation");
      return;
    }
    SamplerOptions options;
    options.numTrees = 75;
    options.cscColumnPointers = fixture.pointers.data();
    options.cscRowIndices = fixture.rows.data();
    options.cscValues = fixture.values.data();
    ClassicSampler sampler(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                           ResponseFamily::gaussian, ySd, 3.0,
                           0.37804942330213542, options, &rng);
    check(sampler.data().hasSparse, "rank-tier sampler holds rank columns");

    const size_t numBurnIn = 200, numSamples = 300;
    std::vector<double> sigmaDraws(numSamples);
    std::vector<double> trainingFits(n * numSamples);
    Results results;
    results.sigma = sigmaDraws.data();
    results.trainingFits = trainingFits.data();
    sampler.run(numBurnIn, numSamples, results);

    std::vector<double> posteriorMean(n, 0.0);
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t i = 0; i < n; ++i)
        posteriorMean[i] +=
          trainingFits[i + s * n] / static_cast<double>(numSamples);
    double sseFit = 0.0, sseMean = 0.0;
    for (size_t i = 0; i < n; ++i) {
      sseFit += (posteriorMean[i] - f[i]) * (posteriorMean[i] - f[i]);
      sseMean += (yMean - f[i]) * (yMean - f[i]);
    }
    check(sseFit < 0.2 * sseMean,
          "sparse end-to-end: fit explains most signal");
    double sigmaPosteriorMean = 0.0;
    for (double s : sigmaDraws)
      sigmaPosteriorMean += s / static_cast<double>(numSamples);
    check(sigmaPosteriorMean > 0.2 && sigmaPosteriorMean < 0.45,
          "sparse end-to-end: sigma near truth (0.3)");
    ext_rng_destroy(rng);
    printf("ok: sparse end-to-end (sigma posterior mean %.3f, sse ratio "
           "%.3f)\n",
           sigmaPosteriorMean, sseFit / sseMean);
  }
}

static void testSparseStateRoundTrip() {
  const size_t n = 240, numChains = 2, numSamples = 5;
  CscFixture fixture;
  fixture.build(n, {0.1, 0.5, 0.08});
  // missing entries on a rank column ride the state round trip too
  for (int k = fixture.pointers[2];
       k < fixture.pointers[2] + 4 && k < fixture.pointers[3]; ++k)
    fixture.values[static_cast<size_t>(k)] = std::nan("");
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = 1.5 * fixture.dense[i] + 0.4 * runif01();

  SamplerOptions options;
  options.numTrees = 15;
  options.numChains = numChains;
  options.cscColumnPointers = fixture.pointers.data();
  options.cscRowIndices = fixture.rows.data();
  options.cscValues = fixture.values.data();

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  };

  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 4242);
  ClassicSampler original(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 1111);
  ClassicSampler restored(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state), "a sparse-store state restores");

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

  check(sigmaA == sigmaB, "restored sparse chains draw identical sigmas");
  check(trainA == trainB, "restored sparse chains draw identical fits");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: sparse state round trip\n");
}

// A mixed dense/CSC predictor set in the overall layout
// [dense0, csc0, dense1, csc1, csc2], dense1 optionally 4-category
// categorical; full is the equivalent all-dense matrix for reference builds.
struct MixedFixture {
  size_t n = 0;
  static const size_t p = 5;
  CscFixture csc;
  std::vector<double> denseSource;  // column-major, n x 2
  std::vector<double> full;         // column-major, n x 5
  std::vector<std::int32_t> sources;
  std::vector<ColumnType> types;

  void build(size_t n_, const std::vector<double>& cscFractions,
             bool categoricalDense1, size_t numMissingCsc2 = 0) {
    n = n_;
    csc.build(n, cscFractions);
    for (int k = csc.pointers[2];
         k < csc.pointers[2] + static_cast<int>(numMissingCsc2) &&
         k < csc.pointers[3]; ++k) {
      csc.values[static_cast<size_t>(k)] = std::nan("");
      csc.dense[static_cast<size_t>(csc.rows[static_cast<size_t>(k)]) +
                2 * n] = std::nan("");
    }
    denseSource.resize(n * 2);
    for (size_t i = 0; i < n; ++i) {
      denseSource[i] = 2.0 * runif01() - 1.0;
      denseSource[i + n] = categoricalDense1
        ? std::floor(runif01() * 4.0) : runif01();
    }
    sources = { 0, ~0, 1, ~1, ~2 };
    types.assign(p, ColumnType::ordinal);
    if (categoricalDense1) types[2] = ColumnType::categorical;
    full.resize(n * p);
    const double* columns[p] = {
      denseSource.data(), csc.dense.data(), denseSource.data() + n,
      csc.dense.data() + n, csc.dense.data() + 2 * n
    };
    for (size_t j = 0; j < p; ++j)
      std::memcpy(full.data() + j * n, columns[j], n * sizeof(double));
  }

  void applyOptions(SamplerOptions& options) {
    options.cscColumnPointers = csc.pointers.data();
    options.cscRowIndices = csc.rows.data();
    options.cscValues = csc.values.data();
    options.mixedDenseValues = denseSource.data();
    options.columnSources = sources.data();
    options.columnTypes = types.data();
  }
};

static void testMixedColumnStore() {
  const size_t n = 300;
  MixedFixture fixture;
  // csc0 rank, csc1 densified, csc2 rank-with-missing; dense1 categorical
  fixture.build(n, {0.05, 0.5, 0.1}, true, 5);

  for (bool useQuantiles : { false, true }) {
    ColumnStore mixed;
    mixed.buildMixed(fixture.denseSource.data(), fixture.csc.pointers.data(),
                     fixture.csc.rows.data(), fixture.csc.values.data(),
                     fixture.sources.data(), n, fixture.p, nullptr, 100,
                     useQuantiles, fixture.types.data());
    ColumnStore reference;
    reference.build(fixture.full.data(), n, fixture.p, 100, useQuantiles,
                    fixture.types.data());

    check(mixed.builtFromCsc && mixed.hasSparse, "mixed store flags itself");
    check(!mixed.columnIsCscBacked(0) && mixed.columnIsCscBacked(1) &&
          !mixed.columnIsCscBacked(2) && mixed.columnIsCscBacked(3) &&
          mixed.columnIsCscBacked(4),
          "the source map splits the column backings");
    check(!mixed.columnIsSparse(0) && mixed.columnIsSparse(1) &&
          !mixed.columnIsSparse(2) && !mixed.columnIsSparse(3) &&
          mixed.columnIsSparse(4),
          "the density threshold tiers the CSC-backed columns");
    check(mixed.rawColumn(0) == fixture.denseSource.data() &&
          mixed.rawColumn(2) == fixture.denseSource.data() + n &&
          mixed.rawColumn(1) == nullptr && mixed.rawColumn(3) == nullptr,
          "raw values are served exactly for dense-backed columns");
    check(mixed.hasMissing[0] == 0 && mixed.hasMissing[1] == 0 &&
          mixed.hasMissing[2] == 0 && mixed.hasMissing[3] == 0 &&
          mixed.hasMissing[4] == 1,
          "stored NaN entries mark only their column missing");

    bool cutsMatch = true;
    for (size_t j = 0; j < fixture.p; ++j) {
      cutsMatch &= mixed.numCuts[j] == reference.numCuts[j];
      cutsMatch &= mixed.cutPoints[j] == reference.cutPoints[j];
    }
    check(cutsMatch, useQuantiles
          ? "mixed quantile cuts bitwise-match the dense builder's"
          : "mixed uniform cuts bitwise-match the dense builder's");
    check(mixed.numCuts[2] == 4, "the dense-backed categorical keeps its levels");

    bool codesMatch = true;
    for (size_t j = 0; j < fixture.p; ++j)
      for (size_t i = 0; i < n; ++i)
        codesMatch &= mixed.codeAt(j, i) == reference.codes[i + j * n];
    check(codesMatch, "mixed codes match the dense builder's at every cell");

    check(mixed.sparseColumn(1).zeroCode == reference.codeFor(1, 0.0),
          "rank zero code is the quantized zero");

    // a view over the mixed parent densifies its codes and gathers raw
    // values only for columns the parent serves raw
    std::vector<size_t> viewRows;
    for (size_t i = 0; i < n; i += 3) viewRows.push_back(i);
    size_t rawColumnsToGather[] = { 0, 1 };
    ColumnStore view;
    view.buildFromParent(mixed, viewRows.data(), viewRows.size(), nullptr, 0,
                         rawColumnsToGather, 2);
    check(!view.hasSparse && !view.builtFromCsc,
          "views over mixed parents densify");
    bool viewMatches = true;
    for (size_t j = 0; j < fixture.p; ++j)
      for (size_t i = 0; i < viewRows.size(); ++i)
        viewMatches &= view.codes[i + j * viewRows.size()] ==
          mixed.codeAt(j, viewRows[i]);
    check(viewMatches, "view gathers parent codes through the mixed layout");
    const double* gathered = view.rawColumn(0);
    bool gatherMatches = gathered != nullptr;
    if (gathered != nullptr)
      for (size_t i = 0; i < viewRows.size(); ++i)
        gatherMatches &= gathered[i] == fixture.denseSource[viewRows[i]];
    check(gatherMatches, "view gathers raw values of a dense-backed column");
    check(view.rawColumn(1) == nullptr,
          "view leaves a CSC-backed column ungathered");
    double mean = 0.0, sd = 0.0;
    double referenceMean = 0.0, referenceSd = 0.0;
    standardizationMomentsForColumn(fixture.denseSource.data(), n,
                                    &referenceMean, &referenceSd);
    check(view.suppliedStandardization(0, &mean, &sd) &&
          mean == referenceMean && sd == referenceSd,
          "view standardization constants come from the parent's full data");
  }
  printf("ok: mixed column store\n");
}

static void testMixedEndToEnd() {
  // all CSC sources densify, so the mixed sampler runs the dense kernels
  // over identical codes: draws bitwise-match the dense-matrix path, the
  // dense-backed categorical column included
  const size_t n = 300;
  MixedFixture fixture;
  fixture.build(n, {0.6, 0.9, 1.0}, true);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = std::sin(3.0 * fixture.denseSource[i]) + fixture.csc.dense[i] +
           0.3 * (fixture.denseSource[i + n] == 2.0 ? 1.0 : 0.0) +
           0.5 * runif01();

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 8471) != 0 ||
      ext_rng_setSeed(rngB, 8471) != 0) {
    check(false, "mixed end-to-end: rng creation");
    return;
  }

  SamplerOptions options;
  options.numTrees = 25;
  options.columnTypes = fixture.types.data();
  ClassicSampler dense(fixture.full.data(), y.data(), n, fixture.p, nullptr,
                       nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                       0.37804942330213542, options, &rngA);

  SamplerOptions mixedOptions;
  mixedOptions.numTrees = 25;
  fixture.applyOptions(mixedOptions);
  ClassicSampler mixed(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                       ResponseFamily::gaussian, 1.0, 3.0,
                       0.37804942330213542, mixedOptions, &rngB);
  check(!mixed.data().hasSparse && mixed.data().builtFromCsc,
        "all-densified mixed sampler holds no rank columns");

  const size_t numBurnIn = 30, numSamples = 40;
  std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
  std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = fitsA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = fitsB.data();
  dense.run(numBurnIn, numSamples, resultsA);
  mixed.run(numBurnIn, numSamples, resultsB);

  check(sigmaA == sigmaB && fitsA == fitsB,
        "mixed sampler bitwise-matches the dense path");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: mixed end-to-end\n");
}

static void testMixedLinearLeaves() {
  // linear leaves designate dense-backed mixed columns (they serve raw
  // values); designating a CSC-backed column refuses at the factory
  const size_t n = 300;
  MixedFixture fixture;
  fixture.build(n, {0.6, 0.9, 1.0}, false);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = (fixture.denseSource[i] > 0.0 ? fixture.denseSource[i + n] : 0.0) +
           0.4 * runif01();

  SamplerOptions mixedOptions;
  mixedOptions.numTrees = 25;
  fixture.applyOptions(mixedOptions);
  size_t covariates[] = { 2 };
  mixedOptions.leafCovariateColumns = covariates;
  mixedOptions.numLeafCovariates = 1;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 555) != 0 ||
      ext_rng_setSeed(rngB, 555) != 0) {
    check(false, "mixed linear leaves: rng creation");
    return;
  }

  {
    SamplerOptions bad = mixedOptions;
    size_t cscBacked[] = { 3 };
    bad.leafCovariateColumns = cscBacked;
    check(createSampler(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                        ResponseFamily::gaussian, 1.0, 3.0,
                        0.37804942330213542, bad, &rngA) == nullptr,
          "CSC-backed leaf covariate refused");
  }

  std::unique_ptr<SamplerBase> mixed = createSampler(
    nullptr, y.data(), n, fixture.p, nullptr, nullptr,
    ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, mixedOptions,
    &rngB);
  check(mixed != nullptr, "dense-backed leaf covariate creates");
  if (mixed == nullptr) {
    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
    return;
  }

  SamplerOptions denseOptions;
  denseOptions.numTrees = 25;
  denseOptions.leafCovariateColumns = covariates;
  denseOptions.numLeafCovariates = 1;
  std::unique_ptr<SamplerBase> dense = createSampler(
    fixture.full.data(), y.data(), n, fixture.p, nullptr, nullptr,
    ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, denseOptions,
    &rngA);
  check(dense != nullptr, "dense reference linear-leaf sampler creates");
  if (dense == nullptr) {
    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
    return;
  }

  const size_t numBurnIn = 30, numSamples = 40;
  std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
  std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = fitsA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = fitsB.data();
  dense->run(numBurnIn, numSamples, resultsA);
  mixed->run(numBurnIn, numSamples, resultsB);
  check(sigmaA == sigmaB && fitsA == fitsB,
        "mixed linear-leaf sampler bitwise-matches the dense path");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: mixed linear leaves\n");
}

static void testMixedStateRoundTrip() {
  const size_t n = 240, numChains = 2, numSamples = 5;
  MixedFixture fixture;
  // both CSC tiers, missing entries, and a dense-backed categorical ride
  // the state round trip
  fixture.build(n, {0.1, 0.5, 0.08}, true, 4);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = 1.5 * fixture.denseSource[i] + fixture.csc.dense[i] +
           0.4 * runif01();

  SamplerOptions options;
  options.numTrees = 15;
  options.numChains = numChains;
  fixture.applyOptions(options);

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], seed + static_cast<std::uint32_t>(c));
    }
  };

  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 6161);
  ClassicSampler original(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 2323);
  ClassicSampler restored(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state), "a mixed-store state restores");

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

  check(sigmaA == sigmaB, "restored mixed chains draw identical sigmas");
  check(trainA == trainB, "restored mixed chains draw identical fits");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: mixed state round trip\n");
}

struct GPLeafFixture {
  static constexpr size_t n = 8;
  std::vector<double> x, z, w;
  ColumnStore store;
  std::vector<size_t> indexBuffer;
  Tree tree;
  double k = 2.0, scale = 0.5 / std::sqrt(10.0), sigmaSq = 0.04;
  double theta[2] = {0.7, 1.3};

  GPLeafFixture() : indexBuffer(n) {
    double x1[] = {0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9};
    double u1[] = {1.3, -0.4, 0.7, 2.2, -1.5, 0.6, 0.1, -0.8};
    double u2[] = {0.2, 1.1, -0.6, 0.4, 0.9, -1.3, 0.5, -0.2};
    x.assign(3 * n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      x[i] = x1[i];
      x[i + n] = u1[i];
      x[i + 2 * n] = u2[i];
    }
    z = {0.35, -0.10, 0.22, 0.55, -0.42, 0.16, 0.03, -0.25};
    w = {1.0, 2.0, 0.5, 1.0, 1.5, 1.0, 0.8, 1.2};
    store.build(x.data(), n, 3, 100);
    tree.initialize(indexBuffer.data(), n);
  }
};

static void testGPLeafMarginal() {
  GPLeafFixture f;

  GPGaussianLeaf leaf;
  leaf.scale = f.scale;
  size_t columns[] = {1};
  leaf.initialize(f.store, columns, 1, f.theta, 256);
  check(leaf.numCovariates() == 1, "gp leaf covariate count");

  // root marginals against the R reference (gp_leaf_reference.R)
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, 0),
            -8.4818528980077694, 1e-9, "gp marginal, weighted root");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), nullptr,
                                                f.k, f.sigmaSq, 0),
            -7.7744742015613255, 1e-9, "gp marginal, unit-weight root");

  // median pairwise-distance lengthscale heuristic matches R's median
  GPGaussianLeaf heuristic;
  heuristic.scale = f.scale;
  heuristic.initialize(f.store, columns, 1, nullptr, 256);
  checkNear(heuristic.lengthscales()[0], 1.0967438055849543, 1e-12,
            "median lengthscale, first covariate");

  // two covariates exercise the multi-column kernel
  GPGaussianLeaf leaf2;
  leaf2.scale = f.scale;
  size_t columns2[] = {1, 2};
  leaf2.initialize(f.store, columns2, 2, f.theta, 256);
  checkNear(leaf2.logIntegratedLikelihoodForNode(f.tree, f.z.data(),
                                                 f.w.data(), f.k, f.sigmaSq, 0),
            -8.7417667358884774, 1e-9, "gp marginal, two covariates");
  checkNear(leaf2.covariateMeans()[0], 0.275, 1e-12, "gp covariate mean");
  checkNear(leaf2.covariateSds()[0], 1.185326959112970, 1e-12,
            "gp covariate sd");
  GPGaussianLeaf heuristic2;
  heuristic2.scale = f.scale;
  heuristic2.initialize(f.store, columns2, 2, nullptr, 256);
  checkNear(heuristic2.lengthscales()[1], 0.94224419672838522, 1e-12,
            "median lengthscale, second covariate");

  // children of a split at x1 ~ 0.5 (cut 50 of the uniform grid)
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(50);
  f.tree.birth(f.store, 0, rule, f.z.data(), f.w.data());
  int32_t leftChild = f.tree.at(0).leftChild;
  check(f.tree.at(leftChild).numObservations() == 4 &&
        f.tree.at(leftChild + 1).numObservations() == 4,
        "gp fixture split partitions 4/4");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild),
            -5.0284693098158764, 1e-9, "gp marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild + 1),
            -3.7990612354161635, 1e-9, "gp marginal, right child");

  // a constant kernel (huge lengthscale) approaches the constant leaf's
  // formula by Sherman-Morrison
  Tree rootTree;
  std::vector<size_t> rootIndices(GPLeafFixture::n);
  rootTree.initialize(rootIndices.data(), GPLeafFixture::n);
  rootTree.setNodeAverages(f.z.data(), f.w.data());
  ConstantGaussianLeaf constant{f.scale};
  double constantScore = constant.logIntegratedLikelihoodForNode(
    rootTree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0);
  GPGaussianLeaf flatKernel;
  flatKernel.scale = f.scale;
  double bigTheta[] = {1.0e6};
  flatKernel.initialize(f.store, columns, 1, bigTheta, 256);
  checkNear(flatKernel.logIntegratedLikelihoodForNode(
              rootTree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0),
            constantScore, 1e-4, "constant-kernel limit matches constant leaf");

  // over the size cap the fallback IS the constant leaf's math, exactly;
  // nodes at or under the cap keep the gp marginal
  GPGaussianLeaf capped;
  capped.scale = f.scale;
  capped.initialize(f.store, columns, 1, f.theta, 4);
  check(capped.logIntegratedLikelihoodForNode(rootTree, f.z.data(),
                                              f.w.data(), f.k, f.sigmaSq, 0) ==
          constantScore,
        "over-cap marginal equals the constant leaf's exactly");
  checkNear(capped.logIntegratedLikelihoodForNode(
              f.tree, f.z.data(), f.w.data(), f.k, f.sigmaSq, leftChild),
            -5.0284693098158764, 1e-9, "at-cap marginal stays gp");

  printf("ok: gp leaf marginal\n");
}

// Zero-weight members have infinite noise variance: they contribute no
// likelihood and their draws follow the conditional law given the
// positive-weight rows. The marginal must equal a leaf holding exactly
// the positive rows, bitwise.
static void testGPLeafZeroWeights(ext_rng* rng) {
  GPLeafFixture f;

  GPGaussianLeaf leaf;
  leaf.scale = f.scale;
  size_t columns[] = {1};
  leaf.initialize(f.store, columns, 1, f.theta, 256);

  std::vector<double> wZero(f.w);
  wZero[2] = 0.0;
  wZero[5] = 0.0;
  double zeroScore = leaf.logIntegratedLikelihoodForNode(
    f.tree, f.z.data(), wZero.data(), f.k, f.sigmaSq, 0);
  check(std::isfinite(zeroScore), "zero-weight gp marginal is finite");

  // a tree whose root holds exactly the positive rows, in the same order
  Tree positiveTree;
  std::vector<size_t> positiveIndices(6);
  positiveTree.initialize(positiveIndices.data(), 6);
  size_t positives[] = {0, 1, 3, 4, 6, 7};
  for (size_t i = 0; i < 6; ++i) positiveIndices[i] = positives[i];
  double positiveScore = leaf.logIntegratedLikelihoodForNode(
    positiveTree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0);
  check(zeroScore == positiveScore,
        "zero-weight marginal equals the positive-subset leaf bitwise");

  // a zero-weight row duplicating a positive row's designated value must
  // draw (nearly) that row's function value every time - the conditional
  // mean at a duplicated input is the input's own value up to the nugget
  const size_t nd = 5;
  double xd1[] = {0.1, 0.3, 0.5, 0.7, 0.9};
  double ud1[] = {1.1, -0.7, 0.4, 1.1, -1.6};  // row 3 duplicates row 0
  std::vector<double> xd(nd * 2);
  for (size_t i = 0; i < nd; ++i) {
    xd[i] = xd1[i];
    xd[i + nd] = ud1[i];
  }
  std::vector<double> zd = {0.4, -0.2, 0.1, 0.0, 0.3};
  std::vector<double> wd = {1.0, 1.0, 2.0, 0.0, 1.5};
  ColumnStore dupStore;
  dupStore.build(xd.data(), nd, 2, 100);
  std::vector<size_t> dupIndices(nd);
  Tree dupTree;
  dupTree.initialize(dupIndices.data(), nd);
  GPGaussianLeaf dupLeaf;
  dupLeaf.scale = f.scale;
  double dupTheta[] = {0.9};
  dupLeaf.initialize(dupStore, columns, 1, dupTheta, 256);

  std::vector<double> fits(nd, 0.0);
  const size_t numDraws = 4000;
  double maxGap = 0.0, minZeroRow = HUGE_VAL, maxZeroRow = -HUGE_VAL;
  bool allFinite = true;
  for (size_t s = 0; s < numDraws; ++s) {
    dupLeaf.beginTreeDraw(dupTree);
    dupLeaf.drawFromPosteriorForNode(rng, dupTree, zd.data(), wd.data(),
                                     f.k, f.sigmaSq, 0, fits.data());
    for (size_t i = 0; i < nd; ++i) allFinite &= std::isfinite(fits[i]);
    double gap = std::fabs(fits[3] - fits[0]);
    if (gap > maxGap) maxGap = gap;
    if (fits[3] < minZeroRow) minZeroRow = fits[3];
    if (fits[3] > maxZeroRow) maxZeroRow = fits[3];
  }
  check(allFinite, "zero-weight gp draws are finite");
  check(maxGap < 1.0e-2,
        "zero-weight duplicate row tracks its positive twin");
  check(maxZeroRow > minZeroRow, "zero-weight row draws vary");

  // no positive members at all: no likelihood - scores like an empty
  // leaf, draws from the prior
  std::vector<double> wAllZero(nd, 0.0);
  check(dupLeaf.logIntegratedLikelihoodForNode(dupTree, zd.data(),
                                               wAllZero.data(), f.k,
                                               f.sigmaSq, 0) == 0.0,
        "all-zero-weight gp marginal scores like an empty leaf");
  dupLeaf.beginTreeDraw(dupTree);
  dupLeaf.drawFromPosteriorForNode(rng, dupTree, zd.data(),
                                   wAllZero.data(), f.k, f.sigmaSq, 0,
                                   fits.data());
  bool priorFinite = true;
  for (size_t i = 0; i < nd; ++i) priorFinite &= std::isfinite(fits[i]);
  check(priorFinite, "all-zero-weight gp draw is finite");

  printf("ok: gp leaf zero weights\n");
}

static void testGPLeafDraw(ext_rng* rng) {
  GPLeafFixture f;
  GPGaussianLeaf leaf;
  leaf.scale = f.scale;
  size_t columns[] = {1};
  leaf.initialize(f.store, columns, 1, f.theta, 256);
  leaf.beginTreeDraw(f.tree);

  // posterior moments of f at the weighted root, from the R reference
  const double expectedMean[8] = {
    0.097801590251712, -0.051930303515435, 0.061351828901789,
    0.097580147953186, -0.100816341185500, 0.052842762505061,
    0.002453754801322, -0.086391000822176};
  const double expectedVar[8] = {
    0.004714667725125, 0.003866729033325, 0.004351767565350,
    0.005197096398532, 0.004634971897899, 0.004301914991551,
    0.004053876726524, 0.003974642704896};

  const int numDraws = 50000;
  std::vector<double> fits(f.n), sum(f.n, 0.0), sumSq(f.n, 0.0);
  double statsSum = 0.0;
  FunctionLeafDrawStats stats = leaf.drawFromPosteriorForNode(
    rng, f.tree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0, fits.data());
  check(stats.numParams == 8.0, "gp draw reports one parameter per member");
  for (int r = 0; r < numDraws; ++r) {
    stats = leaf.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(),
                                          f.k, f.sigmaSq, 0, fits.data());
    statsSum += stats.sumSquaredParams;
    for (size_t i = 0; i < f.n; ++i) {
      sum[i] += fits[i];
      sumSq[i] += fits[i] * fits[i];
    }
  }
  bool meansMatch = true, variancesMatch = true;
  for (size_t i = 0; i < f.n; ++i) {
    double mean = sum[i] / numDraws;
    double variance = sumSq[i] / numDraws - mean * mean;
    meansMatch = meansMatch && std::fabs(mean - expectedMean[i]) < 2e-3;
    variancesMatch =
      variancesMatch && std::fabs(variance - expectedVar[i]) < 5e-4;
  }
  check(meansMatch, "gp posterior draw means");
  check(variancesMatch, "gp posterior draw variances");
  checkNear(statsSum / numDraws, 0.067270287625328065, 5e-3,
            "chi-k contribution matches E[f' C^-1 f]");

  // prior draws are centered with marginal variance (scale / k)^2
  double priorSum = 0.0, priorSumSq = 0.0;
  const int numPriorDraws = 20000;
  leaf.beginTreeDraw(f.tree);
  for (int r = 0; r < numPriorDraws; ++r) {
    leaf.drawFromPriorForNode(rng, f.tree, f.k, 0, fits.data());
    priorSum += fits[0];
    priorSumSq += fits[0] * fits[0];
  }
  double priorMean = priorSum / numPriorDraws;
  double s = f.scale / f.k;
  checkNear(priorMean, 0.0, 2e-3, "gp prior draw mean");
  checkNear(priorSumSq / numPriorDraws - priorMean * priorMean, s * s, 5e-4,
            "gp prior draw variance");

  // a test row duplicating a training row reproduces its drawn value up to
  // the nugget
  const size_t numTest = 5;
  std::vector<double> xTest(numTest * 3);
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < numTest; ++i)
      xTest[i + j * numTest] = f.x[i + j * f.n];
  f.store.buildTest(xTest.data(), numTest);
  leaf.rebuildTestCovariates(f.store);
  leaf.beginTreeDraw(f.tree);
  leaf.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(), f.k,
                                f.sigmaSq, 0, fits.data());
  bool predictionsMatch = true;
  for (size_t i = 0; i < numTest; ++i)
    predictionsMatch = predictionsMatch &&
      std::fabs(leaf.fitForTestObservationForNode(f.tree, 0, i) - fits[i]) <
        1e-3;
  check(predictionsMatch, "duplicated test rows reproduce drawn values");

  // over-cap draws broadcast one constant value and serve it to test rows
  GPGaussianLeaf capped;
  capped.scale = f.scale;
  capped.initialize(f.store, columns, 1, f.theta, 4);
  capped.beginTreeDraw(f.tree);
  f.tree.setNodeAverages(f.z.data(), f.w.data());
  stats = capped.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(),
                                          f.k, f.sigmaSq, 0, fits.data());
  check(stats.numParams == 1.0, "over-cap draw reports one parameter");
  bool allEqual = true;
  for (size_t i = 1; i < f.n; ++i) allEqual = allEqual && fits[i] == fits[0];
  check(allEqual, "over-cap draw broadcasts a constant");
  check(capped.fitForTestObservationForNode(f.tree, 0, 0) == fits[0],
        "over-cap test fit serves the constant");

  // an empty leaf writes no fits, consumes no draws, and predicts zero
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(50);
  f.tree.birth(f.store, 0, rule, f.z.data(), f.w.data());
  int32_t rightChild = f.tree.at(0).leftChild + 1;
  f.tree.at(rightChild).begin = f.tree.at(rightChild).end;
  leaf.beginTreeDraw(f.tree);
  fits.assign(f.n, 123.0);
  leaf.drawFromPosteriorForNode(rng, f.tree, f.z.data(), f.w.data(), f.k,
                                f.sigmaSq, rightChild, fits.data());
  bool untouched = true;
  for (size_t i = 0; i < f.n; ++i) untouched = untouched && fits[i] == 123.0;
  check(untouched, "empty leaf writes no fits");
  check(leaf.fitForTestObservationForNode(f.tree, rightChild, 0) == 0.0,
        "empty leaf predicts zero");

  printf("ok: gp leaf draw\n");
}

static void testGPLeafEndToEnd(ext_rng* rng) {
  const size_t n = 250, p = 1;
  const double pi = 3.141592653589793;
  std::vector<double> x(n), fTrue(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    fTrue[i] = std::sin(4.0 * pi * x[i]);
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * pi * u2);
    y[i] = fTrue[i] + 0.2 * normal;
  }

  double yMean = 0.0, ySd = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= (double) n;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / (double) (n - 1));

  SamplerOptions options;
  options.numTrees = 10;
  size_t covariates[] = {0};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;
  options.gpMaxLeafSize = 100;

  Sampler<GPGaussianLeaf> gpSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);

  const size_t numTest = 5;
  std::vector<double> xTest(x.begin(), x.begin() + numTest);
  gpSampler.setTestPredictors(xTest.data(), numTest);

  const size_t numBurnIn = 150, numSamples = 200;
  std::vector<double> sigmaDraws(numSamples);
  std::vector<double> trainingFits(n * numSamples);
  std::vector<double> testFits(numTest * numSamples);
  Results results;
  results.sigma = sigmaDraws.data();
  results.trainingFits = trainingFits.data();
  results.testFits = testFits.data();
  gpSampler.run(numBurnIn, numSamples, results);

  std::vector<double> posteriorMean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i)
      posteriorMean[i] += trainingFits[i + s * n] / (double) numSamples;

  double sseGp = 0.0, sseMean = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sseGp += (posteriorMean[i] - fTrue[i]) * (posteriorMean[i] - fTrue[i]);
    sseMean += (yMean - fTrue[i]) * (yMean - fTrue[i]);
  }
  check(sseGp < 0.15 * sseMean, "gp end to end explains the smooth signal");

  double sigmaPosteriorMean = 0.0;
  for (double sd : sigmaDraws) sigmaPosteriorMean += sd / (double) numSamples;
  check(sigmaPosteriorMean > 0.15 && sigmaPosteriorMean < 0.3,
        "gp end to end: sigma near truth (0.2)");

  // recorded test fits are the leaves' conditional means at the duplicated
  // rows: equal to the training fits up to the nugget
  double maxTestGap = 0.0;
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < numTest; ++i) {
      double gap = std::fabs(testFits[i + s * numTest] -
                             trainingFits[i + s * n]);
      if (gap > maxTestGap) maxTestGap = gap;
    }
  check(maxTestGap < 0.05, "gp test fits track duplicated training rows");

  // the same run under the constant leaf: the gp fit is closer to the
  // smooth truth
  SamplerOptions constantOptions;
  constantOptions.numTrees = 10;
  Sampler<ConstantGaussianLeaf> constantSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, constantOptions, &rng);
  std::vector<double> constantFits(n * numSamples);
  Results constantResults;
  constantResults.trainingFits = constantFits.data();
  constantSampler.run(numBurnIn, numSamples, constantResults);
  std::vector<double> constantPosteriorMean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i)
      constantPosteriorMean[i] += constantFits[i + s * n] / (double) numSamples;
  double sseConstant = 0.0;
  for (size_t i = 0; i < n; ++i)
    sseConstant += (constantPosteriorMean[i] - fTrue[i]) *
                   (constantPosteriorMean[i] - fTrue[i]);
  check(sseGp < sseConstant, "gp fit beats the constant leaf on smooth truth");

  // sampled k stays sane over the f' C^-1 f accumulation
  SamplerOptions kOptions = options;
  kOptions.updateK = true;
  Sampler<GPGaussianLeaf> kSampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, kOptions, &rng);
  const size_t numKSamples = 50;
  std::vector<double> kDraws(numKSamples), kSigma(numKSamples);
  Results kResults;
  kResults.k = kDraws.data();
  kResults.sigma = kSigma.data();
  kSampler.run(50, numKSamples, kResults);
  bool kSane = true;
  for (double kDraw : kDraws)
    kSane = kSane && std::isfinite(kDraw) && kDraw > 0.05 && kDraw < 50.0;
  check(kSane, "sampled k stays finite and positive under gp leaves");

  // the format surfaces work post-stage 2: a captured state restores into
  // its own sampler and flatten emits records
  SamplerStateData state;
  gpSampler.getState(state);
  check(gpSampler.setState(state), "a gp state restores into its own sampler");
  std::vector<FlatNode> flatNodes;
  std::vector<std::uint32_t> flatCounts;
  gpSampler.flattenTree(0, 0, flatNodes, flatCounts);
  check(!flatNodes.empty() && flatCounts.size() == flatNodes.size(),
        "gp flatten emits records and counts");

  // prior draws compose: fresh structures and prior parameters, then a
  // short continued run stays finite
  gpSampler.sampleTreesFromPrior();
  gpSampler.sampleNodeParametersFromPrior();
  const size_t numContinued = 5;
  std::vector<double> continuedFits(n * numContinued);
  Results continuedResults;
  continuedResults.trainingFits = continuedFits.data();
  gpSampler.run(0, numContinued, continuedResults);
  bool continuedFinite = true;
  for (double fit : continuedFits)
    continuedFinite = continuedFinite && std::isfinite(fit);
  check(continuedFinite, "prior draws continue into finite sampling");

  // the type-erased facade instantiates every member for the function-leaf
  // shape; a fresh sampler's live prediction is the response center (all
  // trees are over-cap roots with zero fits)
  SamplerFacade<GPGaussianLeaf> facade(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);
  std::vector<double> predictions(numTest, 1.0);
  facade.predict(xTest.data(), numTest, predictions.data());
  bool predictionsCentered = true;
  for (double prediction : predictions)
    predictionsCentered = predictionsCentered &&
      std::isfinite(prediction) && prediction == predictions[0];
  check(predictionsCentered, "fresh gp live prediction sits at the center");

  printf("ok: gp leaf end to end (sse ratio %.3f, constant ratio %.3f, "
         "sigma %.3f, test gap %.4f)\n",
         sseGp / sseMean, sseConstant / sseMean, sigmaPosteriorMean,
         maxTestGap);
}

static void testGPLeafFormats(ext_rng* rng) {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = std::sin(3.0 * x[i]) + 0.5 * std::cos(3.0 * x[i + n]) +
           0.2 * normal;
  }
  double yMean = 0.0, ySd = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= (double) n;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / (double) (n - 1));

  SamplerOptions options;
  options.numTrees = 5;
  size_t covariates[] = {0, 1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 2;
  options.gpMaxLeafSize = 60;
  options.keepTrees = true;
  options.numSamplesToStore = 20;

  Sampler<GPGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);

  const size_t numTest = 8;
  std::vector<double> xTest(numTest * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < numTest; ++i)
      xTest[i + j * numTest] = runif01();
  sampler.setTestPredictors(xTest.data(), numTest);

  const size_t numSamples = 20;
  std::vector<double> sigmaDraws(numSamples);
  std::vector<double> testFits(numTest * numSamples);
  Results results;
  results.sigma = sigmaDraws.data();
  results.testFits = testFits.data();
  sampler.run(60, numSamples, results);

  // saved replay reproduces the recorded test fits bitwise: the stored
  // alpha weights and covariate rows are the exact values the live
  // evaluation used, and the replay repeats its arithmetic order
  std::vector<double> replayed(numTest * numSamples, 0.0);
  sampler.predict(xTest.data(), numTest, replayed.data());
  bool replayMatches = true;
  for (size_t s = 0; s < numSamples && replayMatches; ++s)
    for (size_t i = 0; i < numTest && replayMatches; ++i)
      replayMatches = replayed[i + s * numTest] == testFits[i + s * numTest];
  check(replayMatches, "gp saved replay bit-matches recorded test fits");

  // bitwise state round trip into a fresh sampler, then identical
  // continuation of both
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 1234);
  Sampler<GPGaussianLeaf> restored(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng2);
  restored.setTestPredictors(xTest.data(), numTest);

  SamplerStateData state;
  sampler.getState(state);
  check(restored.setState(state), "gp state installs into a fresh sampler");

  // malformed side channels refuse without mutating anything
  SamplerStateData malformed = state;
  malformed.chains[0].savedTreeParams[0].pop_back();
  check(!restored.setState(malformed), "truncated gp side channel refused");
  malformed = state;
  malformed.chains[0].treeParams[0].pop_back();
  check(!restored.setState(malformed), "short gp fits slab refused");

  const size_t numContinued = 10;
  std::vector<double> sigmaA(numContinued), sigmaB(numContinued);
  std::vector<double> trainA(n * numContinued), trainB(n * numContinued);
  std::vector<double> testA(numTest * numContinued),
                      testB(numTest * numContinued);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = trainA.data();
  resultsA.testFits = testA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = trainB.data();
  resultsB.testFits = testB.data();
  sampler.run(0, numContinued, resultsA);
  restored.run(0, numContinued, resultsB);
  check(sigmaA == sigmaB, "restored gp chain draws identical sigmas");
  check(trainA == trainB, "restored gp chain draws identical fits");
  check(testA == testB, "restored gp chain draws identical test fits");

  // flattened live trees emit walkable side-channel blocks
  std::vector<FlatNode> nodes;
  std::vector<std::uint32_t> counts;
  std::vector<double> blocks;
  sampler.flattenTree(0, 0, nodes, counts, &blocks);
  std::vector<size_t> blockOffsets;
  check(!nodes.empty() && counts.size() == nodes.size() &&
          computeFunctionBlockOffsets(blocks.data(), blocks.size(),
                                      (nodes.size() + 1) / 2, 2,
                                      blockOffsets),
        "gp flatten emits records and walkable blocks");

  ext_rng_destroy(rng2);
  printf("ok: gp leaf formats\n");
}

// The cross-sweep kernel cache must be invisible. A cold clone is built
// over the ORIGINAL data (sharing the cut grid), restored from the warm
// sampler's pre-mutation state, and given the identical mutation, so its
// empty cache recomputes what the warm sampler serves from cache. Cycle
// one perturbs the designated column WITHIN its quantization bins: members
// stay identical, so only the regather-clears-the-cache path keeps stale
// kernels from hitting. Cycle two moves the non-designated column across
// bins: observations re-route and the member-list comparison must miss.
static void testGPLeafKernelCache(ext_rng* rng) {
  const size_t n = 150, p = 2;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    x[i + n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = std::sin(3.0 * x[i]) + 0.3 * x[i + n] + 0.2 * normal;
  }
  double yMean = 0.0, ySd = 0.0;
  for (size_t i = 0; i < n; ++i) yMean += y[i];
  yMean /= (double) n;
  for (size_t i = 0; i < n; ++i) ySd += (y[i] - yMean) * (y[i] - yMean);
  ySd = std::sqrt(ySd / (double) (n - 1));

  SamplerOptions options;
  options.numTrees = 5;
  size_t covariates[] = {0};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;
  options.gpMaxLeafSize = 60;

  Sampler<GPGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);

  std::vector<double> sigmaWarm(5);
  Results warm;
  warm.sigma = sigmaWarm.data();
  sampler.run(40, 5, warm);

  const size_t numContinued = 5;
  auto mutateAndCompare = [&](const std::vector<double>& xMutated,
                              const char* acceptLabel,
                              const char* matchLabel) {
    SamplerStateData state;
    sampler.getState(state);

    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngB, 99);
    Sampler<GPGaussianLeaf> cold(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      ySd, 3.0, 0.37804942330213542, options, &rngB);
    bool installed = cold.setState(state);

    bool accepted =
      sampler.setPredictor(xMutated.data(), true, false) ==
        PredictorUpdateResult::accepted &&
      cold.setPredictor(xMutated.data(), true, false) ==
        PredictorUpdateResult::accepted;
    check(accepted, acceptLabel);

    std::vector<double> sigmaA(numContinued), sigmaB(numContinued);
    std::vector<double> trainA(n * numContinued), trainB(n * numContinued);
    Results resultsA, resultsB;
    resultsA.sigma = sigmaA.data();
    resultsA.trainingFits = trainA.data();
    resultsB.sigma = sigmaB.data();
    resultsB.trainingFits = trainB.data();
    sampler.run(0, numContinued, resultsA);
    cold.run(0, numContinued, resultsB);
    check(installed && sigmaA == sigmaB && trainA == trainB, matchLabel);
    ext_rng_destroy(rngB);
  };

  // within-bin designated perturbation: partitions and member lists stay
  // put while the kernels' inputs move
  std::vector<double> xCurrent(x);
  for (size_t i = 0; i < n; ++i)
    xCurrent[i] += (static_cast<double>(i % 3) - 1.0) * 1.0e-9;
  mutateAndCompare(xCurrent, "gp designated-column mutation accepted",
                   "warm cache matches cold clone after designated mutation");

  // cross-bin non-designated change: re-quantization re-routes members
  for (size_t i = 0; i < n; ++i)
    xCurrent[i + n] = 0.05 + 0.9 * xCurrent[i + n];
  mutateAndCompare(xCurrent, "gp non-designated-column mutation accepted",
                   "warm cache matches cold clone after member re-route");

  printf("ok: gp leaf kernel cache\n");
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
  testConstantLeafSuffstatEquivalence();
  testTreeMechanics();
  testTreePriorMath();
  testEndToEndGaussian(rng);
  testRunCancellation(rng);
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
  testWeightedLogistic(rng);
  testCategoricalMechanics();
  testCategoricalPriorMath(rng);
  testEndToEndCategorical(rng);
  testCategoricalMutation(rng);
  testWideCategorical(rng);
  testPooledMaskMechanics(rng);
  testPooledMaskSampler(rng);
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
  testMissingIngestion();
  testMissingMechanics();
  testMissingEndToEnd();
  testLinearLeafMarginal();
  testLinearLeafDraw(rng);
  testLinearLeafEndToEnd(rng);
  testLinearLeafFormats(rng);
  testLinearLeafMutation(rng);
  testLinearLeafViews();
  testGroupedMath(rng);
  testGroupedEndToEnd(rng);
  testGroupedBinary(rng);
  testGroupedStateRoundTrip();
  testSparseKernel();
  testSparseColumnStore();
  testSparseEndToEnd();
  testSparseStateRoundTrip();
  testMixedColumnStore();
  testMixedEndToEnd();
  testMixedLinearLeaves();
  testMixedStateRoundTrip();
  testGPLeafMarginal();
  testGPLeafDraw(rng);
  testGPLeafEndToEnd(rng);
  testGPLeafFormats(rng);
  testGPLeafKernelCache(rng);
  testGPLeafZeroWeights(rng);

  ext_rng_destroy(rng);

  if (failures > 0) {
    printf("%d failure(s)\n", failures);
    return 1;
  }
  printf("all tests passed\n");
  return 0;
}
