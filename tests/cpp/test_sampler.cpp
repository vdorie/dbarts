#include "common.hpp"

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

  ConstantLeafSampler full(x.data(), y.data(), n, p, nullptr, nullptr,
                      ResponseFamily::gaussian, 1.0, 3.0,
                      0.37804942330213542, options, &rngA);

  ColumnStore parent;
  parent.build(x.data(), n, p, options.maxNumCuts);
  std::vector<size_t> rows(n);
  for (size_t i = 0; i < n; ++i) rows[i] = i;
  ColumnStore store;
  store.buildFromParent(parent, rows.data(), n, nullptr, 0);
  ConstantLeafSampler view(std::move(store), y.data(), nullptr, nullptr,
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
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
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
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
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
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
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::probit,
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
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
  ConstantLeafSampler a(x.data(), y.data(), n, p, nullptr, w.data(),
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngA);
  ConstantLeafSampler b(x.data(), y.data(), n, p, nullptr, w.data(),
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
  ConstantLeafSampler wsampler(x.data(), y.data(), n, p, nullptr, w.data(),
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
    ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
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

// Test-fit routing splits its rows across the thread budget above a size
// cutoff; routing draws no rng and each row owns its output slot, so the
// recorded test fits must be bitwise-identical at any thread count. nTest is
// well past the cutoff so numThreads = 8 engages every worker.
static void testTestFitThreadInvariance() {
  const size_t n = 300, nTest = 300000, numSamples = 8;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  std::uint_least32_t st = 20260707u;
  auto u01 = [&]() {
    st = st * 1664525u + 1013904223u;
    return static_cast<double>(st >> 8) * (1.0 / 16777216.0);
  };
  std::vector<double> xTest(nTest * 2);
  for (double& v : xTest) v = u01();

  auto runSampler = [&](size_t numThreads, std::vector<double>& testFits) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 2718);
    SamplerOptions options;
    options.numTrees = 25;
    options.numChains = 1;
    options.numThreads = numThreads;
    ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                           ResponseFamily::gaussian, 1.0, 3.0,
                           0.37804942330213542, options, &rng);
    sampler.setTestPredictors(xTest.data(), nTest);

    std::vector<double> sigma(numSamples);
    testFits.assign(nTest * numSamples, 0.0);
    Results results;
    results.sigma = sigma.data();
    results.testFits = testFits.data();
    sampler.run(50, numSamples, results);
    ext_rng_destroy(rng);
  };

  std::vector<double> serial, threaded2, threaded8;
  runSampler(1, serial);
  runSampler(2, threaded2);
  runSampler(8, threaded8);
  check(serial == threaded2, "2-thread test fits match serial bitwise");
  check(serial == threaded8, "8-thread test fits match serial bitwise");

  bool finite = true;
  for (double f : serial) finite &= std::isfinite(f);
  check(finite, "parallel test fits are finite");
  printf("ok: test-fit thread-count invariance\n");
}

static void testSetData(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

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

static void testSetResponseScaleLock(ext_rng* rng) {
  // setResponse defaults to locking the creation-time transform; updateScale
  // re-anchors it like setData does, holding sigma on the original scale
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  // a linear transform doubling the range: the re-anchored path must double
  // the transform width and halve the internal sigma, the locked path neither
  std::vector<double> yScaled(n);
  for (size_t i = 0; i < n; ++i) yScaled[i] = 2.0 * y[i] + 3.0;

  std::unique_ptr<ConstantLeafSampler> lockedPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& locked(*lockedPtr);
  SamplerStateData before;
  locked.getState(before);
  double widthBefore = before.chains[0].fitMax - before.chains[0].fitMin;
  double sigmaInternalBefore = before.chains[0].sigma;

  locked.setResponse(yScaled.data(), false);
  SamplerStateData afterLocked;
  locked.getState(afterLocked);
  checkNear(afterLocked.chains[0].fitMax - afterLocked.chains[0].fitMin,
            widthBefore, 1e-12, "locked setResponse keeps the transform");
  checkNear(afterLocked.chains[0].sigma, sigmaInternalBefore, 1e-12,
            "locked setResponse keeps internal sigma");

  std::unique_ptr<ConstantLeafSampler> scaledPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& scaled(*scaledPtr);
  SamplerStateData beforeScaled;
  scaled.getState(beforeScaled);
  double widthBeforeScaled =
    beforeScaled.chains[0].fitMax - beforeScaled.chains[0].fitMin;
  double sigmaOriginalBeforeScaled = scaled.sigma(0);
  scaled.setResponse(yScaled.data(), true);
  SamplerStateData afterScaled;
  scaled.getState(afterScaled);
  checkNear(afterScaled.chains[0].fitMax - afterScaled.chains[0].fitMin,
            2.0 * widthBeforeScaled, 1e-12 * widthBeforeScaled,
            "re-anchored setResponse doubles the transform");
  checkNear(scaled.sigma(0), sigmaOriginalBeforeScaled,
            1e-12 * sigmaOriginalBeforeScaled,
            "re-anchored setResponse preserves sigma on the original scale");

  printf("ok: setResponse scale lock\n");
}

static void testSetDataResize(ext_rng* rng) {
  const size_t n = 200, n2 = 320, nTest = 50, nTest2 = 30;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::unique_ptr<ConstantLeafSampler> samplerPtr = makeBurnedInSampler(x, y, n, rng);
  ConstantLeafSampler& sampler(*samplerPtr);

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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::probit,
                         1.0, 3.0, 0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(50, 0, empty);

  std::vector<double> x2, y2Continuous;
  makeMutationData(x2, y2Continuous, n2);
  std::vector<double> y2(n2);
  for (size_t i = 0; i < n2; ++i) y2[i] = y2Continuous[i] > 0.0 ? 1.0 : 0.0;

  sampler.setData(x2.data(), y2.data(), n2, nullptr, nullptr, nullptr, 0);
  check(sampler.numObservations() == n2, "probit setData resizes");

  // latents cold-initialize to 2 y - 1
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr, ResponseFamily::gaussian,
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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

static void testWideCategorical(ext_rng* rng) {
  // an inline categorical column whose masks reach high bit positions (past
  // 31, past the old 53-category double-encoding boundary); the mask rides
  // the flattened node directly
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
  // nonempty assignments; guards against draw schemes that pin low bits
  // for wide masks
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

  // the widest valid mask round-trips through the inline flat payload
  std::vector<FlatNode> flat;
  std::vector<double> paramByNode(tree.nodes.size(), 0.0);
  tree.flatten(store, paramByNode.data(), flat);
  check(flatKindOf(flat[0]) == FlatKind::categoricalInline &&
          flat[0].mask == rule.categoryDirections(),
        "wide mask flattens exactly");
  flat[0].mask = ((1ull << K) - 1ull) & ~1ull;
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
  flat[0].mask = 1ull << K;  // a bit past the K categories
  check(!flatTreeIsWellFormed(store, flat.data(), flat.size()),
        "a mask bit past the categories is rejected");

  // end to end: category codes carry no order, so recovering the two
  // groups requires subset splits over high codes
  SamplerOptions options;
  options.numTrees = 50;
  options.columnTypes = types;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 1, nullptr, nullptr,
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

static void testPooledMaskSampler(ext_rng* rng) {
  // end to end over one pooled column (K = 70, side-channel masks) and one
  // inline column (K = 60, mask in the flat node), with keepTrees, saved-tree
  // replay, live prediction, and a bitwise state round trip mid-run
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
    check(!probe.columnIsPooled(1) && probe.columnIsPooled(0),
          "60 categories stay inline, 70 pool");
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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

  // a stored state restores to an equivalent continuation
  SamplerStateData state;
  sampler.getState(state);
  check(!state.chains[0].forests[0].treeMasks.empty(),
        "wide-column states carry mask channels");
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 4242);
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rng2);
  restored.setTestPredictors(xTest.data(), nTest);
  check(restored.setState(state), "a pooled state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored pooled state reproduces the model");

  // live-tree prediction agrees with the final tree fits
  std::vector<double> livePredict(n);
  {
    ConstantLeafSampler live(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  ConstantLeafSampler samplerA(x.data(), y.data(), n, 2, w2.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngA);
  ext_rng* rngB = makeSeededRng();
  ConstantLeafSampler samplerB(x.data(), y.data(), n, 2, w1.data(), nullptr,
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
  ConstantLeafSampler samplerC(x.data(), y.data(), n, 2, w2.data(), nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rngC);
  samplerC.setTestPredictors(xTest.data(), nTest);
  samplerC.setTestOffset(testOffset.data());

  ext_rng* rngD = makeSeededRng();
  ConstantLeafSampler samplerD(x.data(), y.data(), n, 2, w2.data(), nullptr,
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
  ConstantLeafSampler samplerA(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 5.0, 0.9, optionsA,
                          &rngA);

  SamplerOptions optionsB;
  optionsB.numTrees = 25;
  ext_rng* rngB = makeSeededRng();
  ConstantLeafSampler samplerB(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  ConstantLeafSampler samplerC(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsC, &rngC);

  ext_rng* rngD = makeSeededRng();
  ConstantLeafSampler samplerD(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  ConstantLeafSampler samplerE(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, optionsE, &rngE);
  ext_rng* rngF = makeSeededRng();
  ConstantLeafSampler samplerF(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  ConstantLeafSampler samplerG(x.data(), y.data(), n, 2, nullptr, offset.data(),
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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

  // the state round trip carries the missing directions along
  SamplerStateData state;
  sampler.getState(state);
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 999);
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, &rng2);
  check(restored.setState(state), "a state with missing directions restores");

  checkStructuralRoundTrip(state, restored,
                           "restored state reproduces the missing directions");

  ext_rng_destroy(rng);
  ext_rng_destroy(rng2);
  printf("ok: missing end to end\n");
}

// BCF two-forest sampler: creation, a short run moving both forests, sane
// glue (finite, b0 != b1), setTreatment refresh. Statistical validation is
// benchmarks/R/bcf-exact.R; this is sanity only.
static void testBCFTwoForest(ext_rng* rng) {
  const size_t n = 400, p = 3;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  BCFSpec spec;
  spec.mu.numTrees = 50; spec.mu.base = 0.95; spec.mu.power = 2.0;
  spec.tau.numTrees = 25; spec.tau.base = 0.25; spec.tau.power = 3.0;
  spec.z = z.data();

  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);
  check(sampler.numForests() == 2, "BCF builds two forests");

  // A BCF sampler carries no test treatment vector, so a test blend is
  // ill-defined: the engine must record NaN test fits rather than the bare
  // prognostic forest (bcf-testfits-guard).
  const size_t nTest = 40;
  std::vector<double> xTest(nTest * p);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < nTest; ++i)
      xTest[i + j * nTest] = x[i + j * n];
  sampler.setTestPredictors(xTest.data(), nTest);

  const size_t numBurnIn = 200, numSamples = 200;
  std::vector<double> sigma(numSamples), fits(n * numSamples),
    testFits(nTest * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.trainingFits = fits.data();
  results.testFits = testFits.data();
  sampler.run(numBurnIn, numSamples, results);

  bool testFitsNaN = true;
  for (double f : testFits) testFitsNaN &= std::isnan(f);
  check(testFitsNaN, "BCF test fits come back NaN");

  std::vector<double> muFits(n), tauFits(n);
  sampler.forestTotalFits(0, 0, muFits.data());
  sampler.forestTotalFits(0, 1, tauFits.data());
  double muSS = 0.0, tauSS = 0.0;
  for (size_t i = 0; i < n; ++i) {
    muSS += muFits[i] * muFits[i];
    tauSS += tauFits[i] * tauFits[i];
  }
  check(muSS > 0.0 && tauSS > 0.0, "both BCF forests move off zero");

  double a, b0, b1;
  bool haveGlue = sampler.chain(0).bcfGlue(a, b0, b1);
  check(haveGlue && std::isfinite(a) && std::isfinite(b0) &&
          std::isfinite(b1),
        "BCF glue is finite");
  check(b0 != b1, "BCF treatment scales separate");

  bool sane = true;
  for (size_t i = 0; i < n * numSamples && sane; ++i)
    sane = std::isfinite(fits[i]);
  for (size_t s = 0; s < numSamples && sane; ++s) sane = sigma[s] > 0.0;
  check(sane, "BCF fits finite and sigma positive");

  // The final recorded training fits must be the fitScale * (a * mu + b_z * tau)
  // + shift blend, not the bare mu forest: recover the shift from one row and
  // confirm the linear structure over the rest against the current per-forest
  // fits and glue (the last sample reflects the post-run state).
  double fitScale = sampler.fitScale();
  const double* lastFit = fits.data() + (numSamples - 1) * n;
  auto blend = [&](size_t i) {
    return fitScale * (a * muFits[i] + (z[i] != 0.0 ? b1 : b0) * tauFits[i]);
  };
  double blendShift = lastFit[0] - blend(0);
  bool blendOk = true;
  for (size_t i = 0; i < n && blendOk; ++i)
    blendOk = std::fabs(lastFit[i] - (blend(i) + blendShift)) < 1.0e-8;
  check(blendOk, "BCF training fits are the a*mu + b_z*tau blend");

  std::vector<double> zControl(n, 0.0);
  sampler.setTreatment(zControl.data());
  std::vector<double> sigma2(10), fits2(n * 10);
  Results r2;
  r2.sigma = sigma2.data();
  r2.trainingFits = fits2.data();
  sampler.run(0, 10, r2);
  bool refreshed = true;
  for (size_t i = 0; i < n * 10 && refreshed; ++i)
    refreshed = std::isfinite(fits2[i]);
  check(refreshed, "BCF setTreatment refresh runs");

  printf("ok: BCF two-forest sampler\n");
}

// updateA/updateB false hold the glue fixed at (a, b0, b1) = (1, 0, 1)
// across a run, so the forests fit the a = 1, b_z = z model unchanged.
static void testBCFFixedGlue(ext_rng* rng) {
  const size_t n = 300, p = 2;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + n];
    y[i] = mu + z[i] * tau + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  BCFSpec spec;
  spec.mu.numTrees = 20; spec.mu.base = 0.95; spec.mu.power = 2.0;
  spec.tau.numTrees = 10; spec.tau.base = 0.25; spec.tau.power = 3.0;
  spec.updateA = false; spec.updateB = false;
  spec.z = z.data();

  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);

  Results results;
  sampler.run(50, 50, results);

  double a, b0, b1;
  sampler.chain(0).bcfGlue(a, b0, b1);
  check(a == 1.0 && b0 == 0.0 && b1 == 1.0, "BCF fixed glue holds");

  std::vector<double> tauFits(n);
  sampler.forestTotalFits(0, 1, tauFits.data());
  double tauSS = 0.0;
  for (size_t i = 0; i < n; ++i) tauSS += tauFits[i] * tauFits[i];
  check(tauSS > 0.0, "BCF treatment forest moves under fixed glue");

  printf("ok: BCF fixed glue\n");
}

void runSamplerTests(ext_rng* rng) {
  testBCFTwoForest(rng);
  testBCFFixedGlue(rng);
  testViewSamplerMatchesFull();
  testEndToEndGaussian(rng);
  testRunCancellation(rng);
  testEndToEndProbit(rng);
  testMultiChain();
  testTestFitThreadInvariance();
  testSetData(rng);
  testSetResponseScaleLock(rng);
  testSetDataResize(rng);
  testSetDataQuantileShrink(rng);
  testSetDataProbit(rng);
  testMultiChainSetData();
  testEndToEndLogistic(rng);
  testWeightedLogistic(rng);
  testEndToEndCategorical(rng);
  testWideCategorical(rng);
  testPooledMaskSampler(rng);
  testSetWeightsAndTestOffset();
  testSetControlAndModel();
  testMissingEndToEnd();
}
