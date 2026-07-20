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

// The opt-in fp32 running-residual instantiation (Sampler<..., float>,
// docs/design/reduced-precision-storage.md sec 3b): the fp32 treeY roll and the
// float-input suffstat kernels must run cleanly (this is the ASAN-covered path
// for that new reachable code) and, since the divergence from fp64 sits far
// below the MCMC noise floor, must pass the same recovery bar as the fp64 run.
static void testEndToEndGaussianFp32(ext_rng* rng) {
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
  options.fp32Residual = true;
  Sampler<ConstantGaussianLeaf, float> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    ySd, 3.0, 0.37804942330213542, options, &rng);

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
      posteriorMean[i] += trainingFits[i + s * n] / (double) numSamples;

  double sseFit = 0.0, sseMean = 0.0;
  bool allFinite = true;
  for (size_t i = 0; i < n; ++i) {
    allFinite = allFinite && std::isfinite(posteriorMean[i]);
    sseFit += (posteriorMean[i] - f[i]) * (posteriorMean[i] - f[i]);
    sseMean += (yMean - f[i]) * (yMean - f[i]);
  }
  check(allFinite, "fp32 end-to-end: fits are finite");
  check(sseFit < 0.2 * sseMean, "fp32 end-to-end: fit explains most signal");

  double sigmaPosteriorMean = 0.0;
  for (double s : sigmaDraws) sigmaPosteriorMean += s / (double) numSamples;
  check(sigmaPosteriorMean > 0.2 && sigmaPosteriorMean < 0.45,
        "fp32 end-to-end: sigma near truth (0.3)");

  printf("ok: end-to-end gaussian fp32 residual (sigma posterior mean %.3f, "
         "sse ratio %.3f)\n", sigmaPosteriorMean, sseFit / sseMean);
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

// Integer frequency weights in the logistic family are exactly row
// replication: a weight-w observation's Polya-Gamma latent is PG(w, psi),
// drawn as the sum of w PG(1, psi) variates (docs/plans/weighted-binary.md),
// so its contribution to every leaf sufficient statistic equals that of w
// identical rows. The distinguishing consequence this case pins is that the
// fitted success probability in a region tracks the WEIGHTED empirical rate
// there, sum(w*y)/sum(w) - i.e. the empirical rate of the replicated data -
// not the unweighted rate. Outcome-tied weights (successes counted
// wSuccess-fold) separate the two, so honoring the weights pulls the fit
// materially above the unweighted curve; the loglik precedent (testLogLikelihood
// case 4) and inst/tinytest/test-weighted-logistic.R cover the rest.
static void testWeightedLogistic(ext_rng* rng) {
  const size_t n = 800, p = 3;
  const double wSuccess = 4.0, wFailure = 1.0;
  std::vector<double> x(n * p), y(n), w(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double eta = 4.0 * (x[i] - 0.5);
    double probability = 1.0 / (1.0 + std::exp(-eta));
    y[i] = runif01() < probability ? 1.0 : 0.0;
    w[i] = y[i] != 0.0 ? wSuccess : wFailure;
  }

  SamplerOptions options;
  options.numTrees = 50;
  options.nodeScale = 3.0;

  // (1) determinism: identical weights and seed reproduce the draw bit for
  // bit. This is an identical-construction comparison, so it is immune to the
  // heap-layout / SIMD-reduction-split sensitivity that makes cross-sampler
  // bitwise checks (e.g. weighted vs replicated) unreliable in this harness -
  // which is why the replication equivalence below is asserted statistically.
  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 908) != 0 ||
      ext_rng_setSeed(rngB, 908) != 0) {
    check(false, "weighted logistic: rng creation");
    return;
  }
  ConstantLeafSampler a(x.data(), y.data(), n, p, w.data(), nullptr,
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngA);
  ConstantLeafSampler b(x.data(), y.data(), n, p, w.data(), nullptr,
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngB);
  const size_t detBurnIn = 20, detSamples = 30;
  std::vector<double> fitsA(n * detSamples), fitsB(n * detSamples);
  Results rA, rB;
  rA.trainingFits = fitsA.data();
  rB.trainingFits = fitsB.data();
  a.run(detBurnIn, detSamples, rA);
  b.run(detBurnIn, detSamples, rB);
  bool identical = true;
  for (size_t i = 0; i < n * detSamples && identical; ++i)
    identical = fitsA[i] == fitsB[i];
  check(identical, "weighted logistic is deterministic under a fixed seed");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);

  // (2) a full weighted fit
  ConstantLeafSampler wsampler(x.data(), y.data(), n, p, w.data(), nullptr,
                          ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rng);
  const size_t numBurnIn = 150, numSamples = 200;
  std::vector<double> wfits(n * numSamples);
  Results wres;
  wres.trainingFits = wfits.data();
  wsampler.run(numBurnIn, numSamples, wres);

  // the log-odds fit stays monotone in x1
  double lowSum = 0.0, highSum = 0.0;
  size_t lowCount = 0, highCount = 0;
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < 0.25) { lowSum += wfits[i + s * n]; ++lowCount; }
      if (x[i] > 0.75) { highSum += wfits[i + s * n]; ++highCount; }
    }
  double lowFit = lowSum / (double) lowCount;
  double highFit = highSum / (double) highCount;
  check(highFit > lowFit + 0.5, "weighted logistic recovers monotone signal");

  // the weight-specific claim: each quartile's fitted probability calibrates
  // to that quartile's WEIGHTED success rate sum(w*y)/sum(w), the rate the
  // replicated data would exhibit. Since successes are up-weighted, the
  // unweighted rate sum(y)/count is strictly lower - so calibrating to the
  // weighted rate is exactly what an honored weight does and a dropped or
  // miswired one does not.
  double lowNum = 0.0, lowDen = 0.0, highNum = 0.0, highDen = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (x[i] < 0.25) { lowNum += w[i] * y[i]; lowDen += w[i]; }
    if (x[i] > 0.75) { highNum += w[i] * y[i]; highDen += w[i]; }
  }
  double lowWeightedRate = lowNum / lowDen;
  double highWeightedRate = highNum / highDen;
  double lowProbability = 1.0 / (1.0 + std::exp(-lowFit));
  double highProbability = 1.0 / (1.0 + std::exp(-highFit));
  // 0.12 clears the ~0.1 Monte Carlo calibration error the unweighted analog
  // (testEndToEndLogistic) tolerates at this n, while staying well inside the
  // >=0.2 low-quartile gap that dropping the weights would open.
  checkNear(lowProbability, lowWeightedRate, 0.12,
            "weighted logistic calibrates to the low-quartile weighted rate");
  checkNear(highProbability, highWeightedRate, 0.12,
            "weighted logistic calibrates to the high-quartile weighted rate");

  // omega latents stay positive and finite
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
  std::vector<xint_t> codesBefore(sampler.data().train.codes);

  // identity replacement: same values in new buffers; the rebuilt cuts equal
  // the old ones, every split remaps onto itself, and fits are recovered
  // exactly
  std::vector<double> xCopy(x), yCopy(y);
  sampler.setData(xCopy.data(), yCopy.data(), n, nullptr, nullptr, nullptr, 0);
  check(sampler.data().train.codes == codesBefore, "identity setData preserves codes");
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
    codesInRange &= sampler.data().train.codes[i] <= 3;
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
  std::vector<index_t> indices(n);
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
  std::vector<index_t> rebuiltIndices(n);
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
  check(restored.setState(state, nullptr), "a pooled state restores");

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
  check(restored.setState(state, nullptr), "a state with missing directions restores");

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

// The interweaving glue-ridge rescale (docs/plans/bcf-ridge-interweaving.md):
// after a burn-in that gives the prognostic forest real leaf values, one move
// must (a) leave the combined fit a*mu + b_z*tau invariant to ~1e-12 and (b)
// keep the cached fits self-consistent - a = a0/c, treeFits scaled by exactly
// c, and totalFits still the sum of the tree slabs.
static void testBCFInterweave(ext_rng* rng) {
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
  spec.mu.numTrees = 50;
  spec.tau.numTrees = 25;
  spec.z = z.data();
  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);

  Results results;
  sampler.run(200, 0, results);  // burn in so mu has real leaf values

  size_t muTrees = spec.mu.numTrees;
  std::vector<double> mu0(n), tau0(n), treeFits0(n * muTrees);
  sampler.forestTotalFits(0, 0, mu0.data());
  sampler.forestTotalFits(0, 1, tau0.data());
  sampler.chain(0).forestTreeFits(0, treeFits0.data());
  double a0, b0, b1;
  sampler.chain(0).bcfGlue(a0, b0, b1);

  double c = sampler.chain(0).interweaveGlueRidge();
  check(c > 0.0 && std::isfinite(c), "interweave draws a positive finite c");

  std::vector<double> mu1(n), treeFits1(n * muTrees);
  sampler.forestTotalFits(0, 0, mu1.data());
  sampler.chain(0).forestTreeFits(0, treeFits1.data());
  double a1, b0p, b1p;
  sampler.chain(0).bcfGlue(a1, b0p, b1p);

  // (a) combined fit a*mu + b_z*tau invariant (tau, b0, b1 untouched)
  double maxCombinedDelta = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double before = a0 * mu0[i] + (z[i] != 0.0 ? b1 : b0) * tau0[i];
    double after = a1 * mu1[i] + (z[i] != 0.0 ? b1p : b0p) * tau0[i];
    maxCombinedDelta = std::max(maxCombinedDelta, std::fabs(after - before));
  }
  check(maxCombinedDelta < 1.0e-11, "interweave leaves the combined fit invariant");
  check(b0p == b0 && b1p == b1, "interweave leaves the b glue untouched");

  // (b) a = a0/c and each tree slab scaled by exactly c
  check(std::fabs(a1 - a0 / c) <= 1.0e-13 * std::fabs(a0 / c) + 1.0e-15,
        "interweave sets a = a0 / c");
  double maxScaleDelta = 0.0;
  for (size_t j = 0; j < n * muTrees; ++j)
    maxScaleDelta =
      std::max(maxScaleDelta, std::fabs(treeFits1[j] - c * treeFits0[j]));
  check(maxScaleDelta <= 1.0e-12, "interweave scales tree fits by c");

  // (b) totalFits still equals the sum of the tree slabs
  double maxSumDelta = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double sum = 0.0;
    for (size_t t = 0; t < muTrees; ++t) sum += treeFits1[i + t * n];
    maxSumDelta = std::max(maxSumDelta, std::fabs(sum - mu1[i]));
  }
  check(maxSumDelta < 1.0e-9, "interweave keeps totalFits the sum of tree fits");

  printf("ok: BCF interweave rescale move\n");
}

// The sharp edge (memo section 4): under keepTrees the saved mu leaves are
// flattened before the move, so the move must rescale this sweep's saved slot
// by the same c. Prediction from the saved slot reconstructs the prognostic
// total mu; scale * mu_saved + shift must therefore track scale * mu_live +
// shift, i.e. their difference is a constant shift for every row. An unscaled
// saved slot (mu_saved = mu_live / c) would make that difference row-dependent.
static void testBCFInterweaveKeepTrees(ext_rng* rng) {
  const size_t n = 300, p = 3;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  options.keepTrees = true;
  options.numSamplesToStore = 1;
  BCFSpec spec;
  spec.mu.numTrees = 40;
  spec.tau.numTrees = 20;
  spec.z = z.data();
  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);

  Results results;
  sampler.run(150, 1, results);  // one recorded sweep fills the single slot

  double scale = sampler.fitScale();
  std::vector<double> muLive(n), pred(n);
  sampler.forestTotalFits(0, 0, muLive.data());
  sampler.predict(x.data(), n, pred.data());  // scale * mu_saved + shift

  double d0 = pred[0] - scale * muLive[0];
  double maxSpread = 0.0;
  for (size_t i = 0; i < n; ++i)
    maxSpread =
      std::max(maxSpread, std::fabs((pred[i] - scale * muLive[i]) - d0));
  check(maxSpread < 1.0e-9,
        "keepTrees saved mu slot tracks the rescaled live fit after the move");

  // numThin > 1: the rescale addresses slot (savedSlotBase + sampleNum), with
  // sampleNum = iteration / numThin - numBurnIn - a thinned counter the
  // numThin = 1 case above cannot exercise. Store several slots under thinning;
  // the LAST slot is recorded on the final sweep, so its saved mu must still
  // track the post-run rescaled live fit (a misaddressed rescale would leave it
  // unscaled). A local rng plus a snapshot/restore of the shared runif01 stream
  // keep this block neutral to every downstream test's draw sequence.
  {
    uint64_t savedRngState = rngState;
    const size_t nThin = 300, pThin = 3, numSlots = 4, thin = 3;
    std::vector<double> xThin(nThin * pThin), yThin(nThin), zThin(nThin);
    for (double& v : xThin) v = runif01();
    for (size_t i = 0; i < nThin; ++i) {
      zThin[i] = runif01() < 0.5 ? 1.0 : 0.0;
      double mu = std::sin(3.0 * xThin[i]) + xThin[i + nThin];
      double tau = 1.0 + 2.0 * xThin[i + 2 * nThin];
      yThin[i] = mu + zThin[i] * tau + 0.2 * (runif01() - 0.5);
    }

    SamplerOptions optionsThin;
    optionsThin.keepTrees = true;
    optionsThin.numSamplesToStore = numSlots;
    optionsThin.numThin = thin;
    BCFSpec specThin;
    specThin.mu.numTrees = 40;
    specThin.tau.numTrees = 20;
    specThin.z = zThin.data();
    ext_rng* thinRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(thinRng, 24601u);
    Sampler<ConstantGaussianLeaf> samplerThin(
      xThin.data(), yThin.data(), nThin, pThin, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, optionsThin, specThin, &thinRng);

    Results resultsThin;
    samplerThin.run(50, numSlots, resultsThin);  // numSlots kept, thinned by 3

    double scaleThin = samplerThin.fitScale();
    std::vector<double> muLiveThin(nThin), predThin(nThin * numSlots);
    samplerThin.forestTotalFits(0, 0, muLiveThin.data());
    samplerThin.predict(xThin.data(), nThin, predThin.data());

    const double* lastSlot = predThin.data() + (numSlots - 1) * nThin;
    double dLast = lastSlot[0] - scaleThin * muLiveThin[0];
    double maxSpreadThin = 0.0;
    for (size_t i = 0; i < nThin; ++i)
      maxSpreadThin = std::max(
        maxSpreadThin,
        std::fabs((lastSlot[i] - scaleThin * muLiveThin[i]) - dLast));
    check(maxSpreadThin < 1.0e-9,
          "keepTrees numThin>1 last saved slot tracks the rescaled live fit");

    // every stored slot is a finite prognostic surface; a slot the thinned
    // addressing skipped or double-scaled would show up here
    bool slotsFinite = true;
    for (double v : predThin) slotsFinite &= std::isfinite(v);
    check(slotsFinite, "keepTrees numThin>1 stored slots stay finite");

    ext_rng_destroy(thinRng);
    rngState = savedRngState;
  }

  printf("ok: BCF interweave keepTrees saved slot\n");
}

// growForestFromRoot on a BCF sampler exercises the grow-from-root sweep's own
// two-forest branch (per-forest residual formation, the glue draw, the ridge
// interweave), which no R path reaches - growForestFromRoot is engine-internal.
// This is therefore its only gate: build the two-forest sampler, grow every
// forest from the root, and pin the combined output - both forests finite and
// off zero, glue finite, and a recorded characteristic value of the combined
// internal fit a*mu + b_z*tau. A local rng plus a snapshot/restore of the
// shared runif01 stream keep this test neutral to every downstream test's draws.
static void testBCFGrowForestFromRoot() {
  uint64_t savedRngState = rngState;
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

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 90210u);
  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &localRng);
  check(sampler.numForests() == 2, "BCF grow-from-root builds two forests");

  sampler.chain(0).growForestFromRoot(5);

  std::vector<double> muFits(n), tauFits(n);
  sampler.forestTotalFits(0, 0, muFits.data());
  sampler.forestTotalFits(0, 1, tauFits.data());
  double muSS = 0.0, tauSS = 0.0;
  bool finite = true;
  for (size_t i = 0; i < n; ++i) {
    finite &= std::isfinite(muFits[i]) && std::isfinite(tauFits[i]);
    muSS += muFits[i] * muFits[i];
    tauSS += tauFits[i] * tauFits[i];
  }
  check(finite, "BCF grow-from-root fits are finite");
  check(muSS > 0.0 && tauSS > 0.0, "both grown BCF forests move off zero");

  double a, b0, b1;
  bool haveGlue = sampler.chain(0).bcfGlue(a, b0, b1);
  check(haveGlue && std::isfinite(a) && std::isfinite(b0) &&
          std::isfinite(b1),
        "BCF grow-from-root glue is finite");

  // characteristic value: the mean combined internal fit over the grown
  // two-forest state. Deterministic given localRng seed 90210; a relocation
  // that shifts the grow sweep's draw order or the coupling moves it far past
  // the tolerance, while it survives benign cross-build FP reassociation.
  double combinedMean = 0.0;
  for (size_t i = 0; i < n; ++i)
    combinedMean += a * muFits[i] + (z[i] != 0.0 ? b1 : b0) * tauFits[i];
  combinedMean /= static_cast<double>(n);
  checkNear(combinedMean, -0.060352757243346378, 1e-6,
            "BCF grow-from-root combined fit characteristic value");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: BCF grow-from-root sweep\n");
}

// Build a one-hot n x K category-major count matrix and unit trials from
// single-trial labels - the count-native combiner's single-trial input, whose
// draw stream and working response reduce byte-for-byte to the label path.
static void oneHotCounts(const std::vector<int>& labels, size_t K,
                         std::vector<int>& counts, std::vector<int>& trials) {
  size_t n = labels.size();
  counts.assign(n * K, 0);
  trials.assign(n, 1);
  for (size_t i = 0; i < n; ++i)
    counts[static_cast<size_t>(labels[i]) * n + i] = 1;
}

// The multinomial (softmax) combiner: the interleaved one-vs-rest Polya-Gamma
// draw (positivity, its moment against the current margin, and the K = 2
// reduction to the logistic working-response construction), the K >= 3
// log-sum-exp margin, the softmax combined output and its level-centering
// invariance, the count-native path (PG(n_i) summing, the binomial(n_i, .)
// working response, the byte-identical single-trial reduction), an end-to-end
// monotone K = 3 recovery, and a structural state round-trip.
static void testMultinomial(ext_rng* rng) {
  // ---- isolated combiner: interleaved PG draw and the K = 2 reduction ----
  {
    const size_t n = 6;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);

    // K = 2: the margin C_i0 is exactly forest 1's fit, so category 0's
    // one-vs-rest is binary logistic with eta = f0 - f1 - the reduction pinned.
    std::vector<int> labels = {0, 1, 0, 1, 0, 1};
    std::vector<int> counts, trials;
    oneHotCounts(labels, 2, counts, trials);
    MultinomialSpec spec;
    spec.numCategories = 2;
    spec.counts = counts.data();
    spec.trials = trials.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

    std::vector<Forest<ConstantGaussianLeaf>> forests(2);
    std::vector<double> f0 = {-1.0, 0.5, 1.2, -0.3, 0.0, 0.8};
    std::vector<double> f1 = {0.4, -0.6, 0.2, 1.1, -0.9, 0.3};
    for (size_t k = 0; k < 2; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(n, 0.0);
    }
    forests[0].totalFits = f0;
    forests[1].totalFits = f1;

    combiner.drawForestGlue(0, rng, forests);
    ForestResponse fr =
      combiner.formForestResponse(0, forests, nullptr, nullptr);
    bool positive = true, reduction = true;
    for (size_t i = 0; i < n; ++i) {
      double yi0 = labels[i] == 0 ? 1.0 : 0.0;
      double omega = fr.weights[i];
      positive &= omega > 0.0 && std::isfinite(omega);
      // the logistic path's working response: (y - 1/2)/omega + offset, offset
      // the category margin (f1 for K = 2)
      double expected = (yi0 - 0.5) / omega + f1[i];
      reduction &= std::fabs(fr.response[i] - expected) < 1e-12;
    }
    check(positive, "multinomial interleaved omega positive and finite");
    check(reduction,
          "K = 2 multinomial reproduces the logistic working-response form");

    // the interleaved draw targets PG(1, eta_i0) against the CURRENT margin
    double eta0 = f0[0] - f1[0];
    double pgMean = std::tanh(0.5 * eta0) / (2.0 * eta0);
    const int numDraws = 40000;
    double sum = 0.0;
    for (int d = 0; d < numDraws; ++d) {
      combiner.drawForestGlue(0, rng, forests);
      ForestResponse s =
        combiner.formForestResponse(0, forests, nullptr, nullptr);
      sum += s.weights[0];
    }
    checkNear(sum / numDraws, pgMean, 0.02,
              "interleaved PG(1, eta) mean matches the theoretical moment");

    // combinedTestFits' K = 2 reduction: with the test rows' totalTestFits set
    // to the train totalFits, softmax(f0, f1)_1 = logistic(f1 - f0) and the two
    // channels sum to one per row
    std::vector<double> xTestDummy(n, 0.0);
    data.buildTest(xTestDummy.data(), n);
    forests[0].totalTestFits = f0;
    forests[1].totalTestFits = f1;
    const double* testProbs = combiner.combinedTestFits(forests);
    bool k2reduce = true;
    for (size_t i = 0; i < n; ++i) {
      double p1 = 1.0 / (1.0 + std::exp(f0[i] - f1[i]));
      k2reduce &= std::fabs(testProbs[n + i] - p1) < 1e-12 &&
                  std::fabs(testProbs[i] + testProbs[n + i] - 1.0) < 1e-12;
    }
    check(k2reduce, "combinedTestFits K = 2 reduces to the logistic probability");
  }

  // ---- isolated combiner: K = 3 log-sum-exp margin and softmax output ----
  {
    const size_t n = 4;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);
    std::vector<int> labels = {0, 1, 2, 1};
    std::vector<int> counts, trials;
    oneHotCounts(labels, 3, counts, trials);
    MultinomialSpec spec;
    spec.numCategories = 3;
    spec.counts = counts.data();
    spec.trials = trials.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

    std::vector<Forest<ConstantGaussianLeaf>> forests(3);
    std::vector<std::vector<double>> f = {
      {0.2, -0.5, 1.0, 0.1}, {-0.3, 0.8, -0.2, 0.4}, {0.5, 0.0, 0.6, -0.7}};
    for (size_t k = 0; k < 3; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(n, 0.0);
      forests[k].totalFits = f[k];
    }

    combiner.drawForestGlue(1, rng, forests);
    ForestResponse fr =
      combiner.formForestResponse(1, forests, nullptr, nullptr);
    bool marginOk = true;
    for (size_t i = 0; i < n; ++i) {
      double margin = std::log(std::exp(f[0][i]) + std::exp(f[2][i]));
      double yi1 = labels[i] == 1 ? 1.0 : 0.0;
      double expected = (yi1 - 0.5) / fr.weights[i] + margin;
      marginOk &= std::fabs(fr.response[i] - expected) < 1e-10;
    }
    check(marginOk, "K = 3 log-sum-exp margin enters the working response");

    std::vector<double> before(n * 3);
    std::memcpy(before.data(), combiner.combinedFits(forests),
                n * 3 * sizeof(double));
    bool simplex = true;
    for (size_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (size_t k = 0; k < 3; ++k) {
        double pr = before[k * n + i];
        simplex &= pr > 0.0 && pr < 1.0;
        s += pr;
      }
      simplex &= std::fabs(s - 1.0) < 1e-12;
    }
    check(simplex, "multinomial combined output is a per-observation simplex");

    // the level-centering move leaves the softmax invariant (softmax(f + delta)
    // = softmax(f)), pinning only the non-identified flat direction
    combiner.afterCombine(forests, false, 0, rng);
    const double* after = combiner.combinedFits(forests);
    bool invariant = true;
    for (size_t i = 0; i < n * 3; ++i)
      invariant &= std::fabs(after[i] - before[i]) < 1e-10;
    check(invariant,
          "level-centering leaves the softmax probabilities invariant");

    // combinedTestFits is combinedFits' totalTestFits analog: with the test
    // rows' totalTestFits set to the pre-shift train totalFits `f`, it must be a
    // per-observation simplex AND equal the pre-shift combinedFits (`before`) -
    // the softmax-invariance fact that lets afterCombine leave totalTestFits
    // untouched (the level shift c is common to all K forests)
    std::vector<double> xTestDummy(n, 0.0);
    data.buildTest(xTestDummy.data(), n);
    for (size_t k = 0; k < 3; ++k) forests[k].totalTestFits = f[k];
    const double* testProbs = combiner.combinedTestFits(forests);
    bool testSimplex = true, testMatches = true;
    for (size_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (size_t k = 0; k < 3; ++k) {
        double pr = testProbs[k * n + i];
        testSimplex &= pr > 0.0 && pr < 1.0;
        s += pr;
        testMatches &= std::fabs(pr - before[k * n + i]) < 1e-12;
      }
      testSimplex &= std::fabs(s - 1.0) < 1e-12;
    }
    check(testSimplex, "combinedTestFits is a per-observation test simplex");
    check(testMatches,
          "combinedTestFits equals combinedFits when test rows equal train");
  }

  // ---- count-native combiner: the binomial(n_i, .) working response and the
  //      PG(n_i) = sum-of-n_i-PG(1) mean moment (K = 2, n_i > 1) ----
  {
    const size_t n = 5;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);

    // K = 2 grouped counts: category 0's one-vs-rest is binomial(n_i, .) with
    // margin exactly forest 1's fit (eta = f0 - f1). Category-major counts.
    std::vector<int> trials = {1, 2, 3, 4, 5};
    std::vector<int> y0 = {1, 1, 2, 0, 3};  // category-0 successes, y0 <= trials
    std::vector<int> counts(n * 2);
    for (size_t i = 0; i < n; ++i) {
      counts[i] = y0[i];
      counts[n + i] = trials[i] - y0[i];
    }
    MultinomialSpec spec;
    spec.numCategories = 2;
    spec.counts = counts.data();
    spec.trials = trials.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

    std::vector<Forest<ConstantGaussianLeaf>> forests(2);
    std::vector<double> f0 = {-0.4, 0.3, 0.9, -0.1, 0.5};
    std::vector<double> f1 = {0.2, -0.5, 0.1, 0.7, -0.3};
    for (size_t k = 0; k < 2; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(n, 0.0);
    }
    forests[0].totalFits = f0;
    forests[1].totalFits = f1;

    combiner.drawForestGlue(0, rng, forests);
    ForestResponse fr = combiner.formForestResponse(0, forests, nullptr, nullptr);
    bool positive = true, binomialForm = true;
    for (size_t i = 0; i < n; ++i) {
      double omega = fr.weights[i];
      positive &= omega > 0.0 && std::isfinite(omega);
      // binomial(n_i, .) working response: (y - n_i/2)/omega + margin (f1)
      double expected =
        (static_cast<double>(y0[i]) - trials[i] * 0.5) / omega + f1[i];
      binomialForm &= std::fabs(fr.response[i] - expected) < 1e-12;
    }
    check(positive, "count-native omega positive and finite");
    check(binomialForm,
          "K = 2 count fit forms the binomial(n_i, .) working response");

    // PG(n, eta) is a sum of n PG(1, eta), so its mean is n times the PG(1)
    // moment; match it on a fixed row (weighted-logistic's moment style)
    const size_t row = 4;  // 5 trials
    double eta = f0[row] - f1[row];
    double pgnMean = static_cast<double>(trials[row]) *
                     std::tanh(0.5 * eta) / (2.0 * eta);
    const int numDraws = 40000;
    double sum = 0.0;
    for (int d = 0; d < numDraws; ++d) {
      combiner.drawForestGlue(0, rng, forests);
      ForestResponse s = combiner.formForestResponse(0, forests, nullptr, nullptr);
      sum += s.weights[row];
    }
    checkNear(sum / numDraws, pgnMean, 0.03,
              "PG(n, eta) mean matches n times the PG(1, eta) moment");
  }

  // ---- the byte-identical single-trial reduction: a one-hot count fit draws
  //      the SAME PG stream and forms the SAME working response as the label
  //      path, bit-for-bit at a fixed seed (the in-process neutrality proof) ----
  {
    const size_t n = 8, K = 3;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);
    std::vector<int> labels = {0, 2, 1, 0, 2, 1, 0, 1};
    std::vector<int> counts, trials;
    oneHotCounts(labels, K, counts, trials);
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

    std::vector<Forest<ConstantGaussianLeaf>> forests(K);
    std::vector<std::vector<double>> f = {
      {0.3, -0.4, 0.8, 0.1, -0.6, 0.5, 0.2, -0.1},
      {-0.2, 0.6, -0.3, 0.4, 0.7, -0.5, 0.0, 0.3},
      {0.5, 0.1, -0.7, 0.2, -0.1, 0.6, -0.4, 0.8}};
    for (size_t k = 0; k < K; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(n, 0.0);
      forests[k].totalFits = f[k];
    }

    // Two rngs seeded identically: the combiner (count path) and a hand-rolled
    // reference replaying the label path (one PG(1, psi) per observation,
    // psi = f_if - C_if, working response (indicator - 0.5)/omega + C_if). Same
    // per-forest draw order, so equal draw streams; every omega and response
    // must be bitwise equal, not merely close - the n_i = 1 reduction. The
    // margin mirrors the engine's pairwise log-sum-exp (symmetric, so at K = 3
    // with static fits the prefix/suffix mix is one merge of the two others).
    ext_rng* refRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(refRng, 13579u);
    ext_rng* combRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(combRng, 13579u);

    bool reductionExact = true;
    for (size_t fk = 0; fk < K; ++fk) {
      combiner.drawForestGlue(fk, combRng, forests);
      ForestResponse fr =
        combiner.formForestResponse(fk, forests, nullptr, nullptr);
      for (size_t i = 0; i < n; ++i) {
        double other[2];
        size_t c = 0;
        for (size_t j = 0; j < K; ++j)
          if (j != fk) other[c++] = f[j][i];
        double hi = std::max(other[0], other[1]);
        double lo = std::min(other[0], other[1]);
        double margin = hi + std::log1p(std::exp(lo - hi));
        double omegaRef = ext_rng_simulatePolyaGamma(refRng, f[fk][i] - margin);
        double yif = labels[i] == static_cast<int>(fk) ? 1.0 : 0.0;
        double respRef = (yif - 0.5) / omegaRef + margin;
        reductionExact &=
          fr.weights[i] == omegaRef && fr.response[i] == respRef;
      }
    }
    check(reductionExact,
          "single-trial one-hot count fit reduces to the label path bit-for-bit");
    ext_rng_destroy(refRng);
    ext_rng_destroy(combRng);
  }

  // ---- end-to-end: recover a monotone K = 3 signal ----
  {
    const size_t n = 600, p = 1, K = 3;
    std::vector<double> x(n);
    std::vector<int> labels(n);
    for (size_t i = 0; i < n; ++i) {
      double xi = runif01();
      x[i] = xi;
      double e0 = 2.0 * (0.5 - xi), e2 = 2.0 * (xi - 0.5);  // e1 = 0
      double m = std::max(e0, std::max(0.0, e2));
      double s = std::exp(e0 - m) + std::exp(-m) + std::exp(e2 - m);
      double p0 = std::exp(e0 - m) / s, p1 = std::exp(-m) / s;
      double u = runif01();
      labels[i] = u < p0 ? 0 : (u < p0 + p1 ? 1 : 2);
    }

    std::vector<int> counts, trials;
    oneHotCounts(labels, K, counts, trials);
    SamplerOptions options;
    options.numTrees = 50;
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 50;

    Sampler<ConstantGaussianLeaf> sampler(x.data(), n, p, options, spec, &rng);
    check(sampler.numForests() == K, "multinomial builds K forests");
    check(sampler.numReportedLocations() == K,
          "multinomial reports K locations");
    check(sampler.numVariableCountForests() == K,
          "multinomial reports K variable-count forests");

    const size_t numBurnIn = 200, numSamples = 200;
    std::vector<double> fits(n * K * numSamples);
    std::vector<std::uint32_t> varcounts(p * K * numSamples);
    Results results;
    results.trainingFits = fits.data();
    results.numReportedLocations = K;
    results.variableCounts = varcounts.data();
    results.numVariableCountForests = K;
    sampler.run(numBurnIn, numSamples, results);

    std::vector<double> meanProb(n * K, 0.0);
    bool valid = true;
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t k = 0; k < K; ++k)
        for (size_t i = 0; i < n; ++i) {
          double pr = fits[i + n * (k + K * s)];
          valid &= pr >= 0.0 && pr <= 1.0;
          meanProb[i + n * k] += pr / numSamples;
        }
    check(valid, "multinomial probabilities are in [0, 1]");

    double lowP0 = 0.0, lowP2 = 0.0, highP0 = 0.0, highP2 = 0.0;
    size_t lowN = 0, highN = 0;
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < 0.25) {
        lowP0 += meanProb[i]; lowP2 += meanProb[i + 2 * n]; ++lowN;
      } else if (x[i] > 0.75) {
        highP0 += meanProb[i]; highP2 += meanProb[i + 2 * n]; ++highN;
      }
    }
    check(lowP0 / lowN > lowP2 / lowN, "category 0 dominates at low x");
    check(highP2 / highN > highP0 / highN, "category 2 dominates at high x");

    // per-category varcount: the widened storeSample write lays K forest-major
    // varcount slabs (p each) per sample, slot k = category k's forest. The
    // last sample's slab k must equal forestVariableCounts(k) on the live trees
    // (storeSample is a sweep's final act, so the trees are unchanged since it
    // wrote the counts), the K-forest generalization of the single-forest case.
    bool perCategoryVarcount = true;
    std::vector<std::uint32_t> liveCounts(p);
    for (size_t k = 0; k < K; ++k) {
      sampler.forestVariableCounts(0, k, liveCounts.data());
      const std::uint32_t* slab =
        varcounts.data() + ((numSamples - 1) * K + k) * p;
      for (size_t j = 0; j < p; ++j)
        perCategoryVarcount &= slab[j] == liveCounts[j];
    }
    check(perCategoryVarcount,
          "multinomial storeSample writes each category forest's varcount");

    // structural state round-trip: the K forests serialize through the
    // per-forest list, and the softmax combiner carries no wire state (omega is
    // redrawn on the first restored sweep)
    SamplerStateData state;
    sampler.getState(state);
    check(state.chains[0].forests.size() == K,
          "multinomial state serializes K forests");
    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngB, 4321);
    Sampler<ConstantGaussianLeaf> restored(x.data(), n, p, options, spec, &rngB);
    check(restored.setState(state, nullptr), "multinomial state restores");
    checkStructuralRoundTrip(state, restored,
                             "multinomial restore reproduces the K-forest model");
    ext_rng_destroy(rngB);
  }

  // ---- keepTrees: K-forest saved-tree replay reproduces the recorded test
  // channel bit for bit. predict sums every forest's saved trees at the test
  // rows and softmaxes them; the run accumulated the same forests' totalTestFits
  // and blended them the same way, and neither carries the level-centering shift
  // (afterCombine leaves totalTestFits and the saved leaves alone), so the two
  // agree exactly. A local rng plus a runif01 snapshot keep this neutral to the
  // shared streams the sequenced tests draw from. ----
  {
    uint64_t savedRngState = rngState;
    const size_t n = 200, p = 2, K = 3, numSamples = 6;
    std::vector<double> x(n * p);
    std::vector<int> labels(n);
    for (double& v : x) v = runif01();
    for (size_t i = 0; i < n; ++i) {
      double e0 = 1.5 * (x[i] - 0.5), e1 = x[i + n] - 0.5,
             e2 = 0.5 * (x[i] - x[i + n]);
      double m = std::max(e0, std::max(e1, e2));
      double s = std::exp(e0 - m) + std::exp(e1 - m) + std::exp(e2 - m);
      double p0 = std::exp(e0 - m) / s, p1 = std::exp(e1 - m) / s;
      double u = runif01();
      labels[i] = u < p0 ? 0 : (u < p0 + p1 ? 1 : 2);
    }
    // test rows: a held-out slice of the training design
    const size_t nTest = 25;
    std::vector<double> xTest(nTest * p);
    for (size_t i = 0; i < nTest; ++i)
      for (size_t j = 0; j < p; ++j) xTest[i + j * nTest] = x[i + j * n];

    std::vector<int> counts, trials;
    oneHotCounts(labels, K, counts, trials);
    SamplerOptions options;
    options.numTrees = 30;
    options.keepTrees = true;
    options.numSamplesToStore = numSamples;
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 30;

    ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(localRng, 20260716u);
    Sampler<ConstantGaussianLeaf> sampler(x.data(), n, p, options, spec,
                                          &localRng);
    sampler.setTestPredictors(xTest.data(), nTest);

    std::vector<double> testFits(nTest * K * numSamples);
    Results results;
    results.testFits = testFits.data();
    results.numReportedLocations = K;
    sampler.run(80, numSamples, results);

    std::vector<double> predicted(nTest * K * numSamples);
    sampler.predict(xTest.data(), nTest, predicted.data());
    check(predicted == testFits,
          "multinomial K-forest saved-tree replay equals the recorded test "
          "channel");

    bool simplex = true;
    for (size_t sm = 0; sm < numSamples; ++sm)
      for (size_t i = 0; i < nTest; ++i) {
        double sum = 0.0;
        for (size_t k = 0; k < K; ++k) {
          double pr = predicted[i + nTest * (k + K * sm)];
          simplex &= pr > 0.0 && pr < 1.0;
          sum += pr;
        }
        simplex &= std::fabs(sum - 1.0) < 1e-12;
      }
    check(simplex, "multinomial replayed predictions are per-row simplices");

    ext_rng_destroy(localRng);
    rngState = savedRngState;
  }

  printf("ok: multinomial softmax sampler\n");
}

// growForestFromRoot on a multinomial sampler exercises the grow-from-root
// sweep's own combiner branch (the interleaved per-forest PG draw, the
// per-forest working response, and the level-centering move), which no R path
// reaches - growForestFromRoot is engine-internal (the BCF lesson). A local rng
// plus a snapshot/restore of the shared runif01 stream keep this test neutral.
static void testMultinomialGrowForestFromRoot() {
  uint64_t savedRngState = rngState;
  const size_t n = 400, p = 2, K = 3;
  std::vector<double> x(n * p);
  std::vector<int> labels(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double e0 = 1.5 * x[i], e1 = -1.0 + x[i + n], e2 = 0.5 * (x[i] - x[i + n]);
    double m = std::max(e0, std::max(e1, e2));
    double s = std::exp(e0 - m) + std::exp(e1 - m) + std::exp(e2 - m);
    double p0 = std::exp(e0 - m) / s, p1 = std::exp(e1 - m) / s;
    double u = runif01();
    labels[i] = u < p0 ? 0 : (u < p0 + p1 ? 1 : 2);
  }

  std::vector<int> counts, trials;
  oneHotCounts(labels, K, counts, trials);
  SamplerOptions options;
  options.numTrees = 40;
  MultinomialSpec spec;
  spec.numCategories = K;
  spec.counts = counts.data();
  spec.trials = trials.data();
  spec.forest.numTrees = 40;

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 90210u);
  Sampler<ConstantGaussianLeaf> sampler(x.data(), n, p, options, spec,
                                        &localRng);
  check(sampler.numForests() == K, "multinomial grow-from-root builds K forests");

  sampler.chain(0).growForestFromRoot(5);

  std::vector<std::vector<double>> fits(K, std::vector<double>(n));
  bool finite = true, offZero = true;
  for (size_t k = 0; k < K; ++k) {
    sampler.forestTotalFits(0, k, fits[k].data());
    double ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
      finite &= std::isfinite(fits[k][i]);
      ss += fits[k][i] * fits[k][i];
    }
    offZero &= ss > 0.0;
  }
  check(finite, "multinomial grow-from-root fits are finite");
  check(offZero, "every grown multinomial forest moves off zero");

  bool simplex = true;
  for (size_t i = 0; i < n; ++i) {
    double m = fits[0][i];
    for (size_t k = 1; k < K; ++k) m = std::max(m, fits[k][i]);
    double s = 0.0;
    for (size_t k = 0; k < K; ++k) s += std::exp(fits[k][i] - m);
    double total = 0.0;
    for (size_t k = 0; k < K; ++k) total += std::exp(fits[k][i] - m) / s;
    simplex &= std::fabs(total - 1.0) < 1e-10;
  }
  check(simplex, "grown multinomial softmax is a per-observation simplex");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: multinomial grow-from-root sweep\n");
}

// growForestFromRoot on a GROUPED-COUNT multinomial sampler: the count branch's
// PG(n_i) summing draw and (y - n_i/2) working response are unreachable from R
// through grow-from-root, so the engine test is the only cover (the G4 lesson).
// A local rng plus a snapshot/restore of the shared runif01 stream keep this
// neutral to the sequenced tests.
static void testMultinomialCountGrowForestFromRoot() {
  uint64_t savedRngState = rngState;
  const size_t n = 300, p = 2, K = 3;
  std::vector<double> x(n * p);
  std::vector<int> counts(n * K, 0), trials(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double e0 = 1.5 * x[i], e1 = -1.0 + x[i + n], e2 = 0.5 * (x[i] - x[i + n]);
    double m = std::max(e0, std::max(e1, e2));
    double s = std::exp(e0 - m) + std::exp(e1 - m) + std::exp(e2 - m);
    double pr0 = std::exp(e0 - m) / s, pr1 = std::exp(e1 - m) / s;
    size_t ni = 2 + static_cast<size_t>(runif01() * 4.0);  // 2..5 trials
    trials[i] = static_cast<int>(ni);
    for (size_t t = 0; t < ni; ++t) {
      double u = runif01();
      size_t k = u < pr0 ? 0 : (u < pr0 + pr1 ? 1 : 2);
      ++counts[k * n + i];
    }
  }

  SamplerOptions options;
  options.numTrees = 40;
  MultinomialSpec spec;
  spec.numCategories = K;
  spec.counts = counts.data();
  spec.trials = trials.data();
  spec.forest.numTrees = 40;

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 24680u);
  Sampler<ConstantGaussianLeaf> sampler(x.data(), n, p, options, spec, &localRng);
  check(sampler.numForests() == K,
        "multinomial count grow-from-root builds K forests");

  sampler.chain(0).growForestFromRoot(5);

  std::vector<std::vector<double>> fits(K, std::vector<double>(n));
  bool finite = true, offZero = true, simplex = true;
  for (size_t k = 0; k < K; ++k) {
    sampler.forestTotalFits(0, k, fits[k].data());
    double ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
      finite &= std::isfinite(fits[k][i]);
      ss += fits[k][i] * fits[k][i];
    }
    offZero &= ss > 0.0;
  }
  for (size_t i = 0; i < n; ++i) {
    double m = fits[0][i];
    for (size_t k = 1; k < K; ++k) m = std::max(m, fits[k][i]);
    double se = 0.0;
    for (size_t k = 0; k < K; ++k) se += std::exp(fits[k][i] - m);
    double total = 0.0;
    for (size_t k = 0; k < K; ++k) total += std::exp(fits[k][i] - m) / se;
    simplex &= std::fabs(total - 1.0) < 1e-10;
  }
  check(finite, "multinomial count grow-from-root fits are finite");
  check(offZero, "every grown count multinomial forest moves off zero");
  check(simplex, "grown count multinomial softmax is a per-observation simplex");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: multinomial count grow-from-root sweep\n");
}

// The per-observation log-likelihood channel: requesting it draws no rng and
// mutates no state (computed post-hoc at storeSample), so sigma/train are
// bitwise unchanged, and each family's values equal the closed-form density of
// the recorded fit and sigma on a tiny fixed input.
static void testLogLikelihood() {
  const size_t n = 200, p = 3;
  std::vector<double> x(n * p), yGaussian(n), yBinary(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double eta = 2.0 * (x[i] - 0.5) + x[i + n];
    yGaussian[i] = eta + 0.3 * (runif01() - 0.5);
    yBinary[i] =
      runif01() < 0.5 * std::erfc(-eta / std::sqrt(2.0)) ? 1.0 : 0.0;
  }

  SamplerOptions options;
  options.numTrees = 25;
  const size_t numBurnIn = 30, numSamples = 20;

  // (1) NULL costs nothing: an identical seeded run with the channel off draws
  // the same sigma/train as one with it on
  auto runGaussian = [&](bool withLogLik, std::vector<double>& sigma,
                         std::vector<double>& fits,
                         std::vector<double>& loglik) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 4242);
    ConstantLeafSampler sampler(x.data(), yGaussian.data(), n, p, nullptr,
                                nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, options, &rng);
    sigma.assign(numSamples, 0.0);
    fits.assign(n * numSamples, 0.0);
    loglik.assign(n * numSamples, 0.0);
    Results results;
    results.sigma = sigma.data();
    results.trainingFits = fits.data();
    if (withLogLik) results.logLikelihood = loglik.data();
    sampler.run(numBurnIn, numSamples, results);
    ext_rng_destroy(rng);
  };

  std::vector<double> sigmaOff, fitsOff, llOff, sigmaOn, fitsOn, llOn;
  runGaussian(false, sigmaOff, fitsOff, llOff);
  runGaussian(true, sigmaOn, fitsOn, llOn);
  check(sigmaOff == sigmaOn && fitsOff == fitsOn,
        "requesting logLikelihood leaves sigma and train bitwise unchanged");

  // (2) gaussian: dnorm(y, fit, sigma, log)
  const double logSqrt2Pi = 0.5 * std::log(2.0 * 3.141592653589793);
  bool gaussianOk = true;
  for (size_t s = 0; s < numSamples && gaussianOk; ++s)
    for (size_t i = 0; i < n; ++i) {
      double z = (yGaussian[i] - fitsOn[i + s * n]) / sigmaOn[s];
      double expected = -logSqrt2Pi - std::log(sigmaOn[s]) - 0.5 * z * z;
      if (std::fabs(llOn[i + s * n] - expected) > 1e-9) {
        gaussianOk = false;
        break;
      }
    }
  check(gaussianOk, "gaussian logLikelihood equals the normal log density");

  // (3) probit: log dbinom(y, 1, Phi(eta)); the fit is the latent location
  {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 909);
    ConstantLeafSampler sampler(x.data(), yBinary.data(), n, p, nullptr,
                                nullptr, ResponseFamily::probit, 1.0, 3.0, 1.0,
                                options, &rng);
    std::vector<double> fits(n * numSamples), loglik(n * numSamples);
    Results results;
    results.trainingFits = fits.data();
    results.logLikelihood = loglik.data();
    sampler.run(numBurnIn, numSamples, results);
    ext_rng_destroy(rng);

    const double sqrtHalf = 1.0 / std::sqrt(2.0);
    bool probitOk = true;
    for (size_t s = 0; s < numSamples && probitOk; ++s)
      for (size_t i = 0; i < n; ++i) {
        double eta = fits[i + s * n];
        double expected = yBinary[i] != 0.0
          ? std::log(0.5 * std::erfc(-eta * sqrtHalf))
          : std::log(0.5 * std::erfc(eta * sqrtHalf));
        if (std::fabs(loglik[i + s * n] - expected) > 1e-8) {
          probitOk = false;
          break;
        }
      }
    check(probitOk, "probit logLikelihood equals log dbinom(y, 1, Phi(eta))");
  }

  // (4) weighted logistic: w * log dbinom(y, 1, plogis(eta))
  {
    std::vector<double> w(n);
    for (size_t i = 0; i < n; ++i) w[i] = static_cast<double>(1 + (i % 3));
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 707);
    // the constructor takes (x, y, n, p, weights, offset, family, ...)
    ConstantLeafSampler sampler(x.data(), yBinary.data(), n, p, w.data(),
                                nullptr, ResponseFamily::logistic, 1.0, 3.0,
                                1.0, options, &rng);
    std::vector<double> fits(n * numSamples), loglik(n * numSamples);
    Results results;
    results.trainingFits = fits.data();
    results.logLikelihood = loglik.data();
    sampler.run(numBurnIn, numSamples, results);
    ext_rng_destroy(rng);

    // stable per-trial log-prob: -log(1 + exp(-eta)) for a success,
    // -log(1 + exp(eta)) for a failure; the weighted channel is w times it
    bool logisticOk = true;
    for (size_t s = 0; s < numSamples && logisticOk; ++s)
      for (size_t i = 0; i < n; ++i) {
        double eta = fits[i + s * n];
        if (std::fabs(eta) > 30.0) continue; // avoid exp overflow in the ref
        double perTrial = -std::log1p(std::exp(yBinary[i] != 0.0 ? -eta : eta));
        double expected = w[i] * perTrial;
        if (std::fabs(loglik[i + s * n] - expected) > 1e-8) {
          logisticOk = false;
          break;
        }
      }
    check(logisticOk,
          "weighted logistic logLikelihood equals w * log dbinom(y, 1, p)");
  }

  printf("ok: per-observation log-likelihood channel\n");
}

// Per-forest column restriction: a forest handed a column subset must never
// propose or accept a split outside it over a real MH run.
static void testForestColumnRestriction(ext_rng* rng) {
  const size_t n = 400, p = 5;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  // signal in every column, strongest in the columns we forbid, so an
  // unrestricted forest would certainly split on them
  for (size_t i = 0; i < n; ++i) {
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = 3.0 * x[i + 1 * n] + 3.0 * x[i + 3 * n] + x[i] + x[i + 2 * n] +
           x[i + 4 * n] + 0.2 * normal;
  }

  const std::vector<size_t> allowed = {0, 2};
  std::vector<bool> isAllowed(p, false);
  for (size_t c : allowed) isAllowed[c] = true;

  SamplerOptions options;
  options.numTrees = 40;
  options.forestColumns = allowed.data();
  options.numForestColumns = allowed.size();
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(60, 60, empty);

  SamplerStateData state;
  sampler.getState(state);
  size_t splits = 0;
  bool contained = true;
  for (const std::vector<FlatNode>& tree : state.chains[0].forests[0].trees)
    for (const FlatNode& node : tree)
      if (node.variable != invalidVariable) {
        ++splits;
        if (!isAllowed[static_cast<size_t>(node.variable)]) contained = false;
      }
  check(contained, "restricted forest never splits outside its column subset");
  check(splits > 0, "restricted forest still splits on its allowed columns");

  printf("ok: forest column restriction (%lu splits, all within subset)\n",
         static_cast<unsigned long>(splits));
}

// A restriction naming every column installs a non-null all-ones mask; its
// draws must be byte-identical to the unrestricted (null-mask) default, so the
// masking arithmetic perturbs nothing when it clears no bit. The default path's
// own bitwise identity against the pre-change engine is the equivalence gate.
static void testForestColumnRestrictionAllNeutral() {
  const size_t n = 300, p = 4;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = std::sin(3.0 * x[i]) + x[i + n] + 0.5 * runif01();

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 918) != 0 ||
      ext_rng_setSeed(rngB, 918) != 0) {
    check(false, "column-restriction neutrality: rng creation");
    return;
  }

  SamplerOptions optionsDefault;
  optionsDefault.numTrees = 30;
  ConstantLeafSampler unrestricted(x.data(), y.data(), n, p, nullptr, nullptr,
                                   ResponseFamily::gaussian, 1.0, 3.0,
                                   0.37804942330213542, optionsDefault, &rngA);

  const std::vector<size_t> allColumns = {0, 1, 2, 3};
  SamplerOptions optionsAll = optionsDefault;
  optionsAll.forestColumns = allColumns.data();
  optionsAll.numForestColumns = allColumns.size();
  ConstantLeafSampler fullMask(x.data(), y.data(), n, p, nullptr, nullptr,
                               ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, optionsAll, &rngB);

  const size_t numBurnIn = 40, numSamples = 40;
  std::vector<double> sigmaA(numSamples), sigmaB(numSamples);
  std::vector<double> fitsA(n * numSamples), fitsB(n * numSamples);
  Results resultsA, resultsB;
  resultsA.sigma = sigmaA.data();
  resultsA.trainingFits = fitsA.data();
  resultsB.sigma = sigmaB.data();
  resultsB.trainingFits = fitsB.data();
  unrestricted.run(numBurnIn, numSamples, resultsA);
  fullMask.run(numBurnIn, numSamples, resultsB);

  bool identical = true;
  for (size_t s = 0; s < numSamples && identical; ++s)
    identical = sigmaA[s] == sigmaB[s];
  for (size_t i = 0; i < n * numSamples && identical; ++i)
    identical = fitsA[i] == fitsB[i];
  check(identical, "all-columns mask draws byte-match the unrestricted default");

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  printf("ok: forest column restriction neutral when it clears nothing\n");
}

// A BCF treatment forest handed a moderator subset must never split outside it,
// while the prognostic forest keeps reading the full store. Signal sits in a
// column tau is forbidden to use, so an unrestricted tau would split there.
static void testBCFTauModeratorRestriction(ext_rng* rng) {
  const size_t n = 400, p = 4;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    double mu = 3.0 * x[i + 2 * n] + x[i + 1 * n];  // strongest in col 2
    double tau = 2.0 * x[i] + 3.0 * x[i + 2 * n];   // col 0 allowed, col 2 not
    y[i] = mu + z[i] * tau + 0.2 * (runif01() - 0.5);
  }

  const std::vector<size_t> allowed = {0};
  SamplerOptions options;
  BCFSpec spec;
  spec.mu.numTrees = 50; spec.mu.base = 0.95; spec.mu.power = 2.0;
  spec.tau.numTrees = 40; spec.tau.base = 0.95; spec.tau.power = 2.0;
  spec.tau.columns = allowed.data();
  spec.tau.numColumns = allowed.size();
  spec.z = z.data();

  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);
  Results empty;
  sampler.run(120, 120, empty);

  SamplerStateData state;
  sampler.getState(state);
  size_t tauSplits = 0;
  bool tauContained = true;
  for (const std::vector<FlatNode>& tree : state.chains[0].forests[1].trees)
    for (const FlatNode& node : tree)
      if (node.variable != invalidVariable) {
        ++tauSplits;
        if (node.variable != 0) tauContained = false;
      }
  check(tauContained, "BCF tau never splits outside its moderator subset");
  check(tauSplits > 0, "BCF tau still splits on its allowed moderators");

  bool muUsesForbidden = false;
  for (const std::vector<FlatNode>& tree : state.chains[0].forests[0].trees)
    for (const FlatNode& node : tree)
      if (node.variable != invalidVariable && node.variable != 0)
        muUsesForbidden = true;
  check(muUsesForbidden, "BCF mu reads the full store while tau is restricted");

  printf("ok: BCF tau moderator restriction (%lu tau splits, all in subset)\n",
         static_cast<unsigned long>(tauSplits));
}

void runSamplerTests(ext_rng* rng) {
  testForestColumnRestriction(rng);
  testForestColumnRestrictionAllNeutral();
  testBCFTauModeratorRestriction(rng);
  testBCFTwoForest(rng);
  testBCFFixedGlue(rng);
  testBCFInterweave(rng);
  testBCFInterweaveKeepTrees(rng);
  testBCFGrowForestFromRoot();
  testMultinomial(rng);
  testMultinomialGrowForestFromRoot();
  testMultinomialCountGrowForestFromRoot();
  testViewSamplerMatchesFull();
  testEndToEndGaussian(rng);
  testEndToEndGaussianFp32(rng);
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
  testLogLikelihood();
}
