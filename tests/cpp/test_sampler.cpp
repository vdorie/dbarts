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

// The constant-leaf mu[leafOf] gathers in rollTreeResidual and
// finalizeTotalFits are unrolled by 4 behind an n % 4 prologue, so only a shape
// whose n is not a multiple of 4 runs both halves; every other fixture here uses
// a round n. Two identities pin the loops elementwise, tail included:
//
//   totalFits[i] == sum_t mu_t[leafOf_t[i]]        (finalizeTotalFits)
//   treeY[i]     == y[i] - sum_{t < last} mu_t[leafOf_t[i]]   (the roll)
//
// The right-hand sides are gathered fresh through the dense export, which stays
// rolled, so a prologue that skipped or double-counted an element shows up on
// the left. The roll accumulates incrementally where the reference sums in one
// pass, so the comparison is to within accumulation error, not bitwise.
template <typename SamplerT>
static void checkGatherTailIdentities(SamplerT& sampler, size_t n,
                                      double tolerance, const char* fitLabel,
                                      const char* residualLabel) {
  size_t numTrees = sampler.chain(0).numTrees();
  const std::vector<double>& total = sampler.chain(0).totalFits();
  std::vector<double> fits = sampler.chain(0).treeFits();
  const double* y = sampler.chain(0).workingResponseForTesting();
  const auto& resid = sampler.chain(0).residualForTesting();

  double worstFit = 0.0, worstResidual = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double allButLast = 0.0;
    for (size_t t = 0; t + 1 < numTrees; ++t) allButLast += fits[t * n + i];
    double gathered = allButLast + fits[(numTrees - 1) * n + i];
    worstFit = std::max(worstFit, std::fabs(total[i] - gathered));
    worstResidual =
      std::max(worstResidual,
               std::fabs(static_cast<double>(resid[i]) - (y[i] - allButLast)));
  }
  check(worstFit < tolerance, fitLabel);
  check(worstResidual < tolerance, residualLabel);
  printf("ok: %s (worst fit %.3g, worst residual %.3g)\n", fitLabel, worstFit,
         worstResidual);
}

static void testGatherTailShapes(ext_rng* rng) {
  const size_t n = 1001, p = 4;  // n % 4 == 1: a one-element prologue
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = std::sin(3.0 * x[i]) + 2.0 * x[i + n] * x[i + 2 * n] +
           0.3 * (runif01() - 0.5);

  SamplerOptions options;
  options.numTrees = 25;

  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results results;
  sampler.run(20, 5, results);
  checkGatherTailIdentities(sampler, n, 1e-10,
                            "n % 4 != 0: totalFits matches a fresh mu[leafOf] sum",
                            "n % 4 != 0: rolled residual matches a scalar sum");

  // the fp32 instantiation keeps the rolled loops; it runs the same shape so a
  // future edit that unrolls it too is caught by the same identities, at the
  // looser tolerance a float residual store carries
  options.fp32Residual = true;
  Sampler<ConstantGaussianLeaf, float> single(
    x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
    3.0, 0.37804942330213542, options, &rng);
  Results fp32Results;
  single.run(20, 5, fp32Results);
  checkGatherTailIdentities(single, n, 1e-5,
                            "fp32 n % 4 != 0: totalFits matches a fresh mu[leafOf] sum",
                            "fp32 n % 4 != 0: rolled residual matches a scalar sum");
}

// The fused suffstat's K banks combine strictly left to right,
// ((b0 + b1) + b2) + b3, which is part of the draw law rather than an
// implementation detail. These banks discriminate: 1e-16 is below half an ulp
// of 1.0 (2.22e-16), so left to right each small term rounds away and the sum
// is exactly 1.0, while any pairwise regrouping sums the three small terms
// first and returns 1.0000000000000002. Asserted bitwise, not near.
static void testFusedSuffstatBankCombine() {
  const double contiguous[fusedSuffstatBanks] = {1.0, 1e-16, 1e-16, 1e-16};
  check(combineFusedSuffstatBanks(contiguous, 0, 1) == 1.0,
        "fused bank combine associates left to right");
  // the kernel lays the banks out bank-major: slot s of bank b is at
  // b * stride + s, with stride the tree's node count
  const double strided[2 * fusedSuffstatBanks] = {0.0, 1.0,   0.0, 1e-16,
                                                  0.0, 1e-16, 0.0, 1e-16};
  check(combineFusedSuffstatBanks(strided, 1, 2) == 1.0,
        "fused bank combine reads the banks bank-major");
  printf("ok: fused suffstat bank combine\n");
}

// The fused roll + suffstat pass reproduces rollTreeResidual's per-element
// residual BITWISE - same expressions, same order - which is what makes the
// suffstat's summation association the only draw change the fusion
// introduces. The statistic itself legitimately differs in the last ulps,
// since the fused banks and the stock unroll-by-5 gather associate
// differently, so it is compared to a mixed absolute/relative bound.
//
// n sweeps all four residues of the n % 4 prologue, which puts its elements in
// bank 0 and shifts every later element's bank - the pass's one new failure
// mode, the same reasoning testGatherTailShapes applies to the unrolled
// gathers - plus the short shapes where n < 4 makes the whole vector prologue
// and K collapses to 1.
static void checkFusedSuffstatShape(size_t n, size_t numBurnIn, ext_rng* rng) {
  const size_t p = 3;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i)
    y[i] = std::sin(3.0 * x[i]) + 2.0 * x[i + n] * x[i + 2 * n] +
           0.3 * (runif01() - 0.5);

  SamplerOptions options;
  options.numTrees = 12;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(numBurnIn, 0, empty);

  FusedSuffstatCheck fused =
    sampler.chain(0).checkFusedSuffstatAgainstStockForTesting();
  char label[128];
  std::snprintf(label, sizeof(label),
                "n %% 4 == %lu: fused roll matches the stock roll bitwise",
                static_cast<unsigned long>(n % 4));
  check(fused.residualsAgreeBitwise, label);
  std::snprintf(label, sizeof(label),
                "n %% 4 == %lu: fused sumWeights matches the stock count",
                static_cast<unsigned long>(n % 4));
  check(fused.countsAgreeBitwise, label);
  std::snprintf(label, sizeof(label),
                "n %% 4 == %lu: fused suffstat matches the stock gather",
                static_cast<unsigned long>(n % 4));
  check(fused.worstRelativeError < 1e-9, label);
  std::snprintf(label, sizeof(label),
                "n = %lu: every tree took the fused pass",
                static_cast<unsigned long>(n));
  check(fused.numTreesFused == sampler.chain(0).numTrees(), label);
  printf("ok: fused suffstat at n = %lu (worst statistic gap %.3g)\n",
         static_cast<unsigned long>(n), fused.worstRelativeError);
}

static void testFusedSuffstatMatchesStock(ext_rng* rng) {
  for (size_t n = 1000; n <= 1003; ++n) checkFusedSuffstatShape(n, 20, rng);
  // n < 4: the prologue swallows the whole vector, so the pass is K = 1 and
  // must agree with the stock root kernel, which sums in the same obs order
  for (size_t n = 1; n <= 3; ++n) checkFusedSuffstatShape(n, 2, rng);
}

// Every eligibility clause is otherwise a silent decline, so each one is
// pinned positively through the fused run counter. The clauses are
// compile-time (leaf shape, ResidT) or O(1) (weights, leafOf freshness); the
// stale-map clause is exercised where that state already gets set up, in
// test_moves.cpp's testLeafOfConsistency.
static void testFusedSuffstatDeclines(ext_rng* rng) {
  const size_t n = 300, p = 3;
  std::vector<double> x(n * p), y(n), w(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n] + 0.2 * (runif01() - 0.5);
    w[i] = 0.5 + runif01();
  }

  SamplerOptions options;
  options.numTrees = 20;
  Results empty;

  // covered: the plain gaussian constant leaf, from the very first sweep -
  // there is no gate on a tree whose root is still its only bottom
  {
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                                ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, options, &rng);
    sampler.run(1, 0, empty);
    check(sampler.chain(0).fusedSuffstatRunsForTesting() == options.numTrees,
          "the fused suffstat runs on a fresh sampler's first sweep");
    size_t before = sampler.chain(0).fusedSuffstatRunsForTesting();
    sampler.run(1, 0, empty);
    check(sampler.chain(0).fusedSuffstatRunsForTesting() - before ==
            options.numTrees,
          "the fused suffstat runs for every tree of a burned-in sweep");
  }

  // declined: non-null weights, the one clause that also removes logistic,
  // negbin, t, the heteroscedastic mean forest, BCF and multinomial
  {
    ConstantLeafSampler sampler(x.data(), y.data(), n, p, w.data(), nullptr,
                                ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, options, &rng);
    sampler.run(2, 0, empty);
    check(sampler.chain(0).fusedSuffstatRunsForTesting() == 0,
          "weights decline the fused suffstat");
  }

  // declined: the opt-in fp32 residual, whose roll is left rolled on purpose
  {
    SamplerOptions fp32Options = options;
    fp32Options.fp32Residual = true;
    Sampler<ConstantGaussianLeaf, float> sampler(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, fp32Options, &rng);
    sampler.run(2, 0, empty);
    check(sampler.chain(0).fusedSuffstatRunsForTesting() == 0,
          "the fp32 residual declines the fused suffstat");
  }

  // declined: the vector leaf keeps no node averages at all, and the function
  // leaf keeps the dense slab in place of a leafOf map
  {
    std::vector<size_t> covariates = {2};
    SamplerOptions leafOptions = options;
    leafOptions.leafCovariateColumns = covariates.data();
    leafOptions.numLeafCovariates = 1;
    Sampler<LinearGaussianLeaf> linear(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, leafOptions, &rng);
    linear.run(2, 0, empty);
    check(linear.chain(0).fusedSuffstatRunsForTesting() == 0,
          "the vector leaf declines the fused suffstat");

    leafOptions.gpMaxLeafSize = 100;
    Sampler<GPGaussianLeaf> gp(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
      1.0, 3.0, 0.37804942330213542, leafOptions, &rng);
    gp.run(2, 0, empty);
    check(gp.chain(0).fusedSuffstatRunsForTesting() == 0,
          "the function leaf declines the fused suffstat");
  }

  // declined: both combiners always supply non-null weights, so the two
  // multi-forest models ride the weights clause with no rule of their own
  {
    std::vector<double> z(n), yBcf(n);
    for (size_t i = 0; i < n; ++i) {
      z[i] = runif01() < 0.5 ? 1.0 : 0.0;
      yBcf[i] = 2.0 * x[i] + z[i] * (1.0 + x[i + n]) + 0.2 * (runif01() - 0.5);
    }
    AmplitudeSpec spec;
    spec.mu.numTrees = 20; spec.mu.base = 0.95; spec.mu.power = 2.0;
    spec.tau.numTrees = 15; spec.tau.base = 0.95; spec.tau.power = 2.0;
    spec.z = z.data();
    Sampler<ConstantGaussianLeaf> bcf(
      x.data(), yBcf.data(), n, p, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, options, spec, &rng);
    bcf.run(2, 0, empty);
    check(bcf.chain(0).fusedSuffstatRunsForTesting() == 0,
          "BCF declines the fused suffstat through its combiner weights");
  }
  {
    const size_t K = 3;
    std::vector<int> counts(n * K, 0), trials(n, 1);
    for (size_t i = 0; i < n; ++i) {
      size_t label = static_cast<size_t>(runif01() * static_cast<double>(K));
      counts[(label < K ? label : K - 1) * n + i] = 1;
    }
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = counts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 15;
    Sampler<ConstantGaussianLeaf> multinomial(x.data(), n, p, options, spec,
                                              &rng);
    multinomial.run(2, 0, empty);
    check(multinomial.chain(0).fusedSuffstatRunsForTesting() == 0,
          "multinomial declines the fused suffstat through its omega weights");
  }

  printf("ok: fused suffstat decline set\n");
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
// drawn as the sum of w PG(1, psi) variates
// (docs/design/weighted-logistic.md),
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

  // (3) the counts ARE the PG shape, so replacing them redraws omega on the
  // spot rather than at the next sweep, which would move its trees first. The
  // swap runs the refresh setResponse already runs, from the same stream
  // position, so a swapped sampler is BITWISE one created with the new counts.
  std::vector<double> w2(n);
  for (size_t i = 0; i < n; ++i) w2[i] = 1.0 + (double) (i % 3);
  ext_rng* rngC = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngD = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngC == NULL || rngD == NULL || ext_rng_setSeed(rngC, 5150) != 0 ||
      ext_rng_setSeed(rngD, 5150) != 0) {
    check(false, "weighted logistic: swap rng creation");
    return;
  }
  ConstantLeafSampler built(x.data(), y.data(), n, p, w2.data(), nullptr,
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngC);
  ConstantLeafSampler swapped(x.data(), y.data(), n, p, w.data(), nullptr,
                   ResponseFamily::logistic, 1.0, 3.0, 1.0, options, &rngD);
  built.setResponse(y.data(), false);
  swapped.setWeights(w2.data());
  bool swapMatches = true, swapMoved = false;
  for (size_t i = 0; i < n; ++i) {
    swapMatches &= swapped.latents(0)[i] == built.latents(0)[i];
    if (swapped.latents(0)[i] != 0.25 * w2[i]) swapMoved = true;
  }
  check(swapMatches, "a logistic weight swap equals creation with the counts");
  check(swapMoved, "a logistic weight swap redraws omega off its cold start");
  ext_rng_destroy(rngD);
  ext_rng_destroy(rngC);

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

// The saved-tree replay's fan-out: the thread argument's effect on the
// partition, and the partition's effect on the answer (none). The traversal
// cutoff is overridden for the duration so a fixture this small still fans
// out; without that every count would collapse to one worker and an
// identity-across-counts check would prove nothing.
static void testPredictThreadPartition() {
  const size_t n = 200, nTest = 40, numChains = 2, numSamples = 11;
  const size_t numSlabs = numChains * numSamples;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> xTest(nTest * 2);
  std::uint_least32_t st = 20260824u;
  for (double& v : xTest) {
    st = st * 1664525u + 1013904223u;
    v = static_cast<double>(st >> 8) * (1.0 / 16777216.0);
  }

  std::vector<ext_rng*> rngs(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 4400 + static_cast<uint_least32_t>(c));
  }
  SamplerOptions options;
  options.numTrees = 20;
  options.numChains = numChains;
  options.numThreads = 1;
  options.keepTrees = true;
  options.numSamplesToStore = numSamples;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, rngs.data());
  Results results;
  sampler.run(30, numSamples, results);

  size_t savedCutoff = predictPartition.cutoffOverride;
  predictPartition.cutoffOverride = 1;

  const double sentinel = -9876.5;
  std::vector<double> serial(nTest * numSlabs, sentinel);
  sampler.predict(xTest.data(), nTest, 1, serial.data());
  bool serialWritten = true;
  for (double value : serial) serialWritten &= value != sentinel;
  check(serialWritten && predictPartition.numWorkers == 1,
        "a one-worker replay writes every slab");

  // 3 and 7 divide neither the slab count nor each other, so an
  // even-division assumption in the partition shows up as a gap
  const size_t counts[] = {1, 2, 3, 7, 64};
  bool identical = true, resolvedRight = true, mapRight = true,
       written = true;
  for (size_t requested : counts) {
    std::vector<double> threaded(nTest * numSlabs, sentinel);
    sampler.predict(xTest.data(), nTest, requested, threaded.data());
    identical &= threaded == serial;
    for (double value : threaded) written &= value != sentinel;

    size_t expectedWorkers = requested < numSlabs ? requested : numSlabs;
    resolvedRight &= predictPartition.resolvedThreads == requested;
    resolvedRight &= predictPartition.numWorkers == expectedWorkers;

    // the map must be a contiguous block partition covering every slab: non-
    // decreasing, starting at worker 0, ending at the last worker, with block
    // sizes differing by at most one
    const std::vector<size_t>& map = predictPartition.workerForSlab;
    mapRight &= map.size() == numSlabs;
    if (map.size() != numSlabs) continue;
    std::vector<size_t> owned(expectedWorkers, 0);
    for (size_t slab = 0; slab < numSlabs; ++slab) {
      mapRight &= map[slab] < expectedWorkers;
      if (slab > 0) mapRight &= map[slab] >= map[slab - 1];
      if (map[slab] < expectedWorkers) ++owned[map[slab]];
    }
    size_t smallest = numSlabs, largest = 0;
    for (size_t count : owned) {
      if (count < smallest) smallest = count;
      if (count > largest) largest = count;
    }
    mapRight &= smallest > 0 && largest - smallest <= 1;
    mapRight &= map[0] == 0 && map[numSlabs - 1] == expectedWorkers - 1;
  }
  check(identical && written,
        "the replay is bitwise identical at every worker count");
  check(resolvedRight, "the requested count resolves to that many workers");
  check(mapRight, "every slab is owned by exactly one worker, in blocks");

  // 0 means the sampler's own count
  sampler.setNumThreads(5);
  std::vector<double> viaSampler(nTest * numSlabs, sentinel);
  sampler.predict(xTest.data(), nTest, 0, viaSampler.data());
  check(predictPartition.resolvedThreads == 5 && viaSampler == serial,
        "a zero argument takes the sampler's own count");

  // a sampler whose own count was stored as 0 must still write its output:
  // the floor is what keeps the worker count off zero
  sampler.setNumThreads(0);
  std::vector<double> viaZero(nTest * numSlabs, sentinel);
  sampler.predict(xTest.data(), nTest, 0, viaZero.data());
  check(predictPartition.resolvedThreads == 1 &&
          predictPartition.numWorkers == 1 && viaZero == serial,
        "a sampler count of zero floors to one worker, not to none");

  // with the derived cutoff back in force this fixture is far below it, so
  // the replay stays inline however many threads are asked for
  predictPartition.cutoffOverride = savedCutoff;
  sampler.setNumThreads(1);
  std::vector<double> belowCutoff(nTest * numSlabs, sentinel);
  sampler.predict(xTest.data(), nTest, 8, belowCutoff.data());
  check(predictPartition.resolvedThreads == 8 &&
          predictPartition.numWorkers == 1 && belowCutoff == serial,
        "a replay below the traversal cutoff runs inline");

  for (size_t c = numChains; c > 0; --c) ext_rng_destroy(rngs[c - 1]);
  printf("ok: predict thread partition (%lu slabs)\n",
         static_cast<unsigned long>(numSlabs));
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

// E1 (docs/plans/variance-forest-mutation-routing.md, slice S4): a variance
// forest pins seven n-sized allocations at the creation count, so a
// replacement data set of a different length used to overrun them. setData
// itself returned cleanly and the fault landed on the FOLLOWING run, which
// segfaulted in 4 tries out of 5 and reported a non-finite variance in the
// fifth - which is why this probe belongs in the ASAN leg, where the
// out-of-bounds write is seen on the try that would have survived.
static void testSetDataVarianceForestResize(ext_rng* rng) {
  const size_t n = 200, nGrown = 5000;
  std::vector<double> x, y, xGrown, yGrown;
  makeMutationData(x, y, n);
  makeMutationData(xGrown, yGrown, nGrown);
  SamplerOptions options;
  options.numTrees = 25;
  options.numVarianceTrees = 10;
  ConstantLeafSampler sampler(x.data(), y.data(), n, size_t(2), nullptr,
                              nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  Results empty;
  sampler.run(50, 0, empty);

  auto runAndCheck = [&](size_t count, const char* label) {
    std::vector<double> varianceDraws(count * 5);
    Results results;
    results.varianceFits = varianceDraws.data();
    sampler.run(0, 5, results);
    bool finite = true;
    for (double v : varianceDraws) finite &= std::isfinite(v) && v > 0.0;
    check(finite, label);
  };
  sampler.setData(xGrown.data(), yGrown.data(), nGrown, nullptr, nullptr,
                  nullptr, 0);
  check(sampler.numObservations() == nGrown,
        "setData grows a heteroscedastic sampler");
  runAndCheck(nGrown, "s^2 stays finite and positive after setData grows n");
  sampler.setData(x.data(), y.data(), n, nullptr, nullptr, nullptr, 0);
  check(sampler.numObservations() == n,
        "setData shrinks a heteroscedastic sampler");
  runAndCheck(n, "s^2 stays finite and positive after setData shrinks n");

  printf("ok: setData resize under a variance forest\n");
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
  options.predictors.columnTypes = nullptr;  // set via ctor path below
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  options.predictors.columnTypes = types;
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
  options.predictors.columnTypes = types;
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
  options.predictors.columnTypes = types;
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
  sampler.predict(xTest.data(), nTest, 1, predicted.data());
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
    live.predict(x.data(), n, 1, livePredict.data());
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

// The active-row channel at the sampler surface, where the engine's own
// normalizer and value scan live: an all-ones mask installs NOTHING and is
// bitwise the unmasked run; a fractional or NaN element refuses the whole call
// and perturbs nothing; an all-zeros mask is ACCEPTED and runs with every
// forest at its prior; a family v1 does not build refuses. Closes with the
// vector-leaf arm: installing a mask mid-run must drop the linear leaf's U'WU
// cache, which re-validates on the ordered member list alone and so cannot see
// a weight change that moves no membership.
static void testActiveRows() {
  const size_t n = 200, numSamples = 8;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> weights(n), active(n), composed(n), ones(n, 1.0),
    zeros(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    weights[i] = 0.5 + runif01();
    active[i] = i % 4 == 0 ? 0.0 : 1.0;
    composed[i] = weights[i] * active[i];
  }
  // the bad vectors are otherwise-valid masks, so a partial application would
  // install a real mask and show up in the draws
  std::vector<double> fractional(active), missing(active);
  fractional[n - 3] = 0.5;
  missing[n - 3] = std::nan("");

  SamplerOptions options;
  options.numTrees = 25;
  auto makeSeededRng = []() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 20260816u);
    return rng;
  };
  auto makeSampler = [&](ext_rng*& rng, ResponseFamily family,
                         const double* response, const double* w) {
    rng = makeSeededRng();
    return std::make_unique<ConstantLeafSampler>(
      x.data(), response, n, 2, w, nullptr, family, 1.0, 3.0,
      0.37804942330213542, options, &rng);
  };

  // ---- the normalizer and the refusals, all against one unmasked twin ----
  std::vector<double> sigmaPlain(numSamples), trainPlain(n * numSamples);
  std::vector<double> sigmaOther(numSamples), trainOther(n * numSamples);
  Results resultsPlain, resultsOther;
  resultsPlain.sigma = sigmaPlain.data();
  resultsPlain.trainingFits = trainPlain.data();
  resultsOther.sigma = sigmaOther.data();
  resultsOther.trainingFits = trainOther.data();

  ext_rng* rngPlain;
  auto plain = makeSampler(rngPlain, ResponseFamily::gaussian, y.data(),
                           weights.data());
  plain->run(20, numSamples, resultsPlain);

  struct Degenerate { const double* values; const char* what; bool accepted; };
  const Degenerate degenerates[] = {
    {ones.data(), "an all-ones mask installs nothing", true},
    {fractional.data(), "a fractional element refuses the whole call", false},
    {missing.data(), "an NA element refuses the whole call", false}};
  for (const Degenerate& d : degenerates) {
    ext_rng* rng;
    auto sampler = makeSampler(rng, ResponseFamily::gaussian, y.data(),
                               weights.data());
    // the engine's two false reasons stay SEPARABLE by the probe: a chain
    // that advertises the channel still refuses a non-binary mask, which is
    // what lets a surface answer the capability with its return value while
    // raising on the value
    check(sampler->supportsActiveRows(),
          "the mask channel is advertised whatever the values say");
    check(sampler->setActiveRows(d.values) == d.accepted, d.what);
    sampler->run(20, numSamples, resultsOther);
    check(sigmaPlain == sigmaOther && trainPlain == trainOther,
          "a normalized or refused mask leaves the draws bitwise unchanged");
    ext_rng_destroy(rng);
  }

  // ---- all zeros: accepted, finite, forest at its prior, fits reported ----
  {
    ext_rng* rng;
    auto sampler = makeSampler(rng, ResponseFamily::gaussian, y.data(),
                               weights.data());
    check(sampler->setActiveRows(zeros.data()),
          "an all-zeros mask is accepted, not refused");
    check(sampler->chain(0).sigmaDegreesOfFreedomForTesting() == 3.0,
          "an all-zeros mask leaves the sigma posterior at the prior df");
    sampler->run(20, numSamples, resultsOther);
    bool finite = true;
    for (size_t s = 0; s < numSamples; ++s) finite = finite && sigmaOther[s] > 0.0;
    for (size_t i = 0; i < n * numSamples; ++i)
      finite = finite && std::isfinite(trainOther[i]);
    check(finite, "an all-zeros mask runs finite and still reports every fit");
    ext_rng_destroy(rng);
  }

  // ---- every family accepts, multinomial through its coupling ----
  // The multinomial response holds no precisions to compose a mask into, so it
  // advertises the channel and the combiner carries it (the K interleaved
  // Polya-Gamma precisions); the kernel arm is testActiveRowsMultinomialKernel.
  {
    std::vector<double> binary(n), counts(n);
    for (size_t i = 0; i < n; ++i) {
      binary[i] = y[i] > 0.0 ? 1.0 : 0.0;
      counts[i] = static_cast<double>(i % 5);
    }
    struct Reachable { ResponseFamily family; const double* response; };
    const Reachable reachable[] = {{ResponseFamily::probit, binary.data()},
                                   {ResponseFamily::logistic, binary.data()},
                                   {ResponseFamily::nbinom, counts.data()}};
    for (const Reachable& r : reachable) {
      ext_rng* rng;
      auto sampler = makeSampler(rng, r.family, r.response, nullptr);
      check(sampler->supportsActiveRows() &&
              sampler->setActiveRows(active.data()),
            "a latent-family sampler accepts a mask though it refuses weights");
      ext_rng_destroy(rng);
    }

    // an all-zeros mask on logistic and nbinom also runs finite, exactly as
    // the gaussian arm above
    for (const Reachable& r : reachable) {
      if (r.family == ResponseFamily::probit) continue;
      ext_rng* rng;
      auto sampler = makeSampler(rng, r.family, r.response, nullptr);
      check(sampler->setActiveRows(zeros.data()),
            "an all-zeros mask is accepted on this family too");
      std::vector<double> sigmaZ(numSamples), trainZ(n * numSamples);
      Results resultsZ;
      resultsZ.sigma = sigmaZ.data();
      resultsZ.trainingFits = trainZ.data();
      sampler->run(10, numSamples, resultsZ);
      bool finiteZ = true;
      for (size_t s = 0; s < numSamples; ++s)
        finiteZ = finiteZ && sigmaZ[s] > 0.0;
      for (double v : trainZ) finiteZ = finiteZ && std::isfinite(v);
      check(finiteZ, "an all-zeros mask runs finite and reports every fit");
      ext_rng_destroy(rng);
    }

    // aft needs its own construction (survivalStatus), so it is not in
    // reachable above; same all-zeros claim
    SamplerOptions aftOptions = options;
    std::vector<double> statusAll(n, 1.0);  // no censoring; any finite y works
    aftOptions.survivalStatus = statusAll.data();
    ext_rng* rngAft = makeSeededRng();
    ConstantLeafSampler aftSampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                                   ResponseFamily::aft, 1.0, 3.0,
                                   0.37804942330213542, aftOptions, &rngAft);
    check(aftSampler.setActiveRows(zeros.data()),
          "an all-zeros mask is accepted on an aft sampler");
    std::vector<double> sigmaAftZ(numSamples), trainAftZ(n * numSamples);
    Results resultsAftZ;
    resultsAftZ.sigma = sigmaAftZ.data();
    resultsAftZ.trainingFits = trainAftZ.data();
    aftSampler.run(10, numSamples, resultsAftZ);
    bool finiteAftZ = true;
    for (size_t s = 0; s < numSamples; ++s)
      finiteAftZ = finiteAftZ && sigmaAftZ[s] > 0.0;
    for (double v : trainAftZ) finiteAftZ = finiteAftZ && std::isfinite(v);
    check(finiteAftZ, "an all-zeros mask runs finite on an aft sampler too");
    ext_rng_destroy(rngAft);

    const size_t K = 3;
    std::vector<int> categoryCounts(n * K, 0), trials(n, 1);
    for (size_t i = 0; i < n; ++i) categoryCounts[(i % K) * n + i] = 1;
    MultinomialSpec spec;
    spec.numCategories = K;
    spec.counts = categoryCounts.data();
    spec.trials = trials.data();
    spec.forest.numTrees = 10;
    ext_rng* rngMultinomial = makeSeededRng();
    ConstantLeafSampler multinomial(x.data(), n, 2, options, spec,
                                    &rngMultinomial);
    check(multinomial.supportsActiveRows() &&
            multinomial.setActiveRows(active.data()) &&
            !multinomial.setActiveRows(fractional.data()),
          "a multinomial sampler takes the global mask and refuses a "
          "fractional one");
    // all zeros on a K-forest sampler: every category's forest sits at its
    // prior and the reported softmax probabilities stay a per-row simplex
    check(multinomial.setActiveRows(zeros.data()),
          "an all-zeros mask is accepted on a multinomial sampler");
    std::vector<double> softmax(n * K * numSamples);
    Results resultsMultinomial;
    resultsMultinomial.trainingFits = softmax.data();
    resultsMultinomial.numReportedLocations = K;
    multinomial.run(10, numSamples, resultsMultinomial);
    bool simplex = true;
    for (size_t s = 0; s < numSamples; ++s)
      for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t k = 0; k < K; ++k) sum += softmax[i + n * (k + K * s)];
        simplex = simplex && std::fabs(sum - 1.0) < 1e-12;
      }
    check(simplex,
          "an all-zeros multinomial mask runs and still reports every row's "
          "softmax");
    ext_rng_destroy(rngMultinomial);
  }
  ext_rng_destroy(rngPlain);

  // ---- grouped intercepts delegate, with no edit of their own ----
  // drawGroupEffects already weights its per-group sums by workingWeights(),
  // so an inactive row leaves its group's mean and precision; a group whose
  // every row is inactive falls back to its prior through the same formula,
  // which is coherent and is NOT what deleting the group would do. Pinned
  // BITWISE against setWeights(w * a) rather than for finiteness: a mask that
  // never reached the per-group sums would still run finite.
  {
    std::vector<std::uint32_t> groups(n);
    for (size_t i = 0; i < n; ++i)
      groups[i] = static_cast<std::uint32_t>(i % 5);
    std::vector<double> groupMask(n), groupComposed(n);
    for (size_t i = 0; i < n; ++i) {
      // group 0 entirely inactive, plus a scatter of rows elsewhere
      groupMask[i] = (groups[i] == 0 || i % 7 == 3) ? 0.0 : 1.0;
      groupComposed[i] = weights[i] * groupMask[i];
    }

    SamplerOptions groupedOptions = options;
    groupedOptions.groupIndices = groups.data();
    groupedOptions.numGroups = 5;
    ext_rng* rngMasked = makeSeededRng();
    ext_rng* rngComposed = makeSeededRng();
    ConstantLeafSampler grouped(x.data(), y.data(), n, 2, weights.data(),
                                nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                                0.37804942330213542, groupedOptions,
                                &rngMasked);
    ConstantLeafSampler composedGrouped(
      x.data(), y.data(), n, 2, groupComposed.data(), nullptr,
      ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, groupedOptions,
      &rngComposed);
    check(grouped.supportsActiveRows() &&
            grouped.setActiveRows(groupMask.data()),
          "a grouped sampler forwards the mask to its base family");

    std::vector<double> effects(5 * numSamples),
      effectsComposed(5 * numSamples), trainComposed(n * numSamples),
      sigmaComposed(numSamples);
    Results groupedResults(resultsOther), composedResults;
    groupedResults.groupEffects = effects.data();
    composedResults.sigma = sigmaComposed.data();
    composedResults.trainingFits = trainComposed.data();
    composedResults.groupEffects = effectsComposed.data();
    grouped.run(20, numSamples, groupedResults);
    composedGrouped.run(20, numSamples, composedResults);
    bool finite = true;
    for (double effect : effects) finite = finite && std::isfinite(effect);
    check(finite,
          "an entirely inactive group draws its effect from the prior, finite");
    check(effects == effectsComposed && trainOther == trainComposed &&
            sigmaOther == sigmaComposed,
          "a masked grouped sampler is bitwise setWeights(w * a), effects and "
          "all");
    ext_rng_destroy(rngComposed);
    ext_rng_destroy(rngMasked);
  }

  // ---- the vector leaf: a mid-run install must drop the U'WU cache ----
  // Both arms move the same weights by the same values and both invalidate,
  // so their draws agree bitwise; the mask install moves no MEMBERSHIP, so a
  // cache that was not dropped would serve stale statistics on the leaves at
  // or over minCachedLeafSize and this comparison would part.
  {
    const size_t columns[] = {1};
    SamplerOptions leafOptions = options;
    leafOptions.numTrees = 5;  // few trees, so the leaves clear the cache floor
    leafOptions.leafCovariateColumns = columns;
    leafOptions.numLeafCovariates = 1;

    ext_rng* rngMasked = makeSeededRng();
    ext_rng* rngWeighted = makeSeededRng();
    Sampler<LinearGaussianLeaf> masked(x.data(), y.data(), n, 2, weights.data(),
                                       nullptr, ResponseFamily::gaussian, 1.0,
                                       3.0, 0.37804942330213542, leafOptions,
                                       &rngMasked);
    Sampler<LinearGaussianLeaf> weighted(x.data(), y.data(), n, 2,
                                         weights.data(), nullptr,
                                         ResponseFamily::gaussian, 1.0, 3.0,
                                         0.37804942330213542, leafOptions,
                                         &rngWeighted);
    std::vector<double> sigmaWeighted(numSamples),
      trainWeighted(n * numSamples);
    Results resultsWeighted;
    resultsWeighted.sigma = sigmaWeighted.data();
    resultsWeighted.trainingFits = trainWeighted.data();
    masked.run(20, numSamples, resultsOther);  // warm the cache on both arms
    weighted.run(20, numSamples, resultsWeighted);
    check(masked.setActiveRows(active.data()),
          "a linear-leaf sampler accepts an active-row mask");
    weighted.setWeights(composed.data());
    masked.run(0, numSamples, resultsOther);
    weighted.run(0, numSamples, resultsWeighted);
    check(sigmaOther == sigmaWeighted && trainOther == trainWeighted,
          "a mid-run mask install on a vector leaf equals setWeights(w * a)");
    ext_rng_destroy(rngWeighted);
    ext_rng_destroy(rngMasked);
  }

  printf("ok: active rows at the sampler surface\n");
}

// The mask installed on a GROWN forest, which is where the veto's weight law
// has no gate to hold it: masking rows a tree already split on leaves leaves
// no likelihood term reaches. The forest must keep moving there (the promise
// setActiveRows makes: it sits at its PRIOR, which is a distribution, not a
// freeze), must never install a member-empty leaf while it does - the state
// law export and restore require - and must recover an admissible forest once
// the mask lifts.
static void testActiveRowsOnGrownForest() {
  const size_t n = 300, numSamples = 4;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  std::vector<double> zeros(n, 0.0), partial(n, 1.0);
  for (size_t i = 0; i < n; ++i)
    if (x[i] < 0.5) partial[i] = 0.0;  // a half-space of x0, so leaves fall in

  SamplerOptions options;
  options.numTrees = 25;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20260818u);
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  std::vector<double> sigma(numSamples), train(n * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.trainingFits = train.data();
  sampler.run(100, 0, results);  // grow the forest under no mask at all

  auto forestSignature = [&]() {
    std::uint64_t hash = 0;
    for (size_t t = 0; t < options.numTrees; ++t)
      hash ^= treeStructureSignature(sampler.chain(0).tree(t)) +
              static_cast<std::uint64_t>(t);
    return hash;
  };
  auto countVetoed = [&](const double* weights) {
    size_t vetoed = 0;
    for (size_t t = 0; t < options.numTrees; ++t) {
      const Tree& tree(sampler.chain(0).tree(t));
      std::vector<std::int32_t> bottoms;
      tree.fillBottom(0, bottoms);
      for (std::int32_t b : bottoms)
        vetoed += tree.leafVetoRank(b, weights) != 0 ? 1 : 0;
    }
    return vetoed;
  };
  auto occupied = [&]() {
    bool all = true;
    for (size_t t = 0; t < options.numTrees; ++t)
      all = all && sampler.chain(0).tree(t).bottomNodesAreOccupied();
    return all;
  };

  // ---- the whole forest stranded ----
  check(sampler.setActiveRows(zeros.data()), "an all-zeros mask installs");
  size_t strandedLeaves = countVetoed(zeros.data());
  std::uint64_t before = forestSignature();
  sampler.run(100, numSamples, results);
  check(forestSignature() != before,
        "a fully masked grown forest keeps moving, rather than freezing");
  check(occupied(),
        "no move installs a member-empty leaf while the forest is masked");
  bool finite = true;
  for (double v : train) finite = finite && std::isfinite(v);
  for (double v : sigma) finite = finite && v > 0.0;
  check(finite, "a fully masked grown forest reports finite draws throughout");

  // the state law is what export and restore gate on, and it still holds
  SamplerStateData state;
  sampler.getState(state);
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng2, 4242);
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                               ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, options, &rng2);
  check(restored.setState(state, nullptr),
        "a state drawn under a full mask restores");
  checkStructuralRoundTrip(state, restored,
                           "the masked state reproduces the model");
  ext_rng_destroy(rng2);

  // ---- the mask lifts: every leaf must be admissible again ----
  check(sampler.setActiveRows(nullptr), "the mask clears");
  sampler.run(0, numSamples, results);
  check(countVetoed(nullptr) == 0,
        "the forest is admissible again once the mask lifts");

  // ---- a PARTIAL mask on the grown forest is absorbed back into the set ----
  check(sampler.setActiveRows(partial.data()), "a partial mask installs");
  size_t strandedPartial = countVetoed(partial.data());
  check(strandedPartial > 0,
        "non-vacuity: the partial mask strands leaves of the grown forest");
  sampler.run(200, 0, results);
  check(countVetoed(partial.data()) == 0,
        "the partially masked forest is absorbed back into the admissible set");
  check(occupied(), "and installs no member-empty leaf getting there");

  ext_rng_destroy(rng);
  printf("ok: active rows on a grown forest (%zu leaves stranded by the full "
         "mask, %zu by the partial one)\n",
         strandedLeaves, strandedPartial);
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
  samplerB.predict(xTest.data(), nTest, 1, predicted.data());
  check(predicted == testFits,
        "post-toggle saved predictions equal the recorded test fits");

  samplerB.setTreeStorage(true, 4);
  check(samplerB.savedTreeCapacity() == 4,
        "an unchanged toggle preserves the storage");
  samplerB.setTreeStorage(false, 0);
  check(samplerB.savedTreeCapacity() == 0, "disabling storage frees it");
  std::vector<double> livePredictions(nTest);
  samplerB.predict(xTest.data(), nTest, 1, livePredictions.data());
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
  sampler.predict(xTest, 2, 1, predictions.data());
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
  AmplitudeSpec spec;
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

  // the all-control basis, ROW-major (1, 0) per row: setForestBasis is the sole
  // basis-mutation route, so what was a treatment swap is an install
  std::vector<double> controlBasis(2 * n, 0.0);
  for (size_t i = 0; i < n; ++i) controlBasis[2 * i] = 1.0;
  check(sampler.setForestBasis(1, controlBasis.data(), 2),
        "the sole basis-mutation route installs on the treatment forest");
  std::vector<double> sigma2(10), fits2(n * 10);
  Results r2;
  r2.sigma = sigma2.data();
  r2.trainingFits = fits2.data();
  sampler.run(0, 10, r2);
  bool refreshed = true;
  for (size_t i = 0; i < n * 10 && refreshed; ++i)
    refreshed = std::isfinite(fits2[i]);
  check(refreshed, "BCF setForestBasis refresh runs");

  printf("ok: BCF two-forest sampler\n");
}

// A BCF sampler admits a whole-response swap at updateScale = false: the
// gaussian response re-maps y through the pinned (min_, range_) and touches no
// forest, and the combiner re-derives every per-forest residual from y each
// sweep, so nothing per-forest is stale afterwards. Two arms - the pins hold
// and both forests stay self-consistent, and the swap is bitwise the same
// chain as the already-permitted setOffset(yBuild - yNew), which re-maps to
// the same working response through the same transform.
static void testBCFResponseSwap() {
  const size_t n = 300, p = 3, muTrees = 30, tauTrees = 15;
  std::vector<double> x(n * p), z(n), yBuild(n), yNew(n), delta(n);
  // A local stream, so adding this test does not shift the shared runif01()
  // state the hardcoded characteristic values downstream depend on.
  std::uint64_t state = 20260805u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    yBuild[i] = mu + z[i] * tau + 0.2 * (unif() - 0.5);
    // The offset arm reaches its working response by subtraction, so yNew must
    // be exactly what that subtraction produces; without this round trip the
    // two arms differ in the last ulp and the comparison measures fp rounding.
    delta[i] = yBuild[i] - (0.5 * mu - z[i] * tau + 0.3 * (unif() - 0.5));
    yNew[i] = yBuild[i] - delta[i];
  }

  SamplerOptions options;
  AmplitudeSpec spec;
  spec.mu.numTrees = muTrees; spec.mu.base = 0.95; spec.mu.power = 2.0;
  spec.tau.numTrees = tauTrees; spec.tau.base = 0.25; spec.tau.power = 3.0;
  spec.z = z.data();

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngA, 20260805u);
  ext_rng_setSeed(rngB, 20260805u);
  Sampler<ConstantGaussianLeaf> samplerA(
    x.data(), yBuild.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rngA);
  Sampler<ConstantGaussianLeaf> samplerB(
    x.data(), yBuild.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rngB);
  check(samplerA.supportsResponseMutation(),
        "a BCF sampler opts into multi-forest response mutation");

  const size_t numWarm = 60, numAfter = 40;
  std::vector<double> warmSigma(numWarm), warmFits(n * numWarm);
  Results warm;
  warm.sigma = warmSigma.data();
  warm.trainingFits = warmFits.data();
  samplerA.run(0, numWarm, warm);
  samplerB.run(0, numWarm, warm);

  double fitScale = samplerA.fitScale(),
         leafScale = samplerA.chain(0).leaf().scale, sigma = samplerA.sigma(0);
  samplerA.setResponse(yNew.data(), false);
  samplerB.setOffset(delta.data(), false);
  check(samplerA.fitScale() == fitScale &&
          samplerA.chain(0).leaf().scale == leafScale,
        "BCF response swap pins the response transform and the leaf scale");
  check(samplerA.sigma(0) == sigma,
        "BCF response swap leaves sigma continuous");

  std::vector<double> sigmaA(numAfter), fitsA(n * numAfter), sigmaB(numAfter),
    fitsB(n * numAfter);
  Results rA, rB;
  rA.sigma = sigmaA.data();
  rA.trainingFits = fitsA.data();
  rB.sigma = sigmaB.data();
  rB.trainingFits = fitsB.data();
  samplerA.run(0, numAfter, rA);
  samplerB.run(0, numAfter, rB);

  // Every forest's cached total still sums its own trees' fits, and both
  // forests are still moving off zero. The tolerance is relative rather than
  // ulp-level because finalizeTotalFits reaches the total by cancelling the
  // running residual against the forest response, which BCF divides by b_z:
  // with b0 near its 1e-9 multiplier floor that division carries ~1e-8 of
  // absolute slop into the tau total (present in an unswapped chain too, and
  // unchanged by the swap). A forest left stale would be off by O(1).
  std::vector<double> totals(n), treeFits(n * muTrees), totalsB(n);
  bool consistent = true, moves = true;
  for (size_t f = 0; f < 2; ++f) {
    size_t numTrees = f == 0 ? muTrees : tauTrees;
    samplerA.chain(0).forestTotalFits(f, totals.data());
    samplerA.chain(0).forestTreeFits(f, treeFits.data());
    double worst = 0.0, largest = 0.0, ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double sum = 0.0;
      for (size_t t = 0; t < numTrees; ++t) sum += treeFits[i + t * n];
      worst = std::max(worst, std::fabs(totals[i] - sum));
      largest = std::max(largest, std::fabs(totals[i]));
      ss += totals[i] * totals[i];
    }
    consistent &= worst <= 1.0e-5 * std::max(1.0, largest);
    moves &= ss > 0.0;
  }
  check(consistent, "BCF forest totals still sum their trees after the swap");
  check(moves, "both BCF forests still move after the response swap");

  // The strongest assertion: the swap is the already-permitted offset re-map
  // with a different pointer. Draws, both forests and the glue must agree to
  // the bit. Reported fits are excluded on purpose - storeSample adds the
  // offset back, so the offset arm reports fit + delta by design.
  bool identical = true;
  for (size_t s = 0; s < numAfter; ++s) identical &= sigmaA[s] == sigmaB[s];
  for (size_t f = 0; f < 2; ++f) {
    samplerA.chain(0).forestTotalFits(f, totals.data());
    samplerB.chain(0).forestTotalFits(f, totalsB.data());
    for (size_t i = 0; i < n; ++i) identical &= totals[i] == totalsB[i];
  }
  double a, b0, b1, a2, b02, b12;
  samplerA.chain(0).bcfGlue(a, b0, b1);
  samplerB.chain(0).bcfGlue(a2, b02, b12);
  identical &= a == a2 && b0 == b02 && b1 == b12;
  check(identical, "BCF setResponse(yNew, false) is bitwise the "
                   "setOffset(yBuild - yNew) chain");

  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
  printf("ok: BCF scale-pinned response swap\n");
}

// Regression (heap-use-after-free): a BCF mu-forest interaction constraint must
// not dangle. mu's trees borrow the address of forest.interaction; building the
// tau forest reallocates forests_ and relocates mu, so a value-member constraint
// would leave every mu tree's interaction_ pointing at freed memory - a UAF that
// only ASAN (and non-macOS allocators) expose, which shipped past the single-
// forest tests/cpp gate. The heap-allocated constraint keeps its address across
// the move. Build with a mu order cap, run, and require clean finite fits; under
// -fsanitize=address this aborts on the old bug.
static void testBCFInteractionLifetime() {
  uint64_t savedRngState = rngState;  // RNG-neutral: no downstream snapshot shift
  const size_t n = 300, p = 4;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    z[i] = runif01() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (runif01() - 0.5);
  }

  SamplerOptions options;
  AmplitudeSpec spec;
  spec.mu.numTrees = 40; spec.mu.base = 0.95; spec.mu.power = 2.0;
  spec.mu.interactionMaxOrder = 1;  // constrains the FIRST-built (mu) forest
  spec.tau.numTrees = 20; spec.tau.base = 0.25; spec.tau.power = 3.0;
  spec.z = z.data();

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 90210u);
  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &localRng);

  const size_t numBurnIn = 100, numSamples = 100;
  std::vector<double> sigma(numSamples), fits(n * numSamples);
  Results results;
  results.sigma = sigma.data();
  results.trainingFits = fits.data();
  sampler.run(numBurnIn, numSamples, results);

  std::vector<double> muFits(n);
  sampler.forestTotalFits(0, 0, muFits.data());
  bool clean = true;
  double muSS = 0.0;
  for (size_t i = 0; i < n; ++i) {
    clean &= std::isfinite(muFits[i]);
    muSS += muFits[i] * muFits[i];
  }
  check(clean && muSS > 0.0,
        "BCF mu-interaction constraint survives the tau-forest realloc");

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: BCF mu-interaction constraint lifetime (tau-forest realloc)\n");
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
  AmplitudeSpec spec;
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
  AmplitudeSpec spec;
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
  AmplitudeSpec spec;
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
  sampler.predict(x.data(), n, 1, pred.data());  // scale * mu_saved + shift

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
    AmplitudeSpec specThin;
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
    samplerThin.predict(xThin.data(), nThin, 1, predThin.data());

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
// interweave); the bridge's growFromRoot entry point places no multi-forest
// guard, so this path is R-reachable, not engine-internal only. The fuzz
// mutation suite now walks it too, but only against consistency oracles -
// this is the only gate that pins a hardcoded characteristic value: build the
// two-forest sampler, grow every forest from the root, and pin the combined
// output - both forests finite and off zero, glue finite, and a recorded
// characteristic value of the combined internal fit a*mu + b_z*tau.
static void testBCFGrowForestFromRoot() {
  // A local stream and locally owned generator, so this test neither reads nor
  // shifts the shared runif01() state: its fixture, and so its hardcoded
  // characteristic value, is the same whether or not earlier suites ran.
  std::uint64_t state = 20260819u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  const size_t n = 400, p = 3;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (unif() - 0.5);
  }

  SamplerOptions options;
  AmplitudeSpec spec;
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
  // two-forest state. Deterministic given the local fixture stream and
  // localRng seed 90210, in a filtered run as in a full one; a relocation
  // that shifts the grow sweep's draw order or the coupling moves it far past
  // the tolerance, while it survives benign cross-build FP reassociation.
  double combinedMean = 0.0;
  for (size_t i = 0; i < n; ++i)
    combinedMean += a * muFits[i] + (z[i] != 0.0 ? b1 : b0) * tauFits[i];
  combinedMean /= static_cast<double>(n);
  checkNear(combinedMean, -0.028618738206336595, 1e-6,
            "BCF grow-from-root combined fit characteristic value");

  ext_rng_destroy(localRng);
  printf("ok: BCF grow-from-root sweep\n");
}

// Per-forest replay: predictPerForest at the TRAINING rows reproduces each
// forest's own resident totalFits (forestTotalFits), forest by forest. That is
// the whole contract - RAW, forest-major, no fitScale, no fitShift, no offset -
// and it is stated as a tolerance rather than an equality on purpose: the
// amplitude ridge rescales totalFits as c * sum(mu_t) while it rescales the
// saved leaves themselves, so the replay re-sums sum(c * mu_t) and the two
// associate differently (combiner.hpp's rescale move).
//
// Both replay routes are gated: the saved-slot one under keepTrees, whose LAST
// slot is the run's final recorded sweep and so the live position, and the
// live-tree one without it. A local fixture stream plus a locally owned
// generator keep this neutral to the shared runif01() state.
static void testAmplitudePerForestReplay() {
  std::uint64_t state = 20260820u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  const size_t n = 250, p = 3, numSamples = 4;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (unif() - 0.5);
  }

  AmplitudeSpec spec;
  spec.mu.numTrees = 30;
  spec.tau.numTrees = 15;  // ragged tree counts: each forest drives its own loop
  spec.z = z.data();

  // the two forests' resident totals, and the replayed channels to match them
  std::vector<double> live(n), replayed(n * 2 * numSamples);

  auto maxPerForestGap = [&](Sampler<ConstantGaussianLeaf>& sampler,
                             const double* slab) {
    double worst = 0.0;
    for (size_t f = 0; f < 2; ++f) {
      sampler.forestTotalFits(0, f, live.data());
      for (size_t i = 0; i < n; ++i)
        worst = std::max(worst, std::fabs(slab[i + f * n] - live[i]));
    }
    return worst;
  };

  {
    SamplerOptions options;
    options.keepTrees = true;
    options.numSamplesToStore = numSamples;
    ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(localRng, 5150u);
    Sampler<ConstantGaussianLeaf> sampler(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, options, spec, &localRng);

    Results results;
    sampler.run(120, numSamples, results);
    sampler.predictPerForest(x.data(), n, 1, replayed.data());

    const double* lastSlot = replayed.data() + (numSamples - 1) * n * 2;
    check(maxPerForestGap(sampler, lastSlot) < 1.0e-12,
          "saved-slot per-forest replay reproduces each forest's own total");

    // the RAW convention, pinned against the one channel that does carry the
    // transform: predict replays forest 0 as fitScale * mu + fitShift, so its
    // difference from the raw forest-0 replay is a constant (the shift) for
    // every row. A scale applied inside the per-forest entry would make that
    // difference row-dependent.
    std::vector<double> combined(n * numSamples);
    sampler.predict(x.data(), n, 1, combined.data());
    double scale = sampler.fitScale();
    const double* lastCombined = combined.data() + (numSamples - 1) * n;
    double shift = lastCombined[0] - scale * lastSlot[0];
    double maxShiftSpread = 0.0;
    for (size_t i = 0; i < n; ++i)
      maxShiftSpread = std::max(
        maxShiftSpread,
        std::fabs((lastCombined[i] - scale * lastSlot[i]) - shift));
    check(maxShiftSpread < 1.0e-12,
          "the per-forest channel carries no fit scale of its own");

    bool slotsFinite = true;
    for (double v : replayed) slotsFinite &= std::isfinite(v);
    check(slotsFinite, "every saved per-forest slot is finite");

    ext_rng_destroy(localRng);
  }

  {
    SamplerOptions options;  // keepTrees off: the live-tree route
    ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(localRng, 5151u);
    Sampler<ConstantGaussianLeaf> sampler(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, options, spec, &localRng);

    Results results;
    sampler.run(120, 1, results);
    sampler.predictPerForest(x.data(), n, 1, replayed.data());
    check(maxPerForestGap(sampler, replayed.data()) < 1.0e-12,
          "live-tree per-forest replay reproduces each forest's own total");

    ext_rng_destroy(localRng);
  }

  printf("ok: amplitude per-forest replay\n");
}

// A forest multiplier that is zero - or indistinguishable from zero at the
// tolerance R's own almost-equal comparisons default to - excludes the row from
// that forest's leaf conditionals exactly. A constructed BCF combiner carries
// the neutral glue (a, b0, b1) = (1, 0, 1), so forest 1's control rows have a
// multiplier of exactly zero and must come back with exactly zero weight AND
// exactly zero response; the response half is load-bearing, since the chain
// reads that buffer arithmetically and the node kernels accumulate w * y. The
// snap band is pinned from BOTH sides so the constant cannot drift silently,
// and the treated rows pin that an ordinary multiplier is untouched.
static void testBCFZeroMultiplierSnap() {
  const size_t n = 8;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::vector<double> z(n), y(n), w(n), muFits(n), tauFits(n);
  for (size_t i = 0; i < n; ++i) {
    double t = static_cast<double>(i);
    z[i] = i % 2 == 0 ? 0.0 : 1.0;
    y[i] = 0.5 * t - 1.5;
    w[i] = 0.75 + 0.25 * t;
    muFits[i] = 0.3 * t - 0.7;
    tauFits[i] = 0.9 - 0.2 * t;
  }

  AmplitudeSpec spec;
  spec.z = z.data();
  AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

  std::vector<Forest<ConstantGaussianLeaf>> forests(2);
  for (size_t f = 0; f < 2; ++f) {
    forests[f].numTrees = 1;
    forests[f].leaf.scale = 1.0;
    forests[f].k = 2.0;
    forests[f].treeFits.assign(n, 0.0);
  }
  forests[0].totalFits = muFits;
  forests[1].totalFits = tauFits;

  double a, b0, b1;
  combiner.bcfGlue(a, b0, b1);
  check(a == 1.0 && b0 == 0.0 && b1 == 1.0,
        "a constructed BCF combiner holds the neutral (1, 0, 1) glue");

  ForestResponse fr =
    combiner.formForestResponse(1, forests, y.data(), w.data());
  bool excluded = true, treated = true, capped = true;
  for (size_t i = 0; i < n; ++i) {
    double resid = y[i] - a * muFits[i];
    if (z[i] == 0.0)
      excluded &= fr.response[i] == 0.0 && fr.weights[i] == 0.0;
    else
      treated &= fr.response[i] == resid && fr.weights[i] == w[i];
    // the condition cap the tolerance buys, on every row and every path; NaN
    // from a zero weight meeting an amplified response fails this comparison
    capped &= std::fabs(fr.response[i]) <= std::fabs(resid) * 0x1p26;
  }
  check(excluded,
        "a zero treatment multiplier gives exactly zero weight and response");
  check(treated, "a unit treatment multiplier passes the residual through");
  check(capped, "no reparameterized response exceeds the condition cap");

  // Both band edges, through the glue restore: 2^-27 is inside the band and
  // snaps whichever sign it carries, 2^-26 is the band itself and 2^-25 is
  // outside it, where the reparameterization is the plain division.
  ChainStateData glued;
  glued.hasAmplitudes = true;
  glued.a = 1.0;
  glued.b1 = 1.0;
  for (int sign = -1; sign <= 1; sign += 2) {
    glued.b0 = sign * 0x1p-27;
    combiner.restoreGlue(glued);
    ForestResponse edge =
      combiner.formForestResponse(1, forests, y.data(), w.data());
    bool inside = true;
    for (size_t i = 0; i < n; ++i)
      if (z[i] == 0.0)
        inside &= edge.response[i] == 0.0 && edge.weights[i] == 0.0;
    check(inside, "a multiplier inside the band snaps regardless of sign");
  }

  glued.b0 = 0x1p-25;
  combiner.restoreGlue(glued);
  ForestResponse outside =
    combiner.formForestResponse(1, forests, y.data(), w.data());
  bool divides = true;
  for (size_t i = 0; i < n; ++i) {
    if (z[i] != 0.0) continue;
    double resid = y[i] - muFits[i];
    divides &= outside.response[i] == resid / 0x1p-25 &&
               outside.weights[i] == w[i] * 0x1p-25 * 0x1p-25;
  }
  check(divides, "a multiplier outside the band still divides exactly");

  // The constant-leaf marginal a control-only leaf now feeds, accumulated from
  // the combiner's own output at the zero multiplier rather than from
  // hand-written zeros. Such a leaf carries no weight at all, so its marginal is
  // the honest zero of an uninformative leaf rather than a cancellation of two
  // O(1) terms - and it is zero on BOTH sides of a birth, so a control-only
  // split decision is exactly the prior's instead of the prior's to within a
  // contaminated 1e-15.
  {
    glued.b0 = 0.0;
    combiner.restoreGlue(glued);
    ForestResponse zeroed =
      combiner.formForestResponse(1, forests, y.data(), w.data());
    ConstantGaussianLeaf leaf{1.0 / std::sqrt(25.0)};
    double swP = 0.0, swrP = 0.0, swL = 0.0, swrL = 0.0, swR = 0.0, swrR = 0.0;
    for (size_t i = 0; i < n; ++i) {
      if (z[i] != 0.0) continue;
      double wi = zeroed.weights[i], wri = wi * zeroed.response[i];
      swP += wi;
      swrP += wri;
      if (i < n / 2) { swL += wi; swrL += wri; }
      else           { swR += wi; swrR += wri; }
    }
    double parent = leaf.logIntegratedLikelihood(2.0, 1.0, swP, swrP);
    double left = leaf.logIntegratedLikelihood(2.0, 1.0, swL, swrL);
    double right = leaf.logIntegratedLikelihood(2.0, 1.0, swR, swrR);
    check(parent == 0.0 && left == 0.0 && right == 0.0,
          "a control-only leaf's marginal is exactly zero");
    check(left + right - parent == 0.0,
          "a control-only birth's log MH delta is exactly zero");
  }

  printf("ok: BCF zero multiplier snaps to an exact exclusion\n");
}

// Two generators sit at the same stream position when their serialized states
// agree: the reading that says a pinned draw sequence was consumed exactly,
// with nothing taken, skipped or reordered around it.
static bool sameStreamPosition(const ext_rng* a, const ext_rng* b) {
  size_t length = ext_rng_getSerializedStateLength(a);
  std::vector<unsigned char> stateA(length), stateB(length);
  ext_rng_writeSerializedState(a, stateA.data());
  ext_rng_writeSerializedState(b, stateB.data());
  return stateA == stateB;
}

static ext_rng* makeSeamRng(unsigned int seed = 20260813u) {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, seed);
  return rng;
}

// The combiner at the seam a general per-forest basis family generalizes
// (docs/plans/multiforest-extension-surface.md, M4.0): forestMultiplier's two
// values at both of its call sites, combinedFits' blend, drawGlue's four-draw
// conditional order, and afterCombine's applied scale with its 1.0 skips.
// These are pins, not derivations - the point is that a generalization which
// reaches the same quantity by another route still has to land the same bits.
// Every fixture value is dyadic, so the reference arithmetic is exact and the
// comparisons can be bitwise without a tolerance to hide a wrong wiring.
static void testBCFCombinerSeam() {
  // ---- forestMultiplier's two values, at both call sites ----
  {
    const size_t n = 6;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);

    std::vector<double> z = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<double> mu = {0.5, -1.5, 2.0, 0.25, -0.75, 1.0};
    std::vector<double> tau = {1.0, 0.5, -2.0, 1.5, 0.25, -0.5};
    std::vector<double> y = {2.0, -1.0, 0.5, 3.0, -2.5, 1.25};
    std::vector<double> w = {0.5, 1.0, 2.0, 0.25, 4.0, 1.5};
    // three DISTINCT multipliers, so a wiring that reuses one forest's for the
    // other, or swaps the treated and control coefficients, cannot pass
    const double a = 1.5, b0 = -0.5, b1 = 0.25;

    AmplitudeSpec spec;
    spec.z = z.data();
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    ChainStateData glue;
    glue.hasAmplitudes = true;
    glue.a = a;
    glue.b0 = b0;
    glue.b1 = b1;
    combiner.restoreGlue(glue);

    std::vector<Forest<ConstantGaussianLeaf>> forests(2);
    forests[0].totalFits = mu;
    forests[1].totalFits = tau;

    size_t numTreated = 0;
    for (size_t i = 0; i < n; ++i) numTreated += z[i] != 0.0 ? 1 : 0;
    check(numTreated > 0 && numTreated < n,
          "the seam fixture carries both treated and control rows");

    // forest 0's multiplier is a on EVERY row (row-independent), and the
    // residual it divides is y net of the OTHER forest's b_z tau - the second
    // call site, which a generalization can miswire on its own
    ForestResponse prognostic =
      combiner.formForestResponse(0, forests, y.data(), w.data());
    bool muPair = true;
    for (size_t i = 0; i < n; ++i) {
      double bz = z[i] != 0.0 ? b1 : b0;
      muPair &= prognostic.response[i] == (y[i] - bz * tau[i]) / a &&
                prognostic.weights[i] == w[i] * a * a;
    }
    check(muPair,
          "the prognostic multiplier is a on every row, against a residual net "
          "of b_z tau");

    ForestResponse treatment =
      combiner.formForestResponse(1, forests, y.data(), w.data());
    bool tauPair = true;
    for (size_t i = 0; i < n; ++i) {
      double m = z[i] != 0.0 ? b1 : b0;
      tauPair &= treatment.response[i] == (y[i] - a * mu[i]) / m &&
                 treatment.weights[i] == w[i] * m * m;
    }
    check(tauPair,
          "the treatment multiplier is b1 on treated rows and b0 on control "
          "rows, against a residual net of a mu");

    const double* combined = combiner.combinedFits(forests);
    bool blend = true;
    for (size_t i = 0; i < n; ++i)
      blend &= combined[i] == a * mu[i] + (z[i] != 0.0 ? b1 : b0) * tau[i];
    check(blend, "the combined location is a mu + b_z tau, bitwise");

    // the snap belongs to formForestResponse's REPARAMETERIZATION, not to the
    // model: at a multiplier inside the zero band the treatment pair is
    // exactly zero on the control rows while the blend still carries the exact
    // b0 tau. A generalization that snaps once, in a shared effective-
    // multiplier helper both readers call, fails here and nowhere else.
    glue.b0 = 0x1p-27;
    combiner.restoreGlue(glue);
    ForestResponse snapped =
      combiner.formForestResponse(1, forests, y.data(), w.data());
    const double* exactBlend = combiner.combinedFits(forests);
    bool kept = true;
    for (size_t i = 0; i < n; ++i) {
      if (z[i] != 0.0) continue;
      kept &= snapped.response[i] == 0.0 && snapped.weights[i] == 0.0 &&
              exactBlend[i] == a * mu[i] + 0x1p-27 * tau[i];
    }
    check(kept,
          "the blend keeps the exact multiplier the reparameterization snaps");
  }

  // ---- drawGlue draws a, aVariance, b0, b1 in that order, and nothing else --
  // Zero fits and a unit prior variance leave every conditional mean at 0 and
  // every precision at exactly 1, so each drawn scalar IS its own draw from the
  // stream and the order is read off directly, with no reference arithmetic to
  // go stale. The ridge draw is deliberately absent: it belongs to
  // afterCombine, and a generalization that folds the two methods into one
  // moves this stream.
  {
    const size_t n = 4;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);
    std::vector<double> z = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> y = {0.5, -0.25, 1.0, 0.75};
    std::vector<double> zeros(n, 0.0);

    AmplitudeSpec spec;
    spec.z = z.data();
    spec.bPriorVariance = 1.0;

    auto zeroForests = [&](std::vector<Forest<ConstantGaussianLeaf>>& forests) {
      forests.clear();
      forests.resize(2);
      forests[0].totalFits = zeros;
      forests[1].totalFits = zeros;
    };

    {
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
      std::vector<Forest<ConstantGaussianLeaf>> forests;
      zeroForests(forests);
      ext_rng* live = makeSeamRng();
      ext_rng* reference = makeSeamRng();
      combiner.drawGlue(live, 1.0, y.data(), nullptr, forests);

      double refA = ext_rng_simulateStandardNormal(reference);
      double rate =
        0.5 * refA * refA + 0.5 * spec.aPriorScale * spec.aPriorScale;
      double refAVariance =
        1.0 / ext_rng_simulateGamma(reference, 1.0, 1.0 / rate);
      double refB0 = ext_rng_simulateStandardNormal(reference);
      double refB1 = ext_rng_simulateStandardNormal(reference);

      ChainStateData drawn;
      combiner.serializeGlue(drawn);
      check(drawn.a == refA && drawn.aVariance == refAVariance &&
              drawn.b0 == refB0 && drawn.b1 == refB1,
            "drawGlue draws a, then aVariance, then b0, then b1");
      check(sameStreamPosition(live, reference),
            "drawGlue consumes those four draws and no others");
      ext_rng_destroy(reference);
      ext_rng_destroy(live);
    }

    // the a block's skip takes NO draws, so b0 and b1 are the stream's first
    // two normals; the b block's skip stops the stream after aVariance
    {
      AmplitudeSpec noA = spec;
      noA.updateA = false;
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, noA);
      std::vector<Forest<ConstantGaussianLeaf>> forests;
      zeroForests(forests);
      ext_rng* live = makeSeamRng();
      ext_rng* reference = makeSeamRng();
      combiner.drawGlue(live, 1.0, y.data(), nullptr, forests);
      double refB0 = ext_rng_simulateStandardNormal(reference);
      double refB1 = ext_rng_simulateStandardNormal(reference);
      ChainStateData drawn;
      combiner.serializeGlue(drawn);
      check(drawn.a == 1.0 && drawn.aVariance == 1.0 && drawn.b0 == refB0 &&
              drawn.b1 == refB1 && sameStreamPosition(live, reference),
            "a pinned a block holds its two scalars and takes no draw");
      ext_rng_destroy(reference);
      ext_rng_destroy(live);
    }
    {
      AmplitudeSpec noB = spec;
      noB.updateB = false;
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, noB);
      std::vector<Forest<ConstantGaussianLeaf>> forests;
      zeroForests(forests);
      ext_rng* live = makeSeamRng();
      ext_rng* reference = makeSeamRng();
      combiner.drawGlue(live, 1.0, y.data(), nullptr, forests);
      double refA = ext_rng_simulateStandardNormal(reference);
      double rate =
        0.5 * refA * refA + 0.5 * spec.aPriorScale * spec.aPriorScale;
      (void) ext_rng_simulateGamma(reference, 1.0, 1.0 / rate);
      ChainStateData drawn;
      combiner.serializeGlue(drawn);
      check(drawn.a == refA && drawn.b0 == 0.0 && drawn.b1 == 1.0 &&
              sameStreamPosition(live, reference),
            "a pinned b block holds its two scalars and takes no draw");
      ext_rng_destroy(reference);
      ext_rng_destroy(live);
    }
  }

  // ---- afterCombine's applied scale, its GIG map, and its 1.0 skips ----
  {
    const size_t n = 4;
    std::vector<double> xDummy(n, 0.0);
    ColumnStore data;
    data.build(xDummy.data(), n, 1, 100);
    std::vector<double> z = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> muFits = {0.75, -0.5, 1.25, 0.5};
    std::vector<double> tauFits = {1.0, -0.25, 0.5, 2.0};
    const double leafA = 0.5, leafB = 0.25, tauLeaf = 0.125;
    const double a0 = 1.5, aVariance = 4.0;  // distinct, so a^2/aVariance bites

    AmplitudeSpec spec;
    spec.z = z.data();
    ChainStateData glue;
    glue.hasAmplitudes = true;
    glue.a = a0;
    glue.aVariance = aVariance;

    // one occupied leaf per prognostic tree: L is the leaf count, so a
    // two-tree forest puts the GIG exponent at (L - 1)/2 = 0.5, which the
    // off-by-one L/2 = 1.0 the move map suggests cannot reproduce
    auto build = [&](std::vector<Forest<ConstantGaussianLeaf>>& forests,
                     const std::vector<double>& leaves) {
      forests.clear();
      forests.resize(2);
      for (size_t f = 0; f < 2; ++f) {
        forests[f].leaf.scale = 1.0;
        forests[f].k = 2.0;
      }
      forests[0].numTrees = leaves.size();
      forests[0].totalFits = muFits;
      forests[0].indexBuffer.assign(n * leaves.size(), 0);
      forests[0].trees.resize(leaves.size());
      forests[0].muByTree.assign(leaves.size(), std::vector<double>(1, 0.0));
      for (size_t t = 0; t < leaves.size(); ++t) {
        forests[0].trees[t].initialize(forests[0].indexBuffer.data() + t * n, n);
        forests[0].muByTree[t][0] = leaves[t];
      }
      forests[1].numTrees = 1;
      forests[1].totalFits = tauFits;
      forests[1].indexBuffer.assign(n, 0);
      forests[1].trees.resize(1);
      forests[1].trees[0].initialize(forests[1].indexBuffer.data(), n);
      forests[1].muByTree.assign(1, std::vector<double>(1, tauLeaf));
    };

    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    combiner.restoreGlue(glue);
    std::vector<Forest<ConstantGaussianLeaf>> forests;
    build(forests, {leafA, leafB});

    ext_rng* live = makeSeamRng();
    ext_rng* reference = makeSeamRng();
    double c = combiner.afterCombine(forests, false, 0, live);
    // GIG((L - 1)/2, M / leafVar, a^2 / aVariance), with M the squared sum of
    // the occupied leaves and 1/leafVar the (k / leaf scale)^2 the forest
    // carries
    double squaredSum = leafA * leafA + leafB * leafB;
    double v = ext_rng_simulateGeneralizedInverseGaussian(
      reference, 0.5 * (2.0 - 1.0), squaredSum * (2.0 / 1.0) * (2.0 / 1.0),
      a0 * a0 / aVariance);
    check(c == std::sqrt(v),
          "the applied scale is the square root of the ridge's GIG draw");
    check(sameStreamPosition(live, reference),
          "afterCombine's move is that one draw and no other");

    ChainStateData moved;
    combiner.serializeGlue(moved);
    bool travelled = moved.a == a0 / c && moved.aVariance == aVariance &&
                     forests[0].muByTree[0][0] == leafA * c &&
                     forests[0].muByTree[1][0] == leafB * c;
    for (size_t i = 0; i < n; ++i)
      travelled &= forests[0].totalFits[i] == muFits[i] * c;
    check(travelled,
          "the amplitude travels to a0 / c while every prognostic leaf and fit "
          "scales by exactly c");
    bool untouched = forests[1].muByTree[0][0] == tauLeaf;
    for (size_t i = 0; i < n; ++i)
      untouched &= forests[1].totalFits[i] == tauFits[i];
    check(untouched, "the treatment forest is not on this ridge");
    ext_rng_destroy(reference);
    ext_rng_destroy(live);

    // the three reachable 1.0 returns, each of which must also leave the
    // stream where it found it - that is what keeps a skipped move bitwise.
    // The two non-finite guards below them are unreachable from a fixture.
    struct Skip {
      const char* what;
      std::vector<double> leaves;
      bool updateA;
    };
    const Skip skips[] = {{"a pinned amplitude", {leafA, leafB}, false},
                          {"a single occupied leaf", {leafA}, true},
                          {"an all-zero leaf sum", {0.0, 0.0}, true}};
    for (const Skip& skip : skips) {
      AmplitudeSpec skipSpec = spec;
      skipSpec.updateA = skip.updateA;
      AmplitudeForestCombiner<ConstantGaussianLeaf> skipped(data, skipSpec);
      skipped.restoreGlue(glue);
      std::vector<Forest<ConstantGaussianLeaf>> skipForests;
      build(skipForests, skip.leaves);
      ext_rng* skipRng = makeSeamRng();
      ext_rng* skipReference = makeSeamRng();
      double one = skipped.afterCombine(skipForests, false, 0, skipRng);
      ChainStateData held;
      skipped.serializeGlue(held);
      bool inert = one == 1.0 && held.a == a0 &&
                   sameStreamPosition(skipRng, skipReference);
      for (size_t i = 0; i < n; ++i)
        inert &= skipForests[0].totalFits[i] == muFits[i];
      std::string label =
        std::string("afterCombine returns 1.0 and moves nothing at ") +
        skip.what;
      check(inert, label.c_str());
      ext_rng_destroy(skipReference);
      ext_rng_destroy(skipRng);
    }
  }

  // ---- the same 1.0, through the hook the sweep exposes for these tests ----
  {
    const size_t n = 120, p = 2;
    std::uint64_t state = 20260813u;
    auto unif = [&]() {
      state ^= state << 13; state ^= state >> 7; state ^= state << 17;
      return static_cast<double>(state >> 11) * 0x1.0p-53;
    };
    std::vector<double> x(n * p), y(n), z(n);
    for (double& v : x) v = unif();
    for (size_t i = 0; i < n; ++i) {
      z[i] = unif() < 0.5 ? 1.0 : 0.0;
      y[i] = std::sin(3.0 * x[i]) + z[i] * (1.0 + x[i + n]) +
             0.2 * (unif() - 0.5);
    }

    SamplerOptions options;
    AmplitudeSpec spec;
    spec.mu.numTrees = 10;
    spec.tau.numTrees = 5;
    spec.updateA = false;
    spec.z = z.data();
    ext_rng* rng = makeSeamRng();
    Sampler<ConstantGaussianLeaf> sampler(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
      0.37804942330213542, options, spec, &rng);
    Results results;
    sampler.run(20, 0, results);

    double a0, b0, b1;
    sampler.chain(0).bcfGlue(a0, b0, b1);
    std::vector<double> before(n);
    sampler.forestTotalFits(0, 0, before.data());
    double one = sampler.chain(0).interweaveGlueRidge();
    double a1, b0p, b1p;
    sampler.chain(0).bcfGlue(a1, b0p, b1p);
    std::vector<double> after(n);
    sampler.forestTotalFits(0, 0, after.data());
    bool inert = one == 1.0 && a1 == a0;
    for (size_t i = 0; i < n; ++i) inert &= after[i] == before[i];
    check(inert, "the chain hook reports a skipped ridge move as exactly 1.0");
    ext_rng_destroy(rng);
  }

  printf("ok: BCF combiner seam pins\n");
}

// M4.1's own two facts, beside the seam pins rather than inside them.
//
// (1) The row sum's ASSOCIATION is observable. Under fused multiply-add
// contraction exactly one product in a sum escapes its own rounding - the one
// the closing add absorbs - and combinedFits absorbs FOREST 0's, which is what
// makes the K loop bitwise with the two-forest blend it replaces. The seam
// pin's reference is written as an ordinary expression and so inherits the test
// compiler's own contraction choice; this one spells both candidate
// associations with std::fma, so it states the association instead of assuming
// it, and counts the rows on which the two actually separate.
static void testCombinedFitsAssociation() {
  const size_t n = 64;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::uint64_t state = 20260813u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };
  std::vector<double> z(n), mu(n), tau(n);
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    mu[i] = 4.0 * (unif() - 0.5);
    tau[i] = 4.0 * (unif() - 0.5);
  }
  // deliberately NOT dyadic: a product's low bits have to be live or the two
  // associations agree everywhere and the pin says nothing
  const double a = 1.0834712, b0 = -0.31274, b1 = 0.9124367;

  AmplitudeSpec spec;
  spec.z = z.data();
  AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
  ChainStateData glue;
  glue.hasAmplitudes = true;
  glue.a = a;
  glue.b0 = b0;
  glue.b1 = b1;
  combiner.restoreGlue(glue);

  std::vector<Forest<ConstantGaussianLeaf>> forests(2);
  forests[0].totalFits = mu;
  forests[1].totalFits = tau;

  const double* combined = combiner.combinedFits(forests);
  bool absorbsPrognostic = true;
  size_t discriminating = 0, fusionVisible = 0, twiceRounded = 0,
         absorbsTreatment = 0;
  for (size_t i = 0; i < n; ++i) {
    double bz = z[i] != 0.0 ? b1 : b0;
    // lone multiplies, so each is its own correctly rounded product and the
    // twice-rounded sum below has no multiply left in it to contract
    double prognosticTerm = a * mu[i];
    double treatmentTerm = bz * tau[i];
    double twiceRoundedSum = prognosticTerm + treatmentTerm;
    double absorbsForest0 = std::fma(a, mu[i], treatmentTerm);
    double absorbsForest1 = std::fma(bz, tau[i], prognosticTerm);
    absorbsPrognostic &= combined[i] == absorbsForest0;
    discriminating += absorbsForest0 != absorbsForest1 ? 1 : 0;
    // a row where fusing forest 0's product changes the answer at all: only
    // there can this fixture see whether the target contracts
    if (absorbsForest0 == twiceRoundedSum) continue;
    ++fusionVisible;
    twiceRounded += combined[i] == twiceRoundedSum ? 1 : 0;
    absorbsTreatment +=
      combined[i] == absorbsForest1 && absorbsForest1 != absorbsForest0 ? 1 : 0;
  }
  check(fusionVisible > 0,
        "the association fixture puts forest 0's fused product to work on some "
        "row");

  // The DIRECTION is only a question where the target fuses at all. A build
  // with no fused multiply-add - an x86-64 baseline without FMA, or any
  // -ffp-contract=off build - rounds both products on every row and would fail
  // a direction assertion for a reason that says nothing about this code, so
  // the assertion is skipped rather than reported red. What is NOT skipped is
  // the wiring regression it exists for: absorbing the LAST forest's product.
  if (fusionVisible > 0 && twiceRounded == fusionVisible) {
    printf("ok: combined fits association (no fused multiply-add on this "
           "target; direction assertion skipped, %zu of %zu rows "
           "discriminating)\n", discriminating, n);
    return;
  }

  check(discriminating > 0,
        "the association fixture separates the two roundings on some row");
  check(absorbsTreatment == 0,
        "the combination never absorbs the last forest's product");
  check(absorbsPrognostic,
        "the combination absorbs forest 0's product into the closing add");

  printf("ok: combined fits association (%zu of %zu rows discriminating)\n",
         discriminating, n);
}

// (2) The ORDERING contract on the sole basis-mutation route (M4.3's SUBSUME:
// synthesis is construction-only, setForestBasis is the only operation, so
// LAST INSTALL WINS per forest and the two orderings of a widen and a swap
// commute). Five arms, each red against a different wrong wiring: a route that
// re-synthesizes on install (1, 2), one that rebuilds the amplitudes rather
// than carrying them (2, 3), one that keeps a borrowed z beside the basis (4),
// and one that selects the draw path on widths rather than values (5).
static void testForestBasisOrdering() {
  const size_t n = 6;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::vector<double> z = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  // the COMPLEMENT, so a stale basis is wrong on all six rows rather than some
  std::vector<double> swapped(2 * n);
  for (size_t i = 0; i < n; ++i) {
    swapped[2 * i] = z[i];
    swapped[2 * i + 1] = 1.0 - z[i];
  }
  // a DISCRIMINATING second prognostic column: nonzero and row-varying, so a
  // multiplier that dropped it moves every row (unit values vacate pins)
  std::vector<double> wide(2 * n);
  for (size_t i = 0; i < n; ++i) {
    wide[2 * i] = 1.0;
    wide[2 * i + 1] = 0.5 + 0.25 * static_cast<double>(i);
  }
  std::vector<double> mu = {0.5, -1.5, 2.0, 0.25, -0.75, 1.0};
  std::vector<double> tau = {1.0, 0.5, -2.0, 1.5, 0.25, -0.5};
  std::vector<double> y = {2.0, -1.0, 0.5, 3.0, -2.5, 1.25};
  const double a = 1.5, b0 = -0.5, b1 = 0.25;

  auto build = [&](AmplitudeForestCombiner<ConstantGaussianLeaf>& combiner) {
    ChainStateData glue;
    glue.hasAmplitudes = true;
    glue.a = a;
    glue.b0 = b0;
    glue.b1 = b1;
    combiner.restoreGlue(glue);
  };

  std::vector<Forest<ConstantGaussianLeaf>> forests(2);
  forests[0].totalFits = mu;
  forests[1].totalFits = tau;

  AmplitudeSpec spec;
  spec.z = z.data();

  // --- Arm 1, widen then swap. The widening survives the swap, the prognostic
  // multiplier still reads its second column on every row, and the swapped
  // forest's own block is bitwise what it was. RED under a route that
  // re-synthesizes forest 0 on any install.
  {
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    build(combiner);
    check(combiner.installForestBasis(0, wide.data(), 2),
          "a wider prognostic basis installs");
    ChainStateData widened;
    combiner.serializeGlue(widened);
    // the layout moved, so the added coordinate entered at the neutral 1.0
    check(widened.amplitudeWidths == std::vector<size_t>{2, 2} &&
            widened.amplitudes[0] == a && widened.amplitudes[1] == 1.0 &&
            widened.amplitudes[2] == b0 && widened.amplitudes[3] == b1,
          "a widening carries every block and enters the new coordinate at 1");

    check(combiner.installForestBasis(1, swapped.data(), 2),
          "the swap installs on the widened combiner");
    ChainStateData after;
    combiner.serializeGlue(after);
    check(after.amplitudeWidths == std::vector<size_t>{2, 2} &&
            after.amplitudes == widened.amplitudes,
          "a width-preserving swap leaves every amplitude bitwise");

    const double* combined = combiner.combinedFits(forests);
    bool reads = true;
    for (size_t i = 0; i < n; ++i) {
      double m0 = a * wide[2 * i] + 1.0 * wide[2 * i + 1];
      double bz = swapped[2 * i + 1] != 0.0 ? b1 : b0;
      reads &= combined[i] == m0 * mu[i] + bz * tau[i];
    }
    check(reads, "the widened prognostic basis survives the swap on every row");
  }

  // --- Arm 2, swap then widen. The swapped values survive the layout move and
  // the swapped forest's coordinates are carried bitwise to their new offsets.
  {
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    build(combiner);
    check(combiner.installForestBasis(1, swapped.data(), 2),
          "the swap installs first");
    check(combiner.installForestBasis(0, wide.data(), 2),
          "the widening follows it");
    ChainStateData after;
    combiner.serializeGlue(after);
    check(after.amplitudeWidths == std::vector<size_t>{2, 2} &&
            after.amplitudes[0] == a && after.amplitudes[1] == 1.0 &&
            after.amplitudes[2] == b0 && after.amplitudes[3] == b1,
          "the widening carries the swapped forest's block to its new offset");

    const double* combined = combiner.combinedFits(forests);
    bool reads = true;
    for (size_t i = 0; i < n; ++i) {
      double m0 = a * wide[2 * i] + 1.0 * wide[2 * i + 1];
      double bz = swapped[2 * i + 1] != 0.0 ? b1 : b0;
      reads &= combined[i] == m0 * mu[i] + bz * tau[i];
    }
    check(reads, "both orderings leave the same combination");
  }

  // --- Arm 3, a width-preserving CANONICAL reinstall is the bitwise identity
  // on every amplitude. This is the pin between M4.3 and a bcf-equivalence
  // re-record: the shipped mid-life z swap is exactly this install.
  {
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    build(combiner);
    ChainStateData before;
    combiner.serializeGlue(before);
    check(combiner.installForestBasis(1, swapped.data(), 2),
          "the canonical reinstall runs");
    ChainStateData after;
    combiner.serializeGlue(after);
    check(before.amplitudes == after.amplitudes &&
            before.amplitudeWidths == after.amplitudeWidths &&
            before.amplitudeVariances == after.amplitudeVariances,
          "a width-preserving install is the bitwise identity on the glue");

    // and the reparameterization reads the NEW basis, not the constructed one
    ForestResponse treatment =
      combiner.formForestResponse(1, forests, y.data(), nullptr);
    bool pair = true;
    for (size_t i = 0; i < n; ++i) {
      double m = swapped[2 * i + 1] != 0.0 ? b1 : b0;
      pair &= treatment.response[i] == (y[i] - a * mu[i]) / m &&
              treatment.weights[i] == m * m;
    }
    check(pair, "the reparameterization reads the installed basis");
  }

  // --- Arm 4, the z divergence: the leg Arm 3 structurally cannot see, since
  // Arm 3 asserts amplitudes and the defect is in the NEXT draw. A
  // width-preserving swap to a DIFFERENT complementary pair must change that
  // draw's PARTITION. RED against a retained glue_.z (old z, new basis), green
  // once drawShippedGlue partitions from basis[1].
  {
    // an ASYMMETRIC fit surface, so the two groups' conditionals differ and a
    // mis-grouped draw is visible in the values rather than only in the stream
    std::vector<Forest<ConstantGaussianLeaf>> drawForests(2);
    drawForests[0].totalFits = std::vector<double>(n, 0.0);
    drawForests[1].totalFits = {1.0, 1.0, 1.0, 4.0, 4.0, 4.0};
    std::vector<double> drawY = {1.0, 1.0, 1.0, 8.0, 8.0, 8.0};
    // rows 0-2 control, rows 3-5 treated under the constructed indicator
    std::vector<double> grouped = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<double> groupedBasis(2 * n), flippedBasis(2 * n);
    for (size_t i = 0; i < n; ++i) {
      groupedBasis[2 * i] = 1.0 - grouped[i];
      groupedBasis[2 * i + 1] = grouped[i];
      flippedBasis[2 * i] = grouped[i];
      flippedBasis[2 * i + 1] = 1.0 - grouped[i];
    }
    AmplitudeSpec groupedSpec;
    groupedSpec.z = grouped.data();
    groupedSpec.updateA = false;  // isolate the b block

    auto drawPair = [&](const std::vector<double>& basis, double& out0,
                        double& out1) {
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, groupedSpec);
      build(combiner);
      check(combiner.installForestBasis(1, basis.data(), 2),
            "the grouped basis installs");
      ext_rng* rng = makeSeamRng();
      combiner.drawGlue(rng, 1.0, drawY.data(), nullptr, drawForests);
      ext_rng_destroy(rng);
      double aOut = 0.0;
      check(combiner.bcfGlue(aOut, out0, out1), "the K = 2 reading holds");
    };

    double keptB0 = 0.0, keptB1 = 0.0, flippedB0 = 0.0, flippedB1 = 0.0;
    drawPair(groupedBasis, keptB0, keptB1);
    drawPair(flippedBasis, flippedB0, flippedB1);
    // the SAME rng, the same data, the complementary grouping: the two blocks
    // must trade places rather than reproduce
    check(keptB0 != flippedB0 && keptB1 != flippedB1,
          "a width-preserving swap to a different pair changes the partition");
  }

  // --- Arm 5, draw-path selection is a VALUE predicate, not a width test.
  // Two decisive discriminators, one per half, each keyed on something the
  // two-scalar path structurally CANNOT do.
  //
  // (i) it cannot write a second prognostic coordinate at all: it addresses
  // (a, b0, b1) and nothing else, so under q0 = 2 the added coordinate would
  // stay at the neutral 1.0 the widening entered it at.
  //
  // (ii) it never reads basis[1] beyond the treated indicator - it forms two
  // disjoint group accumulators keyed on it - so changing the CONTROL column
  // alone would leave its draw bitwise unmoved. The general conditional reads
  // the whole row, so it must move. This is the arm a width-only predicate
  // (K == 2, q0 == 1, q1 == 2) passes wrongly: the 0.25/0.75 pair has exactly
  // that shape.
  {
    std::vector<double> weights = {0.5, 1.5, 2.0, 0.25, 1.0, 3.0};
    auto draw = [&](const std::vector<double>* prognostic,
                    const std::vector<double>* treatment) {
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
      build(combiner);
      if (prognostic != nullptr)
        check(combiner.installForestBasis(0, prognostic->data(), 2),
              "the q0 = 2 basis installs");
      if (treatment != nullptr)
        check(combiner.installForestBasis(1, treatment->data(), 2),
              "the non-canonical width-2 basis installs");
      ext_rng* rng = makeSeamRng();
      combiner.drawGlue(rng, 1.0, y.data(), weights.data(), forests);
      ext_rng_destroy(rng);
      ChainStateData drawn;
      combiner.serializeGlue(drawn);
      return drawn.amplitudes;
    };

    std::vector<double> widened = draw(&wide, nullptr);
    check(widened.size() == 4 && widened[1] != 1.0,
          "a q0 = 2 basis routes the draw around the two-scalar path, which "
          "cannot write the second prognostic coordinate");

    std::vector<double> continuousA(2 * n), continuousB(2 * n);
    for (size_t i = 0; i < n; ++i) {
      continuousA[2 * i] = 0.25;
      continuousA[2 * i + 1] = 0.75;
      // the SAME treated column, a different control column: invisible to the
      // two-scalar path, load-bearing for the general one
      continuousB[2 * i] = 3.0 + static_cast<double>(i);
      continuousB[2 * i + 1] = 0.75;
    }
    std::vector<double> drawnA = draw(nullptr, &continuousA);
    std::vector<double> drawnB = draw(nullptr, &continuousB);
    check(drawnA.size() == 3 && drawnB.size() == 3 && drawnA != drawnB,
          "a NON-canonical width-2 basis routes it around too: the control "
          "column moves the draw, which the two-scalar path never reads");

    // (iii) and the canonical shape is still ON the two-scalar path. A
    // canonical pair's control column is DETERMINED by its treated column -
    // 1 - z, and nothing else is canonical - so the discriminating pair holds
    // the treated column FIXED and moves only the control column off that
    // determined value. Both bases then carry the same partition, which is all
    // the two-scalar path reads: it would draw them IDENTICALLY. The general
    // conditional reads the whole row, so it must not. The correct routing
    // therefore separates them - the canonical one keeps the two-scalar path,
    // the perturbed one is moved off it - and a width-only predicate, which
    // would leave BOTH on the two-scalar path, collapses them and turns this
    // red. Unlike the pair it replaces (built from identical expressions, true
    // under any implementation), nothing here is satisfied by construction.
    std::vector<double> canonical(2 * n), controlPerturbed(2 * n);
    for (size_t i = 0; i < n; ++i) {
      canonical[2 * i] = 1.0 - z[i];
      canonical[2 * i + 1] = z[i];
      // the same treated column, a control column a canonical pair could not
      // have: 2 where the canonical one is forced to 1
      controlPerturbed[2 * i] = 2.0 * (1.0 - z[i]);
      controlPerturbed[2 * i + 1] = z[i];
    }
    check(draw(nullptr, &canonical) != draw(nullptr, &controlPerturbed),
          "the canonical shape keeps the two-scalar path, which cannot see "
          "the control column its perturbed twin moves");
  }

  printf("ok: forest basis ordering\n");
}

// (3) The GENERAL q-variate amplitude conditional. Every pin above drives the
// shipped K = 2 shape, which keeps its two-scalar path (drawGlue says why), so
// none of them enters this code at all. Four arms: the conditional recovers a
// known closed-form Gaussian posterior over a CONTINUOUS, non-orthogonal basis
// (the case the factorization exists for); it agrees with the two-scalar path
// on the shipped shape, which is the whole evidence for keeping both; the path
// predicate reads the PRIOR a spec carries and not only its basis shape; and at
// K = 3 every per-forest array is sized by the forest count, which is the read
// that went out of bounds before.
static void testGeneralAmplitudeConditional() {
  const size_t n = 40;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::uint64_t state = 20260814u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };
  std::vector<double> z(n), y(n), w(n), mu(n), tau(n), other(n), wide(2 * n);
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    y[i] = 4.0 * (unif() - 0.5);
    w[i] = 0.5 + unif();
    mu[i] = 2.0 * (unif() - 0.5);
    tau[i] = 2.0 * (unif() - 0.5);
    other[i] = 2.0 * (unif() - 0.5);
    // an intercept and a continuous modifier: the crossproduct's off-diagonal
    // is then nonzero, so the unit triangles are NOT identity and the solve is
    // the general one rather than two scalar draws in disguise
    wide[2 * i] = 1.0;
    wide[2 * i + 1] = 3.0 * unif() - 1.0;
  }
  const double sigma = 0.8, aFixed = 1.25, priorVariance = 0.75;
  const double invSigmaSq = 1.0 / (sigma * sigma);

  // ---- the conditional recovers its own closed form ----
  {
    AmplitudeSpec spec;
    spec.z = z.data();
    spec.updateA = false;  // pin forest 0, so forest 1's block is all that moves
    spec.bPriorVariance = priorVariance;
    spec.generalAmplitudeDraw = true;
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    combiner.installForestBasis(1, wide.data(), 2);
    ChainStateData glue;
    glue.hasAmplitudes = true;
    glue.a = aFixed;
    glue.b0 = 0.0;
    glue.b1 = 0.0;
    combiner.restoreGlue(glue);

    std::vector<Forest<ConstantGaussianLeaf>> forests(2);
    forests[0].totalFits = mu;
    forests[1].totalFits = tau;

    // P = I / priorVariance + sum_i w_i x_i x_i' / sigma^2 and
    // num = sum_i w_i x_i r_i / sigma^2, with x_i the basis row scaled by the
    // forest's own fit and r_i the residual net of a mu
    double p00 = 1.0 / priorVariance, p11 = p00, p01 = 0.0;
    double num0 = 0.0, num1 = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double x0 = wide[2 * i] * tau[i], x1 = wide[2 * i + 1] * tau[i];
      double r = y[i] - aFixed * mu[i];
      p00 += w[i] * x0 * x0 * invSigmaSq;
      p01 += w[i] * x0 * x1 * invSigmaSq;
      p11 += w[i] * x1 * x1 * invSigmaSq;
      num0 += w[i] * x0 * r * invSigmaSq;
      num1 += w[i] * x1 * r * invSigmaSq;
    }
    double det = p00 * p11 - p01 * p01;
    check(det > 0.0 && std::fabs(p01) > 1e-3,
          "the conditional fixture is positive definite and NOT orthogonal");
    double v00 = p11 / det, v11 = p00 / det, v01 = -p01 / det;
    double mean0 = v00 * num0 + v01 * num1, mean1 = v01 * num0 + v11 * num1;

    const size_t draws = 20000;
    double sum0 = 0.0, sum1 = 0.0, sq0 = 0.0, sq1 = 0.0, cross = 0.0;
    ext_rng* rng = makeSeamRng();
    double aOut = 0.0, b0Out = 0.0, b1Out = 0.0;
    for (size_t d = 0; d < draws; ++d) {
      combiner.drawGlue(rng, sigma, y.data(), w.data(), forests);
      combiner.bcfGlue(aOut, b0Out, b1Out);
      sum0 += b0Out; sum1 += b1Out;
      sq0 += b0Out * b0Out; sq1 += b1Out * b1Out; cross += b0Out * b1Out;
    }
    ext_rng_destroy(rng);
    double scale = 1.0 / static_cast<double>(draws);
    double m0 = sum0 * scale, m1 = sum1 * scale;
    double c00 = sq0 * scale - m0 * m0, c11 = sq1 * scale - m1 * m1;
    double c01 = cross * scale - m0 * m1;
    // 4 Monte Carlo standard errors on each mean, 5% on each second moment
    double se0 = std::sqrt(v00 * scale), se1 = std::sqrt(v11 * scale);
    check(std::fabs(m0 - mean0) < 4.0 * se0 &&
            std::fabs(m1 - mean1) < 4.0 * se1,
          "the q-variate draw's mean is the conditional's mean");
    check(std::fabs(c00 / v00 - 1.0) < 0.05 &&
            std::fabs(c11 / v11 - 1.0) < 0.05 &&
            std::fabs(c01 / v01 - 1.0) < 0.05,
          "the q-variate draw's covariance is the conditional's inverse "
          "precision, off-diagonal included");
    check(aOut == aFixed, "a pinned block takes no draw and does not move");
    printf("ok: general amplitude conditional (rho %.3f, mean %.4f/%.4f "
           "against %.4f/%.4f)\n", c01 / std::sqrt(c00 * c11), m0, m1, mean0,
           mean1);
  }

  // ---- the general path and the two-scalar path are the same conditional ----
  // They cannot be BITWISE (drawGlue records the measurement: the b block's
  // per-row products are formed before a branch and accumulated inside it,
  // which fuses unevenly against a straight-line loop), so this asserts what
  // is true - one conditional, one rng stream, agreeing to rounding.
  {
    AmplitudeSpec shipped;
    shipped.z = z.data();
    shipped.bPriorVariance = priorVariance;
    AmplitudeSpec general = shipped;
    general.generalAmplitudeDraw = true;

    double drawn[2][4];
    for (size_t arm = 0; arm < 2; ++arm) {
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(
        data, arm == 0 ? shipped : general);
      std::vector<Forest<ConstantGaussianLeaf>> forests(2);
      forests[0].totalFits = mu;
      forests[1].totalFits = tau;
      ext_rng* rng = makeSeamRng();
      combiner.drawGlue(rng, sigma, y.data(), w.data(), forests);
      ChainStateData out;
      combiner.serializeGlue(out);
      drawn[arm][0] = out.a; drawn[arm][1] = out.aVariance;
      drawn[arm][2] = out.b0; drawn[arm][3] = out.b1;
      ext_rng_destroy(rng);
    }
    bool agrees = true;
    for (size_t j = 0; j < 4; ++j)
      agrees &= std::fabs(drawn[0][j] - drawn[1][j]) <=
                1e-12 * std::fabs(drawn[0][j]);
    check(agrees,
          "the general conditional reproduces the shipped two-scalar draw to "
          "rounding on the same stream");
    // The a block IS bitwise, both accumulators reproducing under weights, so
    // it pins the general path's solve arithmetic end to end: at q = 1 the
    // square-root-free route divides the moment by the pivot ONCE and lands on
    // the shipped n / P + e / sqrt(P) exactly. Swapping in the Cholesky the
    // leaf models use divides twice by sqrt(P) instead and moves this by an
    // ulp - the one assertion in the suite that sees such a substitution
    // through the engine rather than in the helper (test_model.cpp).
    check(drawn[0][0] == drawn[1][0],
          "the general conditional's scalar block is BITWISE the shipped a "
          "draw, which is what the square-root-free solve buys");
  }

  // ---- the path predicate reads the prior, not only the basis shape ----
  // A canonical K = 2 pair whose forest 1 declares a positive half-Cauchy scale
  // is a scale MIXTURE on that block: its prior variance is a live auxiliary.
  // The two-scalar path refreshes forest 0's and no other, so a predicate that
  // read the basis alone would hold forest 1's variance at its declared value
  // for every sweep - a different model, not a different rounding.
  {
    std::vector<double> indicator(2 * n);
    for (size_t i = 0; i < n; ++i) {
      indicator[2 * i] = 1.0 - z[i];
      indicator[2 * i + 1] = z[i];
    }
    const double declared = 0.5;
    AmplitudeSpec routed;
    routed.forests.resize(2);
    routed.forests[0].amplitudePriorScale = 2.0;
    routed.forests[1].basis = indicator.data();
    routed.forests[1].numBasisColumns = 2;
    routed.forests[1].amplitudePriorVariance = declared;
    routed.forests[1].amplitudePriorScale = 1.5;
    // the same spec off the predicate entirely: whatever the routing decides,
    // this arm is the general sweep, so an equality against it is the routing
    // read out rather than a second copy of the arithmetic
    AmplitudeSpec forced = routed;
    forced.generalAmplitudeDraw = true;

    ChainStateData out[2];
    for (size_t arm = 0; arm < 2; ++arm) {
      AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(
        data, arm == 0 ? routed : forced);
      std::vector<Forest<ConstantGaussianLeaf>> forests(2);
      forests[0].totalFits = mu;
      forests[1].totalFits = tau;
      ext_rng* rng = makeSeamRng();
      for (size_t sweep = 0; sweep < 40; ++sweep)
        combiner.drawGlue(rng, sigma, y.data(), w.data(), forests);
      combiner.serializeGlue(out[arm]);
      ext_rng_destroy(rng);
    }
    check(out[0].amplitudeVariances[1] != declared &&
            std::isfinite(out[0].amplitudeVariances[1]) &&
            out[0].amplitudeVariances[1] > 0.0,
          "a scale mixture declared past forest 0 is drawn, not frozen at the "
          "variance it was declared with");
    check(out[0].amplitudes == out[1].amplitudes &&
            out[0].amplitudeVariances == out[1].amplitudeVariances,
          "the spec routes to the general sweep bitwise, which is the path "
          "that refreshes a second scale mixture");
    printf("ok: amplitude path selection reads the prior (forest 1 variance "
           "%.4f against a declared %.4f)\n", out[0].amplitudeVariances[1],
           declared);
  }

  // ---- K = 3: every per-forest array is sized by the forest count ----
  // Before the K-length sizing the basis vector held bcf's two entries however
  // many forests the chain carried, and forest 2's multiplier read past it.
  {
    AmplitudeSpec spec;
    spec.z = z.data();
    spec.bPriorVariance = priorVariance;
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec, 3);
    std::vector<Forest<ConstantGaussianLeaf>> forests(3);
    forests[0].totalFits = mu;
    forests[1].totalFits = tau;
    forests[2].totalFits = other;

    ext_rng* rng = makeSeamRng();
    combiner.drawGlue(rng, sigma, y.data(), w.data(), forests);
    ext_rng_destroy(rng);

    const double* combined = combiner.combinedFits(forests);
    ForestResponse third =
      combiner.formForestResponse(2, forests, y.data(), w.data());
    bool finite = true;
    for (size_t i = 0; i < n; ++i)
      finite &= std::isfinite(combined[i]) && std::isfinite(third.response[i]) &&
                std::isfinite(third.weights[i]);
    check(finite,
          "a third forest combines, reparameterizes and draws its amplitude");
    printf("ok: general amplitude conditional at K = 3\n");
  }
}

// (4) The per-forest ASIS rescale, on the two things only the general move can
// be asked: that it preserves the prior along the orbit at q > 1 (the exponent
// rule, whose off-by-one this arm rejects), and that the combined fit it leaves
// behind is the one it found.
static void testGeneralAmplitudeRidge() {
  const size_t n = 8, numTrees = 3;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::vector<double> z(n), mu(n), tau(n);
  for (size_t i = 0; i < n; ++i) {
    z[i] = i % 2 == 0 ? 0.0 : 1.0;
    mu[i] = 0.5 * static_cast<double>(i) - 1.5;
    tau[i] = 1.25 - 0.25 * static_cast<double>(i);
  }
  const double leafScale = 1.0, k = 2.0;
  const double leafVariance = (leafScale / k) * (leafScale / k);
  const double priorVariance = 0.75;

  AmplitudeSpec spec;
  spec.z = z.data();
  spec.bPriorVariance = priorVariance;
  spec.updateA = false;   // forest 0 holds, so the treatment ridge is alone
  spec.ridgeB = true;     // the b-move: OFF for bcf, ON here
  AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

  std::vector<Forest<ConstantGaussianLeaf>> forests(2);
  std::vector<double> indices(0);
  for (size_t f = 0; f < 2; ++f) {
    forests[f].leaf.scale = leafScale;
    forests[f].k = k;
  }
  forests[0].numTrees = 1;
  forests[0].totalFits = mu;
  forests[0].indexBuffer.assign(n, 0);
  forests[0].trees.resize(1);
  forests[0].trees[0].initialize(forests[0].indexBuffer.data(), n);
  forests[0].muByTree.assign(1, std::vector<double>(1, 0.5));
  forests[1].numTrees = numTrees;
  forests[1].totalFits = tau;
  forests[1].indexBuffer.assign(n * numTrees, 0);
  forests[1].trees.resize(numTrees);
  forests[1].muByTree.assign(numTrees, std::vector<double>(1, 0.0));
  for (size_t t = 0; t < numTrees; ++t)
    forests[1].trees[t].initialize(forests[1].indexBuffer.data() + t * n, n);

  ext_rng* rng = makeSeamRng();

  // ---- the move preserves the prior along the orbit ----
  // The likelihood is constant on the orbit, so the move preserves the
  // posterior IFF it preserves the prior's along-orbit conditional: draw
  // (b0, b1, leaves) from the prior, apply ONE move, and the pushed sample must
  // still be a prior draw. This is the pure-R prototype from
  // docs/design/multiplier-combiner.md, "The exponent rule", run in the
  // engine instead, on second moments rather than KS - at L = 3 and q = 2 the
  // exponent is (L - q)/2 = 0.5, and the off-by-one (L - q + 1)/2 the naive
  // move-map Jacobian gives inflates the leaves by more than 10%.
  const size_t replicates = 20000;
  double amplitudeSquares = 0.0, leafSquares = 0.0;
  ChainStateData drawnGlue;
  drawnGlue.hasAmplitudes = true;
  drawnGlue.a = 1.0;
  for (size_t r = 0; r < replicates; ++r) {
    double priorSd = std::sqrt(priorVariance);
    drawnGlue.b0 = priorSd * ext_rng_simulateStandardNormal(rng);
    drawnGlue.b1 = priorSd * ext_rng_simulateStandardNormal(rng);
    combiner.restoreGlue(drawnGlue);
    for (size_t t = 0; t < numTrees; ++t)
      forests[1].muByTree[t][0] =
        std::sqrt(leafVariance) * ext_rng_simulateStandardNormal(rng);

    combiner.afterCombine(forests, false, 0, rng);

    ChainStateData moved;
    combiner.serializeGlue(moved);
    amplitudeSquares += moved.b0 * moved.b0 + moved.b1 * moved.b1;
    for (size_t t = 0; t < numTrees; ++t)
      leafSquares += forests[1].muByTree[t][0] * forests[1].muByTree[t][0];
  }
  double amplitudeMoment =
    amplitudeSquares / (2.0 * static_cast<double>(replicates)) / priorVariance;
  double leafMoment =
    leafSquares /
    (static_cast<double>(numTrees) * static_cast<double>(replicates)) /
    leafVariance;
  check(std::fabs(amplitudeMoment - 1.0) < 0.04 &&
          std::fabs(leafMoment - 1.0) < 0.04,
        "the q-variate ridge leaves the prior where it found it, at the "
        "(L - q)/2 exponent");
  printf("ok: general amplitude ridge (amplitude %.4f, leaf %.4f, both against "
         "1)\n", amplitudeMoment, leafMoment);

  // ---- the combined fit is invariant on the orbit ----
  {
    forests[1].totalFits = tau;
    for (size_t t = 0; t < numTrees; ++t)
      forests[1].muByTree[t][0] = 0.25 * static_cast<double>(t + 1);
    drawnGlue.b0 = -0.6;
    drawnGlue.b1 = 0.4;
    combiner.restoreGlue(drawnGlue);
    const double* combined = combiner.combinedFits(forests);
    std::vector<double> before(combined, combined + n);
    // the reported scale is the REPORTED forest's, so a move at ANOTHER forest
    // returns 1.0 having moved - the convention the base virtual now states,
    // and the reason this arm's ran-at-all witness is the amplitude itself
    double reported = combiner.afterCombine(forests, false, 0, rng);
    const double* after = combiner.combinedFits(forests);
    ChainStateData moved;
    combiner.serializeGlue(moved);
    bool ran = moved.b0 != drawnGlue.b0 && moved.b1 != drawnGlue.b1 &&
               forests[1].muByTree[0][0] != 0.25;
    bool invariant = true;
    for (size_t i = 0; i < n; ++i)
      invariant &= std::fabs(after[i] - before[i]) <=
                   8.0 * std::numeric_limits<double>::epsilon() *
                     (std::fabs(before[i]) + 1.0);
    check(ran, "the treatment ridge travelled both amplitudes and the leaves");
    check(invariant, "the treatment ridge leaves the combined fit where it "
                     "found it");
    check(reported == 1.0,
          "afterCombine returns 1.0 for a held reported forest while another "
          "forest travels");
  }

  ext_rng_destroy(rng);
}

// (5) The amplitude blocks are addressed THROUGH the per-forest offsets. With a
// two-column prognostic basis the layout is (a0, a1)(b0, b1), so the wire
// format's b0/b1 live at 2 and 3; the retired accessors read 1 and 2 - forest
// 0's second coordinate and forest 1's first. Unconstructible while every
// prognostic basis was one all-ones column, which is why nothing above catches
// it.
static void testAmplitudeOffsetIndexing() {
  const size_t n = 6;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::vector<double> z = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  std::vector<double> mu = {0.5, -1.5, 2.0, 0.25, -0.75, 1.0};
  std::vector<double> tau = {1.0, 0.5, -2.0, 1.5, 0.25, -0.5};
  std::vector<double> wide(2 * n);
  for (size_t i = 0; i < n; ++i) {
    wide[2 * i] = 1.0;
    wide[2 * i + 1] = 0.5 * static_cast<double>(i) - 1.0;
  }
  // four DISTINCT amplitudes, so an accessor off by one coordinate cannot pass
  const double a = 1.5, b0 = -0.5, b1 = 0.25;

  AmplitudeSpec spec;
  spec.z = z.data();
  AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
  combiner.installForestBasis(0, wide.data(), 2);

  // the RAGGED spelling, which is the only one that can name this layout: the
  // widened prognostic block's second coordinate is where the retired
  // three-scalar reading put b0
  ChainStateData glue;
  glue.hasAmplitudes = true;
  glue.amplitudeWidths = {2, 2};
  glue.amplitudes = {a, 1.0, b0, b1};
  glue.amplitudeVariances = {1.0, 1.0};
  combiner.restoreGlue(glue);

  // The round trip alone is NOT the guard: restore and report went through the
  // same accessors, so the retired pair aliased both together and reads its own
  // writes back (measured). The multiplier arm below is what catches it, by
  // reading the blocks through a third path that always used the offsets.
  ChainStateData reported;
  combiner.serializeGlue(reported);
  check(reported.amplitudeWidths == std::vector<size_t>{2, 2} &&
          reported.amplitudes == glue.amplitudes,
        "the glue channel reads the treatment block at its own offset under a "
        "wide prognostic basis");
  // the bcf-shaped reading REFUSES this layout rather than answering for it,
  // which is what keeps a K = 2, q = (1, 2) reader off a wider model
  double aOut = 0.0, b0Out = 0.0, b1Out = 0.0;
  check(!combiner.bcfGlue(aOut, b0Out, b1Out),
        "the three-scalar reading refuses a layout that is not bcf's");

  // the SURFACE route reaches the same install: what the fixture reached
  // through installForestBasis directly, a caller reaches through
  // setForestBasis, and a zero-width or non-finite basis is refused
  check(combiner.setForestBasis(0, wide.data(), 2),
        "the public route installs the same wide basis");
  check(!combiner.setForestBasis(0, wide.data(), 0),
        "a zero-width basis is refused");
  check(!combiner.setForestBasis(2, wide.data(), 1),
        "a forest index past the last is refused");
  std::vector<double> nonFinite(n, std::numeric_limits<double>::infinity());
  check(!combiner.setForestBasis(1, nonFinite.data(), 1),
        "a non-finite basis value is refused");
  combiner.restoreGlue(glue);

  std::vector<Forest<ConstantGaussianLeaf>> forests(2);
  forests[0].totalFits = mu;
  forests[1].totalFits = tau;
  const double* combined = combiner.combinedFits(forests);
  bool blend = true;
  for (size_t i = 0; i < n; ++i) {
    // the widening enters at the neutral amplitude 1.0, so the prognostic
    // multiplier is a + 1.0 * wide_i; an accessor that wrote b0 there instead
    // would move it on every row
    double m0 = a * wide[2 * i] + 1.0 * wide[2 * i + 1];
    double bz = z[i] != 0.0 ? b1 : b0;
    blend &= combined[i] == std::fma(m0, mu[i], bz * tau[i]);
  }
  check(blend,
        "a glue restore under a wide prognostic basis leaves the prognostic "
        "multiplier's own coordinates alone");

  printf("ok: amplitude offset indexing\n");
}

// The caller-supplied per-forest observation weight, on the three questions
// only the engine can answer: the refusal is a returned false and it perturbs
// nothing; the zeroed rows' sufficient statistics are invariant to their own
// responses; and the sigma posterior's degrees of freedom never see the
// channel. The end-to-end null gate, the consumer case and the surface
// refusals are the R suite's (inst/tinytest/test-forest-weights.R).
static void testForestWeights() {
  // A local stream, so adding this test does not shift the shared runif01()
  // state the hardcoded characteristic values downstream depend on.
  std::uint64_t state = 20260810u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  // ---- the zeroed rows' node statistics are invariant to their responses ----
  // With s_i = 0 on a subset S, forest 1's sumWeights and sumWeightedResponse
  // must be BITWISE unchanged when S's RESPONSES are replaced by arbitrary
  // other finite values. NOT stated as "zeroed equals dropped": the suffstat
  // kernel sums 5-wide groups after a length % 5 prologue, so physically
  // removing rows re-associates the sum and moves the last ulp. Every row here
  // is TREATED, whose multiplier is exactly one, so the exclusion is the weight
  // channel's alone - on a control row the zero multiplier would already zero
  // both quantities and the check would be vacuous. The substituted values are
  // enormous, so this doubles as the probe that a zero weight never meets an
  // amplified response and yields NaN.
  {
    const size_t m = 64;
    std::vector<double> xs(m), ys(m), yAlt(m), w(m), sMask(m), muFits(m),
      tauFits(m), zs(m, 1.0);
    for (size_t i = 0; i < m; ++i) {
      xs[i] = unif();
      ys[i] = 0.8 * unif() - 0.4;
      w[i] = 0.5 + unif();
      muFits[i] = 0.3 * unif();
      tauFits[i] = 0.2 * unif();
      sMask[i] = i % 5 == 0 ? 0.0 : 1.0;
      yAlt[i] = sMask[i] == 0.0 ? (i % 2 == 0 ? 1.0e300 : -1.0e300) : ys[i];
    }

    ColumnStore store;
    store.build(xs.data(), m, 1, 20);
    AmplitudeSpec spec;
    spec.z = zs.data();
    AmplitudeForestCombiner<ConstantGaussianLeaf> combiner(store, spec);
    std::vector<Forest<ConstantGaussianLeaf>> forests(2);
    forests[0].totalFits = muFits;
    forests[1].totalFits = tauFits;

    // the combiner's pair, composed with s exactly as the chain composes it the
    // moment formForestResponse returns; the returned pointers alias combiner
    // scratch, so each pair is copied out before the next call
    std::vector<double> yBase, wBase, ySub, wSub;
    auto formPair = [&](const double* response, std::vector<double>& outY,
                        std::vector<double>& outW) {
      ForestResponse fr =
        combiner.formForestResponse(1, forests, response, w.data());
      outY.assign(fr.response, fr.response + m);
      outW.resize(m);
      for (size_t i = 0; i < m; ++i) outW[i] = fr.weights[i] * sMask[i];
    };
    formPair(ys.data(), yBase, wBase);
    formPair(yAlt.data(), ySub, wSub);

    bool substituted = false;
    for (size_t i = 0; i < m; ++i)
      if (sMask[i] == 0.0 && std::fabs(ySub[i]) > 1.0e200) substituted = true;
    check(substituted, "the excluded rows really carry substituted responses");

    std::vector<index_t> indicesA(m), indicesB(m);
    Tree treeA, treeB;
    treeA.initialize(indicesA.data(), m);
    treeB.initialize(indicesB.data(), m);
    Rule rule;
    rule.variableIndex = 0;
    rule.setSplitIndex(static_cast<int32_t>(store.numCuts[0] / 2));
    treeA.birth(store, 0, rule, yBase.data(), wBase.data());
    treeB.birth(store, 0, rule, ySub.data(), wSub.data());
    treeA.setNodeAverages(yBase.data(), wBase.data());
    treeB.setNodeAverages(ySub.data(), wSub.data());
    treeA.bottomScratch.clear();
    treeA.fillBottom(0, treeA.bottomScratch);
    treeB.bottomScratch.clear();
    treeB.fillBottom(0, treeB.bottomScratch);

    bool invariant = treeA.bottomScratch.size() == 2 &&
                     treeB.bottomScratch.size() == 2;
    for (size_t b = 0; b < treeA.bottomScratch.size() && invariant; ++b) {
      const Node& na(treeA.at(treeA.bottomScratch[b]));
      const Node& nb(treeB.at(treeB.bottomScratch[b]));
      invariant = na.numObservations() > 0 &&
                  na.sumWeights == nb.sumWeights &&
                  na.sumWeightedResponse == nb.sumWeightedResponse &&
                  std::isfinite(nb.sumWeightedResponse);
    }
    check(invariant,
          "zero-weight rows leave the node statistics bitwise untouched");
  }

  const size_t n = 240, p = 3;
  std::vector<double> x(n * p), y(n), z(n), unitWeights(n, 1.0), zeros(n, 0.0);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    double mu = std::sin(3.0 * x[i]) + x[i + n];
    double tau = 1.0 + 2.0 * x[i + 2 * n];
    y[i] = mu + z[i] * tau + 0.2 * (unif() - 0.5);
  }

  SamplerOptions options;
  options.numTrees = 20;
  AmplitudeSpec spec;
  spec.mu.numTrees = 20;
  spec.tau.numTrees = 10;
  spec.z = z.data();

  // ---- the refusal is the engine's, and it installs and perturbs nothing ----
  // A single-forest chain handed a non-null per-forest pointer would take the
  // composed path and draw a different chain, so the refusal must land before
  // anything is stored; a multinomial carries K forests, so a forest-count
  // probe in place of the capability probe would let it through.
  ext_rng* rngRefused = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                       nullptr);
  ext_rng* rngPlain = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                     nullptr);
  ext_rng_setSeed(rngRefused, 20260810u);
  ext_rng_setSeed(rngPlain, 20260810u);
  ConstantLeafSampler refused(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rngRefused);
  ConstantLeafSampler plain(x.data(), y.data(), n, p, nullptr, nullptr,
                            ResponseFamily::gaussian, 1.0, 3.0,
                            0.37804942330213542, options, &rngPlain);
  check(!refused.supportsForestWeights(),
        "a single-forest sampler admits no per-forest weight");
  check(!refused.setForestWeights(0, unitWeights.data()),
        "a single-forest sampler refuses a per-forest weight");

  const size_t numSamples = 20;
  std::vector<double> sigmaRefused(numSamples), sigmaPlain(numSamples),
    fitsRefused(n * numSamples), fitsPlain(n * numSamples);
  Results resultsRefused, resultsPlain;
  resultsRefused.sigma = sigmaRefused.data();
  resultsRefused.trainingFits = fitsRefused.data();
  resultsPlain.sigma = sigmaPlain.data();
  resultsPlain.trainingFits = fitsPlain.data();
  refused.run(10, numSamples, resultsRefused);
  plain.run(10, numSamples, resultsPlain);
  bool unperturbed = true;
  for (size_t s = 0; s < numSamples && unperturbed; ++s)
    unperturbed = sigmaRefused[s] == sigmaPlain[s];
  for (size_t i = 0; i < n * numSamples && unperturbed; ++i)
    unperturbed = fitsRefused[i] == fitsPlain[i];
  check(unperturbed, "a refused per-forest weight leaves the draws bitwise");
  ext_rng_destroy(rngPlain);
  ext_rng_destroy(rngRefused);

  {
    const size_t K = 3;
    std::vector<int> counts(n * K, 0), trials(n, 1);
    for (size_t i = 0; i < n; ++i) counts[(i % K) * n + i] = 1;
    MultinomialSpec multinomialSpec;
    multinomialSpec.numCategories = K;
    multinomialSpec.counts = counts.data();
    multinomialSpec.trials = trials.data();
    multinomialSpec.forest.numTrees = 10;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rng, 20260811u);
    ConstantLeafSampler multinomial(x.data(), n, p, options, multinomialSpec,
                                    &rng);
    check(multinomial.numForests() == K &&
            !multinomial.supportsForestWeights(),
          "a multinomial sampler admits no per-forest weight despite K forests");
    check(!multinomial.setForestWeights(1, unitWeights.data()),
          "a multinomial sampler refuses a per-forest weight");
    ext_rng_destroy(rng);
  }

  // ---- BCF installs, refuses an out-of-range forest, and keeps the df ----
  // The sigma posterior's degrees of freedom are nu_0 + #{w_i > 0} over the
  // RESPONSE model's own precisions, so pinning those pins the df; a per-forest
  // weight lives on the chain and must never reach them.
  ext_rng* rngBCF = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngBCF, 20260812u);
  ConstantLeafSampler bcf(x.data(), y.data(), n, p, unitWeights.data(), nullptr,
                          1.0, 3.0, 0.37804942330213542, options, spec,
                          &rngBCF);
  check(bcf.supportsForestWeights(), "a BCF sampler admits a per-forest weight");
  check(!bcf.setForestWeights(2, zeros.data()),
        "a per-forest weight on a forest index out of range is refused");
  check(bcf.setForestWeights(1, zeros.data()) &&
          bcf.setForestWeights(1, nullptr) &&
          bcf.setForestWeights(1, zeros.data()),
        "a per-forest weight installs, clears, and reinstalls");

  // the df itself, not the draw it feeds: nu_0 = 3 plus the count of positive
  // OBSERVATION weights, all n of which are one here
  double dfBefore = bcf.chain(0).sigmaDegreesOfFreedomForTesting();
  std::vector<double> sigmaBCF(numSamples), fitsBCF(n * numSamples);
  Results resultsBCF;
  resultsBCF.sigma = sigmaBCF.data();
  resultsBCF.trainingFits = fitsBCF.data();
  bcf.run(10, numSamples, resultsBCF);
  check(dfBefore == 3.0 + static_cast<double>(n) &&
          bcf.chain(0).sigmaDegreesOfFreedomForTesting() == dfBefore,
        "an all-zero per-forest weight leaves the sigma posterior df at nu + n");
  bool finite = true;
  for (size_t s = 0; s < numSamples && finite; ++s) finite = sigmaBCF[s] > 0.0;
  for (size_t i = 0; i < n * numSamples && finite; ++i)
    finite = std::isfinite(fitsBCF[i]);
  check(finite, "an all-zero per-forest weight keeps the run finite");
  ext_rng_destroy(rngBCF);

  printf("ok: per-forest observation weights\n");
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
    // real single-node trees with sized leaf tables: level-centering accumulates
    // its conditional over OCCUPIED leaves, so a fixture with no trees would make
    // the move a no-op that consumes no draw
    std::vector<double> leafValue = {0.3, -0.4, 0.7};
    for (size_t k = 0; k < 3; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(n, 0.0);
      forests[k].totalFits = f[k];
      forests[k].indexBuffer.assign(n, 0);
      forests[k].trees.resize(1);
      forests[k].trees[0].initialize(forests[k].indexBuffer.data(), n);
      forests[k].muByTree.assign(1, std::vector<double>(1, leafValue[k]));
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

    // and the invariance is not vacuous: one shift c really moved every fit, and
    // the single tree of each forest absorbed all of it (c/m at m = 1)
    double shift = forests[0].totalFits[0] - f[0][0];
    bool shifted = std::fabs(shift) > 1e-8;
    for (size_t k = 0; k < 3; ++k) {
      for (size_t i = 0; i < n; ++i)
        shifted &= std::fabs(forests[k].totalFits[i] - f[k][i] - shift) < 1e-12;
      shifted &=
        std::fabs(forests[k].muByTree[0][0] - leafValue[k] - shift) < 1e-12;
    }
    check(shifted, "level-centering shifts every fit and leaf by the same c");

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
    sampler.predict(xTest.data(), nTest, 1, predicted.data());
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

// K constant-leaf forests for the level-centering seam: forest k gets
// leaves[k].size() stump trees carrying one occupied leaf apiece at the given
// value, a shared leaf scale and leaf k, and the supplied per-row total fits.
// Each tree points into its forest's own index buffer, so the forests are
// filled in place rather than returned by value.
static void buildSeamForests(std::vector<Forest<ConstantGaussianLeaf>>& forests,
                             const std::vector<std::vector<double>>& leaves,
                             const std::vector<std::vector<double>>& fits,
                             double scale, double leafK, size_t n) {
  size_t K = leaves.size();
  forests.clear();
  forests.resize(K);
  for (size_t k = 0; k < K; ++k) {
    size_t numTrees = leaves[k].size();
    forests[k].numTrees = numTrees;
    forests[k].leaf.scale = scale;
    forests[k].k = leafK;
    forests[k].totalFits = fits[k];
    forests[k].indexBuffer.assign(n * numTrees, 0);
    forests[k].trees.resize(numTrees);
    forests[k].muByTree.assign(numTrees, std::vector<double>(1, 0.0));
    for (size_t t = 0; t < numTrees; ++t) {
      forests[k].trees[t].initialize(forests[k].indexBuffer.data() + t * n, n);
      forests[k].muByTree[t][0] = leaves[k][t];
    }
  }
}

// The SECOND live override of the same post-combine virtual (combiner.hpp,
// MultinomialForestCombiner::afterCombine), pinned at the same seam as BCF's
// because M4.2 redefines the virtual for both. Its convention is the opposite
// one: the move is a dataset-wide ADDITIVE shift, drawn from a single standard
// normal, and the returned 1.0 is a constant that says nothing about whether a
// move was made - the base virtual's "the scale its move applied" is already
// false here. Dyadic fixture values, so the reference arithmetic is exact.
//
// Both conditional arms carry a reference the routine cannot supply: arm A the
// documented closed form (docs/design/multinomial.md - at the intercept-only
// configuration the move is an exact independence sampler from the level's
// marginal N(0, tau^2/K)), arm B a precision and numerator computed BY HAND off
// the fixture. Re-deriving the routine's own per-leaf sd here would pin
// nothing: the fixture sets leaf.scale directly, so an oracle written that way
// shares with the routine whatever it does with that scale.
static void testMultinomialCombinerSeam() {
  const size_t n = 4, K = 3;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);
  std::vector<int> labels = {0, 1, 2, 1};
  std::vector<int> counts, trials;
  oneHotCounts(labels, K, counts, trials);
  MultinomialSpec spec;
  spec.numCategories = K;
  spec.counts = counts.data();
  spec.trials = trials.data();
  MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);

  // Arm A, the documented identity. Every forest is a stump ensemble of m = 4
  // trees at leaf scale 1 and leaf k 2, so the per-leaf sd is s = 1/2 and the
  // per-forest total-fit sd is tau = s sqrt(m) = 1: the move must land the
  // grand level at exactly z tau/sqrt(K) for its one standard normal z,
  // WHATEVER the level was before. Each forest's totalFits is seeded to its own
  // leaf sum, which is what a stump ensemble actually fits, because the
  // conditional's mean is built from leaf sums and any other seeding reads a
  // level the identity does not describe; and the pre-move level is
  // deliberately NONZERO, so the single assertion pins the conditional's mean
  // and its sd together.
  {
    const std::vector<std::vector<double>> leaves = {{0.5, -0.25, 0.125, 0.0625},
                                                     {-0.75, 0.25, 0.5, -0.125},
                                                     {1.0, -0.5, 0.25, 0.375}};
    std::vector<std::vector<double>> fits(K);
    double preLevel = 0.0;
    for (size_t k = 0; k < K; ++k) {
      double leafSum = 0.0;
      for (double value : leaves[k]) leafSum += value;
      fits[k].assign(n, leafSum);
      preLevel += leafSum / static_cast<double>(K);
    }
    std::vector<Forest<ConstantGaussianLeaf>> forests;
    buildSeamForests(forests, leaves, fits, 1.0, 2.0, n);
    ext_rng* live = makeSeamRng();
    ext_rng* reference = makeSeamRng();
    combiner.afterCombine(forests, false, 0, live);
    double z = ext_rng_simulateStandardNormal(reference);
    double tau = 1.0;  // (leaf scale / leaf k) * sqrt(m) = (1/2) * 2
    double level = 0.0;
    for (size_t k = 0; k < K; ++k)
      level += forests[k].totalFits[0] / static_cast<double>(K);
    check(std::fabs(preLevel) > 0.1,
          "the identity arm enters at a level the move has to remove");
    checkNear(level, z * tau / std::sqrt(static_cast<double>(K)), 1e-13,
              "the intercept-only move draws the level from N(0, tau^2/K), "
              "independent of the level it found");
    ext_rng_destroy(reference);
    ext_rng_destroy(live);
  }

  // Arm B, the asymmetric conditional. Forest 0 carries TWO trees on purpose:
  // the move's per-tree factors - the conditional's m^2 and m divisors and the
  // uniform c/m absorption - are all invisible at a fixture where every m is 1,
  // and only unequal m makes the MEAN sensitive to a per-forest mis-scaling.
  // HAND-COMPUTED at leaf scale 2 and leaf k 2, so the per-leaf sd is s = 1:
  //   prec = 2/(2^2 * 1) + 1/(1 * 1) + 1/(1 * 1)         = 2.5,
  //   num  = (0.5 - 0.125)/2 + (-0.25)/1 + 0.75/1        = 0.6875,
  //   mean = -num/prec = -0.275,   sd = 1/sqrt(2.5) = 0.6324555.
  // Two seams at two seeds solve the mean and the sd out of the OUTPUT shift,
  // against those decimal targets, in a shape that shares no expression with
  // the routine; the exact reference below then pins the absorption bitwise.
  const std::vector<std::vector<double>> leaves = {{0.5, -0.125},
                                                   {-0.25},
                                                   {0.75}};
  const std::vector<std::vector<double>> fits = {{0.25, -0.5, 1.0, 0.125},
                                                 {0.5, 0.25, -0.75, 1.5},
                                                 {-1.0, 0.75, 0.5, -0.25}};
  check(leaves[0].size() > 1 && leaves[1].size() == 1,
        "the seam fixture carries forests of unequal tree count");
  const unsigned int seamSeeds[2] = {20260813u, 20260901u};
  double shift[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0};
  double returned = 0.0, c = 0.0;
  bool streamPinned = true, additive = true;
  std::vector<Forest<ConstantGaussianLeaf>> forests;
  for (size_t r = 0; r < 2; ++r) {
    buildSeamForests(forests, leaves, fits, 2.0, 2.0, n);
    ext_rng* live = makeSeamRng(seamSeeds[r]);
    ext_rng* reference = makeSeamRng(seamSeeds[r]);
    returned = combiner.afterCombine(forests, false, 0, live);
    normal[r] = ext_rng_simulateStandardNormal(reference);
    streamPinned &= sameStreamPosition(live, reference);
    shift[r] = forests[0].totalFits[0] - fits[0][0];
    c = -0.6875 / 2.5 + normal[r] / std::sqrt(2.5);
    for (size_t k = 0; k < K; ++k) {
      double perTree = c / static_cast<double>(leaves[k].size());
      for (size_t t = 0; t < leaves[k].size(); ++t)
        additive &= forests[k].muByTree[t][0] == leaves[k][t] + perTree;
      for (size_t i = 0; i < n; ++i)
        additive &= forests[k].totalFits[i] == fits[k][i] + c;
    }
    ext_rng_destroy(reference);
    ext_rng_destroy(live);
  }
  check(returned == 1.0,
        "the multinomial post-combine move returns 1.0 having moved");
  check(streamPinned,
        "the level-centering shift is one standard normal and nothing else");
  check(additive,
        "every forest's fits take the whole c while its leaves absorb c/m");
  check(std::fabs(normal[0] - normal[1]) > 0.1,
        "the two seams' standard normals are separated enough to solve on");
  double sd = (shift[0] - shift[1]) / (normal[0] - normal[1]);
  double mean = shift[0] - sd * normal[0];
  checkNear(mean, -0.275, 1e-12,
            "the asymmetric conditional centers at -num/prec");
  checkNear(sd, 0.632455532033676, 1e-12,
            "the asymmetric conditional's sd is 1/sqrt(prec)");

  // the no-move path returns the SAME 1.0 as the move above, which is exactly
  // why the return value cannot be read as this combiner's move indicator
  for (size_t k = 0; k < K; ++k) forests[k].numTrees = 0;
  ext_rng* skipRng = makeSeamRng();
  ext_rng* skipReference = makeSeamRng();
  double one = combiner.afterCombine(forests, false, 0, skipRng);
  bool inert = one == 1.0 && sameStreamPosition(skipRng, skipReference);
  for (size_t k = 0; k < K; ++k)
    for (size_t i = 0; i < n; ++i)
      inert &= forests[k].totalFits[i] == fits[k][i] + c;
  check(inert,
        "a post-combine move with no occupied leaf returns the same 1.0 and "
        "takes no draw");
  ext_rng_destroy(skipReference);
  ext_rng_destroy(skipRng);

  printf("ok: multinomial combiner seam pins\n");
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

// The counts channel on the combiner: a live multinomial combiner's count
// response is replaceable at fixed n and K. Two facts are under test. A
// combiner built over A and swapped to B must form exactly what one built over
// B forms from the same seeded rng - bitwise, since nothing about the swap is
// approximate. And the swap must carry the TRIALS as well as the successes:
// n_i drives the PG draw COUNT and the (y - n_i/2) numerator, so a swap that
// installed only the count matrix would still form the wrong response against
// the wrong-length draw stream.
//
// The burned-in arm is the one with no create-side twin: a combiner several
// glue/form cycles in carries a live omega column set and a running prefix mix,
// and no freshly built combiner can be put in that state. It is checked against
// itself (a self-swap must be inert) and against the count formula (the new
// counts and trials are both in force after the burn).
static void testMultinomialSetCounts() {
  const size_t n = 7, K = 3;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  // two count matrices over the same rows, differing in successes AND in row
  // totals, category-major (column k at k*n)
  std::vector<int> trialsA = {1, 2, 1, 3, 1, 2, 1};
  std::vector<int> trialsB = {2, 1, 3, 1, 4, 1, 2};
  std::vector<int> countsA(n * K, 0), countsB(n * K, 0);
  for (size_t i = 0; i < n; ++i) {
    countsA[(i % K) * n + i] = trialsA[i];
    countsB[((i + 1) % K) * n + i] = trialsB[i];
  }

  // static fits: the combiner is driven directly, so no tree update moves them
  std::vector<Forest<ConstantGaussianLeaf>> forests(K);
  std::vector<std::vector<double>> f = {
    {0.3, -0.4, 0.8, 0.1, -0.6, 0.5, 0.2},
    {-0.2, 0.6, -0.3, 0.4, 0.7, -0.5, 0.0},
    {0.5, 0.1, -0.7, 0.2, -0.1, 0.6, -0.4}};
  for (size_t k = 0; k < K; ++k) {
    forests[k].numTrees = 1;
    forests[k].leaf.scale = 3.0;
    forests[k].k = 2.0;
    forests[k].treeFits.assign(n, 0.0);
    forests[k].totalFits = f[k];
  }

  MultinomialSpec specA, specB;
  specA.numCategories = specB.numCategories = K;
  specA.counts = countsA.data();
  specA.trials = trialsA.data();
  specB.counts = countsB.data();
  specB.trials = trialsB.data();

  // one pass of the sweep's per-forest pair, recording every response and
  // weight when asked; sweeps of them replay the interleaving
  auto cycle = [&](MultinomialForestCombiner<ConstantGaussianLeaf>& combiner,
                   ext_rng* rng, size_t sweeps, std::vector<double>* out) {
    for (size_t s = 0; s < sweeps; ++s)
      for (size_t fk = 0; fk < K; ++fk) {
        combiner.drawForestGlue(fk, rng, forests);
        ForestResponse fr =
          combiner.formForestResponse(fk, forests, nullptr, nullptr);
        if (out == NULL) continue;
        for (size_t i = 0; i < n; ++i) {
          out->push_back(fr.response[i]);
          out->push_back(fr.weights[i]);
        }
      }
  };
  auto seeded = [](std::uint32_t seed) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    return rng;
  };

  // (1) swap-vs-build parity, and its non-vacuity arm
  std::vector<double> recSwap, recBuild, recA;
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specA);
    check(combiner.supportsCountsMutation(),
          "the multinomial combiner owns a replaceable count response");
    combiner.setCounts(countsB.data(), trialsB.data());
    ext_rng* rng = seeded(9001u);
    cycle(combiner, rng, 2, &recSwap);
    ext_rng_destroy(rng);
  }
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specB);
    ext_rng* rng = seeded(9001u);
    cycle(combiner, rng, 2, &recBuild);
    ext_rng_destroy(rng);
  }
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specA);
    ext_rng* rng = seeded(9001u);
    cycle(combiner, rng, 2, &recA);
    ext_rng_destroy(rng);
  }
  check(recSwap.size() == recBuild.size() && recSwap == recBuild,
        "a swapped combiner forms what one built over the new counts forms");
  check(recA != recBuild,
        "and the two count matrices are not the same problem");

  // (2) burned-in: a self-swap is inert, and a swap after the burn installs
  //     both halves of the response
  std::vector<double> recSelf, recControl;
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specA);
    ext_rng* rng = seeded(9002u);
    cycle(combiner, rng, 3, NULL);
    combiner.setCounts(countsA.data(), trialsA.data());
    cycle(combiner, rng, 1, &recSelf);
    ext_rng_destroy(rng);
  }
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specA);
    ext_rng* rng = seeded(9002u);
    cycle(combiner, rng, 3, NULL);
    cycle(combiner, rng, 1, &recControl);
    ext_rng_destroy(rng);
  }
  check(recSelf == recControl,
        "a self-swap on a burned-in combiner changes nothing");

  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, specA);
    ext_rng* rng = seeded(9003u);
    cycle(combiner, rng, 3, NULL);
    combiner.setCounts(countsB.data(), trialsB.data());
    bool formed = true;
    for (size_t fk = 0; fk < K; ++fk) {
      combiner.drawForestGlue(fk, rng, forests);
      ForestResponse fr =
        combiner.formForestResponse(fk, forests, nullptr, nullptr);
      for (size_t i = 0; i < n; ++i) {
        // the margin over the other categories, this fixture's fits being
        // static; the response must be B's (y - n_i/2)/omega + C
        double margin = -HUGE_VAL;
        for (size_t j = 0; j < K; ++j) {
          if (j == fk) continue;
          margin = margin == -HUGE_VAL
            ? f[j][i]
            : std::max(margin, f[j][i]) +
                std::log1p(std::exp(-std::fabs(f[j][i] - margin)));
        }
        double expected =
          (static_cast<double>(countsB[fk * n + i]) - trialsB[i] * 0.5) /
            fr.weights[i] + margin;
        formed &= std::fabs(fr.response[i] - expected) < 1e-12;
      }
    }
    check(formed,
          "a burned-in swap forms the new counts against the new trials");
    ext_rng_destroy(rng);
  }

  printf("ok: multinomial counts channel\n");
}

// The category offset on the combiner: the latent is f_ik + o_ik, so the offset
// enters the margin, is SUBTRACTED back out of the working response, and rides
// the reported blend. Four facts are under test, each of which a plausible
// wrong implementation would break.
//
// The formed response is checked against the formula directly, so a sign or a
// placement error shows here rather than as a shifted posterior. An all-zero
// offset must reproduce the offset-free stream bitwise, and so must clearing
// one, which is what makes the null path structurally today's path rather than
// today's path plus a rounding argument. And the blend must be same-vintage:
// fits written between sweeps - which is what a predictor mutation's tree
// revalidation does, outside the combiner entirely - must appear in the very
// next reported blend, not one sweep later.
static void testMultinomialCategoryOffset() {
  const size_t n = 6, K = 3;
  std::vector<double> xDummy(n, 0.0);
  ColumnStore data;
  data.build(xDummy.data(), n, 1, 100);

  std::vector<int> trials = {1, 2, 1, 3, 2, 1};
  std::vector<int> counts(n * K, 0);
  for (size_t i = 0; i < n; ++i) counts[(i % K) * n + i] = trials[i];

  std::vector<Forest<ConstantGaussianLeaf>> forests(K);
  std::vector<std::vector<double>> f = {
    {0.3, -0.4, 0.8, 0.1, -0.6, 0.5},
    {-0.2, 0.6, -0.3, 0.4, 0.7, -0.5},
    {0.5, 0.1, -0.7, 0.2, -0.1, 0.6}};
  for (size_t k = 0; k < K; ++k) {
    forests[k].numTrees = 1;
    forests[k].leaf.scale = 3.0;
    forests[k].k = 2.0;
    forests[k].treeFits.assign(n, 0.0);
    forests[k].totalFits = f[k];
  }

  // a per-category, per-observation offset with no common row component, so
  // nothing about it lies along the softmax's null direction
  std::vector<double> offset(n * K), zeros(n * K, 0.0);
  for (size_t k = 0; k < K; ++k)
    for (size_t i = 0; i < n; ++i)
      offset[k * n + i] = 0.4 * static_cast<double>(k) -
                          0.15 * static_cast<double>(i % 4);

  MultinomialSpec spec;
  spec.numCategories = K;
  spec.counts = counts.data();
  spec.trials = trials.data();

  auto seeded = [](std::uint32_t seed) {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    return rng;
  };
  auto cycle = [&](MultinomialForestCombiner<ConstantGaussianLeaf>& combiner,
                   ext_rng* rng, size_t sweeps, std::vector<double>* out) {
    for (size_t s = 0; s < sweeps; ++s)
      for (size_t fk = 0; fk < K; ++fk) {
        combiner.drawForestGlue(fk, rng, forests);
        ForestResponse fr =
          combiner.formForestResponse(fk, forests, nullptr, nullptr);
        if (out == NULL) continue;
        for (size_t i = 0; i < n; ++i) {
          out->push_back(fr.response[i]);
          out->push_back(fr.weights[i]);
        }
      }
  };

  // (1) the formula: the margin is the log-sum-exp over the OTHER categories'
  //     offset fits, and the response carries - o_if
  {
    MultinomialSpec withOffset = spec;
    withOffset.offset = offset.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, withOffset);
    ext_rng* rng = seeded(9101u);
    bool formed = true;
    for (size_t fk = 0; fk < K; ++fk) {
      combiner.drawForestGlue(fk, rng, forests);
      ForestResponse fr =
        combiner.formForestResponse(fk, forests, nullptr, nullptr);
      for (size_t i = 0; i < n; ++i) {
        double margin = -HUGE_VAL;
        for (size_t j = 0; j < K; ++j) {
          if (j == fk) continue;
          double raw = f[j][i] + offset[j * n + i];
          margin = margin == -HUGE_VAL
            ? raw
            : std::max(margin, raw) +
                std::log1p(std::exp(-std::fabs(raw - margin)));
        }
        double expected =
          (static_cast<double>(counts[fk * n + i]) - trials[i] * 0.5) /
            fr.weights[i] + margin - offset[fk * n + i];
        formed &= std::fabs(fr.response[i] - expected) < 1e-12;
      }
    }
    check(formed,
          "an offset combiner forms (y - n/2)/omega + C - o against the "
          "offset margin");
    ext_rng_destroy(rng);
  }

  // (2) an all-zero offset is the null path, bitwise, and so is a cleared one;
  //     an installed offset is not (the non-vacuity arm)
  std::vector<double> recNull, recZero, recCleared, recOffset, recSwapped;
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    ext_rng* rng = seeded(9102u);
    cycle(combiner, rng, 2, &recNull);
    ext_rng_destroy(rng);
  }
  {
    MultinomialSpec zeroed = spec;
    zeroed.offset = zeros.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, zeroed);
    ext_rng* rng = seeded(9102u);
    cycle(combiner, rng, 2, &recZero);
    ext_rng_destroy(rng);
  }
  {
    MultinomialSpec withOffset = spec;
    withOffset.offset = offset.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, withOffset);
    ext_rng* rng = seeded(9102u);
    cycle(combiner, rng, 2, &recOffset);
    ext_rng_destroy(rng);
  }
  {
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    combiner.setCategoryOffset(offset.data());
    ext_rng* rng = seeded(9102u);
    cycle(combiner, rng, 2, &recSwapped);
    ext_rng_destroy(rng);
  }
  {
    MultinomialSpec withOffset = spec;
    withOffset.offset = offset.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, withOffset);
    combiner.setCategoryOffset(nullptr);
    ext_rng* rng = seeded(9102u);
    cycle(combiner, rng, 2, &recCleared);
    ext_rng_destroy(rng);
  }
  check(recZero == recNull, "an all-zero category offset is the null path");
  check(recCleared == recNull, "clearing a category offset restores it");
  check(recSwapped == recOffset,
        "installing an offset matches building with it");
  check(recOffset != recNull, "and the offset is not inert");

  // (3) the reported blend is softmax(f + o) and is SAME-VINTAGE: a totalFits
  //     write between sweeps - what a predictor mutation's revalidation does,
  //     never entering the combiner - must show in the next blend
  {
    MultinomialSpec withOffset = spec;
    withOffset.offset = offset.data();
    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, withOffset);
    ext_rng* rng = seeded(9103u);
    cycle(combiner, rng, 1, NULL);
    for (size_t i = 0; i < n; ++i) forests[K - 1].totalFits[i] += 0.75;
    const double* blended = combiner.combinedFits(forests);
    bool current = true;
    for (size_t i = 0; i < n; ++i) {
      double maxRaw = -HUGE_VAL, sumExp = 0.0;
      for (size_t k = 0; k < K; ++k)
        maxRaw = std::max(maxRaw, forests[k].totalFits[i] + offset[k * n + i]);
      for (size_t k = 0; k < K; ++k)
        sumExp +=
          std::exp(forests[k].totalFits[i] + offset[k * n + i] - maxRaw);
      for (size_t k = 0; k < K; ++k) {
        double p =
          std::exp(forests[k].totalFits[i] + offset[k * n + i] - maxRaw) /
          sumExp;
        current &= std::fabs(blended[k * n + i] - p) < 1e-14;
      }
    }
    check(current,
          "the reported blend is the softmax of every category's CURRENT "
          "offset fit");
    for (size_t i = 0; i < n; ++i) forests[K - 1].totalFits[i] -= 0.75;
    ext_rng_destroy(rng);
  }

  // (4) the test blend's own offset: a SEPARATE nTest x K object over the test
  //     rows, entering the reported test softmax where the train offset enters
  //     the train one. Its slab is sized off the test store, so this is also
  //     the arm that walks that buffer under the sanitizers.
  {
    const size_t nTest = 5;
    std::vector<double> xTestDummy(nTest, 0.0);
    data.buildTest(xTestDummy.data(), nTest);
    std::vector<double> testOffset(nTest * K), testZeros(nTest * K, 0.0);
    for (size_t k = 0; k < K; ++k)
      for (size_t i = 0; i < nTest; ++i) {
        forests[k].totalTestFits.resize(nTest);
        forests[k].totalTestFits[i] =
          0.2 * static_cast<double>(k) - 0.35 * static_cast<double>(i % 3);
        testOffset[k * nTest + i] =
          0.5 - 0.3 * static_cast<double>(k) + 0.1 * static_cast<double>(i);
      }

    MultinomialForestCombiner<ConstantGaussianLeaf> combiner(data, spec);
    const double* plain = combiner.combinedTestFits(forests);
    std::vector<double> recPlain(plain, plain + nTest * K);
    combiner.setCategoryTestOffset(testZeros.data());
    const double* zeroed = combiner.combinedTestFits(forests);
    bool zeroIsNull = true;
    for (size_t j = 0; j < nTest * K; ++j)
      zeroIsNull &= zeroed[j] == recPlain[j];
    check(zeroIsNull, "an all-zero category TEST offset is the null path");

    combiner.setCategoryTestOffset(testOffset.data());
    const double* blended = combiner.combinedTestFits(forests);
    bool formed = true, moved = false;
    for (size_t i = 0; i < nTest; ++i) {
      double maxRaw = -HUGE_VAL, sumExp = 0.0;
      for (size_t k = 0; k < K; ++k)
        maxRaw = std::max(maxRaw, forests[k].totalTestFits[i] +
                                    testOffset[k * nTest + i]);
      for (size_t k = 0; k < K; ++k)
        sumExp += std::exp(forests[k].totalTestFits[i] +
                           testOffset[k * nTest + i] - maxRaw);
      for (size_t k = 0; k < K; ++k) {
        double p = std::exp(forests[k].totalTestFits[i] +
                            testOffset[k * nTest + i] - maxRaw) / sumExp;
        formed &= std::fabs(blended[k * nTest + i] - p) < 1e-14;
        moved |= blended[k * nTest + i] != recPlain[k * nTest + i];
      }
    }
    check(formed, "the test blend is the softmax of the offset TEST fits");
    check(moved, "and the test offset is not inert");

    // the test blend follows the test fits, not a snapshot: a totalTestFits
    // write between reports shows in the next one
    for (size_t i = 0; i < nTest; ++i) forests[K - 1].totalTestFits[i] += 0.9;
    const double* after = combiner.combinedTestFits(forests);
    bool current = true;
    for (size_t i = 0; i < nTest; ++i) {
      double maxRaw = -HUGE_VAL, sumExp = 0.0;
      for (size_t k = 0; k < K; ++k)
        maxRaw = std::max(maxRaw, forests[k].totalTestFits[i] +
                                    testOffset[k * nTest + i]);
      for (size_t k = 0; k < K; ++k)
        sumExp += std::exp(forests[k].totalTestFits[i] +
                           testOffset[k * nTest + i] - maxRaw);
      for (size_t k = 0; k < K; ++k)
        current &= std::fabs(after[k * nTest + i] -
                             std::exp(forests[k].totalTestFits[i] +
                                      testOffset[k * nTest + i] - maxRaw) /
                               sumExp) < 1e-14;
    }
    check(current, "the test blend reads every category's CURRENT test fit");

    // clearing returns to the offset-free test path, pointer for pointer
    for (size_t i = 0; i < nTest; ++i) forests[K - 1].totalTestFits[i] -= 0.9;
    combiner.setCategoryTestOffset(nullptr);
    const double* cleared = combiner.combinedTestFits(forests);
    bool restored = true;
    for (size_t j = 0; j < nTest * K; ++j)
      restored &= cleared[j] == recPlain[j];
    check(restored, "clearing a category test offset restores the null path");

    // the TRAIN offset is untouched by any of it
    check(combiner.combinedFits(forests) != nullptr,
          "the train blend survives a test-offset cycle");
  }

  printf("ok: multinomial category offset\n");
}

// The active-row mask on the softmax coupling, at the kernel. The mask is
// GLOBAL, so a masked combiner over n rows is the SAME combiner over the
// compacted active rows: bit for bit in the working response and the composed
// precision at every active row, and variate for variate in the Polya-Gamma
// stream. The comparison is well posed because the margin C_if is row-local (a
// log-sum-exp across the categories AT THAT ROW), so deleting a row changes no
// other row's conditional.
//
// This is the arm that pins the SKIP semantics. The Polya-Gamma sampler is a
// rejection sampler whose consumption depends on psi, so an implementation that
// draws an inactive row's K latents and discards them desynchronizes the two
// streams and every later active row parts. Two further pins ride along: an
// inactive row's composed precision is EXACTLY zero (it leaves every leaf
// sufficient statistic, branch score and leaf draw of every forest), and its
// working response stays finite, because the zero goes into the precision and
// never into omega, which the response divides by.
static void testActiveRowsMultinomialKernel() {
  const size_t n = 12, K = 3;
  std::vector<double> active(n);
  std::vector<size_t> kept;
  for (size_t i = 0; i < n; ++i) {
    active[i] = (i == 0 || i % 5 == 1) ? 0.0 : 1.0;
    if (active[i] != 0.0) kept.push_back(i);
  }
  const size_t m = kept.size();

  // multi-trial rows, so the arm also covers PG(n_i, .)'s variate summing: a
  // skipped row must consume none of its n_i draws
  std::vector<int> counts(n * K, 0), trials(n);
  std::vector<double> fits(n * K);
  for (size_t i = 0; i < n; ++i) {
    trials[i] = 1 + static_cast<int>(i % 3);
    counts[((3 * i + 1) % K) * n + i] = trials[i];
    for (size_t k = 0; k < K; ++k)
      fits[k * n + i] =
        0.4 * static_cast<double>(i) - 1.1 * static_cast<double>(k);
  }
  std::vector<int> countsKept(m * K, 0), trialsKept(m);
  std::vector<double> fitsKept(m * K);
  for (size_t j = 0; j < m; ++j) {
    trialsKept[j] = trials[kept[j]];
    for (size_t k = 0; k < K; ++k) {
      countsKept[k * m + j] = counts[k * n + kept[j]];
      fitsKept[k * m + j] = fits[k * n + kept[j]];
    }
  }

  auto build = [](size_t rows, const int* c, const int* t, const double* f,
                  size_t categories, ColumnStore& data,
                  std::vector<double>& xDummy,
                  std::vector<Forest<ConstantGaussianLeaf>>& forests) {
    xDummy.assign(rows, 0.0);
    data.build(xDummy.data(), rows, 1, 100);
    forests.resize(categories);
    for (size_t k = 0; k < categories; ++k) {
      forests[k].numTrees = 1;
      forests[k].leaf.scale = 3.0;
      forests[k].k = 2.0;
      forests[k].treeFits.assign(rows, 0.0);
      forests[k].totalFits.assign(f + k * rows, f + (k + 1) * rows);
    }
    MultinomialSpec spec;
    spec.numCategories = categories;
    spec.counts = c;
    spec.trials = t;
    return spec;
  };

  ColumnStore dataMasked, dataKept;
  std::vector<double> xMasked, xKept;
  std::vector<Forest<ConstantGaussianLeaf>> forestsMasked, forestsKept;
  MultinomialSpec specMasked = build(n, counts.data(), trials.data(),
                                     fits.data(), K, dataMasked, xMasked,
                                     forestsMasked);
  MultinomialSpec specKept = build(m, countsKept.data(), trialsKept.data(),
                                   fitsKept.data(), K, dataKept, xKept,
                                   forestsKept);
  MultinomialForestCombiner<ConstantGaussianLeaf> masked(dataMasked,
                                                         specMasked);
  MultinomialForestCombiner<ConstantGaussianLeaf> compacted(dataKept, specKept);
  masked.setActiveRows(active.data());

  ext_rng* rngMasked = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngCompacted = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                         NULL);
  ext_rng_setSeed(rngMasked, 20260813u);
  ext_rng_setSeed(rngCompacted, 20260813u);

  bool agrees = true, zeroed = true, finite = true;
  for (size_t f = 0; f < K; ++f) {
    masked.drawForestGlue(f, rngMasked, forestsMasked);
    ForestResponse mf = masked.formForestResponse(f, forestsMasked, nullptr,
                                                  nullptr);
    compacted.drawForestGlue(f, rngCompacted, forestsKept);
    ForestResponse cf = compacted.formForestResponse(f, forestsKept, nullptr,
                                                     nullptr);
    for (size_t j = 0; j < m; ++j)
      agrees = agrees && mf.response[kept[j]] == cf.response[j] &&
               mf.weights[kept[j]] == cf.weights[j];
    for (size_t i = 0; i < n; ++i)
      if (active[i] == 0.0) {
        zeroed = zeroed && mf.weights[i] == 0.0;
        finite = finite && std::isfinite(mf.response[i]);
      }
  }
  check(agrees,
        "a masked multinomial coupling is bitwise the compacted one at every "
        "active row");
  check(zeroed, "an inactive row's composed precision is exactly zero");
  check(finite,
        "an inactive row's working response stays finite (omega is never "
        "zeroed)");
  check(ext_rng_simulateContinuousUniform(rngMasked) ==
          ext_rng_simulateContinuousUniform(rngCompacted),
        "an inactive row's K Polya-Gamma draws are skipped, not discarded");

  // clearing restores the unmasked kernel exactly, which is what makes the
  // engine's all-ones normalization a no-op rather than a near-no-op
  masked.setActiveRows(nullptr);
  ext_rng* rngCleared = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                       NULL);
  ext_rng* rngNever = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngCleared, 20260814u);
  ext_rng_setSeed(rngNever, 20260814u);
  ColumnStore dataNever;
  std::vector<double> xNever;
  std::vector<Forest<ConstantGaussianLeaf>> forestsNever;
  MultinomialSpec specNever = build(n, counts.data(), trials.data(),
                                    fits.data(), K, dataNever, xNever,
                                    forestsNever);
  MultinomialForestCombiner<ConstantGaussianLeaf> never(dataNever, specNever);
  bool restored = true;
  for (size_t f = 0; f < K; ++f) {
    masked.drawForestGlue(f, rngCleared, forestsMasked);
    ForestResponse mf = masked.formForestResponse(f, forestsMasked, nullptr,
                                                  nullptr);
    never.drawForestGlue(f, rngNever, forestsNever);
    ForestResponse nf = never.formForestResponse(f, forestsNever, nullptr,
                                                 nullptr);
    for (size_t i = 0; i < n; ++i)
      restored = restored && mf.response[i] == nf.response[i] &&
                 mf.weights[i] == nf.weights[i];
  }
  check(restored, "clearing the mask restores the unmasked coupling bitwise");

  ext_rng_destroy(rngNever);
  ext_rng_destroy(rngCleared);
  ext_rng_destroy(rngCompacted);
  ext_rng_destroy(rngMasked);
  printf("ok: multinomial active-row mask kernel\n");
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
  AmplitudeSpec spec;
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

// The info dump names its forest. The widened printers must reach forest 1 of
// a BCF sampler on BOTH branches - saved slots under keepTrees, live trees
// without it - while forest 0 renders exactly what it always did. tau is
// restricted to one moderator column, so the forest a dump came from is
// legible in the dump itself: every split tau prints names column 0, and mu
// names others.
static void testPrintTreesForest() {
  // A local stream and locally owned generator, so adding this test shifts
  // neither the shared runif01() state nor the shared rng the hardcoded
  // characteristic values downstream depend on.
  std::uint64_t state = 20260814u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  const size_t n = 300, p = 4;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    y[i] = 3.0 * x[i + 2 * n] + x[i + n] + z[i] * 4.0 * x[i] +
           0.2 * (unif() - 0.5);
  }

  const std::vector<size_t> moderators = {0};
  SamplerOptions options;
  options.keepTrees = true;
  options.numSamplesToStore = 2;
  AmplitudeSpec spec;
  spec.mu.numTrees = 30;
  spec.tau.numTrees = 20;
  spec.tau.columns = moderators.data();
  spec.tau.numColumns = moderators.size();
  spec.z = z.data();

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 9310u);
  Sampler<ConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0,
    0.37804942330213542, options, spec, &rng);
  Results empty;
  sampler.run(150, 2, empty);

  const size_t chain = 0, slot = 0, trees[] = {0, 1, 2};
  auto dump = [&](size_t forestIndex, size_t numSampleIndices) {
    std::string text;
    beginPrintCapture(text);
    sampler.printTrees(&chain, 1, &slot, numSampleIndices, trees, 3,
                       forestIndex, false);
    endPrintCapture();
    return text;
  };
  auto splitVariables = [](const std::string& text) {
    std::vector<long> variables;
    for (size_t at = text.find("var: "); at != std::string::npos;
         at = text.find("var: ", at + 1))
      variables.push_back(std::strtol(text.c_str() + at + 5, nullptr, 10));
    return variables;
  };

  for (int saved = 1; saved >= 0; --saved) {
    if (saved == 0) sampler.setTreeStorage(false, 0);
    std::string branch(saved ? "saved " : "live ");
    std::string mu = dump(0, saved ? 1 : 0), tau = dump(1, saved ? 1 : 0);
    // forest 0's head is the rendering the R surface's format tests pin
    check(mu.compare(0, saved ? 6 : 12, saved ? "TBN: 1" : "n: 300 TBN: ") == 0,
          (branch + "dump of forest 0 opens in its branch's format").c_str());
    check(!tau.empty() && mu != tau,
          (branch + "dump of forest 1 is its own").c_str());
    std::vector<long> muVariables = splitVariables(mu),
                      tauVariables = splitVariables(tau);
    bool tauContained = !tauVariables.empty();
    for (long variable : tauVariables) tauContained &= variable == 0;
    check(tauContained, (branch + "dump of forest 1 is tau's, splitting only "
                                  "on its moderator").c_str());
    bool muBeyond = false;
    for (long variable : muVariables) muBeyond |= variable != 0;
    check(muBeyond,
          (branch + "dump of forest 0 is mu's, reading the full store").c_str());
  }

  ext_rng_destroy(rng);
  printf("ok: printTrees forest index\n");
}

// The mid-chain calibration surface at the engine boundary. Four claims the R
// tests can only see through the bridge: the reported prior scale is the leaf
// scale carried into response units by an independently computed factor; the
// write lands on EVERY chain; a read-then-write is bitwise inert, in the
// internal scale and not merely in the reported one; and a combiner refuses
// the write while the reader still serves each of its forests.
static void testForestCalibration() {
  // A local stream and locally owned generators, so adding this test shifts
  // neither the shared runif01() state nor the shared rng the hardcoded
  // characteristic values downstream depend on.
  std::uint64_t state = 20260813u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  const size_t n = 300, p = 3, numChains = 3, numTrees = 40;
  std::vector<double> x(n * p), y(n), z(n);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    z[i] = unif() < 0.5 ? 1.0 : 0.0;
    y[i] = 7.0 * (x[i] - x[i + n]) + z[i] + 0.4 * (unif() - 0.5);
  }

  ext_rng* rngs[numChains];
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 9200 + static_cast<std::uint32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = numTrees;
  options.numChains = numChains;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, rngs);

  // the reported scale against the factor recomputed from the response
  // transform and the tree count alone
  double factor = sampler.fitScale() * std::sqrt(static_cast<double>(numTrees));
  ForestCalibration calibration = sampler.forestCalibration(0, 0);
  check(calibration.responseScale == sampler.fitScale(),
        "calibration: the reported transform is the response's own");
  check(calibration.priorScale ==
          sampler.chain(0).leaf().scale * factor,
        "calibration: prior scale is the leaf scale in response units");
  check(calibration.priorSd == calibration.priorScale / calibration.k,
        "calibration: prior sd is prior scale over k");
  check(calibration.priorMean == calibration.responseShift,
        "calibration: prior mean is the transform's shift");
  check(!calibration.kHasHyperprior,
        "calibration: a fixed k reports no hyperprior");
  // the default is the family-keyed node scale carried out, so the surface is
  // reading the engine rather than echoing a stored name
  check(std::fabs(calibration.priorScale / (0.5 * sampler.fitScale()) - 1.0) <
          1.0e-12,
        "calibration: an unnamed model reports the node scale in response "
        "units");

  // the write lands on every chain, and a read-then-write does not touch the
  // internal scale on any of them
  check(sampler.setForestPriorScale(0, 2.5),
        "calibration: a single-forest write is accepted");
  std::vector<double> written(numChains);
  for (size_t c = 0; c < numChains; ++c) {
    written[c] = sampler.chain(c).leaf().scale;
    check(std::fabs(sampler.forestCalibration(c, 0).priorScale / 2.5 - 1.0) <
            8.0 * DBL_EPSILON,
          "calibration: the write reaches every chain");
  }
  for (size_t c = 0; c < numChains; ++c)
    sampler.setForestPriorScale(0, sampler.forestCalibration(c, 0).priorScale);
  for (size_t c = 0; c < numChains; ++c)
    check(sampler.chain(c).leaf().scale == written[c],
          "calibration: a read-then-write leaves the internal scale bitwise "
          "untouched");
  // and the skip is not vacuous - a value that is not what is in force writes
  check(sampler.setForestPriorScale(0, 2.5 * (1.0 + 1.0e-9)),
        "calibration: a different value is accepted");
  check(sampler.chain(0).leaf().scale != written[0],
        "calibration: a different value really moves the internal scale");

  // an out-of-range forest is a capability answer, not a raise
  check(!sampler.setForestPriorScale(1, 2.5),
        "calibration: an out-of-range forest refuses the write");
  // and the READER, which has no refusal channel of its own, answers with a
  // default-constructed calibration rather than reading past the last forest -
  // the flat C getter turns that into its own "0, touching nothing"
  ForestCalibration missing = sampler.forestCalibration(0, 1);
  check(missing.priorScale == 0.0 && missing.priorSd == 0.0 &&
          missing.priorMean == 0.0 && missing.k == 0.0 &&
          missing.responseScale == 0.0 && missing.responseShift == 0.0 &&
          !missing.kHasHyperprior,
        "calibration: an out-of-range forest reads as the empty calibration");
  // and the five map quantities are NaN there and on the forest that DOES
  // exist: a single-forest sampler's scale is not map-derived, which is a
  // positive signal rather than a plausible 1.0 a caller would multiply by
  check(std::isnan(missing.amplitudePriorVariance) &&
          std::isnan(missing.amplitudePriorScale) &&
          std::isnan(missing.nodeScaleFactor) &&
          std::isnan(missing.nodeScaleDivisor) &&
          std::isnan(missing.basisRowNorm),
        "calibration: the empty calibration carries no map quantities");
  ForestCalibration offMap = sampler.forestCalibration(0, 0);
  check(std::isnan(offMap.amplitudePriorVariance) &&
          std::isnan(offMap.amplitudePriorScale) &&
          std::isnan(offMap.nodeScaleFactor) &&
          std::isnan(offMap.nodeScaleDivisor) &&
          std::isnan(offMap.basisRowNorm),
        "calibration: a single-forest sampler reports no map quantities");

  {
    ext_rng* bcfRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(bcfRng, 9210);
    SamplerOptions bcfOptions;
    AmplitudeSpec spec;
    spec.mu.numTrees = 30;
    spec.tau.numTrees = 15;
    spec.z = z.data();
    Sampler<ConstantGaussianLeaf> bcf(x.data(), y.data(), n, p, nullptr,
                                      nullptr, 1.0, 3.0, 0.37804942330213542,
                                      bcfOptions, spec, &bcfRng);
    check(bcf.numForests() == 2, "calibration: the BCF fixture has two forests");
    // the reader is total over a combiner's forests, each against its own tree
    // count, so the map's ratio is visible rather than hidden
    ForestCalibration mu = bcf.forestCalibration(0, 0);
    ForestCalibration tau = bcf.forestCalibration(0, 1);
    check(mu.priorScale > 0.0 && tau.priorScale > 0.0,
          "calibration: both BCF forests report a positive prior scale");
    check(mu.k == 1.0 && !mu.kHasHyperprior,
          "calibration: the BCF map pins k at one");
    check(mu.priorScale ==
            bcf.chain(0).leaf().scale * bcf.fitScale() * std::sqrt(30.0),
          "calibration: the BCF prognostic scale carries its own tree count");
    // the writer refuses, because the map owns both halves
    check(!bcf.setForestPriorScale(0, 2.5) &&
            !bcf.setForestPriorScale(1, 2.5),
          "calibration: a combiner refuses the write on every forest");
    ext_rng_destroy(bcfRng);
  }

  for (size_t c = 0; c < numChains; ++c) ext_rng_destroy(rngs[c]);
  printf("ok: forest calibration read/write (prior scale %.4f)\n",
         calibration.priorScale);
}

/// The calibration map's PRODUCT, which nothing in tests/cpp pinned: a forest's
/// node scale is nodeScaleFactor * s / (nodeScaleDivisor * basisRowNorm), and
/// under a latent family forestCalibration reports it exactly - fitScale is 1
/// and priorScale's sqrt(m) cancels the map's, so the reported number IS the
/// map's node scale rather than a transform of it.
///
/// The reader now reports the FACTORS as well as the product
/// (docs/plans/binary-kforest-prior-default.md, leg (c)), so each is pinned
/// against its literal here and the three together against the scale they
/// decompose - a column filled from a neighbouring forest's buffer, or from
/// the other member of the exclusive amplitude pair, cannot pass both.
///
/// Fixture admissibility, in the comment because this is the third set of
/// values: the first left two of the four factors at 1 and the second gave a
/// forest a denominator whose PRODUCT was 1 and a prior scale coincidentally
/// equal to its own factor, so both vacated the pin they were built to close.
/// All nine values are distinct and none is 1 - factors {2.5, 0.4, 1.75},
/// divisors {0.674, 0.8, 0.25}, row norms {3, 5, 6}. The three denominators
/// divisor * norm = {2.022, 4, 1.5} are distinct, none 1, and none equal to any
/// of the nine. The prior scales they induce - probit {1.236400, 0.1,
/// 1.166667}, logistic {2.242581, 0.181380, 2.116099} - are distinct from each
/// other, from s, and each from its own forest's factor, divisor and norm. Tree
/// counts differ per forest, so the sqrt(m) cancellation is read per forest
/// rather than off forest 0's.
///
/// BOTH latent families, because probit ALONE CANNOT SEE A DROPPED ANCHOR:
/// s = 1 there, so factor * s / (divisor * norm) and factor / (divisor * norm)
/// are the same number. Under logistic they differ by 81 percent. That is the
/// near-miss shape M4.4 recorded when a naive anchor missed its logistic arm by
/// 0.865 percent, and it is why this loop is not probit-only.
static void testBCFCalibrationMap() {
  const size_t n = 24, p = 2;
  std::vector<double> x(n * p), y(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) y[i] = i % 2 == 0 ? 1.0 : 0.0;

  // ROW-MAJOR per-forest bases at row norms 3, 5 and 6. Forest 3's rows are a
  // genuine two-column (3.6, 4.8), so a norm formed from one column alone, or
  // from the absolute values, would report 3.6 or 8.4 rather than 6.
  std::vector<double> b0(n, 3.0), b1(2 * n, 0.0), b2(2 * n);
  for (size_t i = 0; i < n; ++i) {
    b1[2 * i + i % 2] = 5.0;
    b2[2 * i] = 3.6;
    b2[2 * i + 1] = 4.8;
  }
  // the row-norm CONVENTION, which no test anywhere exercised: 22 nonzero rows
  // at norms 1..22 (an EVEN count, so the average-of-two-central-order-
  // statistics branch runs, giving 11.5) and two ALL-ZERO rows, excluded rather
  // than counted small. Counting the zeros would give 10.5, and returning the
  // upper central order statistic alone 12. And an ALL-ZERO basis, whose norm
  // falls back to 1 so the forest still reports a finite scale.
  std::vector<double> even(n, 0.0), zeros(n, 0.0);
  for (size_t i = 0; i + 2 < n; ++i) even[i] = static_cast<double>(i + 1);

  const double factors[3] = {2.5, 0.4, 1.75};
  const double divisors[3] = {0.674, 0.8, 0.25};
  const double norms[3] = {3.0, 5.0, 6.0};
  // the amplitude prior each forest carries, which the reader echoes beside
  // the map's own two factors. EXCLUSIVE per forest: forest 0 takes the
  // half-Cauchy scale mixture and reports a NaN variance, forests 1 and 2 the
  // fixed variance and a NaN scale. All three values are distinct from each
  // other and from the nine above, so a reader wired to the wrong forest or to
  // the wrong member of the pair cannot pass.
  const double halfCauchyScale = 1.5;
  const double variances[3] = {0.0, 0.35, 0.125};
  const size_t trees[3] = {30, 15, 40};
  const double* bases[3] = {b0.data(), b1.data(), b2.data()};
  const size_t widths[3] = {1, 2, 2};

  const ResponseFamily families[2] = {ResponseFamily::probit,
                                      ResponseFamily::logistic};
  for (ResponseFamily family : families) {
    const double s = family == ResponseFamily::probit
                       ? 1.0
                       : std::numbers::pi / std::sqrt(3.0);
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 9310);
    SamplerOptions options;
    options.numChains = 1;

    AmplitudeSpec spec;
    spec.family = family;
    spec.forests.resize(3);
    for (size_t f = 0; f < 3; ++f) {
      spec.forests[f].forest.numTrees = trees[f];
      spec.forests[f].nodeScaleFactor = factors[f];
      spec.forests[f].nodeScaleDivisor = divisors[f];
      spec.forests[f].basis = bases[f];
      spec.forests[f].numBasisColumns = widths[f];
      if (f == 0)
        spec.forests[f].amplitudePriorScale = halfCauchyScale;
      else
        spec.forests[f].amplitudePriorVariance = variances[f];
    }
    std::unique_ptr<SamplerBase> sampler =
      createAmplitudeSampler(x.data(), y.data(), n, p, nullptr, nullptr, 1.0,
                             3.0, 0.37804942330213542, options, spec, &rng);
    check(sampler != nullptr, "calibration map: the K = 3 fixture builds");
    if (sampler != nullptr) {
      for (size_t f = 0; f < 3; ++f) {
        ForestCalibration calibration = sampler->forestCalibration(0, f);
        checkNear(calibration.priorScale,
                  factors[f] * s / (divisors[f] * norms[f]), 1.0e-12,
                  "calibration map: the reported scale is the map's product");
        // and what makes that read exact rather than up to a transform
        check(calibration.responseScale == 1.0 && calibration.k == 1.0,
              "calibration map: a latent family reports the identity transform "
              "and k pinned at one");
        // the map's own decomposition, reported per forest: each factor
        // against its literal, and the three together against the scale they
        // decompose, so a column filled from a neighbouring buffer is visible
        check(calibration.nodeScaleFactor == factors[f] &&
                calibration.nodeScaleDivisor == divisors[f],
              "calibration map: the reported factor and divisor are the "
              "forest's own");
        checkNear(calibration.basisRowNorm, norms[f], 1.0e-12,
                  "calibration map: the reported row norm is the forest's "
                  "own");
        checkNear(calibration.nodeScaleFactor * s /
                    (calibration.nodeScaleDivisor * calibration.basisRowNorm),
                  calibration.priorScale, 1.0e-12,
                  "calibration map: the reported decomposition reproduces the "
                  "reported scale");
        // and the amplitude prior, in whichever of its two exclusive spellings
        // the forest carries
        if (f == 0)
          check(calibration.amplitudePriorScale == halfCauchyScale &&
                  std::isnan(calibration.amplitudePriorVariance),
                "calibration map: a scale-mixture forest reports its "
                "half-Cauchy median and no variance");
        else
          check(calibration.amplitudePriorVariance == variances[f] &&
                  std::isnan(calibration.amplitudePriorScale),
                "calibration map: a fixed-variance forest reports its variance "
                "and no half-Cauchy median");
      }
    }

    ext_rng* conventionRng =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(conventionRng, 9311);
    AmplitudeSpec conventionSpec;
    conventionSpec.family = family;
    conventionSpec.forests.resize(2);
    conventionSpec.forests[0].forest.numTrees = 20;
    conventionSpec.forests[0].nodeScaleFactor = 1.25;
    conventionSpec.forests[0].nodeScaleDivisor = 0.5;
    conventionSpec.forests[0].basis = even.data();
    conventionSpec.forests[1].forest.numTrees = 10;
    conventionSpec.forests[1].nodeScaleFactor = 0.9;
    conventionSpec.forests[1].nodeScaleDivisor = 2.0;
    conventionSpec.forests[1].basis = zeros.data();
    std::unique_ptr<SamplerBase> convention = createAmplitudeSampler(
      x.data(), y.data(), n, p, nullptr, nullptr, 1.0, 3.0, 0.37804942330213542,
      options, conventionSpec, &conventionRng);
    check(convention != nullptr, "calibration map: the convention fixture "
                                 "builds");
    if (convention != nullptr) {
      checkNear(convention->forestCalibration(0, 0).priorScale,
                1.25 * s / (0.5 * 11.5), 1.0e-12,
                "calibration map: zero rows are excluded and an even nonzero "
                "count averages the two central order statistics");
      checkNear(convention->forestCalibration(0, 1).priorScale,
                0.9 * s / (2.0 * 1.0), 1.0e-12,
                "calibration map: an all-zero basis falls back to a unit row "
                "norm and still reports a finite scale");
      // and the convention is now READABLE rather than only inferable from the
      // product: the reported norm is the median the map used
      check(convention->forestCalibration(0, 0).basisRowNorm == 11.5 &&
              convention->forestCalibration(0, 1).basisRowNorm == 1.0,
            "calibration map: the reported row norm carries the median and "
            "the all-zero fallback");
    }

    ext_rng_destroy(conventionRng);
    ext_rng_destroy(rng);
  }
  printf("ok: BCF calibration map product\n");
}

// Chain::fitsWithoutOffset against a hand-built fitScale * combined +
// fitShift, and the recorded training write against that plus the offset.
// SINGLE-FOREST necessarily: combinedFits stays private and the public
// forestTotalFits equals it only off a combiner (the BCF leg of this gate is
// the R recombination cell). This is one of the two ARITHMETIC gates - the R
// identity cells cannot see inside the accessor, since the storeSample
// refactor puts both sides of their identity on it - so the fixture carries a
// NON-UNIT response range: at fitScale == 1 an accessor returning the internal
// scale would pass here.
static void testFitsWithoutOffset() {
  // a local stream, so adding this test shifts neither the shared runif01()
  // state nor the shared rng the hardcoded values downstream depend on
  std::uint64_t state = 20260815u;
  auto unif = [&]() {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return static_cast<double>(state >> 11) * 0x1.0p-53;
  };

  const size_t n = 200, p = 3, numChains = 2, numTrees = 25;
  std::vector<double> x(n * p), y(n), offset(n);
  for (double& v : x) v = unif();
  for (size_t i = 0; i < n; ++i) {
    y[i] = 100.0 + 40.0 * x[i] + 20.0 * x[i + n] + 2.0 * (unif() - 0.5);
    offset[i] = 3.0 * (unif() - 0.5);
  }

  ext_rng* rngs[numChains];
  for (size_t c = 0; c < numChains; ++c) {
    rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rngs[c], 4400 + static_cast<std::uint32_t>(c));
  }

  SamplerOptions options;
  options.numTrees = numTrees;
  options.numChains = numChains;
  ConstantLeafSampler sampler(x.data(), y.data(), n, p, nullptr, offset.data(),
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, rngs);

  std::vector<double> trainingFits(n * numChains);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(5, 1, results);

  check(std::fabs(sampler.fitScale() - 1.0) > 0.1,
        "fitsWithoutOffset: the fixture's response scale is not 1");

  std::vector<double> totals(n), actual(n);
  for (size_t c = 0; c < numChains; ++c) {
    ForestCalibration calibration = sampler.forestCalibration(c, 0);
    sampler.forestTotalFits(c, 0, totals.data());
    check(sampler.fitsWithoutOffset(c, actual.data()),
          "fitsWithoutOffset: a single-location sampler reports");
    bool matchesArithmetic = true, matchesRecorded = true;
    for (size_t i = 0; i < n; ++i) {
      double expected =
        calibration.responseScale * totals[i] + calibration.responseShift;
      if (actual[i] != expected) matchesArithmetic = false;
      if (trainingFits[i + c * n] != expected + offset[i])
        matchesRecorded = false;
    }
    check(matchesArithmetic,
          "fitsWithoutOffset: response scale times the combined fit, shifted");
    check(matchesRecorded,
          "fitsWithoutOffset: the recorded training write is it plus offset");
  }

  for (size_t c = 0; c < numChains; ++c) ext_rng_destroy(rngs[c]);
  printf("ok: fitsWithoutOffset (scale %.4f)\n", sampler.fitScale());
}

void runSamplerTests(ext_rng* rng) {
  testFitsWithoutOffset();
  testForestColumnRestriction(rng);
  testForestColumnRestrictionAllNeutral();
  testBCFTauModeratorRestriction(rng);
  testPrintTreesForest();
  testBCFTwoForest(rng);
  testBCFResponseSwap();
  testBCFInteractionLifetime();
  testBCFFixedGlue(rng);
  testBCFInterweave(rng);
  testBCFInterweaveKeepTrees(rng);
  testBCFGrowForestFromRoot();
  testAmplitudePerForestReplay();
  testBCFZeroMultiplierSnap();
  testBCFCombinerSeam();
  testCombinedFitsAssociation();
  testForestBasisOrdering();
  testGeneralAmplitudeConditional();
  testGeneralAmplitudeRidge();
  testAmplitudeOffsetIndexing();
  testForestWeights();
  testMultinomial(rng);
  testMultinomialCombinerSeam();
  testMultinomialGrowForestFromRoot();
  testMultinomialCountGrowForestFromRoot();
  testMultinomialSetCounts();
  testMultinomialCategoryOffset();
  testActiveRowsMultinomialKernel();
  testViewSamplerMatchesFull();
  testEndToEndGaussian(rng);
  testEndToEndGaussianFp32(rng);
  testGatherTailShapes(rng);
  testFusedSuffstatBankCombine();
  testFusedSuffstatMatchesStock(rng);
  testFusedSuffstatDeclines(rng);
  testRunCancellation(rng);
  testEndToEndProbit(rng);
  testMultiChain();
  testTestFitThreadInvariance();
  testPredictThreadPartition();
  testSetData(rng);
  testSetResponseScaleLock(rng);
  testSetDataResize(rng);
  testSetDataVarianceForestResize(rng);
  testSetDataQuantileShrink(rng);
  testSetDataProbit(rng);
  testMultiChainSetData();
  testEndToEndLogistic(rng);
  testWeightedLogistic(rng);
  testEndToEndCategorical(rng);
  testWideCategorical(rng);
  testPooledMaskSampler(rng);
  testActiveRows();
  testActiveRowsOnGrownForest();
  testSetWeightsAndTestOffset();
  testSetControlAndModel();
  testMissingEndToEnd();
  testLogLikelihood();
  testForestCalibration();
  testBCFCalibrationMap();
}
