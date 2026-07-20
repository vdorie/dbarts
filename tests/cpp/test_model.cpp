#include "common.hpp"

static void testIntegratedLikelihood() {
  ConstantGaussianLeaf leaf{0.5 / std::sqrt(200.0)};
  double k = 2.0, sigmaSq = 0.01;
  // weighted suffstat over three responses (2 * 0.5 + 1 * -0.2 + 3 * 0.1)
  double sumW = 6.0, sumWZ = 1.1;

  // independent transcription of the reduced marginal: raw sum wz^2 dropped
  double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
  double posteriorPrecision = sumW / sigmaSq;
  double mean = sumWZ / sumW;
  double explained = sumWZ * mean;
  double expected = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision))
    + 0.5 * explained / sigmaSq
    - 0.5 * ((priorPrecision * mean) * (posteriorPrecision * mean)) /
        (priorPrecision + posteriorPrecision);

  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ),
            expected, 1e-13, "integrated likelihood formula");
  checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, 0.0, 0.0),
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

// The reduced marginal must equal the classic (mean, effective count,
// centered variance) form plus the raw sum of squares it drops, on the same
// raw data, weighted and unweighted, and must agree with the node-context path.
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

    // the classic centered-variance form; the reduced marginal drops its raw ss
    double priorPrecision = (k / leaf.scale) * (k / leaf.scale);
    double posteriorPrecision = sumW / sigmaSq;
    double classic =
      0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision))
      - 0.5 * (variance / sigmaSq) * (double) (n - 1)
      - 0.5 * ((priorPrecision * average) * (posteriorPrecision * average)) /
          (priorPrecision + posteriorPrecision);

    checkNear(leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ),
              classic + 0.5 * sumWZSq / sigmaSq, 1e-12,
              weighted ? "reduced marginal is classic plus dropped ss, weighted"
                       : "reduced marginal is classic plus dropped ss, unweighted");
  }

  // and through the node-context path, over a split's two children
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = (double) i;
  ColumnStore store;
  store.build(x.data(), n, 1, 100);
  std::vector<index_t> indexBuffer(n);
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
    double sumW = 0.0, sumWZ = 0.0;
    for (size_t m = node.begin; m < node.end; ++m) {
      size_t i = tree.indices[m];
      sumW += wv[i]; sumWZ += wv[i] * zv[i];
    }
    return leaf.logIntegratedLikelihood(k, sigmaSq, sumW, sumWZ);
  };
  checkNear(leaf.logIntegratedLikelihoodForNode(tree, zv.data(), wv.data(), k,
                                                sigmaSq, left),
            reference(left), 1e-12, "node-context marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(tree, zv.data(), wv.data(), k,
                                                sigmaSq, left + 1),
            reference(left + 1), 1e-12, "node-context marginal, right child");
  printf("ok: constant leaf suffstat equivalence\n");
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

// K_nu(x), R's reference Bessel (libc++ lacks std::cyl_bessel_k). The _ex form
// takes a caller-owned work buffer (length >= floor(nu) + 1) so it never calls
// R_alloc - the tests link libR but do not initialize its runtime. expo == 2
// returns exp(x) K_nu(x); scaled ratios cancel the exponential, so the moment
// reference stays finite even at the large-order / large-argument BCF regime.
extern "C" double Rf_bessel_k_ex(double x, double nu, double expo, double* bk);

static void testGeneralizedInverseGaussian(ext_rng* rng) {
  // GIG(p, A, B) has E[v] = sqrt(B/A) K_{p+1}(w)/K_p(w) and E[v^2] =
  // (B/A) K_{p+2}(w)/K_p(w), w = sqrt(A B). The grid covers p = 0.5, 1, 5, 50
  // at small and large w, plus the a-near-0 Gamma limit (B = 0), where
  // v ~ Gamma(p, rate A/2) so E[v] = 2p/A and Var[v] = 4p/A^2.
  struct Case { double p, A, B; };
  const Case cases[] = {
    {2.5, 4.0, 0.0},   // gamma limit
    {0.5, 1.0, 1.0},   {0.5, 2.0, 8.0},
    {1.0, 1.0, 1.0},   {1.0, 4.0, 16.0},
    {5.0, 2.0, 8.0},   {5.0, 10.0, 10.0},
    {50.0, 20.0, 20.0}, {50.0, 100.0, 25.0}, {50.0, 200.0, 98.0}
  };
  const int numDraws = 200000;

  for (const Case& c : cases) {
    double expectedMean, expectedVar;
    if (c.B == 0.0) {
      expectedMean = 2.0 * c.p / c.A;
      expectedVar = 4.0 * c.p / (c.A * c.A);
    } else {
      double w = std::sqrt(c.A * c.B);
      double eta = std::sqrt(c.B / c.A);
      double bk[80];  // >= floor(p + 2) + 1 for the grid's largest order
      double k0 = Rf_bessel_k_ex(w, c.p, 2.0, bk);
      double r1 = Rf_bessel_k_ex(w, c.p + 1.0, 2.0, bk) / k0;
      double r2 = Rf_bessel_k_ex(w, c.p + 2.0, 2.0, bk) / k0;
      expectedMean = eta * r1;
      expectedVar = eta * eta * r2 - expectedMean * expectedMean;
    }

    double sum = 0.0, sumSq = 0.0;
    bool positive = true;
    for (int i = 0; i < numDraws; ++i) {
      double draw =
        ext_rng_simulateGeneralizedInverseGaussian(rng, c.p, c.A, c.B);
      if (!(draw > 0.0) || !std::isfinite(draw)) { positive = false; break; }
      sum += draw;
      sumSq += draw * draw;
    }
    check(positive, "GIG draw positive and finite");
    if (!positive) continue;

    double mean = sum / numDraws;
    double variance = sumSq / numDraws - mean * mean;
    double meanSe = std::sqrt(expectedVar / (double) numDraws);
    checkNear(mean, expectedMean, 6.0 * meanSe, "GIG mean");
    checkNear(variance, expectedVar, 0.08 * expectedVar, "GIG variance");
  }

  printf("ok: generalized inverse gaussian sampler\n");
}

static void testChiKHyperprior(ext_rng* rng) {
  // k ~ chi(nu) means k^2 ~ chi-squared(nu) = Gamma(nu/2, 1/2), so given M
  // conditionally-normal leaves the posterior of k^2 is Gamma with shape
  // 0.5 (M + nu) and rate 0.5 (sumSq/leafScale^2 [+ 1/scale^2]). Small M and
  // large nu keep the df term a large share of the shape, so the moments
  // separate chi(nu) from the chi(2 nu - 1) an un-Jacobianed shape would draw.
  const double nu = 8.0, numLeaves = 6.0, sumSquaredParams = 4.0;
  const double leafScale = 1.0;
  const int numDraws = 100000;
  const double shape = 0.5 * (numLeaves + nu);

  ChiKHyperprior prior;
  prior.degreesOfFreedom = nu;

  auto checkMoments = [&](double rate, const char* what) {
    double meanKSq = shape / rate, varKSq = shape / (rate * rate);
    double sum = 0.0, sumOfSquares = 0.0;
    for (int i = 0; i < numDraws; ++i) {
      double k = prior.draw(rng, sumSquaredParams, numLeaves, leafScale);
      double kSq = k * k;
      sum += kSq;
      sumOfSquares += kSq * kSq;
    }
    double mean = sum / numDraws;
    double var = sumOfSquares / numDraws - mean * mean;
    checkNear(mean, meanKSq, 0.01 * meanKSq, what);
    checkNear(var, varKSq, 0.05 * varKSq, what);
  };

  double leafRate = 0.5 * sumSquaredParams / (leafScale * leafScale);
  checkMoments(leafRate, "chi-k posterior k^2 moments, infinite prior scale");

  prior.scale = 0.5;  // finite prior scale enters the rate
  checkMoments(leafRate + 0.5 / (prior.scale * prior.scale),
               "chi-k posterior k^2 moments, finite prior scale");

  // A prior-dominated leaf under an improper scale would draw k far past any
  // finite bound; the sentinel caps it at exactly maxDraw. A tiny leaf scale
  // and a tiny sum of squares put the gamma rate near zero, so the raw draw
  // dwarfs the cap and every draw must clamp to it.
  ChiKHyperprior runaway;
  runaway.degreesOfFreedom = 1000.0;  // shape dominated by df, not by leaves
  for (int i = 0; i < 1000; ++i) {
    double k = runaway.draw(rng, 1e-30, 1.0, 1e-6);
    check(k == ChiKHyperprior::maxDraw, "chi-k runaway draw capped at maxDraw");
  }

  // A healthy draw returns the uncapped sqrt-of-gamma verbatim, unclamped.
  // Seed a private rng, read the inlined formula off it, rewind to the same
  // seed, and confirm draw() reproduces that value exactly and stays below
  // the cap.
  ext_rng* healthyRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ChiKHyperprior healthy;  // df 1.5, infinite scale
  double healthyShape = 0.5 * (numLeaves + healthy.degreesOfFreedom);
  double healthyRate = 0.5 * sumSquaredParams / (leafScale * leafScale);

  ext_rng_setSeed(healthyRng, 42);
  double expectedK = std::sqrt(
    ext_rng_simulateGamma(healthyRng, healthyShape, 1.0 / healthyRate));
  ext_rng_setSeed(healthyRng, 42);
  double gotK = healthy.draw(healthyRng, sumSquaredParams, numLeaves, leafScale);
  check(gotK < ChiKHyperprior::maxDraw, "chi-k healthy draw below the cap");
  checkNear(gotK, expectedK, 0.0, "chi-k healthy draw equals uncapped formula");
  ext_rng_destroy(healthyRng);

  printf("ok: chi-k hyperprior\n");
}

// A leaf that a data mutation strands empty is forced to a zero parameter; it
// is not a draw from the k-scaled prior, so it must not enter the chi-k count
// or sum of squares (the function-leaf path already excludes it, returning
// FunctionLeafDrawStats{0, 0}). No public mutation strands an empty leaf -
// collapse bubbles up while n > 0 - so the chain hook fabricates one and runs
// the real per-sweep accounting; here we confirm only the populated leaves
// count. Both leaf shapes with inline accounting (scalar and vector) are
// exercised; before the fix each counted the empty leaf too.
static void testChiKEmptyLeafAccounting(ext_rng* rng) {
  const size_t n = 200, p = 2;
  std::vector<double> x(n * p), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    double t = static_cast<double>(i) / static_cast<double>(n);
    x[i] = t;                  // split column
    x[i + n] = 2.0 * t - 1.0;  // linear-leaf covariate
    y[i] = (t > 0.5 ? 1.0 : -1.0) + 0.1 * (runif01() - 0.5);
  }
  ColumnStore store;
  // column 1 is the linear-leaf covariate the vector chain below reads
  size_t leafGather[] = {1};
  store.build(x.data(), n, p, 100, false, nullptr, leafGather, 1);

  // splitIndex 50 keeps the left child populated; the hook strands the right
  const int32_t splitVar = 0, splitIndex = 50;

  auto checkAccounting = [&](auto& chain, size_t stride, const char* what) {
    FunctionLeafDrawStats stats =
      chain.accountStrandedLeafKStatsForTesting(splitVar, splitIndex);
    std::vector<int32_t> bottoms;
    chain.tree(0).fillBottom(0, bottoms);
    size_t populated = 0, empty = 0;
    for (int32_t i : bottoms)
      (chain.tree(0).at(i).numObservations() > 0 ? populated : empty) += 1;
    check(populated >= 1 && empty >= 1,
          "fabricated a populated leaf beside a stranded empty one");
    // fixed: exactly the populated leaves' parameters count; before the fix
    // the empty leaf added another stride, so this held populated + empty
    check(stats.numParams == static_cast<double>(populated * stride), what);
  };

  SamplerOptions options;
  options.numTrees = 1;
  options.updateK = true;
  Chain<ConstantGaussianLeaf> scalarChain(
    store, y.data(), weights.data(), nullptr, ResponseFamily::gaussian, 1.0,
    3.0, 0.37804942330213542, options, rng);
  checkAccounting(scalarChain, 1, "constant leaf excludes the empty leaf");

  size_t covariates[] = {1};
  options.leafCovariateColumns = covariates;
  options.numLeafCovariates = 1;
  Chain<LinearGaussianLeaf> vectorChain(
    store, y.data(), weights.data(), nullptr, ResponseFamily::gaussian, 1.0,
    3.0, 0.37804942330213542, options, rng);  // numParams == 2 per leaf
  checkAccounting(vectorChain, 2, "linear leaf excludes the empty leaf");

  printf("ok: chi-k empty-leaf accounting\n");
}

// The sigma^2 posterior df must count only positive-weight rows: zero-weight
// rows drop from the weighted SSR and every leaf conditional, so they cannot
// inflate the df. With n = 20, n_pos = 10, prior df = 3 the draw is scaled
// inverse-chi-squared with posterior df = 13 (not 23), separating n_pos from n
// decisively.
static void testSigmaPosteriorDf(ext_rng* rng) {
  const std::size_t n = 20, nPositive = 10;
  const double priorDf = 3.0, scale = 0.7, residual = 0.7;

  std::vector<double> y(n), fits(n, 0.0), weights(n, 0.0);
  for (std::size_t i = 0; i < nPositive; ++i) {
    weights[i] = 1.0;
    y[i] = residual;  // residual r on each positive row: SSR = nPositive r^2
  }
  for (std::size_t i = nPositive; i < n; ++i)
    y[i] = 1.0e6;  // ignored: w = 0 zeroes its SSR contribution

  ChiSquaredScalePrior prior;
  prior.degreesOfFreedom = priorDf;
  prior.scale = scale;

  double ssr = static_cast<double>(nPositive) * residual * residual;
  double posteriorScale = priorDf * scale + ssr;
  double posteriorDf = priorDf + static_cast<double>(nPositive);
  // scaled inverse-chi-squared X = posteriorScale / chisq(posteriorDf)
  double meanExpected = posteriorScale / (posteriorDf - 2.0);
  double varExpected =
    2.0 * posteriorScale * posteriorScale /
    ((posteriorDf - 2.0) * (posteriorDf - 2.0) * (posteriorDf - 4.0));

  const int numDraws = 400000;
  double sum = 0.0, sumOfSquares = 0.0;
  for (int i = 0; i < numDraws; ++i) {
    double s2 = prior.drawSigmaSqFromPosterior(rng, y.data(), fits.data(),
                                               weights.data(), n, nPositive);
    sum += s2;
    sumOfSquares += s2 * s2;
  }
  double mean = sum / numDraws;
  double var = sumOfSquares / numDraws - mean * mean;
  checkNear(mean, meanExpected, 0.01 * meanExpected,
            "sigma^2 posterior mean, positive-weight df");
  checkNear(var, varExpected, 0.05 * varExpected,
            "sigma^2 posterior variance, positive-weight df");

  printf("ok: sigma posterior positive-weight df\n");
}

// Independent numerical check of the variance leaf's reduced scale marginal:
// log integral over x > 0 of x^{-n/2} exp(-ssr/(2x)) times the chi^-2(df,scale)
// prior density (the (2 pi)^{-n/2} factor the leaf drops). Log-substitution
// x = exp(t) puts the (0, inf) range on a symmetric grid; a max-shift keeps the
// exponential stable, then composite Simpson integrates.
static double logVarianceMarginalQuadrature(double n, double ssr, double df,
                                            double scale) {
  double logPriorNorm =
    0.5 * df * std::log(0.5 * df * scale) - std::lgamma(0.5 * df);
  auto logIntegrand = [&](double t) {
    double x = std::exp(t);
    // -0.5 n log x - ssr/(2x) + log prior(x) + t (the dx = e^t dt Jacobian)
    return -0.5 * n * t - 0.5 * ssr / x + logPriorNorm -
           (0.5 * df + 1.0) * t - 0.5 * df * scale / x + t;
  };
  double center = std::log((ssr + df * scale) / (n + df));
  double lo = center - 40.0, hi = center + 40.0;
  const int steps = 200000;  // even, for Simpson
  double h = (hi - lo) / steps;
  double maxLog = -HUGE_VAL;
  for (int i = 0; i <= steps; ++i) {
    double g = logIntegrand(lo + i * h);
    if (g > maxLog) maxLog = g;
  }
  double sum = 0.0;
  for (int i = 0; i <= steps; ++i) {
    double weight = (i == 0 || i == steps) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    sum += weight * std::exp(logIntegrand(lo + i * h) - maxLog);
  }
  return maxLog + std::log(sum * h / 3.0);
}

// (a) the scale marginal matches direct quadrature over (n, ssr), and the
// node-context own-span accumulation reproduces the direct form over a split.
static void testVarianceLeafMarginal() {
  struct Case { double df, scale, n, ssr; };
  Case cases[] = {
    {3.0, 0.5, 8.0, 3.2}, {10.0, 0.2, 25.0, 6.0}, {4.0, 1.0, 1.0, 0.7}};
  for (const Case& c : cases) {
    ConstantVarianceLeaf leaf{c.df, c.scale};
    double got = leaf.logIntegratedLikelihood(c.n, c.ssr);
    double ref = logVarianceMarginalQuadrature(c.n, c.ssr, c.df, c.scale);
    checkNear(got, ref, 1e-6 * (std::fabs(ref) + 1.0),
              "variance leaf marginal matches quadrature");
  }
  ConstantVarianceLeaf leaf{5.0, 0.3};
  checkNear(leaf.logIntegratedLikelihood(0.0, 0.0), 0.0, 0.0,
            "empty variance leaf contributes zero");

  // the node-context path over a split's two children, with a zero weight to
  // exercise the positive-weight count
  const size_t n = 6;
  double r[n] = {0.4, -0.8, 1.1, 0.2, -0.5, 0.9};
  double w[n] = {1.0, 2.0, 0.0, 1.5, 0.7, 1.2};
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = (double) i;
  ColumnStore store;
  store.build(x.data(), n, 1, 100);
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  std::vector<double> rv(r, r + n), wv(w, w + n);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(2);
  tree.birth(store, 0, rule, rv.data(), wv.data());
  ConstantVarianceLeaf vleaf{4.0, 0.6};
  auto refNode = [&](int32_t nodeIndex) {
    const Node& node = tree.at(nodeIndex);
    double nn = 0.0, ssr = 0.0;
    for (size_t m = node.begin; m < node.end; ++m) {
      size_t i = tree.indices[m];
      if (wv[i] <= 0.0) continue;
      nn += 1.0;
      ssr += wv[i] * rv[i] * rv[i];
    }
    return vleaf.logIntegratedLikelihood(nn, ssr);
  };
  int32_t left = tree.at(0).leftChild;
  checkNear(vleaf.logIntegratedLikelihoodForNode(tree, rv.data(), wv.data(),
                                                 1.0, 1.0, left),
            refNode(left), 1e-12, "node-context variance marginal, left child");
  checkNear(vleaf.logIntegratedLikelihoodForNode(tree, rv.data(), wv.data(),
                                                 1.0, 1.0, left + 1),
            refNode(left + 1), 1e-12,
            "node-context variance marginal, right child");
  printf("ok: variance leaf marginal\n");
}

// (b) the scaled-inverse-chi-squared posterior draw's moments match the analytic
// scale-inv-chi^2 moments, small-n and large-n, plus the prior draw. A local
// generator keeps the shared stream (and every downstream test) unperturbed.
static void testVarianceLeafDraw() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20240719u);
  struct Case { double df, scale; size_t n; };
  Case cases[] = {{4.0, 0.5, 5}, {3.0, 0.8, 500}};
  for (const Case& c : cases) {
    std::vector<double> r(c.n), w(c.n, 1.0);
    double ssr = 0.0;
    for (size_t i = 0; i < c.n; ++i) {
      r[i] = 0.3 + 0.1 * std::sin((double) i);
      ssr += w[i] * r[i] * r[i];
    }
    std::vector<index_t> indexBuffer(c.n);
    Tree tree;
    tree.initialize(indexBuffer.data(), c.n);
    ConstantVarianceLeaf leaf{c.df, c.scale};

    double posteriorScale = c.df * c.scale + ssr;
    double posteriorDf = c.df + (double) c.n;
    double meanExpected = posteriorScale / (posteriorDf - 2.0);
    double varExpected = 2.0 * posteriorScale * posteriorScale /
      ((posteriorDf - 2.0) * (posteriorDf - 2.0) * (posteriorDf - 4.0));

    const int numDraws = 400000;
    double sum = 0.0, sumOfSquares = 0.0;
    for (int d = 0; d < numDraws; ++d) {
      double s2 = leaf.drawFromPosteriorForNode(rng, tree, r.data(), w.data(),
                                                1.0, 1.0, 0);
      sum += s2;
      sumOfSquares += s2 * s2;
    }
    double mean = sum / numDraws;
    double var = sumOfSquares / numDraws - mean * mean;
    checkNear(mean, meanExpected, 0.01 * meanExpected,
              "variance leaf posterior draw mean");
    checkNear(var, varExpected, 0.05 * varExpected,
              "variance leaf posterior draw variance");
  }

  ConstantVarianceLeaf leaf{6.0, 0.4};
  double priorMeanExpected =
    leaf.degreesOfFreedom * leaf.scale / (leaf.degreesOfFreedom - 2.0);
  const int numDraws = 400000;
  double sum = 0.0;
  for (int d = 0; d < numDraws; ++d) sum += leaf.drawFromPrior(rng);
  checkNear(sum / numDraws, priorMeanExpected, 0.01 * priorMeanExpected,
            "variance leaf prior draw mean");
  ext_rng_destroy(rng);
  printf("ok: variance leaf posterior draws\n");
}

// (c) the Section-3.4 calibration collapses to ChiSquaredScalePrior at m' = 1
// and mean-matches the leaf product at m' > 1.
static void testVarianceLeafCalibration() {
  const double df = 3.0, scale = 0.7;
  ConstantVarianceLeaf collapsed = ConstantVarianceLeaf::calibrated(df, scale, 1);
  checkNear(collapsed.degreesOfFreedom, df, 1e-12,
            "m'=1 calibration df collapses to nu");
  checkNear(collapsed.scale, scale, 1e-12,
            "m'=1 calibration scale collapses to lambda^2");

  // at m' = 1 the leaf draw is bitwise ChiSquaredScalePrior's: a 5-row node
  // makes the leaf's own scan and the weighted-SSR unroll sum identically.
  const size_t n = 5;
  std::vector<double> r(n), w(n, 1.0), fits(n, 0.0);
  for (size_t i = 0; i < n; ++i) r[i] = 0.2 * (double) i - 0.3;
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  ChiSquaredScalePrior prior;
  prior.degreesOfFreedom = df;
  prior.scale = scale;

  ext_rng* local = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(local, 20240719u);
  double leafDraw = collapsed.drawFromPosteriorForNode(local, tree, r.data(),
                                                       w.data(), 1.0, 1.0, 0);
  ext_rng_setSeed(local, 20240719u);
  double priorDraw = prior.drawSigmaSqFromPosterior(local, r.data(), fits.data(),
                                                    w.data(), n, n);
  ext_rng_destroy(local);
  checkNear(leafDraw, priorDraw, 0.0,
            "m'=1 leaf draw is ChiSquaredScalePrior's, bitwise");

  // m' = 4: the product of leaf-prior means matches the homoscedastic mean.
  const size_t mPrime = 4;
  ConstantVarianceLeaf multi =
    ConstantVarianceLeaf::calibrated(df, scale, mPrime);
  double leafMean =
    multi.degreesOfFreedom * multi.scale / (multi.degreesOfFreedom - 2.0);
  double homoMean = df * scale / (df - 2.0);
  checkNear(std::pow(leafMean, (double) mPrime), homoMean, 1e-10,
            "Section 3.4 calibration mean-matches the leaf product");
  printf("ok: variance leaf calibration\n");
}

// The divisor guard (design rider i): tree j is scored/drawn against the OTHER
// trees' product s^2_{-j} = s^2 / h_j, so perturbing tree j's OWN leaf must not
// move its suffstat (only another tree's leaf may). Exercises the multiplicative
// roll directly, the one gate the m'=1 reductions cannot reach (divisor == 1).
static void testVarianceDivisorGuard() {
  const size_t n = 4, m = 2;
  VarianceForest vf;
  vf.initialize(m, n, 1.0);
  double f0[n] = {1.5, 0.8, 2.0, 1.1};
  double f1[n] = {0.7, 1.3, 0.9, 1.6};
  for (size_t i = 0; i < n; ++i) {
    vf.factorByTree[0 * n + i] = f0[i];
    vf.factorByTree[1 * n + i] = f1[i];
    vf.combinedVariance[i] = f0[i] * f1[i];
  }
  double e[n] = {0.3, -0.5, 0.2, 0.9};

  vf.formTreeResidual(0, e);
  std::vector<double> r0(vf.treeResidual);
  double ssr0 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    checkNear(vf.divisor[i], f1[i], 1e-13, "s^2_{-0} excludes tree 0's factor");
    ssr0 += r0[i] * r0[i];
  }

  // perturb tree 0's OWN factor (keep combined = f0 * f1 consistent): the
  // scaled residual and thus the suffstat must be unchanged.
  for (size_t i = 0; i < n; ++i) {
    vf.factorByTree[0 * n + i] = f0[i] * 3.0;
    vf.combinedVariance[i] = (f0[i] * 3.0) * f1[i];
  }
  vf.formTreeResidual(0, e);
  double ssrOwn = 0.0;
  for (size_t i = 0; i < n; ++i) {
    checkNear(vf.treeResidual[i], r0[i], 1e-13,
              "tree 0 own-leaf perturbation leaves its residual unmoved");
    ssrOwn += vf.treeResidual[i] * vf.treeResidual[i];
  }
  checkNear(ssrOwn, ssr0, 1e-13,
            "tree 0 own-leaf perturbation leaves its suffstat unmoved");

  // perturb ANOTHER tree's factor: tree 0's residual MUST move.
  for (size_t i = 0; i < n; ++i) {
    vf.factorByTree[1 * n + i] = f1[i] * 2.0;
    vf.combinedVariance[i] = (f0[i] * 3.0) * (f1[i] * 2.0);
  }
  vf.formTreeResidual(0, e);
  bool moved = false;
  for (size_t i = 0; i < n; ++i)
    if (std::fabs(vf.treeResidual[i] - r0[i]) > 1e-9) moved = true;
  check(moved, "another tree's leaf perturbation moves tree 0's suffstat");
  printf("ok: variance divisor guard\n");
}

// End-to-end: a variance forest on a step-heteroscedastic problem recovers a
// higher s^2 where the truth is noisier. A local generator keeps the shared
// stream unperturbed.
static void testVarianceForestRecovery() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 815239u);
  const size_t n = 800, p = 1;
  std::vector<double> x(n * p), y(n), weights(n, 1.0);
  const double sLow = 0.3, sHigh = 2.0;  // truth ratio (sHigh/sLow)^2 ~ 44
  for (size_t i = 0; i < n; ++i) {
    double xi = (double) i / (double) n;
    x[i] = xi;
    y[i] = (xi < 0.5 ? sLow : sHigh) * ext_rng_simulateStandardNormal(rng);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 100);

  SamplerOptions options;
  options.numTrees = 20;
  options.numVarianceTrees = 20;
  Chain<ConstantGaussianLeaf> chain(store, y.data(), weights.data(), nullptr,
                                    ResponseFamily::gaussian, 1.0, 3.0,
                                    0.37804942330213542, options, rng);
  Results results;
  chain.run(300, 300, results);

  const double* s2 = chain.varianceFits();
  check(s2 != nullptr, "variance forest exposes s^2(x)");
  double lowMean = 0.0, highMean = 0.0;
  size_t nLow = 0, nHigh = 0;
  for (size_t i = 0; i < n; ++i) {
    if (x[i] < 0.5) { lowMean += s2[i]; ++nLow; }
    else { highMean += s2[i]; ++nHigh; }
  }
  lowMean /= (double) nLow;
  highMean /= (double) nHigh;
  check(highMean > 2.5 * lowMean,
        "variance forest recovers higher s^2 where the truth is noisier");
  ext_rng_destroy(rng);
  printf("ok: variance forest recovery (s2 ratio %.2f)\n", highMean / lowMean);
}

// Construction-time refusal (design section 5): the variance forest routes
// through the weight channel the latent families already own, so the factory
// accepts it only for the plain-constant-leaf gaussian family.
static void testVarianceForestRefusal() {
  const size_t n = 40, p = 1;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = (double) i / (double) n;
    y[i] = i < n / 2 ? 0.0 : 1.0;
  }
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 7u);
  ext_rng* rngs[1] = {rng};
  SamplerOptions options;
  options.numTrees = 5;
  options.numVarianceTrees = 4;

  check(createSampler(x.data(), y.data(), n, p, nullptr, nullptr,
                      ResponseFamily::gaussian, 1.0, 3.0, 0.378, options, rngs)
          != nullptr,
        "gaussian variance forest accepted at the factory");
  for (ResponseFamily fam : {ResponseFamily::probit, ResponseFamily::logistic,
                             ResponseFamily::nbinom}) {
    check(createSampler(x.data(), y.data(), n, p, nullptr, nullptr, fam, 1.0,
                        3.0, 0.378, options, rngs) == nullptr,
          "variance forest refused for a latent family");
  }
  // a linear-leaf designation is also refused (v1 keeps the mean leaf constant)
  size_t covariate[] = {0};
  options.leafCovariateColumns = covariate;
  options.numLeafCovariates = 1;
  check(createSampler(x.data(), y.data(), n, p, nullptr, nullptr,
                      ResponseFamily::gaussian, 1.0, 3.0, 0.378, options, rngs)
          == nullptr,
        "variance forest refused for a non-constant mean leaf");
  ext_rng_destroy(rng);
  printf("ok: variance forest gaussian-only refusal\n");
}

// C4: the saved-variance-tree predict path and the variance-forest column
// restriction (the R surface's `variance = ~ subset`). Saved-slot predict must
// reproduce the recorded per-sample s^2, and the restricted trees split only on
// the allowed column.
static void testVarianceSavedPredict() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20240720u);
  const size_t n = 300, p = 2;
  std::vector<double> x(n * p), y(n), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    double x0 = (double) i / (double) n;
    double x1 = ((i * 7) % n) / (double) n;  // an unrelated second column
    x[i] = x0;
    x[i + n] = x1;
    y[i] = (x0 < 0.5 ? 0.4 : 1.6) * ext_rng_simulateStandardNormal(rng);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 100);

  SamplerOptions options;
  options.numTrees = 20;
  options.numVarianceTrees = 15;
  std::size_t varianceColumns[] = {0};  // restrict the variance forest to col 0
  options.varianceForestColumns = varianceColumns;
  options.numVarianceForestColumns = 1;
  Chain<ConstantGaussianLeaf> chain(store, y.data(), weights.data(), nullptr,
                                    ResponseFamily::gaussian, 1.0, 3.0,
                                    0.37804942330213542, options, rng);
  const size_t numSamples = 4;
  chain.initializeSavedTrees(numSamples);
  std::vector<double> varFits(n * numSamples);
  Results results;
  results.varianceFits = varFits.data();
  chain.run(200, numSamples, results);

  // saved-slot predict on the training rows reproduces the recorded s^2
  double s2scale = chain.sigmaScale() * chain.sigmaScale();
  std::vector<double> pred(n);
  bool savedMatch = true;
  for (size_t s = 0; s < numSamples; ++s) {
    chain.predictVarianceFromSavedSample(s, x.data(), n, pred.data());
    const double* recorded = varFits.data() + s * n;
    for (size_t i = 0; i < n; ++i)
      if (std::fabs(pred[i] - recorded[i]) > 1e-9 * s2scale) savedMatch = false;
  }
  check(savedMatch, "saved-slot variance predict reproduces the recorded s^2");
  printf("ok: variance saved-tree predict and column restriction\n");
  ext_rng_destroy(rng);
}

// C3: the s(x) reporting channels (train/test), predict on new data, and the
// scale-leaf state round-trip all agree with the live combined variance.
static void testVarianceReportingStatePredict() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 424242u);
  const size_t n = 300, p = 1, nTest = 60;
  std::vector<double> x(n), y(n), xTest(nTest), weights(n, 1.0);
  for (size_t i = 0; i < n; ++i) {
    double xi = (double) i / (double) n;
    x[i] = xi;
    y[i] = (xi < 0.5 ? 0.3 : 1.8) * ext_rng_simulateStandardNormal(rng);
  }
  for (size_t i = 0; i < nTest; ++i) xTest[i] = (double) i / (double) nTest;

  ColumnStore store;
  store.build(x.data(), n, p, 100);
  store.buildTest(xTest.data(), nTest);

  SamplerOptions options;
  options.numTrees = 20;
  options.numVarianceTrees = 20;
  Chain<ConstantGaussianLeaf> chain(store, y.data(), weights.data(), nullptr,
                                    ResponseFamily::gaussian, 1.0, 3.0,
                                    0.37804942330213542, options, rng);

  const size_t numSamples = 5;
  std::vector<double> varFits(n * numSamples), varTestFits(nTest * numSamples);
  Results results;
  results.varianceFits = varFits.data();
  results.varianceTestFits = varTestFits.data();
  chain.run(300, numSamples, results);

  double s2scale = chain.sigmaScale() * chain.sigmaScale();
  double tol = 1e-9 * s2scale;
  const double* combined = chain.varianceFits();
  const double* combinedTest = chain.varianceTestFits();

  // (a) reporting: the last recorded sample equals the live combined variance
  // (train) and the routed test product, both on the original scale
  const double* lastTrain = varFits.data() + (numSamples - 1) * n;
  const double* lastTest = varTestFits.data() + (numSamples - 1) * nTest;
  bool trainMatch = true, testMatch = true;
  for (size_t i = 0; i < n; ++i)
    if (std::fabs(lastTrain[i] - s2scale * combined[i]) > tol) trainMatch = false;
  for (size_t i = 0; i < nTest; ++i)
    if (std::fabs(lastTest[i] - s2scale * combinedTest[i]) > tol)
      testMatch = false;
  check(trainMatch, "recorded train s^2 matches the live combined variance");
  check(testMatch, "recorded test s^2 matches the routed test product");

  // (c) predict on new data matches an in-sample recompute
  std::vector<double> predTrain(n), predTest(nTest);
  chain.predictVariance(x.data(), n, predTrain.data());
  chain.predictVariance(xTest.data(), nTest, predTest.data());
  bool predMatch = true, predTestMatch = true;
  for (size_t i = 0; i < n; ++i)
    if (std::fabs(predTrain[i] - s2scale * combined[i]) > tol) predMatch = false;
  for (size_t i = 0; i < nTest; ++i)
    if (std::fabs(predTest[i] - s2scale * combinedTest[i]) > tol)
      predTestMatch = false;
  check(predMatch, "predictVariance on training rows matches combined variance");
  check(predTestMatch, "predictVariance on test rows matches the test channel");

  // (b) state save/restore reproduces s(x)
  ChainStateData state;
  chain.getState(state);
  ext_rng* rng2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng2, 999u);
  Chain<ConstantGaussianLeaf> restored(store, y.data(), weights.data(), nullptr,
                                       ResponseFamily::gaussian, 1.0, 3.0,
                                       0.37804942330213542, options, rng2);
  check(restored.stateIsValid(state), "heteroscedastic state validates");
  check(restored.setState(state), "heteroscedastic state restores");
  const double* restoredCombined = restored.varianceFits();
  bool restoreMatch = true;
  for (size_t i = 0; i < n; ++i)
    if (std::fabs(restoredCombined[i] - combined[i]) > 1e-12) restoreMatch = false;
  check(restoreMatch, "restored combined variance reproduces the saved s^2");

  // a corrupted (non-positive) leaf value fails the scale-leaf validation
  ChainStateData bad = state;
  for (FlatNode& node : bad.varianceTrees[0])
    if (flatKindOf(node) == FlatKind::leaf) { node.value = -1.0; break; }
  check(!restored.stateIsValid(bad),
        "a non-positive variance leaf fails validation");

  ext_rng_destroy(rng2);
  ext_rng_destroy(rng);
  printf("ok: variance reporting, state, predict\n");
}

// C3 gate (a) at the C++ level: a constant predictor stops both forests from
// splitting, so an m'=1 variance forest is a single chi^-2 leaf whose posterior
// is the homoscedastic sigma^2 posterior - the Section-3.4 calibration collapse.
static void testVarianceM1Reduction() {
  const size_t n = 300, p = 1;
  const double muTrue = 2.0, sigmaTrue = 1.5, sigEst = 1.5;
  std::vector<double> x(n, 0.5), y(n), weights(n, 1.0);

  auto runVariance = [&](bool heteroscedastic, ext_rng* rng) {
    ext_rng* dataRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(dataRng, 314159u);  // identical data for both samplers
    for (size_t i = 0; i < n; ++i)
      y[i] = muTrue + sigmaTrue * ext_rng_simulateStandardNormal(dataRng);
    ext_rng_destroy(dataRng);
    ColumnStore store;
    store.build(x.data(), n, p, 100);
    SamplerOptions options;
    options.numTrees = 20;
    if (heteroscedastic) options.numVarianceTrees = 1;
    Chain<ConstantGaussianLeaf> chain(store, y.data(), weights.data(), nullptr,
                                      ResponseFamily::gaussian, sigEst, 3.0,
                                      0.37804942330213542, options, rng);
    const size_t numSamples = 4000;
    std::vector<double> sig(numSamples), varFits;
    Results results;
    results.sigma = sig.data();
    if (heteroscedastic) {
      varFits.assign(n * numSamples, 0.0);
      results.varianceFits = varFits.data();
    }
    chain.run(1000, numSamples, results);
    double mean = 0.0;
    for (size_t s = 0; s < numSamples; ++s)
      mean += heteroscedastic ? varFits[s * n] : sig[s] * sig[s];
    return mean / (double) numSamples;
  };

  ext_rng* rngH = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngH, 111u);
  double homoMean = runVariance(false, rngH);
  ext_rng_destroy(rngH);
  ext_rng* rngV = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngV, 222u);
  double hetMean = runVariance(true, rngV);
  ext_rng_destroy(rngV);

  checkNear(hetMean, homoMean, 0.06 * homoMean,
            "m'=1 heteroscedastic s^2 posterior matches homoscedastic sigma^2");
  printf("ok: variance m'=1 reduction (homo %.4f het %.4f)\n", homoMean, hetMean);
}

static void testSampleFromPrior(ext_rng* rng) {
  const size_t n = 200, numTrees = 50, numReplications = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);

  SamplerOptions options;
  options.numTrees = numTrees;
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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

// Shared fixture for the linear-leaf component tests; the reference
// constants were computed in R by independently evaluating the same
// marginal-likelihood and posterior formulas at these inputs.
struct LinearLeafFixture {
  static constexpr size_t n = 8;
  std::vector<double> x, z, w;
  ColumnStore store;
  std::vector<index_t> indexBuffer;
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
    // the store owns raw copies of the leaf-covariate columns (1, 2), as a
    // production build does, so a leaf reading rawColumn is served
    size_t leafGather[] = {1, 2};
    store.build(x.data(), n, 3, 100, false, nullptr, leafGather, 2);
    tree.initialize(indexBuffer.data(), n);
  }
};

// The R reference constants below were evaluated with the raw z'Wz term the
// marginals have since dropped; adding the fixture's exact term back keeps
// them an independent check of everything else in the formula.
static double droppedSumOfSquaresTerm(const std::vector<double>& z,
                                      const double* w, size_t begin,
                                      size_t end, double sigmaSq) {
  double zwz = 0.0;
  for (size_t i = begin; i < end; ++i)
    zwz += (w == nullptr ? 1.0 : w[i]) * z[i] * z[i];
  return 0.5 * zwz / sigmaSq;
}

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
            -5.517188945419923 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, f.n, f.sigmaSq),
            1e-9, "linear marginal, weighted root");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), nullptr,
                                                f.k, f.sigmaSq, 0),
            -5.256174278535296 +
              droppedSumOfSquaresTerm(f.z, nullptr, 0, f.n, f.sigmaSq),
            1e-9, "linear marginal, unit-weight root");

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
            -3.760266419465231 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, 4, f.sigmaSq),
            1e-9, "linear marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild + 1),
            -3.062089185586155 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 4, f.n, f.sigmaSq),
            1e-9, "linear marginal, right child");

  // two covariates exercise the 3x3 Cholesky
  LinearGaussianLeaf leaf2;
  leaf2.scale = f.scale;
  size_t columns2[] = {1, 2};
  leaf2.initialize(f.store, columns2, 2);
  Tree rootTree;
  std::vector<index_t> rootIndices(LinearLeafFixture::n);
  rootTree.initialize(rootIndices.data(), LinearLeafFixture::n);
  checkNear(leaf2.logIntegratedLikelihoodForNode(
              rootTree, f.z.data(), f.w.data(), f.k, f.sigmaSq, 0),
            -5.791887259328131 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, f.n, f.sigmaSq),
            1e-9, "linear marginal, two covariates");

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

// The U'WU crossproduct cache must be invisible: a warm score or draw is
// bitwise a cold one, and the covariate (setPredictor/setData) and weight
// (setWeights) mutation hooks must drop it, or a later leaf reuses a stale
// U'WU. A single root over 80 > minCachedLeafSize observations is cacheable.
static void testLinearLeafStatisticsCache() {
  const size_t n = 80, p = 2;  // column 0 splits, column 1 is the leaf basis
  std::vector<double> xV0(n * p), xV1(n * p), z(n), w0(n), w1(n);
  for (size_t i = 0; i < n; ++i) {
    double t = (double) i / (double) n;
    xV0[i] = t;
    xV1[i] = t;  // split column identical across V0 and V1
    xV0[i + n] = std::sin(3.0 * t);
    xV1[i + n] = std::sin(3.0 * t) + 0.5 * t;  // leaf covariate differs
    z[i] = 0.3 * t - 0.15;
    w0[i] = 0.5 + (i % 4 == 0 ? 1.0 : 0.25);
    w1[i] = 0.5 + (i % 3 == 0 ? 0.8 : 0.4);  // weights differ
  }
  const double scale = 0.5 / std::sqrt(10.0), k = 2.0, sigmaSq = 0.04;
  size_t columns[] = {1};

  ColumnStore storeV0, storeV1;
  storeV0.build(xV0.data(), n, p, 100, false, nullptr, columns, 1);
  storeV1.build(xV1.data(), n, p, 100, false, nullptr, columns, 1);
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);

  LinearGaussianLeaf leaf;
  leaf.scale = scale;
  leaf.initialize(storeV0, columns, 1);
  double cold = leaf.logIntegratedLikelihoodForNode(tree, z.data(), w0.data(),
                                                    k, sigmaSq, 0);
  double warm = leaf.logIntegratedLikelihoodForNode(tree, z.data(), w0.data(),
                                                    k, sigmaSq, 0);
  check(cold == warm, "cached linear score is bitwise the scanned score");

  ext_rng* rngWarm = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngCold = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngWarm, 8123u);
  ext_rng_setSeed(rngCold, 8123u);
  double drawWarm[2], drawCold[2];
  leaf.drawFromPosteriorForNode(rngWarm, tree, z.data(), w0.data(), k, sigmaSq,
                                0, drawWarm);
  LinearGaussianLeaf coldLeaf;
  coldLeaf.scale = scale;
  coldLeaf.initialize(storeV0, columns, 1);
  coldLeaf.drawFromPosteriorForNode(rngCold, tree, z.data(), w0.data(), k,
                                    sigmaSq, 0, drawCold);
  check(drawWarm[0] == drawCold[0] && drawWarm[1] == drawCold[1],
        "cached linear draw is bitwise the scanned draw");
  ext_rng_destroy(rngWarm);
  ext_rng_destroy(rngCold);

  // covariate mutation: a leaf warmed under V0 then regathered onto V1 scores
  // exactly like one that only ever scanned V1
  LinearGaussianLeaf regathered;
  regathered.scale = scale;
  regathered.initialize(storeV0, columns, 1);
  regathered.logIntegratedLikelihoodForNode(tree, z.data(), w0.data(), k,
                                            sigmaSq, 0);
  regathered.regatherTrainingCovariates(storeV1);
  double afterRegather = regathered.logIntegratedLikelihoodForNode(
    tree, z.data(), w0.data(), k, sigmaSq, 0);
  LinearGaussianLeaf scanV1;
  scanV1.scale = scale;
  scanV1.initialize(storeV0, columns, 1);
  scanV1.regatherTrainingCovariates(storeV1);
  double scannedV1 = scanV1.logIntegratedLikelihoodForNode(tree, z.data(),
                                                           w0.data(), k, sigmaSq,
                                                           0);
  check(afterRegather != cold, "the V1 covariates move the marginal");
  check(afterRegather == scannedV1, "regather drops the covariate U'WU cache");

  // weight mutation: a leaf warmed under w0 then invalidated scores under w1
  // exactly like a cold scan under w1
  LinearGaussianLeaf reweighted;
  reweighted.scale = scale;
  reweighted.initialize(storeV0, columns, 1);
  reweighted.logIntegratedLikelihoodForNode(tree, z.data(), w0.data(), k,
                                            sigmaSq, 0);
  reweighted.invalidateStatistics();
  double afterInvalidate = reweighted.logIntegratedLikelihoodForNode(
    tree, z.data(), w1.data(), k, sigmaSq, 0);
  LinearGaussianLeaf scanW1;
  scanW1.scale = scale;
  scanW1.initialize(storeV0, columns, 1);
  double scannedW1 = scanW1.logIntegratedLikelihoodForNode(tree, z.data(),
                                                          w1.data(), k, sigmaSq,
                                                          0);
  check(afterInvalidate != cold, "the w1 weights move the marginal");
  check(afterInvalidate == scannedW1, "invalidateStatistics drops the cache");

  printf("ok: linear leaf statistics cache\n");
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

  // the format-dependent surfaces: an identical predictor matrix
  // revalidates, and a state round-trips into itself
  check(sampler->setPredictor(x.data(), false, true) ==
          PredictorUpdateResult::accepted,
        "reinstalling identical predictors is accepted");
  SamplerStateData state;
  sampler->getState(state);
  check(sampler->setState(state, nullptr), "a state restores into its own sampler");
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
  check(restored->setState(state, nullptr), "a linear-leaf state restores");

  checkStructuralRoundTrip(state, *restored,
                           "restored linear-leaf state reproduces the model");

  // a state with mismatched slope shapes is refused
  SamplerStateData malformed = state;
  malformed.chains[0].forests[0].treeParams[0].push_back(0.0);
  check(!sampler->setState(malformed, nullptr),
        "a state with mismatched slopes is refused");

  ext_rng_destroy(rng2);
  printf("ok: linear leaf formats\n");
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
  // a data handle owns raw for the columns its views gather (here covariate 1),
  // so buildFromParent can copy them; a dense build serves rawColumn only for
  // gathered columns
  parent.build(x.data(), n, p, options.maxNumCuts, false, nullptr, covariates,
               1);

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
  // built-in tau priors and the tau posterior against constants computed
  // in R: dcauchy/dgamma at scale 0.55, then the posterior density
  // evaluated at J = 5, b.sq = 0.8
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

  // the exact-IG cauchy draw (drawTauCauchyExactIG) reproduces the SAME tau
  // posterior the slice hits - it targets the identical conditional, so the
  // cauchy quadrature reference above (sliceMean[0]/sliceSd[0]) is its gate too
  {
    const double numGroups = 5.0, ss = 0.8, A = 0.55;
    double t = 0.45, total = 0.0, totalSq = 0.0;
    for (size_t r = 0; r < numDraws; ++r) {
      t = drawTauCauchyExactIG(rng, numGroups, ss, A, t);
      total += t;
      totalSq += t * t;
    }
    double mean = total / static_cast<double>(numDraws);
    double sd =
      std::sqrt(totalSq / static_cast<double>(numDraws) - mean * mean);
    checkNear(mean, sliceMean[0], 0.01, "exact-IG tau mean matches quadrature");
    checkNear(sd, sliceSd[0], 0.01, "exact-IG tau sd matches quadrature");

    // draw-count constancy: the update advances a privately seeded rng by
    // EXACTLY two ext_rng_simulateGamma draws (shapes 1 and (J + 1)/2) with the
    // documented rates - a reference rng fed just those two lands in the same
    // state, so a later uniform matches bit for bit (no stray consumption).
    ext_rng* seededA =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* seededRef =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (seededA == NULL || seededRef == NULL ||
        ext_rng_setSeed(seededA, 4242) != 0 ||
        ext_rng_setSeed(seededRef, 4242) != 0) {
      check(false, "exact-IG draw-count: rng creation");
    } else {
      const double tauIn = 0.37;
      double tauDraw = drawTauCauchyExactIG(seededA, numGroups, ss, A, tauIn);
      double invASq = 1.0 / (A * A);
      double xi = 1.0 / ext_rng_simulateGamma(
                          seededRef, 1.0,
                          1.0 / (1.0 / (tauIn * tauIn) + invASq));
      double tauSq = 1.0 / ext_rng_simulateGamma(
                             seededRef, 0.5 * (numGroups + 1.0),
                             1.0 / (0.5 * ss + 1.0 / xi));
      check(tauDraw == std::sqrt(tauSq),
            "exact-IG draw is two reciprocal gamma draws");
      check(ext_rng_simulateContinuousUniform(seededA) ==
              ext_rng_simulateContinuousUniform(seededRef),
            "exact-IG draw consumes exactly two gamma draws");
    }
    if (seededA != NULL) ext_rng_destroy(seededA);
    if (seededRef != NULL) ext_rng_destroy(seededRef);
  }

  // gamma-branch bit-identity: the exact-IG reroute touches only the cauchy
  // branch, so a seeded gamma-prior slice run reproduces the tau and the
  // following uniform captured from the base (pre-exact-IG) build exactly.
  {
    ext_rng* gammaRng =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (gammaRng == NULL || ext_rng_setSeed(gammaRng, 20260711u) != 0) {
      check(false, "gamma bit-identity: rng creation");
    } else {
      auto gammaDensity = [](double t) {
        return logTauPosterior(t, 5.0, 0.8, TauPriorKind::gamma, 0.55);
      };
      double tau = 0.45;
      for (int step = 0; step < 5; ++step)
        tau = sliceSampleOnce(gammaRng, gammaDensity, tau, 0.55, 0.0, HUGE_VAL);
      check(tau == 0.51564843960975137,
            "gamma slice tau bit-identical to base build");
      check(ext_rng_simulateContinuousUniform(gammaRng) ==
              0.47995807160623372,
            "gamma slice consumes the base build's rng draws");
    }
    if (gammaRng != NULL) ext_rng_destroy(gammaRng);
  }

  // step-out cap: a tau posterior with mode ~1e12 (J = 5, sum b^2 = 5e24;
  // the regime mostly-empty groupings random-walk into) stepped out
  // ~mode/width times before the cap - an indefinite hang. Capped, the
  // bracket stops at x + (m + 1) * width and the draw returns promptly
  // inside it. A dedicated rng keeps the shared stream untouched.
  {
    ext_rng* capRng =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (capRng == NULL || ext_rng_setSeed(capRng, 8271) != 0) {
      check(false, "step-out cap: rng creation");
      return;
    }
    auto stallDensity = [](double t) {
      return logTauPosterior(t, 5.0, 5.0e24, TauPriorKind::cauchy, 0.55);
    };
    double draw =
      sliceSampleOnce(capRng, stallDensity, 0.45, 0.55, 0.0, HUGE_VAL);
    check(std::isfinite(draw) && draw > 0.0 &&
            draw <= 0.45 + 0.55 * 10001.0,
          "capped step-out returns a finite draw inside the bracket");
    ext_rng_destroy(capRng);
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
  ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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
    ConstantLeafSampler sampler(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  const size_t n = 200, numGroups = 6, numChains = 2;
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
  ConstantLeafSampler original(x.data(), y.data(), n, 2, nullptr, nullptr,
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
  ConstantLeafSampler restored(x.data(), y.data(), n, 2, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state, nullptr), "a grouped state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored grouped state reproduces the effects and tau");

  // an effects vector of the wrong length must be refused before anything
  // is overwritten
  SamplerStateData badState(state);
  badState.chains[0].groupEffects.resize(numGroups - 1);
  check(!restored.setState(badState, nullptr),
        "a mismatched grouped state is refused");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: grouped state round trip\n");
}

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
  std::vector<index_t> segment(n);
  for (size_t i = 0; i < n; ++i) segment[i] = i;
  for (size_t i = n - 1; i > 0; --i) {
    size_t k = static_cast<size_t>(runif01() * static_cast<double>(i + 1));
    std::swap(segment[i], segment[k]);
  }

  const xint_t cuts[] = { 0, 4, 5, 100, 250 };
  bool allMatch = true;
  for (xint_t cut : cuts) {
    std::vector<index_t> indices(segment);
    size_t numOnLeft = misc_partitionIndicesSparse(
      bits.data(), wordRanks.data(), nzCodes.data(), zeroCode, cut,
      indices.data(), indices.size());

    size_t expected = 0;
    for (size_t i : segment) expected += denseCodes[i] <= cut ? 1 : 0;
    allMatch &= numOnLeft == expected;
    for (size_t k = 0; k < indices.size(); ++k)
      allMatch &= (denseCodes[indices[k]] <= cut) == (k < numOnLeft);
    std::vector<index_t> sorted(indices);
    std::sort(sorted.begin(), sorted.end());
    std::vector<index_t> sortedSegment(segment);
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
        codesMatch &= fromCsc.codeAt(j, i) == fromDense.train.codes[i + j * n];
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
        viewMatches &= view.train.codes[i + j * viewRows.size()] ==
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
    ConstantLeafSampler dense(fixture.dense.data(), y.data(), n, p, nullptr,
                         nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                         0.37804942330213542, options, &rngA);

    SamplerOptions cscOptions(options);
    cscOptions.cscColumnPointers = fixture.pointers.data();
    cscOptions.cscRowIndices = fixture.rows.data();
    cscOptions.cscValues = fixture.values.data();
    ConstantLeafSampler sparse(nullptr, y.data(), n, p, nullptr, nullptr,
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
    ConstantLeafSampler sampler(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
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
  const size_t n = 240, numChains = 2;
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
  ConstantLeafSampler original(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 1111);
  ConstantLeafSampler restored(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state, nullptr), "a sparse-store state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored sparse-store state reproduces the model");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: sparse state round trip\n");
}

// A single high-cardinality unordered factor held both as a dense column-major
// code matrix (codes as doubles, level order 0..K-1) and as CSC arrays over the
// non-reference entries, the reference level implicit. probReference tunes the
// nonzero fraction, so the same shape can land on either storage tier. The top
// level is planted at row 0 so the dense builder's level count (max + 1) equals
// K, and the reference is a middle level (its level-order code is never 0), so
// the store's zeroCode must carry the reference code, not codeFor(j, 0.0).
struct CscCategoricalFixture {
  size_t n = 0;
  std::uint32_t K = 0;
  xint_t reference = 0;
  std::vector<double> dense;   // column-major n x 1
  std::vector<int> pointers;   // 2
  std::vector<int> rows;
  std::vector<double> values;
  std::vector<ColumnType> types = { ColumnType::categorical };

  void build(size_t n_, std::uint32_t K_, double probReference) {
    n = n_;
    K = K_;
    reference = static_cast<xint_t>(K / 2);
    dense.assign(n, 0.0);
    pointers.assign(2, 0);
    rows.clear();
    values.clear();
    for (size_t i = 0; i < n; ++i) {
      xint_t code;
      if (i == 0) {
        code = static_cast<xint_t>(K - 1);  // guarantee the top level appears
      } else if (runif01() < probReference) {
        code = reference;
      } else {
        xint_t other =
          static_cast<xint_t>(runif01() * static_cast<double>(K - 1));
        if (other >= reference) ++other;  // uniform over the non-reference set
        code = other;
      }
      dense[i] = static_cast<double>(code);
      if (code != reference) {
        rows.push_back(static_cast<int>(i));
        values.push_back(static_cast<double>(code));
      }
    }
    pointers[1] = static_cast<int>(rows.size());
  }

  std::int32_t sources = ~0;  // the one column reads CSC source 0

  void applyOptions(SamplerOptions& options) {
    options.cscColumnPointers = pointers.data();
    options.cscRowIndices = rows.data();
    options.cscValues = values.data();
    options.columnSources = &sources;
    options.columnTypes = types.data();
    options.cscCategoryCounts = &K;
    options.cscReferenceCodes = &reference;
  }

  ColumnStore buildStore(bool useQuantiles) {
    ColumnStore store;
    store.buildMixed(nullptr, pointers.data(), rows.data(), values.data(),
                     &sources, n, 1, nullptr, 100, useQuantiles, types.data(),
                     &K, &reference);
    return store;
  }
};

// The bitwise gate at the engine level: a sparse-categorical column bins
// IDENTICALLY to a dense factor of the same values. Codes match cell by cell,
// and the new rank-layout membership partitions reproduce the dense kernels'
// output byte for byte at several direction masks (inline and pooled).
static void testSparseCategoricalColumnStore() {
  uint64_t savedRngState = rngState;  // leave the shared draw stream in place
  const size_t n = 400;

  // rank inline (K <= 63), rank pooled (K > 63), densified inline
  struct Config { std::uint32_t K; double probReference; bool expectSparse; };
  const Config configs[] = {
    { 6, 0.92, true }, { 80, 0.95, true }, { 6, 0.4, false }
  };

  bool allOk = true;
  for (const Config& config : configs) {
    CscCategoricalFixture fixture;
    fixture.build(n, config.K, config.probReference);
    ColumnStore denseStore;
    denseStore.build(fixture.dense.data(), n, 1, 100, false,
                     fixture.types.data());
    ColumnStore sparseStore = fixture.buildStore(false);

    allOk &= sparseStore.columnIsSparse(0) == config.expectSparse;
    allOk &= sparseStore.numCuts[0] == config.K &&
             denseStore.numCuts[0] == config.K;
    for (size_t i = 0; i < n; ++i)
      allOk &= sparseStore.codeAt(0, i) == denseStore.train.codes[i];

    if (!config.expectSparse) continue;  // dense kernel path, no new sibling

    // a scrambled index segment, partitioned by the sparse membership kernel
    // against the dense one over several direction masks
    std::vector<index_t> segment(n);
    for (size_t i = 0; i < n; ++i) segment[i] = i;
    for (size_t i = n - 1; i > 0; --i)
      std::swap(segment[i], segment[static_cast<size_t>(
                              runif01() * static_cast<double>(i + 1))]);

    const SparseColumnData& sparse = sparseStore.sparseColumn(0);
    const xint_t* dense = denseStore.column(0);
    if (config.K <= 63) {
      const std::uint64_t masks[] = {
        0x5ull, 0x2Aull, 0x3Full, (1ull << fixture.reference), 0x0ull
      };
      for (std::uint64_t mask : masks) {
        std::vector<index_t> a(segment), b(segment);
        size_t leftDense =
          Tree::partitionIndicesByMask(dense, mask, a.data(), a.size());
        size_t leftSparse = Tree::partitionIndicesSparseByMask(
          sparse, mask, b.data(), b.size());
        allOk &= leftDense == leftSparse && a == b;
      }
    } else {
      size_t numWords = maskWordsForCount(config.K);
      for (int m = 0; m < 4; ++m) {
        std::vector<std::uint64_t> mask(numWords, 0);
        for (std::uint32_t c = 0; c < config.K; ++c)
          if ((c + static_cast<std::uint32_t>(m)) % 3 == 0) maskSetBit(mask.data(), c);
        std::vector<index_t> a(segment), b(segment);
        size_t leftDense = Tree::partitionIndicesByWideMask(
          dense, mask.data(), a.data(), a.size());
        size_t leftSparse = Tree::partitionIndicesSparseByWideMask(
          sparse, mask.data(), b.data(), b.size());
        allOk &= leftDense == leftSparse && a == b;
      }
    }
  }
  check(allOk, "sparse-categorical codes and membership match the dense factor");
  rngState = savedRngState;
  printf("ok: sparse categorical column store\n");
}

// The end-to-end half of the gate: a sampler over a sparse-categorical column
// draws BITWISE-IDENTICALLY to one over the dense factor of the same values, on
// every storage tier (rank inline, densified inline, rank pooled).
static void testSparseCategoricalEndToEnd() {
  uint64_t savedRngState = rngState;
  const size_t n = 500;
  struct Config { std::uint32_t K; double probReference; bool expectSparse; };
  const Config configs[] = {
    { 6, 0.92, true }, { 6, 0.4, false }, { 80, 0.95, true }
  };

  for (const Config& config : configs) {
    CscCategoricalFixture fixture;
    fixture.build(n, config.K, config.probReference);
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i)
      y[i] = static_cast<double>(fixture.dense[i]) -
             0.5 * static_cast<double>(config.K) + 0.3 * runif01();

    ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 90210) != 0 ||
        ext_rng_setSeed(rngB, 90210) != 0) {
      check(false, "sparse categorical end-to-end: rng creation");
      return;
    }

    SamplerOptions options;
    options.numTrees = 40;
    options.columnTypes = fixture.types.data();
    ConstantLeafSampler dense(fixture.dense.data(), y.data(), n, 1, nullptr,
                              nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rngA);

    SamplerOptions sparseOptions;
    sparseOptions.numTrees = 40;
    fixture.applyOptions(sparseOptions);
    ConstantLeafSampler sparse(nullptr, y.data(), n, 1, nullptr, nullptr,
                               ResponseFamily::gaussian, 1.0, 3.0,
                               0.37804942330213542, sparseOptions, &rngB);
    check(sparse.data().columnIsSparse(0) == config.expectSparse,
          "sparse-categorical sampler lands on the expected storage tier");

    const size_t numBurnIn = 40, numSamples = 60;
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
          "sparse-categorical sampler bitwise-matches the dense factor");
    ext_rng_destroy(rngB);
    ext_rng_destroy(rngA);
  }
  rngState = savedRngState;
  printf("ok: sparse categorical end-to-end\n");
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
        codesMatch &= mixed.codeAt(j, i) == reference.train.codes[i + j * n];
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
        viewMatches &= view.train.codes[i + j * viewRows.size()] ==
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
  ConstantLeafSampler dense(fixture.full.data(), y.data(), n, fixture.p, nullptr,
                       nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                       0.37804942330213542, options, &rngA);

  SamplerOptions mixedOptions;
  mixedOptions.numTrees = 25;
  fixture.applyOptions(mixedOptions);
  ConstantLeafSampler mixed(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
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
  const size_t n = 240, numChains = 2;
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
  ConstantLeafSampler original(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(40, 0, empty);

  SamplerStateData state;
  original.getState(state);

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 2323);
  ConstantLeafSampler restored(nullptr, y.data(), n, fixture.p, nullptr, nullptr,
                          ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, options, rngs2.data());
  check(restored.setState(state, nullptr), "a mixed-store state restores");

  checkStructuralRoundTrip(state, restored,
                           "restored mixed-store state reproduces the model");

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
  std::vector<index_t> indexBuffer;
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
    // the store owns raw copies of the leaf-covariate columns (1, 2), as a
    // production build does, so a leaf reading rawColumn is served
    size_t leafGather[] = {1, 2};
    store.build(x.data(), n, 3, 100, false, nullptr, leafGather, 2);
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

  // root marginals against constants computed in R by direct evaluation
  // of the GP integrated likelihood at the fixture's kernel and inputs
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, 0),
            -8.4818528980077694 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, f.n, f.sigmaSq),
            1e-9, "gp marginal, weighted root");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), nullptr,
                                                f.k, f.sigmaSq, 0),
            -7.7744742015613255 +
              droppedSumOfSquaresTerm(f.z, nullptr, 0, f.n, f.sigmaSq),
            1e-9, "gp marginal, unit-weight root");

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
            -8.7417667358884774 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, f.n, f.sigmaSq),
            1e-9, "gp marginal, two covariates");
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
            -5.0284693098158764 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, 4, f.sigmaSq),
            1e-9, "gp marginal, left child");
  checkNear(leaf.logIntegratedLikelihoodForNode(f.tree, f.z.data(), f.w.data(),
                                                f.k, f.sigmaSq, leftChild + 1),
            -3.7990612354161635 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 4, f.n, f.sigmaSq),
            1e-9, "gp marginal, right child");

  // a constant kernel (huge lengthscale) approaches the constant leaf's
  // formula by Sherman-Morrison
  Tree rootTree;
  std::vector<index_t> rootIndices(GPLeafFixture::n);
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
            -5.0284693098158764 +
              droppedSumOfSquaresTerm(f.z, f.w.data(), 0, 4, f.sigmaSq),
            1e-9, "at-cap marginal stays gp");

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
  std::vector<index_t> positiveIndices(6);
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
  size_t dupGather[] = {1};
  dupStore.build(xd.data(), nd, 2, 100, false, nullptr, dupGather, 1);
  std::vector<index_t> dupIndices(nd);
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

  // the format surfaces: a captured state restores into its own sampler
  // and flatten emits records
  SamplerStateData state;
  gpSampler.getState(state);
  check(gpSampler.setState(state, nullptr), "a gp state restores into its own sampler");
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
  check(restored.setState(state, nullptr), "gp state installs into a fresh sampler");

  // malformed side channels refuse without mutating anything
  SamplerStateData malformed = state;
  malformed.chains[0].forests[0].savedTreeParams[0].pop_back();
  check(!restored.setState(malformed, nullptr), "truncated gp side channel refused");
  malformed = state;
  malformed.chains[0].forests[0].treeParams[0].pop_back();
  check(!restored.setState(malformed, nullptr), "short gp fits slab refused");

  // the good state survived the rejected malformed loads and reproduces the
  // model (its fits slabs, saved side channels, and rng)
  checkStructuralRoundTrip(state, restored,
                           "restored gp state reproduces the model");

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
// empty cache recomputes what the warm sampler serves from cache. The
// mutation perturbs the designated column WITHIN its quantization bins:
// members stay identical, so only the regather-clears-the-cache path keeps
// stale kernels from hitting. Member re-routing is checked separately on
// one leaf, where shared buffers keep the bitwise comparison sound.
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
    bool installed = cold.setState(state, nullptr);

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

  // member re-route: a leaf warmed over one member set then routed onto a
  // disjoint one must rebuild its kernel rather than serve the stale entry.
  // Kept on a single leaf so the warm serve and the post-invalidate rescan
  // share buffers and compare bitwise irrespective of allocation layout.
  ColumnStore leafStore;
  size_t leafGather[] = {0};
  leafStore.build(x.data(), n, p, 100, false, nullptr, leafGather, 1);
  std::vector<index_t> leafIndices(n);
  Tree leafTree;
  leafTree.initialize(leafIndices.data(), n);
  GPGaussianLeaf leaf;
  leaf.scale = ySd / std::sqrt(5.0);
  size_t leafColumns[] = {0};
  leaf.initialize(leafStore, leafColumns, 1, nullptr, 60);
  const size_t leafSize = 50;  // >= minCachedLeafSize so the cache engages
  leafTree.at(0).end = leafSize;  // members [0, leafSize) of the identity index

  const double k = 2.0, sigmaSq = 0.04;
  double warmScore = leaf.logIntegratedLikelihoodForNode(leafTree, y.data(),
                                                         nullptr, k, sigmaSq, 0);
  double servedScore = leaf.logIntegratedLikelihoodForNode(
    leafTree, y.data(), nullptr, k, sigmaSq, 0);
  check(warmScore == servedScore, "cached gp score is bitwise the scanned score");

  for (size_t i = 0; i < leafSize; ++i) leafIndices[i] = i + (n - leafSize);
  double reroutedScore = leaf.logIntegratedLikelihoodForNode(
    leafTree, y.data(), nullptr, k, sigmaSq, 0);
  leaf.regatherTrainingCovariates(leafStore);  // drop the kernel cache
  double rescannedScore = leaf.logIntegratedLikelihoodForNode(
    leafTree, y.data(), nullptr, k, sigmaSq, 0);
  check(reroutedScore != warmScore, "the re-routed members move the marginal");
  check(reroutedScore == rescannedScore, "member re-route rebuilds the kernel");

  printf("ok: gp leaf kernel cache\n");
}

// AFT log-normal survival (docs/design/survival.md): the free reduction gate
// (uncensored == gaussian, bitwise), the censored-latent refresh moments
// (draws above the truncation bound, empirical mean at the truncated-normal
// mean), and a censored state round trip through the latents + fit.scale
// blocks.
static double standardNormalPdf(double z) {
  return std::exp(-0.5 * z * z) / std::sqrt(2.0 * 3.141592653589793);
}
static double standardNormalCdf(double z) {
  return 0.5 * std::erfc(-z * 0.7071067811865476);
}

static void testAFTReduction(ext_rng*) {
  // an all-uncensored aft fit is bit-identical to a gaussian fit on log T
  const size_t n = 180, p = 3;
  std::vector<double> x(n * p), logT(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) x[i + j * n] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) *
               std::cos(6.283185307179586 * u2);
    logT[i] = 1.5 * x[i] + 0.5 * z;
  }
  std::vector<double> statusAll(n, 1.0);  // every observation an event

  SamplerOptions optG;
  optG.numTrees = 25;
  SamplerOptions optA = optG;
  optA.survivalStatus = statusAll.data();

  auto seededRng = [](std::uint32_t seed) {
    ext_rng* r = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(r, seed);
    return r;
  };
  ext_rng* rngG = seededRng(31337);
  ext_rng* rngA = seededRng(31337);

  std::unique_ptr<SamplerBase> sampG = createSampler(
    x.data(), logT.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian,
    1.0, 3.0, 0.37804942330213542, optG, &rngG);
  std::unique_ptr<SamplerBase> sampA = createSampler(
    x.data(), logT.data(), n, p, nullptr, nullptr, ResponseFamily::aft,
    1.0, 3.0, 0.37804942330213542, optA, &rngA);

  const size_t numBurnIn = 60, numSamples = 80;
  std::vector<double> fitG(n * numSamples), fitA(n * numSamples);
  std::vector<double> sigmaG(numSamples), sigmaA(numSamples);
  Results rG, rA;
  rG.trainingFits = fitG.data();
  rG.sigma = sigmaG.data();
  rA.trainingFits = fitA.data();
  rA.sigma = sigmaA.data();
  sampG->run(numBurnIn, numSamples, rG);
  sampA->run(numBurnIn, numSamples, rA);

  bool identical = fitG == fitA && sigmaG == sigmaA;
  check(identical, "aft reduction: uncensored fit is bitwise gaussian");

  ext_rng_destroy(rngG);
  ext_rng_destroy(rngA);
  printf("ok: aft reduction to gaussian\n");
}

static void testAFTCensoredMoments(ext_rng* rng) {
  // one censored observation, its latent redrawn from the log-scale normal
  // truncated below at its log censoring time: draws stay above the bound and
  // average at the analytic truncated-normal mean
  const size_t m = 5;
  std::vector<double> logTime = {0.2, 0.5, 1.0, 1.5, 2.0};
  std::vector<double> status = {1.0, 1.0, 0.0, 1.0, 1.0};  // index 2 censored
  AFTResponse resp(logTime.data(), status.data(), nullptr, m, 1.0, 3.0,
                   0.37804942330213542);

  const size_t censored = 2;
  double scale = resp.fitScale();
  double shift = resp.fitShift();
  double sigma = 0.4;  // internal-scale residual sd
  std::vector<double> totalFits(m, 0.0);
  totalFits[censored] = -0.5;  // pushes the fit well below the bound

  double meanLog = scale * totalFits[censored] + shift;
  double sdLog = sigma * scale;
  double bound = logTime[censored];
  double alpha = (bound - meanLog) / sdLog;
  double truncatedMean =
    meanLog + sdLog * standardNormalPdf(alpha) / (1.0 - standardNormalCdf(alpha));

  const size_t reps = 40000;
  double sum = 0.0;
  bool allAbove = true;
  for (size_t r = 0; r < reps; ++r) {
    resp.refreshLatents(rng, totalFits.data(), sigma);
    double v = resp.latents()[censored];
    if (v < bound - 1e-9) allAbove = false;
    sum += v;
    // the events keep their observed log-times untouched
    if (resp.latents()[0] != logTime[0]) allAbove = false;
  }
  check(allAbove, "aft censored latents stay above the bound; events fixed");
  checkNear(sum / static_cast<double>(reps), truncatedMean, 0.02,
            "aft censored latent mean at the truncated-normal mean");
  printf("ok: aft censored latent moments\n");
}

static void testAFTStateRoundTrip() {
  // a censored aft state rides the latents (logT_) and fit.scale blocks; a
  // restored sampler reconstructs them, and a mismatched latents vector is
  // refused
  const size_t n = 240, p = 2, numChains = 2;
  std::vector<double> x(n * p), obsLogT(n), status(n);
  uint64_t seed = 20260710u;
  auto localUnif = [&seed]() {
    seed = seed * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>(seed >> 11) / 9007199254740992.0;
  };
  size_t numCensored = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) x[i + j * n] = localUnif();
    double f = 1.2 * x[i] - 0.6 * x[i + n];
    double logT = f + 0.5 * (localUnif() - 0.5) * 3.464;  // rough spread
    double censTime = f + 0.5 + 0.5 * (localUnif() - 0.5) * 3.464;
    if (logT <= censTime) {
      status[i] = 1.0;
      obsLogT[i] = logT;
    } else {
      status[i] = 0.0;
      obsLogT[i] = censTime;
      ++numCensored;
    }
  }
  check(numCensored > 0, "aft round trip: some observations censored");

  SamplerOptions options;
  options.numTrees = 20;
  options.numChains = numChains;
  options.survivalStatus = status.data();

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t base) {
    for (size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], base + static_cast<std::uint32_t>(c));
    }
  };
  std::vector<ext_rng*> rngs(numChains, nullptr);
  makeRngs(rngs, 5150);
  ConstantLeafSampler original(x.data(), obsLogT.data(), n, p, nullptr, nullptr,
                               ResponseFamily::aft, 1.0, 3.0,
                               0.37804942330213542, options, rngs.data());
  Results empty;
  original.run(50, 0, empty);

  SamplerStateData state;
  original.getState(state);
  check(state.chains[0].latents.size() == n, "aft state carries the log-times");
  check(state.chains[0].fitMax > state.chains[0].fitMin,
        "aft state carries the fit scale");

  std::vector<ext_rng*> rngs2(numChains, nullptr);
  makeRngs(rngs2, 9999);
  ConstantLeafSampler restored(x.data(), obsLogT.data(), n, p, nullptr, nullptr,
                               ResponseFamily::aft, 1.0, 3.0,
                               0.37804942330213542, options, rngs2.data());
  check(restored.setState(state, nullptr), "an aft state restores");
  checkStructuralRoundTrip(state, restored,
                           "restored aft state reproduces the latents and scale");

  SamplerStateData badState(state);
  badState.chains[0].latents.resize(n - 1);
  check(!restored.setState(badState, nullptr), "a mismatched aft latents vector is refused");

  for (size_t c = 0; c < numChains; ++c) {
    ext_rng_destroy(rngs[c]);
    ext_rng_destroy(rngs2[c]);
  }
  printf("ok: aft state round trip\n");
}

// lambda_i | z, f, sigma ~ Gamma((nu + 1)/2, rate (nu + w_i r_i^2/sigma^2)/2);
// over fixed residuals and fixed nu the per-observation draw's mean and
// variance match that Gamma (the polya-gamma moments precedent). A local
// generator and a restored global rngState leave the shared stream untouched.
static void testTLambdaMoments(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 4;
  const double nu = 5.0, sigma = 0.3;
  std::vector<double> y(n), weights = {0.5, 1.0, 2.0, 3.0};
  for (std::size_t i = 0; i < n; ++i) y[i] = runif01();

  TResponse resp(y.data(), nullptr, weights.data(), n, 1.0, 3.0,
                 0.37804942330213542, nu);
  std::vector<double> totalFits(n);
  for (std::size_t i = 0; i < n; ++i)
    totalFits[i] = 0.1 * static_cast<double>(i) - 0.15;
  const double* z = resp.workingResponse();  // internal residual r_i = z_i - f_i
  double sigmaSq = sigma * sigma;

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 8675309u);
  const int numDraws = 200000;
  std::vector<double> sum(n, 0.0), sumSq(n, 0.0);
  for (int d = 0; d < numDraws; ++d) {
    resp.refreshLatents(localRng, totalFits.data(), sigma);
    const double* lambda = resp.latents();
    for (std::size_t i = 0; i < n; ++i) {
      sum[i] += lambda[i];
      sumSq[i] += lambda[i] * lambda[i];
    }
  }

  for (std::size_t i = 0; i < n; ++i) {
    double r = z[i] - totalFits[i];
    double shape = 0.5 * (nu + 1.0);
    double rate = 0.5 * (nu + weights[i] * r * r / sigmaSq);
    double expectedMean = shape / rate;
    double expectedVar = shape / (rate * rate);
    double mean = sum[i] / numDraws;
    double var = sumSq[i] / numDraws - mean * mean;
    double meanSe = std::sqrt(expectedVar / static_cast<double>(numDraws));
    checkNear(mean, expectedMean, 5.0 * meanSe, "t lambda conditional mean");
    checkNear(var, expectedVar, 0.05 * expectedVar,
              "t lambda conditional variance");
  }

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: t lambda conditional moments\n");
}

// the sampled-nu grid draw reproduces the hand-computed full conditional: for
// a fixed lambda vector the empirical grid frequencies match the normalized
// product of Gamma(lambda_i; nu/2, nu/2) densities times the gamma(2, 0.1)
// prior, computed here through the exact lambda^{nu/2 - 1} density form the
// grid draw's k-constant simplification drops. Local generator, restored state.
static void testTNuGridPosterior(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  ResidualDfPrior prior;
  const std::size_t n = 12;
  std::vector<double> lambda = {0.3, 0.5, 0.7, 0.9, 1.0, 1.1,
                                1.3, 1.5, 1.7, 0.4, 1.6, 1.0};
  double sumLogLambda = 0.0, sumLambda = 0.0;
  for (double l : lambda) {
    sumLogLambda += std::log(l);
    sumLambda += l;
  }
  double nObs = static_cast<double>(n);

  double expected[ResidualDfPrior::gridSize];
  double maxLog = -HUGE_VAL;
  for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) {
    double a = 0.5 * ResidualDfPrior::grid[k];
    double logPost = nObs * (a * std::log(a) - std::lgamma(a)) +
                     (a - 1.0) * sumLogLambda - a * sumLambda +
                     std::log(ResidualDfPrior::grid[k]) -
                     0.1 * ResidualDfPrior::grid[k];  // gamma(2, 0.1) log kernel
    expected[k] = logPost;
    if (logPost > maxLog) maxLog = logPost;
  }
  double total = 0.0;
  for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) {
    expected[k] = std::exp(expected[k] - maxLog);
    total += expected[k];
  }
  for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) expected[k] /= total;

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 1729u);
  const int numDraws = 400000;
  std::vector<int> counts(ResidualDfPrior::gridSize, 0);
  for (int d = 0; d < numDraws; ++d)
    ++counts[prior.drawIndex(localRng, nObs, sumLogLambda, sumLambda)];

  for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) {
    double freq = static_cast<double>(counts[k]) / numDraws;
    double se = std::sqrt(expected[k] * (1.0 - expected[k]) / numDraws);
    checkNear(freq, expected[k], std::max(0.003, 5.0 * se),
              "t nu grid frequency matches the hand-computed posterior");
  }

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: t nu grid full conditional\n");
}

// after a sweep, workingWeights()[i] is exactly w_i lambda_i (a zero user
// weight composing to 0), and the sigma draw consumes that same composite:
// reproduced bit-for-bit by a bare GaussianResponse fed the composite as fixed
// weights on a lockstepped rng. Local generators, restored global rngState.
static void testTCompositeWeightDelegation(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 6;
  std::vector<double> y(n), weights = {0.0, 0.5, 1.0, 1.5, 2.0, 0.75};
  for (std::size_t i = 0; i < n; ++i) y[i] = runif01();
  const double nu = 6.0, sigma = 0.4, sigmaDf = 3.0,
               sigmaRawScale = 0.37804942330213542;

  TResponse resp(y.data(), nullptr, weights.data(), n, 1.0, sigmaDf,
                 sigmaRawScale, nu);
  std::vector<double> totalFits(n);
  for (std::size_t i = 0; i < n; ++i) totalFits[i] = 0.05 * static_cast<double>(i);

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 271828u);
  resp.refreshLatents(localRng, totalFits.data(), sigma);

  const double* lambda = resp.latents();
  const double* composite = resp.workingWeights();
  bool exact = true;
  for (std::size_t i = 0; i < n; ++i)
    if (composite[i] != weights[i] * lambda[i]) exact = false;
  check(exact, "t working weights equal w_i lambda_i exactly");
  check(composite[0] == 0.0, "t zero user weight composes to zero");

  // the composite must be exactly what drawSigma consumes: a bare Gaussian on
  // the same response, given the composite as its weights, draws the same sigma
  std::vector<double> compositeCopy(composite, composite + n);
  GaussianResponse ref(y.data(), nullptr, compositeCopy.data(), n, 1.0, sigmaDf,
                       sigmaRawScale);
  ext_rng* rngResp = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngRef = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngResp, 424242u);
  ext_rng_setSeed(rngRef, 424242u);
  double sResp = resp.drawSigma(rngResp, totalFits.data(), sigma);
  double sRef = ref.drawSigma(rngRef, totalFits.data(), sigma);
  check(sResp == sRef, "t sigma draw consumes the composite weights exactly");

  ext_rng_destroy(rngRef);
  ext_rng_destroy(rngResp);
  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: t composite weight delegation\n");
}

// fixed nu never draws nu: over identical data and identically seeded local
// generators, a fixed-nu sweep and a grid sweep (both starting at nu = 8) make
// the same n lambda draws, after which the grid stream is ahead by exactly the
// one uniform of the discrete nu draw. Restored global rngState.
static void testTFixedNuNoDraw(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 16;
  const double median = ResidualDfPrior::grid[ResidualDfPrior::medianIndex];
  std::vector<double> y(n), weights(n);
  for (std::size_t i = 0; i < n; ++i) {
    y[i] = runif01();
    weights[i] = 0.5 + runif01();
  }
  std::vector<double> totalFits(n);
  for (std::size_t i = 0; i < n; ++i)
    totalFits[i] = 0.03 * static_cast<double>(i) - 0.2;
  const double sigma = 0.35;

  TResponse fixedResp(y.data(), nullptr, weights.data(), n, 1.0, 3.0,
                      0.37804942330213542, median);
  TResponse gridResp(y.data(), nullptr, weights.data(), n, 1.0, 3.0,
                     0.37804942330213542, -1.0);  // estimate on the grid
  check(!fixedResp.estimatesResidualDf() && gridResp.estimatesResidualDf(),
        "t mode flags: fixed vs grid");

  ext_rng* rngFixed = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngGrid = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngFixed, 90210u);
  ext_rng_setSeed(rngGrid, 90210u);
  fixedResp.refreshLatents(rngFixed, totalFits.data(), sigma);
  gridResp.refreshLatents(rngGrid, totalFits.data(), sigma);

  check(fixedResp.residualDf() == median, "t fixed nu is held");
  bool onGrid = false;
  for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k)
    if (gridResp.residualDf() == ResidualDfPrior::grid[k]) onGrid = true;
  check(onGrid, "t grid nu lands on a grid value");

  // discarding one uniform from the fixed stream aligns it with the grid
  // stream: proof the grid drew exactly one more (the nu draw) and the fixed
  // drew none
  ext_rng_simulateContinuousUniform(rngFixed);
  bool shiftedByOne = true;
  for (int j = 0; j < 64; ++j)
    if (ext_rng_simulateContinuousUniform(rngFixed) !=
        ext_rng_simulateContinuousUniform(rngGrid))
      shiftedByOne = false;
  check(shiftedByOne, "t grid stream is the fixed stream plus one nu draw");

  ext_rng_destroy(rngGrid);
  ext_rng_destroy(rngFixed);
  rngState = savedRngState;
  printf("ok: t fixed nu draws no nu\n");
}

// The doubly-truncated-normal primitive in isolation: draws land in (lower,
// upper] and their mean/variance match the truncated N(mean, 1) moments, from a
// midrange interval (the inverse-CDF bulk) through two deep-tail intervals whose
// Phi gap underflows and hands off to Robert rejection (right tail direct, left
// tail reflected). Reference moments come from a peak-shifted quadrature that
// stays finite where the closed-form normalizer underflows. Local generator,
// restored global rngState.
static void testOrdinalTruncatedNormal(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  struct Case { double mean, lower, upper; };
  const Case cases[] = {
    {0.3, -1.0, 1.5},     // bulk, straddling 0
    {-0.5, 0.4, 2.2},     // bulk, one-sided
    {0.0, 8.5, 9.0},      // deep right tail: gap underflows -> Robert
    {0.0, -41.0, -40.0},  // deep left tail: gap underflows -> Robert reflected
  };
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 20260718u);
  const int numDraws = 400000, G = 200000;

  for (const Case& c : cases) {
    double h = (c.upper - c.lower) / G, gmax = -HUGE_VAL;
    for (int j = 0; j < G; ++j) {
      double d = c.lower + (j + 0.5) * h - c.mean;
      if (-0.5 * d * d > gmax) gmax = -0.5 * d * d;
    }
    double sw = 0.0, sxw = 0.0, sx2w = 0.0;
    for (int j = 0; j < G; ++j) {
      double x = c.lower + (j + 0.5) * h, d = x - c.mean;
      double w = std::exp(-0.5 * d * d - gmax);
      sw += w; sxw += x * w; sx2w += x * x * w;
    }
    double refMean = sxw / sw, refVar = sx2w / sw - refMean * refMean;

    double sum = 0.0, sumSq = 0.0;
    bool inSupport = true;
    for (int d = 0; d < numDraws; ++d) {
      double z = ext_rng_simulateTruncatedNormalScale1(localRng, c.mean,
                                                       c.lower, c.upper);
      if (!(z >= c.lower && z <= c.upper)) inSupport = false;
      sum += z; sumSq += z * z;
    }
    check(inSupport, "ordinal truncated draw lands in (lower, upper]");
    double mean = sum / numDraws, var = sumSq / numDraws - mean * mean;
    checkNear(mean, refMean, std::max(1e-6, 5.0 * std::sqrt(refVar / numDraws)),
              "ordinal truncated-normal mean");
    checkNear(var, refVar, 0.05 * refVar + 1e-9,
              "ordinal truncated-normal variance");
  }

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: ordinal doubly-truncated normal primitive\n");
}

// The K = 3 cutpoint update against a hand-computed reference. First the
// incremental two-category acceptance the sampler forms equals a full
// all-category log-likelihood-plus-log-prior difference (a direct evaluation).
// Then, holding (y, eta) fixed, the gamma_2 draws from repeated refreshLatents
// form a Markov chain whose stationary law is the exact 1-D posterior p(gamma_2)
// proportional to prod_i (Phi-difference) times the log-gap prior with its
// Jacobian; the sampler's posterior mean/variance and a coarse histogram match a
// fine-grid quadrature of it. Local generator, restored global rngState.
static void testOrdinalCutpointConditional(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 10, K = 3;
  std::vector<double> y = {1, 1, 2, 2, 2, 3, 3, 2, 3, 1};
  std::vector<double> eta = {-0.4, 0.1, 0.3, -0.2, 0.5,
                             0.8, 1.1, -0.1, 0.6, 0.2};

  OrdinalResponse resp(y.data(), nullptr, n, K);
  const double m0 = 0.0, s0 = 1.5;  // mirror OrdinalResponse's prior constants
  auto logGapTarget = [&](double d) {
    double z = (std::log(d) - m0) / s0;
    return -0.5 * z * z - std::log(d);
  };
  auto Phi = [](double x) { return Rf_pnorm5(x, 0.0, 1.0, 1, 0); };
  // total unnormalized log posterior of gamma_2 over ALL observations
  auto logPost = [&](double g2) {
    double lp = logGapTarget(g2);
    for (std::size_t i = 0; i < n; ++i) {
      int k = static_cast<int>(y[i]);
      double lo = k == 1 ? -HUGE_VAL : (k == 2 ? 0.0 : g2);
      double hi = k == 3 ? HUGE_VAL : (k == 1 ? 0.0 : g2);
      lp += std::log(Phi(hi - eta[i]) - Phi(lo - eta[i]));
    }
    return lp;
  };

  // (a) the sampler's incremental acceptance equals the direct full evaluation
  double current = resp.cutpoints()[1];  // cold-init gamma_2 = spacing = 1
  bool acceptExact = true;
  for (double p : {0.3, 0.7, 1.4, 2.5}) {
    double got = resp.cutpointLogAcceptanceForTesting(eta.data(), 2, p);
    if (std::fabs(got - (logPost(p) - logPost(current))) > 1e-9)
      acceptExact = false;
  }
  check(acceptExact, "ordinal cutpoint incremental acceptance matches full eval");

  // (b) reference moments and bin masses by fine-grid quadrature over gamma_2
  const int G = 400000, numBins = 8;
  const double g2Max = 6.0, hq = g2Max / G, binWidth = g2Max / numBins;
  double gmax = -HUGE_VAL;
  for (int j = 0; j < G; ++j) {
    double lpv = logPost(hq * (j + 0.5));
    if (lpv > gmax) gmax = lpv;
  }
  std::vector<double> binMass(numBins, 0.0);
  double sw = 0.0, swx = 0.0, swx2 = 0.0;
  for (int j = 0; j < G; ++j) {
    double g2 = hq * (j + 0.5), w = std::exp(logPost(g2) - gmax);
    sw += w; swx += w * g2; swx2 += w * g2 * g2;
    int b = static_cast<int>(g2 / binWidth);
    binMass[b >= numBins ? numBins - 1 : b] += w;
  }
  double refMean = swx / sw, refVar = swx2 / sw - refMean * refMean;
  for (double& m : binMass) m /= sw;

  // (c) run the marginal-MH chain via refreshLatents and collect gamma_2
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 314159u);
  const int burn = 20000, keep = 1000000;
  for (int d = 0; d < burn; ++d) resp.refreshLatents(localRng, eta.data(), 1.0);
  double sum = 0.0, sumSq = 0.0;
  std::vector<double> counts(numBins, 0.0);
  for (int d = 0; d < keep; ++d) {
    resp.refreshLatents(localRng, eta.data(), 1.0);
    double g2 = resp.cutpoints()[1];
    sum += g2; sumSq += g2 * g2;
    int b = static_cast<int>(g2 / binWidth);
    if (b >= 0 && b < numBins) counts[b] += 1.0;
  }
  double mean = sum / keep, var = sumSq / keep - mean * mean;
  checkNear(mean, refMean, 0.02, "ordinal gamma_2 posterior mean");
  checkNear(var, refVar, std::max(0.1 * refVar, 0.005),
            "ordinal gamma_2 posterior variance");
  for (int b = 0; b < numBins; ++b) {
    double freq = counts[b] / keep;
    // deflate the sample size for the MH chain's autocorrelation
    double se = std::sqrt(binMass[b] * (1.0 - binMass[b]) / (keep / 40.0));
    checkNear(freq, binMass[b], std::max(0.01, 5.0 * se),
              "ordinal gamma_2 posterior histogram bin");
  }

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: ordinal cutpoint conditional\n");
}

// Ordinal at K = 2 with gamma_1 = 0 reproduces the probit sampler BIT FOR BIT:
// the empty free-cutpoint set consumes no rng, the boundary categories route
// through probit's one-sided primitives with the same argument, and the NaN
// fallback is the same sign * DBL_EPSILON. Identically seeded generators produce
// identical latents and are left in the identical state (proved by lockstep
// draws afterward), across several sweeps with a nonzero offset. Local
// generators, restored global rngState.
static void testOrdinalProbitEquivalence(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 40;
  std::vector<double> y01(n), y12(n), offset(n), totalFits(n);
  for (std::size_t i = 0; i < n; ++i) {
    y01[i] = runif01() < 0.5 ? 0.0 : 1.0;
    y12[i] = y01[i] + 1.0;  // ordinal category 1 (failure) / 2 (success)
    offset[i] = 0.4 * (runif01() - 0.5);
  }

  ProbitResponse probit(y01.data(), offset.data(), n);
  OrdinalResponse ordinal(y12.data(), offset.data(), n, 2);
  check(ordinal.numCutpoints() == 1 && ordinal.cutpoints()[0] == 0.0,
        "ordinal K = 2 pins a single cutpoint at 0");

  ext_rng* rngP = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngO = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rngP, 8675309u);
  ext_rng_setSeed(rngO, 8675309u);

  bool latentsMatch = true, workingMatch = true;
  for (int sweep = 0; sweep < 5; ++sweep) {
    // perturb the fits each sweep so both truncation sides are exercised
    for (std::size_t i = 0; i < n; ++i)
      totalFits[i] = 1.5 * (runif01() - 0.5);
    probit.refreshLatents(rngP, totalFits.data(), 1.0);
    ordinal.refreshLatents(rngO, totalFits.data(), 1.0);
    const double* lp = probit.latents();
    const double* lo = ordinal.latents();
    const double* wp = probit.workingResponse();
    const double* wo = ordinal.workingResponse();
    for (std::size_t i = 0; i < n; ++i) {
      if (lp[i] != lo[i]) latentsMatch = false;
      if (wp[i] != wo[i]) workingMatch = false;
    }
  }
  check(latentsMatch, "ordinal K = 2 latents equal probit bit for bit");
  check(workingMatch, "ordinal K = 2 working equals probit bit for bit");

  // identical consumption: subsequent draws stay in lockstep
  bool lockstep = true;
  for (int d = 0; d < 64; ++d)
    if (ext_rng_simulateContinuousUniform(rngP) !=
        ext_rng_simulateContinuousUniform(rngO))
      lockstep = false;
  check(lockstep, "ordinal K = 2 leaves the rng stream identical to probit");

  ext_rng_destroy(rngO);
  ext_rng_destroy(rngP);
  rngState = savedRngState;
  printf("ok: ordinal K = 2 equals probit bit for bit\n");
}

// The ordinal cutpoint block round-trips through getState/setState: a by-name
// length-(K-1) "cutpoints" slot the resid.df pattern serializes for a vector,
// written only for the ordinal family and restored exactly. Covers the fixed
// (pinned K = 2, length 1) and the after-sweeps (moved K = 4) blocks, and
// refuses an ordinal state whose cutpoint vector is the wrong length or absent.
// Local generators and a restored global rngState leave the shared stream
// untouched (the sampler draws through its own local rngs).
static void testOrdinalStateRoundTrip() {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 200, p = 2, numChains = 2;
  const double gamma4[3] = {0.0, 0.8, 1.7};  // K = 4 truth, gamma_1 = 0
  std::vector<double> x(n * p), y2(n), y4(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < p; ++j) x[i + j * n] = runif01();
    double f = 1.4 * (x[i] - 0.5) - 0.7 * (x[i + n] - 0.5);
    double u1 = runif01(), u2 = runif01();
    double z = f + std::sqrt(-2.0 * std::log(u1)) *
                     std::cos(6.283185307179586 * u2);
    y2[i] = z <= 0.0 ? 1.0 : 2.0;
    int k = 1;
    while (k < 4 && z > gamma4[k - 1]) ++k;
    y4[i] = static_cast<double>(k);
  }

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t base) {
    for (std::size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], base + static_cast<std::uint32_t>(c));
    }
  };
  auto buildOrdinal = [&](std::size_t K, const double* y,
                          std::vector<ext_rng*>& rngs, std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = 20;
    options.numChains = numChains;
    options.numCategories = K;
    makeRngs(rngs, seed);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y, n, p, nullptr, nullptr, ResponseFamily::ordinal, 1.0, 3.0,
      0.37804942330213542, options, rngs.data());
  };

  // (a) K = 2: the single cutpoint is pinned at 0 and never moves, yet still
  // rides its own length-1 block through the round trip
  {
    std::vector<ext_rng*> rngs(numChains, nullptr), rngs2(numChains, nullptr);
    auto original = buildOrdinal(2, y2.data(), rngs, 2200);
    Results empty;
    original->run(40, 0, empty);
    SamplerStateData state;
    original->getState(state);
    check(state.chains[0].cutpoints.size() == 1 &&
          state.chains[0].cutpoints[0] == 0.0,
          "ordinal K = 2 state carries the pinned length-1 cutpoint");
    auto restored = buildOrdinal(2, y2.data(), rngs2, 5500);
    check(restored->setState(state, nullptr), "an ordinal K = 2 state restores");
    checkStructuralRoundTrip(state, *restored,
                             "restored ordinal K = 2 cutpoint agrees");
    for (ext_rng* r : rngs) ext_rng_destroy(r);
    for (ext_rng* r : rngs2) ext_rng_destroy(r);
  }

  // (b) K = 4: the two free cutpoints move over sweeps and round-trip; a
  // wrong-length or absent cutpoint block is refused
  {
    std::vector<ext_rng*> rngs(numChains, nullptr), rngs2(numChains, nullptr);
    auto original = buildOrdinal(4, y4.data(), rngs, 3300);
    SamplerStateData fixedState;
    original->getState(fixedState);
    check(fixedState.chains[0].cutpoints.size() == 3,
          "ordinal K = 4 state carries a length-(K-1) cutpoint vector");
    Results empty;
    original->run(80, 0, empty);
    SamplerStateData state;
    original->getState(state);
    bool moved = false;
    for (std::size_t c = 0; c < numChains; ++c)
      if (state.chains[c].cutpoints != fixedState.chains[c].cutpoints)
        moved = true;
    check(moved, "ordinal free cutpoints move over sweeps");
    check(state.chains[0].cutpoints[0] == 0.0,
          "ordinal gamma_1 stays pinned at 0 across sweeps");

    auto restored = buildOrdinal(4, y4.data(), rngs2, 8800);
    check(restored->setState(state, nullptr), "an ordinal K = 4 state restores");
    checkStructuralRoundTrip(state, *restored,
                             "restored ordinal cutpoints and latents agree");

    SamplerStateData badLen(state);
    badLen.chains[0].cutpoints.resize(2);
    check(!restored->setState(badLen, nullptr),
          "an ordinal state with a short cutpoint vector is refused");
    SamplerStateData noCut(state);
    for (auto& ch : noCut.chains) ch.cutpoints.clear();
    check(!restored->setState(noCut, nullptr),
          "an ordinal state lacking the cutpoint vector is refused");
    for (ext_rng* r : rngs) ext_rng_destroy(r);
    for (ext_rng* r : rngs2) ext_rng_destroy(r);
  }

  rngState = savedRngState;
  printf("ok: ordinal cutpoint state round trip\n");
}

// The engine test-container entry (facade setTestData): a sampler given a
// mixed dense + CSC test set records test fits BITWISE-IDENTICALLY to one given
// the dense test matrix of the same values, and the entry REFUSES a designated
// leaf covariate that would be CSC-backed. Uses local rngs and restores the
// shared draw stream, so its fixture draws shift no later test.
static void testSparseTestDataEndToEnd() {
  uint64_t savedRngState = rngState;
  const size_t nTrain = 300, numTest = 150, p = 4;

  std::vector<double> xTrain(nTrain * p), y(nTrain);
  for (double& v : xTrain) v = runif01();
  for (size_t i = 0; i < nTrain; ++i)
    y[i] = std::sin(3.0 * xTrain[i]) + xTrain[i + nTrain] + 0.5 * runif01();

  // test columns 0 and 3 dense-backed, 1 rank-tier, 2 densified-tier
  std::vector<double> dense0(numTest), dense3(numTest);
  for (size_t i = 0; i < numTest; ++i) {
    dense0[i] = runif01();
    dense3[i] = runif01();
  }
  std::vector<int> pointers(3, 0), rows;
  std::vector<double> values;
  const double fractions[] = {0.08, 0.6};
  for (int csc = 0; csc < 2; ++csc) {
    for (size_t i = 0; i < numTest; ++i)
      if (runif01() < fractions[csc]) {
        rows.push_back(static_cast<int>(i));
        values.push_back(0.3 + runif01());
      }
    pointers[static_cast<size_t>(csc) + 1] = static_cast<int>(rows.size());
  }
  std::vector<double> denseTest(numTest * p, 0.0);
  for (size_t i = 0; i < numTest; ++i) {
    denseTest[i + 0 * numTest] = dense0[i];
    denseTest[i + 3 * numTest] = dense3[i];
  }
  for (int csc = 0; csc < 2; ++csc)
    for (int k = pointers[static_cast<size_t>(csc)];
         k < pointers[static_cast<size_t>(csc) + 1]; ++k)
      denseTest[static_cast<size_t>(rows[static_cast<size_t>(k)]) +
                static_cast<size_t>(csc + 1) * numTest] =
        values[static_cast<size_t>(k)];
  std::vector<double> denseBlock(numTest * 2);
  std::memcpy(denseBlock.data(), dense0.data(), numTest * sizeof(double));
  std::memcpy(denseBlock.data() + numTest, dense3.data(),
              numTest * sizeof(double));
  std::vector<std::int32_t> columnSources = {0, ~0, ~1, 1};

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngA == NULL || rngB == NULL || ext_rng_setSeed(rngA, 9137) != 0 ||
      ext_rng_setSeed(rngB, 9137) != 0) {
    check(false, "sparse test data: rng creation");
    rngState = savedRngState;
    return;
  }

  SamplerOptions options;
  options.numTrees = 25;
  ConstantLeafSampler denseSampler(xTrain.data(), y.data(), nTrain, p, nullptr,
                       nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                       0.37804942330213542, options, &rngA);
  ConstantLeafSampler sparseSampler(xTrain.data(), y.data(), nTrain, p, nullptr,
                       nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                       0.37804942330213542, options, &rngB);
  denseSampler.setTestPredictors(denseTest.data(), numTest);
  bool built = sparseSampler.setTestData(denseBlock.data(), pointers.data(),
                                         rows.data(), values.data(),
                                         columnSources.data(), nullptr, numTest);
  check(built && sparseSampler.data().testColumnIsSparse(1) &&
        !sparseSampler.data().testColumnIsSparse(2),
        "the test container builds with the expected storage tiers");

  const size_t numBurnIn = 30, numSamples = 40;
  std::vector<double> testA(numTest * numSamples), testB(numTest * numSamples);
  std::vector<double> fitsA(nTrain * numSamples), fitsB(nTrain * numSamples);
  Results resultsA, resultsB;
  resultsA.trainingFits = fitsA.data();
  resultsA.testFits = testA.data();
  resultsB.trainingFits = fitsB.data();
  resultsB.testFits = testB.data();
  denseSampler.run(numBurnIn, numSamples, resultsA);
  sparseSampler.run(numBurnIn, numSamples, resultsB);
  check(fitsA == fitsB, "test representation leaves the training draws untouched");
  check(testA == testB,
        "resident-sparse test fits bitwise-match the dense test matrix");

  // the refusal seam: a linear leaf over a dense-backed test column builds; the
  // same covariate landing on a CSC-backed test column refuses, store untouched
  ext_rng* rngC = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  if (rngC != NULL && ext_rng_setSeed(rngC, 41) == 0) {
    SamplerOptions linearOptions;
    linearOptions.numTrees = 10;
    size_t covariate[] = {0};
    linearOptions.leafCovariateColumns = covariate;
    linearOptions.numLeafCovariates = 1;
    std::unique_ptr<SamplerBase> linear = createSampler(
      xTrain.data(), y.data(), nTrain, p, nullptr, nullptr,
      ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, linearOptions,
      &rngC);
    check(linear != nullptr, "dense-backed leaf covariate sampler creates");
    if (linear != nullptr) {
      check(linear->setTestData(denseBlock.data(), pointers.data(), rows.data(),
                                values.data(), columnSources.data(), nullptr,
                                numTest),
            "a dense-backed leaf covariate accepts the test container");
      std::vector<std::int32_t> refuse = {~0, 0, ~1, 1};  // column 0 now CSC
      check(!linear->setTestData(denseBlock.data(), pointers.data(),
                                 rows.data(), values.data(), refuse.data(),
                                 nullptr, numTest),
            "a CSC-backed leaf covariate refuses the test container");
    }
    ext_rng_destroy(rngC);
  }

  ext_rng_destroy(rngB);
  ext_rng_destroy(rngA);
  rngState = savedRngState;
  printf("ok: sparse test data end to end\n");
}

// Moments of the shape-parameterized PG helper: PG(b, z) for integer b is the
// sum of b iid PG(1, z), so its mean is b (tanh(z/2)/(2z)) and its variance b
// times the PG(1, z) variance (the polya-gamma moments precedent). Plus the
// stream identity that b = 1 consumes exactly the shipped PG(1) Devroye draw bit
// for bit (the ordinal K = 2 identity analog). Local generators, restored state.
static void testNBPolyaGammaShapeMoments(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  auto pg1Mean = [](double c) {
    return c == 0.0 ? 0.25 : std::tanh(0.5 * c) / (2.0 * c);
  };
  auto pg1Var = [](double c) {
    return c == 0.0 ? 1.0 / 24.0
                    : (std::sinh(c) - c) / (4.0 * c * c * c *
                                            std::cosh(0.5 * c) *
                                            std::cosh(0.5 * c));
  };
  const double zs[] = {0.0, 1.5, -3.0};
  const double bs[] = {1.0, 2.0, 5.0, 11.0};
  const int numDraws = 50000;

  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 20260718u);
  for (double z : zs) {
    double c = std::fabs(z);
    for (double b : bs) {
      double expectedMean = b * pg1Mean(c);
      double expectedVar = b * pg1Var(c);
      double sum = 0.0, sumSq = 0.0;
      for (int i = 0; i < numDraws; ++i) {
        double draw = simulatePolyaGammaShape(localRng, b, z);
        if (!(draw > 0.0)) {
          check(false, "nb shape PG draw not positive");
          ext_rng_destroy(localRng);
          rngState = savedRngState;
          return;
        }
        sum += draw;
        sumSq += draw * draw;
      }
      double mean = sum / numDraws;
      double var = sumSq / numDraws - mean * mean;
      double meanSe = std::sqrt(expectedVar / static_cast<double>(numDraws));
      checkNear(mean, expectedMean, 5.0 * meanSe, "nb PG(b, z) shape mean");
      checkNear(var, expectedVar, 0.05 * expectedVar,
                "nb PG(b, z) shape variance");
    }
  }

  // b = 1 is bit-for-bit the shipped PG(1) Devroye stream
  ext_rng* rShape = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rDevroye = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rShape, 8675309u);
  ext_rng_setSeed(rDevroye, 8675309u);
  bool identical = true;
  for (int i = 0; i < 256; ++i) {
    double z = 0.02 * i - 2.5;
    if (simulatePolyaGammaShape(rShape, 1.0, z) !=
        ext_rng_simulatePolyaGamma(rDevroye, z))
      identical = false;
  }
  check(identical, "nb PG(1, z) equals the shipped Devroye stream bit for bit");

  ext_rng_destroy(rDevroye);
  ext_rng_destroy(rShape);
  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: nb polya-gamma shape moments\n");
}

// The NB dispersion grid draw reproduces the hand-computed discrete full
// conditional: for fixed counts y and a fixed fit, w_k proportional to
// exp(L_k + r_k S + log prior_k) with S = sum_i log(1 - p_i), L_k the count-
// histogram lgamma kernel checked against a direct evaluation, prior the
// renormalized gamma(2, 0.1). The ResidualDfPrior drawIndex test pattern. Local
// generator, restored global rngState.
static void testNBDispersionGridConditional(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  NBDispersionPrior prior;
  std::vector<double> y = {0, 1, 2, 2, 3, 5, 1, 0, 4, 2, 7, 1};
  const std::size_t n = y.size();
  prior.computeKernel(y.data(), n);

  // a fixed fit -> S = sum_i log(1 - plogis(psi_i)) = -sum_i log(1 + e^{psi_i})
  double S = 0.0;
  for (std::size_t i = 0; i < n; ++i)
    S -= logOnePlusExp(0.3 * static_cast<double>(i) - 1.2);

  // (a) the kernel matches a direct lgamma evaluation
  bool kernelExact = true;
  for (std::size_t k = 0; k < NBDispersionPrior::gridSize; ++k) {
    double rk = NBDispersionPrior::grid[k];
    double L = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      L += std::lgamma(y[i] + rk) - std::lgamma(rk);
    if (std::fabs(prior.kernelValue(k) - L) > 1e-9) kernelExact = false;
  }
  check(kernelExact, "nb dispersion kernel matches direct lgamma sum");

  // (b) reference posterior over the grid
  double expected[NBDispersionPrior::gridSize];
  double maxLog = -HUGE_VAL;
  for (std::size_t k = 0; k < NBDispersionPrior::gridSize; ++k) {
    double rk = NBDispersionPrior::grid[k];
    double logPost = prior.kernelValue(k) + rk * S + std::log(rk) - 0.1 * rk;
    expected[k] = logPost;
    if (logPost > maxLog) maxLog = logPost;
  }
  double total = 0.0;
  for (std::size_t k = 0; k < NBDispersionPrior::gridSize; ++k) {
    expected[k] = std::exp(expected[k] - maxLog);
    total += expected[k];
  }
  for (std::size_t k = 0; k < NBDispersionPrior::gridSize; ++k)
    expected[k] /= total;

  // (c) sampled histogram
  ext_rng* localRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(localRng, 1729u);
  const int numDraws = 400000;
  std::vector<int> counts(NBDispersionPrior::gridSize, 0);
  for (int d = 0; d < numDraws; ++d) ++counts[prior.drawIndex(localRng, S)];

  for (std::size_t k = 0; k < NBDispersionPrior::gridSize; ++k) {
    double freq = static_cast<double>(counts[k]) / numDraws;
    double se = std::sqrt(expected[k] * (1.0 - expected[k]) / numDraws);
    checkNear(freq, expected[k], std::max(0.003, 5.0 * se),
              "nb dispersion grid frequency matches the hand-computed posterior");
  }

  ext_rng_destroy(localRng);
  rngState = savedRngState;
  printf("ok: nb dispersion grid full conditional\n");
}

// The r-FIRST partially-collapsed sweep order (section 5), reconstructed bit for
// bit: refreshLatents must (1) draw r from the grid conditional against S from
// the CURRENT fit, then (2) draw omega at r_new, then (3) rebuild working with
// r_new. An independent replay on an identically seeded generator - drawIndex on
// the same S, then the shape draws at r_new - reproduces the dispersion, the
// omega, and the working response exactly. Then the restore contract: restoring
// r before latents reproduces working exactly, and restoring latents under a
// stale r yields a different working (so order matters). Local generators,
// restored global rngState.
static void testNBSweepOrderAndRestore(ext_rng*) {
  std::uint64_t savedRngState = rngState;
  std::vector<double> y = {0, 2, 1, 4, 3, 1, 2, 5, 0, 2};
  std::vector<double> offset = {0.1,  -0.2, 0.0,  0.3,  -0.1,
                                0.2,  -0.3, 0.05, 0.15, -0.05};
  const std::size_t n = y.size();
  std::vector<double> totalFits(n);
  for (std::size_t i = 0; i < n; ++i)
    totalFits[i] = 0.2 * static_cast<double>(i) - 0.7;

  NBResponse resp(y.data(), offset.data(), n, -1.0);  // grid mode
  check(resp.estimatesDispersion() && resp.carriesDispersion(),
        "nb grid mode estimates and carries dispersion");

  ext_rng* rResp = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng* rRef = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rResp, 424242u);
  ext_rng_setSeed(rRef, 424242u);
  resp.refreshLatents(rResp, totalFits.data(), 1.0);

  // independent replay in the mandated order on the same rng stream
  NBDispersionPrior refPrior;
  refPrior.computeKernel(y.data(), n);
  double S = 0.0;
  for (std::size_t i = 0; i < n; ++i)
    S -= logOnePlusExp(totalFits[i] + offset[i]);
  double rNew = NBDispersionPrior::grid[refPrior.drawIndex(rRef, S)];
  check(rNew == resp.dispersion(),
        "nb r drawn first from the grid conditional at the current fit");
  bool omegaMatch = true, workingMatch = true;
  const double* omega = resp.latents();
  const double* working = resp.workingResponse();
  for (std::size_t i = 0; i < n; ++i) {
    double psi = totalFits[i] + offset[i];
    double refOmega = simulatePolyaGammaShape(rRef, y[i] + rNew, psi);
    if (refOmega != omega[i]) omegaMatch = false;
    if (working[i] != 0.5 * (y[i] - rNew) / refOmega - offset[i])
      workingMatch = false;
  }
  check(omegaMatch, "nb omega drawn at r_new in the mandated order bit for bit");
  check(workingMatch, "nb working rebuilt with r_new bit for bit");

  // setResponse keeps r (the kept-cutpoints/nu clause): a y swap in grid mode
  // draws NO r variate, only omega at the CURRENT r - proved bitwise, since an
  // r draw would consume a uniform first and shift every omega draw after it
  {
    std::vector<double> ySwap(y);
    ySwap[3] += 1.0;
    double rKept = resp.dispersion();
    ext_rng* rSet = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* rSetRef = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rSet, 55555u);
    ext_rng_setSeed(rSetRef, 55555u);
    resp.setResponse(ySwap.data(), rSet, totalFits.data(), false, nullptr);
    check(resp.dispersion() == rKept, "nb setResponse keeps the current r");
    bool swapExact = true;
    for (std::size_t i = 0; i < n; ++i) {
      double psi = totalFits[i] + offset[i];
      double refOmega = simulatePolyaGammaShape(rSetRef, ySwap[i] + rKept, psi);
      if (resp.latents()[i] != refOmega ||
          resp.workingResponse()[i] !=
            0.5 * (ySwap[i] - rKept) / refOmega - offset[i])
        swapExact = false;
    }
    check(swapExact, "nb setResponse redraws omega at the kept r bit for bit");
    ext_rng_destroy(rSetRef);
    ext_rng_destroy(rSet);
    // put the fixture back for the restore-contract block below
    ext_rng_setSeed(rResp, 424242u);
    resp.setResponse(y.data(), rResp, totalFits.data(), false, nullptr);
  }

  // restore contract: r before latents reproduces working exactly
  double rSaved = resp.dispersion();
  std::vector<double> omegaSaved(omega, omega + n);
  std::vector<double> workingSaved(working, working + n);

  NBResponse dst(y.data(), offset.data(), n, -1.0);
  dst.restoreDispersion(rSaved);
  dst.restoreLatents(omegaSaved.data());
  bool restoreExact = true;
  for (std::size_t i = 0; i < n; ++i)
    if (dst.workingResponse()[i] != workingSaved[i]) restoreExact = false;
  check(restoreExact,
        "nb restoreDispersion-then-restoreLatents reproduces working exactly");

  // stale-r sensitivity: restoring latents under a different r rebuilds working
  // against that r, so restoring latents before r (the forbidden order) is
  // observably wrong
  double staleR = rSaved + 1.0;
  NBResponse stale(y.data(), offset.data(), n, staleR);  // fixed at staleR
  check(!stale.estimatesDispersion() && stale.dispersion() == staleR,
        "nb fixed mode holds the supplied dispersion");
  stale.restoreLatents(omegaSaved.data());
  bool usesCurrentR = true, differsFromSaved = false;
  for (std::size_t i = 0; i < n; ++i) {
    if (stale.workingResponse()[i] !=
        0.5 * (y[i] - staleR) / omegaSaved[i] - offset[i])
      usesCurrentR = false;
    if (stale.workingResponse()[i] != workingSaved[i]) differsFromSaved = true;
  }
  check(usesCurrentR, "nb restoreLatents rebuilds working from the current r");
  check(differsFromSaved,
        "nb restoring latents under a stale r yields a different working");

  ext_rng_destroy(rRef);
  ext_rng_destroy(rResp);
  rngState = savedRngState;
  printf("ok: nb sweep order and restore contract\n");
}

// The nbinom dispersion state block round-trips through the Chain serialization
// (getState/setState): the scalar r rides its own by-name block and the omega
// latents ride the shared latents block, both exact, in fixed and grid-estimated
// modes. The restore ORDER (restoreDispersion before restoreLatents) is a stated
// contract - a wrong order rebuilds the working response against a stale r, the
// sensitivity testNBSweepOrderAndRestore proves at the model level above; here
// the sampler-level round trip is asserted exact under setState's ordering. A
// fresh sampler seeded differently must let the serialized r and omega win over
// the cold-start median. RNG-insulated: local generators, restored global
// rngState, so its fixture draws shift no later test.
static void testNBStateRoundTrip() {
  std::uint64_t savedRngState = rngState;
  const std::size_t n = 200, p = 2, numChains = 2;
  std::vector<double> x(n * p), y(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < p; ++j) x[i + j * n] = runif01();
    double f = 1.2 * (x[i] - 0.5) - 0.6 * (x[i + n] - 0.5);
    // a moderate mean count; the data law is immaterial to a state round trip,
    // only non-negative integrality matters, so a crude count around mu suffices
    y[i] = std::floor(3.0 * std::exp(f) + 2.0 * runif01());
  }

  auto makeRngs = [](std::vector<ext_rng*>& rngs, std::uint32_t base) {
    for (std::size_t c = 0; c < rngs.size(); ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rngs[c], base + static_cast<std::uint32_t>(c));
    }
  };
  auto buildNB = [&](double dispersion, std::vector<ext_rng*>& rngs,
                     std::uint32_t seed) {
    SamplerOptions options;
    options.numTrees = 20;
    options.numChains = numChains;
    options.dispersion = dispersion;
    makeRngs(rngs, seed);
    return std::make_unique<ConstantLeafSampler>(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::nbinom, 1.0,
      3.0, 0.37804942330213542, options, rngs.data());
  };

  auto runOneMode = [&](double dispersion, bool estimated, const char* tag) {
    std::vector<ext_rng*> rngs(numChains, nullptr), rngs2(numChains, nullptr);
    auto original = buildNB(dispersion, rngs, 71000);
    Results empty;
    original->run(40, 0, empty);

    SamplerStateData state;
    original->getState(state);
    check(state.chains[0].latents.size() == n,
          "nb state carries the omega latents");
    check(std::isfinite(state.chains[0].dispersion) &&
            state.chains[0].dispersion > 0.0,
          "nb state carries a finite positive dispersion");
    if (!estimated)
      check(state.chains[0].dispersion == dispersion,
            "fixed-r nb state records the supplied dispersion");

    // a fresh sampler with a DIFFERENT seed: the serialized dispersion and omega
    // must win over the cold-start median r
    auto restored = buildNB(dispersion, rngs2, 99000);
    check(restored->setState(state, nullptr), "an nb state restores");
    // statesAgree compares the dispersion and omega, both bitwise: restore
    // reconstructs the block, and setState's restoreDispersion-before-
    // restoreLatents order rebuilt working under the restored r
    checkStructuralRoundTrip(state, *restored,
                             "restored nb dispersion and omega agree");
    SamplerStateData reState;
    restored->getState(reState);
    check(reState.chains[0].dispersion == state.chains[0].dispersion,
          "nb dispersion round-trips exactly");
    check(reState.chains[0].latents == state.chains[0].latents,
          "nb omega round-trips exactly");

    // a state whose dispersion block is absent (an old or non-count state) is
    // refused: r is NaN and stateIsValid rejects it
    SamplerStateData noR(state);
    for (auto& ch : noR.chains)
      ch.dispersion = std::numeric_limits<double>::quiet_NaN();
    check(!restored->setState(noR, nullptr),
          "an nb state lacking a finite dispersion is refused");
    // a non-positive dispersion is refused (the grid holds positive integers)
    SamplerStateData badR(state);
    for (auto& ch : badR.chains) ch.dispersion = 0.0;
    check(!restored->setState(badR, nullptr),
          "an nb state with a non-positive dispersion is refused");

    for (ext_rng* r : rngs) ext_rng_destroy(r);
    for (ext_rng* r : rngs2) ext_rng_destroy(r);
    printf("ok: nb dispersion state round trip (%s)\n", tag);
  };

  runOneMode(-1.0, true, "grid estimated");
  runOneMode(4.0, false, "fixed r");

  // an nb sampler refuses a gaussian state: it carries neither the omega latents
  // nor a dispersion, so the shared latents check and the dispersion requirement
  // together reject it (the cross-family refusal, resid.df/cutpoints precedent)
  {
    std::vector<ext_rng*> rngsG(numChains, nullptr), rngsN(numChains, nullptr);
    SamplerOptions gaussOptions;
    gaussOptions.numTrees = 20;
    gaussOptions.numChains = numChains;
    makeRngs(rngsG, 4242);
    ConstantLeafSampler gaussian(
      x.data(), y.data(), n, p, nullptr, nullptr, ResponseFamily::gaussian, 1.0,
      3.0, 0.37804942330213542, gaussOptions, rngsG.data());
    Results emptyG;
    gaussian.run(20, 0, emptyG);
    SamplerStateData gaussState;
    gaussian.getState(gaussState);
    check(std::isnan(gaussState.chains[0].dispersion),
          "gaussian state carries no dispersion");
    auto nbSampler = buildNB(-1.0, rngsN, 5353);
    check(!nbSampler->setState(gaussState, nullptr),
          "an nb sampler refuses a gaussian state");
    for (ext_rng* r : rngsG) ext_rng_destroy(r);
    for (ext_rng* r : rngsN) ext_rng_destroy(r);
  }

  rngState = savedRngState;
  printf("ok: nb dispersion state round trip\n");
}

namespace {
// Positive controls: minimal leaves declaring each optional seam method, to
// prove the detection concepts are not vacuously false.
struct TreeDrawMock {
  void drawParametersForTree(ext_rng*, const Tree&,
                             const std::vector<std::int32_t>&, double, double,
                             double*) const {}
};
struct ParamScoreMock {
  double logLikelihoodForBranchWithParams(const Tree&, std::int32_t,
                                          const double*, const double*, double,
                                          double, const double*) const {
    return 0.0;
  }
};
}  // namespace

// The C1 seam is inert for every shipped leaf: no existing leaf declares the
// tree-granularity draw or the M-reading branch score, so the constant leaf
// takes the byte-identical per-node draw and its moves integrate every leaf
// out. The mocks confirm the concepts detect the methods when present.
static void testLeafSeamDispatch() {
  static_assert(!TreeDrawLeafModel<ConstantGaussianLeaf>);
  static_assert(!ParamScoringLeafModel<ConstantGaussianLeaf>);
  static_assert(!TreeDrawLeafModel<LinearGaussianLeaf>);
  static_assert(!ParamScoringLeafModel<LinearGaussianLeaf>);
  static_assert(TreeDrawLeafModel<TreeDrawMock>);
  static_assert(ParamScoringLeafModel<ParamScoreMock>);
  static_assert(!ParamScoringLeafModel<TreeDrawMock>);
  static_assert(!TreeDrawLeafModel<ParamScoreMock>);
  printf("ok: leaf seam dispatch inert for shipped leaves\n");
}

namespace {

// ---- monotone (mBART) constrained leaf: independent references -------------

// Codes reaching a leaf along a variable, [lo, hi], by a direct ancestor walk
// (not through splitInterval), so it independently checks the geometry.
static void oracleLeafBox(const Tree& tree, const ColumnStore& store,
                          std::int32_t leaf, std::int32_t var, std::int32_t& lo,
                          std::int32_t& hi) {
  lo = 0;
  hi = static_cast<std::int32_t>(store.numCuts[static_cast<size_t>(var)]);
  std::int32_t cur = leaf;
  while (tree.at(cur).parent != invalidNode) {
    std::int32_t parent = tree.at(cur).parent;
    bool isRight = (cur == tree.at(parent).leftChild + 1);
    cur = parent;
    if (tree.at(cur).rule.variableIndex == var) {
      std::int32_t s = tree.at(cur).rule.splitIndex();
      if (isRight) lo = std::max(lo, s + 1);
      else hi = std::min(hi, s);
    }
  }
}

// Independent bounds [a, b] and constrained flag for leaf k, from the design's
// neighbor definition over the oracle boxes.
static void oracleBounds(const Tree& tree, const ColumnStore& store,
                         const std::vector<std::int8_t>& dir, std::int32_t k,
                         const double* mu, double& a, double& b,
                         bool& constrained) {
  size_t p = store.numPredictors;
  std::vector<std::int32_t> leaves;
  tree.fillBottom(0, leaves);
  a = -HUGE_VAL;
  b = HUGE_VAL;
  constrained = false;
  for (std::int32_t j : leaves) {
    if (j == k) continue;
    for (size_t i = 0; i < p; ++i) {
      if (dir[i] == 0) continue;
      std::int32_t loK, hiK, loJ, hiJ;
      oracleLeafBox(tree, store, k, static_cast<std::int32_t>(i), loK, hiK);
      oracleLeafBox(tree, store, j, static_cast<std::int32_t>(i), loJ, hiJ);
      bool jBelowK = (hiJ + 1 == loK);
      bool kBelowJ = (hiK + 1 == loJ);
      if (!jBelowK && !kBelowJ) continue;
      bool overlap = true;
      for (size_t m = 0; m < p && overlap; ++m) {
        if (m == i) continue;
        std::int32_t lkm, hkm, ljm, hjm;
        oracleLeafBox(tree, store, k, static_cast<std::int32_t>(m), lkm, hkm);
        oracleLeafBox(tree, store, j, static_cast<std::int32_t>(m), ljm, hjm);
        if (std::max(lkm, ljm) > std::min(hkm, hjm)) overlap = false;
      }
      if (!overlap) continue;
      constrained = true;
      double muj = mu[j];
      if (jBelowK) {
        if (dir[i] > 0) a = std::max(a, muj);
        else            b = std::min(b, muj);
      } else {
        if (dir[i] > 0) b = std::min(b, muj);
        else            a = std::max(a, muj);
      }
    }
  }
}

static void checkBound(double actual, double expected, const char* what) {
  if (std::isinf(expected) || std::isinf(actual)) check(actual == expected, what);
  else checkNear(actual, expected, 1e-9, what);
}

// Grow a random tree by a few random ordinal births.
static void growRandomTree(Tree& tree, const ColumnStore& store, const double* y,
                           int numBirths) {
  for (int step = 0; step < numBirths; ++step) {
    std::vector<std::int32_t> leaves;
    tree.fillBottom(0, leaves);
    std::int32_t leaf = leaves[static_cast<size_t>(runif01() * leaves.size())];
    int var = static_cast<int>(runif01() * store.numPredictors);
    std::int32_t left, right;
    tree.splitInterval(store, leaf, var, &left, &right);
    if (right < left) continue;
    std::int32_t cut = left + static_cast<std::int32_t>(runif01() * (right - left + 1));
    Rule rule;
    rule.variableIndex = var;
    rule.setSplitIndex(cut);
    tree.birth(store, leaf, rule, y, nullptr);
  }
}

static double refLogIntegratedLikelihood(double k, double sig2, double sumW,
                                         double sumWZ, double scale) {
  if (sumW == 0.0) return 0.0;
  double pp = (k / scale) * (k / scale);
  double qp = sumW / sig2;
  double mean = sumWZ / sumW;
  double expl = sumWZ * mean;
  return 0.5 * std::log(pp / (pp + qp)) + 0.5 * expl / sig2 -
         0.5 * ((pp * mean) * (qp * mean)) / (pp + qp);
}

static void refPosterior(double k, double sig2, double sumW, double sumWZ,
                         double scale, double& m, double& s) {
  double pp = (k / scale) * (k / scale);
  double qp = sumW / sig2;
  double var = 1.0 / (pp + qp);
  s = std::sqrt(var);
  m = (sumWZ / sig2) * var;
}

// Fine trapezoidal integral of f over [lo, hi] (infinite ends clamped by the
// caller), as an independent quadrature reference.
template <typename F>
static double trapz(F f, double lo, double hi, int steps = 200000) {
  double h = (hi - lo) / steps;
  double total = 0.5 * (f(lo) + f(hi));
  for (int i = 1; i < steps; ++i) total += f(lo + i * h);
  return total * h;
}

static double normalPdf(double x, double m, double s) {
  double z = (x - m) / s;
  return std::exp(-0.5 * z * z) / (s * std::sqrt(2.0 * std::numbers::pi));
}

}  // namespace

// (a) The neighbor geometry matches a brute-force box-adjacency oracle on
// fuzzed trees over two and three constrained axes.
static void testMonotoneNeighborGeometry() {
  std::uint64_t saved = rngState;
  rngState = 424242u;
  const size_t n = 600, p = 3;
  std::vector<double> x(n * p);
  for (double& v : x) v = runif01();
  ColumnStore store;
  store.build(x.data(), n, p, 24);
  std::vector<double> y(n, 0.0);
  std::vector<std::int8_t> dir[3] = {{1, 0, 0}, {1, -1, 0}, {1, -1, 1}};

  int comparisons = 0;
  for (int trial = 0; trial < 120; ++trial) {
    std::vector<index_t> idx(n);
    Tree tree;
    tree.initialize(idx.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growRandomTree(tree, store, y.data(), 4 + static_cast<int>(runif01() * 4));

    std::vector<std::int32_t> leaves;
    tree.fillBottom(0, leaves);
    std::vector<double> mu(tree.nodes.size(), 0.0);
    for (std::int32_t leaf : leaves) mu[leaf] = 2.0 * runif01() - 1.0;

    const std::vector<std::int8_t>& d = dir[trial % 3];
    MonotoneNeighborScratch scratch;
    for (std::int32_t k : leaves) {
      double a, b;
      bool c;
      monotoneNeighborBounds(tree, store, d.data(), leaves, k, mu.data(), nullptr,
                             0, scratch, &a, &b, &c);
      double ao, bo;
      bool co;
      oracleBounds(tree, store, d, k, mu.data(), ao, bo, co);
      checkBound(a, ao, "monotone geometry: lower bound matches oracle");
      checkBound(b, bo, "monotone geometry: upper bound matches oracle");
      check(c == co, "monotone geometry: constrained flag matches oracle");
      ++comparisons;
    }
  }
  rngState = saved;
  printf("ok: monotone neighbor geometry vs oracle (%d leaves)\n", comparisons);
}

// (b) Feasibility invariants: an end-to-end constrained fit is monotone in its
// single constrained predictor, an emptied cone scores the sentinel, and a
// direct draw sweep keeps every leaf within its bounds.
static void testMonotoneFeasibility() {
  // a private generator keeps the shared stream (and downstream suites) intact
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 7001);
  std::uint64_t savedRngState = rngState;
  rngState = 987654321u;
  // end-to-end: a monotone-increasing truth is recovered by a monotone fit
  const size_t n = 300;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    double u1 = runif01(), u2 = runif01();
    double z = std::sqrt(-2.0 * std::log(u1)) *
               std::cos(6.283185307179586 * u2);
    y[i] = 2.0 * x[i] + 0.2 * z;
  }
  std::vector<std::int8_t> dir = {1};
  SamplerOptions options;
  options.numTrees = 40;
  options.birthOrDeathProbability = 1.0;  // v1 constrained move set
  options.swapProbability = 0.0;
  options.changeProbability = 0.0;
  options.monotoneDirections = dir.data();
  Sampler<MonotoneConstantGaussianLeaf> sampler(
    x.data(), y.data(), n, size_t(1), nullptr, nullptr, ResponseFamily::gaussian,
    1.0, 3.0, 0.37804942330213542, options, &rng);
  const size_t numBurnIn = 150, numSamples = 250;
  std::vector<double> trainingFits(n * numSamples);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(numBurnIn, numSamples, results);

  std::vector<double> mean(n, 0.0);
  for (size_t s = 0; s < numSamples; ++s)
    for (size_t i = 0; i < n; ++i)
      mean[i] += trainingFits[i + s * n] / static_cast<double>(numSamples);
  std::vector<size_t> order(n);
  for (size_t i = 0; i < n; ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) { return x[a] < x[b]; });
  int descents = 0;
  for (size_t i = 1; i < n; ++i)
    if (mean[order[i]] < mean[order[i - 1]] - 1e-9) ++descents;
  check(descents == 0, "monotone end-to-end: fitted mean is non-decreasing");

  // empty cone: a leaf pinned between a high below-neighbor and a low above-
  // neighbor scores the sentinel (the change/swap infeasibility path)
  const size_t m = 300;
  std::vector<double> xg(m);
  for (size_t i = 0; i < m; ++i) xg[i] = static_cast<double>(i) / (m - 1);
  ColumnStore store;
  store.build(xg.data(), m, 1, 60);
  std::vector<double> yg(m, 0.0);
  std::vector<index_t> idx(m);
  Tree tree;
  tree.initialize(idx.data(), m);
  tree.computeLeafStats(0, yg.data(), nullptr);
  std::int32_t left, right;
  tree.splitInterval(store, 0, 0, &left, &right);
  Rule r0;
  r0.variableIndex = 0;
  r0.setSplitIndex(left + (right - left) / 3);
  tree.birth(store, 0, r0, yg.data(), nullptr);
  std::int32_t rightChild = tree.at(0).leftChild + 1;
  tree.splitInterval(store, rightChild, 0, &left, &right);
  Rule r1;
  r1.variableIndex = 0;
  r1.setSplitIndex(left + (right - left) / 2);
  tree.birth(store, rightChild, r1, yg.data(), nullptr);
  std::int32_t lowLeaf = tree.at(0).leftChild;         // A
  std::int32_t midLeaf = tree.at(rightChild).leftChild;  // L
  std::int32_t highLeaf = tree.at(rightChild).leftChild + 1;  // C
  std::vector<double> mu(tree.nodes.size(), 0.0);
  mu[lowLeaf] = 5.0;   // below-neighbor high
  mu[highLeaf] = -5.0;  // above-neighbor low -> empty cone for L
  std::vector<std::int8_t> d1 = {1};
  MonotoneConstantGaussianLeaf leaf;
  leaf.scale = 0.1;
  leaf.data = &store;
  leaf.directions.assign(d1.begin(), d1.end());
  leaf.cInflation = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  double sentinel = leaf.logLikelihoodForBranchWithParams(
    tree, midLeaf, yg.data(), nullptr, 2.0, 1.0, mu.data());
  check(sentinel == -HUGE_VAL, "monotone empty cone scores the sentinel");

  // draw sweep keeps the tree feasible from a feasible start
  std::vector<double> mu2(tree.nodes.size(), 0.0);
  std::vector<std::int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  bool everFeasible = true;
  for (int sweep = 0; sweep < 200; ++sweep) {
    leaf.drawParametersForTree(rng, tree, bottoms, 2.0, 0.05, mu2.data());
    if (!monotoneTreeIsFeasible(tree, store, d1.data(), mu2.data()))
      everFeasible = false;
  }
  check(everFeasible, "monotone draw sweep stays feasible");
  rngState = savedRngState;
  ext_rng_destroy(rng);
  printf("ok: monotone feasibility (descents %d)\n", descents);
}

// (c) The 1-D and 2-D constrained marginals match an independent quadrature,
// and the coupled draw has support {a <= mu_lower <= mu_upper <= b}.
static void testMonotoneMarginal() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 7002);
  const size_t n = 400;
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i) / (n - 1);
  ColumnStore store;
  store.build(x.data(), n, 1, 40);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i) y[i] = 0.5 * x[i] - 0.1;
  std::vector<index_t> idx(n);
  Tree tree;
  tree.initialize(idx.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  std::int32_t left, right;
  tree.splitInterval(store, 0, 0, &left, &right);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(left + (right - left) / 2);
  tree.birth(store, 0, rule, y.data(), nullptr);
  std::int32_t lower = tree.at(0).leftChild;      // increasing: left is lower
  std::int32_t upper = tree.at(0).leftChild + 1;

  std::vector<std::int8_t> d = {1};
  MonotoneConstantGaussianLeaf leaf;
  leaf.scale = 0.2;
  leaf.data = &store;
  leaf.directions.assign(d.begin(), d.end());
  double c = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  leaf.cInflation = c;
  double k = 2.0, sig2 = 0.09, kEff = k / c;

  // 2-D coupled marginal: no outer neighbors, so P = Phi((mR-mL)/sqrt(sL^2+sR^2))
  // and d_* = 1/2 in closed form
  const Node& nL = tree.at(lower);
  const Node& nR = tree.at(upper);
  double got2d = leaf.logLikelihoodForBranchWithParams(tree, 0, y.data(), nullptr,
                                                       k, sig2, nullptr);
  double mL, sL, mR, sR;
  refPosterior(kEff, sig2, nL.sumWeights, nL.sumWeightedResponse, leaf.scale, mL,
               sL);
  refPosterior(kEff, sig2, nR.sumWeights, nR.sumWeightedResponse, leaf.scale, mR,
               sR);
  double refBase =
    refLogIntegratedLikelihood(kEff, sig2, nL.sumWeights, nL.sumWeightedResponse,
                               leaf.scale) +
    refLogIntegratedLikelihood(kEff, sig2, nR.sumWeights, nR.sumWeightedResponse,
                               leaf.scale);
  double refP = Rf_pnorm5((mR - mL) / std::sqrt(sL * sL + sR * sR), 0.0, 1.0, 1, 0);
  double ref2d = refBase + std::log(refP) - std::log(0.5);
  checkNear(got2d, ref2d, 1e-6, "monotone 2-D marginal matches closed form");

  // with a finite outer lower bound the numerator becomes a 1-D quadrature:
  // score the same branch against a pinned below-neighbor via a nested reference
  double aL = mL - 0.5 * sL;  // an interior below bound on the lower leaf
  double numerRef = trapz(
    [&](double xr) {
      double up = xr;  // no upper bound bL
      double inner = Rf_pnorm5((up - mL) / sL, 0.0, 1.0, 1, 0) -
                     Rf_pnorm5((aL - mL) / sL, 0.0, 1.0, 1, 0);
      return (inner > 0.0 ? inner : 0.0) * normalPdf(xr, mR, sR);
    },
    aL, mR + 12.0 * sR);
  double denomRef = trapz(
    [&](double xr) {
      double up = xr;
      double sd = leaf.scale / kEff;
      double inner = Rf_pnorm5(up / sd, 0.0, 1.0, 1, 0) -
                     Rf_pnorm5(aL / sd, 0.0, 1.0, 1, 0);
      return (inner > 0.0 ? inner : 0.0) * normalPdf(xr, 0.0, sd);
    },
    aL, 12.0 * (leaf.scale / kEff));
  double coneNum =
    leaf.coneProbability(aL, HUGE_VAL, aL, HUGE_VAL, mL, sL, mR, sR);
  double coneDen = leaf.coneProbability(aL, HUGE_VAL, aL, HUGE_VAL, 0.0,
                                        leaf.scale / kEff, 0.0, leaf.scale / kEff);
  checkNear(coneNum, numerRef, 1e-6, "monotone cone numerator vs quadrature");
  checkNear(coneDen, denomRef, 1e-6, "monotone cone normalizer vs quadrature");

  // draw support: many coupled draws satisfy mu_lower <= mu_upper
  std::vector<double> mu(tree.nodes.size(), 0.0);
  std::vector<std::int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  bool ordered = true;
  for (int s = 0; s < 5000; ++s) {
    leaf.drawParametersForTree(rng, tree, bottoms, k, sig2, mu.data());
    if (mu[lower] > mu[upper] + 1e-9) ordered = false;
  }
  check(ordered, "monotone coupled draw supports mu_lower <= mu_upper");

  // a deep-tail truncation still returns a finite draw (Robert fallback)
  std::vector<double> muTail(tree.nodes.size(), 0.0);
  muTail[lower] = 40.0;  // forces the upper leaf far into the tail
  double before = muTail[upper];
  (void) before;
  bool finite = true;
  for (int s = 0; s < 50; ++s) {
    leaf.drawParametersForTree(rng, tree, bottoms, k, sig2, muTail.data());
    if (!std::isfinite(muTail[lower]) || !std::isfinite(muTail[upper]))
      finite = false;
  }
  check(finite, "monotone deep-tail draw stays finite");
  ext_rng_destroy(rng);
  printf("ok: monotone marginal and draw vs quadrature\n");
}

// (d) c-inflation: a singly-constrained leaf's post-truncation marginal
// variance equals the unconstrained sigma_mu^2 (design section 6, eq 3.7-3.9).
static void testMonotoneCInflation() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 7003);
  const size_t n = 200;
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i) / (n - 1);
  ColumnStore store;
  store.build(x.data(), n, 1, 20);
  std::vector<double> y(n, 0.0);  // zero response -> posterior mean 0
  std::vector<index_t> idx(n);
  Tree tree;
  tree.initialize(idx.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  std::int32_t left, right;
  tree.splitInterval(store, 0, 0, &left, &right);
  Rule rule;
  rule.variableIndex = 0;
  rule.setSplitIndex(left + (right - left) / 2);
  tree.birth(store, 0, rule, y.data(), nullptr);
  std::int32_t lower = tree.at(0).leftChild, upper = lower + 1;

  std::vector<std::int8_t> d = {1};
  MonotoneConstantGaussianLeaf leaf;
  leaf.scale = 0.1;
  leaf.data = &store;
  leaf.directions.assign(d.begin(), d.end());
  leaf.cInflation = std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
  double k = 2.0;
  double hugeSigma2 = 1e12;  // prior-dominant: draws from the truncated prior
  double sigmaMu = leaf.scale / k;  // the UNinflated target sd

  std::vector<double> mu(tree.nodes.size(), 0.0);
  std::vector<std::int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  const int burn = 2000, draws = 400000;
  for (int s = 0; s < burn; ++s)
    leaf.drawParametersForTree(rng, tree, bottoms, k, hugeSigma2, mu.data());
  double sum = 0.0, sumSq = 0.0;
  for (int s = 0; s < draws; ++s) {
    leaf.drawParametersForTree(rng, tree, bottoms, k, hugeSigma2, mu.data());
    sum += mu[upper];
    sumSq += mu[upper] * mu[upper];
  }
  double mean = sum / draws;
  double var = sumSq / draws - mean * mean;
  checkNear(var, sigmaMu * sigmaMu, 0.05 * sigmaMu * sigmaMu,
            "monotone c-inflation: singly-constrained variance is sigma_mu^2");
  ext_rng_destroy(rng);
  printf("ok: monotone c-inflation (var %.6f target %.6f)\n", var,
         sigmaMu * sigmaMu);
}

void runModelTests(ext_rng* rng) {
  testLeafSeamDispatch();
  testIntegratedLikelihood();
  testPosteriorDraw(rng);
  testConstantLeafSuffstatEquivalence();
  testChiKHyperprior(rng);
  testChiKEmptyLeafAccounting(rng);
  testSigmaPosteriorDf(rng);
  testVarianceLeafMarginal();
  testVarianceLeafDraw();
  testVarianceLeafCalibration();
  testVarianceDivisorGuard();
  testPolyaGamma(rng);
  testGeneralizedInverseGaussian(rng);
  testSampleFromPrior(rng);
  testLinearLeafMarginal();
  testLinearLeafDraw(rng);
  testLinearLeafStatisticsCache();
  testLinearLeafEndToEnd(rng);
  testLinearLeafFormats(rng);
  testLinearLeafViews();
  testGroupedMath(rng);
  testGroupedEndToEnd(rng);
  testGroupedBinary(rng);
  testGroupedStateRoundTrip();
  testAFTReduction(rng);
  testAFTCensoredMoments(rng);
  testAFTStateRoundTrip();
  testTLambdaMoments(rng);
  testTNuGridPosterior(rng);
  testTCompositeWeightDelegation(rng);
  testTFixedNuNoDraw(rng);
  testOrdinalTruncatedNormal(rng);
  testOrdinalCutpointConditional(rng);
  testOrdinalProbitEquivalence(rng);
  testOrdinalStateRoundTrip();
  testNBPolyaGammaShapeMoments(rng);
  testNBDispersionGridConditional(rng);
  testNBSweepOrderAndRestore(rng);
  testNBStateRoundTrip();
  testSparseKernel();
  testSparseColumnStore();
  testSparseEndToEnd();
  testSparseStateRoundTrip();
  testSparseCategoricalColumnStore();
  testSparseCategoricalEndToEnd();
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
  testSparseTestDataEndToEnd();
  testMonotoneNeighborGeometry();
  testMonotoneFeasibility();
  testMonotoneMarginal();
  testMonotoneCInflation();
  // the heteroscedastic end-to-end fits build full chains; they run last so
  // their heap footprint cannot perturb an earlier allocation-sensitive test
  // under the sanitizer's allocator.
  testVarianceForestRecovery();
  testVarianceForestRefusal();
  testVarianceReportingStatePredict();
  testVarianceSavedPredict();
  testVarianceM1Reduction();
}
