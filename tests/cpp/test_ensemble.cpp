// Ensemble-scale sum invariance for the backfit loop. Nearly every exact math
// gate in this suite runs one tree or a handful; the shipped default is 200,
// and the per-tree bookkeeping - the residual rolled from tree to tree, the
// aggregate rebuilt from it once per sweep - is only wrung out at that scale.
// A defect that needs m-tree accumulation to surface (a residual that drifts
// across the loop, a tree's fits retired twice, the last tree's contribution
// lost) passes every small-m gate. Two identities hold here, over EVERY
// observation, after EVERY sweep:
//
//   totalFits[i] == sum_t fits_t[i]                    (finalizeTotalFits)
//   treeY[i]     == y[i] - sum_{t < last} fits_t[i]    (the rolled residual)
//
// They are independent gates on different code even though they are one
// identity algebraically: totalFits is rebuilt as y - treeY + fits_last, so a
// finalize that reads the wrong tree breaks the first alone, while a drifted
// roll breaks both. Both right-hand sides are summed fresh in tree order,
// where the engine reaches the same quantity through an incremental roll, so
// the comparison is to within accumulation error and NOT bitwise: an == would
// fail on legitimate reassociation. The band below sits ~4 orders above the
// worst deviation measured over 200 trees x 30 sweeps (1.2e-15, about sqrt(m)
// ulps of a unit-scale total) and ~9 orders below one leaf value at this
// ensemble size (prior sd nodeScale / sqrt(m) = 0.035), so it discriminates
// against a lost or doubled contribution with room to spare either way.

#include "common.hpp"

namespace {

// n % 4 == 3, so the unrolled constant-leaf gathers run their prologue as well
// as their body; five columns keep 200 trees from all splitting on one.
constexpr size_t ensembleN = 503, ensembleP = 5, ensembleTrees = 200,
                 ensembleSweeps = 30, ensembleBurnIn = 20;
constexpr double ensembleTolerance = 1.0e-11;

/// Both identities for one settled sweep, every observation, no subsampling.
void checkEnsembleSweep(ConstantLeafSampler& sampler, size_t sweep,
                        double& worstFit, double& worstResidual) {
  std::vector<double> total = forestTotals(sampler, 0);
  std::vector<double> fits = TestPeer::treeFits(sampler.chain(0));
  const double* y = TestPeer::workingResponse(sampler.chain(0));
  const std::vector<double>& resid = TestPeer::residual(sampler.chain(0));

  double sweepFit = 0.0, sweepResidual = 0.0;
  for (size_t i = 0; i < ensembleN; ++i) {
    double allButLast = 0.0;
    for (size_t t = 0; t + 1 < ensembleTrees; ++t)
      allButLast += fits[t * ensembleN + i];
    double summed = allButLast + fits[(ensembleTrees - 1) * ensembleN + i];
    sweepFit = std::max(sweepFit, std::fabs(total[i] - summed));
    sweepResidual =
      std::max(sweepResidual, std::fabs(resid[i] - (y[i] - allButLast)));
  }

  char what[128];
  std::snprintf(what, sizeof what,
                "sweep %zu: totalFits sums all %zu trees (worst %.3g)", sweep,
                ensembleTrees, sweepFit);
  check(sweepFit < ensembleTolerance, what);
  std::snprintf(what, sizeof what,
                "sweep %zu: rolled residual retires %zu trees (worst %.3g)",
                sweep, ensembleTrees - 1, sweepResidual);
  check(sweepResidual < ensembleTolerance, what);
  worstFit = std::max(worstFit, sweepFit);
  worstResidual = std::max(worstResidual, sweepResidual);
}

// The level-fibre arm's own shape: small enough that the shift's 10 x 10
// covariance is estimated to a few parts in a thousand over the draw budget,
// large enough that the trees carry several leaves apiece.
constexpr size_t levelN = 300, levelP = 3, levelTrees = 10,
                 levelBurnIn = 60, levelDraws = 200000;
// 10 means plus 55 distinct covariance entries: a two-sided 4.5 leaves the
// whole family under 5e-4 by Bonferroni, and at 2e5 draws the gaussian
// moment approximations below are good to well inside it.
constexpr double levelZThreshold = 4.5;
// the projection is one fused subtraction per tree, so the residual sum is a
// few ulps of a leaf value (~1e-2 here), orders under this
constexpr double levelSumTolerance = 1.0e-12;

/// The closed form of section 1 at the homogeneous constant leaf, read off
/// the frozen forest: v_t = tau^2 / L_t and m_t = -S_t / L_t over tree t's
/// occupied leaves, then the zero-sum conditioning.
struct LevelLaw {
  std::vector<double> mean, variance, covariance;  // covariance is m x m
};

LevelLaw levelLawOf(ConstantLeafSampler& sampler, size_t numTrees) {
  ForestCalibration calibration = sampler.chain(0).forestCalibration(0);
  // priorScale is the response-unit total; the internal per-leaf sd is that
  // divided by the response transform and sqrt(m)
  double tau = calibration.priorSd /
               (calibration.responseScale *
                std::sqrt(static_cast<double>(numTrees)));

  LevelLaw law;
  law.mean.assign(numTrees, 0.0);
  law.variance.assign(numTrees, 0.0);
  std::vector<int32_t> bottoms;
  double sumVariance = 0.0, sumMean = 0.0;
  for (size_t t = 0; t < numTrees; ++t) {
    const Tree& tree = sampler.chain(0).treeInForest(0, t);
    const std::vector<double>& mu = TestPeer::muByTree(sampler.chain(0), t);
    bottoms.clear();
    tree.fillBottom(0, bottoms);
    double leafCount = 0.0, leafSum = 0.0;
    for (int32_t node : bottoms) {
      if (tree.at(node).numObservations() == 0) continue;
      leafCount += 1.0;
      leafSum += mu[static_cast<size_t>(node)];
    }
    law.variance[t] = tau * tau / leafCount;
    law.mean[t] = -leafSum / leafCount;
    sumVariance += law.variance[t];
    sumMean += law.mean[t];
  }
  double ratio = sumMean / sumVariance;
  law.covariance.assign(numTrees * numTrees, 0.0);
  for (size_t t = 0; t < numTrees; ++t) {
    law.mean[t] -= law.variance[t] * ratio;
    for (size_t u = 0; u < numTrees; ++u)
      law.covariance[t * numTrees + u] =
        (t == u ? law.variance[t] : 0.0) -
        law.variance[t] * law.variance[u] / sumVariance;
  }
  return law;
}

/// Empirical first and second moments of a stream of shift vectors, scored
/// against a law as per-coordinate z. The mean's standard error is
/// sqrt(C_tt / N); a gaussian sample covariance entry's is
/// sqrt((C_tt C_uu + C_tu^2) / N).
struct LevelMoments {
  size_t numDraws = 0;
  std::vector<double> sum, sumProducts;

  explicit LevelMoments(size_t numTrees)
    : sum(numTrees, 0.0), sumProducts(numTrees * numTrees, 0.0) {}

  void add(const std::vector<double>& shift) {
    size_t m = sum.size();
    ++numDraws;
    for (size_t t = 0; t < m; ++t) {
      sum[t] += shift[t];
      for (size_t u = t; u < m; ++u)
        sumProducts[t * m + u] += shift[t] * shift[u];
    }
  }

  void score(const LevelLaw& law, double& worstMeanZ, double& worstCovZ) const {
    size_t m = sum.size();
    double n = static_cast<double>(numDraws);
    worstMeanZ = 0.0;
    worstCovZ = 0.0;
    std::vector<double> empiricalMean(m);
    for (size_t t = 0; t < m; ++t) empiricalMean[t] = sum[t] / n;
    for (size_t t = 0; t < m; ++t) {
      double se = std::sqrt(law.covariance[t * m + t] / n);
      worstMeanZ = std::max(worstMeanZ,
                            std::fabs(empiricalMean[t] - law.mean[t]) / se);
      for (size_t u = t; u < m; ++u) {
        double empirical =
          sumProducts[t * m + u] / n - empiricalMean[t] * empiricalMean[u];
        double target = law.covariance[t * m + u];
        double se2 = std::sqrt((law.covariance[t * m + t] *
                                  law.covariance[u * m + u] +
                                target * target) / n);
        worstCovZ = std::max(worstCovZ, std::fabs(empirical - target) / se2);
      }
    }
  }
};

double standardNormal() {
  double u1 = runif01(), u2 = runif01();
  return std::sqrt(-2.0 * std::log(u1 + 1e-300)) *
         std::cos(6.283185307179586 * u2);
}

// ---- the probit rescaling step ---------------------------------------------
//
// A slice step from a fixed point is not a draw from its target, so both
// distributional tests here are in invariance form: start each replication
// from a point drawn exactly from the target, take one step, and score the
// result against the target by Kolmogorov-Smirnov. At 2e5 replications the
// p = 1e-3 bound is 1.95 / sqrt(N).
constexpr size_t rescaleN = 243, rescaleP = 3, rescaleTrees = 8,
                 rescaleBurnIn = 150, rescaleDraws = 200000;
constexpr double rescaleKsBound = 1.95;
constexpr double rescaleMapTolerance = 1.0e-14;

/// A one-dimensional log density tabulated on a uniform grid, its CDF by the
/// trapezoid rule, sampled by inverse CDF and read back by interpolation.
struct GridTarget {
  double lower = 0.0, step = 0.0;
  std::vector<double> cdf;

  template <typename F>
  GridTarget(const F& logDensity, double from, double to, size_t numPoints) {
    // locate the mass on a coarse pass, then tabulate +-30 sd of it finely
    double best = -HUGE_VAL, mode = from;
    for (size_t j = 0; j <= 20000; ++j) {
      double v = from + (to - from) * static_cast<double>(j) / 20000.0;
      double value = logDensity(v);
      if (value > best) {
        best = value;
        mode = v;
      }
    }
    double curvature = -(logDensity(mode + 1e-4) - 2.0 * logDensity(mode) +
                         logDensity(mode - 1e-4)) / 1e-8;
    double sd = curvature > 0.0 ? 1.0 / std::sqrt(curvature) : 1.0;
    lower = mode - 30.0 * sd;
    step = 60.0 * sd / static_cast<double>(numPoints - 1);
    std::vector<double> density(numPoints);
    for (size_t j = 0; j < numPoints; ++j)
      density[j] = std::exp(logDensity(lower + step * static_cast<double>(j)) -
                            best);
    cdf.assign(numPoints, 0.0);
    for (size_t j = 1; j < numPoints; ++j)
      cdf[j] = cdf[j - 1] + 0.5 * step * (density[j - 1] + density[j]);
    double total = cdf.back();
    for (double& c : cdf) c /= total;
  }

  double draw(double u) const {
    size_t j = static_cast<size_t>(
      std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin());
    if (j == 0) return lower;
    if (j >= cdf.size()) return lower + step * static_cast<double>(cdf.size() - 1);
    double width = cdf[j] - cdf[j - 1];
    double fraction = width > 0.0 ? (u - cdf[j - 1]) / width : 0.0;
    return lower + step * (static_cast<double>(j - 1) + fraction);
  }

  double at(double v) const {
    double position = (v - lower) / step;
    if (position <= 0.0) return 0.0;
    size_t j = static_cast<size_t>(position);
    if (j + 1 >= cdf.size()) return 1.0;
    double fraction = position - static_cast<double>(j);
    return cdf[j] + fraction * (cdf[j + 1] - cdf[j]);
  }
};

/// sqrt(N) times the KS distance of a sample from a continuous CDF.
template <typename F>
double scaledKsDistance(std::vector<double> sample, const F& cdf) {
  std::sort(sample.begin(), sample.end());
  double n = static_cast<double>(sample.size()), worst = 0.0;
  for (size_t j = 0; j < sample.size(); ++j) {
    double value = cdf(sample[j]);
    worst = std::max(worst, std::max(static_cast<double>(j + 1) / n - value,
                                     value - static_cast<double>(j) / n));
  }
  return worst * std::sqrt(n);
}

/// The step's four constants read off a chain's raw state: the latents, the
/// offset, the test's own mask and the leaves gathered in tree order, never
/// the step's own arithmetic.
struct RescaleConstants {
  double numActive = 0.0, residualSquares = 0.0, offsetCross = 0.0,
         kTerm = 0.0, degreesOfFreedom = 0.0;

  double logDensity(double v) const {
    double e = std::exp(v);
    return (numActive - degreesOfFreedom) * v - 0.5 * residualSquares * e * e +
           offsetCross * e - kTerm / (e * e);
  }
};

RescaleConstants rescaleConstantsOf(ConstantLeafSampler& sampler,
                                    const std::vector<double>& mask) {
  Chain<ConstantGaussianLeaf>& chain = sampler.chain(0);
  const double* z = TestPeer::latents(chain);
  const double* offset = TestPeer::offset(chain);
  RescaleConstants constants;
  for (size_t i = 0; i < rescaleN; ++i) {
    if (mask[i] == 0.0) continue;
    double fit = 0.0;
    for (size_t t = 0; t < rescaleTrees; ++t)
      fit += TestPeer::muByTree(chain, t)[TestPeer::leafOf(chain, t)[i]];
    double residual = z[i] - fit;
    constants.numActive += 1.0;
    constants.residualSquares += residual * residual;
    constants.offsetCross += offset[i] * residual;
  }
  const ChiKHyperprior& prior = TestPeer::kHyperprior(chain);
  double k = TestPeer::forestK(chain);
  constants.kTerm = k * k / (2.0 * prior.scale * prior.scale);
  constants.degreesOfFreedom = prior.degreesOfFreedom;
  return constants;
}

/// Whether a call declines and leaves the generator where it was.
bool rescaleDeclinesWithoutDraw(Chain<ConstantGaussianLeaf>& chain) {
  ext_rng* rng = TestPeer::rng(chain);
  ext_rng* reference = cloneRng(rng);
  double alpha = 0.0;
  bool taken = TestPeer::drawForestRescale(chain, &alpha);
  bool agree = rngStreamsAgree(reference, rng);
  ext_rng_destroy(reference);
  return !taken && alpha == 1.0 && agree;
}

std::unique_ptr<ConstantLeafSampler>
makeRescaleSampler(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& offset, ResponseFamily family,
                   bool updateK, bool rescale, ext_rng** rngs,
                   size_t numChains = 1, size_t numVarianceTrees = 0) {
  SamplerOptions options;
  options.numTrees = rescaleTrees;
  options.numChains = numChains;
  options.numVarianceTrees = numVarianceTrees;
  options.nodeScale = 3.0;
  options.updateK = updateK;
  options.kHyperprior.degreesOfFreedom = 1.5;
  options.kHyperprior.scale = 2.0;
  options.probitRescaleForest = rescale;
  return std::make_unique<ConstantLeafSampler>(
    x.data(), y.data(), rescaleN, rescaleP, nullptr, offset.data(), family,
    1.0, 3.0, 0.37804942330213542, options, rngs);
}

/// The slice step on a normal target with a closed-form CDF, in invariance
/// form: v0 drawn exactly, one step from it, KS on the result.
double sliceInvarianceKs(ext_rng* rng, double mean, double sd, double width,
                         int stepLimit) {
  std::vector<double> moved(rescaleDraws);
  for (size_t r = 0; r < rescaleDraws; ++r) {
    double start = mean + sd * ext_rng_simulateStandardNormal(rng);
    auto logDensity = [&](double u) {
      double d = (start + u - mean) / sd;
      return -0.5 * d * d;
    };
    moved[r] = start + sliceFromZero(rng, logDensity, width, stepLimit);
  }
  return scaledKsDistance(moved, [&](double v) {
    return 0.5 * std::erfc(-(v - mean) / (sd * std::sqrt(2.0)));
  });
}

void runForestRescaleTests() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20261008u);

  std::vector<double> x(rescaleN * rescaleP), y(rescaleN), offset(rescaleN);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < rescaleN; ++i) {
    double eta = 1.5 * std::sin(3.0 * x[i]) + x[i + rescaleN] - 0.5;
    y[i] = eta + standardNormal() > 0.0 ? 1.0 : 0.0;
    offset[i] = 0.8 * (x[i + 2 * rescaleN] - 0.5);
  }
  // every fifth row inactive, so the leaves hold masked members
  std::vector<double> mask(rescaleN, 1.0);
  for (size_t i = 0; i < rescaleN; i += 5) mask[i] = 0.0;

  std::unique_ptr<ConstantLeafSampler> sampler = makeRescaleSampler(
    x, y, offset, ResponseFamily::probit, true, true, &rng);
  Results results;
  sampler->run(rescaleBurnIn, 0, results);
  check(sampler->setActiveRows(mask.data()), "the probit sampler takes a mask");
  sampler->run(20, 0, results);
  Chain<ConstantGaussianLeaf>& chain = sampler->chain(0);
  // no sweep leaves an empty leaf, so one is stranded; its value is nonzero so
  // that scaling it would show
  check(TestPeer::strandEmptyLeaf(chain, 0.5).second >= 0,
        "rescaling: the fixture strands an empty leaf");
  check(TestPeer::forestRescaleApplies(chain),
        "a burned-in drawn-k probit forest takes the rescaling step");

  // the frozen state
  std::vector<double> frozenZ(TestPeer::latents(chain),
                              TestPeer::latents(chain) + rescaleN);
  std::vector<std::vector<double>> frozenMu(rescaleTrees);
  for (size_t t = 0; t < rescaleTrees; ++t)
    frozenMu[t] = TestPeer::muByTree(chain, t);
  const double frozenK = TestPeer::forestK(chain);
  auto moveFrozenStateTo = [&](double factor) {
    std::vector<double> z(frozenZ);
    for (size_t i = 0; i < rescaleN; ++i)
      if (mask[i] != 0.0) z[i] *= factor;
    TestPeer::restoreLatents(chain, z.data());
    for (size_t t = 0; t < rescaleTrees; ++t) {
      std::vector<double>& mu = TestPeer::muByTree(chain, t);
      for (size_t j = 0; j < mu.size(); ++j) mu[j] = frozenMu[t][j] * factor;
    }
    TestPeer::forestK(chain) = frozenK / factor;
    TestPeer::rebuildTotalFits(chain);
  };
  // the occupied leaves, as (tree, node) pairs
  std::vector<std::pair<size_t, int32_t>> occupied;
  std::vector<int32_t> bottoms;
  size_t numEmpty = 0;
  for (size_t t = 0; t < rescaleTrees; ++t) {
    const Tree& tree = chain.treeInForest(0, t);
    bottoms.clear();
    tree.fillBottom(0, bottoms);
    for (int32_t b : bottoms) {
      if (tree.at(b).numObservations() > 0) occupied.emplace_back(t, b);
      else ++numEmpty;
    }
  }
  check(occupied.size() > rescaleTrees && numEmpty > 0,
        "the frozen trees carry splits and an empty leaf");
  moveFrozenStateTo(1.0);
  RescaleConstants constants = rescaleConstantsOf(*sampler, mask);
  char what[200];

  // 1. The conditional, in invariance form. v0 is the orbit coordinate
  // measured from the frozen state, so its target is the frozen state's own
  // conditional; one step from e^v0 then lands at v1, read off the latents and
  // the leaves, which must agree with each other, with alpha and with k.
  GridTarget target([&](double v) { return constants.logDensity(v); }, -6.0, 6.0,
              200001);
  std::vector<double> landed(rescaleDraws);
  double worstDisagreement = 0.0;
  for (size_t r = 0; r < rescaleDraws; ++r) {
    double start = target.draw(ext_rng_simulateContinuousUniform(rng));
    moveFrozenStateTo(std::exp(start));
    double alpha = 0.0;
    TestPeer::drawForestRescale(chain, &alpha);
    const double* z = TestPeer::latents(chain);
    double zMoved = 0.0, zFrozen = 0.0, muMoved = 0.0, muFrozen = 0.0;
    for (size_t i = 0; i < rescaleN; ++i) {
      if (mask[i] == 0.0) continue;
      zMoved += std::fabs(z[i]);
      zFrozen += std::fabs(frozenZ[i]);
    }
    for (const auto& [t, b] : occupied) {
      muMoved += std::fabs(TestPeer::muByTree(chain, t)[static_cast<size_t>(b)]);
      muFrozen += std::fabs(frozenMu[t][static_cast<size_t>(b)]);
    }
    double fromLatents = std::log(zMoved / zFrozen);
    double fromLeaves = std::log(muMoved / muFrozen);
    double fromK = std::log(frozenK / TestPeer::forestK(chain));
    double fromAlpha = start + std::log(alpha);
    worstDisagreement = std::max(
      {worstDisagreement, std::fabs(fromLeaves - fromLatents),
       std::fabs(fromK - fromLatents), std::fabs(fromAlpha - fromLatents)});
    landed[r] = fromLatents;
  }
  double ks = scaledKsDistance(landed, [&](double v) { return target.at(v); });
  std::snprintf(what, sizeof what,
                "rescaling: one step from the target stays at the target "
                "(sqrt(N) D %.3g against %.3g)", ks, rescaleKsBound);
  check(ks < rescaleKsBound, what);
  std::snprintf(what, sizeof what,
                "rescaling: latents, leaves, k and alpha move by one factor "
                "(worst %.3g)", worstDisagreement);
  check(worstDisagreement < 1.0e-10, what);

  // 2. The mapping, one step from the frozen state.
  moveFrozenStateTo(1.0);
  std::vector<double> totalsBefore = forestTotals(*sampler, 0);
  double alpha = 0.0;
  check(TestPeer::drawForestRescale(chain, &alpha) && alpha != 1.0,
        "rescaling: the frozen state takes the step");
  bool kMapped = std::fabs(TestPeer::forestK(chain) * alpha / frozenK - 1.0) <
                 rescaleMapTolerance;
  bool leavesMapped = true;
  for (size_t t = 0; t < rescaleTrees; ++t) {
    const Tree& tree = chain.treeInForest(0, t);
    bottoms.clear();
    tree.fillBottom(0, bottoms);
    for (int32_t b : bottoms) {
      double now = TestPeer::muByTree(chain, t)[static_cast<size_t>(b)];
      double was = frozenMu[t][static_cast<size_t>(b)];
      leavesMapped = leavesMapped &&
        (tree.at(b).numObservations() > 0
           ? std::fabs(now - alpha * was) <= rescaleMapTolerance * std::fabs(now)
           : now == was);
    }
  }
  const double* z = TestPeer::latents(chain);
  const double* working = TestPeer::workingResponse(chain);
  bool latentsMapped = true, workingMapped = true;
  for (size_t i = 0; i < rescaleN; ++i) {
    latentsMapped = latentsMapped &&
      (mask[i] != 0.0
         ? std::fabs(z[i] - alpha * frozenZ[i]) <= rescaleMapTolerance * std::fabs(z[i])
         : z[i] == frozenZ[i]);
    workingMapped = workingMapped && working[i] == z[i] - offset[i];
  }
  // inside the bound totalFits is scaled in place and the factor accumulates;
  // the scaled cache stays within rounding of the leaves it caches
  auto gatherTotals = [&]() {
    std::vector<double> gathered(rescaleN, 0.0);
    for (size_t t = 0; t < rescaleTrees; ++t)
      for (size_t i = 0; i < rescaleN; ++i)
        gathered[i] +=
          TestPeer::muByTree(chain, t)[TestPeer::leafOf(chain, t)[i]];
    return gathered;
  };
  std::vector<double> totals = forestTotals(*sampler, 0), gathered = gatherTotals();
  bool totalsMapped = TestPeer::totalFitsScale(chain) == alpha;
  for (size_t i = 0; i < rescaleN; ++i)
    totalsMapped = totalsMapped && totals[i] == totalsBefore[i] * alpha &&
      std::fabs(totals[i] - gathered[i]) <=
        rescaleMapTolerance * (1.0 + std::fabs(gathered[i]));
  // past the bound the step re-derives the cache from the scaled leaves in
  // tree order, bitwise, and resets the factor
  moveFrozenStateTo(1.0);
  TestPeer::totalFitsScale(chain) = 1.0e6;
  TestPeer::drawForestRescale(chain, &alpha);
  bool totalsRederived = forestTotals(*sampler, 0) == gatherTotals() &&
                         TestPeer::totalFitsScale(chain) == 1.0;
  check(kMapped, "rescaling: k is divided by alpha");
  check(leavesMapped, "rescaling: every occupied leaf of every tree is "
                      "multiplied by alpha, an empty one left as it was");
  check(latentsMapped, "rescaling: the active latents are multiplied by "
                       "alpha, the inactive ones left as they were");
  check(workingMapped, "rescaling: the working response is alpha z - o");
  check(totalsMapped, "rescaling: inside the bound totalFits is scaled in "
                      "place and the factor accumulates");
  check(totalsRederived, "rescaling: past the bound totalFits is re-derived "
                         "from the leaves in tree order, bitwise, and the "
                         "factor reset");

  // 3. Inertness: each decline consumes no generator draw.
  moveFrozenStateTo(1.0);
  TestPeer::forestK(chain) = HUGE_VAL;
  check(rescaleDeclinesWithoutDraw(chain), "rescaling declines at infinite k");
  moveFrozenStateTo(1.0);
  ChiKHyperprior savedPrior = TestPeer::kHyperprior(chain);
  TestPeer::kHyperprior(chain).scale = HUGE_VAL;
  TestPeer::kHyperprior(chain).degreesOfFreedom = constants.numActive;
  check(rescaleDeclinesWithoutDraw(chain),
        "rescaling declines under an infinite k scale with n <= nu");
  TestPeer::kHyperprior(chain).degreesOfFreedom = constants.numActive - 1.0;
  check(TestPeer::drawForestRescale(chain, &alpha),
        "rescaling is taken under an infinite k scale with n > nu");
  TestPeer::kHyperprior(chain) = savedPrior;
  moveFrozenStateTo(1.0);
  std::vector<double> none(rescaleN, 0.0);
  check(sampler->setActiveRows(none.data()), "an all-inactive mask installs");
  check(rescaleDeclinesWithoutDraw(chain),
        "rescaling declines with every row inactive");
  check(sampler->setActiveRows(mask.data()), "the mask reinstalls");

  ext_rng* otherRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(otherRng, 20261009u);
  std::unique_ptr<ConstantLeafSampler> fixedK = makeRescaleSampler(
    x, y, offset, ResponseFamily::probit, false, true, &otherRng);
  fixedK->run(20, 0, results);
  check(rescaleDeclinesWithoutDraw(fixedK->chain(0)),
        "rescaling declines at a fixed k");
  std::unique_ptr<ConstantLeafSampler> logistic = makeRescaleSampler(
    x, y, offset, ResponseFamily::logistic, true, true, &otherRng);
  logistic->run(20, 0, results);
  check(rescaleDeclinesWithoutDraw(logistic->chain(0)),
        "rescaling declines under the logistic family");
  std::unique_ptr<ConstantLeafSampler> stale = makeRescaleSampler(
    x, y, offset, ResponseFamily::probit, true, true, &otherRng);
  stale->run(20, 0, results);
  stale->sampleTreesFromPrior();
  check(TestPeer::leafOfStale(stale->chain(0), 0) != 0 &&
          rescaleDeclinesWithoutDraw(stale->chain(0)),
        "rescaling declines on a stale tree map");
  // every other family, and a probit fit with a second forest, all drawing k
  std::unique_ptr<ConstantLeafSampler> gaussian = makeRescaleSampler(
    x, y, offset, ResponseFamily::gaussian, true, true, &otherRng);
  gaussian->run(20, 0, results);
  check(rescaleDeclinesWithoutDraw(gaussian->chain(0)),
        "rescaling declines under the gaussian family");
  std::unique_ptr<ConstantLeafSampler> nbinom = makeRescaleSampler(
    x, y, offset, ResponseFamily::nbinom, true, true, &otherRng);
  nbinom->run(20, 0, results);
  check(rescaleDeclinesWithoutDraw(nbinom->chain(0)),
        "rescaling declines under the nbinom family");
  std::unique_ptr<ConstantLeafSampler> variance = makeRescaleSampler(
    x, y, offset, ResponseFamily::gaussian, true, true, &otherRng, 1, 4);
  variance->run(20, 0, results);
  check(variance->chain(0).hasVarianceForest() &&
          rescaleDeclinesWithoutDraw(variance->chain(0)),
        "rescaling declines beside a variance forest");
  {
    std::vector<double> z(rescaleN);
    for (size_t i = 0; i < rescaleN; ++i) z[i] = i % 2 == 0 ? 1.0 : 0.0;
    SamplerOptions options;
    options.nodeScale = 3.0;
    options.updateK = true;
    options.kHyperprior.degreesOfFreedom = 1.5;
    options.kHyperprior.scale = 2.0;
    AmplitudeSpec spec;
    spec.family = ResponseFamily::probit;
    spec.mu.numTrees = rescaleTrees;
    spec.tau.numTrees = 4;
    spec.z = z.data();
    ConstantLeafSampler twoForests(x.data(), y.data(), rescaleN, rescaleP,
                                   nullptr, offset.data(), 1.0, 3.0, 1.0,
                                   options, spec, &otherRng);
    twoForests.run(20, 0, results);
    // a combining chain holds k fixed; marked drawn, only the forest count and
    // the combiner are left to decline on
    Chain<ConstantGaussianLeaf>& twoForestChain = twoForests.chain(0);
    TestPeer::forestUpdateK(twoForestChain) = true;
    TestPeer::kHyperprior(twoForestChain) = options.kHyperprior;
    check(twoForestChain.numForests() == 2 &&
            TestPeer::forestK(twoForestChain) > 0.0 &&
            TestPeer::forestK(twoForestChain) <= DBL_MAX &&
            rescaleDeclinesWithoutDraw(twoForestChain),
          "rescaling declines on a two-forest probit fit");
  }
  ext_rng_destroy(otherRng);

  // the step acts on every chain: two chains, each against its own seed run
  // with the step off, all part
  std::vector<ext_rng*> onRngs(2), offRngs(2);
  for (size_t c = 0; c < 2; ++c) {
    onRngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    offRngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(onRngs[c], 11u + static_cast<std::uint32_t>(c));
    ext_rng_setSeed(offRngs[c], 11u + static_cast<std::uint32_t>(c));
  }
  std::unique_ptr<ConstantLeafSampler> on = makeRescaleSampler(
    x, y, offset, ResponseFamily::probit, true, true, onRngs.data(), 2);
  std::unique_ptr<ConstantLeafSampler> off = makeRescaleSampler(
    x, y, offset, ResponseFamily::probit, true, false, offRngs.data(), 2);
  on->run(30, 0, results);
  off->run(30, 0, results);
  for (size_t c = 0; c < 2; ++c) {
    std::snprintf(what, sizeof what, "rescaling acts on chain %zu", c);
    check(forestTotals(*on, c) != forestTotals(*off, c), what);
  }
  for (size_t c = 0; c < 2; ++c) {
    ext_rng_destroy(onRngs[c]);
    ext_rng_destroy(offRngs[c]);
  }

  // 4. The slice step alone. At the shipped width and limit the limit never
  // binds; at a width of a sd and a limit of 3 it binds on most steps, which
  // is the only run that can tell a randomly split limit from a fixed one.
  double ksShipped = sliceInvarianceKs(rng, 0.3, 0.05, rescaleSliceWidth,
                                       rescaleSliceStepLimit);
  double ksBinding = sliceInvarianceKs(rng, 0.3, 1.0, 1.0, 3);
  std::snprintf(what, sizeof what,
                "slice step: invariant at the shipped width and limit "
                "(sqrt(N) D %.3g) and where the limit binds (%.3g)",
                ksShipped, ksBinding);
  check(ksShipped < rescaleKsBound && ksBinding < rescaleKsBound, what);

  ext_rng_destroy(rng);
  printf("ok: probit rescaling, %zu steps at %zu trees (sqrt(N) D %.3g, worst "
         "factor disagreement %.3g, %zu empty leaves); slice sqrt(N) D %.3g "
         "shipped, %.3g limit binding\n",
         rescaleDraws, rescaleTrees, ks, worstDisagreement, numEmpty,
         ksShipped, ksBinding);
}

}  // namespace

void runEnsembleTests() {
  // own generator and own runif01 stream, restored on the way out: the suite
  // must read the same under a filter as it does in the full run, and must
  // not shift any other suite's hardcoded draws
  std::uint64_t savedRngState = rngState;
  rngState = 0x9e3779b97f4a7c15ull;

  std::vector<double> x(ensembleN * ensembleP), y(ensembleN);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < ensembleN; ++i)
    y[i] = std::sin(3.0 * x[i]) +
           2.0 * x[i + ensembleN] * x[i + 2 * ensembleN] -
           x[i + 3 * ensembleN] + 0.3 * (runif01() - 0.5);

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20260817u);

  SamplerOptions options;
  options.numTrees = ensembleTrees;  // the shipped default, stated not inherited
  options.levelGibbs = LevelGibbsMode::off;
  ConstantLeafSampler sampler(x.data(), y.data(), ensembleN, ensembleP, nullptr,
                              nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);

  double worstFit = 0.0, worstResidual = 0.0;
  Results results;
  for (size_t sweep = 0; sweep < ensembleSweeps; ++sweep) {
    // one sweep per call so the identities see every one, burn-in included:
    // the trees are stumps at sweep 0 and deep by the end, and a drift that
    // needs several sweeps to clear the band has room to show. The tail
    // records, putting the sweep body's recording branch under them too.
    if (sweep < ensembleBurnIn) sampler.run(1, 0, results);
    else sampler.run(0, 1, results);
    checkEnsembleSweep(sampler, sweep, worstFit, worstResidual);
  }

  // a degenerate ensemble (every fit zero, or NaN throughout) satisfies both
  // identities vacuously, so pin that the run actually fit something
  bool allFinite = true;
  double magnitude = 0.0;
  for (double v : forestTotals(sampler, 0)) {
    allFinite = allFinite && std::isfinite(v);
    magnitude = std::max(magnitude, std::fabs(v));
  }
  check(allFinite && magnitude > 0.1,
        "the 200-tree ensemble left a finite, non-trivial fit");

  // The same two identities with the level-fibre step on. What this arm
  // checks is narrow and worth saying: the assertion fires at the END of a
  // sweep, by which point every leaf has been redrawn, so all that reaches it
  // is the projection's own arithmetic residual sum_t c_t riding in through
  // totalFits - of order m eps tau, and the band above sits nine orders over
  // that. It sizes the missing-projection poison, which would put a whole
  // leaf value there, and says nothing about whether the draw's LAW is right;
  // the frozen-forest gate below is what scores that.
  ext_rng* levelRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(levelRng, 20260907u);
  SamplerOptions shiftedOptions;
  shiftedOptions.numTrees = ensembleTrees;
  shiftedOptions.levelGibbs = LevelGibbsMode::on;
  ConstantLeafSampler shifted(x.data(), y.data(), ensembleN, ensembleP,
                              nullptr, nullptr, ResponseFamily::gaussian, 1.0,
                              3.0, 0.37804942330213542, shiftedOptions,
                              &levelRng);
  double shiftedWorstFit = 0.0, shiftedWorstResidual = 0.0;
  for (size_t sweep = 0; sweep < ensembleSweeps; ++sweep) {
    if (sweep < ensembleBurnIn) shifted.run(1, 0, results);
    else shifted.run(0, 1, results);
    checkEnsembleSweep(shifted, sweep, shiftedWorstFit, shiftedWorstResidual);
  }
  ext_rng_destroy(levelRng);

  // The same arm entered from a prior tree draw, which is where a fresh
  // bart2 sampler starts. That entry leaves every obs-to-leaf map marked for
  // rebuild and every leaf at zero against a zeroed aggregate, so a shift
  // taken through the stale maps would put a constant in the residual that
  // totalFits then carries for the rest of the run - and, being carried, it
  // shows on every later sweep rather than only the first.
  ext_rng* priorRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(priorRng, 777u);
  ConstantLeafSampler fromPrior(x.data(), y.data(), ensembleN, ensembleP,
                                nullptr, nullptr, ResponseFamily::gaussian,
                                1.0, 3.0, 0.37804942330213542, shiftedOptions,
                                &priorRng);
  fromPrior.sampleTreesFromPrior();
  double priorWorstFit = 0.0, priorWorstResidual = 0.0;
  for (size_t sweep = 0; sweep < 5; ++sweep) {
    fromPrior.run(1, 0, results);
    checkEnsembleSweep(fromPrior, sweep, priorWorstFit, priorWorstResidual);
  }
  ext_rng_destroy(priorRng);

  // and the step actually moved the leaves it is supposed to: one more sweep
  // whose only difference from the arm above is the flag
  ext_rng* offRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(offRng, 20260907u);
  SamplerOptions offOptions;
  offOptions.numTrees = ensembleTrees;
  offOptions.levelGibbs = LevelGibbsMode::off;
  ConstantLeafSampler off(x.data(), y.data(), ensembleN, ensembleP, nullptr,
                          nullptr, ResponseFamily::gaussian, 1.0, 3.0,
                          0.37804942330213542, offOptions, &offRng);
  off.run(ensembleSweeps, 0, results);
  bool anyDifference = forestTotals(off, 0) != forestTotals(shifted, 0);
  check(anyDifference, "the level-fibre flag moves the sampled path");
  ext_rng_destroy(offRng);

  // ---- the shift's own law, at a frozen forest -------------------------
  //
  // The step is an exact Gibbs draw, so there is no acceptance ratio for a
  // detailed-balance script to score; what replaces one is a direct test of
  // the conditional. Freeze a burned-in forest, restore its leaf tables
  // between calls so every draw scores against the SAME closed form, and read
  // the empirical mean and covariance of c.
  ext_rng* lawRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(lawRng, 3141592u);

  std::vector<double> levelX(levelN * levelP), levelY(levelN);
  for (double& v : levelX) v = runif01();
  for (size_t i = 0; i < levelN; ++i)
    levelY[i] = 2.0 * levelX[i] - levelX[i + levelN] * levelX[i + 2 * levelN] +
                0.3 * (runif01() - 0.5);

  SamplerOptions lawOptions;
  lawOptions.numTrees = levelTrees;
  lawOptions.levelGibbs = LevelGibbsMode::off;
  ConstantLeafSampler lawSampler(levelX.data(), levelY.data(), levelN, levelP,
                                 nullptr, nullptr, ResponseFamily::gaussian,
                                 1.0, 3.0, 0.37804942330213542, lawOptions,
                                 &lawRng);
  lawSampler.run(levelBurnIn, 0, results);

  std::vector<std::vector<double>> frozenMu(levelTrees);
  for (size_t t = 0; t < levelTrees; ++t)
    frozenMu[t] = TestPeer::muByTree(lawSampler.chain(0), t);
  LevelLaw law = levelLawOf(lawSampler, levelTrees);

  LevelMoments moments(levelTrees);
  std::vector<double> shift(levelTrees, 0.0);
  double worstSum = 0.0;
  bool everyDrawEligible = true;
  for (size_t draw = 0; draw < levelDraws; ++draw) {
    for (size_t t = 0; t < levelTrees; ++t)
      TestPeer::muByTree(lawSampler.chain(0), t) = frozenMu[t];
    everyDrawEligible = everyDrawEligible &&
      TestPeer::drawLevelShift(lawSampler.chain(0), shift.data());
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t) total += shift[t];
    worstSum = std::max(worstSum, std::fabs(total));
    moments.add(shift);
  }
  for (size_t t = 0; t < levelTrees; ++t)
    TestPeer::muByTree(lawSampler.chain(0), t) = frozenMu[t];

  check(everyDrawEligible, "the frozen forest is eligible on every draw");
  char what[160];
  std::snprintf(what, sizeof what,
                "the shift sums to zero on every draw (worst %.3g)", worstSum);
  check(worstSum < levelSumTolerance, what);

  double worstMeanZ = 0.0, worstCovZ = 0.0;
  moments.score(law, worstMeanZ, worstCovZ);
  std::snprintf(what, sizeof what,
                "the shift's mean matches its conditional (worst z %.3g)",
                worstMeanZ);
  check(worstMeanZ < levelZThreshold, what);
  std::snprintf(what, sizeof what,
                "the shift's covariance matches its conditional (worst z %.3g)",
                worstCovZ);
  check(worstCovZ < levelZThreshold, what);

  // Poison (i), the projection dropped: u applied directly. sum_t u_t is then
  // N(sum m_t, sum v_t), so the zero-sum check - and with it the ensemble
  // identity above, which is where sum_t c_t reaches totalFits - fails by
  // orders rather than marginally.
  double poisonWorstSum = 0.0;
  for (size_t draw = 0; draw < 1000; ++draw) {
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t)
      total += law.mean[t] + std::sqrt(law.variance[t]) * standardNormal();
    poisonWorstSum = std::max(poisonWorstSum, std::fabs(total));
  }
  std::snprintf(what, sizeof what,
                "poison (i), no projection: the sum moves (%.3g against %.3g)",
                poisonWorstSum, worstSum);
  check(poisonWorstSum > 1.0e6 * levelSumTolerance, what);

  // Poison (ii), the prior-mean term dropped: u centred at zero and projected.
  // The law keeps its covariance and loses its mean, which is the poison that
  // separates a correct conditional from a correctly-shaped one - so the
  // covariance score must still PASS while the mean score fails.
  LevelMoments poisoned(levelTrees);
  double sumVariance = 0.0;
  for (size_t t = 0; t < levelTrees; ++t) sumVariance += law.variance[t];
  for (size_t draw = 0; draw < levelDraws; ++draw) {
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t) {
      shift[t] = std::sqrt(law.variance[t]) * standardNormal();
      total += shift[t];
    }
    double ratio = total / sumVariance;
    for (size_t t = 0; t < levelTrees; ++t)
      shift[t] -= law.variance[t] * ratio;
    poisoned.add(shift);
  }
  double poisonMeanZ = 0.0, poisonCovZ = 0.0;
  poisoned.score(law, poisonMeanZ, poisonCovZ);
  std::snprintf(what, sizeof what,
                "poison (ii), no prior mean: the mean fails (z %.3g)",
                poisonMeanZ);
  check(poisonMeanZ > levelZThreshold, what);
  std::snprintf(what, sizeof what,
                "poison (ii), no prior mean: the covariance still passes "
                "(z %.3g)", poisonCovZ);
  check(poisonCovZ < levelZThreshold, what);

  ext_rng_destroy(lawRng);

  runForestRescaleTests();

  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: ensemble sum invariance, %zu trees x %zu sweeps (worst fit %.3g, "
         "worst residual %.3g)\n", ensembleTrees, ensembleSweeps, worstFit,
         worstResidual);
  printf("ok: level fibre, flag on (worst fit %.3g, worst residual %.3g, "
         "from a prior draw %.3g); "
         "%zu draws at %zu trees (worst |sum| %.3g, mean z %.3g, cov z %.3g); "
         "poisons (no projection |sum| %.3g; no prior mean, z %.3g mean and "
         "%.3g cov)\n",
         shiftedWorstFit, shiftedWorstResidual,
         std::max(priorWorstFit, priorWorstResidual), levelDraws, levelTrees,
         worstSum, worstMeanZ, worstCovZ, poisonWorstSum, poisonMeanZ,
         poisonCovZ);
}
