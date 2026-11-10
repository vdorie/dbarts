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

// The level-fiber arm's own shape: small enough that the shift's 10 x 10
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
struct LevelDistribution {
  std::vector<double> mean, variance, covariance;  // covariance is m x m
};

LevelDistribution levelDistributionOf(ConstantLeafSampler& sampler,
                                      size_t numTrees) {
  ForestCalibration calibration = sampler.chain(0).forestCalibration(0);
  // priorScale is the response-unit total; the internal per-leaf sd is that
  // divided by the response transform and sqrt(m)
  double tau = calibration.priorSd /
               (calibration.responseScale *
                std::sqrt(static_cast<double>(numTrees)));

  LevelDistribution dist;
  dist.mean.assign(numTrees, 0.0);
  dist.variance.assign(numTrees, 0.0);
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
    dist.variance[t] = tau * tau / leafCount;
    dist.mean[t] = -leafSum / leafCount;
    sumVariance += dist.variance[t];
    sumMean += dist.mean[t];
  }
  double ratio = sumMean / sumVariance;
  dist.covariance.assign(numTrees * numTrees, 0.0);
  for (size_t t = 0; t < numTrees; ++t) {
    dist.mean[t] -= dist.variance[t] * ratio;
    for (size_t u = 0; u < numTrees; ++u)
      dist.covariance[t * numTrees + u] =
        (t == u ? dist.variance[t] : 0.0) -
        dist.variance[t] * dist.variance[u] / sumVariance;
  }
  return dist;
}

/// Empirical first and second moments of a stream of shift vectors, scored
/// against a distribution as per-coordinate z. The mean's standard error
/// is sqrt(C_tt / N); a gaussian sample covariance entry's is
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

  void score(const LevelDistribution& dist, double& worstMeanZ,
             double& worstCovZ) const {
    size_t m = sum.size();
    double n = static_cast<double>(numDraws);
    worstMeanZ = 0.0;
    worstCovZ = 0.0;
    std::vector<double> empiricalMean(m);
    for (size_t t = 0; t < m; ++t) empiricalMean[t] = sum[t] / n;
    for (size_t t = 0; t < m; ++t) {
      double se = std::sqrt(dist.covariance[t * m + t] / n);
      worstMeanZ = std::max(worstMeanZ,
                            std::fabs(empiricalMean[t] - dist.mean[t]) / se);
      for (size_t u = t; u < m; ++u) {
        double empirical =
          sumProducts[t * m + u] / n - empiricalMean[t] * empiricalMean[u];
        double target = dist.covariance[t * m + u];
        double se2 = std::sqrt((dist.covariance[t * m + t] *
                                  dist.covariance[u * m + u] +
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

// ---- the moves along the scale orbit: shared fixtures ----------------------
//
// The probit rescaling step and the collapsed move both multiply the occupied
// leaves by a factor and divide k by it, so one frozen-forest fixture, one
// mapping check, one table of declines and one every-chain check serve both.
// A slice step from a fixed point is not a draw from its target, so the
// distributional tests are in invariance form: start each replication from a
// point drawn exactly from the target, take one step, and score the result
// against the target by Kolmogorov-Smirnov. At 2e5 replications the p = 1e-3
// bound is 1.95 / sqrt(N).
constexpr size_t rescaleN = 243, rescaleP = 3, rescaleTrees = 8,
                 rescaleBurnIn = 150, rescaleDraws = 200000;
constexpr double rescaleKsBound = 1.95;
constexpr double rescaleMapTolerance = 1.0e-14;

/// A one-dimensional log density tabulated on a uniform grid, its CDF by the
/// trapezoid rule, sampled by inverse CDF and read back by interpolation.
struct GridTarget {
  double lower = 0.0, step = 0.0;
  std::vector<double> cdf;

  GridTarget() = default;

  /// The same table over a stated interval, for a target whose mass outside
  /// it the caller has bounded.
  template <typename F>
  static GridTarget onInterval(const F& logDensity, double from, double to,
                               size_t numPoints) {
    GridTarget target;
    target.lower = from;
    target.step = (to - from) / static_cast<double>(numPoints - 1);
    std::vector<double> value(numPoints);
    for (size_t j = 0; j < numPoints; ++j)
      value[j] = logDensity(from + target.step * static_cast<double>(j));
    double best = *std::max_element(value.begin(), value.end());
    target.cdf.assign(numPoints, 0.0);
    for (size_t j = 1; j < numPoints; ++j)
      target.cdf[j] = target.cdf[j - 1] +
        0.5 * target.step *
          (std::exp(value[j - 1] - best) + std::exp(value[j] - best));
    double total = target.cdf.back();
    for (double& c : target.cdf) c /= total;
    return target;
  }

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
    if (j >= cdf.size())
      return lower + step * static_cast<double>(cdf.size() - 1);
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

/// One rescaling fixture's data as the references read it, every vector full
/// length whatever the sampler was handed: zeros for no offset, ones for no
/// counts, no trials and no mask.
struct RescaleData {
  std::vector<double> x, y, offset, mask, counts, trials;

  /// The data of rescaleN rows with no offset, counts, trials or mask yet.
  RescaleData()
    : x(rescaleN * rescaleP), y(rescaleN), offset(rescaleN, 0.0),
      mask(rescaleN, 1.0), counts(rescaleN, 1.0), trials(rescaleN, 1.0) {}

  /// Whether row i is in the likelihood: kept by the mask, of a positive
  /// count and of at least one trial.
  bool in(size_t i) const {
    return mask[i] != 0.0 && counts[i] != 0.0 && trials[i] != 0.0;
  }
};

/// check, with the fixture's label ahead of what is checked.
void checkLabeled(bool condition, const char* label, const char* what) {
  char text[240];
  std::snprintf(text, sizeof text, "%s: %s", label, what);
  check(condition, text);
}

/// Whether the probit step declines and leaves the generator where it was.
template <typename C>
bool rescaleDeclinesWithoutDraw(C& chain) {
  ext_rng* rng = TestPeer::rng(chain);
  ext_rng* reference = cloneRng(rng);
  double alpha = 0.0;
  bool taken = TestPeer::drawForestRescale(chain, &alpha);
  bool agree = rngStreamsAgree(reference, rng);
  ext_rng_destroy(reference);
  return !taken && alpha == 1.0 && agree;
}

/// Whether the collapsed move declines, on its schedule and without it, and
/// leaves the generator where it was.
template <typename C>
bool collapsedDeclinesWithoutDraw(C& chain) {
  ext_rng* rng = TestPeer::rng(chain);
  ext_rng* reference = cloneRng(rng);
  double factor = 0.0;
  bool taken = TestPeer::maybeDrawCollapsedRescale(chain);
  taken = TestPeer::drawCollapsedRescale(chain, &factor) || taken;
  bool agree = rngStreamsAgree(reference, rng);
  ext_rng_destroy(reference);
  return !taken && factor == 1.0 && agree;
}

/// A drawn-or-fixed-k sampler on the rescaling fixtures' data. counts, trials
/// and offset may be null. rescale and collapsed are the probit step's and
/// the collapsed move's toggles.
std::unique_ptr<ConstantLeafSampler>
makeRescaleSampler(const std::vector<double>& x, const std::vector<double>& y,
                   const double* counts, const double* trials,
                   const double* offset, ResponseFamily family, bool updateK,
                   bool rescale, bool collapsed, ext_rng** rngs,
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
  options.collapsedRescaleForest = collapsed;
  options.trials = trials;
  return std::make_unique<ConstantLeafSampler>(
    x.data(), y.data(), rescaleN, rescaleP, counts, offset, family, 1.0, 3.0,
    0.37804942330213542, options, rngs);
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

/// A fixture chain's frozen forest and its moves along the scale orbit.
struct FrozenForest {
  Chain<ConstantGaussianLeaf>& chain;
  std::vector<std::vector<double>> mu;
  double k;
  /// The latents and the rows whose latents move with the leaves; both empty
  /// where a move leaves the latents alone.
  std::vector<double> latents, latentMask;
  /// the occupied leaves, as (tree, node) pairs
  std::vector<std::pair<size_t, int32_t>> occupied;
  size_t numEmpty = 0;

  /// scaledRows, where not null, names the rows whose latents the orbit
  /// scales with the leaves, as the probit step's does.
  explicit FrozenForest(Chain<ConstantGaussianLeaf>& c,
                        const std::vector<double>* scaledRows = nullptr)
    : chain(c), mu(rescaleTrees), k(TestPeer::forestK(c)) {
    if (scaledRows != nullptr) {
      latents.assign(TestPeer::latents(c), TestPeer::latents(c) + rescaleN);
      latentMask = *scaledRows;
    }
    std::vector<int32_t> bottoms;
    for (size_t t = 0; t < rescaleTrees; ++t) {
      mu[t] = TestPeer::muByTree(chain, t);
      const Tree& tree = chain.treeInForest(0, t);
      bottoms.clear();
      tree.fillBottom(0, bottoms);
      for (int32_t b : bottoms) {
        if (tree.at(b).numObservations() > 0) occupied.emplace_back(t, b);
        else ++numEmpty;
      }
    }
  }

  /// Puts the chain at the orbit point factor: the occupied leaves and the
  /// held latents of the scaled rows times it, k divided by it, the cache
  /// re-derived.
  void moveTo(double factor) {
    if (!latents.empty()) {
      std::vector<double> z(latents);
      for (size_t i = 0; i < rescaleN; ++i)
        if (latentMask[i] != 0.0) z[i] *= factor;
      TestPeer::restoreLatents(chain, z.data());
    }
    for (size_t t = 0; t < rescaleTrees; ++t)
      TestPeer::muByTree(chain, t) = mu[t];
    for (const auto& [t, b] : occupied)
      TestPeer::muByTree(chain, t)[static_cast<size_t>(b)] *= factor;
    TestPeer::forestK(chain) = k / factor;
    TestPeer::rebuildTotalFits(chain);
  }

  /// The frozen k against zeroed leaves, where every fit is exactly zero.
  void moveToZeroFit() {
    for (size_t t = 0; t < rescaleTrees; ++t)
      std::fill(TestPeer::muByTree(chain, t).begin(),
                TestPeer::muByTree(chain, t).end(), 0.0);
    TestPeer::forestK(chain) = k;
    TestPeer::rebuildTotalFits(chain);
  }

  /// The log factor the chain sits at, read off its leaves, its k and, where
  /// held, its latents.
  double logFactorOfLeaves() const {
    double moved = 0.0, frozen = 0.0;
    for (const auto& [t, b] : occupied) {
      moved += std::fabs(TestPeer::muByTree(chain, t)[static_cast<size_t>(b)]);
      frozen += std::fabs(mu[t][static_cast<size_t>(b)]);
    }
    return std::log(moved / frozen);
  }
  double logFactorOfK() const { return std::log(k / TestPeer::forestK(chain)); }
  double logFactorOfLatents() const {
    const double* z = TestPeer::latents(chain);
    double moved = 0.0, frozen = 0.0;
    for (size_t i = 0; i < rescaleN; ++i) {
      if (latentMask[i] == 0.0) continue;
      moved += std::fabs(z[i]);
      frozen += std::fabs(latents[i]);
    }
    return std::log(moved / frozen);
  }

  /// The chain's leaves gathered per row in tree order.
  std::vector<double> gatherTotals() const {
    std::vector<double> gathered(rescaleN, 0.0);
    for (size_t t = 0; t < rescaleTrees; ++t)
      for (size_t i = 0; i < rescaleN; ++i)
        gathered[i] +=
          TestPeer::muByTree(chain, t)[TestPeer::leafOf(chain, t)[i]];
    return gathered;
  }
};

/// One step from a frozen state at factor 1, and what it mapped.
struct OrbitStep {
  bool taken = false;
  double factor = 0.0;
  bool kMapped = false, leavesMapped = true, totalsMapped = false;
  std::vector<double> latentsBefore, workingBefore, weightsBefore;
};

/// Puts the frozen state at factor 1 and takes step, a callable
/// bool(double* factor) drawing one move, reading off whether k is divided
/// by the factor, every occupied leaf multiplied by it with an empty one
/// left as it was, and totalFits scaled in place with the factor accumulated
/// inside the bound and re-derived past it. The latents, the working
/// response and the precisions before the step are kept for the caller.
template <typename Step>
OrbitStep takeOrbitStep(ConstantLeafSampler& sampler, FrozenForest& frozen,
                        const Step& step) {
  Chain<ConstantGaussianLeaf>& chain = frozen.chain;
  frozen.moveTo(1.0);
  OrbitStep mapped;
  auto copyOf = [](const double* values) {
    return std::vector<double>(values, values + rescaleN);
  };
  std::vector<double> totalsBefore = forestTotals(sampler, 0);
  mapped.latentsBefore = copyOf(TestPeer::latents(chain));
  mapped.workingBefore = copyOf(TestPeer::workingResponse(chain));
  mapped.weightsBefore = copyOf(TestPeer::workingWeights(chain));
  mapped.taken = step(&mapped.factor);
  double factor = mapped.factor;
  mapped.kMapped =
    std::fabs(TestPeer::forestK(chain) * factor / frozen.k - 1.0) <
    rescaleMapTolerance;
  std::vector<int32_t> bottoms;
  for (size_t t = 0; t < rescaleTrees; ++t) {
    const Tree& tree = chain.treeInForest(0, t);
    bottoms.clear();
    tree.fillBottom(0, bottoms);
    for (int32_t b : bottoms) {
      double now = TestPeer::muByTree(chain, t)[static_cast<size_t>(b)];
      double was = frozen.mu[t][static_cast<size_t>(b)];
      mapped.leavesMapped = mapped.leavesMapped &&
        (tree.at(b).numObservations() > 0
           ? std::fabs(now - factor * was) <=
               rescaleMapTolerance * std::fabs(now)
           : now == was);
    }
  }
  std::vector<double> totals = forestTotals(sampler, 0),
                      gathered = frozen.gatherTotals();
  bool inPlace = factor >= 1.0 / totalFitsScaleBound &&
                 factor <= totalFitsScaleBound;
  mapped.totalsMapped =
    TestPeer::totalFitsScale(chain) == (inPlace ? factor : 1.0);
  for (size_t i = 0; i < rescaleN; ++i)
    mapped.totalsMapped = mapped.totalsMapped &&
      (inPlace ? totals[i] == totalsBefore[i] * factor &&
                   std::fabs(totals[i] - gathered[i]) <=
                     rescaleMapTolerance * (1.0 + std::fabs(gathered[i]))
               : totals[i] == gathered[i]);
  return mapped;
}

/// The checks on what takeOrbitStep read off, then one more step with the
/// cache's factor forced past its bound: the step re-derives the cache from
/// the scaled leaves in tree order, bitwise, and resets the factor.
template <typename Step>
void checkOrbitStep(const char* label, const OrbitStep& mapped,
                    ConstantLeafSampler& sampler, FrozenForest& frozen,
                    const Step& step) {
  frozen.moveTo(1.0);
  TestPeer::totalFitsScale(frozen.chain) = 1.0e6;
  double forced = 0.0;
  step(&forced);
  bool totalsRederived = forestTotals(sampler, 0) == frozen.gatherTotals() &&
                         TestPeer::totalFitsScale(frozen.chain) == 1.0;
  checkLabeled(mapped.taken && mapped.factor != 1.0, label,
               "the frozen state takes the step");
  checkLabeled(mapped.kMapped, label, "k is divided by the factor");
  checkLabeled(mapped.leavesMapped, label,
               "every occupied leaf of every tree is multiplied by the "
               "factor, an empty one left as it was");
  checkLabeled(mapped.totalsMapped, label,
               "totalFits is scaled in place with the factor accumulated "
               "inside the bound, re-derived past it");
  checkLabeled(totalsRederived, label,
               "past the bound totalFits is re-derived from the leaves in "
               "tree order, bitwise, and the factor reset");
}

/// The fixtures every move along the orbit declines on with no generator
/// draw, the move's toggles on throughout: its own family at a fixed k and
/// on a stale tree map, every other family it might be built beside, a
/// variance forest, a linear leaf and two forests that share no total.
/// declines is a callable on a chain of either leaf, true where the move
/// declines without a draw.
template <typename Declines>
void checkOrbitDeclines(const char* move, const RescaleData& data,
                        ResponseFamily eligible, ext_rng* rng,
                        const Declines& declines) {
  Results results;
  char what[200];
  auto report = [&](bool declined, const char* where) {
    std::snprintf(what, sizeof what, "%s declines %s", move, where);
    check(declined, what);
  };
  auto build = [&](ResponseFamily family, bool updateK,
                   size_t numVarianceTrees) {
    std::unique_ptr<ConstantLeafSampler> sampler = makeRescaleSampler(
      data.x, data.y, nullptr, nullptr, data.offset.data(), family, updateK,
      true, true, &rng, 1, numVarianceTrees);
    sampler->run(20, 0, results);
    return sampler;
  };
  report(declines(build(eligible, false, 0)->chain(0)), "at a fixed k");
  {
    std::unique_ptr<ConstantLeafSampler> stale = build(eligible, true, 0);
    stale->sampleTreesFromPrior();
    report(TestPeer::leafOfStale(stale->chain(0), 0) != 0 &&
             declines(stale->chain(0)),
           "on a stale tree map");
  }
  const std::pair<ResponseFamily, const char*> families[] = {
    {ResponseFamily::probit, "under the probit family"},
    {ResponseFamily::logistic, "under the logistic family"},
    {ResponseFamily::gaussian, "under the gaussian family"},
    {ResponseFamily::nbinom, "under the nbinom family"}};
  for (const auto& [family, where] : families)
    if (family != eligible)
      report(declines(build(family, true, 0)->chain(0)), where);
  {
    std::unique_ptr<ConstantLeafSampler> variance =
      build(ResponseFamily::gaussian, true, 4);
    report(variance->chain(0).hasVarianceForest() &&
             declines(variance->chain(0)),
           "beside a variance forest");
  }

  // an ordinal fit, a linear leaf and a second forest, all drawing k
  SamplerOptions drawnK;
  drawnK.numTrees = rescaleTrees;
  drawnK.nodeScale = 3.0;
  drawnK.updateK = true;
  drawnK.kHyperprior.degreesOfFreedom = 1.5;
  drawnK.kHyperprior.scale = 2.0;
  std::vector<double> category(rescaleN), z(rescaleN);
  for (size_t i = 0; i < rescaleN; ++i) {
    category[i] = 1.0 + data.y[i] + (data.x[i] > 0.7 ? 1.0 : 0.0);
    z[i] = i % 2 == 0 ? 1.0 : 0.0;
  }
  SamplerOptions options(drawnK);
  options.numCategories = 3;
  ConstantLeafSampler ordinal(data.x.data(), category.data(), rescaleN,
                              rescaleP, nullptr, data.offset.data(),
                              ResponseFamily::ordinal, 1.0, 3.0,
                              0.37804942330213542, options, &rng);
  ordinal.run(20, 0, results);
  report(declines(ordinal.chain(0)), "under the ordinal family");
  std::vector<size_t> covariates = {2};
  options = drawnK;
  options.leafCovariateColumns = covariates.data();
  options.numLeafCovariates = 1;
  Sampler<LinearGaussianLeaf> linear(data.x.data(), data.y.data(), rescaleN,
                                     rescaleP, nullptr, data.offset.data(),
                                     eligible, 1.0, 3.0, 0.37804942330213542,
                                     options, &rng);
  linear.run(20, 0, results);
  report(declines(linear.chain(0)), "on a linear leaf");
  AmplitudeSpec spec;
  spec.family = eligible;
  spec.mu.numTrees = rescaleTrees;
  spec.tau.numTrees = 4;
  spec.z = z.data();
  ConstantLeafSampler twoForests(data.x.data(), data.y.data(), rescaleN,
                                 rescaleP, nullptr, data.offset.data(), 1.0,
                                 3.0, 1.0, drawnK, spec, &rng);
  twoForests.run(20, 0, results);
  // bcf's two forests share no total, mu stating its sd and tau carrying a
  // basis, so its chain holds k fixed; marked drawn, only the forest count
  // and the combiner with no sharers are left to decline on
  Chain<ConstantGaussianLeaf>& twoForestChain = twoForests.chain(0);
  TestPeer::forestUpdateK(twoForestChain) = true;
  TestPeer::kHyperprior(twoForestChain) = drawnK.kHyperprior;
  report(twoForestChain.numForests() == 2 &&
           TestPeer::forestK(twoForestChain) > 0.0 &&
           TestPeer::forestK(twoForestChain) <= DBL_MAX &&
           declines(twoForestChain),
         "on a two-forest fit");
}

/// A move acts on every chain: two chains, each against a twin on its own
/// seed with the move off, each chain's k parting from its twin's after
/// numSweeps. make(on, rngs) builds the two-chain sampler with the move on
/// or off.
template <typename Make>
void checkActsOnEveryChain(const char* move, std::uint32_t seed,
                           size_t numSweeps, const Make& make) {
  std::vector<ext_rng*> onRngs(2), offRngs(2);
  for (size_t c = 0; c < 2; ++c) {
    onRngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    offRngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(onRngs[c], seed + static_cast<std::uint32_t>(c));
    ext_rng_setSeed(offRngs[c], seed + static_cast<std::uint32_t>(c));
  }
  {
    std::unique_ptr<ConstantLeafSampler> on = make(true, onRngs.data());
    std::unique_ptr<ConstantLeafSampler> off = make(false, offRngs.data());
    Results results;
    on->run(numSweeps, 0, results);
    off->run(numSweeps, 0, results);
    char what[200];
    for (size_t c = 0; c < 2; ++c) {
      std::snprintf(what, sizeof what, "%s acts on chain %zu", move, c);
      check(TestPeer::forestK(on->chain(c)) != TestPeer::forestK(off->chain(c)),
            what);
    }
  }
  for (size_t c = 0; c < 2; ++c) {
    ext_rng_destroy(onRngs[c]);
    ext_rng_destroy(offRngs[c]);
  }
}

// ---- the probit rescaling step ---------------------------------------------

void runForestRescaleTests() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20261008u);
  const char* label = "rescaling";

  // every fifth row inactive, so the leaves hold masked members
  RescaleData data;
  for (double& v : data.x) v = runif01();
  for (size_t i = 0; i < rescaleN; ++i) {
    double eta = 1.5 * std::sin(3.0 * data.x[i]) + data.x[i + rescaleN] - 0.5;
    data.y[i] = eta + standardNormal() > 0.0 ? 1.0 : 0.0;
    data.offset[i] = 0.8 * (data.x[i + 2 * rescaleN] - 0.5);
  }
  for (size_t i = 0; i < rescaleN; i += 5) data.mask[i] = 0.0;

  std::unique_ptr<ConstantLeafSampler> sampler = makeRescaleSampler(
    data.x, data.y, nullptr, nullptr, data.offset.data(),
    ResponseFamily::probit, true, true, false, &rng);
  Results results;
  sampler->run(rescaleBurnIn, 0, results);
  check(sampler->setActiveRows(data.mask.data()),
        "the probit sampler takes a mask");
  sampler->run(20, 0, results);
  Chain<ConstantGaussianLeaf>& chain = sampler->chain(0);
  // no sweep leaves an empty leaf, so one is stranded; its value is nonzero so
  // that scaling it would show
  check(TestPeer::strandEmptyLeaf(chain, 0.5).second >= 0,
        "rescaling: the fixture strands an empty leaf");
  check(TestPeer::forestRescaleApplies(chain),
        "a burned-in drawn-k probit forest takes the rescaling step");

  // the frozen state, its latents moving with the leaves on the rows in
  FrozenForest frozen(chain, &data.mask);
  check(frozen.occupied.size() > rescaleTrees && frozen.numEmpty > 0,
        "the frozen trees carry splits and an empty leaf");
  frozen.moveTo(1.0);
  RescaleConstants constants = rescaleConstantsOf(*sampler, data.mask);
  char what[200];

  // 1. The conditional, in invariance form. v0 is the orbit coordinate
  // measured from the frozen state, so its target is the frozen state's own
  // conditional; one step from e^v0 then lands at v1, read off the latents and
  // the leaves, which must agree with each other, with alpha and with k.
  GridTarget target([&](double v) { return constants.logDensity(v); }, -6.0,
                    6.0, 200001);
  std::vector<double> landed(rescaleDraws);
  double worstDisagreement = 0.0;
  for (size_t r = 0; r < rescaleDraws; ++r) {
    double start = target.draw(ext_rng_simulateContinuousUniform(rng));
    frozen.moveTo(std::exp(start));
    double alpha = 0.0;
    TestPeer::drawForestRescale(chain, &alpha);
    double fromLatents = frozen.logFactorOfLatents();
    double fromAlpha = start + std::log(alpha);
    worstDisagreement = std::max(
      {worstDisagreement, std::fabs(frozen.logFactorOfLeaves() - fromLatents),
       std::fabs(frozen.logFactorOfK() - fromLatents),
       std::fabs(fromAlpha - fromLatents)});
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

  // 2. The mapping, one step from the frozen state: the leaves, k and the
  // cache, and the active latents with them.
  auto step = [&](double* alpha) {
    return TestPeer::drawForestRescale(chain, alpha);
  };
  OrbitStep mapped = takeOrbitStep(*sampler, frozen, step);
  const double* z = TestPeer::latents(chain);
  const double* working = TestPeer::workingResponse(chain);
  bool latentsMapped = true, workingMapped = true;
  for (size_t i = 0; i < rescaleN; ++i) {
    latentsMapped = latentsMapped &&
      (data.mask[i] != 0.0
         ? std::fabs(z[i] - mapped.factor * frozen.latents[i]) <=
             rescaleMapTolerance * std::fabs(z[i])
         : z[i] == frozen.latents[i]);
    workingMapped = workingMapped && working[i] == z[i] - data.offset[i];
  }
  checkLabeled(latentsMapped, label,
               "the active latents are multiplied by the factor, the "
               "inactive ones left as they were");
  checkLabeled(workingMapped, label,
               "the working response is the scaled latents less the offset");
  checkOrbitStep(label, mapped, *sampler, frozen, step);

  // 3. Inertness: each decline consumes no generator draw.
  frozen.moveTo(1.0);
  TestPeer::forestK(chain) = HUGE_VAL;
  check(rescaleDeclinesWithoutDraw(chain), "rescaling declines at infinite k");
  frozen.moveTo(1.0);
  ChiKHyperprior savedPrior = TestPeer::kHyperprior(chain);
  TestPeer::kHyperprior(chain).scale = HUGE_VAL;
  TestPeer::kHyperprior(chain).degreesOfFreedom = constants.numActive;
  check(rescaleDeclinesWithoutDraw(chain),
        "rescaling declines under an infinite k scale with n <= nu");
  TestPeer::kHyperprior(chain).degreesOfFreedom = constants.numActive - 1.0;
  double alpha = 0.0;
  check(TestPeer::drawForestRescale(chain, &alpha),
        "rescaling is taken under an infinite k scale with n > nu");
  TestPeer::kHyperprior(chain) = savedPrior;
  frozen.moveTo(1.0);
  std::vector<double> none(rescaleN, 0.0);
  check(sampler->setActiveRows(none.data()), "an all-inactive mask installs");
  check(rescaleDeclinesWithoutDraw(chain),
        "rescaling declines with every row inactive");
  check(sampler->setActiveRows(data.mask.data()), "the mask reinstalls");

  ext_rng* otherRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(otherRng, 20261009u);
  checkOrbitDeclines(label, data, ResponseFamily::probit, otherRng,
                     [](auto& c) { return rescaleDeclinesWithoutDraw(c); });
  ext_rng_destroy(otherRng);
  checkActsOnEveryChain(label, 11u, 30, [&](bool on, ext_rng** rngs) {
    return makeRescaleSampler(data.x, data.y, nullptr, nullptr,
                              data.offset.data(), ResponseFamily::probit, true,
                              on, false, rngs, 2);
  });

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

  sampler.reset();
  ext_rng_destroy(rng);
  printf("ok: probit rescaling, %zu steps at %zu trees (sqrt(N) D %.3g, worst "
         "factor disagreement %.3g, %zu empty leaves); slice sqrt(N) D %.3g "
         "shipped, %.3g limit binding\n",
         rescaleDraws, rescaleTrees, ks, worstDisagreement, frozen.numEmpty,
         ksShipped, ksBinding);
}

// ---- the collapsed rescaling move ------------------------------------------
//
// The move draws its factor with the latents integrated out and leaves them
// to the augmentation draw that follows it in a sweep, so everything it is
// scored against is read off the raw state: the response, the offset, the
// mask and the counts as the test holds them, the leaves gathered in tree
// order and the k prior, never the move's own arithmetic. The distributional
// tests are in invariance form, as above, at a quarter of the replications:
// each one is a slice step over the rows' likelihood, and the p = 1e-3 bound
// is 1.95 / sqrt(N) at any N.
constexpr size_t collapsedDraws = 50000, collapsedRedraws = 2000,
                 collapsedCoinCalls = 4000;
// collapsedRowZBound holds the worst of some two hundred rows' means of a
// right-skewed variate, where one mean's 4.5 would fail a correct sampler on
// one stream in a few hundred.
constexpr double collapsedZBound = 4.5, collapsedRowZBound = 5.5;

/// The move's target at one frozen state: nu, C = k^2 / (2 s^2), and for each
/// row in the likelihood its fit, offset, response, trials and count.
struct CollapsedTarget {
  double degreesOfFreedom = 0.0, kTerm = 0.0;
  std::vector<double> fit, offset, response, trials, count;

  double logDensity(double v) const {
    double factor = std::exp(v), sum = 0.0;
    for (size_t j = 0; j < fit.size(); ++j)
      sum += logisticRowLogLikelihood(response[j], trials[j], count[j],
                                      factor * fit[j] + offset[j]);
    return -degreesOfFreedom * v - kTerm * std::exp(-2.0 * v) + sum;
  }

  /// Its table over an interval outside which the density is below e^-45 of
  /// its value at 0. No row's likelihood passes 1, so left of 0 the log
  /// density is at most -nu v - C e^(-2v), which falls without bound once it
  /// turns, and right of 0 at most -nu v; where the data are separated the
  /// right tail is no thinner than that.
  GridTarget table() const {
    double atZero = logDensity(0.0), lower = 0.0;
    while (2.0 * kTerm * std::exp(-2.0 * lower) <= degreesOfFreedom ||
           -degreesOfFreedom * lower - kTerm * std::exp(-2.0 * lower) >
             atZero - 45.0)
      lower -= 0.25;
    return GridTarget::onInterval([&](double v) { return logDensity(v); },
                                  lower, (45.0 - atZero) / degreesOfFreedom,
                                  100001);
  }
};

/// The move's target at the frozen chain's state under data.
CollapsedTarget collapsedTargetOf(const FrozenForest& frozen,
                                  const RescaleData& data) {
  CollapsedTarget target;
  std::vector<double> fit = frozen.gatherTotals();
  for (size_t i = 0; i < rescaleN; ++i) {
    if (!data.in(i)) continue;
    target.fit.push_back(fit[i]);
    target.offset.push_back(data.offset[i]);
    target.response.push_back(data.y[i]);
    target.trials.push_back(data.trials[i]);
    target.count.push_back(data.counts[i]);
  }
  const ChiKHyperprior& prior = TestPeer::kHyperprior(frozen.chain);
  double chainK = TestPeer::forestK(frozen.chain);
  target.kTerm = chainK * chainK / (2.0 * prior.scale * prior.scale);
  target.degreesOfFreedom = prior.degreesOfFreedom;
  return target;
}

/// The conditional in invariance form, against the target of the frozen
/// state under data: per replication a start drawn from it, the state moved
/// to that orbit point, one move, and the landing point read off the leaves,
/// which k and the reported factor must agree with to within
/// worstDisagreement. Returns sqrt(N) times the KS distance of the landing
/// points from the target, HUGE_VAL if a move declined.
double collapsedInvarianceKs(FrozenForest& frozen, const RescaleData& data,
                             size_t numDraws, double& worstDisagreement) {
  frozen.moveTo(1.0);
  GridTarget target = collapsedTargetOf(frozen, data).table();
  ext_rng* rng = TestPeer::rng(frozen.chain);
  std::vector<double> landed(numDraws);
  bool everyTaken = true;
  for (size_t r = 0; r < numDraws; ++r) {
    double start = target.draw(ext_rng_simulateContinuousUniform(rng));
    frozen.moveTo(std::exp(start));
    double factor = 0.0;
    everyTaken =
      TestPeer::drawCollapsedRescale(frozen.chain, &factor) && everyTaken;
    landed[r] = frozen.logFactorOfLeaves();
    worstDisagreement = std::max(
      {worstDisagreement, std::fabs(frozen.logFactorOfK() - landed[r]),
       std::fabs(start + std::log(factor) - landed[r])});
  }
  if (!everyTaken) return HUGE_VAL;
  return scaledKsDistance(std::move(landed),
                          [&](double v) { return target.at(v); });
}

/// Test 1 at one state: the evaluator against the row-order sum of the
/// reported per-row log-likelihood over the rows in, and against the
/// reference's own sum, at several factors, with nothing held beside the
/// scaled fits and with a part held unscaled, as the sharers' split hands it.
void checkCollapsedEvaluator(const char* label,
                             Chain<ConstantGaussianLeaf>& chain,
                             const std::vector<double>& fits,
                             const RescaleData& data) {
  const ResponseModel& response = TestPeer::response(chain);
  size_t n = fits.size();
  std::vector<double> scaled(n), reported(n), held(n);
  for (size_t i = 0; i < n; ++i)
    held[i] = 0.25 * std::sin(static_cast<double>(i));
  double worstReported = 0.0, worstOwn = 0.0;
  bool rowsOutReported = true;
  for (double factor : {0.25, 1.0, 1.7, 6.0}) {
    for (size_t i = 0; i < n; ++i) scaled[i] = factor * fits[i];
    response.computeLogLikelihood(scaled.data(), 1.0, n, reported.data());
    double sum = 0.0, own = 0.0, ownHeld = 0.0;
    for (size_t i = 0; i < n; ++i) {
      // a row out of the model reports NaN, a binomial row of no trials 0
      if (data.mask[i] == 0.0 || data.counts[i] == 0.0) {
        rowsOutReported = rowsOutReported && std::isnan(reported[i]);
        continue;
      }
      if (data.trials[i] == 0.0) {
        rowsOutReported = rowsOutReported && reported[i] == 0.0;
        continue;
      }
      sum += reported[i];
      own += logisticRowLogLikelihood(data.y[i], data.trials[i], data.counts[i],
                                      factor * fits[i] + data.offset[i]);
      ownHeld += logisticRowLogLikelihood(
        data.y[i], data.trials[i], data.counts[i],
        factor * fits[i] + held[i] + data.offset[i]);
    }
    double got =
      response.logLikelihoodAtScaledFits(fits.data(), nullptr, factor);
    double gotHeld =
      response.logLikelihoodAtScaledFits(fits.data(), held.data(), factor);
    worstReported = std::max(worstReported, std::fabs(got / sum - 1.0));
    worstOwn = std::max({worstOwn, std::fabs(got / own - 1.0),
                         std::fabs(gotHeld / ownHeld - 1.0)});
  }
  char what[200];
  std::snprintf(what, sizeof what,
                "%s: the evaluator is the reported log-likelihood summed in "
                "row order over the rows in (worst %.3g relative)", label,
                worstReported);
  check(worstReported < 1.0e-14 && rowsOutReported, what);
  std::snprintf(what, sizeof what,
                "%s: the evaluator reads the response, offset, mask, "
                "counts and trials in force, and a part held unscaled beside "
                "the fits "
                "(worst %.3g relative)", label, worstOwn);
  check(worstOwn < 1.0e-12, what);
}

/// Test 3, the move as a sweep takes it: over many moves from one state, each
/// followed by the augmentation draw, each row's omega standardized by
/// PG(w_i n_i, psi_i)'s moments at the moved fit. Scored are the mean over
/// every row and move, the same mean signed by whether the factor passed 1,
/// which a draw at the fit before the move shifts one way throughout, and the
/// worst row's mean; a row out keeps its omega, and the working response and
/// the precisions follow the draw. zeroFit runs it from zeroed leaves, where a
/// row with no offset sits at psi = 0 exactly.
void checkAugmentationAfterMove(const char* label, FrozenForest& frozen,
                                const RescaleData& data, bool zeroFit) {
  Chain<ConstantGaussianLeaf>& chain = frozen.chain;
  if (zeroFit) frozen.moveToZeroFit();
  else frozen.moveTo(1.0);
  std::vector<double> fit = frozen.gatherTotals();
  const double* omega = TestPeer::latents(chain);
  std::vector<double> before(omega, omega + rescaleN), rowSum(rescaleN, 0.0);
  double sum = 0.0, signedSum = 0.0, numTerms = 0.0;
  size_t numAtZero = 0;
  bool rowsOutKept = true;
  for (size_t m = 0; m < collapsedRedraws; ++m) {
    if (zeroFit) frozen.moveToZeroFit();
    else frozen.moveTo(1.0);
    double factor = 0.0;
    TestPeer::drawCollapsedRescale(chain, &factor);
    TestPeer::drawAugmentation(chain);
    for (size_t i = 0; i < rescaleN; ++i) {
      if (!data.in(i)) {
        rowsOutKept = rowsOutKept && omega[i] == before[i];
        continue;
      }
      double psi = factor * fit[i] + data.offset[i], mean, variance;
      if (m == 0 && psi == 0.0) ++numAtZero;
      polyaGammaMoments(data.counts[i] * data.trials[i], psi, mean, variance);
      double standardized = (omega[i] - mean) / std::sqrt(variance);
      rowSum[i] += standardized;
      sum += standardized;
      signedSum += factor > 1.0 ? standardized : -standardized;
      numTerms += 1.0;
    }
  }
  double worstRowZ = 0.0;
  for (size_t i = 0; i < rescaleN; ++i)
    worstRowZ = std::max(worstRowZ, std::fabs(rowSum[i]) /
                           std::sqrt(static_cast<double>(collapsedRedraws)));
  double z = sum / std::sqrt(numTerms),
         signedZ = signedSum / std::sqrt(numTerms);
  char what[240];
  std::snprintf(what, sizeof what,
                "%s: after the move the augmentation draw takes each omega as "
                "PG(w n, psi) at the moved fit (z %.3g, signed by the factor "
                "%.3g, worst row %.3g, %zu rows at psi = 0)", label, z, signedZ,
                worstRowZ, numAtZero);
  check(std::fabs(z) < collapsedZBound &&
          std::fabs(signedZ) < collapsedZBound &&
          worstRowZ < collapsedRowZBound && rowsOutKept &&
          (numAtZero > 0) == zeroFit, what);
  const double* working = TestPeer::workingResponse(chain);
  const double* weights = TestPeer::workingWeights(chain);
  bool workingFollows = true;
  for (size_t i = 0; i < rescaleN; ++i) {
    if (!data.in(i)) continue;
    double expected = data.counts[i] * (data.y[i] - 0.5 * data.trials[i]) /
                        omega[i] - data.offset[i];
    workingFollows = workingFollows &&
      std::fabs(working[i] - expected) <=
        1.0e-12 * (1.0 + std::fabs(expected)) &&
      weights[i] == omega[i];
  }
  checkLabeled(workingFollows, label,
               "the working response and precisions follow the draw after "
               "the move");
}

/// Tests 1 to 4 on one burned-in fixture; returns the conditional's
/// sqrt(N) D.
double runCollapsedStateTests(const char* label, ConstantLeafSampler& sampler,
                              const RescaleData& data,
                              double& worstDisagreement) {
  Chain<ConstantGaussianLeaf>& chain = sampler.chain(0);
  FrozenForest frozen(chain);
  checkCollapsedEvaluator(label, chain, forestTotals(sampler, 0), data);
  double ks =
    collapsedInvarianceKs(frozen, data, collapsedDraws, worstDisagreement);
  char what[240];
  std::snprintf(what, sizeof what,
                "%s: one move from the target stays at the target (sqrt(N) D "
                "%.3g against %.3g), leaves, k and the factor moving as one "
                "(worst %.3g)", label, ks, rescaleKsBound, worstDisagreement);
  check(ks < rescaleKsBound && worstDisagreement < 1.0e-10, what);
  checkAugmentationAfterMove(label, frozen, data, false);

  // 4. The mapping, one move from the frozen state: the leaves, k and the
  // cache. The move reads no latent, so it leaves the latents, the working
  // response and the precisions as they were.
  auto step = [&](double* factor) {
    return TestPeer::drawCollapsedRescale(chain, factor);
  };
  OrbitStep mapped = takeOrbitStep(sampler, frozen, step);
  auto unchanged = [](const std::vector<double>& before, const double* now) {
    return std::equal(before.begin(), before.end(), now);
  };
  checkLabeled(unchanged(mapped.latentsBefore, TestPeer::latents(chain)) &&
                 unchanged(mapped.workingBefore,
                           TestPeer::workingResponse(chain)) &&
                 unchanged(mapped.weightsBefore,
                           TestPeer::workingWeights(chain)),
               label,
               "the move leaves the latents, the working response and the "
               "precisions as they were");
  checkOrbitStep(label, mapped, sampler, frozen, step);
  frozen.moveTo(1.0);
  return ks;
}

void runCollapsedRescaleTests() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20261011u);
  const ResponseFamily family = ResponseFamily::logistic;
  const char* link = "collapsed logistic";
  char what[240], label[80];

  // the fixture: an offset that is exactly zero on every seventh row, every
  // fifth row out of the likelihood, counts 1 to 3
  RescaleData data;
  for (double& v : data.x) v = runif01();
  for (size_t i = 0; i < rescaleN; ++i) {
    double eta = 1.5 * std::sin(3.0 * data.x[i]) + data.x[i + rescaleN] - 0.5;
    data.y[i] = eta + standardNormal() > 0.0 ? 1.0 : 0.0;
    data.offset[i] =
      i % 7 == 0 ? 0.0 : 0.8 * (data.x[i + 2 * rescaleN] - 0.5);
    if (i % 5 == 0) data.mask[i] = 0.0;
    data.counts[i] = static_cast<double>(1 + i % 3);
  }
  std::unique_ptr<ConstantLeafSampler> sampler = makeRescaleSampler(
    data.x, data.y, data.counts.data(), nullptr, data.offset.data(), family,
    true, true, true, &rng);
  Results results;
  sampler->run(rescaleBurnIn, 0, results);
  check(sampler->setActiveRows(data.mask.data()),
        "the collapsed fixture takes a mask");
  sampler->run(20, 0, results);
  Chain<ConstantGaussianLeaf>& chain = sampler->chain(0);
  check(TestPeer::strandEmptyLeaf(chain, 0.5).second >= 0,
        "the collapsed fixture strands an empty leaf");
  checkLabeled(TestPeer::collapsedRescaleApplies(chain), link,
               "a burned-in drawn-k forest takes the move");

  // 1 to 4 on the fixture, and the redraw again from a zero fit
  double worstDisagreement = 0.0;
  double ks = runCollapsedStateTests(link, *sampler, data, worstDisagreement);
  FrozenForest frozen(chain);
  check(frozen.occupied.size() > rescaleTrees && frozen.numEmpty > 0,
        "the collapsed fixture's trees carry splits and an empty leaf");
  std::snprintf(label, sizeof label, "%s at a zero fit", link);
  checkAugmentationAfterMove(label, frozen, data, true);
  frozen.moveTo(1.0);

  // 1, continued. A row out given a fit of 1e6 against its response adds
  // nothing; and on a second sampler the evaluator follows setResponse and
  // setData, so it reads the response in force.
  ext_rng* otherRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(otherRng, 20261012u);
  {
    const ResponseModel& response = TestPeer::response(chain);
    std::vector<double> fits = forestTotals(*sampler, 0), poisoned(fits);
    for (size_t i = 0; i < rescaleN; ++i)
      if (data.mask[i] == 0.0) poisoned[i] = data.y[i] != 0.0 ? -1.0e6 : 1.0e6;
    checkLabeled(
      response.logLikelihoodAtScaledFits(poisoned.data(), nullptr, 1.7) ==
        response.logLikelihoodAtScaledFits(fits.data(), nullptr, 1.7),
      link, "a row out of the likelihood adds nothing to the evaluator's sum");

    std::unique_ptr<ConstantLeafSampler> other = makeRescaleSampler(
      data.x, data.y, data.counts.data(), nullptr, data.offset.data(),
      family, true, true, true, &otherRng);
    other->run(10, 0, results);
    RescaleData flipped(data);
    flipped.mask.assign(rescaleN, 1.0);
    for (double& v : flipped.y) v = 1.0 - v;
    check(other->setResponse(flipped.y.data(), false) ==
            ResponseUpdateResult::accepted,
          "the collapsed fixture takes a new response");
    std::snprintf(label, sizeof label, "%s after setResponse", link);
    checkCollapsedEvaluator(label, other->chain(0), forestTotals(*other, 0),
                            flipped);
    const size_t fewer = 120, first = 60;
    RescaleData replaced;
    replaced.x.resize(fewer * rescaleP);
    replaced.y.clear();
    replaced.offset.clear();
    replaced.counts.clear();
    replaced.mask.assign(fewer, 1.0);
    replaced.trials.assign(fewer, 1.0);
    for (size_t i = 0; i < fewer; ++i) {
      for (size_t j = 0; j < rescaleP; ++j)
        replaced.x[j * fewer + i] = data.x[j * rescaleN + first + i];
      replaced.y.push_back(i % 2 == 0 ? data.y[first + i]
                                      : 1.0 - data.y[first + i]);
      replaced.offset.push_back(-data.offset[first + i]);
      replaced.counts.push_back(static_cast<double>(1 + (i + 1) % 3));
    }
    check(other->setData(replaced.x.data(), replaced.y.data(), fewer,
                         replaced.counts.data(), replaced.offset.data(),
                         nullptr, 0),
          "the collapsed fixture takes new rows");
    other->run(5, 0, results);
    std::snprintf(label, sizeof label, "%s after setData", link);
    checkCollapsedEvaluator(label, other->chain(0), forestTotals(*other, 0),
                            replaced);
  }

  // 5. Inertness: each decline consumes no generator draw, on the schedule
  // and off it.
  TestPeer::collapsedRescaleForest(chain) = false;
  check(collapsedDeclinesWithoutDraw(chain),
        "the collapsed move declines with its toggle off");
  TestPeer::collapsedRescaleForest(chain) = true;
  TestPeer::forestK(chain) = HUGE_VAL;
  check(collapsedDeclinesWithoutDraw(chain),
        "the collapsed move declines at infinite k");
  TestPeer::forestK(chain) = 0.0;
  check(collapsedDeclinesWithoutDraw(chain),
        "the collapsed move declines at a k of zero");
  frozen.moveTo(1.0);
  const ChiKHyperprior savedPrior = TestPeer::kHyperprior(chain);
  TestPeer::kHyperprior(chain).scale = HUGE_VAL;
  check(collapsedDeclinesWithoutDraw(chain),
        "the collapsed move declines under an infinite k scale");
  // a scale so small that the prior term overflows at the state itself: the
  // forest is eligible and the move declines with no draw of its own
  TestPeer::kHyperprior(chain).scale = 1.0e-200;
  TestPeer::forestK(chain) = 1.0;
  {
    ext_rng* reference = cloneRng(rng);
    double factor = 0.0;
    bool taken = TestPeer::drawCollapsedRescale(chain, &factor);
    check(TestPeer::collapsedRescaleApplies(chain) && !taken &&
            factor == 1.0 && rngStreamsAgree(reference, rng),
          "the collapsed move declines without a draw where its log density "
          "is not finite at the state");
    ext_rng_destroy(reference);
  }
  TestPeer::kHyperprior(chain) = savedPrior;
  // a k whose square underflows: eligible, and its move finite. From a zero
  // fit the k prior's term alone bounds the orbit on the side k grows on,
  // and the slice is spent before 0.5 (k / s)^2 passes 1.5 log(k / 1e-170)
  // by more than its depth: one move leaves k under 100. A term that had
  // underflowed would leave that side open, and one move in forty past it.
  frozen.moveTo(1.0);
  TestPeer::forestK(chain) = 1.0e-170;
  double tinyFactor = 0.0;
  bool tinyTaken = TestPeer::collapsedRescaleApplies(chain) &&
                   TestPeer::drawCollapsedRescale(chain, &tinyFactor);
  bool tinyFinite = tinyFactor > 0.0 && tinyFactor <= DBL_MAX &&
                    TestPeer::forestK(chain) > 0.0 &&
                    TestPeer::forestK(chain) <= DBL_MAX;
  for (double v : forestTotals(*sampler, 0))
    tinyFinite = tinyFinite && std::isfinite(v);
  double tinyWorstK = 0.0;
  for (size_t m = 0; m < 400; ++m) {
    frozen.moveToZeroFit();
    TestPeer::forestK(chain) = 1.0e-170;
    TestPeer::drawCollapsedRescale(chain, nullptr);
    tinyWorstK = std::max(tinyWorstK, TestPeer::forestK(chain));
  }
  std::snprintf(what, sizeof what,
                "%s: a k whose square underflows takes a finite move, and "
                "from a zero fit no move passes its prior's reach (largest k "
                "%.3g)", link, tinyWorstK);
  check(tinyTaken && tinyFinite && tinyWorstK < 100.0, what);
  frozen.moveTo(1.0);

  checkOrbitDeclines("the collapsed move", data, family, otherRng,
                     [](auto& c) { return collapsedDeclinesWithoutDraw(c); });

  // 6. The coin: an eligible state moves on a quarter of its calls, and a
  // call that does not move consumes exactly the coin.
  size_t numMoved = 0, numDeclinesRead = 0;
  bool oneDrawPerDecline = true;
  for (size_t c = 0; c < collapsedCoinCalls; ++c) {
    ext_rng* reference = numDeclinesRead < 50 ? cloneRng(rng) : nullptr;
    if (TestPeer::maybeDrawCollapsedRescale(chain)) {
      ++numMoved;
      frozen.moveTo(1.0);
    } else if (reference != nullptr) {
      ext_rng_simulateContinuousUniform(reference);
      oneDrawPerDecline = oneDrawPerDecline && rngStreamsAgree(reference, rng);
      ++numDeclinesRead;
    }
    if (reference != nullptr) ext_rng_destroy(reference);
  }
  double calls = static_cast<double>(collapsedCoinCalls);
  double share = static_cast<double>(numMoved) / calls;
  double shareSd = std::sqrt(collapsedRescaleProbability *
                             (1.0 - collapsedRescaleProbability) / calls);
  std::snprintf(what, sizeof what,
                "%s: an eligible sweep takes the move with probability %.2f "
                "(share %.4f, %.2f sd off)", link, collapsedRescaleProbability,
                share, (share - collapsedRescaleProbability) / shareSd);
  check(std::fabs(share - collapsedRescaleProbability) <
          collapsedZBound * shareSd, what);
  check(oneDrawPerDecline && numDeclinesRead == 50,
        "an eligible sweep whose coin says no consumes exactly the coin");

  // 9. Every row out of the likelihood: the move runs, against the target
  // with no likelihood term, k's prior along the orbit.
  RescaleData allOut(data);
  allOut.mask.assign(rescaleN, 0.0);
  check(sampler->setActiveRows(allOut.mask.data()),
        "an all-inactive mask installs");
  double priorDisagreement = 0.0;
  double priorKs =
    collapsedInvarianceKs(frozen, allOut, collapsedDraws, priorDisagreement);
  std::snprintf(what, sizeof what,
                "%s: with every row out the move runs and keeps k's prior "
                "along the orbit (sqrt(N) D %.3g)", link, priorKs);
  check(priorKs < rescaleKsBound && priorDisagreement < 1.0e-10, what);

  // 1 to 4 again on a near-separated fixture with no offset, counts or mask,
  // where the target's right tail is as thick as the prior allows.
  RescaleData separated(data);
  separated.offset.assign(rescaleN, 0.0);
  separated.mask.assign(rescaleN, 1.0);
  separated.counts.assign(rescaleN, 1.0);
  for (size_t i = 0; i < rescaleN; ++i)
    separated.y[i] =
      6.0 * (data.x[i] - 0.5) + 0.3 * standardNormal() > 0.0 ? 1.0 : 0.0;
  std::unique_ptr<ConstantLeafSampler> separatedSampler = makeRescaleSampler(
    separated.x, separated.y, nullptr, nullptr, nullptr, family, true, true,
    true, &otherRng);
  separatedSampler->run(rescaleBurnIn, 0, results);
  std::snprintf(label, sizeof label, "%s, near-separated", link);
  double separatedKs = runCollapsedStateTests(label, *separatedSampler,
                                              separated, worstDisagreement);
  separatedSampler.reset();

  // 1 to 4 again on a binomial response: successes out of 0 to 3 trials, a
  // row of no trials out of the likelihood beside the mask's rows, counts 1
  // and 2, the offset and mask as above.
  RescaleData binomial(data);
  for (size_t i = 0; i < rescaleN; ++i) {
    binomial.trials[i] = static_cast<double>(i % 4);
    binomial.counts[i] = static_cast<double>(1 + i % 2);
    double eta = 1.5 * std::sin(3.0 * data.x[i]) + data.x[i + rescaleN] - 0.5;
    double p = 1.0 / (1.0 + std::exp(-eta));
    binomial.y[i] = 0.0;
    for (size_t j = 0; j < i % 4; ++j)
      binomial.y[i] += runif01() < p ? 1.0 : 0.0;
  }
  std::unique_ptr<ConstantLeafSampler> binomialSampler = makeRescaleSampler(
    binomial.x, binomial.y, binomial.counts.data(), binomial.trials.data(),
    binomial.offset.data(), family, true, true, true, &otherRng);
  binomialSampler->run(rescaleBurnIn, 0, results);
  check(binomialSampler->setActiveRows(binomial.mask.data()),
        "the binomial fixture takes a mask");
  binomialSampler->run(20, 0, results);
  std::snprintf(label, sizeof label, "%s, binomial", link);
  checkLabeled(TestPeer::collapsedRescaleApplies(binomialSampler->chain(0)),
               label, "a burned-in drawn-k forest takes the move");
  double binomialKs = runCollapsedStateTests(label, *binomialSampler, binomial,
                                             worstDisagreement);
  binomialSampler.reset();
  ext_rng_destroy(otherRng);

  // 8. The sweep's opening: a sweep takes the move on its schedule ahead of
  // the augmentation draw, which draws the latents at the scaled fit, and no
  // later step of the sweep moves them. So a twin that takes only those two
  // holds the latents a sweep leaves bit for bit, on a sweep that moves and
  // on one whose coin says no; twins on further seeds until both are seen.
  size_t numOpeningsMoved = 0, numOpeningsStill = 0;
  bool openingsAgree = true;
  for (std::uint32_t seed = 31u; seed < 71u; ++seed) {
    if (numOpeningsMoved > 0 && numOpeningsStill > 0) break;
    ext_rng* sweepRng =
      ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng* twinRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(sweepRng, seed);
    ext_rng_setSeed(twinRng, seed);
    std::unique_ptr<ConstantLeafSampler> swept = makeRescaleSampler(
      data.x, data.y, data.counts.data(), nullptr, data.offset.data(),
      family, true, true, true, &sweepRng);
    std::unique_ptr<ConstantLeafSampler> twin = makeRescaleSampler(
      data.x, data.y, data.counts.data(), nullptr, data.offset.data(),
      family, true, true, true, &twinRng);
    swept->run(20, 0, results);
    twin->run(20, 0, results);
    swept->run(1, 0, results);
    bool moved = TestPeer::maybeDrawCollapsedRescale(twin->chain(0));
    TestPeer::drawAugmentation(twin->chain(0));
    const double* sweptLatents = TestPeer::latents(swept->chain(0));
    openingsAgree = openingsAgree &&
      std::equal(sweptLatents, sweptLatents + rescaleN,
                 TestPeer::latents(twin->chain(0)));
    ++(moved ? numOpeningsMoved : numOpeningsStill);
    swept.reset();
    twin.reset();
    ext_rng_destroy(sweepRng);
    ext_rng_destroy(twinRng);
  }
  std::snprintf(what, sizeof what,
                "%s: a sweep opens with the move on its schedule and then the "
                "augmentation draw (%zu openings that moved, %zu that did not)",
                link, numOpeningsMoved, numOpeningsStill);
  check(openingsAgree && numOpeningsMoved > 0 && numOpeningsStill > 0, what);

  // 7. Every chain.
  checkActsOnEveryChain("the collapsed move", 21u, 50,
                        [&](bool on, ext_rng** rngs) {
    return makeRescaleSampler(data.x, data.y, data.counts.data(), nullptr,
                              data.offset.data(), family, true, true, on,
                              rngs, 2);
  });

  sampler.reset();
  ext_rng_destroy(rng);
  printf("ok: %s, %zu moves at %zu trees (sqrt(N) D %.3g, near-separated "
         "%.3g, binomial %.3g, every row out %.3g; worst factor disagreement "
         "%.3g; coin share %.4f)\n",
         link, collapsedDraws, rescaleTrees, ks, separatedKs, binomialKs,
         priorKs, std::max(worstDisagreement, priorDisagreement), share);
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
  // The shipped default, stated not inherited.
  options.numTrees = ensembleTrees;
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

  // The same two identities with the level-fiber step on. What this arm
  // checks is narrow and worth saying: the assertion fires at the END of a
  // sweep, by which point every leaf has been redrawn, so all that reaches it
  // is the projection's own arithmetic residual sum_t c_t riding in through
  // totalFits - of order m eps tau, and the band above sits nine orders over
  // that. It sizes the missing-projection poison, which would put a whole
  // leaf value there, and says nothing about whether the draw's DISTRIBUTION
  // is right; the frozen-forest gate below is what scores that.
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
  check(anyDifference, "the level-fiber flag moves the sampled path");
  ext_rng_destroy(offRng);

  // ---- the shift's own distribution, at a frozen forest ----------------
  //
  // The step is an exact Gibbs draw, so there is no acceptance ratio for a
  // detailed-balance script to score; what replaces one is a direct test of
  // the conditional. Freeze a burned-in forest, restore its leaf tables
  // between calls so every draw scores against the SAME closed form, and read
  // the empirical mean and covariance of c.
  ext_rng* distRng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(distRng, 3141592u);

  std::vector<double> levelX(levelN * levelP), levelY(levelN);
  for (double& v : levelX) v = runif01();
  for (size_t i = 0; i < levelN; ++i)
    levelY[i] = 2.0 * levelX[i] - levelX[i + levelN] * levelX[i + 2 * levelN] +
                0.3 * (runif01() - 0.5);

  SamplerOptions distOptions;
  distOptions.numTrees = levelTrees;
  distOptions.levelGibbs = LevelGibbsMode::off;
  ConstantLeafSampler distSampler(levelX.data(), levelY.data(), levelN, levelP,
                                  nullptr, nullptr, ResponseFamily::gaussian,
                                  1.0, 3.0, 0.37804942330213542, distOptions,
                                  &distRng);
  distSampler.run(levelBurnIn, 0, results);

  std::vector<std::vector<double>> frozenMu(levelTrees);
  for (size_t t = 0; t < levelTrees; ++t)
    frozenMu[t] = TestPeer::muByTree(distSampler.chain(0), t);
  LevelDistribution dist = levelDistributionOf(distSampler, levelTrees);

  LevelMoments moments(levelTrees);
  std::vector<double> shift(levelTrees, 0.0);
  double worstSum = 0.0;
  bool everyDrawEligible = true;
  for (size_t draw = 0; draw < levelDraws; ++draw) {
    for (size_t t = 0; t < levelTrees; ++t)
      TestPeer::muByTree(distSampler.chain(0), t) = frozenMu[t];
    everyDrawEligible = everyDrawEligible &&
      TestPeer::drawLevelShift(distSampler.chain(0), shift.data());
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t) total += shift[t];
    worstSum = std::max(worstSum, std::fabs(total));
    moments.add(shift);
  }
  for (size_t t = 0; t < levelTrees; ++t)
    TestPeer::muByTree(distSampler.chain(0), t) = frozenMu[t];

  check(everyDrawEligible, "the frozen forest is eligible on every draw");
  char what[160];
  std::snprintf(what, sizeof what,
                "the shift sums to zero on every draw (worst %.3g)", worstSum);
  check(worstSum < levelSumTolerance, what);

  double worstMeanZ = 0.0, worstCovZ = 0.0;
  moments.score(dist, worstMeanZ, worstCovZ);
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
      total += dist.mean[t] + std::sqrt(dist.variance[t]) * standardNormal();
    poisonWorstSum = std::max(poisonWorstSum, std::fabs(total));
  }
  std::snprintf(what, sizeof what,
                "poison (i), no projection: the sum moves (%.3g against %.3g)",
                poisonWorstSum, worstSum);
  check(poisonWorstSum > 1.0e6 * levelSumTolerance, what);

  // Poison (ii), the prior-mean term dropped: u centered at zero and
  // projected. The distribution keeps its covariance and loses its mean, which
  // is the poison that separates a correct conditional from a correctly-shaped
  // one - so the covariance score must still PASS while the mean score fails.
  LevelMoments poisoned(levelTrees);
  double sumVariance = 0.0;
  for (size_t t = 0; t < levelTrees; ++t) sumVariance += dist.variance[t];
  for (size_t draw = 0; draw < levelDraws; ++draw) {
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t) {
      shift[t] = std::sqrt(dist.variance[t]) * standardNormal();
      total += shift[t];
    }
    double ratio = total / sumVariance;
    for (size_t t = 0; t < levelTrees; ++t)
      shift[t] -= dist.variance[t] * ratio;
    poisoned.add(shift);
  }
  double poisonMeanZ = 0.0, poisonCovZ = 0.0;
  poisoned.score(dist, poisonMeanZ, poisonCovZ);
  std::snprintf(what, sizeof what,
                "poison (ii), no prior mean: the mean fails (z %.3g)",
                poisonMeanZ);
  check(poisonMeanZ > levelZThreshold, what);
  std::snprintf(what, sizeof what,
                "poison (ii), no prior mean: the covariance still passes "
                "(z %.3g)", poisonCovZ);
  check(poisonCovZ < levelZThreshold, what);

  ext_rng_destroy(distRng);

  runForestRescaleTests();
  runCollapsedRescaleTests();

  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: ensemble sum invariance, %zu trees x %zu sweeps (worst fit %.3g, "
         "worst residual %.3g)\n", ensembleTrees, ensembleSweeps, worstFit,
         worstResidual);
  printf("ok: level fiber, flag on (worst fit %.3g, worst residual %.3g, "
         "from a prior draw %.3g); "
         "%zu draws at %zu trees (worst |sum| %.3g, mean z %.3g, cov z %.3g); "
         "poisons (no projection |sum| %.3g; no prior mean, z %.3g mean and "
         "%.3g cov)\n",
         shiftedWorstFit, shiftedWorstResidual,
         std::max(priorWorstFit, priorWorstResidual), levelDraws, levelTrees,
         worstSum, worstMeanZ, worstCovZ, poisonWorstSum, poisonMeanZ,
         poisonCovZ);
}
