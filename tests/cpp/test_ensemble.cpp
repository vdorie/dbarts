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
  const std::vector<double>& total = sampler.chain(0).totalFits();
  std::vector<double> fits = sampler.chain(0).treeFits();
  const double* y = sampler.chain(0).workingResponseForTesting();
  const std::vector<double>& resid = sampler.chain(0).residualForTesting();

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
    const std::vector<double>& mu = sampler.chain(0).muByTreeForTesting(t);
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
  for (double v : sampler.chain(0).totalFits()) {
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
  bool anyDifference = false;
  for (size_t i = 0; i < ensembleN; ++i)
    anyDifference = anyDifference ||
      off.chain(0).totalFits()[i] != shifted.chain(0).totalFits()[i];
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
    frozenMu[t] = lawSampler.chain(0).muByTreeForTesting(t);
  LevelLaw law = levelLawOf(lawSampler, levelTrees);

  LevelMoments moments(levelTrees);
  std::vector<double> shift(levelTrees, 0.0);
  double worstSum = 0.0;
  bool everyDrawEligible = true;
  for (size_t draw = 0; draw < levelDraws; ++draw) {
    for (size_t t = 0; t < levelTrees; ++t)
      lawSampler.chain(0).muByTreeForTesting(t) = frozenMu[t];
    everyDrawEligible = everyDrawEligible &&
      lawSampler.chain(0).drawLevelShiftForTesting(shift.data());
    double total = 0.0;
    for (size_t t = 0; t < levelTrees; ++t) total += shift[t];
    worstSum = std::max(worstSum, std::fabs(total));
    moments.add(shift);
  }
  for (size_t t = 0; t < levelTrees; ++t)
    lawSampler.chain(0).muByTreeForTesting(t) = frozenMu[t];

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
