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

  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: ensemble sum invariance, %zu trees x %zu sweeps (worst fit %.3g, "
         "worst residual %.3g)\n", ensembleTrees, ensembleSweeps, worstFit,
         worstResidual);
}
