#ifndef TESTS_CPP_TEST_PEER_HPP
#define TESTS_CPP_TEST_PEER_HPP

// The one definition of bartcore::TestPeer, the friend each engine class
// declares for the state no production path reads. Only the component tests
// include it, so production code cannot name these hooks.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <external/random.h>

#include <bartcore/chain.hpp>
#include <bartcore/model.hpp>

namespace bartcore {

struct TestPeer {
  // Chain

  /// The working response the chain's trees are fitted to.
  template <IntegrableLeafModel L, typename R>
  static const double* workingResponse(Chain<L, R>& chain) {
    return chain.response_->workingResponse();
  }

  /// The dense per-tree fit slab of forest 0, materialized (the constant leaf
  /// gathers its compact tables into the returned buffer).
  template <IntegrableLeafModel L, typename R>
  static std::vector<double> treeFits(const Chain<L, R>& chain) {
    std::vector<double> out(chain.data_.numObservations *
                            chain.forests_[0].numTrees);
    chain.forestTreeFits(0, out.data());
    return out;
  }
  /// Forest f's per-tree fit slabs, tree-major (numObservations x numTrees).
  template <IntegrableLeafModel L, typename R>
  static void forestTreeFits(const Chain<L, R>& chain, std::size_t f,
                             double* out) {
    chain.forestTreeFits(f, out);
  }
  /// Forest f's cached total fits, over a bare chain; a site holding a sampler
  /// reads SamplerBase::forestTotalFits.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<double>& totalFitsInForest(const Chain<L, R>& chain,
                                                      std::size_t f) {
    return chain.forests_[f].totalFits;
  }
  /// Tree t of forest f, and the location the chain hands a latent redraw.
  template <IntegrableLeafModel L, typename R>
  static const Tree& forestTree(const Chain<L, R>& chain, std::size_t f,
                                std::size_t t) {
    return chain.forests_[f].trees[t];
  }
  template <IntegrableLeafModel L, typename R>
  static const double* combinedFits(Chain<L, R>& chain) {
    return chain.combinedFits();
  }
  /// The per-observation precisions the chain's trees and sigma draw read.
  template <IntegrableLeafModel L, typename R>
  static const double* workingWeights(const Chain<L, R>& chain) {
    return chain.response_->workingWeights();
  }
  /// The residual scale on the working scale, which is the one the chain
  /// hands a latent redraw.
  template <IntegrableLeafModel L, typename R>
  static double workingSigma(const Chain<L, R>& chain) {
    return chain.sigma_;
  }
  /// The per-observation weight installed on forest f, or null when none is.
  template <IntegrableLeafModel L, typename R>
  static const double* forestWeights(const Chain<L, R>& chain, std::size_t f) {
    return chain.forestWeights_.empty() ? nullptr : chain.forestWeights_[f];
  }
  /// Tree t's obs-to-leaf map (constant leaf, forest 0), where entry i is the
  /// arena bottom-node index owning observation i.
  template <IntegrableLeafModel L, typename R>
  static const std::uint32_t* leafOf(const Chain<L, R>& chain, std::size_t t) {
    return chain.forests_[0].leafOf.data() + t * chain.data_.numObservations;
  }
  /// Tree t's rebuild mark, the fused suffstat pass's eligibility gate (forest
  /// 0). Non-zero means the map still describes the previous partition, so the
  /// sweep rebuilds before the tree's draw.
  template <IntegrableLeafModel L, typename R>
  static std::uint8_t leafOfStale(const Chain<L, R>& chain, std::size_t t) {
    return chain.forests_[0].leafOfStale[t];
  }
  /// The running residual forest f rolls across a sweep, and the working
  /// response it is rolled against. Together with treeFits they pin the
  /// unrolled mu[leafOf] gathers elementwise, tail included.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<R>& residual(const Chain<L, R>& chain,
                                        std::size_t f = 0) {
    return chain.forests_[f].treeY;
  }
  template <IntegrableLeafModel L, typename R>
  static const double* workingResponse(const Chain<L, R>& chain) {
    return chain.response_->workingResponse();
  }
  /// The current combined variance s^2(x_i), working scale, over the training
  /// or test rows, or null when homoscedastic. Original-scale reporting
  /// multiplies by sigmaScale^2.
  template <IntegrableLeafModel L, typename R>
  static const double* varianceFits(const Chain<L, R>& chain) {
    return chain.varianceForest_
             ? chain.varianceForest_->combinedVariance.data()
             : nullptr;
  }
  template <IntegrableLeafModel L, typename R>
  static const double* varianceTestFits(const Chain<L, R>& chain) {
    return chain.varianceForest_
             ? chain.varianceForest_->combinedVarianceTest.data()
             : nullptr;
  }
  /// The surface the response model holds, which must be varianceFits by
  /// pointer identity after every allocation of the combined-variance storage
  /// (installVarianceSurface is what keeps it so). Null for a family that
  /// keeps no surface.
  template <IntegrableLeafModel L, typename R>
  static const double* installedVarianceSurface(const Chain<L, R>& chain) {
    const ResponseModel* response = chain.response_.get();
    if (auto* gaussian = dynamic_cast<const GaussianResponse*>(response))
      return gaussian->variance_;
    if (auto* aft = dynamic_cast<const AFTResponse*>(response))
      return aft->variance_;
    return nullptr;
  }
  /// The sigma posterior's degrees of freedom, nu_0 plus the count of positive
  /// precisions on the RESPONSE model. A per-forest weight lives on the chain
  /// and must never reach them. Zero for a family that draws no sigma.
  template <IntegrableLeafModel L, typename R>
  static double sigmaDegreesOfFreedom(const Chain<L, R>& chain) {
    const ResponseModel* response = chain.response_.get();
    if (auto* gaussian = dynamic_cast<const GaussianResponse*>(response))
      return sigmaDegreesOfFreedom(*gaussian);
    if (auto* aft = dynamic_cast<const AFTResponse*>(response))
      return sigmaDegreesOfFreedom(*aft);
    return 0.0;
  }
  /// The per-tree factor slab h_j(x_i), tree-major (numVarianceTrees x n),
  /// whose product over j is the combined variance varianceFits reports.
  template <IntegrableLeafModel L, typename R>
  static const double* varianceFactors(const Chain<L, R>& chain) {
    return chain.varianceForest_->factorByTree.data();
  }
  /// The scale leaf's calibration (nu', lambda'^2) in force, which is what a
  /// re-anchoring swap must restate.
  template <IntegrableLeafModel L, typename R>
  static const ConstantVarianceLeaf& varianceLeaf(const Chain<L, R>& chain) {
    return chain.varianceForest_->leaf;
  }
  template <IntegrableLeafModel L, typename R>
  static FunctionLeafDrawStats accountStrandedLeafKStats(
    Chain<L, R>& chain, int32_t variableIndex, int32_t splitIndex) {
    return chain.accountStrandedLeafKStats(variableIndex, splitIndex);
  }
  /// Forest 0's level-fibre shift, drawn once against the leaf tables as they
  /// stand and applied to them. shiftOut receives one c_t per tree, zero where
  /// the tree declined; the return says whether the forest as a whole was
  /// eligible.
  template <IntegrableLeafModel L, typename R>
  static bool drawLevelShift(Chain<L, R>& chain, double* shiftOut) {
    return chain.drawLevelShift(chain.forests_[0], shiftOut);
  }
  /// The probit rescaling step on forest 0, and its predicate; alphaOut
  /// receives the factor applied, 1 on a decline.
  template <IntegrableLeafModel L, typename R>
  static bool drawForestRescale(Chain<L, R>& chain, double* alphaOut) {
    return chain.drawForestRescale(alphaOut);
  }
  template <IntegrableLeafModel L, typename R>
  static bool forestRescaleApplies(const Chain<L, R>& chain) {
    return chain.forestRescaleApplies();
  }
  /// Forest 0's k and its hyperprior, writable, so a gate can move a frozen
  /// state along the step's orbit.
  template <IntegrableLeafModel L, typename R>
  static double& forestK(Chain<L, R>& chain) {
    return chain.forests_[0].k;
  }
  template <IntegrableLeafModel L, typename R>
  static bool& forestUpdateK(Chain<L, R>& chain) {
    return chain.forests_[0].updateK;
  }
  template <IntegrableLeafModel L, typename R>
  static ChiKHyperprior& kHyperprior(Chain<L, R>& chain) {
    return chain.forests_[0].kHyperprior;
  }
  /// The family's latents and offset, and a restore of the former that
  /// rebuilds the working response.
  template <IntegrableLeafModel L, typename R>
  static const double* latents(const Chain<L, R>& chain) {
    return chain.response_->latents();
  }
  template <IntegrableLeafModel L, typename R>
  static const double* offset(const Chain<L, R>& chain) {
    return chain.response_->offset();
  }
  template <IntegrableLeafModel L, typename R>
  static void restoreLatents(Chain<L, R>& chain, const double* latents) {
    chain.response_->restoreLatents(latents);
  }
  /// totalFits re-derived from the leaf tables in tree order, for a test that
  /// writes the tables directly.
  template <IntegrableLeafModel L, typename R>
  static void rebuildTotalFits(Chain<L, R>& chain) {
    chain.rebuildTotalFitsFromTrees();
  }
  /// Forest 0's accumulated in-place factor on totalFits, writable, so a test
  /// can put it at or past the bound that forces a re-derivation.
  template <IntegrableLeafModel L, typename R>
  static double& totalFitsScale(Chain<L, R>& chain) {
    return chain.forests_[0].totalFitsScale;
  }
  /// The chain's generator, so a test can read the stream a kernel left.
  template <IntegrableLeafModel L, typename R>
  static ext_rng* rng(Chain<L, R>& chain) {
    return chain.rng_;
  }
  /// Strands an empty leaf in forest 0, which no public mutation does: the
  /// first leaf whose rows all code at 1 or more on ordinal column 0 is split
  /// just below its lowest code, so every row goes right. The right child
  /// keeps the leaf's value, the empty left one takes emptyValue, and the map
  /// and totalFits are rebuilt, so the fits do not move. Returns the empty
  /// leaf as (tree, node), or (numTrees, -1) when no leaf qualifies.
  template <IntegrableLeafModel L, typename R>
  static std::pair<std::size_t, std::int32_t>
  strandEmptyLeaf(Chain<L, R>& chain, double emptyValue) {
    auto& forest = chain.forests_[0];
    std::vector<std::int32_t> bottoms;
    for (std::size_t t = 0; t < forest.numTrees; ++t) {
      Tree& tree = forest.trees[t];
      bottoms.clear();
      tree.fillBottom(0, bottoms);
      for (std::int32_t b : bottoms) {
        const Node& node = tree.at(b);
        if (node.numObservations() == 0) continue;
        xint_t lowest = naCode;
        for (std::size_t m = node.begin; m < node.end; ++m)
          lowest = std::min(lowest, chain.data_.codeAt(0, tree.indices[m]));
        if (lowest == 0 || lowest == naCode) continue;
        Rule rule;
        rule.variableIndex = 0;
        rule.setSplitIndex(static_cast<std::int32_t>(lowest) - 1);
        tree.birth(chain.data_, b, rule, chain.response_->workingResponse(),
                   chain.response_->workingWeights());
        std::int32_t left = tree.at(b).leftChild;
        std::vector<double>& mu = forest.muByTree[t];
        mu.resize(tree.nodes.size(), 0.0);
        mu[static_cast<std::size_t>(left)] = emptyValue;
        mu[static_cast<std::size_t>(left + 1)] = mu[static_cast<std::size_t>(b)];
        chain.rebuildLeafOf(forest, t);
        chain.rebuildTotalFitsFromTrees();
        return {t, left};
      }
    }
    return {forest.numTrees, -1};
  }
  /// Forest 0's tree t leaf table, writable, so a distributional gate can
  /// restore a frozen leaf state between repeated draws.
  template <IntegrableLeafModel L, typename R>
  static std::vector<double>& muByTree(Chain<L, R>& chain, std::size_t t) {
    return chain.forests_[0].muByTree[t];
  }
  /// (tree, sweep) bodies that took the fused roll + suffstat pass since the
  /// chain was built. Monotone, so tests read differences.
  template <IntegrableLeafModel L, typename R>
  static std::size_t fusedSuffstatRuns(const Chain<L, R>& chain) {
    return chain.fusedSuffstatRuns_;
  }
  template <IntegrableLeafModel L, typename R>
  static FusedSuffstatCheck checkFusedSuffstatAgainstStock(Chain<L, R>& chain) {
    return chain.checkFusedSuffstatAgainstStock();
  }

  /// bcf's (a, b0, b1) reading of the amplitude channel, for the conditionals
  /// and component pins written in its spelling; false off a chain or combiner
  /// whose amplitude layout is not bcf's K = 2, q = (1, 2).
  template <IntegrableLeafModel L, typename R>
  static bool bcfGlue(const Chain<L, R>& chain, double& a, double& b0,
                      double& b1) {
    if (chain.totalAmplitudes() != 3 || chain.numForestAmplitudes(0) != 1)
      return false;
    double out[3];
    chain.combiner_->amplitudes(out);
    a = out[0]; b0 = out[1]; b1 = out[2];
    return true;
  }
  template <IntegrableLeafModel L, typename R>
  static bool bcfGlue(const AmplitudeForestCombiner<L, R>& combiner, double& a,
                      double& b0, double& b1) {
    const auto& glue = combiner.glue_;
    if (glue.amplitudes.size() != 3 || glue.numAmplitudes(0) != 1)
      return false;
    a = glue.a(); b0 = glue.b0(); b1 = glue.b1();
    return true;
  }
  /// Writes bcf's (a, b0, b1) whatever the forests' update switches say,
  /// which restoreGlue does not: a pinned amplitude is model there.
  template <IntegrableLeafModel L, typename R>
  static void setBcfGlue(AmplitudeForestCombiner<L, R>& combiner, double a,
                         double b0, double b1) {
    auto& glue = combiner.glue_;
    glue.a() = a; glue.b0() = b0; glue.b1() = b1;
  }

  // Leaf and response models

  static std::size_t statisticsCacheResidentBytes(
    const LinearGaussianLeaf& leaf) {
    return leaf.statisticsCacheResidentBytes();
  }
  /// The bytes the crossproduct cache prices against its budget.
  static std::size_t statisticsCacheUsedBytes(const LinearGaussianLeaf& leaf) {
    return leaf.statisticsCacheUsedBytes_;
  }
  /// The residual-variance posterior's degrees of freedom, nu_0 + #{w_i > 0}
  /// over the model's own precisions; AFT's is its contained Gaussian's.
  static double sigmaDegreesOfFreedom(const GaussianResponse& response) {
    return response.sigmaSqPrior_.degreesOfFreedom +
           static_cast<double>(response.numPositiveWeights_);
  }
  static double sigmaDegreesOfFreedom(const AFTResponse& response) {
    return sigmaDegreesOfFreedom(*response.gaussian_);
  }
  /// The ordinal log Metropolis acceptance for moving free cutpoint gamma_s to
  /// proposal, and the three per-sweep kernels in isolation (refreshLatents
  /// runs two at once); computeScales returns the free cutpoints' proposal
  /// scales, its only observable.
  static double ordinalThresholdLogAcceptance(const OrdinalResponse& response,
                                              const double* totalFits,
                                              std::size_t s, double proposal) {
    return response.ordinalThresholdLogAcceptance(totalFits, s, proposal);
  }
  static const double* computeScales(OrdinalResponse& response) {
    response.computeScales();
    return response.proposalScale_.data();
  }
  static void updateOrdinalThresholds(OrdinalResponse& response, ext_rng* rng,
                                      const double* totalFits) {
    response.updateOrdinalThresholds(rng, totalFits);
  }
  static void drawLatents(OrdinalResponse& response, ext_rng* rng,
                          const double* totalFits) {
    response.drawLatents(rng, totalFits);
  }
  /// The precomputed shape kernel K_k, the one a negative binomial
  /// response has installed, and the grid probabilities the last drawIndex
  /// normalized in place.
  static double kernelValue(const NBShapePrior& prior, std::size_t k) {
    return prior.kernel_[k];
  }
  static double shapeKernel(const NBResponse& response, std::size_t k) {
    return kernelValue(response.rPrior_, k);
  }
  static double drawnProbability(const NBShapePrior& prior,
                                 std::size_t k) {
    return prior.weight_[k];
  }
};

}  // namespace bartcore

#endif  // TESTS_CPP_TEST_PEER_HPP
