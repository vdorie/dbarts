#ifndef TESTS_CPP_TEST_PEER_HPP
#define TESTS_CPP_TEST_PEER_HPP

// The one definition of bartcore::TestPeer, the friend each engine class
// declares for the state no production path reads. Only the component tests
// include it, so production code cannot name these hooks.

#include <cstddef>
#include <cstdint>
#include <vector>

#include <external/random.h>

#include <bartcore/chain.hpp>
#include <bartcore/model.hpp>

namespace bartcore {

struct TestPeer {
  // Chain

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
  /// The surface the response model holds, which must be varianceFits() by
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
  /// whose product over j is the combined variance varianceFits() reports.
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

  // Leaf and response models

  static std::size_t statisticsCacheResidentBytes(
    const LinearGaussianLeaf& leaf) {
    return leaf.statisticsCacheResidentBytes();
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
  /// The precomputed dispersion kernel L_k, and the one a negative binomial
  /// response has installed, with the collapsed statistic S its grid draw
  /// reads.
  static double kernelValue(const NBDispersionPrior& prior, std::size_t k) {
    return prior.kernel_[k];
  }
  static double dispersionKernel(const NBResponse& response, std::size_t k) {
    return kernelValue(response.rPrior_, k);
  }
  static double collapsedStatistic(const NBResponse& response,
                                   const double* totalFits) {
    return response.collapsedStatistic(totalFits);
  }
};

}  // namespace bartcore

#endif  // TESTS_CPP_TEST_PEER_HPP
