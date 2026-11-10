#ifndef TESTS_CPP_TEST_PEER_HPP
#define TESTS_CPP_TEST_PEER_HPP

// The one definition of bartcore::TestPeer, the friend each engine class
// declares for the state no production path reads. Only the component tests
// include it, so production code cannot name these hooks.

#include <algorithm>
#include <cmath>
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
  /// Whether every mean tree's leaf table, the constant leaf's or the vector
  /// leaf's blocks, covers every node the tree holds, which a state read
  /// from the chain flattens.
  template <IntegrableLeafModel L, typename R>
  static bool leafTablesCoverTrees(const Chain<L, R>& chain) {
    for (const auto& forest : chain.forests_)
      for (std::size_t t = 0; t < forest.numTrees; ++t) {
        std::size_t nodes = forest.trees[t].nodes.size();
        if constexpr (Chain<L, R>::leafIsConstant) {
          if (forest.muByTree[t].size() < nodes) return false;
        } else if constexpr (L::hasVectorParams) {
          if (forest.paramsByTree[t].size() < nodes * forest.leaf.numParams())
            return false;
        }
      }
    return true;
  }
  /// Whether every tree, mean and variance, holds a row in every leaf, and
  /// holds in each leaf the rows the store's codes send there; whether each
  /// mean tree's map (where it is not marked stale) and a vector leaf's fits
  /// are those its leaves give; and whether each variance tree's factor is
  /// one value over each leaf and the variance surface their product. The
  /// chain against its own store, for a step that moved the store.
  template <IntegrableLeafModel L, typename R>
  static bool treesRouteTheirRows(const Chain<L, R>& chain) {
    const ColumnStore& data = chain.data_;
    std::size_t n = data.numObservations;
    auto near = [](double a, double b) {
      return std::fabs(a - b) <= 1e-9 * (1.0 + std::fabs(b));
    };
    auto routes = [&](const Tree& tree, auto&& rowHolds) {
      std::vector<std::int32_t> bottoms;
      tree.fillBottom(0, bottoms);
      std::size_t count = 0;
      for (std::int32_t b : bottoms) {
        const Node& node = tree.at(b);
        if (node.numObservations() == 0) return false;
        for (std::size_t m = node.begin; m < node.end; ++m) {
          std::size_t obs = tree.indices[m];
          if (tree.findBottomNodeForObservation(data, obs) != b ||
              !rowHolds(b, obs, node))
            return false;
          ++count;
        }
      }
      return count == n;
    };
    for (const auto& forest : chain.forests_)
      for (std::size_t t = 0; t < forest.numTrees; ++t) {
        bool holds = routes(forest.trees[t], [&](std::int32_t b,
                                                 std::size_t obs,
                                                 const Node&) {
          if constexpr (Chain<L, R>::leafIsConstant) {
            return forest.leafOfStale[t] != 0 ||
                   forest.leafOf[t * n + obs] == static_cast<std::uint32_t>(b);
          } else if constexpr (L::hasVectorParams) {
            const double* params = forest.paramsByTree[t].data() +
                                   static_cast<std::size_t>(b) *
                                     forest.leaf.numParams();
            return near(forest.treeFits[t * n + obs],
                        forest.leaf.fitForObservation(params, obs));
          } else {
            (void) b;
            (void) obs;
            return true;
          }
        });
        if (!holds) return false;
      }
    if (chain.varianceForest_) {
      const auto& vf = *chain.varianceForest_;
      for (std::size_t j = 0; j < vf.numTrees; ++j) {
        const double* factor = vf.factorByTree.data() + j * n;
        const Tree& tree = vf.trees[j];
        bool holds = routes(tree, [&](std::int32_t, std::size_t obs,
                                      const Node& node) {
          return factor[obs] == factor[tree.indices[node.begin]];
        });
        if (!holds) return false;
      }
      for (std::size_t i = 0; i < n; ++i) {
        double product = 1.0;
        for (std::size_t j = 0; j < vf.numTrees; ++j)
          product *= vf.factorByTree[j * n + i];
        if (!near(vf.combinedVariance[i], product)) return false;
      }
    }
    return true;
  }
  /// Whether each forest's leaf covariates are the store's raw values under
  /// the centers and scales the leaf holds; true for a leaf that takes none.
  template <IntegrableLeafModel L, typename R>
  static bool leafCovariatesAreStandardized(const Chain<L, R>& chain) {
    if constexpr (L::hasVectorParams || L::hasFunctionParams) {
      const ColumnStore& data = chain.data_;
      std::size_t n = data.numObservations;
      for (const auto& forest : chain.forests_) {
        const auto& leaf = forest.leaf;
        std::size_t p = leaf.numCovariates_;
        if (leaf.numObservations_ != n) return false;
        for (std::size_t j = 0; j < p; ++j) {
          const double* raw = data.rawColumn(leaf.columns_[j]);
          if (raw == nullptr) return false;
          for (std::size_t i = 0; i < n; ++i) {
            double expected = std::isnan(raw[i])
              ? 0.0
              : (raw[i] - leaf.means_[j]) / leaf.sds_[j];
            double held = L::hasVectorParams ? leaf.u_[i * p + j]
                                             : leaf.u_[j * n + i];
            if (!(std::fabs(held - expected) <=
                  1e-12 * (1.0 + std::fabs(expected))))
              return false;
          }
        }
      }
    }
    return true;
  }
  /// Whether the chain's test fits, its variance surface and its leaf's
  /// test covariates are sized for the store's test rows.
  template <IntegrableLeafModel L, typename R>
  static bool testStorageCoversTheStore(const Chain<L, R>& chain) {
    std::size_t m = chain.data_.numTestObservations;
    for (const auto& forest : chain.forests_) {
      if (forest.totalTestFits.size() < m || forest.currTestFits.size() < m)
        return false;
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        if (m > 0 && (forest.leaf.numTestObservations_ != m ||
                      forest.leaf.uTest_.size() !=
                        m * forest.leaf.numCovariates_))
          return false;
    }
    return !chain.varianceForest_ ||
           chain.varianceForest_->combinedVarianceTest.size() >= m;
  }
  /// Whether every tree, mean and variance, holds a row in every leaf.
  template <IntegrableLeafModel L, typename R>
  static bool everyLeafHoldsARow(const Chain<L, R>& chain) {
    for (const auto& forest : chain.forests_)
      for (std::size_t t = 0; t < forest.numTrees; ++t)
        if (!forest.trees[t].bottomNodesAreOccupied()) return false;
    if (chain.varianceForest_)
      for (const Tree& tree : chain.varianceForest_->trees)
        if (!tree.bottomNodesAreOccupied()) return false;
    return true;
  }
  /// Whether every monotone tree's leaf values lie in its order; true off the
  /// monotone leaf.
  template <IntegrableLeafModel L, typename R>
  static bool monotoneLeavesInOrder(const Chain<L, R>& chain) {
    return chain.monotoneLeavesInOrder();
  }
  /// Settles the chain with forest f's tree t in flight, as an exception out
  /// of a move on that tree leaves it.
  template <IntegrableLeafModel L, typename R>
  static void settleTreeInFlight(Chain<L, R>& chain, std::size_t f,
                                 std::size_t t) {
    chain.forestInFlight_ = f;
    chain.treeInFlight_ = t;
    chain.settleAfterAllocationFailure();
  }
  /// Forest f's cached total fits, over a bare chain; a site holding a sampler
  /// reads SamplerBase::forestTotalFits.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<double>& totalFitsInForest(const Chain<L, R>& chain,
                                                      std::size_t f) {
    return chain.forests_[f].totalFits;
  }
  /// Forest f's leaf scale, the per-tree prior sd on the working scale.
  template <IntegrableLeafModel L, typename R>
  static double leafScale(const Chain<L, R>& chain, std::size_t f) {
    return chain.forests_[f].leaf.scale;
  }
  /// Mean forest f's tree prior in force.
  template <IntegrableLeafModel L, typename R>
  static const CGMTreePrior& forestTreePrior(const Chain<L, R>& chain,
                                             std::size_t f) {
    return chain.forests_[f].treePrior;
  }
  /// Map forest f's entry, the factor and divisor of its node scale.
  template <IntegrableLeafModel L, typename R>
  static const ForestMapEntry& mapEntry(const Chain<L, R>& chain,
                                        std::size_t f) {
    return chain.forestMap_[f].entry;
  }
  /// Map forest f's statement, and whether it keeps the map's held entry.
  template <IntegrableLeafModel L, typename R>
  static const ForestStatement& mapStatement(const Chain<L, R>& chain,
                                             std::size_t f) {
    return chain.forestMap_[f].statement;
  }
  template <IntegrableLeafModel L, typename R>
  static bool heldAtDefaultK(const Chain<L, R>& chain, std::size_t f) {
    return chain.forestMap_[f].heldAtDefaultK;
  }
  /// Forest f's coefficient prior as the combiner holds it, its variance the
  /// live mixing variance where it mixes, and its block.
  template <IntegrableLeafModel L, typename R>
  static const ForestAmplitudePrior& amplitudePrior(const Chain<L, R>& chain,
                                                    std::size_t f) {
    const auto* combiner =
      dynamic_cast<const AmplitudeForestCombiner<L, R>*>(chain.combiner_.get());
    return combiner->glue_.prior[f];
  }
  template <IntegrableLeafModel L, typename R>
  static std::vector<double> amplitudeBlock(const Chain<L, R>& chain,
                                            std::size_t f) {
    const auto* combiner =
      dynamic_cast<const AmplitudeForestCombiner<L, R>*>(chain.combiner_.get());
    const double* block = combiner->glue_.amplitudesOf(f);
    return std::vector<double>(block, block + combiner->glue_.numAmplitudes(f));
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
  /// The variance forest, for the tree-move probabilities, perturb width and
  /// tree prior in force on it. Dereferences the forest.
  template <IntegrableLeafModel L, typename R>
  static const VarianceForest& varianceForest(const Chain<L, R>& chain) {
    return *chain.varianceForest_;
  }
  /// Mean forest f's move settings in force, read off its fields.
  template <IntegrableLeafModel L, typename R>
  static MoveSettings forestMoveSettings(const Chain<L, R>& chain,
                                         std::size_t f) {
    const Forest<L, R>& forest = chain.forests_[f];
    MoveSettings s;
    s.birthOrDeathProbability = forest.birthOrDeathProbability;
    s.swapProbability = forest.swapProbability;
    s.changeProbability = forest.changeProbability;
    s.perturbProbability = forest.perturbProbability;
    s.perturbWidth = forest.perturbWidth;
    s.ruleGibbsProbability = forest.ruleGibbsProbability;
    s.birthProbability = forest.birthProbability;
    s.levelGibbs = forest.levelGibbs;
    return s;
  }
  /// The variance forest's move settings in force; its mode is unread and
  /// left at the default.
  template <IntegrableLeafModel L, typename R>
  static MoveSettings varianceMoveSettings(const Chain<L, R>& chain) {
    const VarianceForest& vf = *chain.varianceForest_;
    MoveSettings s;
    s.birthOrDeathProbability = vf.birthOrDeathProbability;
    s.swapProbability = vf.swapProbability;
    s.changeProbability = vf.changeProbability;
    s.perturbProbability = vf.perturbProbability;
    s.perturbWidth = vf.perturbWidth;
    s.ruleGibbsProbability = vf.ruleGibbsProbability;
    s.birthProbability = vf.birthProbability;
    return s;
  }
  template <IntegrableLeafModel L, typename R>
  static FunctionLeafDrawStats accountStrandedLeafKStats(
    Chain<L, R>& chain, int32_t variableIndex, int32_t splitIndex) {
    return chain.accountStrandedLeafKStats(variableIndex, splitIndex);
  }
  /// Forest 0's level-fiber shift, drawn once against the leaf tables as they
  /// stand and applied to them. shiftOut receives one c_t per tree, zero where
  /// the tree declined; the return says whether the forest as a whole was
  /// eligible.
  template <IntegrableLeafModel L, typename R>
  static bool drawLevelShift(Chain<L, R>& chain, double* shiftOut) {
    return chain.drawLevelShift(chain.forests_[0], shiftOut);
  }
  /// The same step on forest f.
  template <IntegrableLeafModel L, typename R>
  static bool drawForestLevelShift(Chain<L, R>& chain, std::size_t f) {
    return chain.drawLevelShift(chain.forests_[f]);
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
  /// The collapsed rescaling move on the orbit's forests without its coin,
  /// the move on its schedule, and the predicate; factorOut receives the
  /// factor applied, 1 on a decline.
  template <IntegrableLeafModel L, typename R>
  static bool drawCollapsedRescale(Chain<L, R>& chain, double* factorOut) {
    return chain.drawCollapsedRescale(factorOut);
  }
  template <IntegrableLeafModel L, typename R>
  static bool maybeDrawCollapsedRescale(Chain<L, R>& chain) {
    return chain.maybeDrawCollapsedRescale();
  }
  template <IntegrableLeafModel L, typename R>
  static bool collapsedRescaleApplies(const Chain<L, R>& chain) {
    return chain.collapsedRescaleApplies();
  }
  /// The chain's copy of the engine toggle for that move, writable.
  template <IntegrableLeafModel L, typename R>
  static bool& collapsedRescaleForest(Chain<L, R>& chain) {
    return chain.options_.collapsedRescaleForest;
  }
  /// The augmentation draw a sweep takes after the collapsed move, against
  /// the fits as they stand, with no cancel armed.
  template <IntegrableLeafModel L, typename R>
  static void drawAugmentation(Chain<L, R>& chain) {
    chain.drawAugmentation(false);
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
  /// The family's response model, and the first forest's DART prior.
  template <IntegrableLeafModel L, typename R>
  static const ResponseModel& response(const Chain<L, R>& chain) {
    return *chain.response_;
  }
  template <IntegrableLeafModel L, typename R>
  static DartPrior& dart(Chain<L, R>& chain) {
    return chain.forests_[0].dart;
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
  /// The family's augmentation draw made by hand on the chain's generator,
  /// against its combined fits and working sigma as they stand.
  template <IntegrableLeafModel L, typename R>
  static void refreshLatents(Chain<L, R>& chain) {
    chain.response_->refreshLatents(chain.rng_, chain.combinedFits(),
                                    chain.sigma_);
  }
  /// Every forest's linear-leaf crossproduct cache given no budget, so each
  /// node statistic is a fresh scan.
  template <typename R>
  static void disableStatisticsCache(Chain<LinearGaussianLeaf, R>& chain) {
    for (Forest<LinearGaussianLeaf, R>& forest : chain.forests_) {
      forest.leaf.clearStatisticsCache();
      forest.leaf.statisticsCacheBudget_ = 0;
    }
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
  /// Forest f's tree t leaf table.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<double>& forestMuByTree(const Chain<L, R>& chain,
                                                   std::size_t f,
                                                   std::size_t t) {
    return chain.forests_[f].muByTree[t];
  }
  /// Forest f's tree t parameter blocks (a vector-parameter leaf), one block
  /// of the leaf's numParams per node.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<double>& forestParamsByTree(const Chain<L, R>& chain,
                                                       std::size_t f,
                                                       std::size_t t) {
    return chain.forests_[f].paramsByTree[t];
  }
  /// Forest f's k, whether it is drawn, and its k prior, writable.
  template <IntegrableLeafModel L, typename R>
  static double& forestKAt(Chain<L, R>& chain, std::size_t f) {
    return chain.forests_[f].k;
  }
  template <IntegrableLeafModel L, typename R>
  static bool forestUpdateKAt(const Chain<L, R>& chain, std::size_t f) {
    return chain.forests_[f].updateK;
  }
  template <IntegrableLeafModel L, typename R>
  static const ChiKHyperprior& kHyperpriorAt(const Chain<L, R>& chain,
                                             std::size_t f) {
    return chain.forests_[f].kHyperprior;
  }
  /// Forest f's obs-to-leaf map for tree t.
  template <IntegrableLeafModel L, typename R>
  static const std::uint32_t* forestLeafOf(const Chain<L, R>& chain,
                                           std::size_t f, std::size_t t) {
    return chain.forests_[f].leafOf.data() + t * chain.data_.numObservations;
  }
  /// The forests sharing one total, and the group's k draw against the
  /// statistics the last sweep left.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<std::size_t>& sharers(const Chain<L, R>& chain) {
    return chain.sharers_;
  }
  template <IntegrableLeafModel L, typename R>
  static void drawSharedK(Chain<L, R>& chain) {
    chain.drawSharedK();
  }
  /// Forest f's blocks() group factors on its leaf scale, empty where none is
  /// stored.
  template <IntegrableLeafModel L, typename R>
  static const std::vector<double>& blockScaleFactor(const Chain<L, R>& chain,
                                                     std::size_t f) {
    return chain.forests_[f].blockScaleFactor;
  }
  /// Forest f's k statistic as the last sweep left it: the sum of squares and
  /// the leaf count the k draw reads.
  template <IntegrableLeafModel L, typename R>
  static std::pair<double, double> kStatistic(const Chain<L, R>& chain,
                                              std::size_t f) {
    return {chain.forests_[f].kSumSquaredParams, chain.forests_[f].kNumLeaves};
  }
  /// The combiner's after-sweep step (the multinomial grand-level shift),
  /// drawn from \p rng against the leaf tables as they stand.
  template <IntegrableLeafModel L, typename R>
  static void afterCombine(Chain<L, R>& chain, ext_rng* rng) {
    chain.combiner_->afterCombine(chain.forests_, false, 0, rng);
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
  /// States forest f of a combining chain by hand before its first sweep: its
  /// leaf scale, its coefficient's prior and its block, as no map resolves them.
  template <IntegrableLeafModel L, typename R>
  static void setForestUnits(Chain<L, R>& chain, std::size_t f, double leafScale,
                             const ForestAmplitudePrior& prior,
                             const std::vector<double>& block) {
    chain.forests_[f].leaf.scale = leafScale;
    auto* combiner =
      dynamic_cast<AmplitudeForestCombiner<L, R>*>(chain.combiner_.get());
    combiner->glue_.prior[f] = prior;
    double* amplitudes = combiner->glue_.amplitudesOf(f);
    for (std::size_t j = 0; j < block.size(); ++j) amplitudes[j] = block[j];
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
  /// A gp leaf's member sort: the dispatch, each method on its own, the rule
  /// that picks between them and the bitmap (empty where the rule refuses
  /// every gp leaf of the data), writable so a test can leave a member's bit
  /// set and read which method ran: only the bitmap's walk clears it.
  static const index_t* sortedMembers(const GPGaussianLeaf& leaf,
                                      const Tree& tree, std::int32_t nodeIndex) {
    return leaf.sortedMembers(tree, tree.at(nodeIndex));
  }
  static void sortMembersByBitmap(const GPGaussianLeaf& leaf,
                                  const index_t* members, std::size_t numObs,
                                  index_t* out) {
    leaf.sortMembersByBitmap(members, numObs, out);
  }
  static void sortMembersByComparison(const index_t* members,
                                      std::size_t numObs, index_t* out) {
    GPGaussianLeaf::sortMembersByComparison(members, numObs, out);
  }
  static bool sortsByBitmap(std::size_t numObs, std::size_t numWords) {
    return GPGaussianLeaf::sortsByBitmap(numObs, numWords);
  }
  static std::vector<std::uint64_t>& memberBits(const GPGaussianLeaf& leaf) {
    return leaf.memberBits_;
  }
  /// Whether the kernel cache holds an entry for the tree's node. Read it
  /// before the node's next score or draw, whose lookup rebuilds a dropped
  /// entry.
  static bool holdsKernelEntry(const GPGaussianLeaf& leaf, const Tree& tree,
                               std::int32_t nodeIndex) {
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    for (const auto& cache : leaf.kernelCaches_)
      if (cache.tree == &tree)
        return index < cache.nodes.size() &&
               !cache.nodes[index].members.empty();
    return false;
  }
  /// The residual-variance posterior's degrees of freedom, nu_0 + #{w_i > 0}
  /// over the model's own precisions, nu_0 + sum of the counts over the rows
  /// in under frequency; AFT's is its contained Gaussian's.
  static double sigmaDegreesOfFreedom(const GaussianResponse& response) {
    return response.sigmaSqPrior_.degreesOfFreedom + response.rowCount_;
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

  // Sampler

  /// Sampler::installAcrossGrids with an install step of the test's own.
  template <typename S, typename Install>
  static bool installAcrossGrids(
      S& sampler, const std::vector<std::vector<double>>& storedCutPoints,
      const double* currentPredictors, Install&& install) {
    return sampler.installAcrossGrids(storedCutPoints, currentPredictors,
                                      std::forward<Install>(install));
  }
};

}  // namespace bartcore

#endif  // TESTS_CPP_TEST_PEER_HPP
