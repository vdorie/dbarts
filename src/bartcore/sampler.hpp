#ifndef BARTCORE_SAMPLER_HPP
#define BARTCORE_SAMPLER_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include <external/random.h>
#include <misc/linearAlgebra.h>

#include "data.hpp"
#include "model.hpp"
#include "moves.hpp"
#include "tree.hpp"

namespace bartcore {

struct SamplerOptions {
  size_t numTrees = 200;
  double k = 2.0;
  double nodeScale = 0.5;  // 3.0 for binary responses
  double base = 0.95, power = 2.0;
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  std::uint32_t maxNumCuts = 100;
  // borrowed per-column override of maxNumCuts; copied during construction
  const std::uint32_t* maxNumCutsPerVariable = nullptr;
  bool useQuantiles = false;

  // split-variable selection: fixed weights (borrowed; normalized over
  // available variables at each node) or DART; both null/false = uniform
  const double* splitProbabilities = nullptr;
  bool useDart = false;
  DartPrior dart;

  // k is fixed at .k unless updateK, in which case .k is the starting value
  bool updateK = false;
  ChiKHyperprior kHyperprior;
};

/// Outcome of a transactional predictor change. invalidCutPoints reports a
/// quantile-mode cut refresh whose new column would induce fewer cuts than
/// existing splits require; unlike the reference engine, which errors midway
/// through installation, nothing has been modified.
enum class PredictorUpdateResult { accepted, rolledBack, invalidCutPoints };

/// A sequential per-observation predictor update: stage one observation's
/// leaf moves, test that no leaf empties, then commit or skip, with a single
/// fits rebuild at the end. Type-erased so several samplers sharing an
/// index-aligned column can be swept jointly, committing all-or-none.
class PredictorUpdateSession {
public:
  virtual ~PredictorUpdateSession() = default;
  /// Stages observation i's leaf moves against the running occupancy counts;
  /// true unless some tree's leaf would be left empty.
  virtual bool observationWouldRemainValid(size_t i) = 0;
  /// Installs observation i's staged value and occupancy moves.
  virtual void commitObservation(size_t i) = 0;
  /// Re-routes all observations and rebuilds fits. The sequential guard
  /// admits no empty leaves, so false is an internal invariant violation.
  virtual bool finalize() = 0;
};

/// Posterior draws on the original response scale; caller-owned storage.
struct Results {
  double* sigma = nullptr;          // numSamples
  double* trainingFits = nullptr;   // numObservations x numSamples, or null
  double* testFits = nullptr;       // numTestObservations x numSamples, or null
  std::uint32_t* variableCounts = nullptr;  // numPredictors x numSamples, or null
  double* k = nullptr;              // numSamples, or null; only when k sampled
};

/// Single-chain conjugate backfitting sampler; a faithful port of the classic
/// engine's Gibbs iteration generic over the (integrable) leaf model.
template <IntegrableLeafModel L>
class Sampler {
public:
  Sampler(const double* x, const double* y, size_t numObservations,
          size_t numPredictors, const double* weights, const double* offset,
          bool responseIsBinary, double sigmaEstimate, double sigmaDf,
          double sigmaRawScale, const SamplerOptions& options, ext_rng* rng)
    : options_(options), weights_(weights), rng_(rng) {
    if (options.maxNumCutsPerVariable != nullptr) {
      data_.build(x, numObservations, numPredictors,
                  options.maxNumCutsPerVariable, options.useQuantiles);
    } else {
      data_.build(x, numObservations, numPredictors, options.maxNumCuts,
                  options.useQuantiles);
    }
    options_.maxNumCutsPerVariable = nullptr;  // borrowed; consumed by build

    if (responseIsBinary) {
      response_ = std::make_unique<ProbitResponse>(y, offset, weights,
                                                   numObservations);
    } else {
      response_ = std::make_unique<GaussianResponse>(
        y, offset, weights, numObservations, sigmaEstimate, sigmaDf,
        sigmaRawScale);
    }

    leaf_.scale =
      options.nodeScale / std::sqrt(static_cast<double>(options.numTrees));
    treePrior_.base = options.base;
    treePrior_.power = options.power;
    sigmaIsFixed_ = responseIsBinary;

    if (options.useDart) {
      dart_ = options.dart;
      dart_.initialize(numPredictors);
      treePrior_.splitProbabilities = dart_.probabilities.data();
      splitCounts_.resize(numPredictors);
    } else if (options.splitProbabilities != nullptr) {
      fixedSplitProbabilities_.assign(options.splitProbabilities,
                                      options.splitProbabilities + numPredictors);
      treePrior_.splitProbabilities = fixedSplitProbabilities_.data();
    }

    sigma_ = response_->initialSigma();
    k_ = options.k;

    indexBuffer_.resize(numObservations * options.numTrees);
    trees_.resize(options.numTrees);
    for (size_t t = 0; t < options.numTrees; ++t)
      trees_[t].initialize(indexBuffer_.data() + t * numObservations,
                           numObservations);

    treeFits_.assign(numObservations * options.numTrees, 0.0);
    totalFits_.assign(numObservations, 0.0);
    treeY_.resize(numObservations);
    currFits_.resize(numObservations);
    paramByNode_.clear();
  }

  void setTestPredictors(const double* x_test, size_t numTestObservations) {
    data_.buildTest(x_test, numTestObservations);
    totalTestFits_.assign(numTestObservations, 0.0);
    currTestFits_.resize(numTestObservations);
  }

  /// One thinning-free run; results slots may be null to skip recording.
  void run(size_t numBurnIn, size_t numSamples, Results& results) {
    size_t n = data_.numObservations;
    double* y = response_->workingResponse();

    for (size_t iteration = 0; iteration < numBurnIn + numSamples; ++iteration) {
      bool record = iteration >= numBurnIn;

      MoveContext ctx{data_,
                      treePrior_,
                      options_.birthOrDeathProbability,
                      options_.swapProbability,
                      options_.changeProbability,
                      options_.birthProbability,
                      weights_,
                      k_,
                      scratch_};

      kSumSquaredParams_ = 0.0;
      kNumLeaves_ = 0.0;

      if (record && data_.numTestObservations > 0)
        misc_setVectorToConstant(totalTestFits_.data(),
                                 data_.numTestObservations, 0.0);

      for (size_t t = 0; t < options_.numTrees; ++t) {
        double* oldTreeFits = treeFits_.data() + t * n;

        // treeY = y - (totalFits - oldTreeFits): the residual this tree owns
        std::memcpy(treeY_.data(), y, n * sizeof(double));
        misc_subtractVectorsInPlace(totalFits_.data(), n, treeY_.data());
        misc_addVectorsInPlace(oldTreeFits, n, treeY_.data());

        trees_[t].setNodeAverages(treeY_.data(), weights_);

        bool stepTaken;
        StepType stepType;
        metropolisJumpForTree(ctx, leaf_, rng_, trees_[t], treeY_.data(), sigma_,
                              &stepTaken, &stepType);

        sampleParametersAndSetFits(trees_[t], record);

        misc_subtractVectorsInPlace(oldTreeFits, n, totalFits_.data());
        misc_addVectorsInPlace(currFits_.data(), n, totalFits_.data());
        if (record && data_.numTestObservations > 0)
          misc_addVectorsInPlace(currTestFits_.data(), data_.numTestObservations,
                                 totalTestFits_.data());

        std::memcpy(oldTreeFits, currFits_.data(), n * sizeof(double));
      }

      response_->refreshLatents(rng_, totalFits_.data());
      y = response_->workingResponse();

      if (!sigmaIsFixed_)
        sigma_ = response_->drawSigma(rng_, totalFits_.data(), sigma_);

      if (options_.updateK)
        k_ = options_.kHyperprior.draw(rng_, kSumSquaredParams_, kNumLeaves_,
                                       leaf_.scale);

      if (options_.useDart) {
        std::memset(splitCounts_.data(), 0,
                    splitCounts_.size() * sizeof(std::uint32_t));
        for (size_t t = 0; t < options_.numTrees; ++t)
          trees_[t].countVariableUses(splitCounts_.data());
        dart_.update(rng_, splitCounts_.data());
      }

      if (record) {
        size_t sampleNum = iteration - numBurnIn;
        storeSample(results, sampleNum);
      }
    }
  }

  // Between-sample mutation; new-vector lifetimes are the caller's problem.
  void setOffset(const double* offset, bool updateScale) {
    response_->setOffset(offset, updateScale, &sigma_);
  }
  void setResponse(const double* y) {
    response_->setResponse(y, rng_, totalFits_.data(), &sigma_);
  }
  void setSigma(double sigmaOriginalScale) {
    sigma_ = sigmaOriginalScale / response_->sigmaScale();
  }
  const double* latents() const { return response_->latents(); }

  /// Replace the predictor matrix (borrowed, column-major; the old pointer is
  /// kept on failure). Unless forceUpdate, a leaf that would empty rolls the
  /// whole change back; forceUpdate instead collapses emptied leaves into
  /// their parents.
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) {
    size_t n = data_.numObservations;
    if (updateCutPoints)
      for (size_t j = 0; j < data_.numPredictors; ++j)
        if (!data_.cutsWouldRemainValid(j, newX + j * n))
          return PredictorUpdateResult::invalidCutPoints;

    const double* oldX = data_.x;
    std::vector<xint_t> oldCodes;
    std::vector<std::vector<double>> oldCuts;
    if (!forceUpdate) {
      oldCodes = data_.codes;
      if (updateCutPoints) oldCuts = data_.cutPoints;
    }

    data_.setPredictors(newX, updateCutPoints);

    if (forceUpdate) {
      forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t j = 0; j < data_.numPredictors; ++j)
          data_.quantizeTestColumn(j);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateTreesAndRebuildFits()) {
      data_.x = oldX;
      data_.codes = std::move(oldCodes);
      if (updateCutPoints) data_.cutPoints = std::move(oldCuts);
      for (size_t t = 0; t < options_.numTrees; ++t)
        trees_[t].repartitionSubtree(data_, 0);
      return PredictorUpdateResult::rolledBack;
    }

    if (updateCutPoints && data_.numTestObservations > 0)
      for (size_t j = 0; j < data_.numPredictors; ++j)
        data_.quantizeTestColumn(j);
    return PredictorUpdateResult::accepted;
  }

  /// Overwrite a subset of columns in place; newColumns is column-major,
  /// numObservations x numColumns. Same transaction semantics as
  /// setPredictor.
  PredictorUpdateResult updatePredictor(const double* newColumns,
                                        const size_t* columns,
                                        size_t numColumns, bool forceUpdate,
                                        bool updateCutPoints) {
    size_t n = data_.numObservations;
    if (updateCutPoints)
      for (size_t k = 0; k < numColumns; ++k)
        if (!data_.cutsWouldRemainValid(columns[k], newColumns + k * n))
          return PredictorUpdateResult::invalidCutPoints;

    std::vector<double> oldValues;
    std::vector<xint_t> oldCodes;
    std::vector<std::vector<double>> oldCuts;
    if (!forceUpdate) {
      oldValues.resize(n * numColumns);
      oldCodes.resize(n * numColumns);
      if (updateCutPoints) oldCuts.resize(numColumns);
      for (size_t k = 0; k < numColumns; ++k) {
        std::memcpy(oldValues.data() + k * n, data_.x + columns[k] * n,
                    n * sizeof(double));
        std::memcpy(oldCodes.data() + k * n, data_.codes.data() + columns[k] * n,
                    n * sizeof(xint_t));
        if (updateCutPoints) oldCuts[k] = data_.cutPoints[columns[k]];
      }
    }

    data_.setColumns(newColumns, columns, numColumns, updateCutPoints);

    if (forceUpdate) {
      forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t k = 0; k < numColumns; ++k)
          data_.quantizeTestColumn(columns[k]);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateTreesAndRebuildFits()) {
      double* x_mutable = const_cast<double*>(data_.x);
      for (size_t k = 0; k < numColumns; ++k) {
        size_t j = columns[k];
        std::memcpy(x_mutable + j * n, oldValues.data() + k * n,
                    n * sizeof(double));
        std::memcpy(data_.codes.data() + j * n, oldCodes.data() + k * n,
                    n * sizeof(xint_t));
        if (updateCutPoints) data_.cutPoints[j] = std::move(oldCuts[k]);
      }
      for (size_t t = 0; t < options_.numTrees; ++t)
        trees_[t].repartitionSubtree(data_, 0);
      return PredictorUpdateResult::rolledBack;
    }

    if (updateCutPoints && data_.numTestObservations > 0)
      for (size_t k = 0; k < numColumns; ++k)
        data_.quantizeTestColumn(columns[k]);
    return PredictorUpdateResult::accepted;
  }

  /// Install externally chosen cut points (ascending) for a subset of
  /// columns and unconditionally refresh the trees: splits that fall out of
  /// range or lose their observations collapse into their parents, exactly
  /// as a forced predictor update does.
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints, const size_t* columns,
                    size_t numColumns) {
    for (size_t k = 0; k < numColumns; ++k)
      data_.setCutPointsForColumn(columns[k], newCutPoints[k],
                                  numCutPoints[k]);
    forceRefreshTrees();
  }

  /// Install one column's new values observation-by-observation in random
  /// scan order, rolling back exactly those whose move would empty a leaf;
  /// installed must have room for a flag per observation. Returns finalize()
  /// validity, which the guard makes true by construction.
  bool updatePredictorPerObservation(const double* newColumn, size_t column,
                                     bool* installed) {
    size_t n = data_.numObservations;
    UpdateSessionImpl session(*this, newColumn, column);

    std::vector<size_t> scanOrder(n);
    ext_rng_drawPermutation(rng_, scanOrder.data(), n);

    for (size_t i = 0; i < n; ++i) {
      size_t j = scanOrder[i];
      bool valid = session.observationWouldRemainValid(j);
      installed[j] = valid;
      if (valid) session.commitObservation(j);
    }
    return session.finalize();
  }

  std::unique_ptr<PredictorUpdateSession> beginPredictorUpdate(
    const double* newColumn, size_t column) {
    return std::make_unique<UpdateSessionImpl>(*this, newColumn, column);
  }

  ext_rng* rng() const { return rng_; }

  double sigma() const { return sigma_; }
  double k() const { return k_; }
  bool kIsSampled() const { return options_.updateK; }
  const ColumnStore& data() const { return data_; }
  const Tree& tree(size_t t) const { return trees_[t]; }
  const std::vector<double>& treeFits() const { return treeFits_; }
  const std::vector<double>& totalFits() const { return totalFits_; }
  size_t numObservations() const { return data_.numObservations; }
  size_t numPredictors() const { return data_.numPredictors; }
  size_t numTestObservations() const { return data_.numTestObservations; }

private:
  /// Leaf parameters recovered from a tree's fits, indexed by arena node id;
  /// fits are constant within a leaf, so any member observation's fit is the
  /// parameter. Must run against partitions consistent with the fits, i.e.
  /// before any re-route.
  void recoverParametersFromFits(size_t t, std::vector<double>& paramByNode) {
    Tree& tree(trees_[t]);
    const double* treeFits = treeFits_.data() + t * data_.numObservations;

    paramByNode.assign(tree.nodes.size(), 0.0);
    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      if (node.numObservations() > 0)
        paramByNode[static_cast<size_t>(i)] = treeFits[tree.indices[node.begin]];
    }
  }

  void setTreeFitsFromParameters(size_t t, const std::vector<double>& paramByNode) {
    Tree& tree(trees_[t]);
    double* treeFits = treeFits_.data() + t * data_.numObservations;

    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      double param = paramByNode[static_cast<size_t>(i)];
      if (node.parent == invalidNode) {
        misc_setVectorToConstant(treeFits, node.numObservations(), param);
      } else {
        misc_setIndexedVectorToConstant(treeFits, tree.indices + node.begin,
                                        node.numObservations(), param);
      }
    }
  }

  /// After a predictor change: re-route every tree, failing without fits
  /// changes if any leaf empties (the caller rolls the data back and
  /// repartitions); on success rewrite tree fits from the preserved leaf
  /// parameters. Node averages are left stale; run() recomputes them from
  /// the current residuals before any use.
  bool revalidateTreesAndRebuildFits() {
    size_t n = data_.numObservations;
    std::vector<std::vector<double>> paramsByTree(options_.numTrees);

    bool allValid = true;
    for (size_t t = 0; t < options_.numTrees && allValid; ++t) {
      recoverParametersFromFits(t, paramsByTree[t]);
      trees_[t].repartitionSubtree(data_, 0);
      allValid = trees_[t].bottomNodesAreOccupied();
    }
    if (!allValid) return false;

    for (size_t t = 0; t < options_.numTrees; ++t) {
      double* treeFits = treeFits_.data() + t * n;
      misc_subtractVectorsInPlace(treeFits, n, totalFits_.data());
      setTreeFitsFromParameters(t, paramsByTree[t]);
      misc_addVectorsInPlace(treeFits, n, totalFits_.data());
    }
    return true;
  }

  /// Unconditional refresh: re-route and collapse any node an empty leaf
  /// leaves behind, merging leaf parameters into the collapsed node.
  void forceRefreshTrees() {
    size_t n = data_.numObservations;
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);

    std::vector<double> paramByNode;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      recoverParametersFromFits(t, paramByNode);
      trees_[t].repartitionSubtree(data_, 0);
      trees_[t].collapseEmptyNodes(weights_, paramByNode);
      setTreeFitsFromParameters(t, paramByNode);
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }
  }

  /// Port of the classic engine's InPlacePredictorUpdater: caches each
  /// observation's leaf and per-leaf occupancy by routing codes through the
  /// split structure, so staging a move is two descents and a count check.
  class UpdateSessionImpl final : public PredictorUpdateSession {
  public:
    UpdateSessionImpl(Sampler& sampler, const double* newColumn, size_t column)
      : sampler_(sampler), column_(column), newColumn_(newColumn) {
      size_t n = sampler_.data_.numObservations;
      size_t numTrees = sampler_.options_.numTrees;

      newCodes_.resize(n);
      for (size_t i = 0; i < n; ++i)
        newCodes_[i] = sampler_.data_.codeFor(column_, newColumn_[i]);

      observationLeaf_.resize(numTrees * n);
      leafCounts_.resize(numTrees);
      pendingNewLeaf_.assign(numTrees, invalidNode);
      pendingOldLeaf_.resize(numTrees);

      for (size_t t = 0; t < numTrees; ++t) {
        const Tree& tree(sampler_.trees_[t]);
        leafCounts_[t].assign(tree.nodes.size(), 0);
        for (size_t i = 0; i < n; ++i) {
          int32_t leaf = tree.findBottomNodeForObservation(sampler_.data_, i);
          observationLeaf_[t * n + i] = leaf;
          ++leafCounts_[t][static_cast<size_t>(leaf)];
        }
      }
    }

    bool observationWouldRemainValid(size_t i) override {
      size_t n = sampler_.data_.numObservations;
      bool valid = true;
      for (size_t t = 0; t < sampler_.options_.numTrees && valid; ++t) {
        const Tree& tree(sampler_.trees_[t]);
        int32_t oldLeaf = observationLeaf_[t * n + i];
        int32_t newLeaf = tree.findBottomNodeForObservation(
          sampler_.data_, i, static_cast<int32_t>(column_), newCodes_[i]);

        if (newLeaf == oldLeaf) {
          pendingNewLeaf_[t] = invalidNode;
        } else if (leafCounts_[t][static_cast<size_t>(oldLeaf)] == 1) {
          valid = false;  // last occupant; the move would empty the leaf
        } else {
          pendingOldLeaf_[t] = oldLeaf;
          pendingNewLeaf_[t] = newLeaf;
        }
      }
      return valid;
    }

    void commitObservation(size_t i) override {
      ColumnStore& data(sampler_.data_);
      size_t n = data.numObservations;
      const_cast<double*>(data.x)[i + column_ * n] = newColumn_[i];
      data.codes[i + column_ * n] = newCodes_[i];
      for (size_t t = 0; t < sampler_.options_.numTrees; ++t) {
        if (pendingNewLeaf_[t] != invalidNode) {
          --leafCounts_[t][static_cast<size_t>(pendingOldLeaf_[t])];
          ++leafCounts_[t][static_cast<size_t>(pendingNewLeaf_[t])];
        }
      }
    }

    bool finalize() override { return sampler_.revalidateTreesAndRebuildFits(); }

  private:
    Sampler& sampler_;
    size_t column_;
    const double* newColumn_;
    std::vector<xint_t> newCodes_;
    std::vector<int32_t> observationLeaf_;  // numTrees x numObservations
    std::vector<std::vector<std::uint32_t>> leafCounts_;  // arena-id indexed
    std::vector<int32_t> pendingNewLeaf_;  // invalidNode when no move staged
    std::vector<int32_t> pendingOldLeaf_;
  };

  void sampleParametersAndSetFits(Tree& tree, bool updateTestFits) {
    std::vector<int32_t>& bottoms(tree.bottomScratch);
    bottoms.clear();
    tree.fillBottom(0, bottoms);

    paramByNode_.assign(tree.nodes.size(), 0.0);
    for (int32_t i : bottoms) {
      const Node& node(tree.at(i));
      double param = node.numObservations() == 0
        ? 0.0
        : leaf_.drawFromPosterior(rng_, k_, node.average,
                                  node.numEffectiveObservations, sigma_ * sigma_);
      paramByNode_[static_cast<size_t>(i)] = param;

      if (options_.updateK) {
        kSumSquaredParams_ += param * param;
        kNumLeaves_ += 1.0;
      }

      if (node.parent == invalidNode) {
        misc_setVectorToConstant(currFits_.data(), node.numObservations(), param);
      } else {
        misc_setIndexedVectorToConstant(currFits_.data(),
                                        tree.indices + node.begin,
                                        node.numObservations(), param);
      }
    }

    if (updateTestFits && data_.numTestObservations > 0) {
      for (size_t i = 0; i < data_.numTestObservations; ++i) {
        int32_t leafIndex = tree.findBottomNodeForRow(data_.testRow(i));
        currTestFits_[i] = paramByNode_[static_cast<size_t>(leafIndex)];
      }
    }
  }

  void storeSample(Results& results, size_t sampleNum) {
    size_t n = data_.numObservations;
    double scale = response_->fitScale();
    double shift = response_->fitShift();

    if (results.sigma != nullptr)
      results.sigma[sampleNum] = sigma_ * response_->sigmaScale();

    if (results.k != nullptr) results.k[sampleNum] = k_;

    if (results.trainingFits != nullptr) {
      double* out = results.trainingFits + sampleNum * n;
      for (size_t i = 0; i < n; ++i) out[i] = scale * totalFits_[i] + shift;
      // caller adds any offset back; the engine never sees original-scale y
    }

    if (results.testFits != nullptr && data_.numTestObservations > 0) {
      double* out = results.testFits + sampleNum * data_.numTestObservations;
      for (size_t i = 0; i < data_.numTestObservations; ++i)
        out[i] = scale * totalTestFits_[i] + shift;
    }

    if (results.variableCounts != nullptr) {
      std::uint32_t* out =
        results.variableCounts + sampleNum * data_.numPredictors;
      std::memset(out, 0, data_.numPredictors * sizeof(std::uint32_t));
      for (size_t t = 0; t < options_.numTrees; ++t)
        trees_[t].countVariableUses(out);
    }
  }

  SamplerOptions options_;
  ColumnStore data_;
  const double* weights_;
  ext_rng* rng_;

  L leaf_;
  CGMTreePrior treePrior_;
  DartPrior dart_;
  std::vector<double> fixedSplitProbabilities_;
  std::vector<std::uint32_t> splitCounts_;
  std::unique_ptr<ResponseModel> response_;
  bool sigmaIsFixed_ = false;
  double sigma_ = 1.0;
  double k_ = 2.0;
  double kSumSquaredParams_ = 0.0;
  double kNumLeaves_ = 0.0;

  std::vector<Tree> trees_;
  std::vector<size_t> indexBuffer_;
  std::vector<double> treeFits_;
  std::vector<double> totalFits_, totalTestFits_;
  std::vector<double> treeY_, currFits_, currTestFits_;
  std::vector<double> paramByNode_;
  MoveScratch scratch_;
};

using ClassicSampler = Sampler<ConstantGaussianLeaf>;

}  // namespace bartcore

#endif  // BARTCORE_SAMPLER_HPP
