#ifndef BARTCORE_SAMPLER_HPP
#define BARTCORE_SAMPLER_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <thread>
#include <vector>

#include <external/random.h>

#include "chain.hpp"
#include "data.hpp"
#include "model.hpp"

namespace bartcore {

/// Outcome of a transactional predictor change. invalidCutPoints reports a
/// quantile-mode cut refresh whose new column would induce fewer cuts than
/// existing splits require; unlike the reference engine, which errors midway
/// through installation, nothing has been modified.
enum class PredictorUpdateResult { accepted, rolledBack, invalidCutPoints };

/// A whole sampler's serializable state: per-chain states, the store's cut
/// points (setCutPoints may have replaced the ones creation induces), and
/// the saved-tree write position.
struct SamplerStateData {
  std::vector<ChainStateData> chains;
  std::vector<std::vector<double>> cutPoints;  // empty vector per categorical column
  size_t currentSampleNum = 0;
};

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

/// The sampler proper: a shared column store and one or more chains over it.
/// Chains run independently (optionally on worker threads) and hold all
/// per-chain state; the sampler owns the data and orchestrates transactions,
/// which are all-or-none across every tree of every chain, exactly as the
/// classic engine's were.
template <IntegrableLeafModel L>
class Sampler {
public:
  /// rngs supplies one generator per chain (options.numChains of them). With
  /// numThreads > 1 no rng may call into R (or otherwise share state).
  /// weights apply only to the gaussian family.
  Sampler(const double* x, const double* y, size_t numObservations,
          size_t numPredictors, const double* weights, const double* offset,
          ResponseFamily family, double sigmaEstimate, double sigmaDf,
          double sigmaRawScale, const SamplerOptions& options,
          ext_rng* const* rngs)
    : options_(options), family_(family) {
    if (options.maxNumCutsPerVariable != nullptr) {
      data_.build(x, numObservations, numPredictors,
                  options.maxNumCutsPerVariable, options.useQuantiles,
                  options.columnTypes);
    } else {
      data_.build(x, numObservations, numPredictors, options.maxNumCuts,
                  options.useQuantiles, options.columnTypes);
    }
    options_.maxNumCutsPerVariable = nullptr;  // borrowed; consumed by build
    options_.columnTypes = nullptr;

    chains_.reserve(options.numChains);
    for (size_t c = 0; c < options.numChains; ++c)
      chains_.push_back(std::make_unique<Chain<L>>(
        data_, y, weights, offset, family, sigmaEstimate, sigmaDf,
        sigmaRawScale, options_, rngs[c]));

    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
  }

  // chains reference the store member, so the sampler's address is pinned
  Sampler(const Sampler&) = delete;
  Sampler& operator=(const Sampler&) = delete;

  void setTestPredictors(const double* x_test, size_t numTestObservations) {
    data_.buildTest(x_test, numTestObservations);
    for (auto& chain : chains_) chain->resizeTestStorage();
  }

  /// Runs every chain numBurnIn + numSamples iterations, filling per-chain
  /// slabs of the (chain-major) results arrays; see Results. Chains execute
  /// on up to min(numThreads, numChains) worker threads.
  void run(size_t numBurnIn, size_t numSamples, Results& results) {
    size_t numChains = chains_.size();
    for (auto& chain : chains_) chain->setSavedSlotBase(currentSampleNum_);
    std::vector<Results> chainResults(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      Results& r(chainResults[c]);
      if (results.sigma != nullptr) r.sigma = results.sigma + c * numSamples;
      if (results.k != nullptr) r.k = results.k + c * numSamples;
      if (results.trainingFits != nullptr)
        r.trainingFits =
          results.trainingFits + c * numSamples * data_.numObservations;
      if (results.testFits != nullptr)
        r.testFits =
          results.testFits + c * numSamples * data_.numTestObservations;
      if (results.variableCounts != nullptr)
        r.variableCounts =
          results.variableCounts + c * numSamples * data_.numPredictors;
    }

    size_t numWorkers = options_.numThreads < numChains ? options_.numThreads
                                                        : numChains;
    if (numWorkers <= 1) {
      for (size_t c = 0; c < numChains; ++c)
        chains_[c]->run(numBurnIn, numSamples, chainResults[c]);
    } else {
      std::vector<std::thread> workers;
      workers.reserve(numWorkers);
      for (size_t w = 0; w < numWorkers; ++w) {
        workers.emplace_back([this, w, numWorkers, numChains, numBurnIn,
                              numSamples, &chainResults]() {
          for (size_t c = w; c < numChains; c += numWorkers)
            chains_[c]->run(numBurnIn, numSamples, chainResults[c]);
        });
      }
      for (std::thread& worker : workers) worker.join();
    }

    size_t capacity = savedTreeCapacity();
    if (capacity > 0 && numSamples > 0)
      currentSampleNum_ = (currentSampleNum_ + numSamples) % capacity;
  }

  // Saved trees (keepTrees) and prediction.

  size_t savedTreeCapacity() const {
    return chains_[0]->savedTreeCapacity();
  }
  size_t currentSampleNum() const { return currentSampleNum_; }
  void setCurrentSampleNum(size_t sampleNum) { currentSampleNum_ = sampleNum; }

  const std::vector<FlatNode>& savedTree(size_t chainNum, size_t slot,
                                         size_t treeNum) const {
    return chains_[chainNum]->savedTree(slot, treeNum);
  }
  void flattenTree(size_t chainNum, size_t treeNum,
                   std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts) {
    chains_[chainNum]->flattenTree(treeNum, nodes, counts);
  }

  /// Fits for raw column-major test rows, on the original response scale
  /// (offsets are the caller's problem). With saved trees, out is
  /// numTestObservations x savedTreeCapacity x numChains, chain-major like
  /// Results; without, one slab per chain from the live trees.
  void predict(const double* x_test, size_t numTestObservations, double* out) {
    size_t capacity = savedTreeCapacity();
    for (size_t c = 0; c < chains_.size(); ++c) {
      if (capacity > 0) {
        for (size_t slot = 0; slot < capacity; ++slot)
          chains_[c]->predictFromSavedSample(
            slot, x_test, numTestObservations,
            out + (c * capacity + slot) * numTestObservations);
      } else {
        chains_[c]->predictFromCurrentTrees(x_test, numTestObservations,
                                            out + c * numTestObservations);
      }
    }
  }

  // State serialization: getState captures everything needed to reconstruct
  // the sampler's posterior state in a fresh instance over the same data;
  // setState validates every chain against the state's cut points before
  // touching anything, so failure leaves the sampler unchanged.

  void getState(SamplerStateData& state) {
    state.chains.resize(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->getState(state.chains[c]);
    state.cutPoints = data_.cutPoints;
    state.currentSampleNum = currentSampleNum_;
  }

  bool setState(const SamplerStateData& state) {
    if (state.chains.size() != chains_.size()) return false;
    if (state.cutPoints.size() != data_.numPredictors) return false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.types[j] == ColumnType::categorical) {
        if (!state.cutPoints[j].empty()) return false;
      } else if (state.cutPoints[j].empty() ||
                 state.cutPoints[j].size() > 65535) {
        return false;
      } else {
        for (size_t k = 1; k < state.cutPoints[j].size(); ++k)
          if (state.cutPoints[j][k] <= state.cutPoints[j][k - 1]) return false;
      }
    }

    // install the state's cuts, snapshotting for rollback: tree validity is
    // defined against them
    std::vector<std::vector<double>> oldCutPoints(data_.cutPoints);
    std::vector<std::uint32_t> oldNumCuts(data_.numCuts);
    std::vector<std::uint32_t> oldMaxNumCuts(data_.maxNumCuts);
    std::vector<xint_t> oldCodes(data_.codes);
    std::vector<xint_t> oldTestCodes(data_.testCodes);
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.types[j] == ColumnType::categorical) continue;
      data_.setCutPointsForColumn(j, state.cutPoints[j].data(),
                                  static_cast<std::uint32_t>(
                                    state.cutPoints[j].size()));
    }

    bool allValid = true;
    for (size_t c = 0; c < chains_.size() && allValid; ++c)
      allValid = chains_[c]->stateIsValid(state.chains[c]);

    if (!allValid) {
      data_.cutPoints = std::move(oldCutPoints);
      data_.numCuts = std::move(oldNumCuts);
      data_.maxNumCuts = std::move(oldMaxNumCuts);
      data_.codes = std::move(oldCodes);
      data_.testCodes = std::move(oldTestCodes);
      return false;
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      if (!chains_[c]->setState(state.chains[c])) return false;
    size_t capacity = savedTreeCapacity();
    currentSampleNum_ = capacity > 0 ? state.currentSampleNum % capacity : 0;
    return true;
  }

  // Between-sample mutation, fanned out to every chain; new-vector lifetimes
  // are the caller's problem.
  void setOffset(const double* offset, bool updateScale) {
    for (auto& chain : chains_) chain->setOffset(offset, updateScale);
  }
  void setResponse(const double* y) {
    for (auto& chain : chains_) chain->setResponse(y);
  }
  void setSigma(double sigmaOriginalScale) {
    for (auto& chain : chains_) chain->setSigma(sigmaOriginalScale);
  }
  const double* latents(size_t chainNum = 0) const {
    return chains_[chainNum]->latents();
  }

  /// Replace the entire data set: predictors, response, and optionally
  /// weights, offset, and test predictors, with a possibly different number
  /// of observations (all borrowed; the predictor count is fixed). Not
  /// transactional: cut points are rebuilt from scratch, existing splits are
  /// remapped onto the value-nearest new cuts, and any subtree left invalid
  /// or empty collapses, as in the classic engine. Gaussian chains keep
  /// sigma and the variance prior fixed on the original scale.
  void setData(const double* x, const double* y, size_t numObservations,
               const double* weights, const double* offset,
               const double* x_test, size_t numTestObservations) {
    // recover parameters against the old fits and partitions before anything
    // moves; the old cut values drive the split remap
    std::vector<typename Chain<L>::TreeParameters> params(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->recoverTreeParameters(params[c]);

    std::vector<std::vector<double>> oldCutPoints(data_.cutPoints);

    data_.setData(x, numObservations);
    if (x_test != nullptr && numTestObservations > 0)
      data_.buildTest(x_test, numTestObservations);
    else
      data_.clearTest();

    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->applyNewData(y, weights, offset, oldCutPoints, params[c]);
  }

  /// Replace the predictor matrix (borrowed, column-major; the old pointer is
  /// kept on failure). Unless forceUpdate, a leaf that would empty in any
  /// tree of any chain rolls the whole change back; forceUpdate instead
  /// collapses emptied leaves into their parents.
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) {
    size_t n = data_.numObservations;
    // categorical columns must hold representable codes whether or not cut
    // points refresh
    for (size_t j = 0; j < data_.numPredictors; ++j)
      if ((updateCutPoints ||
           data_.types[j] == ColumnType::categorical) &&
          !data_.cutsWouldRemainValid(j, newX + j * n))
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
      for (auto& chain : chains_) chain->forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t j = 0; j < data_.numPredictors; ++j)
          data_.quantizeTestColumn(j);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateAllChains()) {
      data_.x = oldX;
      data_.codes = std::move(oldCodes);
      if (updateCutPoints) data_.cutPoints = std::move(oldCuts);
      for (auto& chain : chains_) chain->repartitionTrees();
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
    for (size_t k = 0; k < numColumns; ++k)
      if ((updateCutPoints ||
           data_.types[columns[k]] == ColumnType::categorical) &&
          !data_.cutsWouldRemainValid(columns[k], newColumns + k * n))
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
      for (auto& chain : chains_) chain->forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t k = 0; k < numColumns; ++k)
          data_.quantizeTestColumn(columns[k]);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateAllChains()) {
      double* x_mutable = const_cast<double*>(data_.x);
      for (size_t k = 0; k < numColumns; ++k) {
        size_t j = columns[k];
        std::memcpy(x_mutable + j * n, oldValues.data() + k * n,
                    n * sizeof(double));
        std::memcpy(data_.codes.data() + j * n, oldCodes.data() + k * n,
                    n * sizeof(xint_t));
        if (updateCutPoints) data_.cutPoints[j] = std::move(oldCuts[k]);
      }
      for (auto& chain : chains_) chain->repartitionTrees();
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
    for (auto& chain : chains_) chain->forceRefreshTrees();
  }

  /// Install one column's new values observation-by-observation in random
  /// scan order, rolling back exactly those whose move would empty a leaf in
  /// any tree of any chain; installed must have room for a flag per
  /// observation. Returns finalize() validity, which the guard makes true by
  /// construction.
  bool updatePredictorPerObservation(const double* newColumn, size_t column,
                                     bool* installed) {
    size_t n = data_.numObservations;
    UpdateSessionImpl session(*this, newColumn, column);

    std::vector<size_t> scanOrder(n);
    ext_rng_drawPermutation(chains_[0]->rng(), scanOrder.data(), n);

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

  ext_rng* rng() const { return chains_[0]->rng(); }

  ResponseFamily family() const { return family_; }
  double sigma(size_t chainNum = 0) const { return chains_[chainNum]->sigma(); }
  double k(size_t chainNum = 0) const { return chains_[chainNum]->k(); }
  bool kIsSampled() const { return options_.updateK; }
  const ColumnStore& data() const { return data_; }
  const Chain<L>& chain(size_t chainNum) const { return *chains_[chainNum]; }
  size_t numChains() const { return chains_.size(); }
  size_t numTrees() const { return options_.numTrees; }
  size_t numObservations() const { return data_.numObservations; }
  size_t numPredictors() const { return data_.numPredictors; }
  size_t numTestObservations() const { return data_.numTestObservations; }

private:
  /// Two-phase transaction over every chain: validate all trees of all
  /// chains first, then rebuild fits only if everything holds, so a failure
  /// in a late chain never leaves an early chain's fits overwritten.
  bool revalidateAllChains() {
    size_t numChains = chains_.size();
    std::vector<typename Chain<L>::TreeParameters> params(numChains);

    bool allValid = true;
    for (size_t c = 0; c < numChains && allValid; ++c)
      allValid = chains_[c]->revalidateTrees(params[c]);
    if (!allValid) return false;

    for (size_t c = 0; c < numChains; ++c)
      chains_[c]->rebuildFitsFromParameters(params[c]);
    return true;
  }

  /// Port of the classic engine's InPlacePredictorUpdater: caches each
  /// observation's leaf and per-leaf occupancy for every tree of every
  /// chain by routing codes through the split structure, so staging a move
  /// is a descent and a count check per tree.
  class UpdateSessionImpl final : public PredictorUpdateSession {
  public:
    UpdateSessionImpl(Sampler& sampler, const double* newColumn, size_t column)
      : sampler_(sampler), column_(column), newColumn_(newColumn) {
      size_t n = sampler_.data_.numObservations;
      size_t totalNumTrees =
        sampler_.options_.numTrees * sampler_.chains_.size();

      newCodes_.resize(n);
      for (size_t i = 0; i < n; ++i)
        newCodes_[i] = sampler_.data_.codeFor(column_, newColumn_[i]);

      observationLeaf_.resize(totalNumTrees * n);
      leafCounts_.resize(totalNumTrees);
      pendingNewLeaf_.assign(totalNumTrees, invalidNode);
      pendingOldLeaf_.resize(totalNumTrees);

      for (size_t t = 0; t < totalNumTrees; ++t) {
        const Tree& tree(treeAt(t));
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
      for (size_t t = 0; t < leafCounts_.size() && valid; ++t) {
        const Tree& tree(treeAt(t));
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
      for (size_t t = 0; t < leafCounts_.size(); ++t) {
        if (pendingNewLeaf_[t] != invalidNode) {
          --leafCounts_[t][static_cast<size_t>(pendingOldLeaf_[t])];
          ++leafCounts_[t][static_cast<size_t>(pendingNewLeaf_[t])];
        }
      }
    }

    bool finalize() override { return sampler_.revalidateAllChains(); }

  private:
    const Tree& treeAt(size_t t) const {
      size_t numTrees = sampler_.options_.numTrees;
      return sampler_.chains_[t / numTrees]->tree(t % numTrees);
    }

    Sampler& sampler_;
    size_t column_;
    const double* newColumn_;
    std::vector<xint_t> newCodes_;
    std::vector<int32_t> observationLeaf_;  // totalNumTrees x numObservations
    std::vector<std::vector<std::uint32_t>> leafCounts_;  // arena-id indexed
    std::vector<int32_t> pendingNewLeaf_;  // invalidNode when no move staged
    std::vector<int32_t> pendingOldLeaf_;
  };

  SamplerOptions options_;
  ResponseFamily family_;
  ColumnStore data_;
  std::vector<std::unique_ptr<Chain<L>>> chains_;
  size_t currentSampleNum_ = 0;  // next saved-tree slot, wrapping circularly
};

using ClassicSampler = Sampler<ConstantGaussianLeaf>;

}  // namespace bartcore

#endif  // BARTCORE_SAMPLER_HPP
