#ifndef BARTCORE_CHAIN_HPP
#define BARTCORE_CHAIN_HPP

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
  size_t numChains = 1;
  // worker threads for running chains concurrently; only min(numThreads,
  // numChains) are used, and every chain needs its own non-R rng when > 1
  size_t numThreads = 1;
  // every numThin-th iteration is kept; numBurnIn and numSamples count at
  // the kept rate, as in the classic engine
  size_t numThin = 1;
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
  // borrowed per-column types, null for all-ordinal; copied during
  // construction. Categorical columns must hold integer codes 0..K-1 with
  // K <= maxCategories; the caller validates.
  const ColumnType* columnTypes = nullptr;

  // split-variable selection: fixed weights (borrowed; normalized over
  // available variables at each node) or DART; both null/false = uniform
  const double* splitProbabilities = nullptr;
  bool useDart = false;
  DartPrior dart;

  // k is fixed at .k unless updateK, in which case .k is the starting value
  bool updateK = false;
  ChiKHyperprior kHyperprior;

  // when set, every kept sample's trees are flattened into a circular buffer
  // of numSamplesToStore slots (at least 1) per chain, for prediction and
  // reporting after the run; the classic engine's keepTrees
  bool keepTrees = false;
  size_t numSamplesToStore = 0;
};

/// Posterior draws on the original response scale; caller-owned storage.
/// With several chains the arrays hold one slab per chain, chain-major:
/// sigma and k are numSamples x numChains, fits and counts add their leading
/// dimension (e.g. trainingFits is numObservations x numSamples x numChains).
struct Results {
  double* sigma = nullptr;          // numSamples
  double* trainingFits = nullptr;   // numObservations x numSamples, or null
  double* testFits = nullptr;       // numTestObservations x numSamples, or null
  std::uint32_t* variableCounts = nullptr;  // numPredictors x numSamples, or null
  double* k = nullptr;              // numSamples, or null; only when k sampled
};

/// Everything a chain's posterior state comprises, in host-exchangeable
/// form: value-encoded flattened trees, the saved-tree buffer, sigma, k,
/// response latents, DART state, and the serialized rng. Three fields exist
/// so a restored chain continues bitwise identically rather than merely
/// equivalently: sigma is on the engine's internal scale (the original-scale
/// round trip can drop a bit), totalFits carries the running total's
/// floating-point accumulation history, and indices carries each tree's
/// observation ordering, whose within-leaf order the sufficient-statistic
/// sums depend on. totalFits and indices may be left empty; restoration then
/// recomputes them canonically, exact only as far as the last ulp.
struct ChainStateData {
  std::vector<std::vector<FlatNode>> trees;
  std::vector<std::vector<FlatNode>> savedTrees;  // slot-major; empty unless kept
  std::vector<double> totalFits;      // numObservations, or empty
  std::vector<size_t> indices;        // numObservations x numTrees, or empty
  double sigma = 1.0;
  double k = 2.0;
  std::vector<double> latents;            // empty for gaussian
  std::vector<double> dartProbabilities;  // empty when DART is off
  double dartAlpha = 1.0;
  size_t dartNumUpdatesSkipped = 0;
  std::vector<unsigned char> rngState;
};

/// One MCMC chain of the conjugate backfitting sampler: trees, fits, chain
/// parameters, and its own response state and rng, over a shared read-only
/// ColumnStore. Data mutation is orchestrated one level up (Sampler), which
/// calls the tree-refresh methods here; a chain never writes to the store.
template <IntegrableLeafModel L>
class Chain {
public:
  Chain(const ColumnStore& data, const double* y, const double* weights,
        const double* offset, ResponseFamily family, double sigmaEstimate,
        double sigmaDf, double sigmaRawScale, const SamplerOptions& options,
        ext_rng* rng)
    : options_(options), data_(data), weights_(weights), rng_(rng) {
    size_t numObservations = data.numObservations;
    options_.maxNumCutsPerVariable = nullptr;  // consumed by the store build
    options_.columnTypes = nullptr;

    switch (family) {
    case ResponseFamily::probit:
      response_ = std::make_unique<ProbitResponse>(y, offset, numObservations);
      break;
    case ResponseFamily::logistic:
      response_ =
        std::make_unique<LogisticResponse>(y, offset, numObservations);
      break;
    case ResponseFamily::gaussian:
      response_ = std::make_unique<GaussianResponse>(
        y, offset, weights, numObservations, sigmaEstimate, sigmaDf,
        sigmaRawScale);
      break;
    }

    leaf_.scale =
      options.nodeScale / std::sqrt(static_cast<double>(options.numTrees));
    treePrior_.base = options.base;
    treePrior_.power = options.power;
    sigmaIsFixed_ = family != ResponseFamily::gaussian;

    if (options.useDart) {
      dart_ = options.dart;
      dart_.initialize(data.numPredictors);
      treePrior_.splitProbabilities = dart_.probabilities.data();
      splitCounts_.resize(data.numPredictors);
    } else if (options.splitProbabilities != nullptr) {
      fixedSplitProbabilities_.assign(
        options.splitProbabilities,
        options.splitProbabilities + data.numPredictors);
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
    resizeTestStorage();
  }

  /// Called after the shared store's test data changes.
  void resizeTestStorage() {
    totalTestFits_.assign(data_.numTestObservations, 0.0);
    currTestFits_.resize(data_.numTestObservations);
  }

  /// One thinning-free run; results slots may be null to skip recording.
  /// Touches only chain state and the read-only store: safe to run chains
  /// concurrently as long as each has its own rng that never calls into R.
  void run(size_t numBurnIn, size_t numSamples, Results& results) {
    size_t n = data_.numObservations;
    size_t numThin = options_.numThin;
    double* y = response_->workingResponse();
    // per-iteration precisions: user weights for gaussian, the current
    // Polya-Gamma draws for logistic (refreshed with the latents)
    const double* weights = response_->workingWeights();

    size_t totalIterations = (numBurnIn + numSamples) * numThin;
    for (size_t iteration = 0; iteration < totalIterations; ++iteration) {
      bool record = (iteration + 1) % numThin == 0 &&
                    iteration / numThin >= numBurnIn;
      size_t sampleNum = record ? iteration / numThin - numBurnIn : 0;

      MoveContext ctx{data_,
                      treePrior_,
                      options_.birthOrDeathProbability,
                      options_.swapProbability,
                      options_.changeProbability,
                      options_.birthProbability,
                      weights,
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

        trees_[t].setNodeAverages(treeY_.data(), weights);

        bool stepTaken;
        StepType stepType;
        metropolisJumpForTree(ctx, leaf_, rng_, trees_[t], treeY_.data(), sigma_,
                              &stepTaken, &stepType);

        sampleParametersAndSetFits(trees_[t], record);

        // flatten while the freshly drawn parameters are live in paramByNode_
        if (record && savedTreeCapacity_ > 0) {
          size_t slot = (savedSlotBase_ + sampleNum) % savedTreeCapacity_;
          trees_[t].flatten(data_, paramByNode_.data(),
                            savedTrees_[slot * options_.numTrees + t]);
        }

        misc_subtractVectorsInPlace(oldTreeFits, n, totalFits_.data());
        misc_addVectorsInPlace(currFits_.data(), n, totalFits_.data());
        if (record && data_.numTestObservations > 0)
          misc_addVectorsInPlace(currTestFits_.data(), data_.numTestObservations,
                                 totalTestFits_.data());

        std::memcpy(oldTreeFits, currFits_.data(), n * sizeof(double));
      }

      response_->refreshLatents(rng_, totalFits_.data());
      y = response_->workingResponse();
      weights = response_->workingWeights();

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

      if (record) storeSample(results, sampleNum);
    }
  }

  // Between-sample mutation; new-vector lifetimes are the caller's problem.
  void setOffset(const double* offset, bool updateScale) {
    response_->setOffset(offset, updateScale, &sigma_);
  }
  void setWeights(const double* weights) {
    weights_ = weights;
    response_->setWeights(weights);
  }
  void setResponse(const double* y) {
    response_->setResponse(y, rng_, totalFits_.data(), &sigma_);
  }
  void setSigma(double sigmaOriginalScale) {
    sigma_ = sigmaOriginalScale / response_->sigmaScale();
  }
  const double* latents() const { return response_->latents(); }

  /// Replace every tree's structure with a draw from the tree prior over the
  /// current cut grid, empty leaves collapsed, exactly as the reference
  /// engine: fits are left stale, which run() tolerates because totalFits
  /// still sums the per-tree fits.
  void sampleTreesFromPrior() {
    size_t n = data_.numObservations;
    const double* y = response_->workingResponse();
    const double* weights = response_->workingWeights();
    std::vector<double> paramByNode;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      trees_[t].initialize(indexBuffer_.data() + t * n, n);
      growSubtreeFromPrior(trees_[t], 0, y, weights);
      paramByNode.assign(trees_[t].nodes.size(), 0.0);
      trees_[t].collapseEmptyNodes(weights, paramByNode);
    }
  }

  /// Replace every leaf parameter with a draw from the node prior and
  /// rebuild the tree, total, and test fits to match.
  void sampleNodeParametersFromPrior() {
    size_t n = data_.numObservations;
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);
    if (data_.numTestObservations > 0)
      misc_setVectorToConstant(totalTestFits_.data(),
                               data_.numTestObservations, 0.0);

    for (size_t t = 0; t < options_.numTrees; ++t) {
      Tree& tree(trees_[t]);
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      paramByNode_.assign(tree.nodes.size(), 0.0);
      for (int32_t i : tree.bottomScratch)
        paramByNode_[static_cast<size_t>(i)] = leaf_.drawFromPrior(rng_, k_);

      setTreeFitsFromParameters(t, paramByNode_);
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
      for (size_t i = 0; i < data_.numTestObservations; ++i) {
        int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
        totalTestFits_[i] += paramByNode_[static_cast<size_t>(leafIndex)];
      }
    }
  }

  // Tree refreshes after the shared store changed, driven by the sampler.
  // The validate and rebuild phases are split so a transaction can check
  // every chain before any chain's fits are overwritten.

  using TreeParameters = std::vector<std::vector<double>>;

  /// Recover leaf parameters from fits, re-route every tree against the
  /// store's current codes, and report whether all leaves stay occupied.
  bool revalidateTrees(TreeParameters& params) {
    params.resize(options_.numTrees);
    bool allValid = true;
    for (size_t t = 0; t < options_.numTrees && allValid; ++t) {
      recoverParametersFromFits(t, params[t]);
      trees_[t].repartitionSubtree(data_, 0);
      allValid = trees_[t].bottomNodesAreOccupied();
    }
    return allValid;
  }

  /// Second phase of a successful transaction: rewrite tree fits from the
  /// parameters revalidateTrees recovered. Node averages are left stale;
  /// run() recomputes them from current residuals before any use.
  void rebuildFitsFromParameters(const TreeParameters& params) {
    size_t n = data_.numObservations;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      double* treeFits = treeFits_.data() + t * n;
      misc_subtractVectorsInPlace(treeFits, n, totalFits_.data());
      setTreeFitsFromParameters(t, params[t]);
      misc_addVectorsInPlace(treeFits, n, totalFits_.data());
    }
  }

  /// Rollback re-route: restore partitions consistent with the store after
  /// the sampler has restored its old codes.
  void repartitionTrees() {
    for (size_t t = 0; t < options_.numTrees; ++t)
      trees_[t].repartitionSubtree(data_, 0);
  }

  /// First phase of a whole-data replacement: recover every tree's leaf
  /// parameters against the current fits and partitions, before the store or
  /// any per-observation storage moves.
  void recoverTreeParameters(TreeParameters& params) {
    params.resize(options_.numTrees);
    for (size_t t = 0; t < options_.numTrees; ++t)
      recoverParametersFromFits(t, params[t]);
  }

  /// Second phase, after the shared store holds the new predictors and
  /// freshly rebuilt cuts: swap the response state (sigma is preserved on
  /// the original scale), resize per-observation storage, remap split
  /// indices onto the new cut grid, re-route, and collapse anything left
  /// invalid or empty. Node averages are left stale; run() recomputes them.
  void applyNewData(const double* y, const double* weights,
                    const double* offset,
                    const std::vector<std::vector<double>>& oldCutPoints,
                    TreeParameters& params) {
    size_t n = data_.numObservations;
    bool numObservationsChanged = n * options_.numTrees != indexBuffer_.size();

    weights_ = weights;
    response_->setData(y, offset, weights, n, &sigma_);

    if (numObservationsChanged) {
      indexBuffer_.resize(n * options_.numTrees);
      treeFits_.resize(n * options_.numTrees);
      totalFits_.resize(n);
      treeY_.resize(n);
      currFits_.resize(n);
    }
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);

    for (size_t t = 0; t < options_.numTrees; ++t) {
      trees_[t].mapOldCutPointsOntoNew(data_, oldCutPoints, params[t]);
      if (numObservationsChanged)
        trees_[t].resetObservations(indexBuffer_.data() + t * n, n);
      trees_[t].repartitionSubtree(data_, 0);
      trees_[t].collapseEmptyNodes(response_->workingWeights(), params[t]);
      setTreeFitsFromParameters(t, params[t]);
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }

    resizeTestStorage();
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
      trees_[t].collapseEmptyNodes(response_->workingWeights(), paramByNode);
      setTreeFitsFromParameters(t, paramByNode);
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }
  }

  // Saved-tree (keepTrees) storage: a circular buffer of capacity slots,
  // each one kept sample's forest in flattened form. The slot base is set by
  // the sampler before every run so chains write consistent slots without
  // sharing mutable state.

  void initializeSavedTrees(size_t capacity) {
    savedTreeCapacity_ = capacity;
    savedTrees_.assign(capacity * options_.numTrees,
                       std::vector<FlatNode>(1));
  }
  void setSavedSlotBase(size_t base) { savedSlotBase_ = base; }
  size_t savedTreeCapacity() const { return savedTreeCapacity_; }
  const std::vector<FlatNode>& savedTree(size_t slot, size_t t) const {
    return savedTrees_[slot * options_.numTrees + t];
  }

  /// Flatten live tree t with parameters recovered from the current fits;
  /// counts receive the current partition sizes.
  void flattenTree(size_t t, std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts) {
    std::vector<double> params;
    recoverParametersFromFits(t, params);
    trees_[t].flatten(data_, params.data(), nodes, &counts);
  }

  /// Info dump of live tree t, parameters recovered from the current fits.
  void printTree(size_t t, int indentation) {
    std::vector<double> params;
    recoverParametersFromFits(t, params);
    trees_[t].print(data_, params.data(), indentation);
  }

  /// The same for one saved tree, in the reference engine's saved format.
  void printSavedTree(size_t slot, size_t t, int indentation) const {
    const std::vector<FlatNode>& flat(savedTrees_[slot * options_.numTrees + t]);
    printFlatSubtree(data_, flat.data(), indentation);
  }

  /// Fits for raw column-major test rows from one saved sample's trees, on
  /// the original response scale; offsets are the caller's problem.
  void predictFromSavedSample(size_t slot, const double* x_test,
                              size_t numTestObservations, double* out) const {
    std::vector<size_t> indices(numTestObservations);
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    for (size_t t = 0; t < options_.numTrees; ++t) {
      const std::vector<FlatNode>& flat(
        savedTrees_[slot * options_.numTrees + t]);
      for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
      addFlatPredictionsBelow(flat.data(), data_.types.data(), x_test,
                              numTestObservations, indices.data(), 0,
                              numTestObservations, out);
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  /// The same from the live trees, parameters recovered from current fits.
  void predictFromCurrentTrees(const double* x_test,
                               size_t numTestObservations, double* out) {
    std::vector<size_t> indices(numTestObservations);
    std::vector<double> params;
    std::vector<FlatNode> flat;
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    for (size_t t = 0; t < options_.numTrees; ++t) {
      recoverParametersFromFits(t, params);
      trees_[t].flatten(data_, params.data(), flat);
      for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
      addFlatPredictionsBelow(flat.data(), data_.types.data(), x_test,
                              numTestObservations, indices.data(), 0,
                              numTestObservations, out);
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  // Chain state serialization. getState captures everything the posterior
  // state comprises; stateIsValid checks a candidate against the store's
  // current cuts without mutating anything, so a multi-chain restore can be
  // all-or-none; setState trusts that check.

  void getState(ChainStateData& state) {
    state.trees.resize(options_.numTrees);
    std::vector<double> params;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      recoverParametersFromFits(t, params);
      trees_[t].flatten(data_, params.data(), state.trees[t]);
    }
    state.savedTrees = savedTrees_;
    state.totalFits = totalFits_;
    state.indices = indexBuffer_;
    state.sigma = sigma_;
    state.k = k_;
    if (response_->latents() != nullptr) {
      state.latents.assign(response_->latents(),
                           response_->latents() + data_.numObservations);
    } else {
      state.latents.clear();
    }
    if (options_.useDart) {
      state.dartProbabilities = dart_.probabilities;
      state.dartAlpha = dart_.alpha;
      state.dartNumUpdatesSkipped = dart_.numUpdatesSkipped();
    } else {
      state.dartProbabilities.clear();
    }
    state.rngState.resize(ext_rng_getSerializedStateLength(rng_));
    if (!state.rngState.empty())
      ext_rng_writeSerializedState(rng_, state.rngState.data());
  }

  bool stateIsValid(const ChainStateData& state) const {
    if (state.trees.size() != options_.numTrees) return false;
    if (!state.savedTrees.empty() &&
        state.savedTrees.size() != savedTrees_.size())
      return false;
    for (const std::vector<FlatNode>& saved : state.savedTrees)
      if (!flatTreeIsWellFormed(data_, saved.data(), saved.size()))
        return false;
    size_t n = data_.numObservations;
    if (!state.totalFits.empty() && state.totalFits.size() != n) return false;
    if (!state.indices.empty() &&
        state.indices.size() != n * options_.numTrees)
      return false;
    if (!state.latents.empty() &&
        (response_->latents() == nullptr || state.latents.size() != n))
      return false;
    if (options_.useDart && !state.dartProbabilities.empty() &&
        state.dartProbabilities.size() != data_.numPredictors)
      return false;
    if (!state.rngState.empty() &&
        state.rngState.size() != ext_rng_getSerializedStateLength(rng_))
      return false;

    Tree scratch;
    std::vector<size_t> scratchIndices(n);
    std::vector<bool> seen(n);
    std::vector<double> params;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      scratch.initialize(scratchIndices.data(), n);
      if (!scratch.buildFromFlat(data_, state.trees[t].data(),
                                 state.trees[t].size(), params))
        return false;
      if (state.indices.empty()) {
        scratch.repartitionSubtree(data_, 0);
      } else {
        // the stored ordering must permute the observations and respect
        // every node's rule
        const size_t* treeIndices = state.indices.data() + t * n;
        seen.assign(n, false);
        for (size_t i = 0; i < n; ++i) {
          if (treeIndices[i] >= n || seen[treeIndices[i]]) return false;
          seen[treeIndices[i]] = true;
        }
        std::memcpy(scratchIndices.data(), treeIndices, n * sizeof(size_t));
        if (!scratch.setPartitionsFromOrderedIndices(data_, 0)) return false;
      }
      if (!scratch.bottomNodesAreOccupied()) return false;
    }
    return true;
  }

  /// Installs a state stateIsValid accepted; false only on the invariant
  /// violation of a validated tree failing to rebuild.
  bool setState(const ChainStateData& state) {
    size_t n = data_.numObservations;
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);
    std::vector<double> params;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      trees_[t].initialize(indexBuffer_.data() + t * n, n);
      if (!trees_[t].buildFromFlat(data_, state.trees[t].data(),
                                   state.trees[t].size(), params))
        return false;
      if (state.indices.empty()) {
        trees_[t].repartitionSubtree(data_, 0);
      } else {
        std::memcpy(indexBuffer_.data() + t * n, state.indices.data() + t * n,
                    n * sizeof(size_t));
        if (!trees_[t].setPartitionsFromOrderedIndices(data_, 0)) return false;
      }
      setTreeFitsFromParameters(t, params);
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }
    if (!state.totalFits.empty())
      std::memcpy(totalFits_.data(), state.totalFits.data(),
                  n * sizeof(double));
    if (!state.savedTrees.empty()) savedTrees_ = state.savedTrees;
    sigma_ = state.sigma;
    k_ = state.k;
    if (!state.latents.empty())
      response_->restoreLatents(state.latents.data());
    if (options_.useDart && !state.dartProbabilities.empty()) {
      // the tree prior points at this vector's storage; overwrite in place
      std::memcpy(dart_.probabilities.data(), state.dartProbabilities.data(),
                  state.dartProbabilities.size() * sizeof(double));
      dart_.alpha = state.dartAlpha;
      dart_.setNumUpdatesSkipped(state.dartNumUpdatesSkipped);
    }
    if (!state.rngState.empty())
      ext_rng_readSerializedState(rng_, state.rngState.data());
    return true;
  }

  ext_rng* rng() const { return rng_; }
  /// On the original response scale, symmetric with setSigma.
  double sigma() const { return sigma_ * response_->sigmaScale(); }
  double k() const { return k_; }
  size_t numTrees() const { return options_.numTrees; }
  const Tree& tree(size_t t) const { return trees_[t]; }
  const std::vector<double>& treeFits() const { return treeFits_; }
  const std::vector<double>& totalFits() const { return totalFits_; }

private:
  /// The reference engine's recursion: growth is Bernoulli in the
  /// depth-decayed prior probability, rules come from the prior, and empty
  /// children keep growing (availability is rule-based) until the caller
  /// collapses them.
  void growSubtreeFromPrior(Tree& tree, int32_t nodeIndex, const double* y,
                            const double* weights) {
    double growthProbability =
      treePrior_.growthProbability(tree, data_, nodeIndex);
    if (growthProbability <= 0.0 ||
        ext_rng_simulateBernoulli(rng_, growthProbability) == 0)
      return;

    Rule rule = treePrior_.drawRuleAndVariable(tree, data_, rng_, nodeIndex);
    tree.birth(data_, nodeIndex, rule, y, weights);
    int32_t leftChild = tree.at(nodeIndex).leftChild;
    growSubtreeFromPrior(tree, leftChild, y, weights);
    growSubtreeFromPrior(tree, leftChild + 1, y, weights);
  }

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
        int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
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
      if (data_.testOffset != nullptr)
        misc_addVectorsInPlace(data_.testOffset, data_.numTestObservations,
                               out);
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
  const ColumnStore& data_;
  // user weights, held only to forward to the response model; the engine's
  // per-iteration precisions come from response_->workingWeights()
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
  std::vector<std::vector<FlatNode>> savedTrees_;  // slot-major, capacity x numTrees
  size_t savedTreeCapacity_ = 0;
  size_t savedSlotBase_ = 0;
  std::vector<size_t> indexBuffer_;
  std::vector<double> treeFits_;
  std::vector<double> totalFits_, totalTestFits_;
  std::vector<double> treeY_, currFits_, currTestFits_;
  std::vector<double> paramByNode_;
  MoveScratch scratch_;
};

}  // namespace bartcore

#endif  // BARTCORE_CHAIN_HPP
