#ifndef BARTCORE_CHAIN_HPP
#define BARTCORE_CHAIN_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>
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

  // sparse (CSC / dgCMatrix-layout) predictor ingestion: when
  // cscColumnPointers is non-null the x creation argument is ignored and
  // columns build from the triple (all ordinal; borrowed for the sampler's
  // lifetime; the host validates row indices unique and in range). See
  // docs/design/sparse-columns.md.
  const int* cscColumnPointers = nullptr;  // length of the CSC source + 1
  const int* cscRowIndices = nullptr;
  const double* cscValues = nullptr;

  // mixed dense/sparse ingestion: when columnSources is non-null, column j
  // builds from column columnSources[j] of mixedDenseValues (column-major)
  // when nonnegative - categorical allowed, rawColumn served - or from CSC
  // column ~columnSources[j] of the triple above otherwise (ordinal only).
  // Both borrowed for the sampler's lifetime; the host validates. See
  // docs/design/sparse-columns.md.
  const double* mixedDenseValues = nullptr;
  const std::int32_t* columnSources = nullptr;  // length numPredictors

  // linear leaves: the ordinal predictor columns entering every leaf's
  // regression (borrowed; consumed during construction). Empty designates
  // the constant leaf; the factory validates count, range, and type.
  const std::size_t* leafCovariateColumns = nullptr;
  std::size_t numLeafCovariates = 0;

  // GP (function-valued) leaves share the leafCovariateColumns designation;
  // gpLeaves selects them over the linear leaf at the factory (the leaf
  // model type consumes the designation either way). gpLengthscales, when
  // non-null, fixes the kernel lengthscales per designated column on the
  // standardized scale (borrowed; consumed during construction); null
  // selects the median pairwise-distance heuristic. Leaves larger than
  // gpMaxLeafSize score and draw as constant leaves, confining the
  // O(n_leaf^3) kernel math below the cap.
  bool gpLeaves = false;
  const double* gpLengthscales = nullptr;
  std::size_t gpMaxLeafSize = 256;

  // grouped random intercepts (rbart_vi's Gibbs blocks in-core): one group
  // index 0..numGroups-1 per training observation (borrowed; consumed during
  // construction). numGroups == 0 leaves the response ungrouped.
  // tauPriorScale is the original-scale relative scale (sd(y) continuous,
  // 0.5 binary); tauSliceSteps is the per-update slice-step count, the R
  // loop's n.thin coupling.
  const std::uint32_t* groupIndices = nullptr;
  std::size_t numGroups = 0;
  TauPriorKind tauPriorKind = TauPriorKind::cauchy;
  double tauPriorScale = 1.0;
  std::size_t tauSliceSteps = 1;

  // split-variable selection: fixed weights (borrowed; normalized over
  // available variables at each node) or DART; both null/false = uniform
  const double* splitProbabilities = nullptr;
  bool useDart = false;
  DartPrior dart;

  // k is fixed at .k unless updateK, in which case .k is the starting value
  bool updateK = false;
  ChiKHyperprior kHyperprior;

  // gaussian responses only: sigma holds at the constructor's sigmaEstimate
  // and is never drawn (the R-level fixed() residual prior); binary families
  // are always sigma-fixed
  bool sigmaIsFixed = false;

  // when set, every kept sample's trees are flattened into a circular buffer
  // of numSamplesToStore slots (at least 1) per chain, for prediction and
  // reporting after the run; the classic engine's keepTrees
  bool keepTrees = false;
  size_t numSamplesToStore = 0;

  // progress reporting during runs, in the classic engine's format: one
  // "iteration: k (of N)" line every printEvery kept iterations
  bool verbose = false;
  std::uint32_t printEvery = 100;
};

/// Receives chains' formatted progress lines during a run. Runs on worker
/// threads hand lines to a queue the main thread flushes (workers must never
/// call into R); inline runs print directly.
struct ProgressSink {
  virtual ~ProgressSink() = default;
  virtual void report(const char* line) = 0;
};

/// A between-run prior replacement, the classic engine's setModel: every
/// field is applied unconditionally except that DART samplers keep their
/// Dirichlet split machinery (a dbartsModel cannot express DART) and the
/// fixed-sigma binary families ignore the variance prior. Installing a model
/// before any run matches creating with it.
struct ModelParameters {
  double base = 0.95, power = 2.0;
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  double nodeScale = 0.5;
  // k is fixed at .k unless updateK, matching SamplerOptions
  double k = 2.0;
  bool updateK = false;
  ChiKHyperprior kHyperprior;
  // variance prior, anchored to the original-scale sigma estimate exactly
  // as construction is; rawScale is qchisq(1 - quantile, df) / df. With
  // sigmaIsFixed, sigmaEstimate is instead the fixed original-scale sd and
  // sigma is set to it and never drawn again
  bool sigmaIsFixed = false;
  double sigmaEstimate = 1.0;
  double sigmaDf = 3.0;
  double sigmaRawScale = 1.0;
  const double* splitProbabilities = nullptr;  // borrowed; copied on install
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
  // numPredictors x numSamples, or null; filled only under DART
  double* splitProbabilities = nullptr;
  // grouped samplers only, both on the original response scale:
  // tau is numSamples, groupEffects numGroups x numSamples
  double* tau = nullptr;
  double* groupEffects = nullptr;
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
  // vector-parameter leaves only: each tree's slopes, numParams - 1 doubles
  // per leaf in pre-order (the intercept rides in the flattened value);
  // savedTreeParams parallels savedTrees. Empty for scalar leaf models.
  std::vector<std::vector<double>> treeParams;
  std::vector<std::vector<double>> savedTreeParams;
  // wide categorical columns only: each flat tree's mask side channel,
  // parallel to trees/savedTrees. Empty otherwise.
  std::vector<std::vector<std::uint64_t>> treeMasks;
  std::vector<std::vector<std::uint64_t>> savedTreeMasks;
  std::vector<double> totalFits;      // numObservations, or empty
  std::vector<size_t> indices;        // numObservations x numTrees, or empty
  double sigma = 1.0;
  double k = 2.0;
  // the gaussian response transform at capture; max <= min marks scale-free
  double fitMin = 0.0, fitMax = 0.0;
  // the variance prior's internal-scale value; re-anchoring it through the
  // transform is a multiply-divide round trip that can perturb the last
  // bit, so a restore installs this exactly. Negative when not applicable.
  double sigmaPriorScale = -1.0;
  std::vector<double> latents;            // empty for gaussian
  // grouped samplers only, internal scale so restores are exact
  std::vector<double> groupEffects;
  double groupTau = 0.0;
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

    if constexpr (L::hasVectorParams)
      leaf_.initialize(data, options.leafCovariateColumns,
                       options.numLeafCovariates);
    else if constexpr (L::hasFunctionParams)
      leaf_.initialize(data, options.leafCovariateColumns,
                       options.numLeafCovariates, options.gpLengthscales,
                       options.gpMaxLeafSize, options.numChains);
    options_.leafCovariateColumns = nullptr;  // consumed above
    options_.gpLengthscales = nullptr;

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
    // grouped random intercepts decorate the base family; initialization
    // draws b from its prior through this chain's generator
    if (options.numGroups > 0)
      response_ = std::make_unique<GroupedResponse>(
        std::move(response_), numObservations, options.groupIndices,
        options.numGroups, options.tauPriorKind, options.tauPriorScale,
        options.tauSliceSteps, rng);
    options_.groupIndices = nullptr;  // consumed above

    leaf_.scale =
      options.nodeScale / std::sqrt(static_cast<double>(options.numTrees));
    treePrior_.base = options.base;
    treePrior_.power = options.power;
    family_ = family;
    sigmaIsFixed_ = family != ResponseFamily::gaussian || options.sigmaIsFixed;

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
    paramByNode_.clear();
    if constexpr (L::hasVectorParams)
      paramsByTree_.assign(options.numTrees,
                           std::vector<double>(leaf_.numParams(), 0.0));
    resizeTestStorage();
  }

  /// Called after the shared store's test data changes.
  void resizeTestStorage() {
    totalTestFits_.assign(data_.numTestObservations, 0.0);
    currTestFits_.resize(data_.numTestObservations);
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      leaf_.rebuildTestCovariates(data_);
  }

  /// One thinning-free run; results slots may be null to skip recording.
  /// Touches only chain state and the read-only store: safe to run chains
  /// concurrently as long as each has its own rng that never calls into R.
  /// progress, when non-null under verbose, receives one formatted line per
  /// printEvery kept iterations.
  void run(size_t numBurnIn, size_t numSamples, Results& results,
           ProgressSink* progress = nullptr, size_t chainIndex = 0) {
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

      if (options_.verbose && progress != nullptr &&
          (iteration + 1) % numThin == 0 &&
          (iteration / numThin + 1) % options_.printEvery == 0) {
        char line[128];
        if (options_.numChains > 1) {
          std::snprintf(line, sizeof(line),
                        "[%lu] iteration: %lu (of %lu)\n",
                        static_cast<unsigned long>(chainIndex + 1),
                        static_cast<unsigned long>(iteration + 1),
                        static_cast<unsigned long>(totalIterations));
        } else {
          std::snprintf(line, sizeof(line), "iteration: %lu (of %lu)\n",
                        static_cast<unsigned long>(iteration + 1),
                        static_cast<unsigned long>(totalIterations));
        }
        progress->report(line);
      }

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

        // treeY = (y - totalFits) + oldTreeFits is the residual this tree owns,
        // in the prior three-op order so it stays bitwise identical. The same
        // pass strips this tree's old fits from the running total, so the draw
        // can overwrite its slab in place: totalFits now holds every other tree.
        {
          const double* __restrict y_ = y;
          double* __restrict total = totalFits_.data();
          const double* __restrict oldFits = oldTreeFits;
          double* __restrict treeY = treeY_.data();
          for (size_t i = 0; i < n; ++i) {
            treeY[i] = y_[i] - total[i] + oldFits[i];
            total[i] -= oldFits[i];
          }
        }

        // constant-leaf node means, recomputed against this sweep's residual.
        // Linear (vector) leaves accumulate their own per-node statistics and
        // never read node.average, so skip it for them; function (gp) leaves
        // still need it because over-cap nodes delegate to the constant leaf.
        if constexpr (!L::hasVectorParams)
          trees_[t].setNodeAverages(treeY_.data(), weights);

        bool stepTaken;
        StepType stepType;
        metropolisJumpForTree(ctx, leaf_, rng_, trees_[t], treeY_.data(), sigma_,
                              &stepTaken, &stepType);
        // accepted changes and deaths strand pooled mask words; no rule
        // copies are live here, so this is a safe point to reclaim them
        if (data_.hasPooledCategorical)
          trees_[t].compactMaskPoolIfNeeded(data_);

        // the draw writes this tree's new fits straight into its slab
        sampleParametersAndSetFits(t, oldTreeFits, record);

        // flatten while the freshly drawn parameters are live
        if (record && savedTreeCapacity_ > 0) {
          size_t slot = (savedSlotBase_ + sampleNum) % savedTreeCapacity_;
          std::vector<std::uint64_t>* masks = data_.hasWideCategorical
            ? &savedTreeMasks_[slot * options_.numTrees + t] : nullptr;
          if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
            trees_[t].flatten(data_, paramByNode_.data(),
                              savedTrees_[slot * options_.numTrees + t],
                              nullptr, 1, nullptr, masks);
          } else if constexpr (L::hasVectorParams) {
            trees_[t].flatten(data_, paramsByTree_[t].data(),
                              savedTrees_[slot * options_.numTrees + t],
                              nullptr, leaf_.numParams(),
                              &savedTreeParams_[slot * options_.numTrees + t],
                              masks);
          } else {
            // function-valued leaves: records carry per-leaf mean fits for
            // reporting, the side channel carries the draw cache's alpha
            // weights plus covariate rows - the exact values the recorded
            // test fits used, so saved replays bit-match them
            functionLeafValues(trees_[t], oldTreeFits, paramByNode_);
            trees_[t].flatten(data_, paramByNode_.data(),
                              savedTrees_[slot * options_.numTrees + t],
                              nullptr, 1, nullptr, masks);
            std::vector<double>& blocks(
              savedTreeParams_[slot * options_.numTrees + t]);
            blocks.clear();
            for (int32_t i : trees_[t].bottomScratch)
              leaf_.appendLeafBlockFromCache(trees_[t], i, blocks);
          }
        }

        // totalFits += this tree's new fits, restoring the running total. The
        // old fits were already removed above, so this equals the prior
        // (totalFits - oldFits) + newFits exactly: totalFits is bitwise
        // identical, and the draw already left the new fits in the slab (no
        // trailing memcpy).
        {
          double* __restrict total = totalFits_.data();
          const double* __restrict newFits = oldTreeFits;
          for (size_t i = 0; i < n; ++i)
            total[i] += newFits[i];
        }
        if (record && data_.numTestObservations > 0)
          misc_addVectorsInPlace(currTestFits_.data(), data_.numTestObservations,
                                 totalTestFits_.data());
      }

      response_->refreshLatents(rng_, totalFits_.data(), sigma_);
      y = response_->workingResponse();
      weights = response_->workingWeights();

      if (!sigmaIsFixed_)
        sigma_ = response_->drawSigma(rng_, totalFits_.data(), sigma_);

      // a zero sum of squares under an infinite prior scale would make the
      // gamma rate zero and k infinite; with no information, keep k
      if (options_.updateK && kSumSquaredParams_ > 0.0)
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

  void setNumThin(size_t numThin) { options_.numThin = numThin; }
  void setVerbose(bool verbose, std::uint32_t printEvery) {
    options_.verbose = verbose;
    options_.printEvery = printEvery;
  }
  /// The multiplier taking internal-scale fits to the original response
  /// scale: the response range for gaussian, 1 for the binary families.
  double fitScale() const { return response_->fitScale(); }

  /// Install a replacement prior (the classic engine's setModel); see
  /// ModelParameters for the semantics.
  void setModel(const ModelParameters& model) {
    treePrior_.base = model.base;
    treePrior_.power = model.power;
    options_.birthOrDeathProbability = model.birthOrDeathProbability;
    options_.swapProbability = model.swapProbability;
    options_.changeProbability = model.changeProbability;
    options_.birthProbability = model.birthProbability;
    leaf_.scale =
      model.nodeScale / std::sqrt(static_cast<double>(options_.numTrees));

    options_.updateK = model.updateK;
    if (model.updateK) {
      options_.kHyperprior = model.kHyperprior;
    } else {
      k_ = model.k;
    }

    if (family_ == ResponseFamily::gaussian) {
      sigmaIsFixed_ = model.sigmaIsFixed;
      if (model.sigmaIsFixed)
        setSigma(model.sigmaEstimate);
      else
        response_->setSigmaPrior(model.sigmaEstimate, model.sigmaDf,
                                 model.sigmaRawScale);
    }

    if (!options_.useDart) {
      if (model.splitProbabilities == nullptr) {
        fixedSplitProbabilities_.clear();
        treePrior_.splitProbabilities = nullptr;
      } else {
        fixedSplitProbabilities_.assign(
          model.splitProbabilities,
          model.splitProbabilities + data_.numPredictors);
        treePrior_.splitProbabilities = fixedSplitProbabilities_.data();
      }
    }
  }

  /// Sum of squared working residuals, descaled to the original response
  /// scale (binary families report on the latent scale).
  double sumOfSquaredResiduals() {
    double result = misc_htm_computeSumOfSquaredResiduals(
      nullptr, 0, response_->workingResponse(), data_.numObservations,
      totalFits_.data());
    return result * response_->sigmaScale() * response_->sigmaScale();
  }

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
      // fresh structures carry zero parameter blocks until the next draw
      if constexpr (L::hasVectorParams)
        paramsByTree_[t].assign(trees_[t].nodes.size() * leaf_.numParams(),
                                0.0);
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
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
        paramByNode_.assign(tree.nodes.size(), 0.0);
        for (int32_t i : tree.bottomScratch)
          paramByNode_[static_cast<size_t>(i)] = leaf_.drawFromPrior(rng_, k_);

        setTreeFitsFromParameters(t, paramByNode_);
        misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
        for (size_t i = 0; i < data_.numTestObservations; ++i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          totalTestFits_[i] += paramByNode_[static_cast<size_t>(leafIndex)];
        }
      } else if constexpr (L::hasFunctionParams) {
        // prior fits land directly in the tree's fit slab; the per-node
        // prediction cache serves the routed test rows
        double* treeFits = treeFits_.data() + t * n;
        leaf_.beginTreeDraw(tree);
        for (int32_t i : tree.bottomScratch)
          leaf_.drawFromPriorForNode(rng_, tree, k_, i, treeFits);
        misc_addVectorsInPlace(treeFits, n, totalFits_.data());
        for (size_t i = 0; i < data_.numTestObservations; ++i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          totalTestFits_[i] +=
            leaf_.fitForTestObservationForNode(tree, leafIndex, i);
        }
      } else {
        size_t numParams = leaf_.numParams();
        std::vector<double>& treeParams(paramsByTree_[t]);
        treeParams.assign(tree.nodes.size() * numParams, 0.0);
        for (int32_t i : tree.bottomScratch)
          leaf_.drawFromPrior(rng_, k_,
                              treeParams.data() +
                                static_cast<size_t>(i) * numParams);

        setTreeFitsFromParameterBlocks(t, treeParams);
        misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
        for (size_t i = 0; i < data_.numTestObservations; ++i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          totalTestFits_[i] += leaf_.fitForTestObservation(
            treeParams.data() + static_cast<size_t>(leafIndex) * numParams,
            i);
        }
      }
    }
  }

  // Tree refreshes after the shared store changed, driven by the sampler.
  // The validate and rebuild phases are split so a transaction can check
  // every chain before any chain's fits are overwritten.

  using TreeParameters = std::vector<std::vector<double>>;

  /// Recover leaf parameters (from fits for scalar leaves, from the
  /// persisted blocks for vector ones; function-valued leaves keep their
  /// per-observation fits in place, so nothing is recovered), re-route every
  /// tree against the store's current codes, and report whether all leaves
  /// stay occupied.
  bool revalidateTrees(TreeParameters& params) {
    params.resize(options_.numTrees);
    bool allValid = true;
    for (size_t t = 0; t < options_.numTrees && allValid; ++t) {
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams)
        recoverParametersFromFits(t, params[t]);
      else if constexpr (L::hasVectorParams)
        params[t] = paramsByTree_[t];
      else
        params[t].clear();
      trees_[t].repartitionSubtree(data_, 0);
      allValid = trees_[t].bottomNodesAreOccupied();
    }
    return allValid;
  }

  /// Second phase of a successful transaction: rewrite tree fits from the
  /// parameters revalidateTrees recovered. Node averages are left stale;
  /// run() recomputes them from current residuals before any use.
  /// Function-valued leaves only refresh the covariate gather: their
  /// per-observation fits are the parameters and stay in place (the next
  /// sweep's draws replace them under the new values).
  void rebuildFitsFromParameters(const TreeParameters& params) {
    size_t n = data_.numObservations;
    if constexpr (L::hasFunctionParams) {
      (void) n;
      (void) params;
      leaf_.regatherTrainingCovariates(data_);
      return;
    }
    // vector leaves read raw covariate values: pick up the installed ones
    if constexpr (L::hasVectorParams)
      leaf_.regatherTrainingCovariates(data_);
    for (size_t t = 0; t < options_.numTrees; ++t) {
      double* treeFits = treeFits_.data() + t * n;
      misc_subtractVectorsInPlace(treeFits, n, totalFits_.data());
      if constexpr (!L::hasVectorParams)
        setTreeFitsFromParameters(t, params[t]);
      else
        setTreeFitsFromParameterBlocks(t, params[t]);
      misc_addVectorsInPlace(treeFits, n, totalFits_.data());
    }
  }

  /// Rollback re-route: restore partitions consistent with the store after
  /// the sampler has restored its old codes.
  void repartitionTrees() {
    // the restored raw values re-gather to exactly the old covariates
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      leaf_.regatherTrainingCovariates(data_);
    for (size_t t = 0; t < options_.numTrees; ++t)
      trees_[t].repartitionSubtree(data_, 0);
  }

  /// First phase of a whole-data replacement: recover every tree's leaf
  /// parameters against the current fits and partitions, before the store or
  /// any per-observation storage moves.
  void recoverTreeParameters(TreeParameters& params) {
    params.resize(options_.numTrees);
    for (size_t t = 0; t < options_.numTrees; ++t) {
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams)
        recoverParametersFromFits(t, params[t]);
      else if constexpr (L::hasVectorParams)
        params[t] = paramsByTree_[t];
      else
        params[t].clear();
    }
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
    }
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);

    // fresh standardization constants over the replacement data, like the
    // rebuilt cut grid; the persisted parameters carry over as-is, the same
    // approximate continuation the split remap embodies. Function-valued
    // parameters are per-observation fits over the OLD data, so they
    // cold-start at zero instead (the next sweep's draws replace them);
    // structures still carry through the split remap.
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      leaf_.reinitialize(data_);
    if constexpr (L::hasFunctionParams)
      for (size_t t = 0; t < options_.numTrees; ++t)
        params[t].assign(trees_[t].nodes.size(), 0.0);

    size_t paramStride = 1;
    if constexpr (L::hasVectorParams) paramStride = leaf_.numParams();

    for (size_t t = 0; t < options_.numTrees; ++t) {
      trees_[t].mapOldCutPointsOntoNew(data_, oldCutPoints, params[t],
                                       paramStride);
      if (numObservationsChanged)
        trees_[t].resetObservations(indexBuffer_.data() + t * n, n);
      trees_[t].repartitionSubtree(data_, 0);
      trees_[t].collapseEmptyNodes(response_->workingWeights(), params[t],
                                   paramStride);
      if constexpr (!L::hasVectorParams) {
        setTreeFitsFromParameters(t, params[t]);
      } else {
        paramsByTree_[t] = params[t];
        setTreeFitsFromParameterBlocks(t, params[t]);
      }
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }

    resizeTestStorage();
  }

  /// Unconditional refresh: re-route and collapse any node an empty leaf
  /// leaves behind, merging leaf parameters into the collapsed node.
  /// Function-valued leaves keep their per-observation fits (they remain a
  /// coherent state under any partition) and collapse structure only.
  void forceRefreshTrees() {
    size_t n = data_.numObservations;
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);

    if constexpr (L::hasFunctionParams) {
      leaf_.regatherTrainingCovariates(data_);
      std::vector<double> dummyParams;
      for (size_t t = 0; t < options_.numTrees; ++t) {
        trees_[t].repartitionSubtree(data_, 0);
        dummyParams.assign(trees_[t].nodes.size(), 0.0);
        trees_[t].collapseEmptyNodes(response_->workingWeights(), dummyParams);
        misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
      }
    } else if constexpr (!L::hasVectorParams) {
      std::vector<double> paramByNode;
      for (size_t t = 0; t < options_.numTrees; ++t) {
        recoverParametersFromFits(t, paramByNode);
        trees_[t].repartitionSubtree(data_, 0);
        trees_[t].collapseEmptyNodes(response_->workingWeights(), paramByNode);
        setTreeFitsFromParameters(t, paramByNode);
        misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
      }
    } else {
      leaf_.regatherTrainingCovariates(data_);
      size_t numParams = leaf_.numParams();
      for (size_t t = 0; t < options_.numTrees; ++t) {
        trees_[t].repartitionSubtree(data_, 0);
        trees_[t].collapseEmptyNodes(response_->workingWeights(),
                                     paramsByTree_[t], numParams);
        setTreeFitsFromParameterBlocks(t, paramsByTree_[t]);
        misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
      }
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
    // the default slot is a single zero leaf; its slopes are zero too, and
    // a function-valued slot holds one zero-constant block
    if constexpr (L::hasVectorParams)
      savedTreeParams_.assign(
        capacity * options_.numTrees,
        std::vector<double>(leaf_.numParams() - 1, 0.0));
    else if constexpr (L::hasFunctionParams)
      savedTreeParams_.assign(capacity * options_.numTrees,
                              std::vector<double>{0.0, 0.0});
    if (data_.hasWideCategorical)
      savedTreeMasks_.assign(capacity * options_.numTrees,
                             std::vector<std::uint64_t>());
  }
  void setSavedSlotBase(size_t base) { savedSlotBase_ = base; }
  size_t savedTreeCapacity() const { return savedTreeCapacity_; }
  const std::vector<FlatNode>& savedTree(size_t slot, size_t t) const {
    return savedTrees_[slot * options_.numTrees + t];
  }
  /// Slopes of one saved tree (vector-parameter leaves), parallel to
  /// savedTree's pre-order leaves.
  const std::vector<double>& savedTreeSlopes(size_t slot, size_t t) const {
    return savedTreeParams_[slot * options_.numTrees + t];
  }
  /// Flattened mask words of one saved tree (wide categorical columns).
  const std::vector<std::uint64_t>& savedTreeMasks(size_t slot,
                                                   size_t t) const {
    return savedTreeMasks_[slot * options_.numTrees + t];
  }

  /// Flatten live tree t; counts receive the current partition sizes.
  /// Scalar leaves recover parameters from the fits; vector leaves emit
  /// their persisted blocks, the slopes (numParams - 1 per leaf, pre-order)
  /// into slopes when non-null.
  void flattenTree(size_t t, std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes = nullptr,
                   std::vector<std::uint64_t>* masks = nullptr) {
    if constexpr (L::hasFunctionParams) {
      // records carry per-leaf mean fits; slopes, when requested, receives
      // the saved-format side-channel blocks recomputed from the fits
      size_t n = data_.numObservations;
      std::vector<double> values;
      functionLeafValues(trees_[t], treeFits_.data() + t * n, values);
      trees_[t].flatten(data_, values.data(), nodes, &counts, 1, nullptr,
                        masks);
      if (slopes != nullptr) {
        slopes->clear();
        for (int32_t i : trees_[t].bottomScratch)
          leaf_.appendLeafBlock(trees_[t], i, treeFits_.data() + t * n,
                                *slopes);
      }
    } else if constexpr (!L::hasVectorParams) {
      if (slopes != nullptr) slopes->clear();
      std::vector<double> params;
      recoverParametersFromFits(t, params);
      trees_[t].flatten(data_, params.data(), nodes, &counts, 1, nullptr,
                        masks);
    } else {
      trees_[t].flatten(data_, paramsByTree_[t].data(), nodes, &counts,
                        leaf_.numParams(), slopes, masks);
    }
  }

  /// Info dump of live tree t; function-valued leaves print their per-leaf
  /// mean fits.
  void printTree(size_t t, int indentation) {
    if constexpr (L::hasFunctionParams) {
      std::vector<double> values;
      functionLeafValues(trees_[t],
                         treeFits_.data() + t * data_.numObservations,
                         values);
      trees_[t].print(data_, values.data(), indentation);
    } else if constexpr (!L::hasVectorParams) {
      std::vector<double> params;
      recoverParametersFromFits(t, params);
      trees_[t].print(data_, params.data(), indentation);
    } else {
      trees_[t].print(data_, paramsByTree_[t].data(), indentation, 0,
                      leaf_.numParams());
    }
  }

  /// The same for one saved tree, in the reference engine's saved format;
  /// function-valued leaves print their recorded mean values.
  void printSavedTree(size_t slot, size_t t, int indentation) const {
    const std::vector<FlatNode>& flat(savedTrees_[slot * options_.numTrees + t]);
    const std::uint64_t* masks = data_.hasWideCategorical
      ? savedTreeMasks_[slot * options_.numTrees + t].data() : nullptr;
    if constexpr (!L::hasVectorParams) {
      printFlatSubtree(data_, flat.data(), indentation, 0, nullptr, 0, 0,
                       masks);
    } else {
      const std::vector<double>& slopes(
        savedTreeParams_[slot * options_.numTrees + t]);
      printFlatSubtree(data_, flat.data(), indentation, 0, slopes.data(),
                       leaf_.numParams() - 1, 0, masks);
    }
  }

  /// Fits for raw column-major test rows from one saved sample's trees, on
  /// the original response scale; offsets are the caller's problem.
  void predictFromSavedSample(size_t slot, const double* x_test,
                              size_t numTestObservations, double* out) const {
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    const std::uint32_t* numCategories =
      data_.hasWideCategorical ? data_.numCuts.data() : nullptr;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      const std::vector<FlatNode>& flat(
        savedTrees_[slot * options_.numTrees + t]);
      const std::uint64_t* masks = data_.hasWideCategorical
        ? savedTreeMasks_[slot * options_.numTrees + t].data() : nullptr;
      for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
        addFlatPredictionsBelow(flat.data(), data_.types.data(), x_test,
                                numTestObservations, indices.data(), 0,
                                numTestObservations, out, numCategories,
                                masks);
      } else if constexpr (L::hasVectorParams) {
        const std::vector<double>& slopes(
          savedTreeParams_[slot * options_.numTrees + t]);
        addFlatLinearPredictionsBelow(
          flat.data(), data_.types.data(), x_test, numTestObservations,
          indices.data(), 0, numTestObservations, out,
          leaf_.covariateColumns().data(), leaf_.covariateMeans().data(),
          leaf_.covariateSds().data(), leaf_.numParams() - 1, slopes.data(),
          0, numCategories, masks);
      } else {
        const std::vector<double>& blocks(
          savedTreeParams_[slot * options_.numTrees + t]);
        // valid by construction (flatten wrote it) or by stateIsValid
        computeFunctionBlockOffsets(blocks.data(), blocks.size(),
                                    (flat.size() + 1) / 2,
                                    leaf_.numCovariates(), blockOffsets);
        addFlatFunctionPredictionsBelow(
          flat.data(), data_.types.data(), x_test, numTestObservations,
          indices.data(), 0, numTestObservations, out,
          leaf_.covariateColumns().data(), leaf_.covariateMeans().data(),
          leaf_.covariateSds().data(), leaf_.lengthscales().data(),
          leaf_.numCovariates(), blocks.data(), blockOffsets.data(), 0,
          numCategories, masks);
      }
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  /// The same from the live trees; scalar parameters recover from the
  /// current fits, vector ones read their persisted blocks.
  void predictFromCurrentTrees(const double* x_test,
                               size_t numTestObservations, double* out) {
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    std::vector<size_t> indices(numTestObservations);
    std::vector<double> params, slopes;
    std::vector<size_t> blockOffsets;
    std::vector<FlatNode> flat;
    std::vector<std::uint64_t> maskBuffer;
    const std::uint32_t* numCategories =
      data_.hasWideCategorical ? data_.numCuts.data() : nullptr;
    std::vector<std::uint64_t>* masks =
      data_.hasWideCategorical ? &maskBuffer : nullptr;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
        recoverParametersFromFits(t, params);
        trees_[t].flatten(data_, params.data(), flat, nullptr, 1, nullptr,
                          masks);
        addFlatPredictionsBelow(flat.data(), data_.types.data(), x_test,
                                numTestObservations, indices.data(), 0,
                                numTestObservations, out, numCategories,
                                maskBuffer.data());
      } else if constexpr (L::hasFunctionParams) {
        // no draw cache is fresh between runs: blocks recompute from the
        // persisted fits against the current covariates
        size_t n = data_.numObservations;
        functionLeafValues(trees_[t], treeFits_.data() + t * n, params);
        trees_[t].flatten(data_, params.data(), flat, nullptr, 1, nullptr,
                          masks);
        slopes.clear();
        for (int32_t i : trees_[t].bottomScratch)
          leaf_.appendLeafBlock(trees_[t], i, treeFits_.data() + t * n,
                                slopes);
        computeFunctionBlockOffsets(slopes.data(), slopes.size(),
                                    (flat.size() + 1) / 2,
                                    leaf_.numCovariates(), blockOffsets);
        addFlatFunctionPredictionsBelow(
          flat.data(), data_.types.data(), x_test, numTestObservations,
          indices.data(), 0, numTestObservations, out,
          leaf_.covariateColumns().data(), leaf_.covariateMeans().data(),
          leaf_.covariateSds().data(), leaf_.lengthscales().data(),
          leaf_.numCovariates(), slopes.data(), blockOffsets.data(), 0,
          numCategories, maskBuffer.data());
      } else {
        trees_[t].flatten(data_, paramsByTree_[t].data(), flat, nullptr,
                          leaf_.numParams(), &slopes, masks);
        addFlatLinearPredictionsBelow(
          flat.data(), data_.types.data(), x_test, numTestObservations,
          indices.data(), 0, numTestObservations, out,
          leaf_.covariateColumns().data(), leaf_.covariateMeans().data(),
          leaf_.covariateSds().data(), leaf_.numParams() - 1, slopes.data(),
          0, numCategories, maskBuffer.data());
      }
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  // Chain state serialization. getState captures everything the posterior
  // state comprises; stateIsValid checks a candidate against the store's
  // current cuts without mutating anything, so a multi-chain restore can be
  // all-or-none; setState trusts that check. Vector-parameter leaves carry
  // their slopes in treeParams/savedTreeParams alongside the flat trees.

  void getState(ChainStateData& state) {
    state.trees.resize(options_.numTrees);
    if (data_.hasWideCategorical) {
      state.treeMasks.resize(options_.numTrees);
      state.savedTreeMasks = savedTreeMasks_;
    } else {
      state.treeMasks.clear();
      state.savedTreeMasks.clear();
    }
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      std::vector<double> params;
      for (size_t t = 0; t < options_.numTrees; ++t) {
        recoverParametersFromFits(t, params);
        trees_[t].flatten(data_, params.data(), state.trees[t], nullptr, 1,
                          nullptr,
                          data_.hasWideCategorical ? &state.treeMasks[t]
                                                   : nullptr);
      }
      state.treeParams.clear();
      state.savedTreeParams.clear();
    } else if constexpr (L::hasVectorParams) {
      state.treeParams.resize(options_.numTrees);
      for (size_t t = 0; t < options_.numTrees; ++t)
        trees_[t].flatten(data_, paramsByTree_[t].data(), state.trees[t],
                          nullptr, leaf_.numParams(), &state.treeParams[t],
                          data_.hasWideCategorical ? &state.treeMasks[t]
                                                   : nullptr);
      state.savedTreeParams = savedTreeParams_;
    } else {
      // function-valued leaves: records carry reporting means, and each
      // live tree's parameters ARE its fits - one slab per tree in
      // observation order, restored by copy so continuation is bitwise
      size_t n = data_.numObservations;
      std::vector<double> values;
      state.treeParams.resize(options_.numTrees);
      for (size_t t = 0; t < options_.numTrees; ++t) {
        functionLeafValues(trees_[t], treeFits_.data() + t * n, values);
        trees_[t].flatten(data_, values.data(), state.trees[t], nullptr, 1,
                          nullptr,
                          data_.hasWideCategorical ? &state.treeMasks[t]
                                                   : nullptr);
        state.treeParams[t].assign(treeFits_.data() + t * n,
                                   treeFits_.data() + (t + 1) * n);
      }
      state.savedTreeParams = savedTreeParams_;
    }
    state.savedTrees = savedTrees_;
    state.totalFits = totalFits_;
    state.indices = indexBuffer_;
    state.sigma = sigma_;
    state.k = k_;
    response_->getScale(state.fitMin, state.fitMax);
    state.sigmaPriorScale = response_->sigmaPriorScaleInternal();
    if (response_->latents() != nullptr) {
      state.latents.assign(response_->latents(),
                           response_->latents() + data_.numObservations);
    } else {
      state.latents.clear();
    }
    if (response_->numGroupEffects() > 0) {
      state.groupEffects.assign(
        response_->groupEffects(),
        response_->groupEffects() + response_->numGroupEffects());
      state.groupTau = response_->groupTau();
    } else {
      state.groupEffects.clear();
      state.groupTau = 0.0;
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
    // mask channels pair with their flat trees when present; trees holding
    // wide rules without a channel fail the rebuild below
    if (!state.treeMasks.empty() &&
        state.treeMasks.size() != options_.numTrees)
      return false;
    if (!state.savedTreeMasks.empty() &&
        state.savedTreeMasks.size() != state.savedTrees.size())
      return false;
    if constexpr (L::hasVectorParams) {
      // the slope arrays must pair one-to-one with each flat tree's leaves
      // (a well-formed flat tree of m records has (m + 1) / 2 of them)
      size_t numSlopes = leaf_.numParams() - 1;
      if (state.treeParams.size() != options_.numTrees) return false;
      for (size_t t = 0; t < options_.numTrees; ++t)
        if (state.treeParams[t].size() !=
            ((state.trees[t].size() + 1) / 2) * numSlopes)
          return false;
      if (!state.savedTrees.empty()) {
        if (state.savedTreeParams.size() != state.savedTrees.size())
          return false;
        for (size_t s = 0; s < state.savedTrees.size(); ++s)
          if (state.savedTreeParams[s].size() !=
              ((state.savedTrees[s].size() + 1) / 2) * numSlopes)
            return false;
      }
    }
    if constexpr (L::hasFunctionParams) {
      // one fits slab per live tree; saved side channels must walk cleanly
      // against their trees' leaf counts
      if (state.treeParams.size() != options_.numTrees) return false;
      for (size_t t = 0; t < options_.numTrees; ++t)
        if (state.treeParams[t].size() != data_.numObservations) return false;
      if (!state.savedTrees.empty()) {
        if (state.savedTreeParams.size() != state.savedTrees.size())
          return false;
        std::vector<size_t> blockOffsets;
        for (size_t s = 0; s < state.savedTrees.size(); ++s)
          if (!computeFunctionBlockOffsets(
                state.savedTreeParams[s].data(),
                state.savedTreeParams[s].size(),
                (state.savedTrees[s].size() + 1) / 2, leaf_.numCovariates(),
                blockOffsets))
            return false;
      }
    }
    for (size_t s = 0; s < state.savedTrees.size(); ++s) {
      const std::vector<FlatNode>& saved(state.savedTrees[s]);
      const std::uint64_t* masks = state.savedTreeMasks.empty()
        ? nullptr : state.savedTreeMasks[s].data();
      size_t numMaskWords =
        state.savedTreeMasks.empty() ? 0 : state.savedTreeMasks[s].size();
      if (!flatTreeIsWellFormed(data_, saved.data(), saved.size(), masks,
                                numMaskWords))
        return false;
    }
    size_t n = data_.numObservations;
    if (!state.totalFits.empty() && state.totalFits.size() != n) return false;
    if (!state.indices.empty() &&
        state.indices.size() != n * options_.numTrees)
      return false;
    if (!state.latents.empty() &&
        (response_->latents() == nullptr || state.latents.size() != n))
      return false;
    // grouped states must carry a full effects vector for the chain's
    // groups; ungrouped states and chains both hold zero of them
    if (state.groupEffects.size() != response_->numGroupEffects())
      return false;
    if (options_.useDart && !state.dartProbabilities.empty() &&
        state.dartProbabilities.size() != data_.numPredictors)
      return false;
    if (state.fitMax < state.fitMin) return false;

    Tree scratch;
    std::vector<size_t> scratchIndices(n);
    std::vector<bool> seen(n);
    std::vector<double> params;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      scratch.initialize(scratchIndices.data(), n);
      const std::uint64_t* masks =
        state.treeMasks.empty() ? nullptr : state.treeMasks[t].data();
      size_t numMaskWords =
        state.treeMasks.empty() ? 0 : state.treeMasks[t].size();
      if (!scratch.buildFromFlat(data_, state.trees[t].data(),
                                 state.trees[t].size(), params, 1, nullptr,
                                 masks, numMaskWords))
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
    // internal-scale quantities below (tree parameters, fits, sigma) were
    // recorded under this transform; scale-free states leave creation's
    if (state.fitMax > state.fitMin)
      response_->restoreScale(state.fitMin, state.fitMax);
    if (state.sigmaPriorScale >= 0.0)
      response_->restoreSigmaPriorScaleInternal(state.sigmaPriorScale);
    misc_setVectorToConstant(totalFits_.data(), n, 0.0);
    std::vector<double> params;
    for (size_t t = 0; t < options_.numTrees; ++t) {
      trees_[t].initialize(indexBuffer_.data() + t * n, n);
      const std::uint64_t* masks =
        state.treeMasks.empty() ? nullptr : state.treeMasks[t].data();
      size_t numMaskWords =
        state.treeMasks.empty() ? 0 : state.treeMasks[t].size();
      if constexpr (!L::hasVectorParams) {
        if (!trees_[t].buildFromFlat(data_, state.trees[t].data(),
                                     state.trees[t].size(), params, 1,
                                     nullptr, masks, numMaskWords))
          return false;
      } else {
        if (!trees_[t].buildFromFlat(data_, state.trees[t].data(),
                                     state.trees[t].size(), params,
                                     leaf_.numParams(),
                                     state.treeParams[t].data(), masks,
                                     numMaskWords))
          return false;
      }
      if (state.indices.empty()) {
        trees_[t].repartitionSubtree(data_, 0);
      } else {
        std::memcpy(indexBuffer_.data() + t * n, state.indices.data() + t * n,
                    n * sizeof(size_t));
        if (!trees_[t].setPartitionsFromOrderedIndices(data_, 0)) return false;
      }
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
        setTreeFitsFromParameters(t, params);
      } else if constexpr (L::hasVectorParams) {
        paramsByTree_[t] = params;
        setTreeFitsFromParameterBlocks(t, params);
      } else {
        // the recorded slab IS the tree's parameters; copy restores bitwise
        std::memcpy(treeFits_.data() + t * n, state.treeParams[t].data(),
                    n * sizeof(double));
      }
      misc_addVectorsInPlace(treeFits_.data() + t * n, n, totalFits_.data());
    }
    if (!state.totalFits.empty())
      std::memcpy(totalFits_.data(), state.totalFits.data(),
                  n * sizeof(double));
    if (!state.savedTrees.empty()) {
      savedTrees_ = state.savedTrees;
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        savedTreeParams_ = state.savedTreeParams;
      if (data_.hasWideCategorical) {
        if (state.savedTreeMasks.empty())
          savedTreeMasks_.assign(savedTrees_.size(),
                                 std::vector<std::uint64_t>());
        else
          savedTreeMasks_ = state.savedTreeMasks;
      }
    }
    sigma_ = state.sigma;
    k_ = state.k;
    if (!state.latents.empty())
      response_->restoreLatents(state.latents.data());
    if (!state.groupEffects.empty())
      response_->restoreGroupEffects(state.groupEffects.data(),
                                     state.groupTau);
    if (options_.useDart && !state.dartProbabilities.empty()) {
      // the tree prior points at this vector's storage; overwrite in place
      std::memcpy(dart_.probabilities.data(), state.dartProbabilities.data(),
                  state.dartProbabilities.size() * sizeof(double));
      dart_.alpha = state.dartAlpha;
      dart_.setNumUpdatesSkipped(state.dartNumUpdatesSkipped);
    }
    // a serialized generator of a different kind (a single-chain state
    // riding R's stream restored into a dedicated-generator chain, say)
    // cannot be installed; the destination keeps its own stream, which
    // only forfeits bitwise continuation across generator kinds
    if (!state.rngState.empty() &&
        state.rngState.size() == ext_rng_getSerializedStateLength(rng_))
      ext_rng_readSerializedState(rng_, state.rngState.data());
    return true;
  }

  ext_rng* rng() const { return rng_; }
  const L& leaf() const { return leaf_; }
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
  /// before any re-route. Scalar leaves only: every vector-parameter caller
  /// is refused before reaching here (fits are no longer constant per leaf).
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

  /// Function-valued leaves' per-leaf reporting values, indexed by arena
  /// node id: the mean of the leaf's fits (fits indexed by observation),
  /// zero for empty leaves. FlatNode.value carries this for reporting; the
  /// saved side channel serves prediction.
  void functionLeafValues(Tree& tree, const double* fits,
                          std::vector<double>& valueByNode) {
    valueByNode.assign(tree.nodes.size(), 0.0);
    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      if (node.numObservations() == 0) continue;
      double total = 0.0;
      for (size_t m = node.begin; m < node.end; ++m)
        total += fits[tree.indices[m]];
      valueByNode[static_cast<size_t>(i)] =
        total / static_cast<double>(node.numObservations());
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

  /// The vector-parameter sibling of setTreeFitsFromParameters: each member
  /// observation's fit evaluates its leaf's block against the leaf model's
  /// current covariates.
  void setTreeFitsFromParameterBlocks(size_t t,
                                      const std::vector<double>& paramByNode) {
    Tree& tree(trees_[t]);
    double* treeFits = treeFits_.data() + t * data_.numObservations;
    size_t numParams = leaf_.numParams();

    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      const double* params =
        paramByNode.data() + static_cast<size_t>(i) * numParams;
      for (size_t m = node.begin; m < node.end; ++m) {
        size_t obs = tree.indices[m];
        treeFits[obs] = leaf_.fitForObservation(params, obs);
      }
    }
  }

  void sampleParametersAndSetFits(size_t t, double* fits, bool updateTestFits) {
    Tree& tree(trees_[t]);
    std::vector<int32_t>& bottoms(tree.bottomScratch);
    bottoms.clear();
    tree.fillBottom(0, bottoms);

    if constexpr (L::hasFunctionParams) {
      // function-valued draws land one value per member observation directly
      // in `fits` (the fits ARE the parameters); the per-node prediction cache
      // filled by the draws serves the routed test rows while this tree's
      // partitions are unchanged
      leaf_.beginTreeDraw(tree);
      for (int32_t i : bottoms) {
        FunctionLeafDrawStats stats = leaf_.drawFromPosteriorForNode(
          rng_, tree, treeY_.data(), response_->workingWeights(), k_,
          sigma_ * sigma_, i, fits);
        if (options_.updateK) {
          kSumSquaredParams_ += stats.sumSquaredParams;
          kNumLeaves_ += stats.numParams;
        }
      }

      if (updateTestFits && data_.numTestObservations > 0) {
        for (size_t i = 0; i < data_.numTestObservations; ++i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          currTestFits_[i] =
            leaf_.fitForTestObservationForNode(tree, leafIndex, i);
        }
      }
    } else if constexpr (!L::hasVectorParams) {
      paramByNode_.assign(tree.nodes.size(), 0.0);
      for (int32_t i : bottoms) {
        const Node& node(tree.at(i));
        double param = node.numObservations() == 0
          ? 0.0
          : leaf_.drawFromPosteriorForNode(rng_, tree, k_, sigma_ * sigma_, i);
        paramByNode_[static_cast<size_t>(i)] = param;

        if (options_.updateK) {
          kSumSquaredParams_ += param * param;
          kNumLeaves_ += 1.0;
        }

        if (node.parent == invalidNode) {
          misc_setVectorToConstant(fits, node.numObservations(), param);
        } else {
          misc_setIndexedVectorToConstant(fits,
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
    } else {
      size_t numParams = leaf_.numParams();
      // draws land in the tree's persistent blocks: they are the source of
      // truth for flatten, state, prediction, and the mutation flows
      std::vector<double>& treeParams(paramsByTree_[t]);
      treeParams.assign(tree.nodes.size() * numParams, 0.0);
      for (int32_t i : bottoms) {
        const Node& node(tree.at(i));
        double* params = treeParams.data() +
                         static_cast<size_t>(i) * numParams;
        // empty leaves keep the zero block without consuming draws, matching
        // the scalar path's zero parameter
        if (node.numObservations() > 0)
          leaf_.drawFromPosteriorForNode(rng_, tree, treeY_.data(),
                                         response_->workingWeights(), k_,
                                         sigma_ * sigma_, i, params);

        if (options_.updateK) {
          // every coordinate shares the scale / k prior sd, so the scaled-chi
          // posterior accumulates them all, the leaf count scaled to match
          for (size_t j = 0; j < numParams; ++j)
            kSumSquaredParams_ += params[j] * params[j];
          kNumLeaves_ += static_cast<double>(numParams);
        }

        for (size_t m = node.begin; m < node.end; ++m) {
          size_t obs = tree.indices[m];
          fits[obs] = leaf_.fitForObservation(params, obs);
        }
      }

      if (updateTestFits && data_.numTestObservations > 0) {
        for (size_t i = 0; i < data_.numTestObservations; ++i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          currTestFits_[i] = leaf_.fitForTestObservation(
            treeParams.data() + static_cast<size_t>(leafIndex) * numParams,
            i);
        }
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
      // original-scale convention, matching the classic engine and the
      // recorded test fits: any offset is part of the fit
      const double* offset = response_->offset();
      if (offset != nullptr) misc_addVectorsInPlace(offset, n, out);
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

    if (results.splitProbabilities != nullptr && options_.useDart) {
      double* out =
        results.splitProbabilities + sampleNum * data_.numPredictors;
      std::memcpy(out, dart_.probabilities.data(),
                  data_.numPredictors * sizeof(double));
    }

    // grouped channels de-scale like sigma: b is a deviation, so no shift
    if (results.tau != nullptr)
      results.tau[sampleNum] = response_->groupTau() * response_->sigmaScale();

    if (results.groupEffects != nullptr) {
      std::size_t numGroups = response_->numGroupEffects();
      double* out = results.groupEffects + sampleNum * numGroups;
      const double* effects = response_->groupEffects();
      double sigmaScale = response_->sigmaScale();
      for (std::size_t j = 0; j < numGroups; ++j)
        out[j] = effects[j] * sigmaScale;
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
  ResponseFamily family_ = ResponseFamily::gaussian;
  bool sigmaIsFixed_ = false;
  double sigma_ = 1.0;
  double k_ = 2.0;
  double kSumSquaredParams_ = 0.0;
  double kNumLeaves_ = 0.0;

  std::vector<Tree> trees_;
  std::vector<std::vector<FlatNode>> savedTrees_;  // slot-major, capacity x numTrees
  // vector-parameter leaves: each live tree's parameter blocks, arena-
  // indexed with stride numParams and kept consistent with the tree's
  // structure between sweeps (fits are no longer constant per leaf, so
  // parameters persist instead of being recovered); savedTreeParams_ holds
  // the slopes of each saved tree, parallel to savedTrees_. Both stay empty
  // for scalar leaf models.
  std::vector<std::vector<double>> paramsByTree_;
  std::vector<std::vector<double>> savedTreeParams_;
  // wide categorical columns (> 53 levels): each saved tree's flattened
  // mask words, parallel to savedTrees_; empty otherwise
  std::vector<std::vector<std::uint64_t>> savedTreeMasks_;
  size_t savedTreeCapacity_ = 0;
  size_t savedSlotBase_ = 0;
  std::vector<size_t> indexBuffer_;
  std::vector<double> treeFits_;
  std::vector<double> totalFits_, totalTestFits_;
  std::vector<double> treeY_, currTestFits_;
  std::vector<double> paramByNode_;
  MoveScratch scratch_;
};

}  // namespace bartcore

#endif  // BARTCORE_CHAIN_HPP
