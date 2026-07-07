#ifndef BARTCORE_CHAIN_HPP
#define BARTCORE_CHAIN_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <vector>

#include <external/random.h>
#include <misc/linearAlgebra.h>
#include <misc/thread.h>

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
/// form: value-encoded flattened trees, the saved-tree buffer, sigma
/// (original scale), k, response latents, DART state, and the serialized
/// rng. Restore rebuilds the rest canonically - partitions from the tree
/// structure and cut points, totalFits by summing the tree fits, the
/// variance prior by re-anchoring through the transform - so a restored
/// chain continues equivalently, not bitwise: the last ulp of the dropped
/// accumulation history is not reproduced.
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
  double sigma = 1.0;  // original response scale
  double k = 2.0;
  // the gaussian response transform at capture; max <= min marks scale-free
  double fitMin = 0.0, fitMax = 0.0;
  std::vector<double> latents;            // empty for gaussian
  // grouped samplers only, internal scale so restores are exact
  std::vector<double> groupEffects;
  double groupTau = 0.0;
  std::vector<double> dartProbabilities;  // empty when DART is off
  double dartAlpha = 1.0;
  size_t dartNumUpdatesSkipped = 0;
  std::vector<unsigned char> rngState;
};

/// One ensemble of the backfitting sampler: its trees, their fits and working
/// residual, the leaf model, the split selector, and the per-forest prior and
/// options. A Chain holds one-or-more (BCF combines a prognostic and a
/// treatment forest); everything a Forest omits - the rng, the response model,
/// sigma, and the shared store - is chain-level and passed in by the methods
/// that drive it.
template <IntegrableLeafModel L>
struct Forest {
  // per-forest options: tree count, tree-move probabilities, the k
  // hyperprior, and whether the split selector is DART
  size_t numTrees = 200;
  double birthOrDeathProbability = 0.5;
  double swapProbability = 0.1;
  double changeProbability = 0.4;
  double birthProbability = 0.5;
  bool updateK = false;
  ChiKHyperprior kHyperprior;
  bool useDart = false;

  L leaf;
  CGMTreePrior treePrior;
  DartPrior dart;
  std::vector<double> fixedSplitProbabilities;
  std::vector<std::uint32_t> splitCounts;

  // k is fixed unless updateK; the two accumulators gather the leaf sum of
  // squares and count over a sweep, feeding the k hyperprior draw
  double k = 2.0;
  double kSumSquaredParams = 0.0;
  double kNumLeaves = 0.0;

  std::vector<Tree> trees;
  std::vector<size_t> indexBuffer;
  std::vector<double> treeFits;
  std::vector<double> totalFits, totalTestFits;
  std::vector<double> treeY, currTestFits;
  std::vector<double> paramByNode;
  // vector-parameter leaves: each live tree's parameter blocks, arena-
  // indexed with stride numParams and kept consistent with the tree's
  // structure between sweeps (fits are no longer constant per leaf, so
  // parameters persist instead of being recovered); savedTreeParams holds
  // the slopes of each saved tree, parallel to savedTrees. Both stay empty
  // for scalar leaf models.
  std::vector<std::vector<double>> paramsByTree;
  MoveScratch scratch;

  // Saved-tree (keepTrees) storage: a circular buffer of capacity slots,
  // each one kept sample's forest in flattened form (slot-major, capacity x
  // numTrees). The slot base is set by the sampler before every run so
  // chains write consistent slots without sharing mutable state.
  std::vector<std::vector<FlatNode>> savedTrees;
  std::vector<std::vector<double>> savedTreeParams;
  // pooled categorical columns (> 63 levels): each saved tree's flattened
  // mask words, parallel to savedTrees; empty otherwise
  std::vector<std::vector<std::uint64_t>> savedTreeMasks;
  size_t savedTreeCapacity = 0;
  size_t savedSlotBase = 0;
};

/// Per-forest calibration for a BCF sampler: the prognostic (mu) and
/// treatment (tau) forests carry bcf's distinct tree counts and priors
/// (docs/design/bcf.md). Node scales are fixed - the adaptive magnitude
/// lives in the glue (a for mu, b0/b1 for tau), so no k hyperprior runs.
struct BCFForestSpec {
  std::size_t numTrees = 200;
  double base = 0.95, power = 2.0;
  double nodeScale = 0.5;
  double k = 1.0;
  double birthOrDeathProbability = 0.5, swapProbability = 0.1,
         changeProbability = 0.4, birthProbability = 0.5;
};

/// The two model specs plus the treatment vector a BCF chain is built from.
struct BCFSpec {
  BCFForestSpec mu, tau;
  const double* z = nullptr;    // borrowed 0/1 treatment indicator per obs
  double aPriorScale = 1.0;     // half-Cauchy median for the mu scalar a
  double bPriorVariance = 0.5;  // N(0, .) prior variance for b0, b1
};

/// The combining response's glue (docs/design/bcf.md): the prognostic scalar
/// a (half-Cauchy via the scale-mixture auxiliary aVariance) and the
/// treatment scales b0/b1, plus the sweep's per-forest scratch. y = a mu +
/// b_z tau + eps; a Chain holds one only in BCF mode.
struct BCFState {
  const double* z = nullptr;
  double a = 1.0, b0 = 0.0, b1 = 1.0;
  double aVariance = 1.0;
  double aPriorScale = 1.0;
  double bPriorVariance = 0.5;
  std::vector<double> combined, forestResponse, forestWeights;
};

/// One MCMC chain of the conjugate backfitting sampler: one-or-more forests
/// (trees, fits, per-forest prior), and its own response state and rng, over
/// a shared read-only ColumnStore. Data mutation is orchestrated one level up
/// (Sampler), which calls the tree-refresh methods here; a chain never writes
/// to the store.
///
/// The three leaf shapes (scalar, vector, function-valued parameters) each
/// discriminate in exactly one method per concern. Adding a leaf shape
/// touches: the constructor's leaf initialization, resizeTestStorage,
/// leafTracksNodeAverages, storeSavedTreeRecord, sampleParametersAndSetFits,
/// sampleNodeParametersFromPrior, recoverLeafParameters, setTreeFits,
/// rebuildFitsFromParameters, applyNewData, forceRefreshTrees,
/// initializeSavedTrees, flattenTree, printTree, printSavedTree,
/// addFlatPredictions, getState, setState, and stateIsValid.
template <IntegrableLeafModel L>
class Chain {
public:
  /// Scalar and function-valued leaves read per-node weighted means
  /// (function leaves because their over-cap nodes delegate to the constant
  /// leaf); vector leaves accumulate their own statistics.
  static constexpr bool leafTracksNodeAverages = !L::hasVectorParams;

  Chain(const ColumnStore& data, const double* y, const double* weights,
        const double* offset, ResponseFamily family, double sigmaEstimate,
        double sigmaDf, double sigmaRawScale, const SamplerOptions& options,
        ext_rng* rng)
    : options_(options), data_(data), weights_(weights), rng_(rng) {
    size_t numObservations = data.numObservations;
    options_.maxNumCutsPerVariable = nullptr;  // consumed by the store build
    options_.columnTypes = nullptr;

    forests_.emplace_back();
    Forest<L>& forest = forests_.back();
    forest.numTrees = options.numTrees;
    forest.birthOrDeathProbability = options.birthOrDeathProbability;
    forest.swapProbability = options.swapProbability;
    forest.changeProbability = options.changeProbability;
    forest.birthProbability = options.birthProbability;
    forest.updateK = options.updateK;
    forest.kHyperprior = options.kHyperprior;
    forest.useDart = options.useDart;

    if constexpr (L::hasVectorParams)
      forest.leaf.initialize(data, options.leafCovariateColumns,
                             options.numLeafCovariates);
    else if constexpr (L::hasFunctionParams)
      forest.leaf.initialize(data, options.leafCovariateColumns,
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
        std::make_unique<LogisticResponse>(y, offset, weights, numObservations);
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

    forest.leaf.scale =
      options.nodeScale / std::sqrt(static_cast<double>(forest.numTrees));
    forest.treePrior.base = options.base;
    forest.treePrior.power = options.power;
    family_ = family;
    sigmaIsFixed_ = family != ResponseFamily::gaussian || options.sigmaIsFixed;

    if (options.useDart) {
      forest.dart = options.dart;
      forest.dart.initialize(data.numPredictors);
      forest.treePrior.splitProbabilities = forest.dart.probabilities.data();
      forest.splitCounts.resize(data.numPredictors);
    } else if (options.splitProbabilities != nullptr) {
      forest.fixedSplitProbabilities.assign(
        options.splitProbabilities,
        options.splitProbabilities + data.numPredictors);
      forest.treePrior.splitProbabilities = forest.fixedSplitProbabilities.data();
    }

    sigma_ = response_->initialSigma();
    forest.k = options.k;

    forest.indexBuffer.resize(numObservations * forest.numTrees);
    forest.trees.resize(forest.numTrees);
    for (size_t t = 0; t < forest.numTrees; ++t)
      forest.trees[t].initialize(forest.indexBuffer.data() + t * numObservations,
                                 numObservations);

    forest.treeFits.assign(numObservations * forest.numTrees, 0.0);
    forest.totalFits.assign(numObservations, 0.0);
    forest.treeY.resize(numObservations);
    forest.paramByNode.clear();
    if constexpr (L::hasVectorParams)
      forest.paramsByTree.assign(forest.numTrees,
                                 std::vector<double>(forest.leaf.numParams(), 0.0));
    resizeTestStorage();
  }

  /// Two-forest BCF chain (docs/design/bcf.md): a prognostic forest mu
  /// (forest 0) and a treatment forest tau (forest 1) combined on a gaussian
  /// response as y = a mu + b_z tau + eps. Constant leaves only; both forests
  /// read the full store (the moderator subset arrives with data ownership).
  Chain(const ColumnStore& data, const double* y, const double* weights,
        const double* offset, double sigmaEstimate, double sigmaDf,
        double sigmaRawScale, const SamplerOptions& options,
        const BCFSpec& spec, ext_rng* rng)
    : options_(options), data_(data), weights_(weights), rng_(rng) {
    static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                  "BCF is a constant-leaf model");
    options_.maxNumCutsPerVariable = nullptr;
    options_.columnTypes = nullptr;
    response_ = std::make_unique<GaussianResponse>(
      y, offset, weights, data.numObservations, sigmaEstimate, sigmaDf,
      sigmaRawScale);
    family_ = ResponseFamily::gaussian;
    sigmaIsFixed_ = options.sigmaIsFixed;
    sigma_ = response_->initialSigma();

    buildBCFForest(spec.mu);
    buildBCFForest(spec.tau);

    bcf_ = std::make_unique<BCFState>();
    bcf_->z = spec.z;
    bcf_->aPriorScale = spec.aPriorScale;
    bcf_->bPriorVariance = spec.bPriorVariance;
    resizeTestStorage();
  }

  ~Chain() {
    if (testFitPool_ != nullptr) misc_mt_destroy(testFitPool_);
  }

  std::size_t numForests() const { return forests_.size(); }
  /// Re-forms b_{z_i} and both residuals on the next sweep; z is borrowed.
  void setTreatment(const double* z) { if (bcf_) bcf_->z = z; }
  /// BCF glue on the combining response; false for a non-BCF chain.
  bool bcfGlue(double& a, double& b0, double& b1) const {
    if (!bcf_) return false;
    a = bcf_->a; b0 = bcf_->b0; b1 = bcf_->b1;
    return true;
  }
  /// The forest's constant-leaf function values on the internal scale (mu for
  /// forest 0, tau for forest 1); numObservations doubles.
  void forestTotalFits(std::size_t f, double* out) const {
    std::memcpy(out, forests_[f].totalFits.data(),
                data_.numObservations * sizeof(double));
  }

  /// Between-run reconfiguration; the test-fit pool is rebuilt lazily to
  /// the new share of the budget on the next routing.
  void setNumThreads(size_t numThreads) { options_.numThreads = numThreads; }

  /// Called after the shared store's test data changes.
  void resizeTestStorage() {
    for (Forest<L>& forest : forests_) {
      forest.totalTestFits.assign(data_.numTestObservations, 0.0);
      forest.currTestFits.resize(data_.numTestObservations);
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        forest.leaf.rebuildTestCovariates(data_);
    }
  }

  /// One thinning-free run; results slots may be null to skip recording.
  /// Touches only chain state and the read-only store: safe to run chains
  /// concurrently as long as each has its own rng that never calls into R.
  /// progress, when non-null under verbose, receives one formatted line per
  /// printEvery kept iterations.
  /// Runs the chain, returning true if it stopped early because shouldCancel
  /// (polled once per sweep, called only on the thread that owns this chain)
  /// asked it to. The check touches no sampled state, so a run that is not
  /// cancelled is bitwise identical to one without it.
  bool run(size_t numBurnIn, size_t numSamples, Results& results,
           ProgressSink* progress = nullptr, size_t chainIndex = 0,
           const std::function<bool()>* shouldCancel = nullptr) {
    size_t n = data_.numObservations;
    size_t numThin = options_.numThin;
    double* y = response_->workingResponse();
    // per-iteration precisions: user weights for gaussian, the current
    // Polya-Gamma draws for logistic (refreshed with the latents)
    const double* weights = response_->workingWeights();

    size_t totalIterations = (numBurnIn + numSamples) * numThin;
    for (size_t iteration = 0; iteration < totalIterations; ++iteration) {
      // cooperative cancellation: stop between sweeps so no draw is left
      // half-applied. The results filled so far are discarded by the caller.
      if (shouldCancel != nullptr && (*shouldCancel)()) return true;

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

      for (size_t f = 0; f < forests_.size(); ++f) {
        Forest<L>& forest = forests_[f];
        // single-forest samplers backfit against the shared working response
        // and weights unchanged; a BCF forest sees its residual net of the
        // other forest's scaled contribution, divided by its own multiplier
        const double* forestY = y;
        const double* forestWeights = weights;
        if (bcf_) {
          formForestResponse(f, y, weights);
          forestY = bcf_->forestResponse.data();
          forestWeights = bcf_->forestWeights.data();
        }
        MoveContext ctx{data_,
                        forest.treePrior,
                        forest.birthOrDeathProbability,
                        forest.swapProbability,
                        forest.changeProbability,
                        forest.birthProbability,
                        forestWeights,
                        forest.k,
                        forest.scratch};

        forest.kSumSquaredParams = 0.0;
        forest.kNumLeaves = 0.0;

        if (record && data_.numTestObservations > 0)
          misc_setVectorToConstant(forest.totalTestFits.data(),
                                   data_.numTestObservations, 0.0);

        // treeY carries a running residual across the sweep: entering tree
        // t's body it holds y minus every other tree's current fits (new for
        // trees already drawn, old for the rest) - exactly the residual tree
        // t owns - so the draw can overwrite the slab in place. One fused
        // pass per tree retires the previous tree's new fits and admits this
        // tree's old ones; totalFits is stale until rebuilt after the loop.
        for (size_t t = 0; t < forest.numTrees; ++t) {
          double* treeFits = forest.treeFits.data() + t * n;

          if (t == 0) {
            const double* __restrict y_ = forestY;
            const double* __restrict total = forest.totalFits.data();
            const double* __restrict oldFits = treeFits;
            double* __restrict resid = forest.treeY.data();
            for (size_t i = 0; i < n; ++i)
              resid[i] = y_[i] - total[i] + oldFits[i];
          } else {
            const double* __restrict prevFits = treeFits - n;
            const double* __restrict oldFits = treeFits;
            double* __restrict resid = forest.treeY.data();
            for (size_t i = 0; i < n; ++i)
              resid[i] += oldFits[i] - prevFits[i];
          }

          // constant-leaf node means, recomputed against this sweep's residual
          if constexpr (leafTracksNodeAverages)
            forest.trees[t].setNodeAverages(forest.treeY.data(), forestWeights);

          bool stepTaken;
          StepType stepType;
          metropolisJumpForTree(ctx, forest.leaf, rng_, forest.trees[t],
                                forest.treeY.data(), sigma_, &stepTaken,
                                &stepType);
          // accepted changes and deaths strand pooled mask words; no rule
          // copies are live here, so this is a safe point to reclaim them
          if (data_.hasPooledCategorical)
            forest.trees[t].compactMaskPoolIfNeeded(data_);

          // the draw writes this tree's new fits straight into its slab
          sampleParametersAndSetFits(forest, t, treeFits, record);

          // flatten while the freshly drawn parameters are live
          if (record && forest.savedTreeCapacity > 0)
            storeSavedTreeRecord(
              forest, t,
              (forest.savedSlotBase + sampleNum) % forest.savedTreeCapacity,
              treeFits);

          if (record && data_.numTestObservations > 0)
            misc_addVectorsInPlace(forest.currTestFits.data(),
                                   data_.numTestObservations,
                                   forest.totalTestFits.data());
        }

        // rebuild the running total for the latent/sigma updates and
        // recording: the residual still includes the last tree's slab, whose
        // new fits retire here instead of in a pass of their own
        if (forest.numTrees > 0) {
          const size_t last = forest.numTrees - 1;
          const double* __restrict y_ = forestY;
          const double* __restrict resid = forest.treeY.data();
          const double* __restrict lastFits = forest.treeFits.data() + last * n;
          double* __restrict total = forest.totalFits.data();
          for (size_t i = 0; i < n; ++i)
            total[i] = y_[i] - resid[i] + lastFits[i];
        }
      }

      // a single forest reports its own fits; BCF the a mu + b_z tau blend
      const double* combined = combinedFits();
      response_->refreshLatents(rng_, combined, sigma_);
      y = response_->workingResponse();
      weights = response_->workingWeights();

      if (!sigmaIsFixed_)
        sigma_ = response_->drawSigma(rng_, combined, sigma_);

      if (bcf_) drawGlue(y, weights);

      for (Forest<L>& forest : forests_) {
        // a zero sum of squares under an infinite prior scale would make the
        // gamma rate zero and k infinite; with no information, keep k
        if (forest.updateK && forest.kSumSquaredParams > 0.0)
          forest.k = forest.kHyperprior.draw(rng_, forest.kSumSquaredParams,
                                             forest.kNumLeaves, forest.leaf.scale);

        if (forest.useDart) {
          std::memset(forest.splitCounts.data(), 0,
                      forest.splitCounts.size() * sizeof(std::uint32_t));
          for (size_t t = 0; t < forest.numTrees; ++t)
            forest.trees[t].countVariableUses(forest.splitCounts.data());
          forest.dart.update(rng_, forest.splitCounts.data());
        }
      }

      if (record) storeSample(results, sampleNum);
    }
    return false;
  }

  // Between-sample mutation; new-vector lifetimes are the caller's problem.
  void setOffset(const double* offset, bool updateScale) {
    response_->setOffset(offset, updateScale, &sigma_);
  }
  void setWeights(const double* weights) {
    weights_ = weights;
    response_->setWeights(weights);
  }
  void setResponse(const double* y, bool updateScale) {
    response_->setResponse(y, rng_, forests_[0].totalFits.data(), updateScale,
                           &sigma_);
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
    Forest<L>& forest = forests_[0];
    forest.treePrior.base = model.base;
    forest.treePrior.power = model.power;
    forest.birthOrDeathProbability = model.birthOrDeathProbability;
    forest.swapProbability = model.swapProbability;
    forest.changeProbability = model.changeProbability;
    forest.birthProbability = model.birthProbability;
    forest.leaf.scale =
      model.nodeScale / std::sqrt(static_cast<double>(forest.numTrees));

    forest.updateK = model.updateK;
    if (model.updateK) {
      forest.kHyperprior = model.kHyperprior;
    } else {
      forest.k = model.k;
    }

    if (family_ == ResponseFamily::gaussian) {
      sigmaIsFixed_ = model.sigmaIsFixed;
      if (model.sigmaIsFixed)
        setSigma(model.sigmaEstimate);
      else
        response_->setSigmaPrior(model.sigmaEstimate, model.sigmaDf,
                                 model.sigmaRawScale);
    }

    if (!forest.useDart) {
      if (model.splitProbabilities == nullptr) {
        forest.fixedSplitProbabilities.clear();
        forest.treePrior.splitProbabilities = nullptr;
      } else {
        forest.fixedSplitProbabilities.assign(
          model.splitProbabilities,
          model.splitProbabilities + data_.numPredictors);
        forest.treePrior.splitProbabilities =
          forest.fixedSplitProbabilities.data();
      }
    }
  }

  /// Sum of squared working residuals, descaled to the original response
  /// scale (binary families report on the latent scale).
  double sumOfSquaredResiduals() {
    double result = misc_computeSumOfSquaredResiduals(
      response_->workingResponse(), data_.numObservations,
      forests_[0].totalFits.data());
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
    for (Forest<L>& forest : forests_)
      for (size_t t = 0; t < forest.numTrees; ++t) {
        forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
        growSubtreeFromPrior(forest, forest.trees[t], 0, y, weights);
        paramByNode.assign(forest.trees[t].nodes.size(), 0.0);
        forest.trees[t].collapseEmptyNodes(data_, weights, paramByNode);
        // fresh structures carry zero parameter blocks until the next draw
        if constexpr (L::hasVectorParams)
          forest.paramsByTree[t].assign(
            forest.trees[t].nodes.size() * forest.leaf.numParams(), 0.0);
      }
  }

  /// Replace every leaf parameter with a draw from the node prior and
  /// rebuild the tree, total, and test fits to match.
  void sampleNodeParametersFromPrior() {
    size_t n = data_.numObservations;
    for (Forest<L>& forest : forests_) {
      misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);
      if (data_.numTestObservations > 0)
        misc_setVectorToConstant(forest.totalTestFits.data(),
                                 data_.numTestObservations, 0.0);

      for (size_t t = 0; t < forest.numTrees; ++t) {
        Tree& tree(forest.trees[t]);
        tree.bottomScratch.clear();
        tree.fillBottom(0, tree.bottomScratch);
        if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
          forest.paramByNode.assign(tree.nodes.size(), 0.0);
          for (int32_t i : tree.bottomScratch)
            forest.paramByNode[static_cast<size_t>(i)] =
              forest.leaf.drawFromPrior(rng_, forest.k);

          setTreeFitsFromParameters(forest, t, forest.paramByNode);
          misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                                 forest.totalFits.data());
          routeTestRows(data_.numTestObservations, [&](size_t i) {
            int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
            forest.totalTestFits[i] +=
              forest.paramByNode[static_cast<size_t>(leafIndex)];
          });
        } else if constexpr (L::hasFunctionParams) {
          // prior fits land directly in the tree's fit slab; the per-node
          // prediction cache serves the routed test rows
          double* treeFits = forest.treeFits.data() + t * n;
          forest.leaf.beginTreeDraw(tree);
          for (int32_t i : tree.bottomScratch)
            forest.leaf.drawFromPriorForNode(rng_, tree, forest.k, i, treeFits);
          misc_addVectorsInPlace(treeFits, n, forest.totalFits.data());
          routeTestRows(data_.numTestObservations, [&](size_t i) {
            int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
            forest.totalTestFits[i] +=
              forest.leaf.fitForTestObservationForNode(tree, leafIndex, i);
          });
        } else {
          size_t numParams = forest.leaf.numParams();
          std::vector<double>& treeParams(forest.paramsByTree[t]);
          treeParams.assign(tree.nodes.size() * numParams, 0.0);
          for (int32_t i : tree.bottomScratch)
            forest.leaf.drawFromPrior(rng_, forest.k,
                                      treeParams.data() +
                                        static_cast<size_t>(i) * numParams);

          setTreeFitsFromParameterBlocks(forest, t, treeParams);
          misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                                 forest.totalFits.data());
          routeTestRows(data_.numTestObservations, [&](size_t i) {
            int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
            forest.totalTestFits[i] += forest.leaf.fitForTestObservation(
              treeParams.data() + static_cast<size_t>(leafIndex) * numParams,
              i);
          });
        }
      }
    }
  }

  // Tree refreshes after the shared store changed, driven by the sampler.
  // The validate and rebuild phases are split so a transaction can check
  // every chain before any chain's fits are overwritten.

  using TreeParameters = std::vector<std::vector<double>>;

  /// Leaf parameters of tree t in transferable form: recovered from the
  /// fits for scalar leaves, copied from the persisted blocks for vector
  /// ones; function-valued leaves keep their per-observation fits in place,
  /// so nothing is recovered.
  void recoverLeafParameters(Forest<L>& forest, size_t t,
                             std::vector<double>& params) {
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams)
      recoverParametersFromFits(forest, t, params);
    else if constexpr (L::hasVectorParams)
      params = forest.paramsByTree[t];
    else
      params.clear();
  }

  /// Rewrites tree t's fit slab from recovered parameters, persisting them
  /// for vector leaves; not a function-leaf flow (their fits ARE the
  /// parameters and every caller handles them separately).
  void setTreeFits(Forest<L>& forest, size_t t,
                   const std::vector<double>& params) {
    static_assert(!L::hasFunctionParams);
    if constexpr (!L::hasVectorParams) {
      setTreeFitsFromParameters(forest, t, params);
    } else {
      forest.paramsByTree[t] = params;
      setTreeFitsFromParameterBlocks(forest, t, params);
    }
  }

  /// Recover leaf parameters (from fits for scalar leaves, from the
  /// persisted blocks for vector ones; function-valued leaves keep their
  /// per-observation fits in place, so nothing is recovered), re-route every
  /// tree against the store's current codes, and report whether all leaves
  /// stay occupied.
  bool revalidateTrees(TreeParameters& params) {
    Forest<L>& forest = forests_[0];
    params.resize(forest.numTrees);
    bool allValid = true;
    for (size_t t = 0; t < forest.numTrees && allValid; ++t) {
      recoverLeafParameters(forest, t, params[t]);
      forest.trees[t].repartitionSubtree(data_, 0);
      allValid = forest.trees[t].bottomNodesAreOccupied();
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
    Forest<L>& forest = forests_[0];
    dropStaleMissingDirections();
    if constexpr (L::hasFunctionParams) {
      (void) params;
      forest.leaf.regatherTrainingCovariates(data_);
    } else {
      // vector leaves read raw covariate values: pick up the installed ones
      if constexpr (L::hasVectorParams)
        forest.leaf.regatherTrainingCovariates(data_);
      size_t n = data_.numObservations;
      for (size_t t = 0; t < forest.numTrees; ++t) {
        double* treeFits = forest.treeFits.data() + t * n;
        misc_subtractVectorsInPlace(treeFits, n, forest.totalFits.data());
        setTreeFits(forest, t, params[t]);
        misc_addVectorsInPlace(treeFits, n, forest.totalFits.data());
      }
    }
  }

  /// Rollback re-route: restore partitions consistent with the store after
  /// the sampler has restored its old codes.
  void repartitionTrees() {
    for (Forest<L>& forest : forests_) {
      // the restored raw values re-gather to exactly the old covariates
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        forest.leaf.regatherTrainingCovariates(data_);
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].repartitionSubtree(data_, 0);
    }
  }

  /// First phase of a whole-data replacement: recover every tree's leaf
  /// parameters against the current fits and partitions, before the store or
  /// any per-observation storage moves.
  void recoverTreeParameters(TreeParameters& params) {
    Forest<L>& forest = forests_[0];
    params.resize(forest.numTrees);
    for (size_t t = 0; t < forest.numTrees; ++t)
      recoverLeafParameters(forest, t, params[t]);
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
    Forest<L>& forest = forests_[0];
    size_t n = data_.numObservations;
    bool numObservationsChanged =
      n * forest.numTrees != forest.indexBuffer.size();

    weights_ = weights;
    response_->setData(y, offset, weights, n, &sigma_);

    if (numObservationsChanged) {
      forest.indexBuffer.resize(n * forest.numTrees);
      forest.treeFits.resize(n * forest.numTrees);
      forest.totalFits.resize(n);
      forest.treeY.resize(n);
    }
    misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);

    // fresh standardization constants over the replacement data, like the
    // rebuilt cut grid; the persisted parameters carry over as-is, the same
    // approximate continuation the split remap embodies. Function-valued
    // parameters are per-observation fits over the OLD data, so they
    // cold-start at zero instead (the next sweep's draws replace them);
    // structures still carry through the split remap.
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      forest.leaf.reinitialize(data_);
    if constexpr (L::hasFunctionParams)
      for (size_t t = 0; t < forest.numTrees; ++t)
        params[t].assign(forest.trees[t].nodes.size(), 0.0);

    size_t paramStride = 1;
    if constexpr (L::hasVectorParams) paramStride = forest.leaf.numParams();

    dropStaleMissingDirections();
    for (size_t t = 0; t < forest.numTrees; ++t) {
      forest.trees[t].mapOldCutPointsOntoNew(data_, oldCutPoints, params[t],
                                             paramStride);
      if (numObservationsChanged)
        forest.trees[t].resetObservations(forest.indexBuffer.data() + t * n, n);
      forest.trees[t].repartitionSubtree(data_, 0);
      forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                         params[t], paramStride);
      if constexpr (!L::hasVectorParams) {
        setTreeFitsFromParameters(forest, t, params[t]);
      } else {
        forest.paramsByTree[t] = params[t];
        setTreeFitsFromParameterBlocks(forest, t, params[t]);
      }
      misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                             forest.totalFits.data());
    }

    resizeTestStorage();
  }

  /// After a data mutation re-quantizes the store, drop every tree's stale
  /// missing directions so the live masks stay within the reachable gauge.
  void dropStaleMissingDirections() {
    for (Forest<L>& forest : forests_)
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].dropStaleMissingDirections(data_);
  }

  /// Unconditional refresh: re-route and collapse any node an empty leaf
  /// leaves behind, merging leaf parameters into the collapsed node.
  /// Function-valued leaves keep their per-observation fits (they remain a
  /// coherent state under any partition) and collapse structure only.
  void forceRefreshTrees() {
    size_t n = data_.numObservations;
    dropStaleMissingDirections();

    for (Forest<L>& forest : forests_) {
      misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);

      if constexpr (L::hasFunctionParams) {
        forest.leaf.regatherTrainingCovariates(data_);
        std::vector<double> dummyParams;
        for (size_t t = 0; t < forest.numTrees; ++t) {
          forest.trees[t].repartitionSubtree(data_, 0);
          dummyParams.assign(forest.trees[t].nodes.size(), 0.0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             dummyParams);
          misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                                 forest.totalFits.data());
        }
      } else if constexpr (!L::hasVectorParams) {
        std::vector<double> paramByNode;
        for (size_t t = 0; t < forest.numTrees; ++t) {
          recoverParametersFromFits(forest, t, paramByNode);
          forest.trees[t].repartitionSubtree(data_, 0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             paramByNode);
          setTreeFitsFromParameters(forest, t, paramByNode);
          misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                                 forest.totalFits.data());
        }
      } else {
        forest.leaf.regatherTrainingCovariates(data_);
        size_t numParams = forest.leaf.numParams();
        for (size_t t = 0; t < forest.numTrees; ++t) {
          forest.trees[t].repartitionSubtree(data_, 0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             forest.paramsByTree[t], numParams);
          setTreeFitsFromParameterBlocks(forest, t, forest.paramsByTree[t]);
          misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                                 forest.totalFits.data());
        }
      }
    }
  }

  // Saved-tree (keepTrees) storage: a circular buffer of capacity slots,
  // each one kept sample's forest in flattened form. The slot base is set by
  // the sampler before every run so chains write consistent slots without
  // sharing mutable state.

  void initializeSavedTrees(size_t capacity) {
    for (Forest<L>& forest : forests_) {
      forest.savedTreeCapacity = capacity;
      forest.savedTrees.assign(capacity * forest.numTrees,
                               std::vector<FlatNode>(1));
      // the default slot is a single zero leaf; its slopes are zero too, and
      // a function-valued slot holds one zero-constant block
      if constexpr (L::hasVectorParams)
        forest.savedTreeParams.assign(
          capacity * forest.numTrees,
          std::vector<double>(forest.leaf.numParams() - 1, 0.0));
      else if constexpr (L::hasFunctionParams)
        forest.savedTreeParams.assign(capacity * forest.numTrees,
                                      std::vector<double>{0.0, 0.0});
      if (data_.hasPooledCategorical)
        forest.savedTreeMasks.assign(capacity * forest.numTrees,
                                     std::vector<std::uint64_t>());
    }
  }
  void setSavedSlotBase(size_t base) {
    for (Forest<L>& forest : forests_) forest.savedSlotBase = base;
  }
  size_t savedTreeCapacity() const { return forests_[0].savedTreeCapacity; }
  const std::vector<FlatNode>& savedTree(size_t slot, size_t t) const {
    const Forest<L>& forest = forests_[0];
    return forest.savedTrees[slot * forest.numTrees + t];
  }
  /// Slopes of one saved tree (vector-parameter leaves), parallel to
  /// savedTree's pre-order leaves.
  const std::vector<double>& savedTreeSlopes(size_t slot, size_t t) const {
    const Forest<L>& forest = forests_[0];
    return forest.savedTreeParams[slot * forest.numTrees + t];
  }
  /// Flattened mask words of one saved tree (wide categorical columns).
  const std::vector<std::uint64_t>& savedTreeMasks(size_t slot,
                                                   size_t t) const {
    const Forest<L>& forest = forests_[0];
    return forest.savedTreeMasks[slot * forest.numTrees + t];
  }

  /// Flatten live tree t into saved slot `slot` while its freshly drawn
  /// parameters are live; treeFits is the tree's slab. Function-valued
  /// leaves' records carry per-leaf mean fits for reporting, and their side
  /// channel the draw cache's alpha weights plus covariate rows - the exact
  /// values the recorded test fits used, so saved replays bit-match them.
  void storeSavedTreeRecord(Forest<L>& forest, size_t t, size_t slot,
                            const double* treeFits) {
    std::vector<std::uint64_t>* masks = data_.hasPooledCategorical
      ? &forest.savedTreeMasks[slot * forest.numTrees + t] : nullptr;
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      forest.trees[t].flatten(data_, forest.paramByNode.data(),
                              forest.savedTrees[slot * forest.numTrees + t],
                              nullptr, 1, nullptr, masks);
    } else if constexpr (L::hasVectorParams) {
      forest.trees[t].flatten(data_, forest.paramsByTree[t].data(),
                              forest.savedTrees[slot * forest.numTrees + t],
                              nullptr, forest.leaf.numParams(),
                              &forest.savedTreeParams[slot * forest.numTrees + t],
                              masks);
    } else {
      functionLeafValues(forest.trees[t], treeFits, forest.paramByNode);
      forest.trees[t].flatten(data_, forest.paramByNode.data(),
                              forest.savedTrees[slot * forest.numTrees + t],
                              nullptr, 1, nullptr, masks);
      std::vector<double>& blocks(
        forest.savedTreeParams[slot * forest.numTrees + t]);
      blocks.clear();
      for (int32_t i : forest.trees[t].bottomScratch)
        forest.leaf.appendLeafBlockFromCache(forest.trees[t], i, blocks);
    }
  }

  /// Flatten live tree t; counts receive the current partition sizes.
  /// Scalar leaves recover parameters from the fits; vector leaves emit
  /// their persisted blocks, the slopes (numParams - 1 per leaf, pre-order)
  /// into slopes when non-null.
  void flattenTree(size_t t, std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes = nullptr,
                   std::vector<std::uint64_t>* masks = nullptr) {
    Forest<L>& forest = forests_[0];
    if constexpr (L::hasFunctionParams) {
      // records carry per-leaf mean fits; slopes, when requested, receives
      // the saved-format side-channel blocks recomputed from the fits
      size_t n = data_.numObservations;
      std::vector<double> values;
      functionLeafValues(forest.trees[t], forest.treeFits.data() + t * n,
                         values);
      forest.trees[t].flatten(data_, values.data(), nodes, &counts, 1, nullptr,
                              masks);
      if (slopes != nullptr) {
        slopes->clear();
        for (int32_t i : forest.trees[t].bottomScratch)
          forest.leaf.appendLeafBlock(forest.trees[t], i,
                                      forest.treeFits.data() + t * n, *slopes);
      }
    } else if constexpr (!L::hasVectorParams) {
      if (slopes != nullptr) slopes->clear();
      std::vector<double> params;
      recoverParametersFromFits(forest, t, params);
      forest.trees[t].flatten(data_, params.data(), nodes, &counts, 1, nullptr,
                              masks);
    } else {
      forest.trees[t].flatten(data_, forest.paramsByTree[t].data(), nodes,
                              &counts, forest.leaf.numParams(), slopes, masks);
    }
  }

  /// Info dump of live tree t; function-valued leaves print their per-leaf
  /// mean fits.
  void printTree(size_t t, int indentation) {
    Forest<L>& forest = forests_[0];
    if constexpr (L::hasFunctionParams) {
      std::vector<double> values;
      functionLeafValues(forest.trees[t],
                         forest.treeFits.data() + t * data_.numObservations,
                         values);
      forest.trees[t].print(data_, values.data(), indentation);
    } else if constexpr (!L::hasVectorParams) {
      std::vector<double> params;
      recoverParametersFromFits(forest, t, params);
      forest.trees[t].print(data_, params.data(), indentation);
    } else {
      forest.trees[t].print(data_, forest.paramsByTree[t].data(), indentation, 0,
                            forest.leaf.numParams());
    }
  }

  /// The same for one saved tree, in the reference engine's saved format;
  /// function-valued leaves print their recorded mean values.
  void printSavedTree(size_t slot, size_t t, int indentation) const {
    const Forest<L>& forest = forests_[0];
    const std::vector<FlatNode>& flat(
      forest.savedTrees[slot * forest.numTrees + t]);
    const std::uint64_t* masks = data_.hasPooledCategorical
      ? forest.savedTreeMasks[slot * forest.numTrees + t].data() : nullptr;
    if constexpr (!L::hasVectorParams) {
      printFlatSubtree(data_, flat.data(), indentation, 0, nullptr, 0, 0,
                       masks);
    } else {
      const std::vector<double>& slopes(
        forest.savedTreeParams[slot * forest.numTrees + t]);
      printFlatSubtree(data_, flat.data(), indentation, 0, slopes.data(),
                       forest.leaf.numParams() - 1, 0, masks);
    }
  }

  /// Adds one flattened tree's predictions for raw column-major test rows
  /// to out, dispatching on the leaf shape's record format: plain leaf
  /// values, slope blocks, or function blocks (whose offsets are valid by
  /// construction or by stateIsValid). sideChannel is null for scalar
  /// leaves; indices and blockOffsets are caller scratch.
  void addFlatPredictions(const std::vector<FlatNode>& flat,
                          const std::vector<double>* sideChannel,
                          const std::uint64_t* masks, const double* x_test,
                          size_t numTestObservations,
                          std::vector<size_t>& indices,
                          std::vector<size_t>& blockOffsets,
                          double* out) const {
    const L& leaf = forests_[0].leaf;
    for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      addFlatPredictionsBelow(flat.data(), x_test, numTestObservations,
                              indices.data(), 0, numTestObservations, out,
                              masks);
    } else if constexpr (L::hasVectorParams) {
      addFlatLinearPredictionsBelow(
        flat.data(), x_test, numTestObservations,
        indices.data(), 0, numTestObservations, out,
        leaf.covariateColumns().data(), leaf.covariateMeans().data(),
        leaf.covariateSds().data(), leaf.numParams() - 1,
        sideChannel->data(), 0, masks);
    } else {
      computeFunctionBlockOffsets(sideChannel->data(), sideChannel->size(),
                                  (flat.size() + 1) / 2,
                                  leaf.numCovariates(), blockOffsets);
      addFlatFunctionPredictionsBelow(
        flat.data(), x_test, numTestObservations,
        indices.data(), 0, numTestObservations, out,
        leaf.covariateColumns().data(), leaf.covariateMeans().data(),
        leaf.covariateSds().data(), leaf.lengthscales().data(),
        leaf.numCovariates(), sideChannel->data(), blockOffsets.data(), 0,
        masks);
    }
  }

  /// Fits for raw column-major test rows from one saved sample's trees, on
  /// the original response scale; offsets are the caller's problem.
  void predictFromSavedSample(size_t slot, const double* x_test,
                              size_t numTestObservations, double* out) const {
    const Forest<L>& forest = forests_[0];
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    for (size_t t = 0; t < forest.numTrees; ++t) {
      const std::uint64_t* masks = data_.hasPooledCategorical
        ? forest.savedTreeMasks[slot * forest.numTrees + t].data() : nullptr;
      const std::vector<double>* sideChannel = forest.savedTreeParams.empty()
        ? nullptr : &forest.savedTreeParams[slot * forest.numTrees + t];
      addFlatPredictions(forest.savedTrees[slot * forest.numTrees + t],
                         sideChannel, masks, x_test, numTestObservations,
                         indices, blockOffsets, out);
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  /// The same from the live trees, flattened on the fly; function-valued
  /// leaves recompute their blocks from the persisted fits against the
  /// current covariates (no draw cache is fresh between runs).
  void predictFromCurrentTrees(const double* x_test,
                               size_t numTestObservations, double* out) {
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    std::vector<double> slopes;
    std::vector<std::uint32_t> counts;
    std::vector<FlatNode> flat;
    std::vector<std::uint64_t> maskBuffer;
    std::vector<std::uint64_t>* masks =
      data_.hasPooledCategorical ? &maskBuffer : nullptr;
    for (size_t t = 0; t < forests_[0].numTrees; ++t) {
      flattenTree(t, flat, counts, &slopes, masks);
      addFlatPredictions(flat, &slopes, maskBuffer.data(), x_test,
                         numTestObservations, indices, blockOffsets, out);
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
    Forest<L>& forest = forests_[0];
    state.trees.resize(forest.numTrees);
    if (data_.hasPooledCategorical) {
      state.treeMasks.resize(forest.numTrees);
      state.savedTreeMasks = forest.savedTreeMasks;
    } else {
      state.treeMasks.clear();
      state.savedTreeMasks.clear();
    }
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      std::vector<double> params;
      for (size_t t = 0; t < forest.numTrees; ++t) {
        recoverParametersFromFits(forest, t, params);
        forest.trees[t].flatten(data_, params.data(), state.trees[t], nullptr, 1,
                                nullptr,
                                data_.hasPooledCategorical ? &state.treeMasks[t]
                                                         : nullptr);
      }
      state.treeParams.clear();
      state.savedTreeParams.clear();
    } else if constexpr (L::hasVectorParams) {
      state.treeParams.resize(forest.numTrees);
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].flatten(data_, forest.paramsByTree[t].data(),
                                state.trees[t], nullptr, forest.leaf.numParams(),
                                &state.treeParams[t],
                                data_.hasPooledCategorical ? &state.treeMasks[t]
                                                         : nullptr);
      state.savedTreeParams = forest.savedTreeParams;
    } else {
      // function-valued leaves: records carry reporting means, and each
      // live tree's parameters ARE its fits - one slab per tree in
      // observation order, restored by copy
      size_t n = data_.numObservations;
      std::vector<double> values;
      state.treeParams.resize(forest.numTrees);
      for (size_t t = 0; t < forest.numTrees; ++t) {
        functionLeafValues(forest.trees[t], forest.treeFits.data() + t * n,
                           values);
        forest.trees[t].flatten(data_, values.data(), state.trees[t], nullptr, 1,
                                nullptr,
                                data_.hasPooledCategorical ? &state.treeMasks[t]
                                                         : nullptr);
        state.treeParams[t].assign(forest.treeFits.data() + t * n,
                                   forest.treeFits.data() + (t + 1) * n);
      }
      state.savedTreeParams = forest.savedTreeParams;
    }
    state.savedTrees = forest.savedTrees;
    state.sigma = sigma();
    state.k = forest.k;
    response_->getScale(state.fitMin, state.fitMax);
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
    if (forest.useDart) {
      state.dartProbabilities = forest.dart.probabilities;
      state.dartAlpha = forest.dart.alpha;
      state.dartNumUpdatesSkipped = forest.dart.numUpdatesSkipped();
    } else {
      state.dartProbabilities.clear();
    }
    state.rngState.resize(ext_rng_getSerializedStateLength(rng_));
    if (!state.rngState.empty())
      ext_rng_writeSerializedState(rng_, state.rngState.data());
  }

  bool stateIsValid(const ChainStateData& state) const {
    const Forest<L>& forest = forests_[0];
    if (state.trees.size() != forest.numTrees) return false;
    if (!state.savedTrees.empty() &&
        state.savedTrees.size() != forest.savedTrees.size())
      return false;
    // mask channels pair with their flat trees when present; trees holding
    // wide rules without a channel fail the rebuild below
    if (!state.treeMasks.empty() &&
        state.treeMasks.size() != forest.numTrees)
      return false;
    if (!state.savedTreeMasks.empty() &&
        state.savedTreeMasks.size() != state.savedTrees.size())
      return false;
    if constexpr (L::hasVectorParams) {
      // the slope arrays must pair one-to-one with each flat tree's leaves
      // (a well-formed flat tree of m records has (m + 1) / 2 of them)
      size_t numSlopes = forest.leaf.numParams() - 1;
      if (state.treeParams.size() != forest.numTrees) return false;
      for (size_t t = 0; t < forest.numTrees; ++t)
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
      if (state.treeParams.size() != forest.numTrees) return false;
      for (size_t t = 0; t < forest.numTrees; ++t)
        if (state.treeParams[t].size() != data_.numObservations) return false;
      if (!state.savedTrees.empty()) {
        if (state.savedTreeParams.size() != state.savedTrees.size())
          return false;
        std::vector<size_t> blockOffsets;
        for (size_t s = 0; s < state.savedTrees.size(); ++s)
          if (!computeFunctionBlockOffsets(
                state.savedTreeParams[s].data(),
                state.savedTreeParams[s].size(),
                (state.savedTrees[s].size() + 1) / 2, forest.leaf.numCovariates(),
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
    if (!state.latents.empty() &&
        (response_->latents() == nullptr || state.latents.size() != n))
      return false;
    // grouped states must carry a full effects vector for the chain's
    // groups; ungrouped states and chains both hold zero of them
    if (state.groupEffects.size() != response_->numGroupEffects())
      return false;
    if (forest.useDart && !state.dartProbabilities.empty() &&
        state.dartProbabilities.size() != data_.numPredictors)
      return false;
    if (state.fitMax < state.fitMin) return false;

    Tree scratch;
    std::vector<size_t> scratchIndices(n);
    std::vector<double> params;
    for (size_t t = 0; t < forest.numTrees; ++t) {
      scratch.initialize(scratchIndices.data(), n);
      const std::uint64_t* masks =
        state.treeMasks.empty() ? nullptr : state.treeMasks[t].data();
      size_t numMaskWords =
        state.treeMasks.empty() ? 0 : state.treeMasks[t].size();
      if (!scratch.buildFromFlat(data_, state.trees[t].data(),
                                 state.trees[t].size(), params, 1, nullptr,
                                 masks, numMaskWords))
        return false;
      scratch.repartitionSubtree(data_, 0);
      if (!scratch.bottomNodesAreOccupied()) return false;
    }
    return true;
  }

  /// Installs a state stateIsValid accepted; false only on the invariant
  /// violation of a validated tree failing to rebuild.
  bool setState(const ChainStateData& state) {
    Forest<L>& forest = forests_[0];
    size_t n = data_.numObservations;
    // the internal-scale tree parameters and fits below were recorded under
    // this transform; scale-free states leave creation's. restoreScale
    // re-anchors the variance prior through it.
    if (state.fitMax > state.fitMin)
      response_->restoreScale(state.fitMin, state.fitMax);
    misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);
    std::vector<double> params;
    for (size_t t = 0; t < forest.numTrees; ++t) {
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
      const std::uint64_t* masks =
        state.treeMasks.empty() ? nullptr : state.treeMasks[t].data();
      size_t numMaskWords =
        state.treeMasks.empty() ? 0 : state.treeMasks[t].size();
      if constexpr (!L::hasVectorParams) {
        if (!forest.trees[t].buildFromFlat(data_, state.trees[t].data(),
                                           state.trees[t].size(), params, 1,
                                           nullptr, masks, numMaskWords))
          return false;
      } else {
        if (!forest.trees[t].buildFromFlat(data_, state.trees[t].data(),
                                           state.trees[t].size(), params,
                                           forest.leaf.numParams(),
                                           state.treeParams[t].data(), masks,
                                           numMaskWords))
          return false;
      }
      forest.trees[t].repartitionSubtree(data_, 0);
      if constexpr (L::hasFunctionParams) {
        // the recorded slab IS the tree's parameters; copy restores bitwise
        std::memcpy(forest.treeFits.data() + t * n, state.treeParams[t].data(),
                    n * sizeof(double));
      } else {
        setTreeFits(forest, t, params);
      }
      misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                             forest.totalFits.data());
    }
    if (!state.savedTrees.empty()) {
      forest.savedTrees = state.savedTrees;
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        forest.savedTreeParams = state.savedTreeParams;
      if (data_.hasPooledCategorical) {
        if (state.savedTreeMasks.empty())
          forest.savedTreeMasks.assign(forest.savedTrees.size(),
                                       std::vector<std::uint64_t>());
        else
          forest.savedTreeMasks = state.savedTreeMasks;
      }
    }
    setSigma(state.sigma);
    forest.k = state.k;
    if (!state.latents.empty())
      response_->restoreLatents(state.latents.data());
    if (!state.groupEffects.empty())
      response_->restoreGroupEffects(state.groupEffects.data(),
                                     state.groupTau);
    if (forest.useDart && !state.dartProbabilities.empty()) {
      // the tree prior points at this vector's storage; overwrite in place
      std::memcpy(forest.dart.probabilities.data(),
                  state.dartProbabilities.data(),
                  state.dartProbabilities.size() * sizeof(double));
      forest.dart.alpha = state.dartAlpha;
      forest.dart.setNumUpdatesSkipped(state.dartNumUpdatesSkipped);
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
  const L& leaf() const { return forests_[0].leaf; }
  /// On the original response scale, symmetric with setSigma.
  double sigma() const { return sigma_ * response_->sigmaScale(); }
  double k() const { return forests_[0].k; }
  size_t numTrees() const { return forests_[0].numTrees; }
  const Tree& tree(size_t t) const { return forests_[0].trees[t]; }
  const std::vector<double>& treeFits() const { return forests_[0].treeFits; }
  const std::vector<double>& totalFits() const { return forests_[0].totalFits; }

private:
  template <typename F> struct TestFitRange { size_t begin, end; F* fn; };
  template <typename F> static void runTestFitRange(void* data) {
    TestFitRange<F>* r = static_cast<TestFitRange<F>*>(data);
    for (size_t i = r->begin; i < r->end; ++i) (*r->fn)(i);
  }

  /// Apply fn to each test row. Routing draws no rng and each row writes its
  /// own output slot, so splitting the range across this chain's share of the
  /// thread budget yields byte-identical results at any thread count. Serial
  /// below the cutoff, where dispatch outweighs the routing it saves.
  template <typename F>
  void routeTestRows(size_t numTest, F fn) {
    size_t chains = options_.numChains > 0 ? options_.numChains : 1;
    size_t budget = options_.numThreads / chains;
    if (budget >= 2 && numTest >= testFitParallelCutoff) {
      if (testFitPool_ == nullptr ||
          misc_mt_getNumThreads(testFitPool_) != budget) {
        if (testFitPool_ != nullptr) misc_mt_destroy(testFitPool_);
        misc_mt_create(&testFitPool_, budget);
      }
      size_t numThreads, perThread, offByOne;
      misc_mt_getNumThreadsForJob(testFitPool_, numTest,
                                  testFitParallelCutoff / 2, &numThreads,
                                  &perThread, &offByOne);
      if (numThreads > 1) {
        std::vector<TestFitRange<F>> ranges(numThreads);
        std::vector<void*> ptrs(numThreads);
        for (size_t w = 0, start = 0; w < numThreads; ++w) {
          size_t count = perThread - (w < offByOne ? 0 : 1);
          ranges[w] = TestFitRange<F>{start, start + count, &fn};
          ptrs[w] = &ranges[w];
          start += count;
        }
        misc_mt_runTasks(testFitPool_, &runTestFitRange<F>, ptrs.data(),
                         numThreads);
        return;
      }
    }
    for (size_t i = 0; i < numTest; ++i) fn(i);
  }

  /// The reference engine's recursion: growth is Bernoulli in the
  /// depth-decayed prior probability, rules come from the prior, and empty
  /// children keep growing (availability is rule-based) until the caller
  /// collapses them.
  void growSubtreeFromPrior(Forest<L>& forest, Tree& tree, int32_t nodeIndex,
                            const double* y, const double* weights) {
    double growthProbability =
      forest.treePrior.growthProbability(tree, data_, nodeIndex);
    if (growthProbability <= 0.0 ||
        ext_rng_simulateBernoulli(rng_, growthProbability) == 0)
      return;

    Rule rule =
      forest.treePrior.drawRuleAndVariable(tree, data_, rng_, nodeIndex);
    tree.birth(data_, nodeIndex, rule, y, weights);
    int32_t leftChild = tree.at(nodeIndex).leftChild;
    growSubtreeFromPrior(forest, tree, leftChild, y, weights);
    growSubtreeFromPrior(forest, tree, leftChild + 1, y, weights);
  }

  /// Leaf parameters recovered from a tree's fits, indexed by arena node id;
  /// fits are constant within a leaf, so any member observation's fit is the
  /// parameter. Must run against partitions consistent with the fits, i.e.
  /// before any re-route. Scalar leaves only: every vector-parameter caller
  /// is refused before reaching here (fits are no longer constant per leaf).
  void recoverParametersFromFits(Forest<L>& forest, size_t t,
                                 std::vector<double>& paramByNode) {
    Tree& tree(forest.trees[t]);
    const double* treeFits =
      forest.treeFits.data() + t * data_.numObservations;

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

  void setTreeFitsFromParameters(Forest<L>& forest, size_t t,
                                 const std::vector<double>& paramByNode) {
    Tree& tree(forest.trees[t]);
    double* treeFits = forest.treeFits.data() + t * data_.numObservations;

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
  void setTreeFitsFromParameterBlocks(Forest<L>& forest, size_t t,
                                      const std::vector<double>& paramByNode) {
    Tree& tree(forest.trees[t]);
    double* treeFits = forest.treeFits.data() + t * data_.numObservations;
    size_t numParams = forest.leaf.numParams();

    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      const double* params =
        paramByNode.data() + static_cast<size_t>(i) * numParams;
      for (size_t m = node.begin; m < node.end; ++m) {
        size_t obs = tree.indices[m];
        treeFits[obs] = forest.leaf.fitForObservation(params, obs);
      }
    }
  }

  void sampleParametersAndSetFits(Forest<L>& forest, size_t t, double* fits,
                                  bool updateTestFits) {
    Tree& tree(forest.trees[t]);
    std::vector<int32_t>& bottoms(tree.bottomScratch);
    bottoms.clear();
    tree.fillBottom(0, bottoms);

    if constexpr (L::hasFunctionParams) {
      // function-valued draws land one value per member observation directly
      // in `fits` (the fits ARE the parameters); the per-node prediction cache
      // filled by the draws serves the routed test rows while this tree's
      // partitions are unchanged
      forest.leaf.beginTreeDraw(tree);
      for (int32_t i : bottoms) {
        FunctionLeafDrawStats stats = forest.leaf.drawFromPosteriorForNode(
          rng_, tree, forest.treeY.data(), response_->workingWeights(),
          forest.k, sigma_ * sigma_, i, fits);
        if (forest.updateK) {
          forest.kSumSquaredParams += stats.sumSquaredParams;
          forest.kNumLeaves += stats.numParams;
        }
      }

      if (updateTestFits && data_.numTestObservations > 0)
        routeTestRows(data_.numTestObservations, [&](size_t i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          forest.currTestFits[i] =
            forest.leaf.fitForTestObservationForNode(tree, leafIndex, i);
        });
    } else if constexpr (!L::hasVectorParams) {
      forest.paramByNode.assign(tree.nodes.size(), 0.0);
      for (int32_t i : bottoms) {
        const Node& node(tree.at(i));
        double param = node.numObservations() == 0
          ? 0.0
          : forest.leaf.drawFromPosteriorForNode(rng_, tree, forest.k,
                                                 sigma_ * sigma_, i);
        forest.paramByNode[static_cast<size_t>(i)] = param;

        if (forest.updateK) {
          forest.kSumSquaredParams += param * param;
          forest.kNumLeaves += 1.0;
        }

        if (node.parent == invalidNode) {
          misc_setVectorToConstant(fits, node.numObservations(), param);
        } else {
          misc_setIndexedVectorToConstant(fits,
                                          tree.indices + node.begin,
                                          node.numObservations(), param);
        }
      }

      if (updateTestFits && data_.numTestObservations > 0)
        routeTestRows(data_.numTestObservations, [&](size_t i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          forest.currTestFits[i] =
            forest.paramByNode[static_cast<size_t>(leafIndex)];
        });
    } else {
      size_t numParams = forest.leaf.numParams();
      // draws land in the tree's persistent blocks: they are the source of
      // truth for flatten, state, prediction, and the mutation flows
      std::vector<double>& treeParams(forest.paramsByTree[t]);
      treeParams.assign(tree.nodes.size() * numParams, 0.0);
      for (int32_t i : bottoms) {
        const Node& node(tree.at(i));
        double* params = treeParams.data() +
                         static_cast<size_t>(i) * numParams;
        // empty leaves keep the zero block without consuming draws, matching
        // the scalar path's zero parameter
        if (node.numObservations() > 0)
          forest.leaf.drawFromPosteriorForNode(rng_, tree, forest.treeY.data(),
                                               response_->workingWeights(),
                                               forest.k, sigma_ * sigma_, i,
                                               params);

        if (forest.updateK) {
          // every coordinate shares the scale / k prior sd, so the scaled-chi
          // posterior accumulates them all, the leaf count scaled to match
          for (size_t j = 0; j < numParams; ++j)
            forest.kSumSquaredParams += params[j] * params[j];
          forest.kNumLeaves += static_cast<double>(numParams);
        }

        for (size_t m = node.begin; m < node.end; ++m) {
          size_t obs = tree.indices[m];
          fits[obs] = forest.leaf.fitForObservation(params, obs);
        }
      }

      if (updateTestFits && data_.numTestObservations > 0)
        routeTestRows(data_.numTestObservations, [&](size_t i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, data_.testRow(i));
          forest.currTestFits[i] = forest.leaf.fitForTestObservation(
            treeParams.data() + static_cast<size_t>(leafIndex) * numParams,
            i);
        });
    }
  }

  /// A BCF forest, built self-contained so the single-forest constructor
  /// (the bitwise-gated path) is untouched: constant leaf, fixed k, no DART.
  void buildBCFForest(const BCFForestSpec& spec) {
    std::size_t n = data_.numObservations;
    forests_.emplace_back();
    Forest<L>& forest = forests_.back();
    forest.numTrees = spec.numTrees;
    forest.birthOrDeathProbability = spec.birthOrDeathProbability;
    forest.swapProbability = spec.swapProbability;
    forest.changeProbability = spec.changeProbability;
    forest.birthProbability = spec.birthProbability;
    forest.updateK = false;
    forest.useDart = false;
    forest.k = spec.k;
    forest.leaf.scale =
      spec.nodeScale / std::sqrt(static_cast<double>(spec.numTrees));
    forest.treePrior.base = spec.base;
    forest.treePrior.power = spec.power;
    forest.indexBuffer.resize(n * spec.numTrees);
    forest.trees.resize(spec.numTrees);
    for (std::size_t t = 0; t < spec.numTrees; ++t)
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
    forest.treeFits.assign(n * spec.numTrees, 0.0);
    forest.totalFits.assign(n, 0.0);
    forest.treeY.resize(n);
  }

  double forestMultiplier(std::size_t f, std::size_t i) const {
    if (f == 0) return bcf_->a;
    return bcf_->z[i] != 0.0 ? bcf_->b1 : bcf_->b0;
  }

  /// Effective response and precision forest f's constant-leaf draws see so
  /// that m_f * f_f explains the residual (y minus the other forests' scaled
  /// contributions): response r_i / m_f, weight w_i m_f^2, which reproduce
  /// the leaf's node sums without touching the leaf math. |m_f| is floored to
  /// keep the division finite in the pathological near-zero-scale case.
  void formForestResponse(std::size_t f, const double* y, const double* w) {
    std::size_t n = data_.numObservations;
    bcf_->forestResponse.resize(n);
    bcf_->forestWeights.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      double resid = y[i];
      for (std::size_t g = 0; g < forests_.size(); ++g)
        if (g != f) resid -= forestMultiplier(g, i) * forests_[g].totalFits[i];
      double m = forestMultiplier(f, i);
      if (std::fabs(m) < 1.0e-9) m = m < 0.0 ? -1.0e-9 : 1.0e-9;
      bcf_->forestResponse[i] = resid / m;
      bcf_->forestWeights[i] = (w == nullptr ? 1.0 : w[i]) * m * m;
    }
  }

  const double* combinedFits() {
    if (!bcf_) return forests_[0].totalFits.data();
    std::size_t n = data_.numObservations;
    bcf_->combined.resize(n);
    const double* mu = forests_[0].totalFits.data();
    const double* tau = forests_[1].totalFits.data();
    for (std::size_t i = 0; i < n; ++i)
      bcf_->combined[i] = bcf_->a * mu[i] +
        (bcf_->z[i] != 0.0 ? bcf_->b1 : bcf_->b0) * tau[i];
    return bcf_->combined.data();
  }

  /// The glue's Gaussian full conditionals (docs/design/bcf.md): a as the mu
  /// coefficient (prior N(0, aVariance), whose half-Cauchy scale mixture is
  /// refreshed after via an inverse-gamma auxiliary), b0/b1 as the tau
  /// coefficients over control/treated (prior N(0, bPriorVariance)).
  void drawGlue(const double* y, const double* w) {
    std::size_t n = data_.numObservations;
    const double* mu = forests_[0].totalFits.data();
    const double* tau = forests_[1].totalFits.data();
    double invSigmaSq = 1.0 / (sigma_ * sigma_);

    double aPrec = 1.0 / bcf_->aVariance, aNum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double wi = w == nullptr ? 1.0 : w[i];
      double bz = bcf_->z[i] != 0.0 ? bcf_->b1 : bcf_->b0;
      double r = y[i] - bz * tau[i];
      aPrec += wi * mu[i] * mu[i] * invSigmaSq;
      aNum += wi * mu[i] * r * invSigmaSq;
    }
    bcf_->a =
      aNum / aPrec + ext_rng_simulateStandardNormal(rng_) / std::sqrt(aPrec);

    // t_1 scale mixture: aVariance ~ IG(1/2, scale^2/2) mixes N(0, aVariance)
    // to Cauchy(0, scale), so the conditional's rate carries scale^2, not its
    // inverse
    double rate = 0.5 * bcf_->a * bcf_->a +
                  0.5 * bcf_->aPriorScale * bcf_->aPriorScale;
    bcf_->aVariance = 1.0 / ext_rng_simulateGamma(rng_, 1.0, 1.0 / rate);

    double bPrec = 1.0 / bcf_->bPriorVariance;
    double p0 = bPrec, n0 = 0.0, p1 = bPrec, n1 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double wi = w == nullptr ? 1.0 : w[i];
      double r = y[i] - bcf_->a * mu[i];
      double prec = wi * tau[i] * tau[i] * invSigmaSq;
      double num = wi * tau[i] * r * invSigmaSq;
      if (bcf_->z[i] != 0.0) { p1 += prec; n1 += num; }
      else { p0 += prec; n0 += num; }
    }
    bcf_->b0 = n0 / p0 + ext_rng_simulateStandardNormal(rng_) / std::sqrt(p0);
    bcf_->b1 = n1 / p1 + ext_rng_simulateStandardNormal(rng_) / std::sqrt(p1);
  }

  void storeSample(Results& results, size_t sampleNum) {
    Forest<L>& forest = forests_[0];
    size_t n = data_.numObservations;
    double scale = response_->fitScale();
    double shift = response_->fitShift();

    if (results.sigma != nullptr)
      results.sigma[sampleNum] = sigma_ * response_->sigmaScale();

    if (results.k != nullptr) results.k[sampleNum] = forest.k;

    if (results.trainingFits != nullptr) {
      double* out = results.trainingFits + sampleNum * n;
      if (bcf_) {
        const double* mu = forests_[0].totalFits.data();
        const double* tau = forests_[1].totalFits.data();
        for (size_t i = 0; i < n; ++i)
          out[i] = scale * (bcf_->a * mu[i] +
                            (bcf_->z[i] != 0.0 ? bcf_->b1 : bcf_->b0) * tau[i]) +
                   shift;
      } else {
        for (size_t i = 0; i < n; ++i)
          out[i] = scale * forest.totalFits[i] + shift;
      }
      // original-scale convention, matching the classic engine and the
      // recorded test fits: any offset is part of the fit
      const double* offset = response_->offset();
      if (offset != nullptr) misc_addVectorsInPlace(offset, n, out);
    }

    if (results.testFits != nullptr && data_.numTestObservations > 0) {
      double* out = results.testFits + sampleNum * data_.numTestObservations;
      for (size_t i = 0; i < data_.numTestObservations; ++i)
        out[i] = scale * forest.totalTestFits[i] + shift;
      if (data_.testOffset != nullptr)
        misc_addVectorsInPlace(data_.testOffset, data_.numTestObservations,
                               out);
    }

    if (results.variableCounts != nullptr) {
      std::uint32_t* out =
        results.variableCounts + sampleNum * data_.numPredictors;
      std::memset(out, 0, data_.numPredictors * sizeof(std::uint32_t));
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].countVariableUses(out);
    }

    if (results.splitProbabilities != nullptr && forest.useDart) {
      double* out =
        results.splitProbabilities + sampleNum * data_.numPredictors;
      std::memcpy(out, forest.dart.probabilities.data(),
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

  std::unique_ptr<ResponseModel> response_;
  ResponseFamily family_ = ResponseFamily::gaussian;
  bool sigmaIsFixed_ = false;
  double sigma_ = 1.0;

  // one-or-more forests (size 1 for every non-BCF sampler); the sweep,
  // mutation fan-out, and prediction loop over them
  std::vector<Forest<L>> forests_;

  // the BCF combining response's glue and sweep scratch; null off BCF, and
  // the whole two-forest sweep collapses to the single-forest path when so
  std::unique_ptr<BCFState> bcf_;

  // Persistent pool for parallel test-fit routing, sized to this chain's
  // share of the thread budget; created lazily, never below the cutoff. The
  // forests borrow it through routeTestRows.
  misc_mt_manager_t testFitPool_ = nullptr;
  static constexpr size_t testFitParallelCutoff = 65536;
};

}  // namespace bartcore

#endif  // BARTCORE_CHAIN_HPP
