#ifndef BARTCORE_CHAIN_HPP
#define BARTCORE_CHAIN_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <type_traits>
#include <vector>

#include <external/io.h> // ext_throwError
#include <external/random.h>
#include <misc/linearAlgebra.h>
#include <misc/thread.h>

#include "combiner.hpp"
#include "data.hpp"
#include "grow.hpp"
#include "model.hpp"
#include "moves.hpp"
#include "scan.hpp"
#include "tree.hpp"

namespace bartcore {

/// A residual variance prior on the ORIGINAL response scale. A heteroscedastic
/// chain retains its creation value because these three numbers - and nothing
/// else - calibrate the scale leaf, once, at construction.
struct ResidualPrior {
  double sigmaEstimate = 0.0, sigmaDf = 0.0, sigmaRawScale = 0.0;
};

struct SamplerOptions {
  size_t numTrees = 200;
  size_t numChains = 1;
  // worker threads for running chains concurrently; only min(numThreads,
  // numChains) are used, and every chain needs its own non-R rng when > 1
  size_t numThreads = 1;
  // every numThin-th iteration is kept; numBurnIn and numSamples count at
  // the kept rate
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

  // Every predictor value the store ingests, in one borrowed view: the dense
  // block and/or the CSC (dgCMatrix-layout) triple, the per-column source map
  // that mixes them, and the typing channel (column types, declared level
  // counts, CSC reference codes). Unmapped, the x creation argument supplies
  // the dense block and the build retains nothing; mapped, the view's own
  // block is copied into the store and its CSC slices stay borrowed for the
  // store's lifetime. Consumed during construction, which clears it. See
  // PredictorSource (data.hpp) and docs/design/sparse-columns.md.
  PredictorSource predictors;

  // linear leaves: the ordinal predictor columns entering every leaf's
  // regression (borrowed; consumed during construction). Empty designates
  // the constant leaf; the factory validates count, range, and type.
  const std::size_t* leafCovariateColumns = nullptr;
  std::size_t numLeafCovariates = 0;

  // monotone (mBART) constraints: a borrowed per-predictor direction in
  // {-1, 0, +1} (length numPredictors), consumed at construction. A nonzero
  // entry selects the constrained constant-leaf instantiation at the factory;
  // null - or all zero - keeps the unchanged constant-leaf path
  // (docs/design/monotone.md). Not yet exposed through the R surface (C++/tests
  // only).
  const std::int8_t* monotoneDirections = nullptr;

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

  // AFT survival: per-observation status (1 = uncensored event, 0 = right-
  // censored) required when family is aft, ignored otherwise. Borrowed; the
  // response copies it during construction (see AFTResponse). The y creation
  // argument then holds the log survival/censoring times.
  const double* survivalStatus = nullptr;

  // split-variable selection: fixed weights (borrowed; normalized over
  // available variables at each node) or DART; both null/false = uniform
  const double* splitProbabilities = nullptr;
  bool useDart = false;
  DartPrior dart;

  // per-forest split-variable restriction: the column indices this forest may
  // split on (borrowed; consumed during construction). Null, or count 0,
  // leaves every column available - the default, nullptr-guarded so the
  // availability path is byte-for-byte unchanged.
  const std::size_t* forestColumns = nullptr;
  std::size_t numForestColumns = 0;

  // per-forest interaction constraint (docs/design/interaction-constraints.md):
  // interactionMaxOrder caps the DISTINCT split variables on
  // any root-to-leaf path (0 = uncapped); interactionForbiddenPairs lists
  // forbidden co-occurrence pairs as 2 * interactionNumForbiddenPairs column
  // indices (borrowed, consumed at construction). Both defaults leave every
  // path unconstrained - the availability path is byte-for-byte unchanged. Not
  // yet exposed through the R surface (C++/tests only; single-forest path).
  std::size_t interactionMaxOrder = 0;
  const std::size_t* interactionForbiddenPairs = nullptr;
  std::size_t interactionNumForbiddenPairs = 0;

  // per-forest block-additive constraint
  // (docs/design/interaction-constraints.md): confine each WHOLE tree to one
  // declared group of predictors so the
  // ensemble is exactly f = sum_G f_G. blockOfColumn (length numPredictors,
  // borrowed) gives each column's 0-based group (negative = in no block, only for a
  // restricted forest); blockTreeCounts (length numBlocks, borrowed) is the
  // deterministic contiguous per-group tree capacity, summing to numTrees. numBlocks
  // 0 (default) leaves every tree unrestricted - byte-for-byte unchanged. Consumed
  // at construction.
  std::size_t numBlocks = 0;
  const std::int32_t* blockOfColumn = nullptr;
  const std::size_t* blockTreeCounts = nullptr;

  // heteroscedastic variance forest (HBART, docs/design/heteroscedastic.md):
  // numVarianceTrees > 0 adds a SECOND forest modeling s^2(x) as a product of
  // scaled-inverse-chi-squared leaves, coupled to the mean forest through the
  // precision (weight) channel w_i^mean = user_w_i / s^2(x_i). 0 (default) is
  // homoscedastic - byte-for-byte unchanged, no variance forest built. Gaussian
  // family + plain constant leaf only; the factory refuses the combination
  // otherwise. varianceBase/variancePower are the variance trees' CGM prior.
  // Not yet exposed through the R surface (C++/tests only).
  std::size_t numVarianceTrees = 0;
  double varianceBase = 0.95, variancePower = 2.0;
  // optional split-variable restriction for the variance forest (borrowed
  // 0-based column indices, consumed at construction): the columns s^2(x) may
  // split on. Null or count 0 leaves every column available (all mean
  // predictors), the `variance = ~ x1 + x2` selector's default.
  const std::size_t* varianceForestColumns = nullptr;
  std::size_t numVarianceForestColumns = 0;

  // k is fixed at .k unless updateK, in which case .k is the starting value
  bool updateK = false;
  ChiKHyperprior kHyperprior;

  // gaussian responses only: sigma holds at the constructor's sigmaEstimate
  // and is never drawn (the R-level fixed() residual prior); binary families
  // are always sigma-fixed
  bool sigmaIsFixed = false;

  // continuous (gaussian family) responses only: a finite residualDf selects
  // Student-t errors via the scale-mixture augmentation (TResponse,
  // docs/design/robust-errors.md) - a positive value fixes nu, <= 0 estimates
  // it on the grid; NaN (default) and +Inf keep the Gaussian law. Non-gaussian
  // families ignore it (the bridge refuses the combination host-side).
  double residualDf = std::numeric_limits<double>::quiet_NaN();

  // ordinal (cumulative-probit) responses only: the number of ordered category
  // levels K (docs/design/ordinal.md), selecting OrdinalResponse with a K-1
  // cutpoint vector. 0 (default) is a non-ordinal response; the bridge sets it
  // from the ordered-factor level count and refuses K < 2.
  std::size_t numCategories = 0;

  // negative-binomial counts (nbinom family) only: the dispersion spec
  // (docs/design/negative-binomial.md sections 4, 5), the residualDf sign
  // convention - a positive value fixes r there (an integer), a non-positive
  // value estimates r on the capped grid. Only the nbinom construction reads it;
  // NaN (default) is ignored by every other family.
  double dispersion = std::numeric_limits<double>::quiet_NaN();

  // when set, every kept sample's trees are flattened into a circular buffer
  // of numSamplesToStore slots (at least 1) per chain, for prediction and
  // reporting after the run
  bool keepTrees = false;
  size_t numSamplesToStore = 0;

  // progress reporting during runs: one "iteration: k (of N)" line every
  // printEvery kept iterations
  bool verbose = false;
  std::uint32_t printEvery = 100;

  // opt-in fp32 running residual (docs/design/reduced-precision-storage.md):
  // stores Forest::treeY in fp32 with fp64-accumulated reductions,
  // halving the dominant suffstat gather's memory traffic at large n. The
  // factory mints the fp32 instantiation ONLY for the gaussian constant-leaf
  // path; the default (false) instantiation is byte-for-byte the fp64 engine.
  bool fp32Residual = false;
};

/// Receives chains' formatted progress lines during a run. Worker-thread runs
/// hand lines to a queue the main thread flushes (workers must never call into
/// R); inline runs print directly.
struct ProgressSink {
  virtual ~ProgressSink() = default;
  virtual void report(const char* line) = 0;
};

/// A between-run prior replacement: every field is applied unconditionally
/// except that DART samplers keep their Dirichlet split machinery (a
/// dbartsModel cannot express DART) and the fixed-sigma binary families
/// ignore the variance prior. Installing a model before any run matches
/// creating with it.
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
///
/// Under BCF (two forests) trainingFits carries the a * mu + b_z * tau blend,
/// but testFits is filled with NaN (no test treatment vector to blend), and
/// k, variableCounts, and splitProbabilities report the PROGNOSTIC (mu) forest
/// only; the treatment forest is reached through the per-forest channels
/// (forestTotalFits + bcfGlue, docs/design/bcf.md). logLikelihood is likewise
/// NaN-filled under BCF (the blended per-observation location is not visible
/// to the response model).
struct Results {
  double* sigma = nullptr;          // numSamples
  double* trainingFits = nullptr;   // numObservations x numSamples, or null
  double* testFits = nullptr;       // numTestObservations x numSamples, or null
  // numPredictors x numVariableCountForests x numSamples, or null; the forest
  // dimension is 1 for every additive model (the exact numPredictors x
  // numSamples layout), numCategories for multinomial's per-category channels
  std::uint32_t* variableCounts = nullptr;
  double* k = nullptr;              // numSamples, or null; only when k sampled
  // numPredictors x numSamples, or null; filled only under DART
  double* splitProbabilities = nullptr;
  // grouped samplers only, both on the original response scale:
  // tau is numSamples, groupEffects numGroups x numSamples
  double* tau = nullptr;
  double* groupEffects = nullptr;
  // per-draw training log-likelihood, numObservations x numSamples, or null;
  // gaussian and binary families, NaN under BCF
  double* logLikelihood = nullptr;
  // per-sample cutpoints, numCutpoints x numSamples, or null; filled only when
  // the response family carriesCutpoints() (ordinal's K-1 thresholds). Every
  // other family leaves it null and allocates nothing.
  double* cutpoints = nullptr;
  // heteroscedastic variance surface s^2(x), original response scale, or null:
  // varianceFits is numObservations x numSamples, varianceTestFits is
  // numTestObservations x numSamples. A SEPARATELY-typed-forest channel (not a
  // numReportedLocations widening of the response location); filled only when a
  // variance forest is present, left untouched otherwise.
  double* varianceFits = nullptr;
  double* varianceTestFits = nullptr;
  // per-draw per-forest reporting, both null unless the forest coupling defines
  // them (ForestCombiner::forestReportingIsDefined; BCF alone today), so every
  // other model allocates and computes nothing here: forestFits is
  // numObservations x numForests x numSamples, forest-major within a sample, on
  // the INTERNAL scale (the forests' own function values, as forestTotalFits
  // reports them); glue is 3 x numSamples, each column the (a, b0, b1) that
  // recombines them into that draw's location. A SEPARATELY-typed channel pair,
  // as the variance surface is, not a numReportedLocations widening.
  double* forestFits = nullptr;
  double* glue = nullptr;
  // per-observation channels the trainingFits/testFits arrays carry: 1 for
  // every additive model (the exact current layout), more for a multi-location
  // combiner. The run bridge sizes the fits buffers by it and Sampler strides
  // per chain by it; storeSample reads the count from the combiner directly.
  std::size_t numReportedLocations = 1;
  // per-sample forests the variableCounts array carries: 1 for every additive
  // model (the exact current layout), numCategories for multinomial. The run
  // bridge sizes variableCounts by it and Sampler strides per chain by it;
  // storeSample reads the count from the combiner directly.
  std::size_t numVariableCountForests = 1;
  // per-sample cutpoints the cutpoints array carries: 0 for every family but
  // ordinal (K-1). The run bridge sizes cutpoints by it and Sampler strides per
  // chain by it; storeSample reads the count from the response directly.
  std::size_t numCutpoints = 0;
};

/// A host's per-sweep conditioning hook, invoked before every sweep on the
/// running thread with the chain index, the 0-based sweep counter, and whether
/// the sweep is discarded burn-in; returns true to stop the run early, exactly
/// as shouldCancel does. Only set when chains run inline: worker-thread chains
/// must not call into R.
using SweepCallback =
  std::function<bool(std::size_t chainIndex, std::size_t sweepIndex,
                     bool isBurnIn)>;

/// The heteroscedastic variance ensemble (HBART; docs/design/heteroscedastic.md
/// sections 5-6): a second forest of ConstantVarianceLeaf trees whose product
/// s^2(x_i) = prod_j h_j(x_i) modulates the mean forest's precision through the
/// weight channel (w_i^mean = user_w_i / s^2(x_i)). It is distinctly typed from
/// the mean Forest<L, ResidT> - a scale leaf, and a MULTIPLICATIVE running residual
/// (divide by the OTHER trees' product s^2_{-j}, not the additive roll) - so a
/// Chain holds it as a nullable member, not another Forest<L, ResidT>. Gaussian only.
struct VarianceForest {
  std::size_t numTrees = 0;
  double birthOrDeathProbability = 0.5, swapProbability = 0.1,
         changeProbability = 0.4, birthProbability = 0.5;
  ConstantVarianceLeaf leaf;
  CGMTreePrior treePrior;
  std::vector<Tree> trees;
  std::vector<index_t> indexBuffer;
  // per-tree multiplicative factor h_j(x_i), tree-major (numTrees x n); the
  // combined variance s^2(x_i) = prod_j factorByTree[j * n + i]
  std::vector<double> factorByTree;
  std::vector<double> combinedVariance;      // s^2(x_i), length n
  std::vector<double> combinedVarianceTest;  // s^2(x_test_i), length nTest
  std::vector<double> meanResidual;      // scratch e_i = y_i - f(x_i), per sweep
  std::vector<double> divisor;           // scratch s^2_{-j}(x_i), per tree
  std::vector<double> treeResidual;      // scratch e_i / sqrt(s^2_{-j}), per tree
  // per-forest split-variable restriction (0/1 byte per predictor, 1 =
  // splittable); empty leaves every column available - the default path
  std::vector<std::uint8_t> columnMask;
  // saved-sample flattened variance trees (keepTrees), a circular buffer of
  // capacity slots x numTrees, slot-major, for posterior predict on new data
  std::vector<std::vector<FlatNode>> savedTrees;
  std::size_t savedTreeCapacity = 0, savedSlotBase = 0;
  MoveScratch scratch;

  std::size_t numObservations() const { return combinedVariance.size(); }

  /// Build numTrees root trees over n observations and seed every leaf factor
  /// so the product equals initialVariance (each root h = initialVariance^(1/m')).
  void initialize(std::size_t numTrees_, std::size_t n, double initialVariance) {
    numTrees = numTrees_;
    indexBuffer.assign(n * numTrees, 0);
    trees.resize(numTrees);
    for (std::size_t t = 0; t < numTrees; ++t)
      trees[t].initialize(indexBuffer.data() + t * n, n);
    double rootFactor =
      std::pow(initialVariance, 1.0 / static_cast<double>(numTrees));
    factorByTree.assign(numTrees * n, rootFactor);
    combinedVariance.assign(n, initialVariance);
    meanResidual.assign(n, 0.0);
    divisor.assign(n, 1.0);
    treeResidual.assign(n, 0.0);
  }

  /// The multiplicative roll for tree j: divisor holds s^2_{-j}(x_i) =
  /// s^2(x_i) / h_j(x_i), the OTHER trees' product EXCLUDING tree j, and
  /// treeResidual holds the scaled mean residual e_i / sqrt(s^2_{-j}). Excluding
  /// h_j is the divisor guard: perturbing tree j's own leaf
  /// must not move its own suffstat, only another tree's leaf may.
  void formTreeResidual(std::size_t j, const double* meanResidualIn) {
    std::size_t n = numObservations();
    const double* __restrict factor = factorByTree.data() + j * n;
    for (std::size_t i = 0; i < n; ++i) {
      divisor[i] = combinedVariance[i] / factor[i];
      treeResidual[i] = meanResidualIn[i] / std::sqrt(divisor[i]);
    }
  }

  /// Scatter tree j leaf b's drawn factor h to every member observation and
  /// fold it back into the combined variance (s^2 = s^2_{-j} * h). divisor must
  /// hold this tree's current s^2_{-j} (from formTreeResidual).
  void applyLeafFactor(std::size_t j, std::int32_t node, double h) {
    std::size_t n = numObservations();
    double* __restrict factor = factorByTree.data() + j * n;
    const Tree& tree = trees[j];
    const Node& b = tree.at(node);
    for (std::size_t m = b.begin; m < b.end; ++m) {
      std::size_t i = tree.indices[m];
      factor[i] = h;
      combinedVariance[i] = divisor[i] * h;
    }
  }
};

/// Bank count of the fused roll + node-average suffstat pass. K fixes the
/// summation order and so is part of the draw law: a knob here would mean the
/// package has no single law, with every baseline and every bug report
/// conditional on it. Measured: K = 1 loses (the single-accumulator FP
/// dependency chain), K = 4 is green on both architectures
/// (docs/design/memory-wall-frontier.md secs 11-12).
constexpr std::size_t fusedSuffstatBanks = 4;

/// Fixed-order combine of one node's banks: ((b0 + b1) + b2) + b3, strictly
/// left to right, banks laid out bank-major with `stride` slots each. Any
/// regrouping (pairwise, SIMD reduce) is a different sum and a different draw.
inline double combineFusedSuffstatBanks(const double* acc, std::size_t slot,
                                        std::size_t stride) {
  double sum = acc[slot];
  for (std::size_t b = 1; b < fusedSuffstatBanks; ++b)
    sum += acc[b * stride + slot];
  return sum;
}

/// What Chain::checkFusedSuffstatAgainstStockForTesting reports.
struct FusedSuffstatCheck {
  /// Trees whose eligibility predicate accepted; the rest declined to stock.
  std::size_t numTreesFused = 0;
  /// The fused roll wrote treeY bit-for-bit as rollTreeResidual does.
  bool residualsAgreeBitwise = true;
  /// sumWeights agreed bitwise (both sides report the member count).
  bool countsAgreeBitwise = true;
  /// Worst gap in sumWeightedResponse, relative to max(|stock|, 1): the two
  /// legitimately differ there, since the fused banks and the stock
  /// unroll-by-5 associate differently, and a well-fit leaf's sum passes
  /// through zero, where a purely relative measure is meaningless.
  double worstRelativeError = 0.0;
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
template <IntegrableLeafModel L, typename ResidT = double>
class Chain {
public:
  static_assert(std::is_same_v<ResidT, double> || std::is_same_v<ResidT, float>,
                "the residual storage type is fp64 (default) or opt-in fp32");
  // ResidT is the running-residual (treeY) element type: double by default,
  // byte-identical to the fp64 engine; float only for the opt-in gaussian
  // constant-leaf path the factory mints (reduced-precision-storage.md sec 3b).

  /// Scalar and function-valued leaves read per-node weighted means
  /// (function leaves because their over-cap nodes delegate to the constant
  /// leaf); vector leaves accumulate their own statistics.
  static constexpr bool leafTracksNodeAverages = !L::hasVectorParams;

  /// The constant leaf stores compact node-indexed mu tables plus a per-tree
  /// obs-to-leaf map in place of the dense per-tree fit slab; vector and
  /// function leaves keep the slab.
  static constexpr bool leafIsConstant =
    !L::hasVectorParams && !L::hasFunctionParams;

  Chain(const ColumnStore& data, const double* y, const double* weights,
        const double* offset, ResponseFamily family, double sigmaEstimate,
        double sigmaDf, double sigmaRawScale, const SamplerOptions& options,
        ext_rng* rng)
    : options_(options), data_(data), weights_(weights), rng_(rng) {
    size_t numObservations = data.numObservations;
    options_.maxNumCutsPerVariable = nullptr;  // consumed by the store build
    options_.predictors = {};

    forests_.emplace_back();
    Forest<L, ResidT>& forest = forests_.back();
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
                             options.numLeafCovariates, options.numChains);
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
      // a finite residualDf selects the Student-t scale-mixture error law
      // (TResponse: > 0 fixes nu, <= 0 estimates it on the grid); NaN and Inf
      // keep the Gaussian law (docs/design/robust-errors.md)
      if (std::isfinite(options.residualDf))
        response_ = std::make_unique<TResponse>(
          y, offset, weights, numObservations, sigmaEstimate, sigmaDf,
          sigmaRawScale, options.residualDf);
      else
        response_ = std::make_unique<GaussianResponse>(
          y, offset, weights, numObservations, sigmaEstimate, sigmaDf,
          sigmaRawScale);
      break;
    case ResponseFamily::aft:
      // y holds the log survival/censoring times; status marks the censored
      response_ = std::make_unique<AFTResponse>(
        y, options.survivalStatus, offset, numObservations, sigmaEstimate,
        sigmaDf, sigmaRawScale);
      break;
    case ResponseFamily::ordinal:
      // y holds one-based category indices in {1..K}; a ConstantGaussianLeaf
      // single-forest model like probit, sigma fixed at 1 (docs/design/ordinal.md)
      response_ = std::make_unique<OrdinalResponse>(y, offset, numObservations,
                                                    options.numCategories);
      break;
    case ResponseFamily::nbinom:
      // y holds non-negative counts; the forest fits the log-odds latent under
      // the Polya-Gamma augmentation, sigma fixed at 1, with dispersion r fixed
      // (options.dispersion > 0) or grid-estimated (docs/design/negative-binomial.md)
      response_ = std::make_unique<NBResponse>(y, offset, numObservations,
                                               options.dispersion);
      break;
    }
    options_.survivalStatus = nullptr;  // consumed above
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
    // the constrained constant leaf reads box geometry from the store and its
    // per-column directions from the borrowed spec; the c-inflation matches a
    // constrained leaf's post-truncation variance to the unconstrained prior
    if constexpr (TreeDrawLeafModel<L>) {
      forest.leaf.data = &data_;
      forest.leaf.directions.assign(
        options.monotoneDirections,
        options.monotoneDirections + data.numPredictors);
      forest.leaf.cInflation =
        std::sqrt(std::numbers::pi / (std::numbers::pi - 1.0));
    }
    forest.treePrior.base = options.base;
    forest.treePrior.power = options.power;
    family_ = family;
    // aft draws sigma conjugately like gaussian; only the binary families fix it
    sigmaIsFixed_ = (family != ResponseFamily::gaussian &&
                     family != ResponseFamily::aft) ||
                    options.sigmaIsFixed;

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

    // A restricted forest clears the availability of every unlisted column; an
    // empty list leaves the mask empty and the trees unrestricted (the default
    // availability path, byte-for-byte).
    if (options.forestColumns != nullptr && options.numForestColumns > 0) {
      forest.columnMask.assign(data.numPredictors, 0);
      for (size_t c = 0; c < options.numForestColumns; ++c)
        forest.columnMask[options.forestColumns[c]] = 1;
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].setColumnMask(forest.columnMask.data());
    }
    options_.forestColumns = nullptr;  // consumed above

    // Per-forest interaction constraint, installed like the column mask: an
    // unset (or inactive) constraint leaves every tree's pointer null and the
    // availability path byte-for-byte unchanged.
    if (options.interactionMaxOrder > 0 || options.interactionNumForbiddenPairs > 0) {
      forest.interaction = std::make_unique<InteractionConstraint>();
      forest.interaction->build(data.numPredictors, options.interactionMaxOrder,
                                options.interactionForbiddenPairs,
                                options.interactionNumForbiddenPairs);
      if (forest.interaction->active())
        for (size_t t = 0; t < forest.numTrees; ++t)
          forest.trees[t].setInteractionConstraint(forest.interaction.get());
      else
        forest.interaction.reset();
    }
    options_.interactionForbiddenPairs = nullptr;  // consumed above

    // Block-additive constraint: confine each tree to one group. The
    // single-forest path carries no base columnMask, so the block row is the
    // group membership directly; numBlocks 0 leaves every tree unrestricted.
    installBlockMasks(forest, data.numPredictors, options.numBlocks,
                      options.blockOfColumn, options.blockTreeCounts);
    options_.blockOfColumn = nullptr;    // consumed above
    options_.blockTreeCounts = nullptr;  // consumed above

    initForestFitStorage(forest, numObservations);
    forest.totalFits.assign(numObservations, 0.0);
    forest.treeY.resize(numObservations);
    forest.paramByNode.clear();
    if constexpr (L::hasVectorParams)
      forest.paramsByTree.assign(forest.numTrees,
                                 std::vector<double>(forest.leaf.numParams(), 0.0));

    // heteroscedastic: a second, distinctly typed variance forest coupled to
    // this constant-leaf gaussian mean forest through the weight channel. Built
    // only for the plain constant leaf and gaussian family (the factory refuses
    // every other combination); homoscedastic leaves varianceForest_ null and
    // this whole path compiled/branched out, so the mean sweep is byte-identical.
    if constexpr (std::is_same_v<L, ConstantGaussianLeaf>)
      if (family == ResponseFamily::gaussian && options.numVarianceTrees > 0)
        buildVarianceForest(options, sigmaEstimate, sigmaDf, sigmaRawScale);

    resizeTestStorage();
  }

  /// Two-forest BCF chain (docs/design/bcf.md): a prognostic forest mu
  /// (forest 0) and a treatment forest tau (forest 1) combined on a gaussian
  /// response as y = a mu + b_z tau + eps. Constant leaves only; both forests
  /// read the full store unless tau carries a moderator subset in
  /// BCFForestSpec.columns (the default, no list, restricts nothing).
  Chain(const ColumnStore& data, const double* y, const double* weights,
        const double* offset, double sigmaEstimate, double sigmaDf,
        double sigmaRawScale, const SamplerOptions& options,
        const BCFSpec& spec, ext_rng* rng)
    : options_(options), data_(data), weights_(weights), rng_(rng) {
    static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                  "BCF is a constant-leaf model");
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
    options_.forestColumns = nullptr;  // BCF restriction arrives via BCFForestSpec
    response_ = std::make_unique<GaussianResponse>(
      y, offset, weights, data.numObservations, sigmaEstimate, sigmaDf,
      sigmaRawScale);
    family_ = ResponseFamily::gaussian;
    sigmaIsFixed_ = options.sigmaIsFixed;
    sigma_ = response_->initialSigma();

    // Calibration map (docs/design/bcf.md): s is the sample sd of the
    // range-scaled response (y mapped to [-0.5, 0.5]). The prognostic total
    // mu ~ N(0, s^2) so the half-Cauchy a (median aPriorScale) puts it at
    // aPriorScale sd(y); the treatment total tau ~ N(0, (sdModerate s /
    // 0.674)^2) so with b1 - b0 ~ N(0, 2 bPriorVariance) and half-normal
    // median 0.674 the effect (b1 - b0) tau sits at sdModerate sd(y). The
    // map fixes k at 1 and overrides the host node prior for both forests.
    constexpr double kHalfNormalMedian = 0.674;
    double s = scaledResponseSd();
    buildBCFForest(spec.mu, s);
    buildBCFForest(spec.tau, spec.sdModerate * s / kHalfNormalMedian);

    combiner_ = std::make_unique<BCFForestCombiner<L, ResidT>>(data, spec);
    resizeTestStorage();
  }

  /// K-forest multinomial (softmax) chain (docs/design/multinomial.md): K
  /// symmetric constant-leaf category forests coupled through a softmax
  /// likelihood with an interleaved one-vs-rest Polya-Gamma augmentation and a
  /// likelihood-invariant level-centering move, all owned by
  /// MultinomialForestCombiner. The grouped-count response (an n x K count matrix
  /// and per-observation trials n_i, single-trial labels entering as one-hot with
  /// n_i = 1) rides the spec; there is no sigma (fixed, like the binary families).
  Chain(const ColumnStore& data, const SamplerOptions& options,
        const MultinomialSpec& spec, ext_rng* rng)
    : options_(options), data_(data), weights_(nullptr), rng_(rng) {
    static_assert(!L::hasVectorParams && !L::hasFunctionParams,
                  "multinomial is a constant-leaf model");
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
    options_.forestColumns = nullptr;
    response_ = std::make_unique<MultinomialResponse>(data.numObservations);
    // logistic marks the binary-family sigma semantics (fixed at 1); the
    // softmax has no sigma of its own, and family() is not read on this path.
    family_ = ResponseFamily::logistic;
    sigmaIsFixed_ = true;
    sigma_ = response_->initialSigma();

    for (std::size_t k = 0; k < spec.numCategories; ++k)
      buildMultinomialForest(spec.forest, spec.nodeScale, spec.k);

    combiner_ =
      std::make_unique<MultinomialForestCombiner<L, ResidT>>(data, spec);
    resizeTestStorage();
  }

  ~Chain() {
    if (testFitPool_ != nullptr) misc_mt_destroy(testFitPool_);
  }

  std::size_t numForests() const { return forests_.size(); }
  std::size_t numTreesInForest(std::size_t f) const {
    return forests_[f].numTrees;
  }
  /// The current combined variance s^2(x_i), working scale, over the training
  /// (varianceFits) or test (varianceTestFits) rows, or null when homoscedastic.
  /// Original-scale reporting multiplies by sigmaScale^2 (storeSample, predict).
  const double* varianceFits() const {
    return varianceForest_ ? varianceForest_->combinedVariance.data() : nullptr;
  }
  const double* varianceTestFits() const {
    return varianceForest_ ? varianceForest_->combinedVarianceTest.data()
                           : nullptr;
  }
  bool hasVarianceForest() const { return varianceForest_ != nullptr; }
  std::size_t numVarianceTrees() const {
    return varianceForest_ ? varianceForest_->numTrees : 0;
  }
  /// The residual prior the scale leaf was calibrated from at creation;
  /// meaningful only under a variance forest, where nothing recalibrates it.
  const ResidualPrior& varianceLeafPrior() const { return varianceLeafPrior_; }

  /// s^2(x) on the ORIGINAL scale for new rows of a Columns predictor source
  /// from one saved sample's variance trees; the per-tree factors MULTIPLY
  /// into the product.
  template <typename Columns>
  void predictVarianceFromSavedSample(std::size_t slot, const Columns& columns,
                                      std::size_t numTest, double* out) {
    VarianceForest& vf = *varianceForest_;
    for (std::size_t i = 0; i < numTest; ++i) out[i] = 1.0;
    std::vector<std::size_t> indices(numTest);
    std::vector<std::size_t> blockOffsets;
    std::vector<double> treeFit(numTest);
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      misc_setVectorToConstant(treeFit.data(), numTest, 0.0);
      addFlatPredictions(vf.savedTrees[slot * vf.numTrees + j], nullptr,
                         nullptr, columns, numTest, indices, blockOffsets,
                         treeFit.data());
      for (std::size_t i = 0; i < numTest; ++i) out[i] *= treeFit[i];
    }
    double s = response_->sigmaScale();
    for (std::size_t i = 0; i < numTest; ++i) out[i] *= s * s;
  }
  /// The response's original-scale factor (range for gaussian); a variance is
  /// reported on the original scale as the working s^2 times its square.
  double sigmaScale() const { return response_->sigmaScale(); }

  /// The variance surface s^2(x) on the ORIGINAL response scale for new rows
  /// of a Columns predictor source, from the current variance trees; null-safe
  /// no-op off a variance forest. Each tree is flattened and replayed like the
  /// mean forest's predict, but the per-tree factors MULTIPLY into the product.
  template <typename Columns>
  void predictVariance(const Columns& columns, std::size_t numTest,
                       double* out) {
    if (!varianceForest_) return;
    VarianceForest& vf = *varianceForest_;
    for (std::size_t i = 0; i < numTest; ++i) out[i] = 1.0;
    std::vector<std::size_t> indices(numTest);
    std::vector<std::size_t> blockOffsets;
    std::vector<double> leafValues, treeFit(numTest);
    std::vector<FlatNode> flat;
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      recoverVarianceLeafValues(vf, j, leafValues);
      vf.trees[j].flatten(data_, leafValues.data(), flat);
      misc_setVectorToConstant(treeFit.data(), numTest, 0.0);
      addFlatPredictions(flat, nullptr, nullptr, columns, numTest, indices,
                         blockOffsets, treeFit.data());
      for (std::size_t i = 0; i < numTest; ++i) out[i] *= treeFit[i];
    }
    double s = response_->sigmaScale();
    for (std::size_t i = 0; i < numTest; ++i) out[i] *= s * s;
  }
  /// Per-observation channels the recorded fits carry: the combiner's location
  /// count (1 for every additive combiner, BCF included), 1 off any combiner.
  std::size_t numReportedLocations() const {
    return combiner_ ? combiner_->numReportedLocations() : 1;
  }
  /// Per-sample forests the recorded split-usage channel carries: the
  /// combiner's forest count (1 for every additive combiner, BCF included),
  /// 1 off any combiner. The run bridge reads it to size the varcount array.
  std::size_t numVariableCountForests() const {
    return combiner_ ? combiner_->numVariableCountForests() : 1;
  }
  /// Per-sample cutpoints the recorded cutpoint channel carries: the response's
  /// K-1 for a cutpoint-carrying family (ordinal), 0 for every other.
  std::size_t numCutpoints() const { return response_->numCutpoints(); }
  /// Whether the recorded test-fit channel carries a defined value: true off any
  /// combiner (single forest) and for a combiner that blends a test surface
  /// (multinomial softmax), false for BCF (no test treatment vector). The bridge
  /// gates its test-surface refusal on this so BCF stays refused while
  /// multinomial and single-forest samplers are allowed.
  bool testFitsAreDefined() const {
    return combiner_ ? combiner_->testFitsAreDefined() : true;
  }
  /// Whether the recorded per-forest fits and glue channels carry defined
  /// values: the coupling's own answer (true for BCF), false off any combiner.
  /// The run bridge reads it to decide whether to allocate them at all.
  bool forestReportingIsDefined() const {
    return combiner_ != nullptr && combiner_->forestReportingIsDefined();
  }
  bool usesDart() const { return forests_[0].useDart; }
  /// Re-forms b_{z_i} and both residuals on the next sweep; z is borrowed.
  void setTreatment(const double* z) {
    if (combiner_) combiner_->setTreatment(z);
  }
  /// Whether this chain's forest coupling permits the response-side conduit -
  /// setResponse and setOffset at updateScale = false, and setWeights, which
  /// has no scale to pin at all; false off any combiner, since a single-forest
  /// chain never consults it. The gaussian conjunct is
  /// structural, not advisory: a latent family reads forests_[0].totalFits as
  /// though it were the combined location, which is false on a coupling.
  bool supportsResponseMutation() const {
    return combiner_ && combiner_->supportsResponseMutation() &&
           family_ == ResponseFamily::gaussian;
  }
  /// Whether this chain admits a caller-supplied per-forest weight. The
  /// setter's refusal and SamplerShape::supportsForestWeights both read this
  /// one predicate, so the advertised capability and the refusal cannot
  /// disagree.
  bool supportsForestWeights() const {
    return combiner_ != nullptr && combiner_->supportsForestWeights();
  }

  /// Installs a BORROWED per-observation weight s on forest f, clearing it at a
  /// null pointer; returns false, installing nothing, when the coupling admits
  /// no such weight or f names no forest.
  ///
  /// s is a multiplicative PRECISION factor on forest f's own leaf
  /// conditionals, composing with the observation weight so that forest f's
  /// draws see w_i m_f^2 s_i. It does NOT remove the row from occupancy, from
  /// the empty-leaf veto, from the combination (the row still receives
  /// m_f f_f(x_i)), or from the residual sigma degrees of freedom, which count
  /// positive OBSERVATION weights; s_i = 0 says only that row i carries no
  /// information about forest f, and its leaves stay well-defined prior draws.
  ///
  /// Two edges a consumer is misled without: at s_i = 0 with a nonzero
  /// multiplier only the WEIGHT is zeroed - the response stays the
  /// reparameterized residual r_i / m_f - so the reported-fit exactness an
  /// exactly zero multiplier buys does NOT follow this channel; and the weight
  /// lives on the chain rather than in the serialized state, so a pipeline that
  /// REBUILDS a sampler and restores its state silently drops the weight and
  /// draws a different model while the states still agree.
  bool setForestWeights(std::size_t f, const double* s) {
    if (!supportsForestWeights() || f >= forests_.size()) return false;
    if (forestWeights_.empty()) forestWeights_.assign(forests_.size(), nullptr);
    forestWeights_[f] = s;
    return true;
  }

  /// BCF glue on the combining response; false for a non-BCF chain.
  bool bcfGlue(double& a, double& b0, double& b1) const {
    return combiner_ ? combiner_->bcfGlue(a, b0, b1) : false;
  }
  /// The forest's constant-leaf function values on the internal scale (mu for
  /// forest 0, tau for forest 1); numObservations doubles.
  void forestTotalFits(std::size_t f, double* out) const {
    std::memcpy(out, forests_[f].totalFits.data(),
                data_.numObservations * sizeof(double));
  }

  /// Forest f's per-predictor split usage, accumulated across its trees into
  /// out (numPredictors entries, zeroed here); the per-forest analog of the
  /// reported-forest variable-count channel storeSample records, addressing an
  /// arbitrary forest so a multi-forest model can report each forest's splits.
  void forestVariableCounts(std::size_t f, std::uint32_t* out) const {
    const Forest<L, ResidT>& forest = forests_[f];
    std::memset(out, 0, data_.numPredictors * sizeof(std::uint32_t));
    for (std::size_t t = 0; t < forest.numTrees; ++t)
      forest.trees[t].countVariableUses(out);
  }

  /// Forest f's per-tree fit slabs, tree-major (numObservations x numTrees); a
  /// consistency read of the cached fits for tests.
  void forestTreeFits(std::size_t f, double* out) const {
    const Forest<L, ResidT>& forest = forests_[f];
    size_t n = data_.numObservations;
    if constexpr (leafIsConstant) {
      // materialize the compact fits by gather (identical bytes to the slab)
      for (size_t t = 0; t < forest.numTrees; ++t) {
        const double* mu = forest.muByTree[t].data();
        const std::uint32_t* leaf = forest.leafOf.data() + t * n;
        double* o = out + t * n;
        for (size_t i = 0; i < n; ++i) o[i] = mu[leaf[i]];
      }
    } else {
      std::memcpy(out, forest.treeFits.data(),
                  n * forest.numTrees * sizeof(double));
    }
  }

  /// Fires the combiner's post-combine move (BCF: the interweaving glue-ridge
  /// rescale, BCFForestCombiner<L>::afterCombine) outside a sweep, for the
  /// component tests; returns the applied scale (1.0 off BCF or when the move
  /// is skipped).
  double interweaveGlueRidge(bool record = false, std::size_t sampleNum = 0) {
    return combiner_
      ? combiner_->afterCombine(forests_, record, sampleNum, rng_)
      : 1.0;
  }

  /// Between-run reconfiguration; the test-fit pool is rebuilt lazily to
  /// the new share of the budget on the next routing.
  void setNumThreads(size_t numThreads) { options_.numThreads = numThreads; }

  /// Called after the shared store's test data changes.
  void resizeTestStorage() {
    for (Forest<L, ResidT>& forest : forests_) {
      forest.totalTestFits.assign(data_.numTestObservations, 0.0);
      forest.currTestFits.resize(data_.numTestObservations);
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        forest.leaf.rebuildTestCovariates(data_);
    }
    if (varianceForest_)
      varianceForest_->combinedVarianceTest.assign(data_.numTestObservations,
                                                   1.0);
  }

  /// One run of (numBurnIn + numSamples) * numThin sweeps, recording the
  /// numSamples post-burn-in kept draws; results slots may be null to skip
  /// recording.
  /// Touches only chain state and the read-only store: safe to run chains
  /// concurrently as long as each has its own rng that never calls into R.
  /// progress, when non-null under verbose, receives one formatted line per
  /// printEvery kept iterations.
  /// Runs the chain, returning true if it stopped early because shouldCancel
  /// (polled once per sweep, called only on the thread that owns this chain)
  /// or onSweep asked it to. Both touch no sampled state, so a run with neither
  /// set is bitwise identical to one without them.
  bool run(size_t numBurnIn, size_t numSamples, Results& results,
           ProgressSink* progress = nullptr, size_t chainIndex = 0,
           const std::function<bool()>* shouldCancel = nullptr,
           const SweepCallback* onSweep = nullptr) {
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
      // the host's conditioning hook fires unthrottled before the sweep it
      // conditions; returning true stops as shouldCancel does
      if (onSweep != nullptr &&
          (*onSweep)(chainIndex, iteration, iteration / numThin < numBurnIn))
        return true;

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

      // heteroscedastic: the mean forest backfits against precisions divided by
      // the current variance surface, w_i^mean = user_w_i / s^2(x_i), with the
      // global sigma fixed at 1 (the variance forest IS the residual variance).
      // Homoscedastic skips this and the mean sweep sees the shared weights.
      if (varianceForest_) {
        formMeanWeights();
        weights = meanWeights_.data();
      }

      for (size_t f = 0; f < forests_.size(); ++f) {
        Forest<L, ResidT>& forest = forests_[f];
        // single-forest samplers backfit against the shared working response
        // and weights unchanged; a BCF forest sees its residual net of the
        // other forest's scaled contribution, divided by its own multiplier
        const double* forestY = y;
        const double* forestWeights = weights;
        if (combiner_) {
          // the interleaved coupling draws forest f's latents against the
          // current margins here, immediately before formForestResponse reads
          // them (a no-op for BCF); base no-op keeps every additive path bitwise
          combiner_->drawForestGlue(f, rng_, forests_);
          ForestResponse fr = combiner_->formForestResponse(f, forests_, y,
                                                            weights);
          forestY = fr.response;
          forestWeights = composeForestWeights(f, fr.weights);
        }
        MoveContext ctx{data_,
                        forest.treePrior,
                        forest.birthOrDeathProbability,
                        forest.swapProbability,
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
          double* treeFits = nullptr;
          if constexpr (!leafIsConstant)
            treeFits = forest.treeFits.data() + t * n;

          // the roll and the node-average suffstat go in one obs-order pass
          // where the fusion is eligible; everywhere else the stock pair runs
          if (!rollAndSetNodeAveragesFused(forest, t, forestY, forestWeights)) {
            rollTreeResidual(forest, t, forestY);

            // constant-leaf node means, recomputed against this sweep's residual
            if constexpr (leafTracksNodeAverages)
              forest.trees[t].setNodeAverages(forest.treeY.data(),
                                              forestWeights);
          }

          bool stepTaken;
          StepType stepType;
          int32_t changedNode = invalidNode;
          // hand a constrained-conjugate leaf its persistent mu block so a
          // branch score can read frozen neighbor values; compiled out for the
          // conjugate leaves, which read no leaf parameters
          if constexpr (leafIsConstant && ParamScoringLeafModel<L>)
            ctx.leafParams = forest.muByTree[t].data();
          metropolisJumpForTree(ctx, forest.leaf, rng_, forest.trees[t],
                                forest.treeY.data(), sigma_, &stepTaken,
                                &stepType, &changedNode);
          // accepted changes and deaths strand pooled mask words; no rule
          // copies are live here, so this is a safe point to reclaim them
          if (data_.hasPooledCategorical)
            forest.trees[t].compactMaskPoolIfNeeded(data_);

          // leafOf catches up with the settled move: rejections restore the
          // partition exactly and write nothing, an accepted move patches only
          // its repartitioned subtree, and a tree marked stale (wholesale
          // structure reset with fits left cached) rebuilds in full
          if constexpr (leafIsConstant) {
            bool wasStale = forest.leafOfStale[t] != 0;
            if (wasStale) {
              rebuildLeafOf(forest, t);
            } else if (stepTaken) {
              updateLeafOfBelow(forest, t, changedNode);
            }
            // the constrained leaf reads muByTree through the move phase (frozen
            // neighbors) and its draw updates the surviving block in place, so it
            // must stay sized and feasible across an accepted birth/death rather
            // than be zeroed and refilled. A birth seeds its two children with
            // the parent's (feasible) value; a death seeds the merged leaf inside
            // its neighbor bounds; the draw then redraws every leaf.
            if constexpr (TreeDrawLeafModel<L>)
              maintainMonotoneLeafStore(forest, t, wasStale, stepTaken, stepType,
                                        changedNode);
          }

          // the draw writes this tree's new leaf values: the constant leaf's mu
          // table, or the dense slab
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

        // rebuild the running total for the latent/sigma updates and recording;
        // the last tree's new fits retire here instead of in a pass of their own
        finalizeTotalFits(forest, forestY);
      }

      // a single forest reports its own fits; BCF the a mu + b_z tau blend
      const double* combined = combinedFits();
      response_->refreshLatents(rng_, combined, sigma_);
      y = response_->workingResponse();
      weights = response_->workingWeights();
      // a latent family's refresh changes the weights U'WU is cached against
      if constexpr (L::hasVectorParams)
        if (response_->workingWeightsVaryPerSweep())
          for (Forest<L, ResidT>& forest : forests_)
            forest.leaf.invalidateStatistics();

      if (!sigmaIsFixed_)
        sigma_ = response_->drawSigma(rng_, combined, sigma_);

      if (combiner_) {
        combiner_->drawGlue(rng_, sigma_, y, weights, forests_);
        combiner_->afterCombine(forests_, record, sampleNum, rng_);
      }

      // heteroscedastic: after the mean forest settles (sigma stays fixed at 1),
      // backfit the variance forest against the mean residual and refresh the
      // combined s^2(x) that next sweep's mean weights divide by.
      if (varianceForest_) {
        sweepVarianceForest(y, combined);
        if (record && data_.numTestObservations > 0) refreshVarianceTestFits();
        if (record && varianceForest_->savedTreeCapacity > 0)
          storeVarianceSavedTrees(sampleNum);
      }

      for (Forest<L, ResidT>& forest : forests_) {
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
    if constexpr (L::hasVectorParams)
      forests_[0].leaf.invalidateStatistics();
  }
  void setResponse(const double* y, bool updateScale) {
    response_->setResponse(y, rng_, forests_[0].totalFits.data(), updateScale,
                           &sigma_);
    // a latent family's setResponse refreshes the Polya-Gamma weights U'WU
    // depends on; a gaussian one moves only the residual
    if constexpr (L::hasVectorParams)
      if (response_->workingWeightsVaryPerSweep())
        forests_[0].leaf.invalidateStatistics();
  }
  /// Unguarded for the structurally pinned binary families: their restore
  /// paths (setState and installForest) reinstall the donor's own pinned
  /// value, so the write is benign, and the user-facing change is refused at
  /// the bridge (refusePinnedSigmaChange). Under a variance forest sigma is
  /// not a parameter at all - buildVarianceForest pins it at 1 on the working
  /// scale and the variance surface carries the residual scale from there - so
  /// every write is dropped, which is what closes the internal callers that
  /// have no error channel.
  void setSigma(double sigmaOriginalScale) {
    if (varianceForest_) return;
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

  /// Install a replacement prior; see ModelParameters for the semantics.
  void setModel(const ModelParameters& model) {
    Forest<L, ResidT>& forest = forests_[0];
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

    // both arms move a variance forest's pin - the fixed one by installing an
    // estimate over the working-scale 1, the other by re-prioring and
    // unpinning - so the whole clause is skipped there; the incoming residual
    // prior addresses the scale leaf instead (below)
    if (!varianceForest_ &&
        (family_ == ResponseFamily::gaussian ||
         family_ == ResponseFamily::aft)) {
      sigmaIsFixed_ = model.sigmaIsFixed;
      if (model.sigmaIsFixed)
        setSigma(model.sigmaEstimate);
      else
        response_->setSigmaPrior(model.sigmaEstimate, model.sigmaDf,
                                 model.sigmaRawScale);
    }

    // Under a variance forest the residual prior calibrates the scale leaf
    // rather than sigma - the same three numbers buildVarianceForest consumed,
    // through the same conversion to the working scale (initialSigma =
    // sigmaEstimate / sigmaScale). The retained triple moves with the leaf, so
    // it always names the prior the live calibration came from. The next sweep
    // redraws every factor under the new prior; the current surface is left
    // alone, as a homoscedastic setSigmaPrior leaves the current sigma.
    // KNOWN GAP: sigmaScale() is read here and at creation, so an intervening
    // updateScale response/offset swap recalibrates onto the NEW scale while an
    // untouched leaf keeps the old one.
    if (varianceForest_) {
      double workingSigma = model.sigmaEstimate / response_->sigmaScale();
      varianceForest_->leaf = ConstantVarianceLeaf::calibrated(
        model.sigmaDf, workingSigma * workingSigma * model.sigmaRawScale,
        varianceForest_->numTrees);
      varianceLeafPrior_ = {model.sigmaEstimate, model.sigmaDf,
                            model.sigmaRawScale};
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
  /// current cut grid, empty leaves collapsed, and return the forest to the
  /// zero-fit state a freshly constructed chain carries: every tree's leaf
  /// parameters are zero, hence every tree's fit is identically zero,
  /// totalFits and totalTestFits are zero, and the constant leaf's obs-to-leaf
  /// map is all-root. The cached-fits evaluator the next sweep's residual roll
  /// reads - mu[leafOf] for the constant leaf, the dense treeFits row
  /// otherwise - therefore agrees with totalFits, and every map entry is in
  /// bounds for the reset arena. Leaving the pre-reset fits behind instead is
  /// what broke once: zeroing mu alone left totalFits summing fits the cache
  /// could no longer reproduce, and the displacement that injects is carried
  /// forward by every later sweep rather than healing.
  ///
  /// The reset is FOREST-ONLY. sigma, k, the DART probabilities, the BCF /
  /// multinomial glue, the variance forest, the saved-tree buffers and the
  /// response's latent block (z, omega, cutpoints, dispersion) are untouched,
  /// exactly as sampleNodeParametersFromPrior leaves them; a caller wanting a
  /// true restart on a latent family follows with setResponse.
  ///
  /// Two deviations from a freshly constructed chain are deliberate, and both
  /// are load-bearing: leafOfStale[t] is left at 1 where a fresh chain has 0,
  /// and muByTree[t] is sized to the drawn tree where a fresh chain holds a
  /// single root.
  ///
  /// Do NOT rebuild the map with rebuildLeafOf (or installLeafOfAndAddToTotal)
  /// here. Either clears leafOfStale, which makes the fused roll+suffstat pass
  /// eligible on the first post-reset sweep; that pass is bitwise-identical on
  /// the residual but deliberately NOT on the suffstat association, so it
  /// would move the draws of every call site - bart2's default init is this
  /// entry with no parameter draw after it, so its sweep 1 depends on the
  /// decline. The all-root memset keeps the decline and is in bounds for every
  /// tree while gathering mu[0] == 0, so the roll's cached term is exactly
  /// zero against the zeroed totalFits.
  void sampleTreesFromPrior() {
    size_t n = data_.numObservations;
    const double* y = response_->workingResponse();
    const double* weights = response_->workingWeights();
    std::vector<double> paramByNode;
    for (Forest<L, ResidT>& forest : forests_) {
      for (size_t t = 0; t < forest.numTrees; ++t) {
        forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
        growSubtreeFromPrior(forest, forest.trees[t], 0, y, weights);
        paramByNode.assign(forest.trees[t].nodes.size(), 0.0);
        forest.trees[t].collapseEmptyNodes(data_, weights, paramByNode);
        // fresh structures carry zero parameter blocks until the next draw
        if constexpr (L::hasVectorParams)
          forest.paramsByTree[t].assign(
            forest.trees[t].nodes.size() * forest.leaf.numParams(), 0.0);
        if constexpr (leafIsConstant) {
          // mu must be resized to the grown tree here: a ParamScoringLeafModel
          // (monotone) reads muByTree for frozen neighbors DURING the first
          // move, before the stale rebuild in maintainMonotoneLeafStore runs,
          // so a block left at its size-1 root would be read out of bounds.
          // The all-equal zero seed is monotone-feasible (every neighbor bound
          // holds with equality) and is exactly what the first draw reassigns,
          // so the conjugate constant leaf - which never reads mu in the move
          // - is byte-unchanged.
          forest.muByTree[t].assign(forest.trees[t].nodes.size(), 0.0);
          // the map is the evaluator's other half and moves with the values:
          // assign never shrinks capacity, so a draw that SHRINKS the arena
          // would leave leafOf naming node ids past mu's .size() and the
          // gather would read stale bytes rather than fault. All-root is in
          // bounds for every tree (nodes.size() >= 1) and gathers mu[0] == 0.
          std::memset(forest.leafOf.data() + t * n, 0,
                      n * sizeof(std::uint32_t));
          // still marked for rebuild at the tree's own draw, which is also
          // what declines the fused suffstat pass for this sweep
          forest.leafOfStale[t] = 1;
        } else if constexpr (L::hasVectorParams || L::hasFunctionParams) {
          // vector and function leaves carry the dense slab in place of the
          // (mu, leafOf) pair; for the function leaf the row IS the parameters
          misc_setVectorToConstant(forest.treeFits.data() + t * n, n, 0.0);
        }
      }
      misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);
      if (data_.numTestObservations > 0)
        misc_setVectorToConstant(forest.totalTestFits.data(),
                                 data_.numTestObservations, 0.0);
    }
  }

  /// Warm-start producer: run numSweeps of XBART-style grow-from-root in place
  /// of the MH move, leaving a legal forest the exact sweeps then own. This is
  /// a SEPARATE sweep loop that DUPLICATES run()'s per-iteration body (the
  /// residual roll, node averages, sampleParametersAndSetFits, and the
  /// latent/sigma/k/DART updates) rather than branching inside it, so the
  /// default run path stays byte-identical; the only substitution is the
  /// per-tree step, where growTreeFromRoot rebuilds the tree from a fresh root
  /// against the tree's residual. Constant leaf only (the scan's marginal);
  /// vector and function leaves are refused on the R surface and no-op here.
  void growForestFromRoot(size_t numSweeps) {
    if constexpr (L::hasVectorParams || L::hasFunctionParams) {
      (void) numSweeps;  // non-constant leaves fall back to prior-grown init
      return;
    } else {
      size_t n = data_.numObservations;
      double* y = response_->workingResponse();
      const double* weights = response_->workingWeights();
      GrowScratch growScratch;

      for (size_t sweep = 0; sweep < numSweeps; ++sweep) {
        // heteroscedastic: the scan reads precisions, and under a variance
        // forest the global sigma is pinned at 1 on the working scale, so it
        // must see w_i / s^2(x_i) - run()'s own pre-step. Scanning against unit
        // precisions instead would price every split against a residual
        // variance of 1 where the data's is s^2(x_i). The variance forest is
        // NOT swept here: grow-from-root initializes the mean forest, and the
        // following run's first sweep fits the variance surface against it.
        if (varianceForest_) {
          formMeanWeights();
          weights = meanWeights_.data();
        }

        for (size_t f = 0; f < forests_.size(); ++f) {
          Forest<L, ResidT>& forest = forests_[f];
          const double* forestY = y;
          const double* forestWeights = weights;
          if (combiner_) {
            combiner_->drawForestGlue(f, rng_, forests_);
            ForestResponse fr = combiner_->formForestResponse(f, forests_, y,
                                                              weights);
            forestY = fr.response;
            forestWeights = composeForestWeights(f, fr.weights);
          }

          forest.kSumSquaredParams = 0.0;
          forest.kNumLeaves = 0.0;

          for (size_t t = 0; t < forest.numTrees; ++t) {
            rollTreeResidual(forest, t, forestY);

            // grow a fresh tree from the root against tree t's residual, in
            // place of metropolisJumpForTree; the reset returns it to a single
            // root over the full index buffer, then setNodeAverages primes the
            // root statistic the scan's no-split term reads
            forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
            forest.trees[t].setNodeAverages(forest.treeY.data(), forestWeights);
            growTreeFromRoot(data_, forest.treePrior, forest.leaf, rng_,
                             forest.trees[t], 0, forest.treeY.data(),
                             forestWeights, forest.k, sigma_, growScratch);
            if (data_.hasPooledCategorical)
              forest.trees[t].compactMaskPoolIfNeeded(data_);

            // regrowth replaces the partition wholesale
            rebuildLeafOf(forest, t);
            sampleParametersAndSetFits(forest, t, nullptr, false);
          }

          finalizeTotalFits(forest, forestY);
        }

        const double* combined = combinedFits();
        response_->refreshLatents(rng_, combined, sigma_);
        y = response_->workingResponse();
        weights = response_->workingWeights();

        if (!sigmaIsFixed_)
          sigma_ = response_->drawSigma(rng_, combined, sigma_);

        if (combiner_) {
          combiner_->drawGlue(rng_, sigma_, y, weights, forests_);
          combiner_->afterCombine(forests_, false, 0, rng_);
        }

        for (Forest<L, ResidT>& forest : forests_) {
          if (forest.updateK && forest.kSumSquaredParams > 0.0)
            forest.k = forest.kHyperprior.draw(rng_, forest.kSumSquaredParams,
                                               forest.kNumLeaves,
                                               forest.leaf.scale);
          if (forest.useDart) {
            std::memset(forest.splitCounts.data(), 0,
                        forest.splitCounts.size() * sizeof(std::uint32_t));
            for (size_t t = 0; t < forest.numTrees; ++t)
              forest.trees[t].countVariableUses(forest.splitCounts.data());
            forest.dart.update(rng_, forest.splitCounts.data());
          }
        }
      }
    }
  }

  /// Replace every leaf parameter with a draw from the node prior and
  /// rebuild the tree, total, and test fits to match.
  void sampleNodeParametersFromPrior() {
    size_t n = data_.numObservations;
    for (Forest<L, ResidT>& forest : forests_) {
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
          // the draws become the live leaf block, so under a constraint they
          // must come from the prior truncated to the monotone cone - drawn
          // per tree, since the cone couples a tree's leaves
          if constexpr (TreeDrawLeafModel<L>) {
            if (!forest.leaf.drawFromPriorForTree(rng_, tree, tree.bottomScratch,
                                                  forest.k,
                                                  forest.paramByNode.data()))
              ext_throwError("monotone prior draw: no feasible leaf vector in "
                             "%d attempts", L::priorDrawMaxAttempts);
          } else {
            for (int32_t i : tree.bottomScratch)
              forest.paramByNode[static_cast<size_t>(i)] =
                forest.leaf.drawFromPrior(rng_, forest.k);
          }

          setTreeFitsFromParameters(forest, t, forest.paramByNode);
          if (forest.leafOfStale[t] != 0) rebuildLeafOf(forest, t);
          addTreeFitsToTotal(forest, t);
          routeTestRows(data_.numTestObservations, [&](size_t i) {
            int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
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
            int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
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
            int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
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
  /// One entry per forest of forests_, in forest order.
  using ForestParameters = std::vector<TreeParameters>;

  /// What the validate phase of a predictor transaction hands the rebuild
  /// phase: every forest's recovered leaf parameters, plus the per-forest list
  /// of trees the two phases agree to visit. Forest 0's list is always every
  /// tree in order - its subtract-then-add round trip is the recorded
  /// single-forest arithmetic, and (T - f) + f != T at the ULP, so pruning it
  /// would move a baseline - while a forest past it lists only the trees that
  /// split on a column the transaction touched. Skipping the rest is exact,
  /// not approximate: an untouched tree's partition, parameters and fits are
  /// unchanged, and leaving its contribution to forest.totalFits alone
  /// preserves it bitwise where a round trip would perturb it
  /// (docs/plans/multiforest-predictor-mutation.md, "Pruning"). An empty
  /// params[f][t] is NOT a legal skip marker - a function-valued leaf
  /// legitimately produces one - so the list is explicit and shared.
  ///
  /// The variance forest rides the same handoff in its own pair of members: it
  /// lives outside forests_, carries node-indexed leaf FACTORS rather than mean
  /// parameters, and prunes on the same predicate a forest past the first does.
  struct ForestRevalidation {
    ForestParameters params;
    std::vector<std::vector<std::size_t>> survivors;
    TreeParameters varianceParams;
    std::vector<std::size_t> varianceSurvivors;
  };

  /// Leaf parameters of tree t in transferable form: recovered from the
  /// fits for scalar leaves, copied from the persisted blocks for vector
  /// ones; function-valued leaves keep their per-observation fits in place,
  /// so nothing is recovered.
  void recoverLeafParameters(Forest<L, ResidT>& forest, size_t t,
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
  void setTreeFits(Forest<L, ResidT>& forest, size_t t,
                   const std::vector<double>& params) {
    static_assert(!L::hasFunctionParams);
    if constexpr (!L::hasVectorParams) {
      setTreeFitsFromParameters(forest, t, params);
    } else {
      forest.paramsByTree[t] = params;
      setTreeFitsFromParameterBlocks(forest, t, params);
    }
  }

  /// THE pruning predicate: forest's trees that split on at least one of the
  /// numTouched columns, in tree order. Both consumers share it - the
  /// transaction's revalidate/rebuild pair through collectSurvivors, which
  /// exempts forest 0, and the per-observation session's cell guard through
  /// treesSplittingOnColumn, which does not. scratch is a numPredictors-sized
  /// census buffer the caller owns so the census does not allocate per tree.
  /// Templated on the ensemble so the variance forest, which is not a
  /// Forest<L, ResidT>, prunes on the same predicate rather than a copy of it.
  template <typename Ensemble>
  void collectSplittingTrees(const Ensemble& forest,
                             const std::size_t* touched, std::size_t numTouched,
                             std::vector<std::uint32_t>& scratch,
                             std::vector<std::size_t>& out) const {
    out.clear();
    scratch.resize(data_.numPredictors);
    for (size_t t = 0; t < forest.numTrees; ++t) {
      std::memset(scratch.data(), 0,
                  data_.numPredictors * sizeof(std::uint32_t));
      forest.trees[t].countVariableUses(scratch.data());
      for (size_t k = 0; k < numTouched; ++k)
        if (scratch[touched[k]] != 0) {
          out.push_back(t);
          break;
        }
    }
  }

  /// Forest f's survivor list for a transaction touching numTouched columns
  /// (a null list means every column, where the predicate is vacuous):
  /// forest 0 keeps every tree, a forest past it keeps the trees that split on
  /// a touched column.
  void collectSurvivors(const Forest<L, ResidT>& forest, std::size_t f,
                        const std::size_t* touched, std::size_t numTouched,
                        std::vector<std::uint32_t>& scratch,
                        std::vector<std::size_t>& out) const {
    if (f == 0 || touched == nullptr) {
      out.clear();
      out.resize(forest.numTrees);
      for (size_t t = 0; t < forest.numTrees; ++t) out[t] = t;
      return;
    }
    collectSplittingTrees(forest, touched, numTouched, scratch, out);
  }

  /// Forest f's trees that split on column: what a per-observation session
  /// caches for that forest. The same predicate the rebuild prunes with,
  /// WITHOUT forest 0's exemption - that exemption preserves forest 0's
  /// recorded subtract-then-add arithmetic and is no part of the predicate.
  /// So the guarded set is a subset of the revalidated set in every forest,
  /// and a tree in the difference is one whose partition the revalidation
  /// reproduces unchanged from an occupied pre-state: that is why the session's
  /// finalize cannot fail (docs/plans/multiforest-predictor-mutation.md,
  /// "Pruning"). Dropping a tree from the cache cannot move the returned mask
  /// either: with no node splitting on column the override never fires during
  /// the descent, so the new leaf is the old one and no move is ever staged.
  void treesSplittingOnColumn(std::size_t f, std::size_t column,
                              std::vector<std::uint32_t>& scratch,
                              std::vector<std::size_t>& out) const {
    collectSplittingTrees(forests_[f], &column, 1, scratch, out);
  }

  /// The variance forest's own pair. It sits outside forests_ and has no
  /// numTreesInForest entry, so it is addressed by numVarianceTrees() and
  /// reached through these rather than through a forest index. The pruning is
  /// the f >= 1 rule (no forest-0 exemption): the variance path carries no
  /// recorded subtract-then-add arithmetic to preserve, so skipping a tree no
  /// touched column reaches is the exact implementation, not an optimization.
  void varianceSurvivors(const std::size_t* touched, std::size_t numTouched,
                         std::vector<std::uint32_t>& scratch,
                         std::vector<std::size_t>& out) const {
    const VarianceForest& vf = *varianceForest_;
    if (touched == nullptr) {
      out.clear();
      out.resize(vf.numTrees);
      for (std::size_t j = 0; j < vf.numTrees; ++j) out[j] = j;
      return;
    }
    collectSplittingTrees(vf, touched, numTouched, scratch, out);
  }
  void varianceTreesSplittingOnColumn(std::size_t column,
                                      std::vector<std::uint32_t>& scratch,
                                      std::vector<std::size_t>& out) const {
    collectSplittingTrees(*varianceForest_, &column, 1, scratch, out);
  }

  /// Recover leaf parameters (from fits for scalar leaves, from the
  /// persisted blocks for vector ones; function-valued leaves keep their
  /// per-observation fits in place, so nothing is recovered), re-route every
  /// tree of every forest against the store's current codes, and report
  /// whether all leaves stay occupied. touched names the columns the
  /// transaction moved (null for a whole-matrix swap); it drives the j-split
  /// pruning of forests past the first, which forest 0 does not take. A
  /// single-forest chain therefore runs exactly the loop it always ran.
  bool revalidateTrees(ForestRevalidation& state, const std::size_t* touched,
                       std::size_t numTouched) {
    state.params.resize(forests_.size());
    state.survivors.resize(forests_.size());
    std::vector<std::uint32_t> census;
    bool allValid = true;
    for (size_t f = 0; f < forests_.size() && allValid; ++f) {
      Forest<L, ResidT>& forest = forests_[f];
      collectSurvivors(forest, f, touched, numTouched, census,
                       state.survivors[f]);
      TreeParameters& params = state.params[f];
      params.resize(forest.numTrees);
      const std::vector<std::size_t>& survivors = state.survivors[f];
      for (size_t k = 0; k < survivors.size() && allValid; ++k) {
        size_t t = survivors[k];
        recoverLeafParameters(forest, t, params[t]);
        forest.trees[t].repartitionSubtree(data_, 0);
        allValid = forest.trees[t].bottomNodesAreOccupied();
      }
    }
    // the variance forest, appended after the forests_ body: a homoscedastic
    // chain never takes the branch, so its path is the one it always ran
    if (allValid && varianceForest_ != nullptr)
      allValid = revalidateVarianceTrees(state, touched, numTouched, census);
    return allValid;
  }

  /// The variance forest's validate phase: recover each surviving tree's
  /// node-indexed leaf factors through its LIVE partition FIRST, then
  /// repartition against the store's current codes and report occupancy. The
  /// order is the one refreshVarianceForest documents - the recovery reads the
  /// per-observation slab through each leaf's current members, so recovering
  /// after the repartition would read the new partition's members out of the
  /// old leaves' slots.
  ///
  /// Collapses nothing, drops no missing directions and scatters nothing: this
  /// half must be undoable by a repartition alone, and factorByTree and
  /// combinedVariance stay untouched so a rollback restores the state exactly.
  bool revalidateVarianceTrees(ForestRevalidation& state,
                               const std::size_t* touched,
                               std::size_t numTouched,
                               std::vector<std::uint32_t>& census) {
    VarianceForest& vf = *varianceForest_;
    varianceSurvivors(touched, numTouched, census, state.varianceSurvivors);
    TreeParameters& params = state.varianceParams;
    params.resize(vf.numTrees);
    bool allValid = true;
    const std::vector<std::size_t>& survivors = state.varianceSurvivors;
    for (std::size_t k = 0; k < survivors.size() && allValid; ++k) {
      std::size_t j = survivors[k];
      recoverVarianceLeafValues(vf, j, params[j]);
      vf.trees[j].repartitionSubtree(data_, 0);
      allValid = vf.trees[j].bottomNodesAreOccupied();
    }
    return allValid;
  }

  /// Second phase of a successful transaction: rewrite tree fits from the
  /// parameters revalidateTrees recovered, over the same survivor lists.
  /// Node averages are left stale; run() recomputes them from current
  /// residuals before any use. Function-valued leaves only refresh the
  /// covariate gather: their per-observation fits are the parameters and stay
  /// in place (the next sweep's draws replace them under the new values).
  void rebuildFitsFromParameters(const ForestRevalidation& state) {
    dropStaleMissingDirections();
    for (size_t f = 0; f < forests_.size(); ++f) {
      Forest<L, ResidT>& forest = forests_[f];
      const TreeParameters& params = state.params[f];
      const std::vector<std::size_t>& survivors = state.survivors[f];
      if constexpr (L::hasFunctionParams) {
        (void) params;
        (void) survivors;
        forest.leaf.regatherTrainingCovariates(data_);
      } else {
        // vector leaves read raw covariate values: pick up the installed ones
        if constexpr (L::hasVectorParams)
          forest.leaf.regatherTrainingCovariates(data_);
        for (size_t k = 0; k < survivors.size(); ++k) {
          size_t t = survivors[k];
          // the subtract reads the pre-reroute (mu, leafOf) pair, exactly the
          // cached fits totalFits still sums; the map then tracks the
          // repartition revalidateTrees performed
          subtractTreeFitsFromTotal(forest, t);
          setTreeFits(forest, t, params[t]);
          if constexpr (leafIsConstant) {
            installLeafOfAndAddToTotal(forest, t);
          } else {
            addTreeFitsToTotal(forest, t);
          }
        }
      }
    }
    if (varianceForest_) rebuildVarianceFactors(state);
  }

  /// The variance forest's rebuild phase: drop the missing directions the new
  /// codes stranded, scatter each surviving tree's recovered factors through
  /// the partition the validate phase installed, and recompute s^2(x) as the
  /// fresh product. A pruned tree keeps its slab entries, which its unchanged
  /// partition still describes, so the product is the same one it would have
  /// been round tripped to.
  ///
  /// The direction drop follows the repartition here rather than preceding it
  /// (refreshVarianceForest's order): with hasMissing false for the column, the
  /// bit routes nothing, so clearing it before or after the routing is the same
  /// partition (tree.hpp, dropStaleMissingDirections). It runs over every tree,
  /// as the mean side's chain-level drop does, since a stale bit outside the
  /// reachable gauge would fail a later flatten wherever it sits.
  void rebuildVarianceFactors(const ForestRevalidation& state) {
    VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    for (std::size_t j = 0; j < vf.numTrees; ++j)
      vf.trees[j].dropStaleMissingDirections(data_);
    for (std::size_t k = 0; k < state.varianceSurvivors.size(); ++k) {
      std::size_t j = state.varianceSurvivors[k];
      Tree& tree = vf.trees[j];
      const std::vector<double>& leafValues = state.varianceParams[j];
      double* factor = vf.factorByTree.data() + j * n;
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      for (std::int32_t b : tree.bottomScratch) {
        double h = leafValues[static_cast<std::size_t>(b)];
        const Node& node = tree.at(b);
        for (std::size_t m = node.begin; m < node.end; ++m)
          factor[tree.indices[m]] = h;
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      double product = 1.0;
      for (std::size_t j = 0; j < vf.numTrees; ++j)
        product *= vf.factorByTree[j * n + i];
      vf.combinedVariance[i] = product;
    }
  }

  /// Rollback re-route: restore partitions consistent with the store after
  /// the sampler has restored its old codes.
  void repartitionTrees() {
    for (Forest<L, ResidT>& forest : forests_) {
      // the restored raw values re-gather to exactly the old covariates
      if constexpr (L::hasVectorParams || L::hasFunctionParams)
        forest.leaf.regatherTrainingCovariates(data_);
      for (size_t t = 0; t < forest.numTrees; ++t)
        forest.trees[t].repartitionSubtree(data_, 0);
    }
    // the variance trees the validate phase re-routed. Load-bearing: without
    // it a rejected transaction leaves s^2(x) routed by the proposal. Nothing
    // else has to be put back - the validate phase leaves factorByTree and
    // combinedVariance alone, so restoring the partition restores the state.
    if (varianceForest_)
      for (Tree& tree : varianceForest_->trees)
        tree.repartitionSubtree(data_, 0);
  }

  /// First phase of a whole-data replacement: recover every tree's leaf
  /// parameters against the current fits and partitions, before the store or
  /// any per-observation storage moves.
  void recoverTreeParameters(TreeParameters& params) {
    Forest<L, ResidT>& forest = forests_[0];
    params.resize(forest.numTrees);
    for (size_t t = 0; t < forest.numTrees; ++t)
      recoverLeafParameters(forest, t, params[t]);
  }

  /// Companion first phase for the variance forest: the node-indexed leaf
  /// factors of every variance tree, read through the LIVE partition. The
  /// recovery strides the per-observation slab by data_.numObservations, so it
  /// cannot wait until applyNewData - by then the store carries the
  /// replacement count and the stride is wrong. Left empty (and unread) off a
  /// variance forest.
  void recoverVarianceParameters(TreeParameters& params) {
    if (varianceForest_ == nullptr) return;
    VarianceForest& vf = *varianceForest_;
    params.resize(vf.numTrees);
    for (std::size_t j = 0; j < vf.numTrees; ++j)
      recoverVarianceLeafValues(vf, j, params[j]);
  }

  /// Second phase, after the shared store holds the new predictors and
  /// freshly rebuilt cuts: swap the response state (sigma is preserved on
  /// the original scale), resize per-observation storage, remap split
  /// indices onto the new cut grid, re-route, and collapse anything left
  /// invalid or empty. Node averages are left stale; run() recomputes them.
  /// varianceParams carries recoverVarianceParameters' factors and is read
  /// only under a variance forest.
  void applyNewData(const double* y, const double* weights,
                    const double* offset,
                    const std::vector<std::vector<double>>& oldCutPoints,
                    TreeParameters& params,
                    const TreeParameters& varianceParams) {
    Forest<L, ResidT>& forest = forests_[0];
    size_t n = data_.numObservations;
    bool numObservationsChanged =
      n * forest.numTrees != forest.indexBuffer.size();

    weights_ = weights;
    response_->setData(y, offset, weights, n, &sigma_);

    if (numObservationsChanged) {
      forest.indexBuffer.resize(n * forest.numTrees);
      if constexpr (leafIsConstant)
        forest.leafOf.assign(n * forest.numTrees, 0);
      else
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
      if constexpr (leafIsConstant) rebuildLeafOf(forest, t);
      addTreeFitsToTotal(forest, t);
    }

    resizeTestStorage();

    // the variance forest lives outside forests_ and pins seven n-sized
    // allocations at the creation count, so a replacement data set of a
    // different length overruns them (meanWeights_ first, in formMeanWeights).
    // Resize, then re-anchor through the rebuilt grid with the factors
    // recovered before the store moved. Appended after the forest-0 body:
    // nothing above it moves, so a homoscedastic chain draws what it drew.
    if (varianceForest_) {
      bool varianceCountChanged = varianceForest_->numObservations() != n;
      if (varianceCountChanged) resizeVarianceStorage(n);
      refreshVarianceForest(&oldCutPoints, &varianceParams,
                            varianceCountChanged);
    }
  }

  /// After a data mutation re-quantizes the store, drop every tree's stale
  /// missing directions so the live masks stay within the reachable gauge.
  void dropStaleMissingDirections() {
    for (Forest<L, ResidT>& forest : forests_)
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

    for (Forest<L, ResidT>& forest : forests_) {
      misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);

      if constexpr (L::hasFunctionParams) {
        forest.leaf.regatherTrainingCovariates(data_);
        std::vector<double> dummyParams;
        for (size_t t = 0; t < forest.numTrees; ++t) {
          forest.trees[t].repartitionSubtree(data_, 0);
          dummyParams.assign(forest.trees[t].nodes.size(), 0.0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             dummyParams);
          addTreeFitsToTotal(forest, t);
        }
      } else if constexpr (!L::hasVectorParams) {
        std::vector<double> paramByNode;
        for (size_t t = 0; t < forest.numTrees; ++t) {
          recoverParametersFromFits(forest, t, paramByNode);
          forest.trees[t].repartitionSubtree(data_, 0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             paramByNode);
          setTreeFitsFromParameters(forest, t, paramByNode);
          if constexpr (leafIsConstant) rebuildLeafOf(forest, t);
          addTreeFitsToTotal(forest, t);
        }
      } else {
        forest.leaf.regatherTrainingCovariates(data_);
        size_t numParams = forest.leaf.numParams();
        for (size_t t = 0; t < forest.numTrees; ++t) {
          forest.trees[t].repartitionSubtree(data_, 0);
          forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                             forest.paramsByTree[t], numParams);
          setTreeFitsFromParameterBlocks(forest, t, forest.paramsByTree[t]);
          addTreeFitsToTotal(forest, t);
        }
      }
    }

    // the variance forest lives outside forests_ and needs the same re-route;
    // the grid is installed rather than rebuilt on every path that reaches
    // here (setCutPoints, a forced predictor swap), so no remap. Appended, and
    // called per site rather than from dropStaleMissingDirections, which three
    // paths share (refreshVarianceForest drops its own directions).
    if (varianceForest_) refreshVarianceForest(nullptr);
  }

  // Saved-tree (keepTrees) storage: a circular buffer of capacity slots,
  // each one kept sample's forest in flattened form. The slot base is set by
  // the sampler before every run so chains write consistent slots without
  // sharing mutable state.

  void initializeSavedTrees(size_t capacity) {
    if (varianceForest_) {
      varianceForest_->savedTreeCapacity = capacity;
      varianceForest_->savedTrees.assign(
        capacity * varianceForest_->numTrees, std::vector<FlatNode>(1));
    }
    for (Forest<L, ResidT>& forest : forests_) {
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
    for (Forest<L, ResidT>& forest : forests_) forest.savedSlotBase = base;
    if (varianceForest_) varianceForest_->savedSlotBase = base;
  }
  size_t savedTreeCapacity() const { return forests_[0].savedTreeCapacity; }
  const std::vector<FlatNode>& savedTree(size_t slot, size_t t,
                                         size_t forestIndex = 0) const {
    const Forest<L, ResidT>& forest = forests_[forestIndex];
    return forest.savedTrees[slot * forest.numTrees + t];
  }
  /// Slopes of one saved tree (vector-parameter leaves), parallel to
  /// savedTree's pre-order leaves.
  const std::vector<double>& savedTreeSlopes(size_t slot, size_t t,
                                             size_t forestIndex = 0) const {
    const Forest<L, ResidT>& forest = forests_[forestIndex];
    return forest.savedTreeParams[slot * forest.numTrees + t];
  }
  /// Flattened mask words of one saved tree (wide categorical columns).
  const std::vector<std::uint64_t>& savedTreeMasks(
    size_t slot, size_t t, size_t forestIndex = 0) const {
    const Forest<L, ResidT>& forest = forests_[forestIndex];
    return forest.savedTreeMasks[slot * forest.numTrees + t];
  }

  /// Flatten live tree t into saved slot `slot` while its freshly drawn
  /// parameters are live; treeFits is the tree's slab. Function-valued
  /// leaves' records carry per-leaf mean fits for reporting, and their side
  /// channel the draw cache's alpha weights plus covariate rows - the exact
  /// values the recorded test fits used, so saved replays bit-match them.
  void storeSavedTreeRecord(Forest<L, ResidT>& forest, size_t t, size_t slot,
                            const double* treeFits) {
    std::vector<std::uint64_t>* masks = data_.hasPooledCategorical
      ? &forest.savedTreeMasks[slot * forest.numTrees + t] : nullptr;
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      forest.trees[t].flatten(data_, forest.muByTree[t].data(),
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
                   std::vector<std::uint64_t>* masks = nullptr,
                   size_t forestIndex = 0) {
    Forest<L, ResidT>& forest = forests_[forestIndex];
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
    Forest<L, ResidT>& forest = forests_[0];
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

  /// The same for one saved tree; function-valued leaves print their
  /// recorded mean values.
  void printSavedTree(size_t slot, size_t t, int indentation) const {
    const Forest<L, ResidT>& forest = forests_[0];
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

  /// Adds one flattened tree's predictions for numTestObservations rows of a
  /// Columns predictor source to out, dispatching on the leaf shape's record
  /// format: plain leaf values, slope blocks, or function blocks (whose
  /// offsets are valid by construction or by stateIsValid). sideChannel is
  /// null for scalar leaves; indices and blockOffsets are caller scratch.
  template <typename Columns>
  void addFlatPredictions(const std::vector<FlatNode>& flat,
                          const std::vector<double>* sideChannel,
                          const std::uint64_t* masks, const Columns& x_test,
                          size_t numTestObservations,
                          std::vector<size_t>& indices,
                          std::vector<size_t>& blockOffsets,
                          double* out) const {
    const L& leaf = forests_[0].leaf;
    for (size_t i = 0; i < numTestObservations; ++i) indices[i] = i;
    if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
      addFlatPredictionsBelow(flat.data(), x_test, indices.data(), 0,
                              numTestObservations, out, masks);
    } else if constexpr (L::hasVectorParams) {
      addFlatLinearPredictionsBelow(
        flat.data(), x_test,
        indices.data(), 0, numTestObservations, out,
        leaf.covariateColumns().data(), leaf.covariateMeans().data(),
        leaf.covariateSds().data(), leaf.numParams() - 1,
        sideChannel->data(), 0, masks);
    } else {
      computeFunctionBlockOffsets(sideChannel->data(), sideChannel->size(),
                                  (flat.size() + 1) / 2,
                                  leaf.numCovariates(), blockOffsets);
      addFlatFunctionPredictionsBelow(
        flat.data(), x_test,
        indices.data(), 0, numTestObservations, out,
        leaf.covariateColumns().data(), leaf.covariateMeans().data(),
        leaf.covariateSds().data(), leaf.lengthscales().data(),
        leaf.numCovariates(), sideChannel->data(), blockOffsets.data(), 0,
        masks);
    }
  }

  /// Fits for the test rows of a Columns predictor source from one saved
  /// sample's trees, on the original response scale; offsets are the caller's
  /// problem.
  template <typename Columns>
  void predictFromSavedSample(size_t slot, const Columns& columns,
                              size_t numTestObservations, double* out) const {
    const Forest<L, ResidT>& forest = forests_[0];
    misc_setVectorToConstant(out, numTestObservations, 0.0);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    for (size_t t = 0; t < forest.numTrees; ++t) {
      const std::uint64_t* masks = data_.hasPooledCategorical
        ? forest.savedTreeMasks[slot * forest.numTrees + t].data() : nullptr;
      const std::vector<double>* sideChannel = forest.savedTreeParams.empty()
        ? nullptr : &forest.savedTreeParams[slot * forest.numTrees + t];
      addFlatPredictions(forest.savedTrees[slot * forest.numTrees + t],
                         sideChannel, masks, columns, numTestObservations,
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
  template <typename Columns>
  void predictFromCurrentTrees(const Columns& columns,
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
      addFlatPredictions(flat, &slopes, maskBuffer.data(), columns,
                         numTestObservations, indices, blockOffsets, out);
    }
    double scale = response_->fitScale();
    double shift = response_->fitShift();
    for (size_t i = 0; i < numTestObservations; ++i)
      out[i] = scale * out[i] + shift;
  }

  /// K-forest softmax replay of one saved sample (multinomial): sum every one
  /// of the K forests' saved trees at the new rows into a location-major raw
  /// slab, then softmax per row into out (nTest x K, out[k*nTest + i]) through
  /// the same map storeSample's test channel applies to totalTestFits. The
  /// level-centering grand shift is absent from the saved (flattened) leaves,
  /// but softmax is invariant to a shift common to all K categories, so the
  /// replayed probabilities are the identified ones (docs/design/multinomial.md).
  /// The per-forest total is on the internal (softmax log-odds) scale, not the
  /// fitScale-shifted response scale, exactly as totalTestFits is - fitScale is
  /// the identity for the multinomial response, so no conversion applies.
  template <typename Columns>
  void predictFromSavedSampleMulti(size_t slot, const Columns& columns,
                                   size_t numTestObservations,
                                   double* out) const {
    size_t K = forests_.size();
    std::vector<double> raw(numTestObservations * K);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    for (size_t f = 0; f < K; ++f) {
      const Forest<L, ResidT>& forest = forests_[f];
      double* forestRaw = raw.data() + f * numTestObservations;
      misc_setVectorToConstant(forestRaw, numTestObservations, 0.0);
      for (size_t t = 0; t < forest.numTrees; ++t) {
        const std::uint64_t* masks = data_.hasPooledCategorical
          ? forest.savedTreeMasks[slot * forest.numTrees + t].data() : nullptr;
        const std::vector<double>* sideChannel = forest.savedTreeParams.empty()
          ? nullptr : &forest.savedTreeParams[slot * forest.numTrees + t];
        addFlatPredictions(forest.savedTrees[slot * forest.numTrees + t],
                           sideChannel, masks, columns, numTestObservations,
                           indices, blockOffsets, forestRaw);
      }
    }
    softmaxLocationMajor(raw.data(), numTestObservations, K, out);
  }

  /// The same K-forest softmax replay from the live trees, flattened on the
  /// fly; reached only when keepTrees is off, which the R surface refuses for a
  /// multinomial predict, so it exists for completeness and the engine tests.
  template <typename Columns>
  void predictFromCurrentTreesMulti(const Columns& columns,
                                    size_t numTestObservations, double* out) {
    size_t K = forests_.size();
    std::vector<double> raw(numTestObservations * K);
    std::vector<size_t> indices(numTestObservations);
    std::vector<size_t> blockOffsets;
    std::vector<double> slopes;
    std::vector<std::uint32_t> counts;
    std::vector<FlatNode> flat;
    std::vector<std::uint64_t> maskBuffer;
    std::vector<std::uint64_t>* masks =
      data_.hasPooledCategorical ? &maskBuffer : nullptr;
    for (size_t f = 0; f < K; ++f) {
      double* forestRaw = raw.data() + f * numTestObservations;
      misc_setVectorToConstant(forestRaw, numTestObservations, 0.0);
      for (size_t t = 0; t < forests_[f].numTrees; ++t) {
        flattenTree(t, flat, counts, &slopes, masks, f);
        addFlatPredictions(flat, &slopes, maskBuffer.data(), columns,
                           numTestObservations, indices, blockOffsets,
                           forestRaw);
      }
    }
    softmaxLocationMajor(raw.data(), numTestObservations, K, out);
  }

  // Chain state serialization. getState captures everything the posterior
  // state comprises; stateIsValid checks a candidate against the store's
  // current cuts without mutating anything, so a multi-chain restore can be
  // all-or-none; setState trusts that check. Vector-parameter leaves carry
  // their slopes in treeParams/savedTreeParams alongside the flat trees.

  void getState(ChainStateData& state) {
    state.forests.resize(forests_.size());
    for (size_t f = 0; f < forests_.size(); ++f) {
      Forest<L, ResidT>& forest = forests_[f];
      ForestStateData& fs = state.forests[f];
      fs.trees.resize(forest.numTrees);
      if (data_.hasPooledCategorical) {
        fs.treeMasks.resize(forest.numTrees);
        fs.savedTreeMasks = forest.savedTreeMasks;
      } else {
        fs.treeMasks.clear();
        fs.savedTreeMasks.clear();
      }
      if constexpr (!L::hasVectorParams && !L::hasFunctionParams) {
        std::vector<double> params;
        for (size_t t = 0; t < forest.numTrees; ++t) {
          recoverParametersFromFits(forest, t, params);
          forest.trees[t].flatten(data_, params.data(), fs.trees[t], nullptr, 1,
                                  nullptr,
                                  data_.hasPooledCategorical ? &fs.treeMasks[t]
                                                           : nullptr);
        }
        fs.treeParams.clear();
        fs.savedTreeParams.clear();
      } else if constexpr (L::hasVectorParams) {
        fs.treeParams.resize(forest.numTrees);
        for (size_t t = 0; t < forest.numTrees; ++t)
          forest.trees[t].flatten(data_, forest.paramsByTree[t].data(),
                                  fs.trees[t], nullptr, forest.leaf.numParams(),
                                  &fs.treeParams[t],
                                  data_.hasPooledCategorical ? &fs.treeMasks[t]
                                                           : nullptr);
        fs.savedTreeParams = forest.savedTreeParams;
      } else {
        // function-valued leaves: records carry reporting means, and each
        // live tree's parameters ARE its fits - one slab per tree in
        // observation order, restored by copy
        size_t n = data_.numObservations;
        std::vector<double> values;
        fs.treeParams.resize(forest.numTrees);
        for (size_t t = 0; t < forest.numTrees; ++t) {
          functionLeafValues(forest.trees[t], forest.treeFits.data() + t * n,
                             values);
          forest.trees[t].flatten(data_, values.data(), fs.trees[t], nullptr, 1,
                                  nullptr,
                                  data_.hasPooledCategorical ? &fs.treeMasks[t]
                                                           : nullptr);
          fs.treeParams[t].assign(forest.treeFits.data() + t * n,
                                  forest.treeFits.data() + (t + 1) * n);
        }
        fs.savedTreeParams = forest.savedTreeParams;
      }
      fs.savedTrees = forest.savedTrees;
      fs.k = forest.k;
      // written for EVERY forest, not just the response-derived ones (BCF's):
      // the block is self-describing, and a data-independent scale simply
      // records the value a same-spec destination already constructed
      fs.leafScale = forest.leaf.scale;
    }
    Forest<L, ResidT>& forest = forests_[0];
    state.sigma = sigma();
    response_->getScale(state.fitMin, state.fitMax);
    if (response_->latents() != nullptr) {
      state.latents.assign(response_->latents(),
                           response_->latents() + data_.numObservations);
    } else {
      state.latents.clear();
    }
    // a t response's lambda rides latents above; its nu is the companion
    // scalar block, absent (NaN) for every other family
    state.residualDf = response_->carriesResidualDf()
                         ? response_->residualDf()
                         : std::numeric_limits<double>::quiet_NaN();
    // the ordinal-only cutpoint vector (length K-1); z rides latents above. A
    // non-ordinal chain carries none and writes no block.
    if (response_->carriesCutpoints()) {
      state.cutpoints.assign(
        response_->cutpoints(),
        response_->cutpoints() + response_->numCutpoints());
    } else {
      state.cutpoints.clear();
    }
    // an NB response's dispersion r is a companion scalar block (the resid.df
    // pattern); omega rides latents above. Absent (NaN) for every other family.
    state.dispersion = response_->carriesDispersion()
                         ? response_->dispersion()
                         : std::numeric_limits<double>::quiet_NaN();
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
    // the combiner fills the BCF-shaped glue wire format (hasBCF included); a
    // single-forest chain carries no combiner and leaves it off
    state.hasBCF = false;
    if (combiner_) combiner_->serializeGlue(state);

    // heteroscedastic: flatten every variance tree with its recovered positive
    // leaf factors (working scale); empty off a variance forest
    if (varianceForest_) {
      VarianceForest& vf = *varianceForest_;
      state.varianceTrees.resize(vf.numTrees);
      std::vector<double> leafValues;
      for (std::size_t j = 0; j < vf.numTrees; ++j) {
        recoverVarianceLeafValues(vf, j, leafValues);
        vf.trees[j].flatten(data_, leafValues.data(), state.varianceTrees[j]);
      }
    } else {
      state.varianceTrees.clear();
    }
  }

  bool stateIsValid(const ChainStateData& state) const {
    if (state.forests.size() != forests_.size()) return false;
    size_t n = data_.numObservations;
    Tree scratch;
    std::vector<index_t> scratchIndices(n);
    std::vector<double> params;
    for (size_t f = 0; f < forests_.size(); ++f) {
      const Forest<L, ResidT>& forest = forests_[f];
      const ForestStateData& fs = state.forests[f];
      if (fs.trees.size() != forest.numTrees) return false;
      if (!fs.savedTrees.empty() &&
          fs.savedTrees.size() != forest.savedTrees.size())
        return false;
      // mask channels pair with their flat trees when present; trees holding
      // wide rules without a channel fail the rebuild below
      if (!fs.treeMasks.empty() && fs.treeMasks.size() != forest.numTrees)
        return false;
      if (!fs.savedTreeMasks.empty() &&
          fs.savedTreeMasks.size() != fs.savedTrees.size())
        return false;
      if constexpr (L::hasVectorParams) {
        // the slope arrays must pair one-to-one with each flat tree's leaves
        // (a well-formed flat tree of m records has (m + 1) / 2 of them)
        size_t numSlopes = forest.leaf.numParams() - 1;
        if (fs.treeParams.size() != forest.numTrees) return false;
        for (size_t t = 0; t < forest.numTrees; ++t)
          if (fs.treeParams[t].size() !=
              ((fs.trees[t].size() + 1) / 2) * numSlopes)
            return false;
        if (!fs.savedTrees.empty()) {
          if (fs.savedTreeParams.size() != fs.savedTrees.size()) return false;
          for (size_t s = 0; s < fs.savedTrees.size(); ++s)
            if (fs.savedTreeParams[s].size() !=
                ((fs.savedTrees[s].size() + 1) / 2) * numSlopes)
              return false;
        }
      }
      if constexpr (L::hasFunctionParams) {
        // one fits slab per live tree; saved side channels must walk cleanly
        // against their trees' leaf counts
        if (fs.treeParams.size() != forest.numTrees) return false;
        for (size_t t = 0; t < forest.numTrees; ++t)
          if (fs.treeParams[t].size() != n) return false;
        if (!fs.savedTrees.empty()) {
          if (fs.savedTreeParams.size() != fs.savedTrees.size()) return false;
          std::vector<size_t> blockOffsets;
          for (size_t s = 0; s < fs.savedTrees.size(); ++s)
            if (!computeFunctionBlockOffsets(
                  fs.savedTreeParams[s].data(), fs.savedTreeParams[s].size(),
                  (fs.savedTrees[s].size() + 1) / 2, forest.leaf.numCovariates(),
                  blockOffsets))
              return false;
        }
      }
      for (size_t s = 0; s < fs.savedTrees.size(); ++s) {
        const std::vector<FlatNode>& saved(fs.savedTrees[s]);
        const std::uint64_t* masks = fs.savedTreeMasks.empty()
          ? nullptr : fs.savedTreeMasks[s].data();
        size_t numMaskWords =
          fs.savedTreeMasks.empty() ? 0 : fs.savedTreeMasks[s].size();
        if (!flatTreeIsWellFormed(data_, saved.data(), saved.size(), masks,
                                  numMaskWords))
          return false;
      }
      for (size_t t = 0; t < forest.numTrees; ++t) {
        scratch.initialize(scratchIndices.data(), n);
        // carry this forest's interaction constraint and column mask (both null
        // when unrestricted, so the walks short-circuit and the default path is
        // unchanged)
        scratch.setInteractionConstraint(forest.interaction.get());
        // the effective per-tree mask: a block row (already intersected with any
        // base mask at install) when blocks() is active, else the shared columnMask
        scratch.setColumnMask(
          !forest.blockMasks.empty()
            ? forest.blockMasks.data() + forest.blockOfTree[t] * data_.numPredictors
            : (forest.columnMask.empty() ? nullptr : forest.columnMask.data()));
        const std::uint64_t* masks =
          fs.treeMasks.empty() ? nullptr : fs.treeMasks[t].data();
        size_t numMaskWords =
          fs.treeMasks.empty() ? 0 : fs.treeMasks[t].size();
        if (!scratch.buildFromFlat(data_, fs.trees[t].data(),
                                   fs.trees[t].size(), params, 1, nullptr,
                                   masks, numMaskWords))
          return false;
        scratch.repartitionSubtree(data_, 0);
        if (!scratch.bottomNodesAreOccupied()) return false;
        // a state install must not admit a tree that violates the constraint
        // (design "Containment"): the availability predicate is not self-
        // checking, so treeLogProbability would mis-score a donor grown
        // unconstrained. Trivially passes for an unconstrained forest.
        if (!scratch.interactionSubtreeIsValid(0)) return false;
        // the same containment reasoning for the split-variable restriction (BCF
        // moderators, a column-restricted variance forest): a donor splitting on
        // a forbidden column would be mis-scored against an availability menu
        // that excludes it. Trivially passes for an unrestricted forest.
        if (!scratch.columnMaskSubtreeIsValid(0)) return false;
      }
    }
    if (!state.latents.empty() &&
        (response_->latents() == nullptr || state.latents.size() != n))
      return false;
    // a t sampler needs both its mixing precisions (lambda, in latents) and a
    // positive residual df; an old gaussian state, or one from a gaussian
    // sampler, carries neither and cannot continue the mixture
    if (response_->carriesResidualDf() &&
        (state.latents.size() != n || !(state.residualDf > 0.0)))
      return false;
    // an ordinal sampler needs its full length-(K-1) cutpoint vector; an old
    // state, or one from another family, carries none and cannot continue
    if (response_->carriesCutpoints() &&
        state.cutpoints.size() != response_->numCutpoints())
      return false;
    // an NB sampler needs both its omega latents (in latents) and a finite
    // positive dispersion r; an old state, or one from another family, carries
    // neither and cannot continue the augmentation (docs/design/negative-binomial.md)
    if (response_->carriesDispersion() &&
        (state.latents.size() != n || !(state.dispersion > 0.0) ||
         !std::isfinite(state.dispersion)))
      return false;
    // grouped states must carry a full effects vector for the chain's
    // groups; ungrouped states and chains both hold zero of them
    if (state.groupEffects.size() != response_->numGroupEffects())
      return false;
    if (forests_[0].useDart && !state.dartProbabilities.empty() &&
        state.dartProbabilities.size() != data_.numPredictors)
      return false;
    if (state.fitMax < state.fitMin) return false;
    // heteroscedastic: a variance state must carry one flat tree per variance
    // tree, each well-formed AND with every leaf a strictly positive scale (a
    // variance, unlike a Gaussian mean leaf) - the scale-leaf validation - AND
    // with every bottom occupied against this sampler's data, the criterion the
    // mean branch above imposes tree by tree. Without the occupancy pass an
    // installed variance tree could report a scale no row supports, which is
    // exactly what the transactional veto refuses to create
    // (docs/plans/multiforest-predictor-mutation.md, S3 item 5).
    if (varianceForest_) {
      if (state.varianceTrees.size() != varianceForest_->numTrees) return false;
      const std::uint8_t* varianceMask = varianceForest_->columnMask.empty()
        ? nullptr : varianceForest_->columnMask.data();
      for (const std::vector<FlatNode>& tree : state.varianceTrees) {
        if (!flatTreeIsWellFormed(data_, tree.data(), tree.size(), nullptr, 0))
          return false;
        for (const FlatNode& node : tree)
          if (flatKindOf(node) == FlatKind::leaf && !(node.value > 0.0))
            return false;
        scratch.initialize(scratchIndices.data(), n);
        // the variance forest carries no interaction constraint; clear the
        // mean loop's, which the shared scratch would otherwise still hold
        scratch.setInteractionConstraint(nullptr);
        scratch.setColumnMask(varianceMask);
        if (!scratch.buildFromFlat(data_, tree.data(), tree.size(), params))
          return false;
        scratch.repartitionSubtree(data_, 0);
        if (!scratch.bottomNodesAreOccupied()) return false;
      }
    } else if (!state.varianceTrees.empty()) {
      return false;
    }
    return true;
  }

  /// Whether every live tree in `state` satisfies its forest's interaction
  /// constraint - the warm-start containment gate (design "Containment"): a
  /// donor grown under a different (or no) constraint may hold a tree this
  /// sampler's constraint forbids, which treeLogProbability would mis-score.
  /// installForests calls this before touching live state so it can report a
  /// clear refusal; setState reaches the same guarantee through stateIsValid.
  /// Trivially true for an unconstrained forest, so the default warm start is
  /// byte-for-byte unchanged. Mirrors stateIsValid's structural scratch build
  /// (paramStride 1: only the rule structure, not leaf params, is examined).
  bool interactionStateFeasible(const ChainStateData& state) const {
    if (state.forests.size() != forests_.size()) return true;  // shape gate elsewhere
    size_t n = data_.numObservations;
    Tree scratch;
    std::vector<index_t> scratchIndices(n);
    std::vector<double> params;
    for (size_t f = 0; f < forests_.size(); ++f) {
      const Forest<L, ResidT>& forest = forests_[f];
      if (!forest.interaction || !forest.interaction->active()) continue;  // unconstrained: nothing to check
      const ForestStateData& fs = state.forests[f];
      if (fs.trees.size() != forest.numTrees) return true;  // shape gate elsewhere
      for (size_t t = 0; t < forest.numTrees; ++t) {
        scratch.initialize(scratchIndices.data(), n);
        scratch.setInteractionConstraint(forest.interaction.get());
        const std::uint64_t* masks =
          fs.treeMasks.empty() ? nullptr : fs.treeMasks[t].data();
        size_t numMaskWords =
          fs.treeMasks.empty() ? 0 : fs.treeMasks[t].size();
        // a malformed tree is the caller's shape concern (installForest fails
        // it); here we only judge the interaction feasibility of a buildable one
        if (!scratch.buildFromFlat(data_, fs.trees[t].data(), fs.trees[t].size(),
                                   params, 1, nullptr, masks, numMaskWords))
          continue;
        if (!scratch.interactionSubtreeIsValid(0)) return false;
      }
    }
    return true;
  }

  /// Whether every live tree in `state` respects its forest's column mask - the
  /// warm-start containment gate for the split-variable restriction (BCF
  /// moderators, a column-restricted variance forest): a donor grown under a
  /// different (or no) restriction may hold a tree that splits on a column this
  /// sampler's forest forbids, which splitVariableLogProbability would mis-score
  /// against an availability menu (collectAvailableVariables) that excludes it.
  /// installForests calls this before touching live state so it can report a
  /// clear refusal; setState reaches the same guarantee through stateIsValid.
  /// Trivially true for an unrestricted forest, so the default warm start is
  /// byte-for-byte unchanged. Mirrors interactionStateFeasible's scratch build
  /// (paramStride 1: only the rule structure, not leaf params, is examined).
  bool columnMaskStateFeasible(const ChainStateData& state) const {
    if (state.forests.size() != forests_.size()) return true;  // shape gate elsewhere
    size_t n = data_.numObservations;
    Tree scratch;
    std::vector<index_t> scratchIndices(n);
    std::vector<double> params;
    for (size_t f = 0; f < forests_.size(); ++f) {
      const Forest<L, ResidT>& forest = forests_[f];
      // a blocks() forest carries its per-tree restriction in blockMasks, not the
      // shared columnMask; either (or an intersected tau block row) must be gated
      if (forest.columnMask.empty() && forest.blockMasks.empty()) continue;  // unrestricted
      const ForestStateData& fs = state.forests[f];
      if (fs.trees.size() != forest.numTrees) return true;  // shape gate elsewhere
      for (size_t t = 0; t < forest.numTrees; ++t) {
        scratch.initialize(scratchIndices.data(), n);
        // the effective per-tree mask: a block row (already intersected with any
        // base moderator/variance mask at install) when blocks() is active, else
        // the shared forest columnMask
        scratch.setColumnMask(forest.blockMasks.empty()
          ? forest.columnMask.data()
          : forest.blockMasks.data() + forest.blockOfTree[t] * data_.numPredictors);
        const std::uint64_t* masks =
          fs.treeMasks.empty() ? nullptr : fs.treeMasks[t].data();
        size_t numMaskWords =
          fs.treeMasks.empty() ? 0 : fs.treeMasks[t].size();
        // a malformed tree is the caller's shape concern (installForest fails
        // it); here we only judge the column feasibility of a buildable one
        if (!scratch.buildFromFlat(data_, fs.trees[t].data(), fs.trees[t].size(),
                                   params, 1, nullptr, masks, numMaskWords))
          continue;
        if (!scratch.columnMaskSubtreeIsValid(0)) return false;
      }
    }
    // the variance forest sits outside forests_, so a `variance = ~ subset`
    // restriction needs its own pass over the state's variance trees - the
    // same scratch build, against varianceForest_->columnMask. A state
    // carrying no variance trees skips, as does an unrestricted variance
    // forest (no mask).
    if (varianceForest_ && !varianceForest_->columnMask.empty() &&
        state.varianceTrees.size() == varianceForest_->numTrees) {
      for (const std::vector<FlatNode>& flat : state.varianceTrees) {
        scratch.initialize(scratchIndices.data(), n);
        scratch.setColumnMask(varianceForest_->columnMask.data());
        if (!scratch.buildFromFlat(data_, flat.data(), flat.size(), params))
          continue;
        if (!scratch.columnMaskSubtreeIsValid(0)) return false;
      }
    }
    return true;
  }

  /// Rebuilds forest f's live trees, partitions, and fits from a flat state's
  /// live channel against the current cut grid, zeroing and re-accumulating
  /// totalFits. False if a flat tree fails to rebuild. Shared by setState and
  /// the warm-start installForest.
  bool rebuildLiveForest(size_t f, const ForestStateData& fs,
                         std::vector<double>& params) {
    Forest<L, ResidT>& forest = forests_[f];
    size_t n = data_.numObservations;
    misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);
    for (size_t t = 0; t < forest.numTrees; ++t) {
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
      const std::uint64_t* masks =
        fs.treeMasks.empty() ? nullptr : fs.treeMasks[t].data();
      size_t numMaskWords =
        fs.treeMasks.empty() ? 0 : fs.treeMasks[t].size();
      if constexpr (!L::hasVectorParams) {
        if (!forest.trees[t].buildFromFlat(data_, fs.trees[t].data(),
                                           fs.trees[t].size(), params, 1,
                                           nullptr, masks, numMaskWords))
          return false;
      } else {
        if (!forest.trees[t].buildFromFlat(data_, fs.trees[t].data(),
                                           fs.trees[t].size(), params,
                                           forest.leaf.numParams(),
                                           fs.treeParams[t].data(), masks,
                                           numMaskWords))
          return false;
      }
      forest.trees[t].repartitionSubtree(data_, 0);
      // containment backstop (design): the live tree carries this forest's
      // constraint and column mask, so a warm-start donor grown unconstrained is
      // caught before treeLogProbability can mis-score it. installForests
      // pre-checks for the clear error; this guarantees the invariant on any
      // live-install path. Trivially passes for an unrestricted forest (null
      // short-circuit).
      if (!forest.trees[t].interactionSubtreeIsValid(0)) return false;
      if (!forest.trees[t].columnMaskSubtreeIsValid(0)) return false;
      if constexpr (L::hasFunctionParams) {
        // the recorded slab IS the tree's parameters; copy restores bitwise
        std::memcpy(forest.treeFits.data() + t * n, fs.treeParams[t].data(),
                    n * sizeof(double));
      } else {
        setTreeFits(forest, t, params);
      }
      if constexpr (leafIsConstant) rebuildLeafOf(forest, t);
      addTreeFitsToTotal(forest, t);
    }
    return true;
  }

  /// Cross-grid warm start's per-forest rebuild: the donor's flat trees were
  /// grown on donorCutPoints, a grid this sampler no longer carries, so build
  /// each tree's structure against that grid - where its ordinal split values
  /// resolve exactly - then remap the recovered split indices onto the live
  /// grid, collapsing splits the coarser or shifted grid starves, exactly as
  /// setData's applyNewData does. The donor grid is installed over the shared
  /// store only for the structural buildFromFlat (ScopedCutGrid; the live
  /// observation codes are never touched) and reverts before the remap and
  /// repartition. Standardization constants and, for vector or function leaves,
  /// per-observation leaf state re-anchor to the live data the same way
  /// applyNewData's do (node parameters carry through the remap; function fits,
  /// being per-donor-observation, cold-start). False if a flat tree fails to
  /// rebuild or the remapped tree violates this forest's constraint. store is
  /// this chain's own shared data (the const data_ aliases it), handed in
  /// mutably because only the sampler owns it: the donor grid is swapped over it
  /// for the structural build alone.
  bool rebuildLiveForestRemapped(
      size_t f, const ForestStateData& fs,
      const std::vector<std::vector<double>>& donorCutPoints,
      std::vector<double>& params, ColumnStore& store) {
    Forest<L, ResidT>& forest = forests_[f];
    size_t n = data_.numObservations;
    misc_setVectorToConstant(forest.totalFits.data(), n, 0.0);

    // fresh standardization / per-observation leaf state over the live data,
    // the same re-anchoring applyNewData performs after a data replacement
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      forest.leaf.reinitialize(data_);

    size_t paramStride = 1;
    if constexpr (L::hasVectorParams) paramStride = forest.leaf.numParams();

    for (size_t t = 0; t < forest.numTrees; ++t) {
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
      const std::uint64_t* masks =
        fs.treeMasks.empty() ? nullptr : fs.treeMasks[t].data();
      size_t numMaskWords =
        fs.treeMasks.empty() ? 0 : fs.treeMasks[t].size();

      // recover the donor split indices by building against the donor grid; the
      // swap is structural (codes untouched) and reverts before the remap even
      // on the failure return (ScopedCutGrid's destructor)
      {
        ScopedCutGrid donorGrid(store, donorCutPoints);
        if constexpr (!L::hasVectorParams) {
          if (!forest.trees[t].buildFromFlat(store, fs.trees[t].data(),
                                             fs.trees[t].size(), params, 1,
                                             nullptr, masks, numMaskWords))
            return false;
        } else {
          if (!forest.trees[t].buildFromFlat(store, fs.trees[t].data(),
                                             fs.trees[t].size(), params,
                                             forest.leaf.numParams(),
                                             fs.treeParams[t].data(), masks,
                                             numMaskWords))
            return false;
        }
      }

      // function-valued fits are per-donor-observation; cold-start them, as
      // applyNewData does (the next sweep's draws replace them)
      if constexpr (L::hasFunctionParams)
        params.assign(forest.trees[t].nodes.size(), 0.0);

      // remap onto the live grid, then re-route and collapse starved or emptied
      // subtrees - the applyNewData body, without the response/store rebuild
      forest.trees[t].mapOldCutPointsOntoNew(data_, donorCutPoints, params,
                                             paramStride);
      forest.trees[t].repartitionSubtree(data_, 0);
      forest.trees[t].collapseEmptyNodes(data_, response_->workingWeights(),
                                         params, paramStride);
      // containment backstop (as rebuildLiveForest): a donor grown unconstrained
      // cannot leave an infeasible live tree, even after the remap collapses
      if (!forest.trees[t].interactionSubtreeIsValid(0)) return false;
      if (!forest.trees[t].columnMaskSubtreeIsValid(0)) return false;
      if constexpr (!L::hasVectorParams) {
        setTreeFitsFromParameters(forest, t, params);
      } else {
        forest.paramsByTree[t] = params;
        setTreeFitsFromParameterBlocks(forest, t, params);
      }
      if constexpr (leafIsConstant) rebuildLeafOf(forest, t);
      addTreeFitsToTotal(forest, t);
    }
    return true;
  }

  /// Warm start: seed the live forest(s), sigma, and k from a donor's flat
  /// trees, leaving this chain's rng, latents, group effects, and saved-tree
  /// buffer untouched - the donor supplies a starting position, not a
  /// continuation. Callers guarantee shape compatibility; false signals only a
  /// flat tree that failed to rebuild. donorCutPoints null installs the donor's
  /// splits verbatim (the donor shares this sampler's grid); non-null remaps
  /// them onto the live grid, collapsing starved splits (a cross-grid start),
  /// and requires store - this chain's shared data, mutable so the donor grid
  /// can be swapped in for the structural build.
  bool installForest(const ChainStateData& state,
                     const std::vector<std::vector<double>>* donorCutPoints =
                       nullptr,
                     ColumnStore* store = nullptr) {
    if (state.forests.size() != forests_.size()) return false;
    if (state.fitMax > state.fitMin)
      response_->restoreScale(state.fitMin, state.fitMax);
    std::vector<double> params;
    for (size_t f = 0; f < forests_.size(); ++f) {
      const ForestStateData& fs = state.forests[f];
      if (fs.trees.size() != forests_[f].numTrees) return false;
      bool rebuilt = donorCutPoints == nullptr
        ? rebuildLiveForest(f, fs, params)
        : rebuildLiveForestRemapped(f, fs, *donorCutPoints, params, *store);
      if (!rebuilt) return false;
      forests_[f].k = fs.k;
      // the leaf prior's other half, adopted like k, sigma, the transform, DART
      // and the glue already are: a donor's trees were drawn under its scale, so
      // installing the trees without it leaves a hybrid (donor units,
      // destination calibration). An absent block (0.0, or any non-positive or
      // non-finite value - k's posture, no new refusal) leaves construction's.
      if (fs.leafScale > 0.0) forests_[f].leaf.scale = fs.leafScale;
    }
    setSigma(state.sigma);
    Forest<L, ResidT>& forest = forests_[0];
    if (forest.useDart && !state.dartProbabilities.empty()) {
      std::memcpy(forest.dart.probabilities.data(),
                  state.dartProbabilities.data(),
                  state.dartProbabilities.size() * sizeof(double));
      forest.dart.alpha = state.dartAlpha;
      forest.dart.setNumUpdatesSkipped(state.dartNumUpdatesSkipped);
    }
    if (combiner_) combiner_->restoreGlue(state);
    return true;
  }

  /// The warm start's variance half, run by installForests right after
  /// installForest so a refusal can name the variance forest rather than the
  /// shape. Same grid: the donor's flat trees resolve against the live grid and
  /// rebuildVarianceForest installs them verbatim, as setState does. Cross
  /// grid: build each tree's structure against the donor grid - where its
  /// ordinal split values resolve exactly - then hand it to
  /// refreshVarianceForest's remap arm, mirroring rebuildLiveForestRemapped
  /// (the recovered factors come from the flat trees, not from the live slab,
  /// which still holds the destination's own surface). False when a flat tree
  /// fails to rebuild, or when a rebuilt tree leaves a bottom unoccupied: an
  /// empty bottom carries no drawn factor, and while recoverVarianceLeafValues
  /// keeps that safe by abstaining (1.0), installing one would report a scale
  /// this data never supported. Only the same-grid arm can produce one - the
  /// cross-grid remap collapses empty nodes.
  bool installVarianceForest(
      const std::vector<std::vector<FlatNode>>& trees,
      const std::vector<std::vector<double>>* donorCutPoints,
      ColumnStore* store) {
    VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    if (trees.size() != vf.numTrees) return false;
    if (donorCutPoints == nullptr) {
      if (!rebuildVarianceForest(trees)) return false;
    } else {
      TreeParameters recovered(vf.numTrees);
      {
        ScopedCutGrid donorGrid(*store, *donorCutPoints);
        for (std::size_t j = 0; j < vf.numTrees; ++j) {
          Tree& tree = vf.trees[j];
          tree.initialize(vf.indexBuffer.data() + j * n, n);
          if (!tree.buildFromFlat(*store, trees[j].data(), trees[j].size(),
                                  recovered[j]))
            return false;
        }
      }
      refreshVarianceForest(donorCutPoints, &recovered);
    }
    for (std::size_t j = 0; j < vf.numTrees; ++j)
      if (!vf.trees[j].bottomNodesAreOccupied()) return false;
    return true;
  }

  /// Installs a state stateIsValid accepted; false only on the invariant
  /// violation of a validated tree failing to rebuild.
  bool setState(const ChainStateData& state) {
    if (state.forests.size() != forests_.size()) return false;
    // the internal-scale tree parameters and fits below were recorded under
    // this transform; scale-free states leave creation's. restoreScale
    // re-anchors the variance prior through it.
    if (state.fitMax > state.fitMin)
      response_->restoreScale(state.fitMin, state.fitMax);
    std::vector<double> params;
    for (size_t f = 0; f < forests_.size(); ++f) {
      Forest<L, ResidT>& forest = forests_[f];
      const ForestStateData& fs = state.forests[f];
      if (!rebuildLiveForest(f, fs, params)) return false;
      if (!fs.savedTrees.empty()) {
        forest.savedTrees = fs.savedTrees;
        if constexpr (L::hasVectorParams || L::hasFunctionParams)
          forest.savedTreeParams = fs.savedTreeParams;
        if (data_.hasPooledCategorical) {
          if (fs.savedTreeMasks.empty())
            forest.savedTreeMasks.assign(forest.savedTrees.size(),
                                         std::vector<std::uint64_t>());
          else
            forest.savedTreeMasks = fs.savedTreeMasks;
        }
      }
      forest.k = fs.k;
      // as in installForest above. CONSEQUENCE: a setModel(node.scale) issued
      // AFTER the last storeState no longer survives a save/load re-creation,
      // since the state's scale now wins - exactly the wart k already had,
      // applied consistently to both halves of the leaf prior.
      if (fs.leafScale > 0.0) forest.leaf.scale = fs.leafScale;
    }
    setSigma(state.sigma);
    // RESTORE CONTRACT (docs/design/negative-binomial.md section 5): an NB
    // response's restoreLatents rebuilds the working response from omega AND r,
    // so the dispersion MUST be reinstalled before the latents - a restore that
    // installs omega first would rebuild working against the stale r. stateIsValid
    // guaranteed a finite positive r for an NB sampler.
    if (response_->carriesDispersion())
      response_->restoreDispersion(state.dispersion);
    if (!state.latents.empty())
      response_->restoreLatents(state.latents.data());
    // stateIsValid guaranteed a positive df for a t sampler; fixed mode
    // reinstalls its constant, estimated mode its last grid draw
    if (response_->carriesResidualDf())
      response_->restoreResidualDf(state.residualDf);
    // stateIsValid guaranteed a full length-(K-1) cutpoint vector for an
    // ordinal sampler; z was restored above under these same cutpoints
    if (response_->carriesCutpoints())
      response_->restoreCutpoints(state.cutpoints.data());
    if (!state.groupEffects.empty())
      response_->restoreGroupEffects(state.groupEffects.data(),
                                     state.groupTau);
    Forest<L, ResidT>& forest = forests_[0];
    if (forest.useDart && !state.dartProbabilities.empty()) {
      // the tree prior points at this vector's storage; overwrite in place
      std::memcpy(forest.dart.probabilities.data(),
                  state.dartProbabilities.data(),
                  state.dartProbabilities.size() * sizeof(double));
      forest.dart.alpha = state.dartAlpha;
      forest.dart.setNumUpdatesSkipped(state.dartNumUpdatesSkipped);
    }
    if (combiner_) combiner_->restoreGlue(state);
    // heteroscedastic: rebuild the variance trees and recompute s^2(x) from the
    // restored positive factors (stateIsValid checked count, form, positivity)
    if (varianceForest_ && !rebuildVarianceForest(state.varianceTrees))
      return false;
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
  /// The dense per-tree fit slab, materialized (the constant leaf gathers its
  /// compact tables into the returned buffer); a consistency read for tests.
  std::vector<double> treeFits() const {
    std::vector<double> out(data_.numObservations * forests_[0].numTrees);
    forestTreeFits(0, out.data());
    return out;
  }
  const std::vector<double>& totalFits() const { return forests_[0].totalFits; }
  /// The per-forest sibling of totalFits(), which addresses forest 0 alone;
  /// a whole-sampler consistency read (the fuzz snapshot) needs every forest's.
  const std::vector<double>& totalFitsInForest(std::size_t f) const {
    return forests_[f].totalFits;
  }
  /// Test hook: tree t's obs-to-leaf map (constant leaf, forest 0), where entry
  /// i is the arena bottom-node index owning observation i.
  const std::uint32_t* leafOfForTesting(size_t t) const {
    return forests_[0].leafOf.data() + t * data_.numObservations;
  }
  /// Test hook: tree t's rebuild mark, the fused suffstat pass's eligibility
  /// gate (forest 0). Non-zero means the map still describes the previous
  /// partition, so the sweep rebuilds before the tree's draw.
  std::uint8_t leafOfStaleForTesting(size_t t) const {
    return forests_[0].leafOfStale[t];
  }
  /// Test hooks: the running residual forest 0 rolls across a sweep, and the
  /// working response it is rolled against. Together with treeFits() they pin
  /// the unrolled mu[leafOf] gathers elementwise, tail included.
  const std::vector<ResidT>& residualForTesting() const {
    return forests_[0].treeY;
  }
  const double* workingResponseForTesting() const {
    return response_->workingResponse();
  }
  /// Test hook: the sigma posterior's degrees of freedom, nu_0 plus the count
  /// of positive precisions on the RESPONSE model. A per-forest weight lives on
  /// the chain and must never reach them.
  double sigmaDegreesOfFreedomForTesting() const {
    return response_->sigmaDegreesOfFreedomForTesting();
  }
  /// Forest f's tree t. Chain::tree reaches forest 0 alone and keeps that
  /// meaning for its own callers; this is what a whole-sampler walk resolves
  /// through - the per-observation session's survivor table, and the
  /// structural checks (routing agreement) that must see a BCF or multinomial
  /// forest at all.
  const Tree& treeInForest(std::size_t f, std::size_t t) const {
    return forests_[f].trees[t];
  }
  /// Variance tree j. The scale twin of treeInForest: the variance forest sits
  /// outside forests_ and has no numTreesInForest entry, so a whole-sampler
  /// walk - the per-observation session's survivor table, the structural
  /// checks - reaches its trees through this and counts them with
  /// numVarianceTrees(). Dereferences the forest; a caller gates on
  /// hasVarianceForest().
  const Tree& varianceTree(std::size_t j) const {
    return varianceForest_->trees[j];
  }
  /// Test hook: the per-tree factor slab h_j(x_i), tree-major
  /// (numVarianceTrees x n), whose product over j is the combined variance
  /// varianceFits() reports. Nothing else exposes it to a test.
  const double* varianceFactorsForTesting() const {
    return varianceForest_->factorByTree.data();
  }

  /// Test hook: split forest 0's tree 0 at (variableIndex, splitIndex) and
  /// strand its right child empty, then run tree 0's parameter draw exactly
  /// as a sweep does (updateK forced on) and report the k accounting the
  /// chi-k draw consumes. The stranded empty leaf must contribute nothing to
  /// the leaf count or the sum of squares, matching the function-leaf path.
  /// No public mutation strands an empty leaf, so this fabricates one.
  FunctionLeafDrawStats accountStrandedLeafKStatsForTesting(int32_t variableIndex,
                                                            int32_t splitIndex) {
    Forest<L, ResidT>& forest = forests_[0];
    Tree& tree = forest.trees[0];
    const double* residual = response_->workingResponse();
    const double* weights = response_->workingWeights();
    size_t n = data_.numObservations;

    Rule rule;
    rule.variableIndex = variableIndex;
    rule.setSplitIndex(splitIndex);
    tree.birth(data_, 0, rule, residual, weights);
    int32_t rightChild = tree.at(0).leftChild + 1;
    tree.at(rightChild).begin = tree.at(rightChild).end;
    if constexpr (leafIsConstant) rebuildLeafOf(forest, 0);

    forest.updateK = true;
    forest.kSumSquaredParams = 0.0;
    forest.kNumLeaves = 0.0;
    if constexpr (leafTracksNodeAverages)
      tree.setNodeAverages(residual, weights);
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      std::copy(residual, residual + n, forest.treeY.begin());

    std::vector<double> fits(n, 0.0);
    sampleParametersAndSetFits(forest, 0, fits.data(), false);
    return FunctionLeafDrawStats{forest.kSumSquaredParams, forest.kNumLeaves};
  }

  /// Diagnostic: (tree, sweep) bodies that took the fused roll + suffstat
  /// pass since this chain was built. Monotone, so tests read differences.
  size_t fusedSuffstatRunsForTesting() const { return fusedSuffstatRuns_; }

  /// Test hook: for each of forest 0's trees, run the stock
  /// rollTreeResidual + setNodeAverages pair and the fused pass over the SAME
  /// entering residual, and report whether they agree. The two are specified
  /// to differ in exactly one way - the suffstat's summation association - and
  /// nothing else pins that; in particular the bitwise residual identity is
  /// what keeps this landing out of the shifting class on treeY.
  ///
  /// Drives forest 0 against the response model's own working response and
  /// weights, so it is meaningful only for a single-forest chain with no
  /// variance forest (elsewhere run() supplies combiner or mean weights).
  /// Draws no rng and makes no move, so the chain's draw sequence is
  /// untouched; it does advance the fused run counter, and it leaves treeY
  /// and the node statistics where a roll-only sweep would - both are
  /// recomputed at the top of every tree body.
  FusedSuffstatCheck checkFusedSuffstatAgainstStockForTesting() {
    FusedSuffstatCheck result;
    if constexpr (!leafIsConstant) {
      return result;
    } else {
      Forest<L, ResidT>& forest = forests_[0];
      const double* forestY = response_->workingResponse();
      const double* forestWeights = response_->workingWeights();
      size_t n = data_.numObservations;
      std::vector<ResidT> entering, stockResid;
      std::vector<double> stockSumWeights, stockSumResponse;

      for (size_t t = 0; t < forest.numTrees; ++t) {
        Tree& tree(forest.trees[t]);
        entering.assign(forest.treeY.begin(), forest.treeY.end());

        rollTreeResidual(forest, t, forestY);
        tree.setNodeAverages(forest.treeY.data(), forestWeights);
        stockResid.assign(forest.treeY.begin(), forest.treeY.end());
        stockSumWeights.clear();
        stockSumResponse.clear();
        tree.bottomScratch.clear();
        tree.fillBottom(0, tree.bottomScratch);
        for (int32_t b : tree.bottomScratch) {
          stockSumWeights.push_back(tree.at(b).sumWeights);
          stockSumResponse.push_back(tree.at(b).sumWeightedResponse);
        }

        std::copy(entering.begin(), entering.end(), forest.treeY.begin());
        if (!rollAndSetNodeAveragesFused(forest, t, forestY, forestWeights)) {
          std::copy(stockResid.begin(), stockResid.end(),
                    forest.treeY.begin());
          continue;
        }
        ++result.numTreesFused;
        for (size_t i = 0; i < n; ++i)
          if (forest.treeY[i] != stockResid[i])
            result.residualsAgreeBitwise = false;
        size_t b = 0;
        for (int32_t nodeIndex : tree.bottomScratch) {
          const Node& node(tree.at(nodeIndex));
          if (node.sumWeights != stockSumWeights[b])
            result.countsAgreeBitwise = false;
          double reference = stockSumResponse[b];
          double scale = std::max(std::fabs(reference), 1.0);
          result.worstRelativeError =
            std::max(result.worstRelativeError,
                     std::fabs(node.sumWeightedResponse - reference) / scale);
          ++b;
        }
      }
      return result;
    }
  }

private:
  /// Composes forest f's installed per-observation weight into the precisions
  /// the combiner just formed, returning a chain-owned scratch pointer; the
  /// combiner's own pointer passes through with no allocation, no copy and no
  /// arithmetic when forest f carries no weight, which is the whole cost of the
  /// channel on every configuration that installs none. Called the moment
  /// formForestResponse returns, before the tree loop reads the precisions -
  /// grow-from-root consumes them immediately, in setNodeAverages. The
  /// combiner's weights are never null, so neither is the composed vector.
  const double* composeForestWeights(std::size_t f, const double* weights) {
    if (forestWeights_.empty() || forestWeights_[f] == nullptr) return weights;
    std::size_t n = data_.numObservations;
    forestWeightScratch_.resize(n);
    const double* s = forestWeights_[f];
    for (std::size_t i = 0; i < n; ++i)
      forestWeightScratch_[i] = weights[i] * s[i];
    return forestWeightScratch_.data();
  }

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

  /// Recursive growth from the prior: growth is Bernoulli in the
  /// depth-decayed prior probability, rules come from the prior, and empty
  /// children keep growing (availability is rule-based) until the caller
  /// collapses them.
  void growSubtreeFromPrior(Forest<L, ResidT>& forest, Tree& tree, int32_t nodeIndex,
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

  /// Constant-leaf fit storage sizing: node-indexed mu tables (one zero-value
  /// root per tree) and the tree-major obs-to-leaf map (all-root is all zeros).
  /// Vector and function leaves size the dense slab instead.
  void initForestFitStorage(Forest<L, ResidT>& forest, size_t n) {
    if constexpr (leafIsConstant) {
      forest.muByTree.assign(forest.numTrees, std::vector<double>(1, 0.0));
      forest.leafOf.assign(n * forest.numTrees, 0);
      forest.leafOfStale.assign(forest.numTrees, 0);
    } else {
      forest.treeFits.assign(n * forest.numTrees, 0.0);
    }
  }

  /// Rewrite tree t's obs-to-leaf entries for the members of nodeIndex's
  /// subtree (constant leaf): every member of a bottom node points at that
  /// node. Entries outside the subtree's segment are untouched, so an accepted
  /// move's patch costs its repartitioned members, not n.
  void updateLeafOfBelow(Forest<L, ResidT>& forest, size_t t, int32_t nodeIndex) {
    Tree& tree(forest.trees[t]);
    std::uint32_t* leaf = forest.leafOf.data() + t * data_.numObservations;
    tree.bottomScratch.clear();
    tree.fillBottom(nodeIndex, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      for (size_t m = node.begin; m < node.end; ++m)
        leaf[tree.indices[m]] = static_cast<std::uint32_t>(i);
    }
  }

  /// Full obs-to-leaf rebuild for a wholesale partition change; clears the
  /// tree's staleness mark.
  void rebuildLeafOf(Forest<L, ResidT>& forest, size_t t) {
    updateLeafOfBelow(forest, t, 0);
    forest.leafOfStale[t] = 0;
  }

  /// rebuildLeafOf + addTreeFitsToTotal in one walk over the new bottoms
  /// (each observation receives one map write and one add either way, so the
  /// fusion is value-identical); the data-mutation transaction runs this per
  /// tree, where the separate passes measurably cost.
  void installLeafOfAndAddToTotal(Forest<L, ResidT>& forest, size_t t) {
    Tree& tree(forest.trees[t]);
    const std::vector<double>& mu(forest.muByTree[t]);
    std::uint32_t* leaf = forest.leafOf.data() + t * data_.numObservations;
    double* total = forest.totalFits.data();
    tree.bottomScratch.clear();
    tree.fillBottom(0, tree.bottomScratch);
    for (int32_t i : tree.bottomScratch) {
      const Node& node(tree.at(i));
      double value = mu[static_cast<size_t>(i)];
      for (size_t m = node.begin; m < node.end; ++m) {
        size_t obs = tree.indices[m];
        leaf[obs] = static_cast<std::uint32_t>(i);
        total[obs] += value;
      }
    }
    forest.leafOfStale[t] = 0;
  }

  /// Build the nullable variance forest (heteroscedastic): calibrate its scale
  /// leaf from THIS chain's homoscedastic sigma prior on the working scale
  /// (df = sigmaDf, scale = initialSigma^2 * rawScale, the GaussianResponse
  /// derivation), seed s^2(x) at the initial variance, then fix the global
  /// sigma at 1 - the variance forest carries the residual variance from here.
  /// At numVarianceTrees == 1 the calibration reproduces the sigma prior exactly.
  void buildVarianceForest(const SamplerOptions& options, double sigmaEstimate,
                           double sigmaDf, double sigmaRawScale) {
    varianceLeafPrior_ = {sigmaEstimate, sigmaDf, sigmaRawScale};
    std::size_t n = data_.numObservations;
    double initialVariance = sigma_ * sigma_;  // sigma_ still holds initialSigma
    double priorScale = initialVariance * sigmaRawScale;
    varianceForest_ = std::make_unique<VarianceForest>();
    VarianceForest& vf = *varianceForest_;
    vf.birthOrDeathProbability = options.birthOrDeathProbability;
    vf.swapProbability = options.swapProbability;
    vf.changeProbability = options.changeProbability;
    vf.birthProbability = options.birthProbability;
    vf.treePrior.base = options.varianceBase;
    vf.treePrior.power = options.variancePower;
    vf.leaf = ConstantVarianceLeaf::calibrated(sigmaDf, priorScale,
                                               options.numVarianceTrees);
    vf.initialize(options.numVarianceTrees, n, initialVariance);
    // restrict the variance trees to the `variance = ~ subset` columns; empty
    // leaves every column available, byte-for-byte the unrestricted path
    if (options.varianceForestColumns != nullptr &&
        options.numVarianceForestColumns > 0) {
      vf.columnMask.assign(data_.numPredictors, 0);
      for (std::size_t c = 0; c < options.numVarianceForestColumns; ++c)
        vf.columnMask[options.varianceForestColumns[c]] = 1;
      for (std::size_t t = 0; t < vf.numTrees; ++t)
        vf.trees[t].setColumnMask(vf.columnMask.data());
    }
    sigmaIsFixed_ = true;
    sigma_ = 1.0;
    meanWeights_.assign(n, 0.0);
  }

  /// Divide the user precisions by the current variance surface for the mean
  /// forest's next sweep: meanWeights_[i] = user_w_i / s^2(x_i) (unit weight
  /// where the user supplied none). One of the two heteroscedastic n-scratch
  /// vectors (the other is the variance forest's per-tree divisor s^2_{-j}).
  void formMeanWeights() {
    const VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    const double* userWeights = response_->workingWeights();
    for (std::size_t i = 0; i < n; ++i)
      meanWeights_[i] =
        (userWeights == nullptr ? 1.0 : userWeights[i]) / vf.combinedVariance[i];
  }

  /// One variance-forest sweep: backfit each variance tree against the mean
  /// residual scaled by the OTHER trees' product s^2_{-j} (the multiplicative
  /// roll), score/move it with the existing conjugate machinery over the scale
  /// leaf, redraw its leaf factors, and fold them into s^2(x). The scale leaf's
  /// suffstat carries the user case weights (sum_i w_i e_i^2 / s^2_{-j}).
  void sweepVarianceForest(const double* y, const double* meanFits) {
    VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    const double* userWeights = response_->workingWeights();
    for (std::size_t i = 0; i < n; ++i) vf.meanResidual[i] = y[i] - meanFits[i];

    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      vf.formTreeResidual(j, vf.meanResidual.data());
      MoveContext ctx{data_,
                      vf.treePrior,
                      vf.birthOrDeathProbability,
                      vf.swapProbability,
                      vf.birthProbability,
                      userWeights,
                      1.0,  // k: unread by the scale leaf's marginal
                      vf.scratch};
      bool stepTaken;
      StepType stepType;
      int32_t changedNode = invalidNode;
      metropolisJumpForTree(ctx, vf.leaf, rng_, vf.trees[j],
                            vf.treeResidual.data(), 1.0, &stepTaken, &stepType,
                            &changedNode);
      if (data_.hasPooledCategorical)
        vf.trees[j].compactMaskPoolIfNeeded(data_);

      std::vector<int32_t>& bottoms(vf.trees[j].bottomScratch);
      bottoms.clear();
      vf.trees[j].fillBottom(0, bottoms);
      for (int32_t b : bottoms) {
        double h = vf.leaf.drawFromPosteriorForNode(
          rng_, vf.trees[j], vf.treeResidual.data(), userWeights, 1.0, 1.0, b);
        vf.applyLeafFactor(j, b, h);
      }
    }
    // recompute s^2(x) as the fresh product (drift reset, so it equals the
    // reporting/predict recompute exactly); the running per-tree update above
    // stays correct within the sweep, this only retires its accumulated rounding
    for (std::size_t i = 0; i < n; ++i) {
      double product = 1.0;
      for (std::size_t j = 0; j < vf.numTrees; ++j)
        product *= vf.factorByTree[j * n + i];
      vf.combinedVariance[i] = product;
    }
  }

  /// Node-indexed (arena id) leaf factors for variance tree j, recovered from
  /// the per-observation slab (every member of a leaf carries its value): the
  /// scale analogue of recoverParametersFromFits, for flattening and predict.
  ///
  /// A bottom with no members carries no drawn factor and reads the
  /// multiplicative identity 1.0 - that tree abstains where it has no training
  /// support - so an unsupported test row reads the product over the trees that
  /// DO have support instead of zero, and a flattened state stays inside its
  /// own strict-positivity check. The empty-leaf veto keeps a live tree's
  /// bottoms occupied, so only an installed tree (setState today, a warm start
  /// once installForest carries variance trees) can need the fallback;
  /// refreshVarianceForest asserts the live invariant at its recover step.
  /// Internal slots are never read - flatten and the merges take bottoms only.
  void recoverVarianceLeafValues(const VarianceForest& vf, std::size_t j,
                                 std::vector<double>& out) const {
    const Tree& tree = vf.trees[j];
    const double* factor = vf.factorByTree.data() + j * data_.numObservations;
    out.assign(tree.nodes.size(), 1.0);
    recoverVarianceLeafValuesBelow(tree, 0, factor, out);
  }

  /// The recursion, over the LIVE tree rather than the node arena: a released
  /// pair stays in `nodes` on the free list, reads as a bottom, and keeps the
  /// index range it held when it was live - which a setData that SHRINKS n
  /// puts past the end of the resized index buffer. Walking from the root is
  /// also what the mean side's functionLeafValues does.
  void recoverVarianceLeafValuesBelow(const Tree& tree, std::int32_t nodeIndex,
                                      const double* factor,
                                      std::vector<double>& out) const {
    const Node& node = tree.at(nodeIndex);
    if (node.isBottom()) {
      if (node.numObservations() > 0)
        out[static_cast<std::size_t>(nodeIndex)] =
          factor[tree.indices[node.begin]];
      return;
    }
    recoverVarianceLeafValuesBelow(tree, node.leftChild, factor, out);
    recoverVarianceLeafValuesBelow(tree, node.leftChild + 1, factor, out);
  }

  /// Flatten every variance tree into saved slot `slot` for this kept sample,
  /// so predict can replay the posterior variance surface on new data.
  void storeVarianceSavedTrees(std::size_t sampleNum) {
    VarianceForest& vf = *varianceForest_;
    std::size_t slot =
      (vf.savedSlotBase + sampleNum) % vf.savedTreeCapacity;
    std::vector<double> leafValues;
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      recoverVarianceLeafValues(vf, j, leafValues);
      vf.trees[j].flatten(data_, leafValues.data(),
                          vf.savedTrees[slot * vf.numTrees + j]);
    }
  }

  /// Route the test rows through every variance tree and multiply their leaf
  /// factors into s^2(x_test), the test analogue of the maintained combined
  /// variance; called only when a recorded sweep needs the test channel.
  void refreshVarianceTestFits() {
    VarianceForest& vf = *varianceForest_;
    std::size_t nTest = data_.numTestObservations;
    for (std::size_t i = 0; i < nTest; ++i) vf.combinedVarianceTest[i] = 1.0;
    std::vector<double> leafValues;
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      recoverVarianceLeafValues(vf, j, leafValues);
      const Tree& tree = vf.trees[j];
      routeTestRows(nTest, [&](std::size_t i) {
        int32_t leaf = tree.findBottomNodeForRow(data_, i);
        vf.combinedVarianceTest[i] *=
          leafValues[static_cast<std::size_t>(leaf)];
      });
    }
  }

  /// Rebuild the variance forest from a flat state's variance trees: rebuild
  /// each tree, scatter its positive leaf factors to the per-observation slab
  /// through the restored partition, then recompute s^2(x) as the product.
  bool rebuildVarianceForest(const std::vector<std::vector<FlatNode>>& trees) {
    VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    if (trees.size() != vf.numTrees) return false;
    std::vector<double> leafValues;
    for (std::size_t i = 0; i < n; ++i) vf.combinedVariance[i] = 1.0;
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      Tree& tree = vf.trees[j];
      tree.initialize(vf.indexBuffer.data() + j * n, n);
      if (!tree.buildFromFlat(data_, trees[j].data(), trees[j].size(),
                              leafValues))
        return false;
      tree.repartitionSubtree(data_, 0);
      double* factor = vf.factorByTree.data() + j * n;
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      for (int32_t b : tree.bottomScratch) {
        double h = leafValues[static_cast<std::size_t>(b)];
        const Node& node = tree.at(b);
        for (std::size_t m = node.begin; m < node.end; ++m) {
          std::size_t i = tree.indices[m];
          factor[i] = h;
          vf.combinedVariance[i] *= h;
        }
      }
    }
    return true;
  }

  /// Re-length the seven n-sized allocations a variance forest pins at the
  /// creation count - meanWeights_ and the forest's indexBuffer, factorByTree,
  /// combinedVariance, meanResidual, divisor, treeResidual - for a replacement
  /// data set. Every one is scratch or re-derived: refreshVarianceForest,
  /// which the caller must run next, scatters factorByTree through the new
  /// partition and rebuilds combinedVariance from it, and each sweep fills the
  /// other four before reading them. combinedVarianceTest belongs to
  /// resizeTestStorage.
  void resizeVarianceStorage(std::size_t n) {
    VarianceForest& vf = *varianceForest_;
    vf.indexBuffer.assign(n * vf.numTrees, 0);
    vf.factorByTree.assign(n * vf.numTrees, 1.0);
    vf.combinedVariance.assign(n, 1.0);
    vf.meanResidual.assign(n, 0.0);
    vf.divisor.assign(n, 1.0);
    vf.treeResidual.assign(n, 0.0);
    meanWeights_.assign(n, 0.0);
  }

  /// Re-anchor the variance forest to the store after a mutation moved the
  /// predictors, the cut grid, or the observation count: the scale analogue of
  /// the forests_ body of forceRefreshTrees and applyNewData, which loop
  /// forests_ and so never reach it. oldCutPoints remaps split indices onto a
  /// rebuilt grid (null when the grid is installed rather than rebuilt, as
  /// setCutPoints and a forced predictor swap leave it).
  ///
  /// The order is load-bearing. Recovery reads the per-observation slab through
  /// each leaf's CURRENT members, so it must precede every partition change;
  /// recovering afterwards reads the new partition's members out of the old
  /// leaves' slots. Everything downstream then works on node-indexed factors,
  /// exactly as the mean forest's parameters do.
  ///
  /// recoveredFactors supplies those node-indexed factors from an EARLIER
  /// recovery, which is what a whole-data replacement needs: the slab's stride
  /// is the old observation count and the store no longer reports it. Null
  /// recovers here instead. numObservationsChanged says the caller has already
  /// resized the storage (resizeVarianceStorage), so the trees must be
  /// re-pointed at the moved index buffer; the paths that install a grid over
  /// unchanged data leave it false.
  void refreshVarianceForest(
      const std::vector<std::vector<double>>* oldCutPoints,
      const TreeParameters* recoveredFactors = nullptr,
      bool numObservationsChanged = false) {
    VarianceForest& vf = *varianceForest_;
    std::size_t n = data_.numObservations;
    std::vector<double> leafValues;
    for (std::size_t j = 0; j < vf.numTrees; ++j) {
      Tree& tree = vf.trees[j];
      // the empty-leaf veto admits no unoccupied bottom into live state, so
      // every factor recovered HERE is a drawn (positive) scale. A
      // caller supplying recoveredFactors read them itself and owns that
      // invariant: setData's from its own live partition, the warm start's
      // from a donor's flat trees, whose partition is not built until the
      // repartition below.
      assert(recoveredFactors != nullptr || tree.bottomNodesAreOccupied());
      if (recoveredFactors != nullptr)
        leafValues = (*recoveredFactors)[j];
      else
        recoverVarianceLeafValues(vf, j, leafValues);
      tree.dropStaleMissingDirections(data_);
      if (oldCutPoints != nullptr)
        tree.mapOldCutPointsOntoNew<GeometricMerge>(data_, *oldCutPoints,
                                                    leafValues);
      if (numObservationsChanged)
        tree.resetObservations(vf.indexBuffer.data() + j * n, n);
      tree.repartitionSubtree(data_, 0);
      tree.collapseEmptyNodes<GeometricMerge>(
        data_, response_->workingWeights(), leafValues);
      // scatter the surviving factors back through the new partition
      double* factor = vf.factorByTree.data() + j * n;
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      for (int32_t b : tree.bottomScratch) {
        double h = leafValues[static_cast<std::size_t>(b)];
        const Node& node = tree.at(b);
        for (std::size_t m = node.begin; m < node.end; ++m)
          factor[tree.indices[m]] = h;
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      double product = 1.0;
      for (std::size_t j = 0; j < vf.numTrees; ++j)
        product *= vf.factorByTree[j * n + i];
      vf.combinedVariance[i] = product;
    }
  }

  /// treeY <- the residual tree t owns, admitting tree t's old fits and (t > 0)
  /// retiring tree t-1's new fits in observation order. The constant leaf
  /// gathers each tree's fit through mu[leafOf]; the dense slab reads it direct.
  ///
  /// The constant-leaf gathers are unrolled by 4 (n % 4 prologue first) for
  /// ResidT = double and MUST STAY THAT WAY: at -O2 the rolled loop compiles
  /// scalar on every ISA - no compiler emits a hardware gather at any flag
  /// level - and the unroll is what lets gcc assemble a 128-bit (256-bit under
  /// -mavx) software gather plus a packed add, and lets clang pair the index
  /// loads. On x86 the t > 0 body vectorizes only its add; its subtractions
  /// stay scalar. Unrolling reorders iterations, not the operations within an
  /// iteration, so every element is bit-for-bit what the rolled form writes.
  /// ResidT = float deliberately keeps the rolled loop: clang already
  /// vectorizes it with NEON lane-insert gathers, and the unroll would only
  /// bolt a scalar prologue onto that. Measured x86 run -3.0 to -3.9%; see
  /// docs/plans/setpredictor-leafof-rebuild.md.
  ///
  /// leaf and leafPrev are two rows of one allocation; __restrict is sound on
  /// them because both are read-only through the block, not because the rows
  /// are disjoint.
  void rollTreeResidual(Forest<L, ResidT>& forest, size_t t, const double* forestY) {
    size_t n = data_.numObservations;
    // the delta/total is formed in double from double inputs (mu, y, total);
    // the fp32 rounding, when ResidT = float, happens ONLY on the store here -
    // static_cast<double> is the identity default path (byte-identical)
    ResidT* __restrict resid = forest.treeY.data();
    if constexpr (leafIsConstant) {
      const double* __restrict mu = forest.muByTree[t].data();
      const std::uint32_t* __restrict leaf = forest.leafOf.data() + t * n;
#ifndef NDEBUG
      // mu and leafOf are two halves of one cached-fits evaluator; resizing
      // one without the other reads past mu's .size() but inside its capacity,
      // which returns stale values instead of faulting and which the local
      // ASAN build cannot see (detect_container_overflow is off). Mirrors the
      // fused pass's map assert. R's build defines NDEBUG, so this is live
      // only in tests/cpp.
      for (size_t j = 0; j < n; ++j)
        assert(leaf[j] < forest.muByTree[t].size());
#endif
      if (t == 0) {
        const double* __restrict y_ = forestY;
        const double* __restrict total = forest.totalFits.data();
        if constexpr (std::is_same_v<ResidT, double>) {
          size_t i = 0, nMod4 = n % 4;
          for ( ; i < nMod4; ++i)
            resid[i] = static_cast<ResidT>(y_[i] - total[i] + mu[leaf[i]]);
          for ( ; i < n; i += 4) {
            resid[i] = static_cast<ResidT>(y_[i] - total[i] + mu[leaf[i]]);
            resid[i + 1] =
              static_cast<ResidT>(y_[i + 1] - total[i + 1] + mu[leaf[i + 1]]);
            resid[i + 2] =
              static_cast<ResidT>(y_[i + 2] - total[i + 2] + mu[leaf[i + 2]]);
            resid[i + 3] =
              static_cast<ResidT>(y_[i + 3] - total[i + 3] + mu[leaf[i + 3]]);
          }
        } else {
          for (size_t i = 0; i < n; ++i)
            resid[i] = static_cast<ResidT>(y_[i] - total[i] + mu[leaf[i]]);
        }
      } else {
        const double* __restrict muPrev = forest.muByTree[t - 1].data();
        const std::uint32_t* __restrict leafPrev =
          forest.leafOf.data() + (t - 1) * n;
#ifndef NDEBUG
        for (size_t j = 0; j < n; ++j)
          assert(leafPrev[j] < forest.muByTree[t - 1].size());
#endif
        if constexpr (std::is_same_v<ResidT, double>) {
          size_t i = 0, nMod4 = n % 4;
          for ( ; i < nMod4; ++i)
            resid[i] = static_cast<ResidT>(static_cast<double>(resid[i]) +
                                           (mu[leaf[i]] - muPrev[leafPrev[i]]));
          for ( ; i < n; i += 4) {
            resid[i] = static_cast<ResidT>(static_cast<double>(resid[i]) +
                                           (mu[leaf[i]] - muPrev[leafPrev[i]]));
            resid[i + 1] =
              static_cast<ResidT>(static_cast<double>(resid[i + 1]) +
                                  (mu[leaf[i + 1]] - muPrev[leafPrev[i + 1]]));
            resid[i + 2] =
              static_cast<ResidT>(static_cast<double>(resid[i + 2]) +
                                  (mu[leaf[i + 2]] - muPrev[leafPrev[i + 2]]));
            resid[i + 3] =
              static_cast<ResidT>(static_cast<double>(resid[i + 3]) +
                                  (mu[leaf[i + 3]] - muPrev[leafPrev[i + 3]]));
          }
        } else {
          for (size_t i = 0; i < n; ++i)
            resid[i] = static_cast<ResidT>(static_cast<double>(resid[i]) +
                                           (mu[leaf[i]] - muPrev[leafPrev[i]]));
        }
      }
    } else {
      const double* __restrict treeFits = forest.treeFits.data() + t * n;
      if (t == 0) {
        const double* __restrict y_ = forestY;
        const double* __restrict total = forest.totalFits.data();
        for (size_t i = 0; i < n; ++i)
          resid[i] = static_cast<ResidT>(y_[i] - total[i] + treeFits[i]);
      } else {
        const double* __restrict prevFits = treeFits - n;
        for (size_t i = 0; i < n; ++i)
          resid[i] = static_cast<ResidT>(static_cast<double>(resid[i]) +
                                         (treeFits[i] - prevFits[i]));
      }
    }
  }

  /// The residual roll and the pre-move node-average suffstat in ONE
  /// observation-order pass: roll resid[i] exactly as rollTreeResidual does,
  /// then scatter-add it into acc[leafOf[i]], a node-indexed accumulator small
  /// enough to stay L1-resident, so setNodeAverages' random gather over
  /// indices[] never runs. That gather is 39-40% of a sweep
  /// (docs/design/memory-wall-frontier.md sec 10); the fusion measures 1.41x
  /// to 1.54x on Zen2 and 1.07x to 1.25x on M1 Max (secs 11-12).
  ///
  /// Returns false, having written nothing, when the fusion is not eligible;
  /// the caller then runs the stock rollTreeResidual + Tree::setNodeAverages
  /// pair, which stays the only path for every other caller of either.
  /// Eligible iff every clause holds, each compile-time or O(1):
  ///   leafIsConstant     - the scatter is into acc[leafOf[i]], and only the
  ///                        constant leaf keeps a leafOf map
  ///   ResidT == double   - the fp32 roll is deliberately left rolled, for the
  ///                        codegen reason rollTreeResidual documents
  ///   weights == nullptr - a weighted statistic needs a second bank set (sum
  ///                        w as well as sum w r), unmeasured. This one clause
  ///                        is what declines BCF, multinomial, logistic,
  ///                        negbin, t and the heteroscedastic mean forest,
  ///                        whose combiners and families always supply weights
  ///   !leafOfStale[t]    - a stale map still describes the PREVIOUS
  ///                        partition, so acc would be indexed by node ids the
  ///                        current tree need not have
  ///
  /// Exactness contract. resid[i] is formed with the SAME expressions in the
  /// SAME order as rollTreeResidual, so treeY comes out bit-for-bit identical
  /// and the suffstat association is the only draw change anywhere - do not
  /// "simplify" the arithmetic below. The association itself is: bank
  /// assignment positional (the n % 4 prologue accumulates into bank 0, and
  /// thereafter element i goes to bank (i - n % 4) mod 4), the combine
  /// ((b0 + b1) + b2) + b3 left to right always, and the pass scalar
  /// throughout - no SIMD reduce, no reassociation, and structurally no FMA
  /// contraction on an accumulate that holds no multiply. The order therefore
  /// depends only on n and the partition, not on ISA, lane width, or worker
  /// count: within a host, draws stay bitwise identical across every SIMD
  /// dispatch level and thread count.
  bool rollAndSetNodeAveragesFused(Forest<L, ResidT>& forest, size_t t,
                                   const double* forestY,
                                   const double* forestWeights) {
    if constexpr (!leafIsConstant || !std::is_same_v<ResidT, double>) {
      (void) forest; (void) t; (void) forestY; (void) forestWeights;
      return false;
    } else {
      if (forestWeights != nullptr || forest.leafOfStale[t] != 0) return false;

      static_assert(fusedSuffstatBanks == 4,
                    "the unrolled body below assigns four banks by hand");
      size_t n = data_.numObservations;
      Tree& tree(forest.trees[t]);
      size_t numNodes = tree.nodes.size();
      fusedAcc_.assign(numNodes * fusedSuffstatBanks, 0.0);
      double* __restrict acc = fusedAcc_.data();

      double* __restrict resid = forest.treeY.data();
      const double* __restrict mu = forest.muByTree[t].data();
      const std::uint32_t* __restrict leaf = forest.leafOf.data() + t * n;
      size_t i = 0, nMod4 = n % 4;

#ifndef NDEBUG
      // The map addresses acc and then the bottoms; a map that is fresh but
      // wrong scatters into arena slots no bottom reads back. R's build
      // defines NDEBUG, so this is live only in tests/cpp.
      for (size_t j = 0; j < n; ++j)
        assert(leaf[j] < numNodes &&
               tree.at(static_cast<int32_t>(leaf[j])).isBottom());
#endif

      if (t == 0) {
        const double* __restrict y_ = forestY;
        const double* __restrict total = forest.totalFits.data();
        for ( ; i < nMod4; ++i) {
          double r = y_[i] - total[i] + mu[leaf[i]];
          resid[i] = r;
          acc[leaf[i]] += r;
        }
        for ( ; i < n; i += 4) {
          double r0 = y_[i] - total[i] + mu[leaf[i]];
          double r1 = y_[i + 1] - total[i + 1] + mu[leaf[i + 1]];
          double r2 = y_[i + 2] - total[i + 2] + mu[leaf[i + 2]];
          double r3 = y_[i + 3] - total[i + 3] + mu[leaf[i + 3]];
          resid[i] = r0;
          resid[i + 1] = r1;
          resid[i + 2] = r2;
          resid[i + 3] = r3;
          acc[leaf[i]] += r0;
          acc[numNodes + leaf[i + 1]] += r1;
          acc[2 * numNodes + leaf[i + 2]] += r2;
          acc[3 * numNodes + leaf[i + 3]] += r3;
        }
      } else {
        const double* __restrict muPrev = forest.muByTree[t - 1].data();
        const std::uint32_t* __restrict leafPrev =
          forest.leafOf.data() + (t - 1) * n;
        for ( ; i < nMod4; ++i) {
          double r = resid[i] + (mu[leaf[i]] - muPrev[leafPrev[i]]);
          resid[i] = r;
          acc[leaf[i]] += r;
        }
        for ( ; i < n; i += 4) {
          double r0 = resid[i] + (mu[leaf[i]] - muPrev[leafPrev[i]]);
          double r1 =
            resid[i + 1] + (mu[leaf[i + 1]] - muPrev[leafPrev[i + 1]]);
          double r2 =
            resid[i + 2] + (mu[leaf[i + 2]] - muPrev[leafPrev[i + 2]]);
          double r3 =
            resid[i + 3] + (mu[leaf[i + 3]] - muPrev[leafPrev[i + 3]]);
          resid[i] = r0;
          resid[i + 1] = r1;
          resid[i + 2] = r2;
          resid[i + 3] = r3;
          acc[leaf[i]] += r0;
          acc[numNodes + leaf[i + 1]] += r1;
          acc[2 * numNodes + leaf[i + 2]] += r2;
          acc[3 * numNodes + leaf[i + 3]] += r3;
        }
      }

      // the bottoms take the accumulated sums in place of the gather;
      // unweighted, so sumWeights is the count the partition already knows -
      // which is what misc_computeIndexedSufficientStatisticsFast reports too
      tree.bottomScratch.clear();
      tree.fillBottom(0, tree.bottomScratch);
      for (int32_t b : tree.bottomScratch) {
        Node& node(tree.at(b));
        node.sumWeights = static_cast<double>(node.numObservations());
        node.sumWeightedResponse =
          combineFusedSuffstatBanks(acc, static_cast<size_t>(b), numNodes);
      }
      ++fusedSuffstatRuns_;
      return true;
    }
  }

  /// Rebuild totalFits after the sweep loop: the last tree's new fits retire
  /// here instead of in a pass of their own.
  ///
  /// The constant-leaf gather is unrolled by 4 (n % 4 prologue first) for
  /// ResidT = double for the reason rollTreeResidual documents at length: the
  /// rolled form never vectorizes at -O2 on any ISA, the unroll is what buys
  /// the packed software gather, and it is elementwise so it stays
  /// bit-for-bit identical to the rolled form. Do not roll it back up.
  /// ResidT = float keeps the rolled loop so clang's own codegen there is
  /// left exactly as it is.
  void finalizeTotalFits(Forest<L, ResidT>& forest, const double* forestY) {
    if (forest.numTrees == 0) return;
    size_t n = data_.numObservations;
    const size_t last = forest.numTrees - 1;
    const double* __restrict y_ = forestY;
    const ResidT* __restrict resid = forest.treeY.data();
    double* __restrict total = forest.totalFits.data();
    if constexpr (leafIsConstant) {
      const double* __restrict mu = forest.muByTree[last].data();
      const std::uint32_t* __restrict leaf = forest.leafOf.data() + last * n;
      if constexpr (std::is_same_v<ResidT, double>) {
        size_t i = 0, nMod4 = n % 4;
        for ( ; i < nMod4; ++i) total[i] = y_[i] - resid[i] + mu[leaf[i]];
        for ( ; i < n; i += 4) {
          total[i] = y_[i] - resid[i] + mu[leaf[i]];
          total[i + 1] = y_[i + 1] - resid[i + 1] + mu[leaf[i + 1]];
          total[i + 2] = y_[i + 2] - resid[i + 2] + mu[leaf[i + 2]];
          total[i + 3] = y_[i + 3] - resid[i + 3] + mu[leaf[i + 3]];
        }
      } else {
        for (size_t i = 0; i < n; ++i) total[i] = y_[i] - resid[i] + mu[leaf[i]];
      }
    } else {
      const double* __restrict lastFits = forest.treeFits.data() + last * n;
      for (size_t i = 0; i < n; ++i)
        total[i] = y_[i] - resid[i] + lastFits[i];
    }
  }

  /// totalFits += tree t's fits (constant: gathered through mu[leafOf]).
  void addTreeFitsToTotal(Forest<L, ResidT>& forest, size_t t) {
    size_t n = data_.numObservations;
    if constexpr (leafIsConstant) {
      const double* mu = forest.muByTree[t].data();
      const std::uint32_t* leaf = forest.leafOf.data() + t * n;
      double* total = forest.totalFits.data();
      for (size_t i = 0; i < n; ++i) total[i] += mu[leaf[i]];
    } else {
      misc_addVectorsInPlace(forest.treeFits.data() + t * n, n,
                             forest.totalFits.data());
    }
  }

  /// totalFits -= tree t's current fits, the inverse of addTreeFitsToTotal.
  void subtractTreeFitsFromTotal(Forest<L, ResidT>& forest, size_t t) {
    size_t n = data_.numObservations;
    if constexpr (leafIsConstant) {
      const double* mu = forest.muByTree[t].data();
      const std::uint32_t* leaf = forest.leafOf.data() + t * n;
      double* total = forest.totalFits.data();
      for (size_t i = 0; i < n; ++i) total[i] -= mu[leaf[i]];
    } else {
      misc_subtractVectorsInPlace(forest.treeFits.data() + t * n, n,
                                  forest.totalFits.data());
    }
  }

  /// Leaf parameters recovered for tree t, indexed by arena node id. The
  /// constant leaf's persistent mu table already holds them; the resize is a
  /// no-op fit that keeps the arena-length invariant callers rely on.
  void recoverParametersFromFits(Forest<L, ResidT>& forest, size_t t,
                                 std::vector<double>& paramByNode) {
    paramByNode = forest.muByTree[t];
    paramByNode.resize(forest.trees[t].nodes.size(), 0.0);
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

  /// Constant leaves install the mu table only: callers whose partitions
  /// changed rebuild leafOf themselves, and the rest keep the current map.
  void setTreeFitsFromParameters(Forest<L, ResidT>& forest, size_t t,
                                 const std::vector<double>& paramByNode) {
    if constexpr (leafIsConstant) {
      forest.muByTree[t] = paramByNode;
    } else {
      // function-leaf cold-start over a fresh partition: scatter the per-node
      // value into the dense slab
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
  }

  /// The vector-parameter sibling of setTreeFitsFromParameters: each member
  /// observation's fit evaluates its leaf's block against the leaf model's
  /// current covariates.
  void setTreeFitsFromParameterBlocks(Forest<L, ResidT>& forest, size_t t,
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

  /// Keep the constrained leaf's mu block consistent with an accepted move so
  /// the following draw sweeps a feasible, correctly sized state. A stale tree
  /// resets to the all-equal (feasible) zero vector; an accepted birth grows
  /// the block and seeds the two children with the split parent's value
  /// (feasible for both and ordered); an accepted death seeds the merged leaf
  /// with a point inside its neighbor bounds. Rejections leave the block, which
  /// still matches the restored structure.
  void maintainMonotoneLeafStore(Forest<L, ResidT>& forest, size_t t, bool wasStale,
                                 bool stepTaken, StepType stepType,
                                 int32_t changedNode) {
    Tree& tree(forest.trees[t]);
    std::vector<double>& mu(forest.muByTree[t]);
    if (wasStale) {
      mu.assign(tree.nodes.size(), 0.0);
      return;
    }
    if (!stepTaken) return;
    if (stepType == StepType::birth) {
      // grow the block to the new node count, then draw the two children from
      // their exact constrained conditional (docs/design/monotone.md eq. 4.17);
      // the sibling slots are
      // read but not their stale mu, so growing with zeros before the draw is
      // safe
      mu.resize(tree.nodes.size(), 0.0);
      forest.leaf.redrawAfterBirth(rng_, tree, changedNode, forest.k,
                                   sigma_ * sigma_, mu.data());
    } else if (stepType == StepType::death) {
      mu.resize(tree.nodes.size(), 0.0);
      forest.leaf.redrawAfterDeath(rng_, tree, changedNode, forest.k,
                                   sigma_ * sigma_, mu.data());
    }
  }

  void sampleParametersAndSetFits(Forest<L, ResidT>& forest, size_t t, double* fits,
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
          int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
          forest.currTestFits[i] =
            forest.leaf.fitForTestObservationForNode(tree, leafIndex, i);
        });
    } else if constexpr (!L::hasVectorParams) {
      // the constant leaf fills its mu table only; leafOf is already current
      // (patched at move acceptance, or rebuilt when marked stale)
      (void) fits;
      std::vector<double>& mu(forest.muByTree[t]);
      if constexpr (TreeDrawLeafModel<L>) {
        // the surviving block carries the previous sweep's feasible draw (the
        // move phase kept it sized and feasible); grow it for any node the tree
        // gained without disturbing existing leaves, then the coupled Gibbs
        // sweep updates every leaf in place. Fixed k only, the truncated draw
        // carries no clean chi-k statistic (docs/design/monotone.md section 6)
        mu.resize(tree.nodes.size(), 0.0);
        forest.leaf.drawParametersForTree(rng_, tree, bottoms, forest.k,
                                          sigma_ * sigma_, mu.data());
      } else {
        mu.assign(tree.nodes.size(), 0.0);
        for (int32_t i : bottoms) {
          const Node& node(tree.at(i));
          double param = node.numObservations() == 0
            ? 0.0
            : forest.leaf.drawFromPosteriorForNode(rng_, tree, forest.k,
                                                   sigma_ * sigma_, i);
          mu[static_cast<size_t>(i)] = param;

          // a forced-zero empty leaf is not a draw from the k-scaled prior, so
          // it carries no information about k; skip it as the function path does
          if (forest.updateK && node.numObservations() > 0) {
            forest.kSumSquaredParams += param * param;
            forest.kNumLeaves += 1.0;
          }
        }
      }

      if (updateTestFits && data_.numTestObservations > 0)
        routeTestRows(data_.numTestObservations, [&](size_t i) {
          int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
          forest.currTestFits[i] = mu[static_cast<size_t>(leafIndex)];
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

        // a forced-zero empty leaf is not a draw from the k-scaled prior, so
        // it carries no information about k; skip it as the function path does
        if (forest.updateK && node.numObservations() > 0) {
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
          int32_t leafIndex = tree.findBottomNodeForRow(data_, i);
          forest.currTestFits[i] = forest.leaf.fitForTestObservation(
            treeParams.data() + static_cast<size_t>(leafIndex) * numParams,
            i);
        });
    }
  }

  /// Sample sd of the range-scaled working response, the anchor the BCF
  /// calibration map (docs/design/bcf.md) states its per-forest leaf scales
  /// against.
  double scaledResponseSd() const {
    std::size_t n = data_.numObservations;
    const double* yScaled = response_->workingResponse();
    double mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) mean += yScaled[i];
    mean /= static_cast<double>(n);
    double sumSquares = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double d = yScaled[i] - mean;
      sumSquares += d * d;
    }
    return std::sqrt(sumSquares / static_cast<double>(n - 1));
  }

  /// Variant A block-additive install (docs/design/interaction-constraints.md):
  /// confine each whole tree to one declared group of predictors so the ensemble
  /// is exactly f = sum_G f_G. Build the numBlocks group-membership rows (each
  /// intersected with any base columnMask already on the forest - the F3 BCF-tau
  /// case, so the moderator restriction is never lost), assign trees to groups by
  /// the deterministic contiguous capacity (trees [0, c0) -> group 0, [c0, c0+c1)
  /// -> group 1, ...; NO rng), and point every tree's column mask at its group's
  /// row. A no-op when numBlocks is 0, so the default path is byte-for-byte
  /// unchanged. blockOfColumn[j] < 0 marks a column in no block (masked out).
  static void installBlockMasks(Forest<L, ResidT>& forest,
                                std::size_t numPredictors, std::size_t numBlocks,
                                const std::int32_t* blockOfColumn,
                                const std::size_t* blockTreeCounts) {
    if (numBlocks == 0) return;
    const bool haveBase = !forest.columnMask.empty();
    forest.blockMasks.assign(numBlocks * numPredictors, 0);
    for (std::size_t g = 0; g < numBlocks; ++g)
      for (std::size_t j = 0; j < numPredictors; ++j) {
        const bool inGroup = blockOfColumn[j] >= 0 &&
          static_cast<std::size_t>(blockOfColumn[j]) == g;
        const bool baseAllows = !haveBase || forest.columnMask[j] != 0;
        forest.blockMasks[g * numPredictors + j] =
          (inGroup && baseAllows) ? 1 : 0;
      }
    forest.blockOfTree.assign(forest.numTrees, 0);
    std::size_t t = 0;
    for (std::size_t g = 0; g < numBlocks && t < forest.numTrees; ++g)
      for (std::size_t c = 0; c < blockTreeCounts[g] && t < forest.numTrees; ++c)
        forest.blockOfTree[t++] = g;
    for (std::size_t tt = 0; tt < forest.numTrees; ++tt)
      forest.trees[tt].setColumnMask(
        forest.blockMasks.data() + forest.blockOfTree[tt] * numPredictors);
  }

  /// A BCF forest, built self-contained so the single-forest constructor
  /// (covered by the draw-for-draw equivalence benchmark) is untouched:
  /// constant leaf, fixed k = 1 (the map's convention), no DART. nodeScale is
  /// the map-derived total.
  void buildBCFForest(const BCFForestSpec& spec, double nodeScale) {
    std::size_t n = data_.numObservations;
    forests_.emplace_back();
    Forest<L, ResidT>& forest = forests_.back();
    forest.numTrees = spec.numTrees;
    forest.birthOrDeathProbability = spec.birthOrDeathProbability;
    forest.swapProbability = spec.swapProbability;
    forest.changeProbability = spec.changeProbability;
    forest.birthProbability = spec.birthProbability;
    forest.updateK = false;
    forest.useDart = false;
    forest.k = 1.0;
    forest.leaf.scale =
      nodeScale / std::sqrt(static_cast<double>(spec.numTrees));
    forest.treePrior.base = spec.base;
    forest.treePrior.power = spec.power;
    forest.indexBuffer.resize(n * spec.numTrees);
    forest.trees.resize(spec.numTrees);
    for (std::size_t t = 0; t < spec.numTrees; ++t)
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);

    // A restricted forest clears the availability of every unlisted column,
    // installed on its trees exactly as the single-forest ctor does; an empty
    // list leaves the mask empty and the trees unrestricted (the default
    // availability path, byte-for-byte). Only tau carries a list; mu reads the
    // full store.
    if (spec.columns != nullptr && spec.numColumns > 0) {
      forest.columnMask.assign(data_.numPredictors, 0);
      for (std::size_t c = 0; c < spec.numColumns; ++c)
        forest.columnMask[spec.columns[c]] = 1;
      for (std::size_t t = 0; t < spec.numTrees; ++t)
        forest.trees[t].setColumnMask(forest.columnMask.data());
    }

    // this forest's optional interaction constraint (installed like the single-
    // forest ctor): mu and tau carry independent caps. An unset (or inactive)
    // constraint leaves every tree's pointer null and the path unchanged.
    if (spec.interactionMaxOrder > 0 || spec.interactionNumForbiddenPairs > 0) {
      forest.interaction = std::make_unique<InteractionConstraint>();
      forest.interaction->build(data_.numPredictors, spec.interactionMaxOrder,
                                spec.interactionForbiddenPairs,
                                spec.interactionNumForbiddenPairs);
      if (forest.interaction->active())
        for (std::size_t t = 0; t < spec.numTrees; ++t)
          forest.trees[t].setInteractionConstraint(forest.interaction.get());
      else
        forest.interaction.reset();
    }

    // Variant A block-additive constraint (mu / tau independent): the block rows
    // intersect tau's moderator columnMask installed above (F3), so a restricted
    // tau's per-tree mask is group AND moderators. numBlocks 0 leaves it unchanged.
    installBlockMasks(forest, data_.numPredictors, spec.numBlocks,
                      spec.blockOfColumn, spec.blockTreeCounts);

    initForestFitStorage(forest, n);
    forest.totalFits.assign(n, 0.0);
    forest.treeY.resize(n);
  }

  /// One symmetric category forest for a multinomial chain: constant leaf, no
  /// DART, no split restriction, fixed k, leaf scale nodeScale/sqrt(numTrees)
  /// (the pi*sqrt(3)/sqrt(2) anchor, docs/design/multinomial.md). Every category
  /// forest is identical.
  void buildMultinomialForest(const MultinomialForestSpec& spec,
                              double nodeScale, double k) {
    std::size_t n = data_.numObservations;
    forests_.emplace_back();
    Forest<L, ResidT>& forest = forests_.back();
    forest.numTrees = spec.numTrees;
    forest.birthOrDeathProbability = spec.birthOrDeathProbability;
    forest.swapProbability = spec.swapProbability;
    forest.changeProbability = spec.changeProbability;
    forest.birthProbability = spec.birthProbability;
    forest.updateK = false;
    forest.useDart = false;
    forest.k = k;
    forest.leaf.scale =
      nodeScale / std::sqrt(static_cast<double>(spec.numTrees));
    forest.treePrior.base = spec.base;
    forest.treePrior.power = spec.power;
    forest.indexBuffer.resize(n * spec.numTrees);
    forest.trees.resize(spec.numTrees);
    for (std::size_t t = 0; t < spec.numTrees; ++t)
      forest.trees[t].initialize(forest.indexBuffer.data() + t * n, n);
    initForestFitStorage(forest, n);
    forest.totalFits.assign(n, 0.0);
    forest.treeY.resize(n);
  }

  /// A single forest reports its own total fits directly; a multi-forest chain
  /// asks the combiner for the blended per-observation location. The null test
  /// stays chain-side so the single-forest path returns a bare pointer with no
  /// virtual call and no copy.
  const double* combinedFits() {
    return combiner_ ? combiner_->combinedFits(forests_)
                     : forests_[0].totalFits.data();
  }

  void storeSample(Results& results, size_t sampleNum) {
    // the scalar channels (k, variable counts, split probabilities) and the
    // single-forest fit paths address the reported forest; the combiner names
    // it (BCF: the prognostic mu, forest 0), a single-forest chain is forest 0
    std::size_t reportedIndex = combiner_ ? combiner_->reportedForest() : 0;
    Forest<L, ResidT>& forest = forests_[reportedIndex];
    size_t n = data_.numObservations;
    double scale = response_->fitScale();
    double shift = response_->fitShift();

    if (results.sigma != nullptr)
      results.sigma[sampleNum] = sigma_ * response_->sigmaScale();

    if (results.k != nullptr) results.k[sampleNum] = forest.k;

    // heteroscedastic: the variance surface s^2(x) on the original scale (the
    // working product times sigmaScale^2), train and test, its own channel
    if (varianceForest_) {
      double varScale = response_->sigmaScale() * response_->sigmaScale();
      const VarianceForest& vf = *varianceForest_;
      if (results.varianceFits != nullptr) {
        double* out = results.varianceFits + sampleNum * n;
        for (size_t i = 0; i < n; ++i) out[i] = varScale * vf.combinedVariance[i];
      }
      if (results.varianceTestFits != nullptr && data_.numTestObservations > 0) {
        size_t nTest = data_.numTestObservations;
        double* out = results.varianceTestFits + sampleNum * nTest;
        for (size_t i = 0; i < nTest; ++i)
          out[i] = varScale * vf.combinedVarianceTest[i];
      }
    }

    // the combined output carries one per-observation channel per reported
    // location (one everywhere but a multi-location combiner); the writes stride
    // location-major within a sample so L = 1 is the exact current byte layout
    size_t numLocations = combiner_ ? combiner_->numReportedLocations() : 1;

    if (results.trainingFits != nullptr) {
      double* out = results.trainingFits + sampleNum * n * numLocations;
      // the combiner owns the blend; recompute it against the post-glue scalars
      // this sweep settled on (combinedFits at the sweep top ran before the
      // glue draw). Off any combiner the reported forest's own fits are it.
      const double* combined =
        combiner_ ? combiner_->combinedFits(forests_) : forest.totalFits.data();
      const double* offset = response_->offset();
      for (size_t loc = 0; loc < numLocations; ++loc) {
        double* dst = out + loc * n;
        const double* src = combined + loc * n;
        for (size_t i = 0; i < n; ++i)
          dst[i] = scale * src[i] + shift;
        // original-scale convention, matching the recorded test fits: any
        // offset is part of the fit. A per-observation offset add would be
        // wrong for a multi-location combiner writing probability channels
        // (softmax) rather than an additive location, but multinomial creation
        // carries no offset, so this stays null on that path.
        if (offset != nullptr) misc_addVectorsInPlace(offset, n, dst);
      }
    }

    if (results.testFits != nullptr && data_.numTestObservations > 0) {
      size_t nTest = data_.numTestObservations;
      double* out = results.testFits + sampleNum * nTest * numLocations;
      if (combiner_ && !combiner_->testFitsAreDefined()) {
        // A BCF test blend a * mu + b_z * tau is ill-defined here: the API
        // carries no test treatment vector, so only the bare prognostic
        // forest could be recorded, which silently misreports the fit.
        // Flag the channel as unusable; BCF consumers recombine per forest
        // via forestTotalFits + the bcfGlue coefficients (docs/design/bcf.md).
        for (size_t i = 0; i < nTest * numLocations; ++i)
          out[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        // A single forest (or single-location combiner) reports its own test
        // fits directly; a multi-location combiner (multinomial) asks it to
        // blend the K forests' totalTestFits into the K softmax test
        // probabilities. The level-centering grand shift is common to all K
        // forests and absent from totalTestFits, but softmax invariance to a
        // common shift makes the blend correct without it (afterCombine leaves
        // totalTestFits untouched; docs/design/multinomial.md). At L = 1 testSrc
        // is the bare totalTestFits and the write is byte-identical.
        const double* testSrc = (combiner_ && numLocations > 1)
          ? combiner_->combinedTestFits(forests_)
          : forest.totalTestFits.data();
        for (size_t loc = 0; loc < numLocations; ++loc) {
          double* dst = out + loc * nTest;
          for (size_t i = 0; i < nTest; ++i)
            dst[i] = scale * testSrc[loc * nTest + i] + shift;
          if (data_.testOffset != nullptr)
            misc_addVectorsInPlace(data_.testOffset, nTest, dst);
        }
      }
    }

    // per-forest reporting, for a coupling that defines it (BCF): each forest's
    // own internal-scale function values, forest-major within a sample, and the
    // scalars that recombine them. Both are reads of state this sweep already
    // settled on - the same values forestTotalFits and bcfGlue hand a caller
    // that drives one sweep at a time - so the channels consume no rng and
    // mutate no state.
    if (combiner_ && combiner_->forestReportingIsDefined()) {
      if (results.forestFits != nullptr) {
        double* out = results.forestFits + sampleNum * n * forests_.size();
        for (std::size_t f = 0; f < forests_.size(); ++f)
          std::memcpy(out + f * n, forests_[f].totalFits.data(),
                      n * sizeof(double));
      }
      if (results.glue != nullptr) {
        double* out = results.glue + sampleNum * 3;
        combiner_->bcfGlue(out[0], out[1], out[2]);
      }
    }

    if (results.variableCounts != nullptr) {
      // one varcount slab per reported forest, forest-major within a sample so
      // count 1 (single forest and BCF) is the exact current byte layout; a
      // multi-forest combiner (multinomial) records each category forest's
      // splits into slot j = variableCountForest(j)
      std::size_t numVCForests =
        combiner_ ? combiner_->numVariableCountForests() : 1;
      for (std::size_t j = 0; j < numVCForests; ++j) {
        std::size_t f =
          combiner_ ? combiner_->variableCountForest(j) : reportedIndex;
        forestVariableCounts(
          f, results.variableCounts +
               (sampleNum * numVCForests + j) * data_.numPredictors);
      }
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

    // the K-1 ordinal thresholds, aligned with this sweep's latent draw; the
    // marginal cutpoint MH ran inside refreshLatents before this store, so
    // gamma is the accepted value the drawn latents are consistent with
    if (results.cutpoints != nullptr && response_->carriesCutpoints()) {
      std::size_t numCutpoints = response_->numCutpoints();
      std::memcpy(results.cutpoints + sampleNum * numCutpoints,
                  response_->cutpoints(), numCutpoints * sizeof(double));
    }

    if (results.logLikelihood != nullptr) {
      double* out = results.logLikelihood + sampleNum * n;
      if (combiner_ && !combiner_->logLikelihoodIsDefined()) {
        // The BCF per-observation location blends two forests through the
        // glue coefficients, which the response model cannot see; scoring
        // forest 0 alone would misreport it. NaN-flag as testFits does.
        for (size_t i = 0; i < n; ++i)
          out[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        response_->computeLogLikelihood(forest.totalFits.data(), sigma_, n,
                                        out);
      }
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
  std::vector<Forest<L, ResidT>> forests_;

  // the multi-forest coupling (BCF today); null for every single-forest
  // sampler, so the sweep, reporting, and state paths collapse to the direct
  // forest-0 path when so and pay no virtual call
  std::unique_ptr<ForestCombiner<L, ResidT>> combiner_;

  // caller-supplied per-forest observation weights, one BORROWED pointer per
  // forest (null = none). Left EMPTY until the first install, which is the
  // pass-through gate every configuration that installs none takes. The scratch
  // holds the composed precisions of whichever forest is being formed, so one
  // buffer serves them all: nothing reads a forest's precisions past its own
  // tree loop.
  std::vector<const double*> forestWeights_;
  std::vector<double> forestWeightScratch_;

  // heteroscedastic variance forest (docs/design/heteroscedastic.md); null for
  // every homoscedastic sampler, so its sweep and the weight division are
  // branched out and the mean path is byte-identical. meanWeights_ is the mean
  // forest's per-sweep divided precisions user_w_i / s^2(x_i), sized only when
  // the variance forest is built.
  std::unique_ptr<VarianceForest> varianceForest_;
  std::vector<double> meanWeights_;
  ResidualPrior varianceLeafPrior_;

  // Persistent pool for parallel test-fit routing, sized to this chain's
  // share of the thread budget; created lazily, never below the cutoff. The
  // forests borrow it through routeTestRows.
  misc_mt_manager_t testFitPool_ = nullptr;
  static constexpr size_t testFitParallelCutoff = 65536;

  // Node-indexed scatter-add accumulator for the fused roll, fusedSuffstatBanks
  // copies laid out bank-major. Sized per tree (a handful of nodes), so it
  // stays L1-resident; it grows to the largest tree this chain has swept and
  // then stops reallocating. Chain-owned mutable scratch: it is safe only
  // because a chain's sweep is sequential, and any future in-chain parallelism
  // has to privatize it.
  std::vector<double> fusedAcc_;
  // Diagnostic only, never read by the sampler: how many (tree, sweep) bodies
  // took the fused pass. Every eligibility clause is otherwise a silent
  // decline, which would let a refactor give back the gather unnoticed.
  size_t fusedSuffstatRuns_ = 0;
};

}  // namespace bartcore

#endif  // BARTCORE_CHAIN_HPP
