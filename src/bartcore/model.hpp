#ifndef BARTCORE_MODEL_HPP
#define BARTCORE_MODEL_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <cfloat>
#include <memory>
#include <numbers>
#include <vector>

#include <external/random.h>
#include <misc/stats.h>
#include <misc/linearAlgebra.h>

#include "data.hpp"
#include "tree.hpp"

namespace bartcore {

/// Leaf models the conjugate engine can run must be integrable: they expose
/// a closed-form log marginal over their parameters, scored through node
/// context - (tree, working response, weights, node) - which is all the
/// moves ever see. Parameter shape splits the hierarchy below that: scalar
/// models (one parameter per leaf) keep the scalar draw interface,
/// vector models write numParams() doubles per leaf and evaluate fits per
/// observation, function-valued models draw one value per member observation
/// directly into the fit vector (the fits ARE the parameters). Chain code
/// branches on the shape traits at compile time, so a scalar leaf model
/// never instantiates the vector or function paths.
template <typename L>
concept LeafModelCore =
  requires(const L leaf, const Tree& tree, const double* v, double d,
           std::int32_t node) {
    { L::hasVectorParams } -> std::convertible_to<bool>;
    { L::hasFunctionParams } -> std::convertible_to<bool>;
    { leaf.logIntegratedLikelihoodForNode(tree, v, v, d, d, node) }
      -> std::same_as<double>;
  };

template <typename L>
concept ScalarLeafModel = LeafModelCore<L> && !L::hasVectorParams &&
  !L::hasFunctionParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, double d,
           std::int32_t node) {
    { leaf.drawFromPosteriorForNode(rng, tree, d, d, node) }
      -> std::same_as<double>;
    { leaf.drawFromPrior(rng, d) } -> std::same_as<double>;
  };

template <typename L>
concept VectorLeafModel = LeafModelCore<L> && L::hasVectorParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, const double* v,
           double d, std::int32_t node, double* out, std::size_t i) {
    { leaf.numParams() } -> std::same_as<std::size_t>;
    { leaf.drawFromPosteriorForNode(rng, tree, v, v, d, d, node, out) }
      -> std::same_as<void>;
    { leaf.drawFromPrior(rng, d, out) } -> std::same_as<void>;
    { leaf.fitForObservation(v, i) } -> std::same_as<double>;
    { leaf.fitForTestObservation(v, i) } -> std::same_as<double>;
  };

/// Chi-k accumulation from one function-valued leaf draw: the standardized
/// sum of squares (f' C^-1 f for a GP draw, the squared value for a
/// constant-fallback draw) and its parameter count.
struct FunctionLeafDrawStats {
  double sumSquaredParams = 0.0;
  double numParams = 0.0;
};

/// Function-valued leaves: draws write one value per member observation into
/// the caller's fit vector, so parameters need no per-tree storage. Test
/// rows evaluate through a per-draw cache the chain resets with
/// beginTreeDraw before each tree's parameter sweep.
template <typename L>
concept FunctionLeafModel = LeafModelCore<L> && L::hasFunctionParams &&
  requires(const L leaf, ext_rng* rng, const Tree& tree, const double* v,
           double d, std::int32_t node, double* fits, std::size_t i) {
    { leaf.beginTreeDraw(tree) } -> std::same_as<void>;
    { leaf.drawFromPosteriorForNode(rng, tree, v, v, d, d, node, fits) }
      -> std::same_as<FunctionLeafDrawStats>;
    { leaf.drawFromPriorForNode(rng, tree, d, node, fits) }
      -> std::same_as<FunctionLeafDrawStats>;
    { leaf.fitForTestObservationForNode(tree, node, i) }
      -> std::same_as<double>;
  };

template <typename L>
concept IntegrableLeafModel =
  ScalarLeafModel<L> || VectorLeafModel<L> || FunctionLeafModel<L>;

/// Constant Gaussian leaf: mu ~ N(0, (scale / k)^2), Gaussian likelihood.
struct ConstantGaussianLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = false;

  double scale;  // nodeScale / sqrt(numTrees)

  // The (sum w, sum wz, sum wz^2) sufficient statistic drives the marginal
  // and the draw directly. This is the crossproduct form of the classic CGM
  // formula: its centered response sum of squares is
  // sumWZSq - sumWZ^2 / sumW and its posterior precision is sumW / sigma^2,
  // so the two forms are algebraically equal and differ only in rounding.
  double logIntegratedLikelihood(double k, double residualVariance,
                                 double sumWeights, double sumWeightedResponse,
                                 double sumWeightedResponseSq) const {
    if (sumWeights == 0.0) return 0.0;

    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = sumWeights / residualVariance;
    double mean = sumWeightedResponse / sumWeights;
    double centeredSumOfSquares =
      sumWeightedResponseSq - sumWeightedResponse * mean;

    double result;
    result  = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision));
    result -= 0.5 * centeredSumOfSquares / residualVariance;
    result -= 0.5 * ((priorPrecision * mean) * (posteriorPrecision * mean)) /
              (priorPrecision + posteriorPrecision);
    return result;
  }

  double drawFromPosterior(ext_rng* rng, double k, double sumWeights,
                           double sumWeightedResponse,
                           double residualVariance) const {
    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = sumWeights / residualVariance;
    double posteriorMean = (sumWeightedResponse / residualVariance) /
                           (priorPrecision + posteriorPrecision);
    double posteriorSd = 1.0 / std::sqrt(priorPrecision + posteriorPrecision);
    return posteriorMean + posteriorSd * ext_rng_simulateStandardNormal(rng);
  }

  double drawFromPrior(ext_rng* rng, double k) const {
    return (scale / k) * ext_rng_simulateStandardNormal(rng);
  }

  // The node-context interface the engine drives; delegates to the scalar
  // math above against the node's cached sufficient statistic. y and weights
  // are unused now that the node caches the full triple.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double*,
                                        const double*, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    return logIntegratedLikelihood(k, residualVariance, node.sumWeights,
                                   node.sumWeightedResponse,
                                   node.sumWeightedResponseSq);
  }

  double drawFromPosteriorForNode(ext_rng* rng, const Tree& tree, double k,
                                  double residualVariance,
                                  int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    return drawFromPosterior(rng, k, node.sumWeights,
                             node.sumWeightedResponse, residualVariance);
  }
};

static_assert(ScalarLeafModel<ConstantGaussianLeaf>);

/// Linear Gaussian leaf: each bottom node fits an intercept plus a linear
/// term in q designated ordinal predictor columns, standardized internally
/// to the training mean and sample sd. All q + 1 coefficients are iid
/// N(0, (scale / k)^2), so the marginal over them is a closed-form ridge
/// regression and the conjugate MH moves apply unchanged. The sufficient
/// statistics (U'WU, U'Wz, z'Wz) accumulate per call over the node's index
/// segment - the same O(n_leaf) walk as the constant leaf's variance,
/// widened by (q + 1)^2 - so move rollbacks need no cache invalidation.
/// Missing covariate values enter at the standardized mean (zero); rules on
/// the same column still route the missingness itself. Linear leaves read
/// raw predictor values through ColumnStore::rawColumn: the borrowed matrix
/// normally, gathered copies (with parent-derived standardization constants)
/// on a view built with the designation; a view without them is refused
/// upstream.
struct LinearGaussianLeaf {
  static constexpr bool hasVectorParams = true;
  static constexpr bool hasFunctionParams = false;
  /// Sufficient-statistic scratch lives on the stack; the factory validates.
  static constexpr std::size_t maxNumCovariates = 8;

  double scale = 1.0;  // nodeScale / sqrt(numTrees)

  std::size_t numParams() const { return numCovariates_ + 1; }
  std::size_t numCovariates() const { return numCovariates_; }
  const std::vector<std::size_t>& covariateColumns() const { return columns_; }
  const std::vector<double>& covariateMeans() const { return means_; }
  const std::vector<double>& covariateSds() const { return sds_; }

  /// Gather and standardize the designated columns from the store's raw
  /// training values. Standardization constants come from the observed
  /// values alone; a constant (or all-missing) column keeps sd 1 and
  /// degrades to an extra intercept the ridge prior absorbs.
  void initialize(const ColumnStore& data, const std::size_t* columns,
                  std::size_t numColumns, std::size_t numChains = 1) {
    numCovariates_ = numColumns;
    columns_.assign(columns, columns + numColumns);
    statisticsCacheBudget_ =
      statisticsCacheTotalBudgetBytes / (numChains == 0 ? 1 : numChains);
    reinitialize(data);
  }

  /// Re-derive the standardization constants and covariates for the stored
  /// designation, for whole-data replacement (possibly a new number of
  /// observations): the analogue of setData rebuilding the cut grid. On a
  /// view the raw values come from the gathered copies and the constants
  /// from the parent's full data - the same calibration inheritance as the
  /// copied cut grid, so every fold runs the prior a full-data fit would.
  void reinitialize(const ColumnStore& data) {
    clearStatisticsCache();
    std::size_t numColumns = numCovariates_;
    numObservations_ = data.numObservations;
    means_.assign(numColumns, 0.0);
    sds_.assign(numColumns, 1.0);
    u_.resize(numObservations_ * numColumns);
    for (std::size_t j = 0; j < numColumns; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double mean, sd;
      if (!data.suppliedStandardization(columns_[j], &mean, &sd))
        standardizationMomentsForColumn(column, numObservations_, &mean, &sd);
      means_[j] = mean;
      sds_[j] = sd;
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - mean) / sds_[j];
    }
    rebuildTestCovariates(data);
  }

  /// Regather the training covariates under the existing standardization
  /// constants, for the predictor-mutation surface (including rollback,
  /// which restores the old raw values): the prior's calibration is sticky
  /// across in-place value changes, the way refreshCutsForColumn keeps the
  /// cut count. Whole-data replacement re-initializes instead, refreshing
  /// the constants the way setData rebuilds the cut grid.
  void regatherTrainingCovariates(const ColumnStore& data) {
    clearStatisticsCache();
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Drop the crossproduct cache when U'WU's other inputs change: the case
  /// weights (setWeights) or the per-sweep Polya-Gamma refresh of a latent
  /// family. Covariate and whole-data mutations clear it through the two
  /// regather paths above.
  void invalidateStatistics() const { clearStatisticsCache(); }

  /// Regather the test covariates under the training standardization; called
  /// whenever the store's test data changes.
  void rebuildTestCovariates(const ColumnStore& data) {
    numTestObservations_ = data.numTestObservations;
    uTest_.assign(numTestObservations_ * numCovariates_, 0.0);
    if (numTestObservations_ == 0) return;
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawTestColumn(columns_[j]);
      if (column == nullptr) continue;
      double* u = uTest_.data() + j * numTestObservations_;
      for (std::size_t i = 0; i < numTestObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  double fitForObservation(const double* params, std::size_t i) const {
    double result = params[0];
    for (std::size_t j = 0; j < numCovariates_; ++j)
      result += params[j + 1] * u_[i + j * numObservations_];
    return result;
  }

  double fitForTestObservation(const double* params, std::size_t i) const {
    double result = params[0];
    for (std::size_t j = 0; j < numCovariates_; ++j)
      result += params[j + 1] * uTest_[i + j * numTestObservations_];
    return result;
  }

  /// log integral of prod N(z_i; b0 + b'u_i, sigma^2 / w_i) against the
  /// N(0, (scale / k)^2 I) prior on (b0, b), dropping the same terms the
  /// constant leaf drops (they depend only on n and w, which the moves'
  /// before/after comparisons share). With M = U'WU + tau sigma^2 I and
  /// tau = (k / scale)^2:
  ///   0.5 (q+1) log(tau sigma^2) - 0.5 log det M
  ///     - 0.5 (z'Wz - b'M^-1 b) / sigma^2,  b = U'Wz,
  /// which reduces exactly to the constant leaf's formula at q = 0.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    if (tree.at(nodeIndex).numObservations() == 0) return 0.0;

    std::size_t p = numParams();
    double crossproduct[maxStatisticSize], projection[maxNumCovariates + 1];
    double responseSumOfSquares;
    accumulateNodeStatistics(tree, y, weights, nodeIndex, crossproduct,
                             projection, &responseSumOfSquares);

    double ridge = (k / scale) * (k / scale) * residualVariance;
    for (std::size_t a = 0; a < p; ++a) crossproduct[a * p + a] += ridge;
    choleskyDecompose(crossproduct, p);

    double halfLogDet = 0.0;
    for (std::size_t a = 0; a < p; ++a)
      halfLogDet += std::log(crossproduct[a * p + a]);

    // b' M^-1 b through the forward half-solve alone
    solveLowerTriangular(crossproduct, p, projection);
    double quadraticForm = 0.0;
    for (std::size_t a = 0; a < p; ++a)
      quadraticForm += projection[a] * projection[a];

    return 0.5 * static_cast<double>(p) * std::log(ridge) - halfLogDet -
           0.5 * (responseSumOfSquares - quadraticForm) / residualVariance;
  }

  /// (b0, b) | z ~ N(M^-1 b, sigma^2 M^-1): with M = LL' and v = L^-1 b, the
  /// draw is L'^-1 (v + sigma eps) for eps standard normal, drawn in
  /// coordinate order (intercept first). Empty leaves zero the block without
  /// consuming generator draws, like the constant leaf's empty-leaf zero.
  void drawFromPosteriorForNode(ext_rng* rng, const Tree& tree,
                                const double* y, const double* weights,
                                double k, double residualVariance,
                                int32_t nodeIndex, double* out) const {
    std::size_t p = numParams();
    if (tree.at(nodeIndex).numObservations() == 0) {
      for (std::size_t a = 0; a < p; ++a) out[a] = 0.0;
      return;
    }

    double crossproduct[maxStatisticSize];
    double responseSumOfSquares;
    accumulateNodeStatistics(tree, y, weights, nodeIndex, crossproduct, out,
                             &responseSumOfSquares);

    double ridge = (k / scale) * (k / scale) * residualVariance;
    for (std::size_t a = 0; a < p; ++a) crossproduct[a * p + a] += ridge;
    choleskyDecompose(crossproduct, p);

    solveLowerTriangular(crossproduct, p, out);
    double sigma = std::sqrt(residualVariance);
    for (std::size_t a = 0; a < p; ++a)
      out[a] += sigma * ext_rng_simulateStandardNormal(rng);
    solveLowerTriangularTransposed(crossproduct, p, out);
  }

  void drawFromPrior(ext_rng* rng, double k, double* out) const {
    std::size_t p = numParams();
    for (std::size_t a = 0; a < p; ++a)
      out[a] = (scale / k) * ext_rng_simulateStandardNormal(rng);
  }

private:
  static constexpr std::size_t maxStatisticSize =
    (maxNumCovariates + 1) * (maxNumCovariates + 1);

  /// One pass over the node's index segment: crossproduct receives the full
  /// symmetric U'WU (row-major (q+1) x (q+1), leading intercept column),
  /// projection U'Wz, and responseSumOfSquares z'Wz. U'WU is served from the
  /// crossproduct cache when the node's member list matches an entry; the
  /// residual-dependent projection and z'Wz always rescan. A served value is
  /// bitwise the fresh scan's, since the entry was built by the identical
  /// fused loop over the same members, covariates, and weights.
  void accumulateNodeStatistics(const Tree& tree, const double* y,
                                const double* weights, int32_t nodeIndex,
                                double* crossproduct, double* projection,
                                double* responseSumOfSquares) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t p = numParams();
    double row[maxNumCovariates + 1];
    row[0] = 1.0;

    const double* cached = lookupCrossproduct(tree, node, nodeIndex);
    if (cached != nullptr) {
      for (std::size_t a = 0; a < p; ++a) projection[a] = 0.0;
      *responseSumOfSquares = 0.0;
      for (std::size_t m = node.begin; m < node.end; ++m) {
        std::size_t i = tree.indices[m];
        double w = weights == nullptr ? 1.0 : weights[i];
        double z = y[i];
        for (std::size_t j = 0; j < numCovariates_; ++j)
          row[j + 1] = u_[i + j * numObservations_];
        for (std::size_t a = 0; a < p; ++a) {
          double scaled = w * row[a];
          projection[a] += scaled * z;
        }
        *responseSumOfSquares += w * z * z;
      }
      std::memcpy(crossproduct, cached, p * p * sizeof(double));
      return;
    }

    for (std::size_t a = 0; a < p * p; ++a) crossproduct[a] = 0.0;
    for (std::size_t a = 0; a < p; ++a) projection[a] = 0.0;
    *responseSumOfSquares = 0.0;
    for (std::size_t m = node.begin; m < node.end; ++m) {
      std::size_t i = tree.indices[m];
      double w = weights == nullptr ? 1.0 : weights[i];
      double z = y[i];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        row[j + 1] = u_[i + j * numObservations_];
      for (std::size_t a = 0; a < p; ++a) {
        double scaled = w * row[a];
        projection[a] += scaled * z;
        for (std::size_t b = a; b < p; ++b)
          crossproduct[a * p + b] += scaled * row[b];
      }
      *responseSumOfSquares += w * z * z;
    }
    for (std::size_t a = 0; a < p; ++a)
      for (std::size_t b = a + 1; b < p; ++b)
        crossproduct[b * p + a] = crossproduct[a * p + b];
    storeCrossproduct(tree, node, nodeIndex, crossproduct);
  }

  /// In-place lower Cholesky of a symmetric positive-definite matrix; the
  /// ridge guarantees definiteness, so there is no failure path.
  static void choleskyDecompose(double* m, std::size_t p) {
    for (std::size_t j = 0; j < p; ++j) {
      double diagonal = m[j * p + j];
      for (std::size_t a = 0; a < j; ++a)
        diagonal -= m[j * p + a] * m[j * p + a];
      diagonal = std::sqrt(diagonal);
      m[j * p + j] = diagonal;
      for (std::size_t i = j + 1; i < p; ++i) {
        double value = m[i * p + j];
        for (std::size_t a = 0; a < j; ++a)
          value -= m[i * p + a] * m[j * p + a];
        m[i * p + j] = value / diagonal;
      }
    }
  }

  static void solveLowerTriangular(const double* l, std::size_t p, double* x) {
    for (std::size_t i = 0; i < p; ++i) {
      double value = x[i];
      for (std::size_t a = 0; a < i; ++a) value -= l[i * p + a] * x[a];
      x[i] = value / l[i * p + i];
    }
  }

  static void solveLowerTriangularTransposed(const double* l, std::size_t p,
                                             double* x) {
    for (std::size_t i = p; i > 0; --i) {
      double value = x[i - 1];
      for (std::size_t a = i; a < p; ++a) value -= l[a * p + (i - 1)] * x[a];
      x[i - 1] = value / l[(i - 1) * p + (i - 1)];
    }
  }

  /// One leaf's cached U'WU (row-major p x p), tagged with the exact ordered
  /// member list that built it. Coherence: every lookup re-validates by
  /// comparing that list against tree.indices[begin..end], so any structural
  /// move or rejected-move rollback that alters membership fails the compare
  /// and rebuilds - the key is rollback-stable with no per-move hook. U'WU's
  /// other inputs (covariates, weights) are held fixed within a cache
  /// lifetime by clearing wholesale from the two covariate regathers and from
  /// invalidateStatistics (weights); within a sweep the score and draw phases
  /// see identical weights, so a served value is bitwise the fresh scan's.
  struct CachedNodeStatistics {
    std::vector<std::size_t> members;
    double crossproduct[maxStatisticSize];
  };
  struct TreeStatisticsCache {
    const Tree* tree = nullptr;
    std::vector<CachedNodeStatistics> nodes;  // arena-indexed
  };

  /// Leaves below this rescan; their U'WU is cheap and churns fast.
  static constexpr std::size_t minCachedLeafSize = 32;
  /// Byte ceiling over cached member lists, split across chains at initialize;
  /// when spent, further leaves rescan (still correct, just uncached).
  static constexpr std::size_t statisticsCacheTotalBudgetBytes =
    static_cast<std::size_t>(256) << 20;

  static std::size_t statisticsEntryBytes(std::size_t numObs) {
    return numObs * sizeof(std::size_t);
  }

  TreeStatisticsCache& statisticsCacheForTree(const Tree& tree) const {
    for (TreeStatisticsCache& cache : statisticsCaches_)
      if (cache.tree == &tree) return cache;
    statisticsCaches_.emplace_back();
    statisticsCaches_.back().tree = &tree;
    return statisticsCaches_.back();
  }

  void clearStatisticsCache() const {
    statisticsCaches_.clear();
    statisticsCacheUsedBytes_ = 0;
  }

  /// The node's cached U'WU when its ordered member list matches, else null.
  const double* lookupCrossproduct(const Tree& tree, const Node& node,
                                   int32_t nodeIndex) const {
    std::size_t numObs = node.numObservations();
    if (numObs < minCachedLeafSize) return nullptr;
    TreeStatisticsCache& cache = statisticsCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) return nullptr;
    const CachedNodeStatistics& entry = cache.nodes[index];
    if (entry.members.size() == numObs &&
        std::memcmp(entry.members.data(), tree.indices + node.begin,
                    numObs * sizeof(std::size_t)) == 0)
      return entry.crossproduct;
    return nullptr;
  }

  /// Record the freshly scanned U'WU for the node, subject to the byte budget.
  void storeCrossproduct(const Tree& tree, const Node& node, int32_t nodeIndex,
                         const double* crossproduct) const {
    std::size_t numObs = node.numObservations();
    if (numObs < minCachedLeafSize) return;
    std::size_t p = numParams();
    TreeStatisticsCache& cache = statisticsCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) cache.nodes.resize(index + 1);
    CachedNodeStatistics& entry = cache.nodes[index];
    std::size_t oldBytes = statisticsEntryBytes(entry.members.size());
    std::size_t newBytes = statisticsEntryBytes(numObs);
    if (statisticsCacheUsedBytes_ - oldBytes + newBytes >
        statisticsCacheBudget_) {
      if (!entry.members.empty()) {
        statisticsCacheUsedBytes_ -= oldBytes;
        entry.members.clear();
      }
      return;
    }
    statisticsCacheUsedBytes_ += newBytes - oldBytes;
    entry.members.assign(tree.indices + node.begin, tree.indices + node.end);
    std::memcpy(entry.crossproduct, crossproduct, p * p * sizeof(double));
  }

  std::size_t numCovariates_ = 0;
  std::size_t numObservations_ = 0;
  std::size_t numTestObservations_ = 0;
  std::vector<std::size_t> columns_;
  std::vector<double> means_, sds_;
  std::vector<double> u_;      // standardized, column-major n x q
  std::vector<double> uTest_;  // standardized, column-major numTest x q
  mutable std::vector<TreeStatisticsCache> statisticsCaches_;
  mutable std::size_t statisticsCacheUsedBytes_ = 0;
  std::size_t statisticsCacheBudget_ = statisticsCacheTotalBudgetBytes;
};

static_assert(VectorLeafModel<LinearGaussianLeaf>);

/// Gaussian-process leaf: each bottom node fits a smooth function of q
/// designated ordinal predictor columns, standardized internally exactly as
/// the linear leaf's. The leaf function over the standardized covariates u
/// is
///   f ~ GP(0, (scale / k)^2 C),
/// with C a squared-exponential correlation kernel (per-column lengthscales
/// theta plus a small fixed nugget for conditioning), so the prior variance
/// ties to k as the constant and linear leaves do and a constant kernel
/// recovers the constant leaf. The marginal over f is the closed form
/// log N(z; 0, sigma^2 W^-1 + (scale / k)^2 C) - one O(n_leaf^3) Cholesky -
/// so the conjugate MH moves apply unchanged. The drawn function values at
/// the leaf's training rows ARE the parameters and land directly in the
/// per-tree fits; test rows evaluate the conditional mean c(x*)' C^-1 f
/// through per-draw cached weights.
///
/// Leaves larger than maxLeafSize score and draw as constant leaves instead
/// (delegating to ConstantGaussianLeaf's exact math). The proposal's veto
/// guardrail would deadlock: every tree starts as a single root leaf holding
/// all observations, and a birth splitting one over-cap leaf into two
/// over-cap children can never be accepted against a vetoed-vs-doubly-vetoed
/// comparison. The fallback instead makes oversized regions behave exactly
/// like constant-leaf BART - data-informed splits throughout - with the GP
/// refinement switching on once a leaf falls under the cap; the scoring rule
/// is a deterministic function of leaf membership, so the MH comparisons
/// stay coherent.
struct GPGaussianLeaf {
  static constexpr bool hasVectorParams = false;
  static constexpr bool hasFunctionParams = true;
  /// The replay predictor's stack scratch bound; the factory validates.
  static constexpr std::size_t maxNumCovariates = maxFunctionLeafCovariates;
  /// Conditioning jitter on the correlation diagonal; not a modeling knob.
  static constexpr double nugget = 1.0e-6;

  double scale = 1.0;  // nodeScale / sqrt(numTrees)

  std::size_t numCovariates() const { return numCovariates_; }
  const std::vector<std::size_t>& covariateColumns() const { return columns_; }
  const std::vector<double>& covariateMeans() const { return means_; }
  const std::vector<double>& covariateSds() const { return sds_; }
  const std::vector<double>& lengthscales() const { return lengthscales_; }
  std::size_t maxLeafSize() const { return maxLeafSize_; }

  /// Gather and standardize the designated columns and fix the kernel
  /// lengthscales: the supplied per-column values (standardized scale) when
  /// non-null, otherwise the median pairwise-distance heuristic over the
  /// standardized training values.
  void initialize(const ColumnStore& data, const std::size_t* columns,
                  std::size_t numColumns, const double* lengthscales,
                  std::size_t maxLeafSize, std::size_t numChains = 1) {
    numCovariates_ = numColumns;
    columns_.assign(columns, columns + numColumns);
    maxLeafSize_ = maxLeafSize;
    kernelCacheBudget_ =
      kernelCacheTotalBudgetBytes / (numChains == 0 ? 1 : numChains);
    if (lengthscales != nullptr)
      suppliedLengthscales_.assign(lengthscales, lengthscales + numColumns);
    else
      suppliedLengthscales_.clear();
    reinitialize(data);
  }

  /// Re-derive standardization constants, covariates, and (when not
  /// supplied) lengthscales for the stored designation; the whole-data
  /// analogue of the rebuilt cut grid, sharing the linear leaf's calibration
  /// inheritance on views.
  void reinitialize(const ColumnStore& data) {
    clearKernelCaches();
    std::size_t numColumns = numCovariates_;
    numObservations_ = data.numObservations;
    means_.assign(numColumns, 0.0);
    sds_.assign(numColumns, 1.0);
    u_.resize(numObservations_ * numColumns);
    for (std::size_t j = 0; j < numColumns; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double mean, sd;
      if (!data.suppliedStandardization(columns_[j], &mean, &sd))
        standardizationMomentsForColumn(column, numObservations_, &mean, &sd);
      means_[j] = mean;
      sds_[j] = sd;
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - mean) / sds_[j];
    }
    lengthscales_.assign(numColumns, 1.0);
    for (std::size_t j = 0; j < numColumns; ++j)
      lengthscales_[j] = suppliedLengthscales_.empty()
        ? medianPairwiseDistance(u_.data() + j * numObservations_,
                                 numObservations_)
        : suppliedLengthscales_[j];
    rebuildTestCovariates(data);
  }

  /// Regather the training covariates under the existing standardization
  /// constants and lengthscales, for the predictor-mutation surface; the
  /// prior's calibration is sticky across in-place value changes.
  void regatherTrainingCovariates(const ColumnStore& data) {
    clearKernelCaches();
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double* u = u_.data() + j * numObservations_;
      for (std::size_t i = 0; i < numObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Regather the test covariates under the training standardization; called
  /// whenever the store's test data changes.
  void rebuildTestCovariates(const ColumnStore& data) {
    numTestObservations_ = data.numTestObservations;
    uTest_.assign(numTestObservations_ * numCovariates_, 0.0);
    if (numTestObservations_ == 0) return;
    for (std::size_t j = 0; j < numCovariates_; ++j) {
      const double* column = data.rawTestColumn(columns_[j]);
      if (column == nullptr) continue;
      double* u = uTest_.data() + j * numTestObservations_;
      for (std::size_t i = 0; i < numTestObservations_; ++i)
        u[i] = isNA(column[i]) ? 0.0 : (column[i] - means_[j]) / sds_[j];
    }
  }

  /// Reset the per-draw prediction cache before one tree's parameter sweep;
  /// draws fill it node by node and the test-fit loop reads it while the
  /// tree's partitions are unchanged.
  void beginTreeDraw(const Tree& tree) const {
    nodeConstant_.assign(tree.nodes.size(), 0.0);
    nodeAlphaOffset_.assign(tree.nodes.size(), -1);
    alphaBuffer_.clear();
    evictStaleKernelEntries(tree);
  }

  /// score = -0.5 log det V + 0.5 log det(sigma^2 W^-1) - 0.5 z' V^-1 z with
  /// V = (scale / k)^2 C + sigma^2 W^-1, dropping the same per-observation
  /// terms the constant leaf drops; a constant kernel reduces this to the
  /// constant leaf's formula by Sherman-Morrison.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) return 0.0;
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      return fallback.logIntegratedLikelihoodForNode(tree, y, weights, k,
                                                     residualVariance,
                                                     nodeIndex);
    }
    if (weights != nullptr && anyZeroWeight(tree, node, weights))
      return logIntegratedLikelihoodOverPositiveWeights(
        tree, node, y, weights, k, residualVariance, nodeIndex, numObs);

    double s2 = (scale / k) * (scale / k);

    const CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    if (entry != nullptr) {
      cholV_.resize(numObs * numObs);
      for (std::size_t a = 0; a < numObs * numObs; ++a)
        cholV_[a] = entry->kernel[a] * s2;
    } else {
      gatherLeafCovariates(tree, node, numObs);
      buildKernel(numObs, cholV_);
      for (std::size_t a = 0; a < numObs * numObs; ++a) cholV_[a] *= s2;
    }
    double logDetNoise = 0.0;
    for (std::size_t r = 0; r < numObs; ++r) {
      double w =
        weights == nullptr ? 1.0 : weights[tree.indices[node.begin + r]];
      double noise = residualVariance / w;
      logDetNoise += std::log(noise);
      cholV_[r * numObs + r] += noise;
    }
    choleskyDecomposeLeaf(cholV_.data(), numObs);

    double halfLogDetV = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      halfLogDetV += std::log(cholV_[r * numObs + r]);

    vectorScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      vectorScratch_[r] = y[tree.indices[node.begin + r]];
    solveLowerLeaf(cholV_.data(), numObs, vectorScratch_.data());
    double quadraticForm = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      quadraticForm += vectorScratch_[r] * vectorScratch_[r];

    return -halfLogDetV + 0.5 * logDetNoise - 0.5 * quadraticForm;
  }

  /// Draws f | z into fits at the leaf's member observations by Matheron's
  /// rule - f = f0 + s^2 C V^-1 (z - f0 - e0) with f0 ~ N(0, s^2 C) and
  /// e0 ~ N(0, sigma^2 W^-1) - which has exactly the posterior law and needs
  /// only the two factorizations already in hand. Consumes 2 n_leaf standard
  /// normals (n_leaf + n_positive when zero-weight members are present);
  /// empty leaves consume none and cache a zero constant. Caches the
  /// prediction weights alpha = C^-1 f for the node and returns the chi-k
  /// contribution (f' alpha over n_leaf coordinates).
  FunctionLeafDrawStats drawFromPosteriorForNode(
    ext_rng* rng, const Tree& tree, const double* y, const double* weights,
    double k, double residualVariance, int32_t nodeIndex,
    double* fits) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      setConstantCache(nodeIndex, 0.0);
      return FunctionLeafDrawStats();
    }
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      double value = fallback.drawFromPosteriorForNode(rng, tree, k,
                                                       residualVariance,
                                                       nodeIndex);
      for (std::size_t m = node.begin; m < node.end; ++m)
        fits[tree.indices[m]] = value;
      setConstantCache(nodeIndex, value);
      return FunctionLeafDrawStats{value * value, 1.0};
    }
    if (weights != nullptr && anyZeroWeight(tree, node, weights))
      return drawFromPosteriorOverPositiveWeights(
        rng, tree, node, y, weights, k, residualVariance, nodeIndex, fits);

    double s = scale / k, s2 = s * s;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);
    cholV_.resize(numObs * numObs);
    for (std::size_t a = 0; a < numObs * numObs; ++a)
      cholV_[a] = s2 * kernel[a];
    for (std::size_t r = 0; r < numObs; ++r) {
      double w =
        weights == nullptr ? 1.0 : weights[tree.indices[node.begin + r]];
      cholV_[r * numObs + r] += residualVariance / w;
    }
    choleskyDecomposeLeaf(cholV_.data(), numObs);

    // f0 = s L_C eps, drawn first in row order, then e0 row by row
    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
    }
    vectorScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double w =
        weights == nullptr ? 1.0 : weights[tree.indices[node.begin + r]];
      vectorScratch_[r] = y[tree.indices[node.begin + r]] - fScratch_[r] -
                          std::sqrt(residualVariance / w) *
                            ext_rng_simulateStandardNormal(rng);
    }
    solveLowerLeaf(cholV_.data(), numObs, vectorScratch_.data());
    solveLowerTransposedLeaf(cholV_.data(), numObs, vectorScratch_.data());
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a < numObs; ++a)
        value += kernel[r * numObs + a] * vectorScratch_[a];
      fScratch_[r] += s2 * value;
      fits[tree.indices[node.begin + r]] = fScratch_[r];
    }

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Prior draw f = s L_C eps into fits, with the same prediction cache and
  /// chi-k accounting as the posterior draw; over-cap leaves draw the
  /// constant leaf's prior. Consumes n_leaf standard normals (one over cap,
  /// none when empty).
  FunctionLeafDrawStats drawFromPriorForNode(ext_rng* rng, const Tree& tree,
                                             double k, int32_t nodeIndex,
                                             double* fits) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      setConstantCache(nodeIndex, 0.0);
      return FunctionLeafDrawStats();
    }
    if (numObs > maxLeafSize_) {
      ConstantGaussianLeaf fallback{scale};
      double value = fallback.drawFromPrior(rng, k);
      for (std::size_t m = node.begin; m < node.end; ++m)
        fits[tree.indices[m]] = value;
      setConstantCache(nodeIndex, value);
      return FunctionLeafDrawStats{value * value, 1.0};
    }

    double s = scale / k;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);

    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
      fits[tree.indices[node.begin + r]] = fScratch_[r];
    }

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Append one leaf's saved side-channel block (the
  /// computeFunctionBlockOffsets layout) from the live draw cache: the
  /// cached alpha plus the members' plain standardized rows, in member
  /// order; constant-valued nodes (over-cap, empty) append [0, constant].
  /// Valid under the same freshness condition as the test-fit evaluation.
  void appendLeafBlockFromCache(const Tree& tree, int32_t nodeIndex,
                                std::vector<double>& blocks) const {
    std::ptrdiff_t offset =
      nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)];
    if (offset < 0) {
      blocks.push_back(0.0);
      blocks.push_back(nodeConstant_[static_cast<std::size_t>(nodeIndex)]);
      return;
    }
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    blocks.push_back(static_cast<double>(numObs));
    const double* alpha = alphaBuffer_.data() + offset;
    blocks.insert(blocks.end(), alpha, alpha + numObs);
    appendLeafRows(tree, node, numObs, blocks);
  }

  /// The same block recomputed from a tree's persisted fits (the fits ARE
  /// the parameters): over-cap leaves append their constant (fits are
  /// uniform there), gp leaves refactor alpha = C^-1 f against the current
  /// covariates. Serves flatten and live prediction between runs, when no
  /// draw cache is fresh.
  void appendLeafBlock(const Tree& tree, int32_t nodeIndex,
                       const double* fits, std::vector<double>& blocks) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    if (numObs == 0) {
      blocks.push_back(0.0);
      blocks.push_back(0.0);
      return;
    }
    if (numObs > maxLeafSize_) {
      blocks.push_back(0.0);
      blocks.push_back(fits[tree.indices[node.begin]]);
      return;
    }

    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, kernel_);
    cholK_.assign(kernel_.begin(),
                  kernel_.begin() +
                    static_cast<std::ptrdiff_t>(numObs * numObs));
    choleskyDecomposeLeaf(cholK_.data(), numObs);

    blocks.push_back(static_cast<double>(numObs));
    std::size_t alphaStart = blocks.size();
    for (std::size_t r = 0; r < numObs; ++r)
      blocks.push_back(fits[tree.indices[node.begin + r]]);
    solveLowerLeaf(cholK_.data(), numObs, blocks.data() + alphaStart);
    solveLowerTransposedLeaf(cholK_.data(), numObs,
                             blocks.data() + alphaStart);
    appendLeafRows(tree, node, numObs, blocks);
  }

  /// Conditional mean c(x*)' alpha at one test row from the node's cached
  /// draw; constant-valued nodes (over-cap, empty) return their constant.
  /// Valid only between the node's draw and the next structure change - the
  /// alpha weights pair with the member ordering at draw time.
  double fitForTestObservationForNode(const Tree& tree, int32_t nodeIndex,
                                      std::size_t testIndex) const {
    std::ptrdiff_t offset = nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)];
    if (offset < 0) return nodeConstant_[static_cast<std::size_t>(nodeIndex)];
    const Node& node(tree.at(nodeIndex));
    std::size_t numObs = node.numObservations();
    const double* alpha = alphaBuffer_.data() + offset;
    double result = 0.0;
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      double distanceSq = 0.0;
      for (std::size_t j = 0; j < numCovariates_; ++j) {
        double difference = (uTest_[testIndex + j * numTestObservations_] -
                             u_[obs + j * numObservations_]) /
                            lengthscales_[j];
        distanceSq += difference * difference;
      }
      result += std::exp(-0.5 * distanceSq) * alpha[r];
    }
    return result;
  }

private:
  /// Median of pairwise absolute differences on the standardized scale,
  /// deterministically subsampled to at most 512 rows; degenerate columns
  /// fall back to 1.
  static double medianPairwiseDistance(const double* u, std::size_t n) {
    std::size_t stride = 1 + (n - 1) / 512;
    std::vector<double> values;
    for (std::size_t i = 0; i < n; i += stride) values.push_back(u[i]);
    if (values.size() < 2) return 1.0;
    std::vector<double> distances;
    distances.reserve(values.size() * (values.size() - 1) / 2);
    for (std::size_t a = 0; a < values.size(); ++a)
      for (std::size_t b = a + 1; b < values.size(); ++b)
        distances.push_back(std::fabs(values[a] - values[b]));
    std::sort(distances.begin(), distances.end());
    std::size_t count = distances.size();
    double median = count % 2 == 1
      ? distances[count / 2]
      : 0.5 * (distances[count / 2 - 1] + distances[count / 2]);
    return median > 0.0 ? median : 1.0;
  }

  /// Leaf rows' standardized covariates pre-divided by the lengthscales,
  /// row-major numObs x q, so the kernel is exp(-0.5 ||row_r - row_c||^2).
  void gatherLeafCovariates(const Tree& tree, const Node& node,
                            std::size_t numObs) const {
    leafU_.resize(numObs * numCovariates_);
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        leafU_[r * numCovariates_ + j] =
          u_[obs + j * numObservations_] / lengthscales_[j];
    }
  }

  /// Squared-exponential correlation plus the nugget diagonal over the
  /// gathered leaf rows.
  void buildKernel(std::size_t numObs, std::vector<double>& out) const {
    out.resize(numObs * numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      out[r * numObs + r] = 1.0 + nugget;
      for (std::size_t c = r + 1; c < numObs; ++c) {
        double distanceSq = 0.0;
        for (std::size_t j = 0; j < numCovariates_; ++j) {
          double difference = leafU_[r * numCovariates_ + j] -
                              leafU_[c * numCovariates_ + j];
          distanceSq += difference * difference;
        }
        double value = std::exp(-0.5 * distanceSq);
        out[r * numObs + c] = value;
        out[c * numObs + r] = value;
      }
    }
  }

  void setConstantCache(int32_t nodeIndex, double value) const {
    nodeConstant_[static_cast<std::size_t>(nodeIndex)] = value;
    nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)] = -1;
  }

  /// Append the members' plain standardized covariate rows (row-major
  /// numObs x q, member order; lengthscales NOT baked in - replays divide)
  /// to a saved side-channel block.
  void appendLeafRows(const Tree& tree, const Node& node, std::size_t numObs,
                      std::vector<double>& blocks) const {
    for (std::size_t r = 0; r < numObs; ++r) {
      std::size_t obs = tree.indices[node.begin + r];
      for (std::size_t j = 0; j < numCovariates_; ++j)
        blocks.push_back(u_[obs + j * numObservations_]);
    }
  }

  /// alpha = C^-1 f through the correlation Cholesky, appended to the
  /// per-tree buffer; returns the chi-k contribution f' alpha.
  FunctionLeafDrawStats cacheAlphaForNode(int32_t nodeIndex,
                                          std::size_t numObs,
                                          const double* cholK) const {
    std::size_t offset = alphaBuffer_.size();
    alphaBuffer_.insert(alphaBuffer_.end(), fScratch_.begin(),
                        fScratch_.begin() +
                          static_cast<std::ptrdiff_t>(numObs));
    double* alpha = alphaBuffer_.data() + offset;
    solveLowerLeaf(cholK, numObs, alpha);
    solveLowerTransposedLeaf(cholK, numObs, alpha);
    nodeAlphaOffset_[static_cast<std::size_t>(nodeIndex)] =
      static_cast<std::ptrdiff_t>(offset);
    double sumSquared = 0.0;
    for (std::size_t r = 0; r < numObs; ++r)
      sumSquared += fScratch_[r] * alpha[r];
    return FunctionLeafDrawStats{sumSquared,
                                 static_cast<double>(numObs)};
  }

  /// In-place lower Cholesky; the nugget and noise diagonals guarantee
  /// definiteness, so there is no failure path.
  static void choleskyDecomposeLeaf(double* m, std::size_t p) {
    for (std::size_t j = 0; j < p; ++j) {
      double diagonal = m[j * p + j];
      for (std::size_t a = 0; a < j; ++a)
        diagonal -= m[j * p + a] * m[j * p + a];
      diagonal = std::sqrt(diagonal);
      m[j * p + j] = diagonal;
      for (std::size_t i = j + 1; i < p; ++i) {
        double value = m[i * p + j];
        for (std::size_t a = 0; a < j; ++a)
          value -= m[i * p + a] * m[j * p + a];
        m[i * p + j] = value / diagonal;
      }
    }
  }

  static void solveLowerLeaf(const double* l, std::size_t p, double* x) {
    for (std::size_t i = 0; i < p; ++i) {
      double value = x[i];
      for (std::size_t a = 0; a < i; ++a) value -= l[i * p + a] * x[a];
      x[i] = value / l[i * p + i];
    }
  }

  static void solveLowerTransposedLeaf(const double* l, std::size_t p,
                                       double* x) {
    for (std::size_t i = p; i > 0; --i) {
      double value = x[i - 1];
      for (std::size_t a = i; a < p; ++a) value -= l[a * p + (i - 1)] * x[a];
      x[i - 1] = value / l[(i - 1) * p + (i - 1)];
    }
  }

  /// One leaf's cached correlation kernel (nugget included) and its lower
  /// Cholesky factor, tagged with the exact member list that built it. The
  /// kernel depends only on membership order and the fixed lengthscales -
  /// not sigma, k, or the response - so entries survive across sweeps until
  /// an accepted move re-routes observations; every lookup re-validates by
  /// comparing the member list, making a hit bitwise identical to a fresh
  /// build with no invalidation hooks to miss. A vacant slot has an empty
  /// member list.
  struct CachedLeafKernel {
    std::vector<std::size_t> members;
    std::vector<double> kernel;
    std::vector<double> cholK;  // filled lazily by the first draw
  };
  struct TreeKernelCache {
    const Tree* tree = nullptr;
    std::vector<CachedLeafKernel> nodes;  // arena-indexed
  };

  /// Leaves below this recompute; their kernels are cheap and churn fast.
  static constexpr std::size_t minCachedLeafSize = 32;
  /// Hard byte ceiling over all cached kernels, split evenly across chains
  /// at initialize; when spent, further leaves recompute (still correct,
  /// just uncached).
  static constexpr std::size_t kernelCacheTotalBudgetBytes =
    static_cast<std::size_t>(256) << 20;

  static std::size_t kernelEntryBytes(std::size_t numObs) {
    return numObs * sizeof(std::size_t) +
           2 * numObs * numObs * sizeof(double);
  }

  TreeKernelCache& kernelCacheForTree(const Tree& tree) const {
    for (TreeKernelCache& cache : kernelCaches_)
      if (cache.tree == &tree) return cache;
    kernelCaches_.emplace_back();
    kernelCaches_.back().tree = &tree;
    return kernelCaches_.back();
  }

  void clearKernelCaches() const {
    kernelCaches_.clear();
    kernelCacheUsedBytes_ = 0;
  }

  /// The validated cache entry for a leaf, building the kernel into a fresh
  /// entry when none matches and the budget allows; null when the leaf is
  /// too small to bother or the budget refuses, in which case the caller
  /// computes into scratch. The returned pointer is valid until the next
  /// cache mutation, which no caller spans.
  CachedLeafKernel* cachedKernelForNode(const Tree& tree, const Node& node,
                                        int32_t nodeIndex,
                                        std::size_t numObs) const {
    if (numObs < minCachedLeafSize) return nullptr;
    TreeKernelCache& cache = kernelCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(nodeIndex);
    if (index >= cache.nodes.size()) cache.nodes.resize(index + 1);
    CachedLeafKernel& entry = cache.nodes[index];
    if (entry.members.size() == numObs &&
        std::memcmp(entry.members.data(), tree.indices + node.begin,
                    numObs * sizeof(std::size_t)) == 0)
      return &entry;

    std::size_t oldBytes = kernelEntryBytes(entry.members.size());
    std::size_t newBytes = kernelEntryBytes(numObs);
    if (kernelCacheUsedBytes_ - oldBytes + newBytes > kernelCacheBudget_) {
      if (!entry.members.empty()) {
        kernelCacheUsedBytes_ -= oldBytes;
        entry = CachedLeafKernel();
      }
      return nullptr;
    }
    kernelCacheUsedBytes_ += newBytes - oldBytes;
    entry = CachedLeafKernel();  // release any stale capacity
    entry.members.assign(tree.indices + node.begin, tree.indices + node.end);
    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, entry.kernel);
    return &entry;
  }

  /// The leaf's correlation kernel and lower Cholesky factor for a draw,
  /// served from the cache when possible and computed into the scratch
  /// buffers otherwise; either way the values are bitwise those of a fresh
  /// build.
  void kernelAndFactorForNode(const Tree& tree, const Node& node,
                              int32_t nodeIndex, std::size_t numObs,
                              const double** kernel,
                              const double** cholK) const {
    CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    if (entry != nullptr) {
      if (entry->cholK.empty()) {
        entry->cholK.assign(entry->kernel.begin(), entry->kernel.end());
        choleskyDecomposeLeaf(entry->cholK.data(), numObs);
      }
      *kernel = entry->kernel.data();
      *cholK = entry->cholK.data();
      return;
    }
    gatherLeafCovariates(tree, node, numObs);
    buildKernel(numObs, kernel_);
    cholK_.assign(kernel_.begin(),
                  kernel_.begin() +
                    static_cast<std::ptrdiff_t>(numObs * numObs));
    choleskyDecomposeLeaf(cholK_.data(), numObs);
    *kernel = kernel_.data();
    *cholK = cholK_.data();
  }

  /// True when any member observation carries a zero weight. Those rows
  /// have infinite noise variance, so they contribute no likelihood and
  /// the leaf takes the positive-subset paths below; the constant and
  /// linear leaves get the same behavior for free by multiplying by w.
  bool anyZeroWeight(const Tree& tree, const Node& node,
                     const double* weights) const {
    for (std::size_t m = node.begin; m < node.end; ++m)
      if (weights[tree.indices[m]] == 0.0) return true;
    return false;
  }

  /// The marginal over only the positive-weight members - the exact limit
  /// of infinite noise variance on the rest. Still a deterministic
  /// function of leaf membership and the weights, so MH comparisons stay
  /// coherent; a leaf with no positive-weight members scores 0 like an
  /// empty leaf.
  double logIntegratedLikelihoodOverPositiveWeights(
    const Tree& tree, const Node& node, const double* y,
    const double* weights, double k, double residualVariance,
    int32_t nodeIndex, std::size_t numObs) const {
    positiveScratch_.clear();
    for (std::size_t r = 0; r < numObs; ++r)
      if (weights[tree.indices[node.begin + r]] > 0.0)
        positiveScratch_.push_back(r);
    std::size_t numPos = positiveScratch_.size();
    if (numPos == 0) return 0.0;

    double s2 = (scale / k) * (scale / k);

    const CachedLeafKernel* entry =
      cachedKernelForNode(tree, node, nodeIndex, numObs);
    const double* kernel;
    if (entry != nullptr) {
      kernel = entry->kernel.data();
    } else {
      gatherLeafCovariates(tree, node, numObs);
      buildKernel(numObs, kernel_);
      kernel = kernel_.data();
    }
    cholV_.resize(numPos * numPos);
    for (std::size_t a = 0; a < numPos; ++a)
      for (std::size_t b = 0; b < numPos; ++b)
        cholV_[a * numPos + b] =
          kernel[positiveScratch_[a] * numObs + positiveScratch_[b]] * s2;
    double logDetNoise = 0.0;
    for (std::size_t a = 0; a < numPos; ++a) {
      double w = weights[tree.indices[node.begin + positiveScratch_[a]]];
      double noise = residualVariance / w;
      logDetNoise += std::log(noise);
      cholV_[a * numPos + a] += noise;
    }
    choleskyDecomposeLeaf(cholV_.data(), numPos);

    double halfLogDetV = 0.0;
    for (std::size_t a = 0; a < numPos; ++a)
      halfLogDetV += std::log(cholV_[a * numPos + a]);

    vectorScratch_.resize(numPos);
    for (std::size_t a = 0; a < numPos; ++a)
      vectorScratch_[a] = y[tree.indices[node.begin + positiveScratch_[a]]];
    solveLowerLeaf(cholV_.data(), numPos, vectorScratch_.data());
    double quadraticForm = 0.0;
    for (std::size_t a = 0; a < numPos; ++a)
      quadraticForm += vectorScratch_[a] * vectorScratch_[a];

    return -halfLogDetV + 0.5 * logDetNoise - 0.5 * quadraticForm;
  }

  /// Matheron's rule conditioning only on the positive-weight members: f0
  /// covers every member (numObs normals, row order, as always), e0 and
  /// the correction only the positive rows (numPos more normals), so
  /// zero-weight rows draw from the correct conditional law and the
  /// prediction cache stays defined over the full member list. With no
  /// positive members the draw is the prior draw.
  FunctionLeafDrawStats drawFromPosteriorOverPositiveWeights(
    ext_rng* rng, const Tree& tree, const Node& node, const double* y,
    const double* weights, double k, double residualVariance,
    int32_t nodeIndex, double* fits) const {
    std::size_t numObs = node.numObservations();
    positiveScratch_.clear();
    for (std::size_t r = 0; r < numObs; ++r)
      if (weights[tree.indices[node.begin + r]] > 0.0)
        positiveScratch_.push_back(r);
    std::size_t numPos = positiveScratch_.size();

    double s = scale / k, s2 = s * s;

    const double* kernel;
    const double* cholK;
    kernelAndFactorForNode(tree, node, nodeIndex, numObs, &kernel, &cholK);

    epsScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r)
      epsScratch_[r] = ext_rng_simulateStandardNormal(rng);
    fScratch_.resize(numObs);
    for (std::size_t r = 0; r < numObs; ++r) {
      double value = 0.0;
      for (std::size_t a = 0; a <= r; ++a)
        value += cholK[r * numObs + a] * epsScratch_[a];
      fScratch_[r] = s * value;
    }

    if (numPos > 0) {
      cholV_.resize(numPos * numPos);
      for (std::size_t a = 0; a < numPos; ++a)
        for (std::size_t b = 0; b < numPos; ++b)
          cholV_[a * numPos + b] =
            s2 * kernel[positiveScratch_[a] * numObs + positiveScratch_[b]];
      for (std::size_t a = 0; a < numPos; ++a) {
        double w = weights[tree.indices[node.begin + positiveScratch_[a]]];
        cholV_[a * numPos + a] += residualVariance / w;
      }
      choleskyDecomposeLeaf(cholV_.data(), numPos);

      vectorScratch_.resize(numPos);
      for (std::size_t a = 0; a < numPos; ++a) {
        std::size_t r = positiveScratch_[a];
        double w = weights[tree.indices[node.begin + r]];
        vectorScratch_[a] = y[tree.indices[node.begin + r]] - fScratch_[r] -
                            std::sqrt(residualVariance / w) *
                              ext_rng_simulateStandardNormal(rng);
      }
      solveLowerLeaf(cholV_.data(), numPos, vectorScratch_.data());
      solveLowerTransposedLeaf(cholV_.data(), numPos,
                               vectorScratch_.data());
      for (std::size_t r = 0; r < numObs; ++r) {
        double value = 0.0;
        for (std::size_t a = 0; a < numPos; ++a)
          value += kernel[r * numObs + positiveScratch_[a]] *
                   vectorScratch_[a];
        fScratch_[r] += s2 * value;
      }
    }
    for (std::size_t r = 0; r < numObs; ++r)
      fits[tree.indices[node.begin + r]] = fScratch_[r];

    return cacheAlphaForNode(nodeIndex, numObs, cholK);
  }

  /// Drop entries whose node no longer exists, is no longer a bottom node,
  /// or whose membership changed, keeping the budget for live leaves.
  /// Lookups re-validate regardless; this is hygiene, not correctness.
  void evictStaleKernelEntries(const Tree& tree) const {
    for (TreeKernelCache& cache : kernelCaches_) {
      if (cache.tree != &tree) continue;
      for (std::size_t index = 0; index < cache.nodes.size(); ++index) {
        CachedLeafKernel& entry = cache.nodes[index];
        if (entry.members.empty()) continue;
        bool live = index < tree.nodes.size();
        if (live) {
          const Node& node(tree.at(static_cast<int32_t>(index)));
          live = node.isBottom() &&
                 node.numObservations() == entry.members.size() &&
                 std::memcmp(entry.members.data(),
                             tree.indices + node.begin,
                             entry.members.size() * sizeof(std::size_t)) == 0;
        }
        if (!live) {
          kernelCacheUsedBytes_ -= kernelEntryBytes(entry.members.size());
          entry = CachedLeafKernel();
        }
      }
      break;
    }
  }

  std::size_t numCovariates_ = 0;
  std::size_t numObservations_ = 0;
  std::size_t numTestObservations_ = 0;
  std::size_t maxLeafSize_ = 256;
  std::vector<std::size_t> columns_;
  std::vector<double> means_, sds_;
  std::vector<double> lengthscales_, suppliedLengthscales_;
  std::vector<double> u_;      // standardized, column-major n x q
  std::vector<double> uTest_;  // standardized, column-major numTest x q
  // per-call scratch and the per-tree-draw prediction cache; mutable because
  // scoring and drawing are logically const, and safe because a leaf model
  // instance belongs to a single chain (the CGMTreePrior precedent)
  mutable std::vector<double> leafU_;    // row-major numObs x q, / theta
  mutable std::vector<double> kernel_, cholK_, cholV_;
  mutable std::vector<double> epsScratch_, fScratch_, vectorScratch_;
  mutable std::vector<std::size_t> positiveScratch_;  // w > 0 member offsets
  mutable std::vector<double> alphaBuffer_;
  mutable std::vector<double> nodeConstant_;         // arena-indexed
  mutable std::vector<std::ptrdiff_t> nodeAlphaOffset_;  // -1 = constant
  // the cross-sweep kernel cache, keyed by tree address (stable for a
  // chain's lifetime) and node arena index
  mutable std::vector<TreeKernelCache> kernelCaches_;
  mutable std::size_t kernelCacheUsedBytes_ = 0;
  std::size_t kernelCacheBudget_ = kernelCacheTotalBudgetBytes;
};

static_assert(FunctionLeafModel<GPGaussianLeaf>);

/// Chipman-George-McCulloch tree structure prior. Split-variable selection
/// is uniform over available variables when splitProbabilities is null;
/// otherwise proportional to the supplied weights restricted to available
/// variables. DART points this at its Gibbs-updated probability vector.
struct CGMTreePrior {
  double base = 0.95;
  double power = 2.0;
  const double* splitProbabilities = nullptr;  // length numPredictors, or null
  // pooled-column scratch (reachable set, pattern draw); mutable because
  // rule scoring and drawing are logically const, and safe because a prior
  // instance belongs to a single chain
  mutable std::vector<std::uint64_t> reachableScratch_;
  mutable std::vector<std::uint64_t> patternScratch_;
  // per-variable availability from Tree::collectAvailableVariables
  mutable std::vector<std::uint8_t> availableScratch_;

  double growthProbability(const Tree& tree, const ColumnStore& data,
                           int32_t nodeIndex) const {
    if (!tree.hasAnyAvailableVariable(data, nodeIndex)) return 0.0;
    return base / std::pow(1.0 + static_cast<double>(tree.depthOf(nodeIndex)), power);
  }

  double splitVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                     int32_t nodeIndex) const {
    availableScratch_.resize(data.numPredictors);
    std::size_t numAvailable =
      tree.collectAvailableVariables(data, nodeIndex, availableScratch_.data());
    if (splitProbabilities == nullptr)
      return -std::log(static_cast<double>(numAvailable));

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (availableScratch_[j]) totalProbability += splitProbabilities[j];
    return std::log(
      splitProbabilities[tree.at(nodeIndex).rule.variableIndex] / totalProbability);
  }

  /// Categorical rules are uniform over the direction assignments of the
  /// categories reachable at the node that leave neither side empty, both
  /// orientations counted: 2^R - 2 of them. Unreachable categories carry no
  /// probability (their direction bits are pinned to zero by the moves).
  double ruleForVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                       int32_t nodeIndex) const {
    int32_t variableIndex = tree.at(nodeIndex).rule.variableIndex;
    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical) {
      size_t j = static_cast<size_t>(variableIndex);
      size_t numReachable;
      if (data.columnIsPooled(j)) {
        size_t numWords = maskWordsForCount(data.numCuts[j]);
        reachableScratch_.resize(numWords);
        tree.reachableCategoriesWide(data, nodeIndex, variableIndex,
                                     reachableScratch_.data());
        numReachable = maskPopcount(reachableScratch_.data(), numWords);
      } else {
        numReachable = static_cast<size_t>(std::popcount(
          tree.reachableCategories(data, nodeIndex, variableIndex)));
      }
      // past 54 reachable pow(2, R) - 2 is no longer exact; the closed form
      // is a deterministic function of R, so MH ratios stay consistent
      if (numReachable > 54)
        return -(static_cast<double>(numReachable) * std::log(2.0) +
                 std::log1p(-std::pow(2.0, 1.0 -
                                             static_cast<double>(numReachable))));
      return -std::log(std::pow(2.0, static_cast<double>(numReachable)) - 2.0);
    }
    int32_t left, right;
    tree.splitInterval(data, nodeIndex, variableIndex, &left, &right);
    double logNumRules = std::log(static_cast<double>(right - left + 1));
    // a column with missing values widens every rule by its two-way
    // missing direction (the categorical count absorbs the missing
    // category through the reachable mask instead)
    if (data.hasMissing[static_cast<size_t>(variableIndex)])
      logNumRules += std::log(2.0);
    return -logNumRules;
  }

  double treeLogProbability(const Tree& tree, const ColumnStore& data,
                            int32_t nodeIndex = 0) const {
    double growth = growthProbability(tree, data, nodeIndex);
    if (tree.at(nodeIndex).isBottom()) return std::log(1.0 - growth);

    double result = std::log(growth);
    result += splitVariableLogProbability(tree, data, nodeIndex);
    result += ruleForVariableLogProbability(tree, data, nodeIndex);
    result += treeLogProbability(tree, data, tree.at(nodeIndex).leftChild);
    result += treeLogProbability(tree, data, tree.at(nodeIndex).leftChild + 1);
    return result;
  }

  int32_t drawSplitVariable(const Tree& tree, const ColumnStore& data,
                            ext_rng* rng, int32_t nodeIndex) const {
    availableScratch_.resize(data.numPredictors);
    std::size_t numGood =
      tree.collectAvailableVariables(data, nodeIndex, availableScratch_.data());

    if (splitProbabilities == nullptr) {
      std::size_t variableNumber =
        ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, numGood);
      std::size_t count = 0;
      for (std::size_t j = 0; j < data.numPredictors; ++j)
        if (availableScratch_[j]) {
          if (count == variableNumber) return static_cast<int32_t>(j);
          ++count;
        }
      return invalidVariable;
    }

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (availableScratch_[j]) totalProbability += splitProbabilities[j];

    double cutoff = ext_rng_simulateContinuousUniform(rng) * totalProbability;
    double runningProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j) {
      if (availableScratch_[j]) {
        runningProbability += splitProbabilities[j];
        if (runningProbability >= cutoff) return static_cast<int32_t>(j);
      }
    }
    return invalidVariable;  // unreachable with valid probabilities
  }

  /// Build a categorical direction mask from a pattern whose bit k gives
  /// the side of the k-th reachable category (ascending); unreachable
  /// categories stay zero (left), the canonical gauge.
  static std::uint64_t categoryDirectionsForPattern(std::uint64_t reachable,
                                                    std::uint64_t pattern) {
    std::uint64_t directions = 0;
    int bit = 0;
    while (reachable != 0) {
      std::uint32_t category = static_cast<std::uint32_t>(
        std::countr_zero(reachable));
      if ((pattern >> bit) & 1u) directions |= 1ull << category;
      reachable &= reachable - 1;
      ++bit;
    }
    return directions;
  }

  /// Uniform over the direction patterns of numReachable categories that
  /// leave neither side empty, drawn bit by bit with the two all-same
  /// patterns rejected: a single range draw of u * (2^R - 2) has only the
  /// generator's granularity, which for wide masks pins the low pattern
  /// bits to functions of the high ones.
  static std::uint64_t drawCategoryPattern(ext_rng* rng, int numReachable) {
    // 63 categories plus the missing pseudo-category can fill the word
    std::uint64_t allRight =
      numReachable >= 64 ? ~0ull : (1ull << numReachable) - 1ull;
    std::uint64_t pattern;
    do {
      pattern = 0;
      for (int bit = 0; bit < numReachable; ++bit)
        if (ext_rng_simulateBernoulli(rng, 0.5) == 1) pattern |= 1ull << bit;
    } while (pattern == 0 || pattern == allRight);
    return pattern;
  }

  /// The pooled-column extension of drawCategoryPattern: the identical
  /// bit-by-bit scheme, one draw per reachable position ascending, the two
  /// all-same patterns rejected.
  static void drawCategoryPatternWide(ext_rng* rng, size_t numReachable,
                                      std::uint64_t* pattern,
                                      size_t numWords) {
    bool degenerate = true;
    while (degenerate) {
      for (size_t w = 0; w < numWords; ++w) pattern[w] = 0;
      for (size_t bit = 0; bit < numReachable; ++bit)
        if (ext_rng_simulateBernoulli(rng, 0.5) == 1)
          pattern[bit >> 6] |= 1ull << (bit & 63u);
      degenerate = maskIsZero(pattern, numWords) ||
                   maskPopcount(pattern, numWords) == numReachable;
    }
  }

  /// The pooled-column extension of categoryDirectionsForPattern.
  static void categoryDirectionsForPatternWide(const std::uint64_t* reachable,
                                               const std::uint64_t* pattern,
                                               std::uint64_t* directions,
                                               size_t numWords) {
    for (size_t w = 0; w < numWords; ++w) directions[w] = 0;
    size_t bit = 0;
    for (size_t w = 0; w < numWords; ++w) {
      std::uint64_t remaining = reachable[w];
      while (remaining != 0) {
        std::uint32_t category =
          static_cast<std::uint32_t>(std::countr_zero(remaining)) +
          64u * static_cast<std::uint32_t>(w);
        if (((pattern[bit >> 6] >> (bit & 63u)) & 1u) != 0)
          maskSetBit(directions, category);
        remaining &= remaining - 1;
        ++bit;
      }
    }
  }

  /// The tree is non-const because a pooled categorical draw allocates its
  /// mask in the tree's pool; the caller truncates back on rejection.
  Rule drawRuleForVariable(Tree& tree, const ColumnStore& data,
                           ext_rng* rng, int32_t nodeIndex,
                           int32_t variableIndex) const {
    Rule result;
    result.variableIndex = variableIndex;

    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical) {
      size_t j = static_cast<size_t>(variableIndex);
      if (data.columnIsPooled(j)) {
        size_t numWords = maskWordsForCount(data.numCuts[j]);
        reachableScratch_.resize(numWords);
        tree.reachableCategoriesWide(data, nodeIndex, variableIndex,
                                     reachableScratch_.data());
        size_t numReachable =
          maskPopcount(reachableScratch_.data(), numWords);
        patternScratch_.resize(numWords);
        drawCategoryPatternWide(rng, numReachable, patternScratch_.data(),
                                numWords);
        size_t offset = tree.allocateMask(numWords);
        categoryDirectionsForPatternWide(reachableScratch_.data(),
                                         patternScratch_.data(),
                                         tree.mutableMaskWordsFor(offset),
                                         numWords);
        result.setMaskOffset(offset);
        return result;
      }
      std::uint64_t reachable =
        tree.reachableCategories(data, nodeIndex, variableIndex);
      int numReachable = std::popcount(reachable);
      std::uint64_t pattern = drawCategoryPattern(rng, numReachable);
      result.setCategoryDirections(
        categoryDirectionsForPattern(reachable, pattern));
      return result;
    }

    int32_t left, right;
    tree.splitInterval(data, nodeIndex, variableIndex, &left, &right);
    result.setSplitIndex(static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, left, right + 1)));
    // the missing direction is part of the rule, a symmetric coin drawn
    // only when the column can route a missing value - NA-free data spends
    // no extra generator draws
    if (data.hasMissing[static_cast<size_t>(variableIndex)])
      result.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);
    return result;
  }

  Rule drawRuleAndVariable(Tree& tree, const ColumnStore& data,
                           ext_rng* rng, int32_t nodeIndex) const {
    int32_t variableIndex = drawSplitVariable(tree, data, rng, nodeIndex);
    return drawRuleForVariable(tree, data, rng, nodeIndex, variableIndex);
  }
};

/// DART (Linero 2018): Dirichlet prior over split-variable probabilities.
/// s | counts ~ Dirichlet(alpha/p + m_1, ..., alpha/p + m_p), sampled via
/// normalized gammas. The concentration alpha is optionally sampled on a
/// fixed grid of lambda = alpha / (alpha + rho) with lambda ~ Beta(a, b);
/// all lgamma terms of the grid are precomputed, so the per-iteration alpha
/// update is O(gridSize) multiply-adds.
struct DartPrior {
  double alpha = 1.0;
  bool updateAlpha = true;
  double betaA = 0.5, betaB = 1.0, rho = 0.0;  // rho <= 0 means numPredictors
  // Iterations to hold s uniform before updates begin, so the forest is
  // likelihood-informed when counts first enter the Dirichlet; starting cold
  // can lock s away from the signal variables. The BART package convention
  // is half of burn-in.
  std::size_t updateDelay = 0;
  std::vector<double> probabilities;

  void initialize(std::size_t numPredictors) {
    probabilities.assign(numPredictors, 1.0 / static_cast<double>(numPredictors));
    if (rho <= 0.0) rho = static_cast<double>(numPredictors);
    if (updateAlpha) {
      gridAlpha_.resize(gridSize);
      gridConstant_.resize(gridSize);
      gridWeight_.resize(gridSize);
      double p = static_cast<double>(numPredictors);
      for (std::size_t k = 0; k < gridSize; ++k) {
        double lambda = static_cast<double>(k + 1) / static_cast<double>(gridSize + 1);
        double alphaK = lambda * rho / (1.0 - lambda);
        gridAlpha_[k] = alphaK;
        gridConstant_[k] = std::lgamma(alphaK) - p * std::lgamma(alphaK / p) +
                           (betaA - 1.0) * std::log(lambda) +
                           (betaB - 1.0) * std::log(1.0 - lambda);
      }
    }
  }

  void update(ext_rng* rng, const std::uint32_t* splitCounts) {
    if (numUpdatesSkipped_ < updateDelay) {
      ++numUpdatesSkipped_;
      return;
    }

    std::size_t numPredictors = probabilities.size();
    double p = static_cast<double>(numPredictors);

    double total = 0.0;
    for (std::size_t j = 0; j < numPredictors; ++j) {
      double draw = ext_rng_simulateGamma(
        rng, alpha / p + static_cast<double>(splitCounts[j]), 1.0);
      probabilities[j] = draw > 1.0e-300 ? draw : 1.0e-300;
      total += probabilities[j];
    }
    for (std::size_t j = 0; j < numPredictors; ++j) probabilities[j] /= total;

    if (!updateAlpha) return;

    double sumLogProbabilities = 0.0;
    for (std::size_t j = 0; j < numPredictors; ++j)
      sumLogProbabilities += std::log(probabilities[j]);

    double maxLogPosterior = -HUGE_VAL;
    for (std::size_t k = 0; k < gridSize; ++k) {
      gridWeight_[k] = gridConstant_[k] + (gridAlpha_[k] / p) * sumLogProbabilities;
      if (gridWeight_[k] > maxLogPosterior) maxLogPosterior = gridWeight_[k];
    }
    double totalWeight = 0.0;
    for (std::size_t k = 0; k < gridSize; ++k) {
      gridWeight_[k] = std::exp(gridWeight_[k] - maxLogPosterior);
      totalWeight += gridWeight_[k];
    }
    for (std::size_t k = 0; k < gridSize; ++k) gridWeight_[k] /= totalWeight;

    alpha = gridAlpha_[ext_rng_drawFromDiscreteDistribution(rng, gridWeight_.data(),
                                                            gridSize)];
  }

  static constexpr std::size_t gridSize = 1000;

  // State serialization: the delay counter is the only hidden mutable state
  // (alpha and probabilities are public members).
  std::size_t numUpdatesSkipped() const { return numUpdatesSkipped_; }
  void setNumUpdatesSkipped(std::size_t value) { numUpdatesSkipped_ = value; }

private:
  std::vector<double> gridAlpha_, gridConstant_, gridWeight_;
  std::size_t numUpdatesSkipped_ = 0;
};

/// Chi hyperprior on the end-node precision parameter k. Since k ~ chi(nu)
/// gives k^2 ~ Gamma(nu/2, 1/2), the posterior of k^2 is gamma with shape
/// 0.5 (M + nu); a finite prior scale adds 0.5 / scale^2 to the rate.
struct ChiKHyperprior {
  double degreesOfFreedom = 1.5;
  double scale = HUGE_VAL;  // infinite = flat in the rate term

  /// Sentinel ceiling on the sampled k. An improper or weak prior scale leaves
  /// the prior-dominated k Gibbs fixed point with a growth factor above 1, so
  /// k can run away toward a leaf sd that is already statistically zero on the
  /// standardized [-0.5, 0.5] response scale; capping the draw bounds the
  /// excursion. Behavior-neutral outside that runaway regime: healthy draws
  /// sit far below it, so the cap never engages.
  static constexpr double maxDraw = 1.0e6;

  double draw(ext_rng* rng, double sumSquaredParams, double totalNumLeaves,
              double leafScale) const {
    double shape = 0.5 * (totalNumLeaves + degreesOfFreedom);
    // classic form: numTrees * s_sq / nodeScale^2 == s_sq / leafScale^2
    double rate = 0.5 * sumSquaredParams / (leafScale * leafScale);
    if (std::fabs(scale) <= DBL_MAX) rate += 0.5 / (scale * scale);
    double k = std::sqrt(ext_rng_simulateGamma(rng, shape, 1.0 / rate));
    return k > maxDraw ? maxDraw : k;
  }
};

/// Conjugate chi-squared residual variance prior; scale arrives already
/// derived (qchisq(1 - quantile, df) / df, then multiplied by the initial
/// sigma^2 estimate at initialization).
struct ChiSquaredScalePrior {
  double degreesOfFreedom = 3.0;
  double scale = 1.0;

  /// numPositiveWeights is the count of rows with weight > 0 (== numObservations
  /// when unweighted): zero-weight rows drop from the weighted SSR and every
  /// leaf conditional, so they must not inflate the posterior df either.
  double drawSigmaSqFromPosterior(ext_rng* rng, const double* y,
                                  const double* totalFits, const double* weights,
                                  std::size_t numObservations,
                                  std::size_t numPositiveWeights) const {
    double sumOfSquaredResiduals = weights == nullptr
      ? misc_computeSumOfSquaredResiduals(y, numObservations, totalFits)
      : misc_computeWeightedSumOfSquaredResiduals(y, numObservations, weights,
                                                  totalFits);
    double posteriorDegreesOfFreedom =
      degreesOfFreedom + static_cast<double>(numPositiveWeights);
    double posteriorScale = degreesOfFreedom * scale + sumOfSquaredResiduals;
    return posteriorScale /
           ext_rng_simulateChiSquared(rng, posteriorDegreesOfFreedom);
  }
};

/// Response families the sampler can run; gaussian fits the response
/// directly, the binary families fit a latent working response.
enum class ResponseFamily { gaussian, probit, logistic };

/// Response models own the working response the backfitting engine sees and
/// any per-iteration latent refresh; concrete classes own their O(n) loops.
/// Response/offset pointers are borrowed: the caller keeps them alive.
class ResponseModel {
public:
  virtual ~ResponseModel() = default;

  /// The engine's working response; contents may change across iterations
  /// for latent-variable models.
  virtual double* workingResponse() = 0;

  /// Per-observation precisions the engine weights tree updates by: the
  /// fixed user weights for gaussian, the current Polya-Gamma draws for
  /// logistic, null (unit weights) for probit. Refreshed alongside the
  /// working response.
  virtual const double* workingWeights() const = 0;

  /// Whether workingWeights() changes across iterations (the latent
  /// Polya-Gamma refresh). Cross-sweep sufficient-statistic caches drop on
  /// each refresh when true; false families weight by fixed user values.
  virtual bool workingWeightsVaryPerSweep() const { return false; }

  /// Called once per iteration after all trees update; sigma is the chain's
  /// current residual sd on the internal scale (1 for the binary families),
  /// which only the grouped decorator's conjugate update consumes.
  virtual void refreshLatents(ext_rng* rng, const double* totalFits,
                              double sigma) = 0;

  /// Residual standard deviation update; fixed-sigma models return sigma.
  virtual double drawSigma(ext_rng* rng, const double* totalFits,
                           double sigma) = 0;

  // Between-sample mutation (the embedded-Gibbs API). sigmaInOut is the
  // chain's current sigma on the internal scale; models that rescale keep it
  // fixed on the original scale across the change. updateScale re-anchors the
  // transform to the new response, as setOffset does; false locks it.
  virtual void setResponse(const double* y, ext_rng* rng,
                           const double* totalFits, bool updateScale,
                           double* sigmaInOut) = 0;
  virtual void setOffset(const double* offset, bool updateScale,
                         double* sigmaInOut) = 0;

  /// Whole-data replacement, possibly changing the number of observations;
  /// latent states are cold-initialized since the fits they condition on are
  /// stale. Models that rescale keep sigmaInOut fixed on the original scale.
  virtual void setData(const double* y, const double* offset,
                       const double* weights, std::size_t numObservations,
                       double* sigmaInOut) = 0;

  /// Replace the case weights (borrowed). A bare pointer swap: nothing
  /// rescales, and the weighted residuals enter the next iteration's node
  /// statistics and sigma draw. Only gaussian responses carry weights;
  /// elsewhere a no-op (the host rejects earlier).
  virtual void setWeights(const double*) {}

  /// Replace the residual-variance prior: re-anchors to the supplied
  /// original-scale sigma estimate exactly as construction does, so a swap
  /// before any run matches creating with the new prior. A no-op for the
  /// fixed-sigma binary families.
  virtual void setSigmaPrior(double /*sigmaEstimate*/, double /*degreesOfFreedom*/,
                             double /*rawScale*/) {}

  virtual const double* latents() const { return nullptr; }

  /// The current training offset (borrowed), or null. Recorded training
  /// fits add it back, keeping them on the original scale.
  virtual const double* offset() const = 0;

  /// State restoration: overwrite the latent state with previously stored
  /// latents() values and rebuild the working response from them. A no-op
  /// for models without latents.
  virtual void restoreLatents(const double*) {}

  /// State serialization of the response transform: gaussian's offset-
  /// adjusted (min, max), whose evolution under setOffset(updateScale) is
  /// otherwise unrecoverable from the data; scale-free families report
  /// (0, 0) and ignore restoration.
  virtual void getScale(double& min, double& max) const { min = max = 0.0; }
  virtual void restoreScale(double /*min*/, double /*max*/) {}

  virtual double initialSigma() const = 0;

  // Sample de-scaling to the original response scale.
  virtual double fitScale() const = 0;
  virtual double fitShift() const = 0;
  virtual double sigmaScale() const = 0;

  /// Grouped random effects, when the model carries them (the GroupedResponse
  /// decorator): the per-group intercepts and tau, on the internal scale.
  /// Ungrouped models report none.
  virtual const double* groupEffects() const { return nullptr; }
  virtual std::size_t numGroupEffects() const { return 0; }
  virtual double groupTau() const { return 0.0; }
  /// State restoration of the grouped channel; values are internal-scale
  /// exactly as reported, so restores are exact. A no-op when ungrouped.
  virtual void restoreGroupEffects(const double* /*effects*/,
                                   double /*tau*/) {}
};

class GaussianResponse final : public ResponseModel {
public:
  /// sigmaRawScale is qchisq(1 - quantile, df) / df; offset may be null.
  GaussianResponse(const double* y, const double* offset, const double* weights,
                   std::size_t numObservations, double sigmaEstimate,
                   double sigmaDf, double sigmaRawScale)
    : y_(y), offset_(offset), weights_(weights),
      numObservations_(numObservations),
      numPositiveWeights_(countPositiveWeights(weights, numObservations)) {
    yRescaled_.resize(numObservations);
    rescale();

    initialSigma_ = sigmaEstimate / range_;
    sigmaSqPrior_.degreesOfFreedom = sigmaDf;
    sigmaSqPrior_.scale = initialSigma_ * initialSigma_ * sigmaRawScale;
  }

  double* workingResponse() override { return yRescaled_.data(); }
  const double* workingWeights() const override { return weights_; }
  const double* offset() const override { return offset_; }
  void refreshLatents(ext_rng*, const double*, double) override {}

  double drawSigma(ext_rng* rng, const double* totalFits, double) override {
    return std::sqrt(sigmaSqPrior_.drawSigmaSqFromPosterior(
      rng, yRescaled_.data(), totalFits, weights_, numObservations_,
      numPositiveWeights_));
  }

  void setResponse(const double* y, ext_rng*, const double*, bool updateScale,
                   double* sigmaInOut) override {
    if (updateScale) {
      // sigma and the variance prior stay fixed on the original scale
      double sigmaUnscaled = *sigmaInOut * range_;
      double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

      y_ = y;
      rescale();

      sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
      *sigmaInOut = sigmaUnscaled / range_;
    } else {
      // reuse the existing data scale
      y_ = y;
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
      if (offset_ != nullptr)
        misc_subtractVectorsInPlace(offset_, numObservations_, yRescaled_.data());
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
      misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                       1.0 / range_);
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
    }
  }

  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double* sigmaInOut) override {
    // sigma and the variance prior stay fixed on the original scale
    double sigmaUnscaled = *sigmaInOut * range_;
    double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

    y_ = y;
    offset_ = offset;
    weights_ = weights;
    numObservations_ = numObservations;
    numPositiveWeights_ = countPositiveWeights(weights, numObservations);
    yRescaled_.resize(numObservations);
    rescale();

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
    *sigmaInOut = sigmaUnscaled / range_;
  }

  void setWeights(const double* weights) override {
    weights_ = weights;
    numPositiveWeights_ = countPositiveWeights(weights, numObservations_);
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    double sigmaInternal = sigmaEstimate / range_;
    sigmaSqPrior_.degreesOfFreedom = degreesOfFreedom;
    sigmaSqPrior_.scale = sigmaInternal * sigmaInternal * rawScale;
  }

  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    if (updateScale) {
      double sigmaUnscaled = *sigmaInOut * range_;
      double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

      offset_ = offset;
      rescale();

      sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
      *sigmaInOut = sigmaUnscaled / range_;
    } else {
      // reuse the existing data scale
      offset_ = offset;
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
      if (offset_ != nullptr)
        misc_subtractVectorsInPlace(offset_, numObservations_, yRescaled_.data());
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
      misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                       1.0 / range_);
      misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
    }
  }

  double initialSigma() const override { return initialSigma_; }
  double fitScale() const override { return range_; }
  double fitShift() const override { return range_ * 0.5 + min_; }
  double sigmaScale() const override { return range_; }

  void getScale(double& min, double& max) const override {
    min = min_;
    max = max_;
  }

  /// Installs a stored transform, rebuilding the working response under it
  /// and re-anchoring the variance prior on the original scale, exactly as
  /// the mutation methods do; the chain's internal sigma is restored from
  /// the same state afterwards, so no sigmaInOut here.
  void restoreScale(double min, double max) override {
    double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

    min_ = min;
    max_ = max;
    range_ = max_ - min_;
    if (range_ == 0.0) range_ = 1.0;

    std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
    if (offset_ != nullptr)
      misc_subtractVectorsInPlace(offset_, numObservations_,
                                  yRescaled_.data());
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
    misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                     1.0 / range_);
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
  }

private:
  /// Zero-weight rows are documented as ignored; the sigma posterior df counts
  /// only positive-weight rows (all n when unweighted).
  static std::size_t countPositiveWeights(const double* weights,
                                          std::size_t numObservations) {
    if (weights == nullptr) return numObservations;
    std::size_t count = 0;
    for (std::size_t i = 0; i < numObservations; ++i)
      if (weights[i] > 0.0) ++count;
    return count;
  }

  void rescale() {
    if (offset_ != nullptr) {
      misc_subtractVectors(offset_, numObservations_, y_, yRescaled_.data());
    } else {
      std::memcpy(yRescaled_.data(), y_, numObservations_ * sizeof(double));
    }

    min_ = yRescaled_[0];
    max_ = yRescaled_[0];
    for (std::size_t i = 1; i < numObservations_; ++i) {
      if (yRescaled_[i] < min_) min_ = yRescaled_[i];
      if (yRescaled_[i] > max_) max_ = yRescaled_[i];
    }
    range_ = max_ - min_;
    if (range_ == 0.0) range_ = 1.0;

    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -min_);
    misc_scalarMultiplyVectorInPlace(yRescaled_.data(), numObservations_,
                                     1.0 / range_);
    misc_addScalarToVectorInPlace(yRescaled_.data(), numObservations_, -0.5);
  }

  const double* y_;
  const double* offset_;
  const double* weights_;
  std::size_t numObservations_;
  std::size_t numPositiveWeights_;
  std::vector<double> yRescaled_;
  double min_ = 0.0, max_ = 0.0, range_ = 1.0;
  double initialSigma_ = 1.0;
  ChiSquaredScalePrior sigmaSqPrior_;
};

/// Weights are deliberately unsupported: weighted probit has no exact
/// latent-variable form, so latent draws are unweighted by design. The
/// pre-1.0 engine instead scaled the latents by 1 / sqrt(w), which is not a
/// coherent weighted-likelihood model; that scaling was dropped.
class ProbitResponse final : public ResponseModel {
public:
  ProbitResponse(const double* y, const double* offset,
                 std::size_t numObservations)
    : y_(y), offset_(offset), numObservations_(numObservations) {
    latents_.resize(numObservations);
    working_.resize(numObservations);
    // z = 2 y - 1: -1 for failures, 1 for successes.
    misc_setVectorToConstant(latents_.data(), numObservations, -1.0);
    misc_addVectorsInPlaceWithMultiplier(y, numObservations, 2.0,
                                         latents_.data());
    rebuildWorking();
  }

  double* workingResponse() override { return working_.data(); }
  const double* workingWeights() const override { return nullptr; }
  const double* offset() const override { return offset_; }

  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double mean = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double sign = 2.0 * y_[i] - 1.0;
      double z =
        sign * ext_rng_simulateLowerTruncatedNormalScale1(rng, sign * mean, 0.0);
      latents_[i] = !std::isnan(z) ? z : sign * DBL_EPSILON;
    }
    rebuildWorking();
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    refreshLatents(rng, totalFits, 1.0);
  }

  void setOffset(const double* offset, bool, double*) override {
    offset_ = offset;
    rebuildWorking();
  }

  void setData(const double* y, const double* offset, const double*,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    latents_.resize(numObservations);
    working_.resize(numObservations);
    // cold init, z = 2 y - 1
    misc_setVectorToConstant(latents_.data(), numObservations, -1.0);
    misc_addVectorsInPlaceWithMultiplier(y, numObservations, 2.0,
                                         latents_.data());
    rebuildWorking();
  }

  const double* latents() const override { return latents_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(latents_.data(), latents, numObservations_ * sizeof(double));
    rebuildWorking();
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

private:
  void rebuildWorking() {
    std::memcpy(working_.data(), latents_.data(),
                numObservations_ * sizeof(double));
    if (offset_ != nullptr)
      misc_subtractVectorsInPlace(offset_, numObservations_, working_.data());
  }

  const double* y_;
  const double* offset_;
  std::size_t numObservations_;
  std::vector<double> latents_;
  std::vector<double> working_;
};

/// Logistic via Polya-Gamma augmentation (Polson, Scott & Windle 2013):
/// given the fit eta_i = f(x_i) + offset_i, omega_i ~ PG(1, eta_i), and
/// kappa_i / omega_i with kappa_i = y_i - 0.5 is conditionally
/// N(eta_i, 1 / omega_i). The backfitting engine therefore sees working
/// response kappa_i / omega_i - offset_i under per-iteration precision
/// weights omega_i, running the same weighted conjugate updates as a
/// weighted gaussian with sigma fixed at 1. Case weights are unsupported.
/// latents() exposes the omega draws.
class LogisticResponse final : public ResponseModel {
public:
  /// weights are observation counts (positive integers) or null for unit
  /// weights; the host validates integrality. A count-w observation's
  /// Polya-Gamma latent is PG(w, psi), drawn exactly as the sum of w
  /// independent PG(1, psi) variates, so null weights reproduce the
  /// unweighted stream bit for bit.
  LogisticResponse(const double* y, const double* offset,
                   const double* weights, std::size_t numObservations)
    : y_(y), offset_(offset), weights_(weights),
      numObservations_(numObservations) {
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  double* workingResponse() override { return working_.data(); }
  const double* workingWeights() const override { return omega_.data(); }
  bool workingWeightsVaryPerSweep() const override { return true; }
  const double* offset() const override { return offset_; }

  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double offset = offset_ != nullptr ? offset_[i] : 0.0;
      double psi = totalFits[i] + offset;
      long reps = weights_ != nullptr ? std::lround(weights_[i]) : 1L;
      double omega = ext_rng_simulatePolyaGamma(rng, psi);
      for (long c = 1; c < reps; ++c)
        omega += ext_rng_simulatePolyaGamma(rng, psi);
      omega_[i] = omega;
      double weight = weights_ != nullptr ? weights_[i] : 1.0;
      working_[i] = weight * (y_[i] - 0.5) / omega - offset;
    }
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool, double*) override {
    y_ = y;
    refreshLatents(rng, totalFits, 1.0);
  }

  void setOffset(const double* offset, bool, double*) override {
    // omega and kappa / omega stand; only the shift into the working
    // response moves
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double unshifted = working_[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      working_[i] = unshifted - (offset != nullptr ? offset[i] : 0.0);
    }
    offset_ = offset;
  }

  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    weights_ = weights;
    numObservations_ = numObservations;
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  const double* latents() const override { return omega_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(omega_.data(), latents, numObservations_ * sizeof(double));
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double weight = weights_ != nullptr ? weights_[i] : 1.0;
      working_[i] = weight * (y_[i] - 0.5) / omega_[i] -
                    (offset_ != nullptr ? offset_[i] : 0.0);
    }
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

private:
  /// Deterministic cold start, the analogue of probit's z = 2 y - 1: omega
  /// at PG(w, 0)'s mean of w/4, so the working response starts at
  /// 4 (y - 1/2) - offset independent of the weight; real draws replace it
  /// after the first tree sweep.
  void coldStart() {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double weight = weights_ != nullptr ? weights_[i] : 1.0;
      omega_[i] = 0.25 * weight;
      working_[i] =
        4.0 * (y_[i] - 0.5) - (offset_ != nullptr ? offset_[i] : 0.0);
    }
  }

  const double* y_;
  const double* offset_;
  const double* weights_;
  std::size_t numObservations_;
  std::vector<double> omega_;
  std::vector<double> working_;
};

/// Tau priors the in-core grouped sampler supports (rbart.priors); custom
/// R-function priors keep rbart_vi's R loop instead.
enum class TauPriorKind { cauchy, gamma };

/// log p(tau) for the built-in priors, matching rbart.priors evaluated with
/// scale already converted to the engine's internal response scale:
/// dcauchy(tau; 0, scale) or dgamma(tau; shape 2.5, scale), both log.
/// Normalizing constants are kept so component tests compare against R
/// directly; both are scale families, so the internal-scale conversion is
/// exactly a division of the scale parameter.
inline double logTauPrior(TauPriorKind kind, double scale, double tau) {
  if (kind == TauPriorKind::cauchy) {
    double z = tau / scale;
    return -std::log(std::numbers::pi * scale) - std::log1p(z * z);
  }
  return 1.5 * std::log(tau) - tau / scale - 2.5 * std::log(scale) -
         std::lgamma(2.5);
}

/// The R loop's tau log posterior given the group effects:
///   -J log(tau) - 0.5 sum(b^2) / tau^2 + log p(tau)
/// on (0, Inf); anything outside scores -inf.
inline double logTauPosterior(double tau, double numGroups,
                              double sumSquaredEffects, TauPriorKind kind,
                              double priorScale) {
  if (tau <= 0.0 || std::isinf(tau)) return -HUGE_VAL;
  return -numGroups * std::log(tau) -
         0.5 * sumSquaredEffects / (tau * tau) +
         logTauPrior(kind, priorScale, tau);
}

/// One slice-sampling step (Neal 2003: stepping out, then shrinkage) of a
/// log density on (lower, upper); width is the stepping-out unit and the
/// interval endpoints clamp to the boundary, as R's sliceSample does.
/// Stepping out is bounded by Neal's m at 10^4 expansions per side: healthy
/// runs measure under 100 expansions, while an unbounded search chasing a
/// distant mode costs ~mode/width density evaluations (an effective hang
/// once tau has random-walked past ~1e8, reachable through mostly-empty
/// groupings). A capped bracket still contains the current point, so
/// shrinkage samples correctly inside it. Shrinkage terminates in exact
/// arithmetic; its iteration cap turns numeric pathology into keeping the
/// current point.
template <typename LogDensity>
double sliceSampleOnce(ext_rng* rng, const LogDensity& logDensity, double x,
                       double width, double lower, double upper) {
  constexpr int maxStepOut = 10000;
  double logHeight = logDensity(x) - ext_rng_simulateExponential(rng, 1.0);
  double u = ext_rng_simulateContinuousUniform(rng);
  double left = x - u * width;
  double right = left + width;
  int steps = maxStepOut;
  while (steps-- > 0 && left > lower && logDensity(left) > logHeight)
    left -= width;
  if (left < lower) left = lower;
  steps = maxStepOut;
  while (steps-- > 0 && right < upper && logDensity(right) > logHeight)
    right += width;
  if (right > upper) right = upper;

  for (int iteration = 0; iteration < 1000; ++iteration) {
    double proposal =
      left + ext_rng_simulateContinuousUniform(rng) * (right - left);
    if (logDensity(proposal) > logHeight) return proposal;
    if (proposal < x) left = proposal; else right = proposal;
  }
  return x;
}

/// One conjugate draw of grouped random intercepts: with working response
/// z_i ~ N(F_i + b_g(i), sigma^2 / w_i) and prior b_j ~ N(0, tau^2), the
/// groups are independent normals
///   prec_j = (sum_{i in j} w_i) / sigma^2 + 1 / tau^2
///   mean_j = (sum_{i in j} w_i (z_i - F_i)) / (sigma^2 prec_j).
/// Null weights are unit weights; a group with no observations draws from
/// its prior through the same formula. Draws land in effects (also used as
/// the residual accumulator) in group order, one standard normal each;
/// weightScratch supplies numGroups doubles of scratch.
inline void drawGroupEffects(ext_rng* rng, const double* z,
                             const double* weights, const double* totalFits,
                             const std::uint32_t* groupIndex,
                             std::size_t numObservations,
                             std::size_t numGroups, double sigma, double tau,
                             double* effects, double* weightScratch) {
  for (std::size_t j = 0; j < numGroups; ++j) {
    effects[j] = 0.0;
    weightScratch[j] = 0.0;
  }
  for (std::size_t i = 0; i < numObservations; ++i) {
    double w = weights == nullptr ? 1.0 : weights[i];
    std::uint32_t j = groupIndex[i];
    weightScratch[j] += w;
    effects[j] += w * (z[i] - totalFits[i]);
  }
  double sigmaSq = sigma * sigma;
  for (std::size_t j = 0; j < numGroups; ++j) {
    double precision = weightScratch[j] / sigmaSq + 1.0 / (tau * tau);
    double mean = effects[j] / (sigmaSq * precision);
    effects[j] =
      mean + ext_rng_simulateStandardNormal(rng) / std::sqrt(precision);
  }
}

/// Grouped random intercepts as a response-model decorator, running
/// rbart_vi's Gibbs blocks in-core over any base family: the trees condition
/// on workingResponse() = base z - b_g(i), so recorded fits stay f(x)-only
/// and the user offset and gaussian scale anchoring belong to the base model
/// untouched. b and tau live on the base model's internal scale; the tau
/// prior's relative scale converts once at construction (both priors are
/// scale families) and reported draws de-scale by sigmaScale(), like sigma.
/// The weighted conjugate update deliberately replaces the R loop's
/// unweighted group means, which is what makes Polya-Gamma logistic weights
/// compose. Group structure is fixed at creation; the bridge refuses
/// setResponse/setData on grouped samplers (the overrides below delegate
/// faithfully for same-length replacements regardless).
class GroupedResponse final : public ResponseModel {
public:
  /// tauPriorScale is the original-scale relative scale (sd(y) continuous,
  /// 0.5 binary); the prior's 2.5 multiplier and the internal-scale
  /// conversion are applied here. Initialization matches the R loop: tau at
  /// a fifth of the relative scale, b drawn from its prior.
  GroupedResponse(std::unique_ptr<ResponseModel> base,
                  std::size_t numObservations,
                  const std::uint32_t* groupIndices, std::size_t numGroups,
                  TauPriorKind tauPriorKind, double tauPriorScale,
                  std::size_t tauSliceSteps, ext_rng* rng)
    : base_(std::move(base)), numObservations_(numObservations),
      groupIndex_(groupIndices, groupIndices + numObservations),
      priorKind_(tauPriorKind), sliceSteps_(tauSliceSteps) {
    double scale = base_->sigmaScale();
    priorScale_ = 2.5 * tauPriorScale / scale;
    tau_ = (tauPriorScale / 5.0) / scale;
    groupEffects_.resize(numGroups);
    for (std::size_t j = 0; j < numGroups; ++j)
      groupEffects_[j] = tau_ * ext_rng_simulateStandardNormal(rng);
    weightScratch_.resize(numGroups);
    workingResponse_.resize(numObservations);
    fitScratch_.resize(numObservations);
    rebuildWorking();
  }

  double* workingResponse() override { return workingResponse_.data(); }
  const double* workingWeights() const override {
    return base_->workingWeights();
  }
  bool workingWeightsVaryPerSweep() const override {
    return base_->workingWeightsVaryPerSweep();
  }
  const double* offset() const override { return base_->offset(); }

  /// rbart_vi's blocks in their order relative to the tree sweep: the base
  /// refreshes its latents against f + b, then b draws conjugately against
  /// the fresh working response and weights, then tau by slice sampling,
  /// keeping the last of sliceSteps_ steps as the R loop does. sigma is the
  /// previous draw, exactly the value the R loop's b update conditions on.
  void refreshLatents(ext_rng* rng, const double* totalFits,
                      double sigma) override {
    shiftFits(totalFits);
    base_->refreshLatents(rng, fitScratch_.data(), sigma);
    drawGroupEffects(rng, base_->workingResponse(), base_->workingWeights(),
                     totalFits, groupIndex_.data(), numObservations_,
                     groupEffects_.size(), sigma, tau_, groupEffects_.data(),
                     weightScratch_.data());

    double sumSquaredEffects = 0.0;
    for (double effect : groupEffects_)
      sumSquaredEffects += effect * effect;
    double numGroups = static_cast<double>(groupEffects_.size());
    TauPriorKind kind = priorKind_;
    double priorScale = priorScale_;
    auto logDensity = [=](double tau) {
      return logTauPosterior(tau, numGroups, sumSquaredEffects, kind,
                             priorScale);
    };
    for (std::size_t step = 0; step < sliceSteps_; ++step)
      tau_ = sliceSampleOnce(rng, logDensity, tau_, priorScale_, 0.0,
                             HUGE_VAL);

    rebuildWorking();
  }

  /// sigma conditions on f + b through the same shifted scratch.
  double drawSigma(ext_rng* rng, const double* totalFits,
                   double sigma) override {
    shiftFits(totalFits);
    return base_->drawSigma(rng, fitScratch_.data(), sigma);
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   bool updateScale, double* sigmaInOut) override {
    shiftFits(totalFits);
    base_->setResponse(y, rng, fitScratch_.data(), updateScale, sigmaInOut);
    rebuildWorking();
  }

  /// The base handles any scale re-anchoring; b and tau keep their values
  /// against the new transform, so an updateScale offset change shifts what
  /// they mean on the original scale. rbart_vi's in-core path never calls
  /// this.
  void setOffset(const double* offset, bool updateScale,
                 double* sigmaInOut) override {
    base_->setOffset(offset, updateScale, sigmaInOut);
    rebuildWorking();
  }

  /// Group indices are per-observation and fixed at creation, so only a
  /// same-length replacement is coherent; the bridge refuses the call.
  void setData(const double* y, const double* offset, const double* weights,
               std::size_t numObservations, double* sigmaInOut) override {
    base_->setData(y, offset, weights, numObservations, sigmaInOut);
    rebuildWorking();
  }

  void setWeights(const double* weights) override {
    base_->setWeights(weights);
  }

  void setSigmaPrior(double sigmaEstimate, double degreesOfFreedom,
                     double rawScale) override {
    base_->setSigmaPrior(sigmaEstimate, degreesOfFreedom, rawScale);
  }

  const double* latents() const override { return base_->latents(); }

  void restoreLatents(const double* latents) override {
    base_->restoreLatents(latents);
    rebuildWorking();
  }

  void getScale(double& min, double& max) const override {
    base_->getScale(min, max);
  }
  void restoreScale(double min, double max) override {
    base_->restoreScale(min, max);
    rebuildWorking();
  }

  double initialSigma() const override { return base_->initialSigma(); }
  double fitScale() const override { return base_->fitScale(); }
  double fitShift() const override { return base_->fitShift(); }
  double sigmaScale() const override { return base_->sigmaScale(); }

  const double* groupEffects() const override {
    return groupEffects_.data();
  }
  std::size_t numGroupEffects() const override {
    return groupEffects_.size();
  }
  double groupTau() const override { return tau_; }
  void restoreGroupEffects(const double* effects, double tau) override {
    std::memcpy(groupEffects_.data(), effects,
                groupEffects_.size() * sizeof(double));
    tau_ = tau;
    rebuildWorking();
  }

private:
  void shiftFits(const double* totalFits) {
    for (std::size_t i = 0; i < numObservations_; ++i)
      fitScratch_[i] = totalFits[i] + groupEffects_[groupIndex_[i]];
  }

  void rebuildWorking() {
    const double* z = base_->workingResponse();
    for (std::size_t i = 0; i < numObservations_; ++i)
      workingResponse_[i] = z[i] - groupEffects_[groupIndex_[i]];
  }

  std::unique_ptr<ResponseModel> base_;
  std::size_t numObservations_;
  std::vector<std::uint32_t> groupIndex_;  // per observation, 0..J-1
  std::vector<double> groupEffects_;       // b, internal scale
  double tau_ = 1.0;                       // internal scale
  TauPriorKind priorKind_;
  double priorScale_ = 1.0;                // internal scale, 2.5 applied
  std::size_t sliceSteps_ = 1;
  std::vector<double> workingResponse_;    // base z minus b_g(i)
  std::vector<double> fitScratch_;         // f + b for base updates
  std::vector<double> weightScratch_;      // per-group weight sums
};

}  // namespace bartcore

#endif  // BARTCORE_MODEL_HPP
