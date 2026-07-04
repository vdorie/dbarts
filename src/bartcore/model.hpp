#ifndef BARTCORE_MODEL_HPP
#define BARTCORE_MODEL_HPP

#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <cfloat>
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
/// models (one parameter per leaf) keep the classic scalar draw interface,
/// vector models write numParams() doubles per leaf and evaluate fits per
/// observation. Chain code branches on hasVectorParams at compile time, so
/// the scalar paths compile exactly as they always have.
template <typename L>
concept LeafModelCore =
  requires(const L leaf, const Tree& tree, const double* v, double d,
           std::int32_t node) {
    { L::hasVectorParams } -> std::convertible_to<bool>;
    { leaf.logIntegratedLikelihoodForNode(tree, v, v, d, d, node) }
      -> std::same_as<double>;
  };

template <typename L>
concept ScalarLeafModel = LeafModelCore<L> && !L::hasVectorParams &&
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

template <typename L>
concept IntegrableLeafModel = ScalarLeafModel<L> || VectorLeafModel<L>;

/// Constant Gaussian leaf: mu ~ N(0, (scale / k)^2), Gaussian likelihood.
struct ConstantGaussianLeaf {
  static constexpr bool hasVectorParams = false;

  double scale;  // nodeScale / sqrt(numTrees)

  double logIntegratedLikelihood(double k, double residualVariance,
                                 double average, double numEffectiveObservations,
                                 double variance, std::size_t numObservations) const {
    if (numObservations == 0) return 0.0;

    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = numEffectiveObservations / residualVariance;

    double result;
    result  = 0.5 * std::log(priorPrecision / (priorPrecision + posteriorPrecision));
    result -= 0.5 * (variance / residualVariance) *
              static_cast<double>(numObservations - 1);
    result -= 0.5 * ((priorPrecision * average) * (posteriorPrecision * average)) /
              (priorPrecision + posteriorPrecision);
    return result;
  }

  double drawFromPosterior(ext_rng* rng, double k, double average,
                           double numEffectiveObservations,
                           double residualVariance) const {
    double priorPrecision = (k / scale) * (k / scale);
    double posteriorPrecision = numEffectiveObservations / residualVariance;
    double posteriorMean =
      posteriorPrecision * average / (priorPrecision + posteriorPrecision);
    double posteriorSd = 1.0 / std::sqrt(priorPrecision + posteriorPrecision);
    return posteriorMean + posteriorSd * ext_rng_simulateStandardNormal(rng);
  }

  double drawFromPrior(ext_rng* rng, double k) const {
    return (scale / k) * ext_rng_simulateStandardNormal(rng);
  }

  // The node-context interface the engine drives; delegates to the scalar
  // math above against the node's cached statistics.
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* weights, double k,
                                        double residualVariance,
                                        int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    double variance = tree.computeVariance(nodeIndex, y, weights);
    return logIntegratedLikelihood(k, residualVariance, node.average,
                                   node.numEffectiveObservations, variance,
                                   node.numObservations());
  }

  double drawFromPosteriorForNode(ext_rng* rng, const Tree& tree, double k,
                                  double residualVariance,
                                  int32_t nodeIndex) const {
    const Node& node(tree.at(nodeIndex));
    return drawFromPosterior(rng, k, node.average,
                             node.numEffectiveObservations, residualVariance);
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
                  std::size_t numColumns) {
    numCovariates_ = numColumns;
    columns_.assign(columns, columns + numColumns);
    reinitialize(data);
  }

  /// Re-derive the standardization constants and covariates for the stored
  /// designation, for whole-data replacement (possibly a new number of
  /// observations): the analogue of setData rebuilding the cut grid. On a
  /// view the raw values come from the gathered copies and the constants
  /// from the parent's full data - the same calibration inheritance as the
  /// copied cut grid, so every fold runs the prior a full-data fit would.
  void reinitialize(const ColumnStore& data) {
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
  /// projection U'Wz, and responseSumOfSquares z'Wz.
  void accumulateNodeStatistics(const Tree& tree, const double* y,
                                const double* weights, int32_t nodeIndex,
                                double* crossproduct, double* projection,
                                double* responseSumOfSquares) const {
    const Node& node(tree.at(nodeIndex));
    std::size_t p = numParams();
    for (std::size_t a = 0; a < p * p; ++a) crossproduct[a] = 0.0;
    for (std::size_t a = 0; a < p; ++a) projection[a] = 0.0;
    *responseSumOfSquares = 0.0;

    double row[maxNumCovariates + 1];
    row[0] = 1.0;
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

  std::size_t numCovariates_ = 0;
  std::size_t numObservations_ = 0;
  std::size_t numTestObservations_ = 0;
  std::vector<std::size_t> columns_;
  std::vector<double> means_, sds_;
  std::vector<double> u_;      // standardized, column-major n x q
  std::vector<double> uTest_;  // standardized, column-major numTest x q
};

static_assert(VectorLeafModel<LinearGaussianLeaf>);

/// Chipman-George-McCulloch tree structure prior. Split-variable selection
/// is uniform over available variables when splitProbabilities is null;
/// otherwise proportional to the supplied weights restricted to available
/// variables (the classic engine's splitprobs semantics). DART points this
/// at its Gibbs-updated probability vector.
struct CGMTreePrior {
  double base = 0.95;
  double power = 2.0;
  const double* splitProbabilities = nullptr;  // length numPredictors, or null
  // pooled-column scratch (reachable set, pattern draw); mutable because
  // rule scoring and drawing are logically const, and safe because a prior
  // instance belongs to a single chain
  mutable std::vector<std::uint64_t> reachableScratch_;
  mutable std::vector<std::uint64_t> patternScratch_;

  double growthProbability(const Tree& tree, const ColumnStore& data,
                           int32_t nodeIndex) const {
    if (tree.numVariablesAvailable(data, nodeIndex) == 0) return 0.0;
    return base / std::pow(1.0 + static_cast<double>(tree.depthOf(nodeIndex)), power);
  }

  double splitVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                     int32_t nodeIndex) const {
    if (splitProbabilities == nullptr)
      return -std::log(static_cast<double>(tree.numVariablesAvailable(data, nodeIndex)));

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (tree.variableAvailable(data, nodeIndex, static_cast<int32_t>(j)))
        totalProbability += splitProbabilities[j];
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
    if (splitProbabilities == nullptr) {
      std::size_t numGood = tree.numVariablesAvailable(data, nodeIndex);
      std::size_t variableNumber =
        ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, numGood);
      return tree.findIthAvailableVariable(data, nodeIndex, variableNumber);
    }

    double totalProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (tree.variableAvailable(data, nodeIndex, static_cast<int32_t>(j)))
        totalProbability += splitProbabilities[j];

    double cutoff = ext_rng_simulateContinuousUniform(rng) * totalProbability;
    double runningProbability = 0.0;
    for (std::size_t j = 0; j < data.numPredictors; ++j) {
      if (tree.variableAvailable(data, nodeIndex, static_cast<int32_t>(j))) {
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

/// Chi hyperprior on the end-node precision parameter k. The posterior of
/// k^2 is gamma in the sum of squared leaf parameters across the forest;
/// a finite prior scale adds 0.5 / scale^2 to the rate (classic
/// ChiHyperprior::drawFromPosterior).
struct ChiKHyperprior {
  double degreesOfFreedom = 1.25;
  double scale = HUGE_VAL;  // infinite = flat in the rate term

  double draw(ext_rng* rng, double sumSquaredParams, double totalNumLeaves,
              double leafScale) const {
    double shape = 0.5 * (totalNumLeaves + 2.0 * degreesOfFreedom - 1.0);
    // classic form: numTrees * s_sq / nodeScale^2 == s_sq / leafScale^2
    double rate = 0.5 * sumSquaredParams / (leafScale * leafScale);
    if (std::fabs(scale) <= DBL_MAX) rate += 0.5 / (scale * scale);
    return std::sqrt(ext_rng_simulateGamma(rng, shape, 1.0 / rate));
  }
};

/// Conjugate chi-squared residual variance prior; scale arrives already
/// derived (qchisq(1 - quantile, df) / df, then multiplied by the initial
/// sigma^2 estimate at initialization, as in the classic engine).
struct ChiSquaredScalePrior {
  double degreesOfFreedom = 3.0;
  double scale = 1.0;

  double drawSigmaSqFromPosterior(ext_rng* rng, const double* y,
                                  const double* totalFits, const double* weights,
                                  std::size_t numObservations) const {
    double sumOfSquaredResiduals = weights == nullptr
      ? misc_htm_computeSumOfSquaredResiduals(nullptr, 0, y, numObservations,
                                              totalFits)
      : misc_htm_computeWeightedSumOfSquaredResiduals(nullptr, 0, y,
                                                      numObservations, weights,
                                                      totalFits);
    double posteriorDegreesOfFreedom =
      degreesOfFreedom + static_cast<double>(numObservations);
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

  /// Called once per iteration after all trees update.
  virtual void refreshLatents(ext_rng* rng, const double* totalFits) = 0;

  /// Residual standard deviation update; fixed-sigma models return sigma.
  virtual double drawSigma(ext_rng* rng, const double* totalFits,
                           double sigma) = 0;

  // Between-sample mutation (the embedded-Gibbs API). sigmaInOut is the
  // chain's current sigma on the internal scale; models that rescale keep it
  // fixed on the original scale across the change.
  virtual void setResponse(const double* y, ext_rng* rng,
                           const double* totalFits, double* sigmaInOut) = 0;
  virtual void setOffset(const double* offset, bool updateScale,
                         double* sigmaInOut) = 0;

  /// Whole-data replacement, possibly changing the number of observations;
  /// latent states are cold-initialized since the fits they condition on are
  /// stale. Models that rescale keep sigmaInOut fixed on the original scale.
  virtual void setData(const double* y, const double* offset,
                       const double* weights, std::size_t numObservations,
                       double* sigmaInOut) = 0;

  /// Replace the case weights (borrowed). Like the classic engine's, a bare
  /// pointer swap: nothing rescales, and the weighted residuals enter the
  /// next iteration's node statistics and sigma draw. Only gaussian
  /// responses carry weights; elsewhere a no-op (the host rejects earlier).
  virtual void setWeights(const double*) {}

  /// Replace the residual-variance prior (the classic engine's setModel):
  /// re-anchors to the supplied original-scale sigma estimate exactly as
  /// construction does, so a swap before any run matches creating with the
  /// new prior. A no-op for the fixed-sigma binary families.
  virtual void setSigmaPrior(double /*sigmaEstimate*/, double /*degreesOfFreedom*/,
                             double /*rawScale*/) {}

  virtual const double* latents() const { return nullptr; }

  /// The current training offset (borrowed), or null. Recorded training
  /// fits add it back, matching the classic engine's original-scale
  /// convention.
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
  /// The variance prior's internal-scale value, for exact state carry:
  /// re-anchoring the prior through the transform (restoreScale) is a
  /// multiply-divide round trip that can perturb the last bit, so restores
  /// install the recorded value directly. Negative when not applicable.
  virtual double sigmaPriorScaleInternal() const { return -1.0; }
  virtual void restoreSigmaPriorScaleInternal(double /*scale*/) {}

  virtual double initialSigma() const = 0;

  // Sample de-scaling to the original response scale.
  virtual double fitScale() const = 0;
  virtual double fitShift() const = 0;
  virtual double sigmaScale() const = 0;
};

class GaussianResponse final : public ResponseModel {
public:
  /// sigmaRawScale is qchisq(1 - quantile, df) / df; offset may be null.
  GaussianResponse(const double* y, const double* offset, const double* weights,
                   std::size_t numObservations, double sigmaEstimate,
                   double sigmaDf, double sigmaRawScale)
    : y_(y), offset_(offset), weights_(weights),
      numObservations_(numObservations) {
    yRescaled_.resize(numObservations);
    rescale();

    initialSigma_ = sigmaEstimate / range_;
    sigmaSqPrior_.degreesOfFreedom = sigmaDf;
    sigmaSqPrior_.scale = initialSigma_ * initialSigma_ * sigmaRawScale;
  }

  double* workingResponse() override { return yRescaled_.data(); }
  const double* workingWeights() const override { return weights_; }
  const double* offset() const override { return offset_; }
  void refreshLatents(ext_rng*, const double*) override {}

  double drawSigma(ext_rng* rng, const double* totalFits, double) override {
    return std::sqrt(sigmaSqPrior_.drawSigmaSqFromPosterior(
      rng, yRescaled_.data(), totalFits, weights_, numObservations_));
  }

  void setResponse(const double* y, ext_rng*, const double*,
                   double* sigmaInOut) override {
    // sigma and the variance prior stay fixed on the original scale
    double sigmaUnscaled = *sigmaInOut * range_;
    double priorUnscaled = sigmaSqPrior_.scale * range_ * range_;

    y_ = y;
    rescale();

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
    *sigmaInOut = sigmaUnscaled / range_;
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
    yRescaled_.resize(numObservations);
    rescale();

    sigmaSqPrior_.scale = priorUnscaled / (range_ * range_);
    *sigmaInOut = sigmaUnscaled / range_;
  }

  void setWeights(const double* weights) override { weights_ = weights; }

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

  double sigmaPriorScaleInternal() const override {
    return sigmaSqPrior_.scale;
  }
  void restoreSigmaPriorScaleInternal(double scale) override {
    sigmaSqPrior_.scale = scale;
  }

private:
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
  std::vector<double> yRescaled_;
  double min_ = 0.0, max_ = 0.0, range_ = 1.0;
  double initialSigma_ = 1.0;
  ChiSquaredScalePrior sigmaSqPrior_;
};

/// Weights are deliberately unsupported: the reference engine scales the
/// latent draws by 1 / sqrt(w), which does not correspond to a coherent
/// weighted-likelihood model, so the behavior was stripped rather than
/// ported.
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

  void refreshLatents(ext_rng* rng, const double* totalFits) override {
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
                   double*) override {
    y_ = y;
    refreshLatents(rng, totalFits);
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
    // cold init, z = 2 y - 1, as the reference engine's initializeLatents
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
  LogisticResponse(const double* y, const double* offset,
                   std::size_t numObservations)
    : y_(y), offset_(offset), numObservations_(numObservations) {
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  double* workingResponse() override { return working_.data(); }
  const double* workingWeights() const override { return omega_.data(); }
  const double* offset() const override { return offset_; }

  void refreshLatents(ext_rng* rng, const double* totalFits) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double offset = offset_ != nullptr ? offset_[i] : 0.0;
      omega_[i] = ext_rng_simulatePolyaGamma(rng, totalFits[i] + offset);
      working_[i] = (y_[i] - 0.5) / omega_[i] - offset;
    }
  }

  double drawSigma(ext_rng*, const double*, double sigma) override {
    return sigma;
  }

  void setResponse(const double* y, ext_rng* rng, const double* totalFits,
                   double*) override {
    y_ = y;
    refreshLatents(rng, totalFits);
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

  void setData(const double* y, const double* offset, const double*,
               std::size_t numObservations, double*) override {
    y_ = y;
    offset_ = offset;
    numObservations_ = numObservations;
    omega_.resize(numObservations);
    working_.resize(numObservations);
    coldStart();
  }

  const double* latents() const override { return omega_.data(); }

  void restoreLatents(const double* latents) override {
    std::memcpy(omega_.data(), latents, numObservations_ * sizeof(double));
    for (std::size_t i = 0; i < numObservations_; ++i)
      working_[i] = (y_[i] - 0.5) / omega_[i] -
                    (offset_ != nullptr ? offset_[i] : 0.0);
  }

  double initialSigma() const override { return 1.0; }
  double fitScale() const override { return 1.0; }
  double fitShift() const override { return 0.0; }
  double sigmaScale() const override { return 1.0; }

private:
  /// Deterministic cold start, the analogue of probit's z = 2 y - 1:
  /// omega at its PG(1, 0) mean of 1/4, so the working response starts at
  /// 4 kappa = +/- 2; real draws replace it after the first tree sweep.
  void coldStart() {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      omega_[i] = 0.25;
      working_[i] =
        4.0 * (y_[i] - 0.5) - (offset_ != nullptr ? offset_[i] : 0.0);
    }
  }

  const double* y_;
  const double* offset_;
  std::size_t numObservations_;
  std::vector<double> omega_;
  std::vector<double> working_;
};

}  // namespace bartcore

#endif  // BARTCORE_MODEL_HPP
