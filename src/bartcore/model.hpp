#ifndef BARTCORE_MODEL_HPP
#define BARTCORE_MODEL_HPP

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

/// Leaf models the conjugate engine can run must be integrable: they expose a
/// closed-form log marginal over their parameter. The suffstat shape is fixed
/// to (average, effective count, variance) for the constant leaf; the concept
/// widens when non-constant leaves land.
template <typename L>
concept IntegrableLeafModel = requires(const L leaf, ext_rng* rng, double d,
                                       std::size_t n) {
  { leaf.logIntegratedLikelihood(d, d, d, d, d, n) } -> std::same_as<double>;
  { leaf.drawFromPosterior(rng, d, d, d, d) } -> std::same_as<double>;
};

/// Constant Gaussian leaf: mu ~ N(0, (scale / k)^2), Gaussian likelihood.
struct ConstantGaussianLeaf {
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
};

static_assert(IntegrableLeafModel<ConstantGaussianLeaf>);

/// Chipman-George-McCulloch tree structure prior. Split-variable selection
/// is uniform over available variables when splitProbabilities is null;
/// otherwise proportional to the supplied weights restricted to available
/// variables (the classic engine's splitprobs semantics). DART points this
/// at its Gibbs-updated probability vector.
struct CGMTreePrior {
  double base = 0.95;
  double power = 2.0;
  const double* splitProbabilities = nullptr;  // length numPredictors, or null

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

  double ruleForVariableLogProbability(const Tree& tree, const ColumnStore& data,
                                       int32_t nodeIndex) const {
    int32_t left, right;
    tree.splitInterval(data, nodeIndex, tree.at(nodeIndex).rule.variableIndex,
                       &left, &right);
    return -std::log(static_cast<double>(right - left + 1));
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

  Rule drawRuleForVariable(const Tree& tree, const ColumnStore& data,
                           ext_rng* rng, int32_t nodeIndex,
                           int32_t variableIndex) const {
    int32_t left, right;
    tree.splitInterval(data, nodeIndex, variableIndex, &left, &right);
    Rule result;
    result.variableIndex = variableIndex;
    result.splitIndex = static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, left, right + 1));
    return result;
  }

  Rule drawRuleAndVariable(const Tree& tree, const ColumnStore& data,
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

private:
  std::vector<double> gridAlpha_, gridConstant_, gridWeight_;
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

/// Response models own the working response the backfitting engine sees and
/// any per-iteration latent refresh; concrete classes own their O(n) loops.
/// Response/offset pointers are borrowed: the caller keeps them alive.
class ResponseModel {
public:
  virtual ~ResponseModel() = default;

  /// The engine's working response; contents may change across iterations
  /// for latent-variable models.
  virtual double* workingResponse() = 0;

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

  virtual const double* latents() const { return nullptr; }

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

class ProbitResponse final : public ResponseModel {
public:
  ProbitResponse(const double* y, const double* offset, const double* weights,
                 std::size_t numObservations)
    : y_(y), offset_(offset), weights_(weights),
      numObservations_(numObservations) {
    latents_.resize(numObservations);
    working_.resize(numObservations);
    // z = 2 y - 1: -1 for failures, 1 for successes.
    misc_setVectorToConstant(latents_.data(), numObservations, -1.0);
    misc_addVectorsInPlaceWithMultiplier(y, numObservations, 2.0,
                                         latents_.data());
    rebuildWorking();
  }

  double* workingResponse() override { return working_.data(); }

  void refreshLatents(ext_rng* rng, const double* totalFits) override {
    for (std::size_t i = 0; i < numObservations_; ++i) {
      double mean = totalFits[i] + (offset_ != nullptr ? offset_[i] : 0.0);
      double sign = 2.0 * y_[i] - 1.0;
      double z = weights_ == nullptr
        ? sign * ext_rng_simulateLowerTruncatedNormalScale1(rng, sign * mean, 0.0)
        : sign * ext_rng_simulateLowerTruncatedNormal(
                   rng, sign * mean, 1.0 / std::sqrt(weights_[i]), 0.0);
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

  const double* latents() const override { return latents_.data(); }

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
  const double* weights_;
  std::size_t numObservations_;
  std::vector<double> latents_;
  std::vector<double> working_;
};

}  // namespace bartcore

#endif  // BARTCORE_MODEL_HPP
