#ifndef BARTCORE_FACADE_HPP
#define BARTCORE_FACADE_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

#include "sampler.hpp"

namespace bartcore {

/// Type-erased boundary for hosts (the R bridge): all type information is
/// resolved once, in the factory; every later call is one virtual hop into
/// fully typed code.
class SamplerBase {
public:
  virtual ~SamplerBase() = default;

  virtual bool run(std::size_t numBurnIn, std::size_t numSamples,
                   Results& results,
                   const std::function<bool()>& pollInterrupt = {},
                   const SweepCallback& onSweep = {}) = 0;
  virtual void setOffset(const double* offset, bool updateScale) = 0;
  virtual void setResponse(const double* y, bool updateScale) = 0;
  virtual void setWeights(const double* weights) = 0;
  virtual void setSigma(double sigmaOriginalScale) = 0;
  virtual void setTestPredictors(const double* x_test,
                                 std::size_t numTestObservations) = 0;
  virtual void setTestOffset(const double* testOffset) = 0;
  virtual void setData(const double* x, const double* y,
                       std::size_t numObservations, const double* weights,
                       const double* offset, const double* x_test,
                       std::size_t numTestObservations,
                       const double* testOffset) = 0;
  virtual PredictorUpdateResult setPredictor(const double* newX,
                                             bool forceUpdate,
                                             bool updateCutPoints) = 0;
  virtual PredictorUpdateResult updatePredictor(const double* newColumns,
                                                const std::size_t* columns,
                                                std::size_t numColumns,
                                                bool forceUpdate,
                                                bool updateCutPoints) = 0;
  virtual void setCutPoints(const double* const* newCutPoints,
                            const std::uint32_t* numCutPoints,
                            const std::size_t* columns,
                            std::size_t numColumns) = 0;
  virtual bool updatePredictorPerObservation(const double* newColumn,
                                             std::size_t column,
                                             bool* installed) = 0;
  virtual std::unique_ptr<PredictorUpdateSession> beginPredictorUpdate(
    const double* newColumn, std::size_t column) = 0;
  // Saved trees, prediction, and state serialization.
  virtual std::size_t savedTreeCapacity() const = 0;
  virtual std::size_t currentSampleNum() const = 0;
  virtual const std::vector<FlatNode>& savedTree(std::size_t chainNum,
                                                 std::size_t slot,
                                                 std::size_t treeNum) const = 0;
  /// Slopes of one saved tree, parallel to savedTree's pre-order leaves;
  /// meaningful only for vector-parameter leaf models.
  virtual const std::vector<double>& savedTreeSlopes(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum) const = 0;
  /// Flattened mask words of one saved tree; meaningful only when the store
  /// has wide categorical columns.
  virtual const std::vector<std::uint64_t>& savedTreeMasks(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum) const = 0;
  /// slopes, when non-null, receives a vector-parameter tree's slopes
  /// (one block per leaf in pre-order); cleared for scalar leaf models.
  /// masks, when non-null, receives the wide categorical side channel.
  virtual void flattenTree(std::size_t chainNum, std::size_t treeNum,
                           std::vector<FlatNode>& nodes,
                           std::vector<std::uint32_t>& counts,
                           std::vector<double>* slopes = nullptr,
                           std::vector<std::uint64_t>* masks = nullptr) = 0;
  virtual void predict(const double* x_test,
                       std::size_t numTestObservations, double* out) = 0;
  virtual void getState(SamplerStateData& state) = 0;
  virtual bool setState(const SamplerStateData& state) = 0;
  virtual void sampleTreesFromPrior() = 0;
  virtual void sampleNodeParametersFromPrior() = 0;
  virtual void setNumThreads(std::size_t numThreads) = 0;
  virtual void setNumThin(std::size_t numThin) = 0;
  virtual void setVerbose(bool verbose, std::uint32_t printEvery) = 0;
  virtual double fitScale() const = 0;
  virtual void setTreeStorage(bool keepTrees, std::size_t numSamplesToStore) = 0;
  virtual void setModel(const ModelParameters& model) = 0;
  virtual double sumOfSquaredResiduals(std::size_t chainNum) = 0;
  virtual void printTrees(const std::size_t* chainIndices,
                          std::size_t numChainIndices,
                          const std::size_t* sampleIndices,
                          std::size_t numSampleIndices,
                          const std::size_t* treeIndices,
                          std::size_t numTreeIndices) = 0;

  virtual ext_rng* rng() const = 0;
  /// The leaf covariate designation (linear and GP leaves); 0/null for
  /// scalar leaf models.
  virtual std::size_t numLeafCovariates() const = 0;
  virtual const std::size_t* leafCovariateColumns() const = 0;
  /// True for function-valued (GP) leaf models, whose state and reporting
  /// layouts differ from the vector-parameter ones.
  virtual bool usesFunctionLeaves() const = 0;
  virtual ResponseFamily family() const = 0;
  /// Grouped random intercepts: the group count, 0 when ungrouped.
  virtual std::size_t numGroups() const = 0;
  virtual const ColumnStore& data() const = 0;
  virtual const double* latents(std::size_t chainNum) const = 0;
  virtual double sigma(std::size_t chainNum) const = 0;
  virtual bool kIsSampled() const = 0;
  virtual bool usesDart() const = 0;
  virtual std::size_t numChains() const = 0;
  virtual std::size_t numThreads() const = 0;
  /// Forest count: 1 for every non-BCF sampler, 2 for BCF (prognostic +
  /// treatment). The bridge refuses state serialization above 1 (step 4).
  virtual std::size_t numForests() const = 0;
  /// BCF surface (docs/design/bcf.md); no-op/false off BCF. out receives
  /// {a, b0, b1}; forestTotalFits writes numObservations internal-scale fits.
  virtual void setTreatment(const double* z) = 0;
  virtual bool bcfGlue(std::size_t chainNum, double* out) const = 0;
  virtual void forestTotalFits(std::size_t chainNum, std::size_t forestIndex,
                               double* out) const = 0;
  virtual std::size_t numTrees() const = 0;
  virtual std::size_t numObservations() const = 0;
  virtual std::size_t numPredictors() const = 0;
  virtual std::size_t numTestObservations() const = 0;
};

template <IntegrableLeafModel L>
class SamplerFacade final : public SamplerBase {
public:
  template <typename... Args>
  explicit SamplerFacade(Args&&... args) : impl_(std::forward<Args>(args)...) {}

  bool run(std::size_t numBurnIn, std::size_t numSamples, Results& results,
           const std::function<bool()>& pollInterrupt = {},
           const SweepCallback& onSweep = {}) override {
    return impl_.run(numBurnIn, numSamples, results, pollInterrupt, onSweep);
  }
  void setOffset(const double* offset, bool updateScale) override {
    impl_.setOffset(offset, updateScale);
  }
  void setResponse(const double* y, bool updateScale) override {
    impl_.setResponse(y, updateScale);
  }
  void setWeights(const double* weights) override {
    impl_.setWeights(weights);
  }
  void setSigma(double sigmaOriginalScale) override {
    impl_.setSigma(sigmaOriginalScale);
  }
  void setTestPredictors(const double* x_test,
                         std::size_t numTestObservations) override {
    impl_.setTestPredictors(x_test, numTestObservations);
  }
  void setTestOffset(const double* testOffset) override {
    impl_.setTestOffset(testOffset);
  }
  void setData(const double* x, const double* y, std::size_t numObservations,
               const double* weights, const double* offset,
               const double* x_test, std::size_t numTestObservations,
               const double* testOffset) override {
    impl_.setData(x, y, numObservations, weights, offset, x_test,
                  numTestObservations, testOffset);
  }
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) override {
    return impl_.setPredictor(newX, forceUpdate, updateCutPoints);
  }
  PredictorUpdateResult updatePredictor(const double* newColumns,
                                        const std::size_t* columns,
                                        std::size_t numColumns,
                                        bool forceUpdate,
                                        bool updateCutPoints) override {
    return impl_.updatePredictor(newColumns, columns, numColumns, forceUpdate,
                                 updateCutPoints);
  }
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints,
                    const std::size_t* columns,
                    std::size_t numColumns) override {
    impl_.setCutPoints(newCutPoints, numCutPoints, columns, numColumns);
  }
  bool updatePredictorPerObservation(const double* newColumn,
                                     std::size_t column,
                                     bool* installed) override {
    return impl_.updatePredictorPerObservation(newColumn, column, installed);
  }
  std::unique_ptr<PredictorUpdateSession> beginPredictorUpdate(
    const double* newColumn, std::size_t column) override {
    return impl_.beginPredictorUpdate(newColumn, column);
  }
  std::size_t savedTreeCapacity() const override {
    return impl_.savedTreeCapacity();
  }
  std::size_t currentSampleNum() const override {
    return impl_.currentSampleNum();
  }
  const std::vector<FlatNode>& savedTree(std::size_t chainNum,
                                         std::size_t slot,
                                         std::size_t treeNum) const override {
    return impl_.savedTree(chainNum, slot, treeNum);
  }
  const std::vector<double>& savedTreeSlopes(
    std::size_t chainNum, std::size_t slot,
    std::size_t treeNum) const override {
    return impl_.savedTreeSlopes(chainNum, slot, treeNum);
  }
  const std::vector<std::uint64_t>& savedTreeMasks(
    std::size_t chainNum, std::size_t slot,
    std::size_t treeNum) const override {
    return impl_.savedTreeMasks(chainNum, slot, treeNum);
  }
  void flattenTree(std::size_t chainNum, std::size_t treeNum,
                   std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes,
                   std::vector<std::uint64_t>* masks) override {
    impl_.flattenTree(chainNum, treeNum, nodes, counts, slopes, masks);
  }
  void predict(const double* x_test, std::size_t numTestObservations,
               double* out) override {
    impl_.predict(x_test, numTestObservations, out);
  }
  void getState(SamplerStateData& state) override { impl_.getState(state); }
  bool setState(const SamplerStateData& state) override {
    return impl_.setState(state);
  }
  void sampleTreesFromPrior() override { impl_.sampleTreesFromPrior(); }
  void sampleNodeParametersFromPrior() override {
    impl_.sampleNodeParametersFromPrior();
  }
  void setNumThreads(std::size_t numThreads) override {
    impl_.setNumThreads(numThreads);
  }
  void setNumThin(std::size_t numThin) override { impl_.setNumThin(numThin); }
  void setVerbose(bool verbose, std::uint32_t printEvery) override {
    impl_.setVerbose(verbose, printEvery);
  }
  double fitScale() const override { return impl_.fitScale(); }
  void setTreeStorage(bool keepTrees, std::size_t numSamplesToStore) override {
    impl_.setTreeStorage(keepTrees, numSamplesToStore);
  }
  void setModel(const ModelParameters& model) override {
    impl_.setModel(model);
  }
  double sumOfSquaredResiduals(std::size_t chainNum) override {
    return impl_.sumOfSquaredResiduals(chainNum);
  }
  void printTrees(const std::size_t* chainIndices,
                  std::size_t numChainIndices,
                  const std::size_t* sampleIndices,
                  std::size_t numSampleIndices, const std::size_t* treeIndices,
                  std::size_t numTreeIndices) override {
    impl_.printTrees(chainIndices, numChainIndices, sampleIndices,
                     numSampleIndices, treeIndices, numTreeIndices);
  }

  ext_rng* rng() const override { return impl_.rng(); }
  std::size_t numLeafCovariates() const override {
    return impl_.numLeafCovariates();
  }
  const std::size_t* leafCovariateColumns() const override {
    return impl_.leafCovariateColumns();
  }
  bool usesFunctionLeaves() const override { return L::hasFunctionParams; }
  ResponseFamily family() const override { return impl_.family(); }
  std::size_t numGroups() const override { return impl_.numGroups(); }
  const ColumnStore& data() const override { return impl_.data(); }
  const double* latents(std::size_t chainNum) const override {
    return impl_.latents(chainNum);
  }
  double sigma(std::size_t chainNum) const override {
    return impl_.sigma(chainNum);
  }
  bool kIsSampled() const override { return impl_.kIsSampled(); }
  bool usesDart() const override { return impl_.usesDart(); }
  std::size_t numChains() const override { return impl_.numChains(); }
  std::size_t numThreads() const override { return impl_.numThreads(); }
  std::size_t numForests() const override { return impl_.numForests(); }
  void setTreatment(const double* z) override { impl_.setTreatment(z); }
  bool bcfGlue(std::size_t chainNum, double* out) const override {
    return impl_.bcfGlue(chainNum, out);
  }
  void forestTotalFits(std::size_t chainNum, std::size_t forestIndex,
                       double* out) const override {
    impl_.forestTotalFits(chainNum, forestIndex, out);
  }
  std::size_t numTrees() const override { return impl_.numTrees(); }
  std::size_t numObservations() const override { return impl_.numObservations(); }
  std::size_t numPredictors() const override { return impl_.numPredictors(); }
  std::size_t numTestObservations() const override {
    return impl_.numTestObservations();
  }

  Sampler<L>& impl() { return impl_; }

private:
  Sampler<L> impl_;
};

/// One sequential sweep over samplers sharing an index-aligned predictor
/// column: each observation is installed in every sampler or in none, so the
/// fits never diverge. installed receives one flag per observation. Returns
/// the conjunction of the finalize() validities, true by construction of the
/// guard.
inline bool updatePredictorPerObservationJointly(
  SamplerBase* const* samplers, std::size_t numSamplers,
  const double* newColumn, const std::size_t* columns, bool* installed) {
  if (numSamplers == 0) return true;

  std::size_t numObservations = samplers[0]->numObservations();

  std::vector<std::unique_ptr<PredictorUpdateSession>> sessions;
  sessions.reserve(numSamplers);
  for (std::size_t k = 0; k < numSamplers; ++k)
    sessions.push_back(samplers[k]->beginPredictorUpdate(newColumn, columns[k]));

  std::vector<std::size_t> scanOrder(numObservations);
  ext_rng_drawPermutation(samplers[0]->rng(), scanOrder.data(), numObservations);

  for (std::size_t i = 0; i < numObservations; ++i) {
    std::size_t j = scanOrder[i];
    bool valid = true;
    for (std::size_t k = 0; k < numSamplers && valid; ++k)
      valid = sessions[k]->observationWouldRemainValid(j);
    installed[j] = valid;
    if (valid)
      for (std::size_t k = 0; k < numSamplers; ++k)
        sessions[k]->commitObservation(j);
  }

  bool allValid = true;
  for (std::size_t k = 0; k < numSamplers; ++k)
    allValid &= sessions[k]->finalize();
  return allValid;
}

/// The constant-leaf instantiation over the response families. rngs supplies
/// one generator per chain (options.numChains of them).
inline std::unique_ptr<SamplerBase> createConstantLeafSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  ResponseFamily family, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* const* rngs) {
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, family,
    sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
}

/// Shared designation check for the leaf-covariate dispatchers: every
/// designated column must be in range, non-categorical, and able to serve
/// contiguous raw values - the callers differ only in how they answer those
/// last two questions.
template <typename IsCategorical, typename ServesRawValues>
inline bool leafCovariateDesignationIsValid(const SamplerOptions& options,
                                            std::size_t numPredictors,
                                            IsCategorical isCategorical,
                                            ServesRawValues servesRawValues) {
  if (options.leafCovariateColumns == nullptr ||
      options.numLeafCovariates > LinearGaussianLeaf::maxNumCovariates)
    return false;
  for (std::size_t k = 0; k < options.numLeafCovariates; ++k) {
    std::size_t j = options.leafCovariateColumns[k];
    if (j >= numPredictors || isCategorical(j) || !servesRawValues(j))
      return false;
  }
  return true;
}

/// Dispatch on the leaf model: designated leaf covariates select the
/// linear-leaf instantiation - or the GP one under options.gpLeaves -
/// anything else the classic constant leaf. Returns null on an invalid
/// designation - more than maxNumCovariates columns (both leaf models share
/// the bound), a column out of range, or a categorical column (category
/// codes are unordered; interact through splits instead) - which the host
/// turns into its own error.
inline std::unique_ptr<SamplerBase> createSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  ResponseFamily family, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* const* rngs) {
  if (options.numLeafCovariates == 0)
    return createConstantLeafSampler(x, y, numObservations, numPredictors,
                                     weights, offset, family, sigmaEstimate,
                                     sigmaDf, sigmaRawScale, options, rngs);

  // CSC-backed columns hold no contiguous raw values for a leaf model:
  // pure-CSC stores are refused outright, mixed builds per designated column
  if (options.cscColumnPointers != nullptr && options.columnSources == nullptr)
    return nullptr;
  if (!leafCovariateDesignationIsValid(
        options, numPredictors,
        [&](std::size_t j) {
          return options.columnTypes != nullptr &&
                 options.columnTypes[j] == ColumnType::categorical;
        },
        [&](std::size_t j) {
          return options.columnSources == nullptr ||
                 options.columnSources[j] >= 0;
        }))
    return nullptr;
  if (options.gpLeaves)
    return std::make_unique<SamplerFacade<GPGaussianLeaf>>(
      x, y, numObservations, numPredictors, weights, offset, family,
      sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
  return std::make_unique<SamplerFacade<LinearGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, family,
    sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
}

/// As createSampler, but over a pre-built store - typically a row-subset
/// view (ColumnStore::buildFromParent) sharing a parent's cut grid.
/// y/weights/offset match the store's observations and stay alive for the
/// sampler's lifetime; the raw-x mutation surface must not be called on the
/// result when the store holds no raw values. Designated leaf covariates
/// select the linear-leaf instantiation exactly as createSampler does, with
/// one extra requirement: the store must serve raw values for every
/// designated column (a view gathers them when built with the designation).
/// Returns null on an invalid designation.
inline std::unique_ptr<SamplerBase> createSamplerOverStore(
  ColumnStore&& store, const double* y, const double* weights,
  const double* offset, ResponseFamily family, double sigmaEstimate,
  double sigmaDf, double sigmaRawScale, const SamplerOptions& options,
  ext_rng* const* rngs) {
  if (options.numLeafCovariates == 0)
    return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
      std::move(store), y, weights, offset, family, sigmaEstimate, sigmaDf,
      sigmaRawScale, options, rngs);

  if (!leafCovariateDesignationIsValid(
        options, store.numPredictors,
        [&](std::size_t j) { return store.types[j] == ColumnType::categorical; },
        [&](std::size_t j) { return store.rawColumn(j) != nullptr; }))
    return nullptr;
  if (options.gpLeaves)
    return std::make_unique<SamplerFacade<GPGaussianLeaf>>(
      std::move(store), y, weights, offset, family, sigmaEstimate, sigmaDf,
      sigmaRawScale, options, rngs);
  return std::make_unique<SamplerFacade<LinearGaussianLeaf>>(
    std::move(store), y, weights, offset, family, sigmaEstimate, sigmaDf,
    sigmaRawScale, options, rngs);
}

/// A BCF two-forest sampler (docs/design/bcf.md): constant-leaf and gaussian
/// only, so the single instantiation. rngs supplies one generator per chain.
inline std::unique_ptr<SamplerBase> createBCFSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  double sigmaEstimate, double sigmaDf, double sigmaRawScale,
  const SamplerOptions& options, const BCFSpec& spec, ext_rng* const* rngs) {
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, sigmaEstimate,
    sigmaDf, sigmaRawScale, options, spec, rngs);
}

}  // namespace bartcore

#endif  // BARTCORE_FACADE_HPP
