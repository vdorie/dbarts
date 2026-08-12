#ifndef BARTCORE_FACADE_HPP
#define BARTCORE_FACADE_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "sampler.hpp"

namespace bartcore {

/// Every nullary count and capability of a sampler in one by-value POD, so a
/// new one costs a field and a fill line rather than a declaration, a forward,
/// and an override. Read through SamplerBase::shape(), which fills it on
/// demand: a cached copy would go stale behind any mutator.
///
/// Live state (the rng, the data store, sigma, the latents, the current sample
/// number, the fit scale) stays on its own accessor; only quantities fixed
/// between mutations belong here.
struct SamplerShape {
  std::size_t numObservations;
  std::size_t numPredictors;
  std::size_t numTestObservations;
  std::size_t numChains;
  std::size_t numThreads;
  std::size_t numTrees;
  /// Forest count: 1 for single-forest samplers, 2 for BCF (prognostic +
  /// treatment), K for a K-category multinomial.
  std::size_t numForests;
  /// Grouped random intercepts: the group count, 0 when ungrouped.
  std::size_t numGroups;
  /// The leaf covariate designation (linear and GP leaves); 0/null for
  /// scalar leaf models. The columns are borrowed from the leaf model and
  /// stay valid for the sampler's lifetime.
  std::size_t numLeafCovariates;
  const std::size_t* leafCovariateColumns;
  /// Per-observation channels the recorded fits carry: 1 for every additive
  /// model (BCF included), more for a multi-location combiner. The run bridge
  /// reads it to size trainingFits/testFits; internal, invisible to dbarts.h.
  std::size_t numReportedLocations;
  /// Per-sample forests the recorded variable-count channel carries: 1 for
  /// every additive model (BCF included), numCategories for multinomial. The
  /// run bridge reads it to size the varcount array; internal, invisible to
  /// dbarts.h.
  std::size_t numVariableCountForests;
  /// Per-sample cutpoints the recorded cutpoint channel carries: 0 for every
  /// family but ordinal (K-1). The run bridge reads it to size and name the
  /// cutpoints channel, present only when nonzero; internal, invisible to
  /// dbarts.h.
  std::size_t numCutpoints;
  /// Saved samples the tree store holds, 0 when keepTrees is off.
  std::size_t savedTreeCapacity;
  ResponseFamily family;
  /// Whether a heteroscedastic variance forest is present; the s^2(x)
  /// channels of run and predictVariance gate on it.
  bool hasVarianceForest;
  /// The residual prior the variance forest's scale leaf was calibrated from
  /// at creation, original scale; all zero when homoscedastic.
  ResidualPrior varianceLeafPrior;
  /// True for function-valued (GP) leaf models, whose state and reporting
  /// layouts differ from the vector-parameter ones.
  bool usesFunctionLeaves;
  bool kIsSampled;
  bool usesDart;
  /// Whether the forest coupling permits the response-side conduit -
  /// setResponse and setOffset at updateScale = false, and setWeights; false
  /// off any combiner, and false for a non-gaussian response. The bridge gates
  /// its multi-forest refusals on it.
  bool supportsResponseMutation;
  /// Whether the forest coupling admits a caller-supplied per-forest,
  /// per-observation weight (BCF alone today). Derived from the same predicate
  /// Chain::setForestWeights refuses on, so the bridge's capability probe and
  /// the engine's refusal cannot disagree. Internal, invisible to dbarts.h.
  bool supportsForestWeights;
  /// Whether the forest coupling owns a count-matrix response that setCounts
  /// can replace (the multinomial softmax alone today). The category count K
  /// rides numReportedLocations, which already is K for that coupling and 1 for
  /// every additive model, so the bridge reads the pair rather than a second
  /// name for the same number. Internal, invisible to dbarts.h.
  bool supportsCountsMutation;
  /// Whether the recorded test-fit channel carries a defined value: false only
  /// for BCF (no test treatment vector to blend off-sample). The bridge's
  /// test-surface refusal gates on it, so multinomial and single-forest
  /// samplers are allowed while BCF stays refused. Internal, invisible to
  /// dbarts.h.
  bool testFitsAreDefined;
  /// Whether the recorded per-forest fits and glue channels carry defined
  /// values: true only for a coupling that composes its forests through scalar
  /// glue (BCF today), so a run over any other model allocates neither. The run
  /// bridge reads it to decide whether the two channels exist; internal,
  /// invisible to dbarts.h.
  bool forestReportingIsDefined;
};

// Bridge entry points longjmp out of Rf_error past every destructor, so no
// shape field may own storage.
static_assert(std::is_trivially_destructible_v<SamplerShape>,
              "SamplerShape must not own storage: Rf_error skips its dtor");

/// Type-erased boundary for hosts (the R bridge): all type information is
/// resolved once, in the factory; every later call is one virtual hop into
/// fully typed code.
class SamplerBase {
public:
  virtual ~SamplerBase() = default;

  /// The sampler's counts and capabilities, assembled on every call. Fetch it
  /// once at the top of an entry point rather than per use.
  virtual SamplerShape shape() const = 0;

  virtual bool run(std::size_t numBurnIn, std::size_t numSamples,
                   Results& results,
                   const std::function<bool()>& pollInterrupt = {},
                   const SweepCallback& onSweep = {}) = 0;
  virtual void setOffset(const double* offset, bool updateScale) = 0;
  virtual void setResponse(const double* y, bool updateScale) = 0;
  virtual void setWeights(const double* weights) = 0;
  virtual void setSigma(double sigmaOriginalScale) = 0;
  /// Build the typed test store from a borrowed predictor view against the
  /// training cut grid, owning its raw. Returns false without touching the
  /// test store when a designated leaf covariate column would be CSC-backed.
  virtual bool setTestData(const PredictorSource& source) = 0;
  /// Dense convenience spelling (the dbarts.h shape): a plain column-major
  /// test matrix. Every column is dense, so the CSC-backed leaf-covariate
  /// refusal cannot fire.
  void setTestPredictors(const double* x_test,
                         std::size_t numTestObservations) {
    setTestData(densePredictorSource(x_test, numTestObservations,
                                     data().numPredictors));
  }
  virtual void setTestOffset(const double* testOffset) = 0;
  virtual void setData(const double* x, const double* y,
                       std::size_t numObservations, const double* weights,
                       const double* offset, const double* x_test,
                       std::size_t numTestObservations,
                       const double* testOffset) = 0;
  /// Replace the predictors from a borrowed view; only a dense block is
  /// consumable (PredictorSource::isDenseBlock), since every mutation kernel
  /// indexes the values column-major.
  virtual PredictorUpdateResult setPredictor(const PredictorSource& newX,
                                             bool forceUpdate,
                                             bool updateCutPoints) = 0;
  virtual PredictorUpdateResult updatePredictor(
    const PredictorSource& newColumns, const std::size_t* columns,
    std::size_t numColumns, bool forceUpdate, bool updateCutPoints) = 0;
  /// Dense convenience spellings (the dbarts.h shape): plain column-major
  /// blocks over the store's own row count.
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) {
    const ColumnStore& store = data();
    return setPredictor(densePredictorSource(newX, store.numObservations,
                                             store.numPredictors),
                        forceUpdate, updateCutPoints);
  }
  PredictorUpdateResult updatePredictor(const double* newColumns,
                                        const std::size_t* columns,
                                        std::size_t numColumns,
                                        bool forceUpdate,
                                        bool updateCutPoints) {
    return updatePredictor(
      densePredictorSource(newColumns, data().numObservations, numColumns),
      columns, numColumns, forceUpdate, updateCutPoints);
  }
  virtual void setCutPoints(const double* const* newCutPoints,
                            const std::uint32_t* numCutPoints,
                            const std::size_t* columns,
                            std::size_t numColumns,
                            const double* currentPredictors) = 0;
  virtual bool updatePredictorPerObservation(const double* newColumn,
                                             std::size_t column,
                                             bool* installed) = 0;
  virtual std::unique_ptr<PredictorUpdateSession> beginPredictorUpdate(
    const double* newColumn, std::size_t column) = 0;
  virtual std::size_t currentSampleNum() const = 0;
  /// forestIndex selects the forest to read (0 for every non-BCF sampler,
  /// 1 for the BCF treatment forest).
  virtual const std::vector<FlatNode>& savedTree(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex = 0) const = 0;
  /// Slopes of one saved tree, parallel to savedTree's pre-order leaves;
  /// meaningful only for vector-parameter leaf models.
  virtual const std::vector<double>& savedTreeSlopes(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex = 0) const = 0;
  /// Flattened mask words of one saved tree; meaningful only when the store
  /// has wide categorical columns.
  virtual const std::vector<std::uint64_t>& savedTreeMasks(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex = 0) const = 0;
  /// slopes, when non-null, receives a vector-parameter tree's slopes
  /// (one block per leaf in pre-order); cleared for scalar leaf models.
  /// masks, when non-null, receives the wide categorical side channel.
  virtual void flattenTree(std::size_t chainNum, std::size_t treeNum,
                           std::vector<FlatNode>& nodes,
                           std::vector<std::uint32_t>& counts,
                           std::vector<double>* slopes = nullptr,
                           std::vector<std::uint64_t>* masks = nullptr,
                           std::size_t forestIndex = 0) = 0;
  /// Fits for new rows of a borrowed predictor view; a CSC-backed view routes
  /// its rows resident, without a dense n x p materialization.
  virtual void predict(const PredictorSource& source,
                       std::size_t numTestObservations, double* out) = 0;
  /// Heteroscedastic variance surface s^2(x) on new rows (original scale);
  /// SamplerShape::hasVarianceForest gates it.
  virtual void predictVariance(const PredictorSource& source,
                               std::size_t numTestObservations,
                               double* out) = 0;
  /// Dense convenience spellings (the dbarts.h shape): a plain column-major
  /// block of new rows, as setTestPredictors takes.
  void predict(const double* x_test, std::size_t numTestObservations,
               double* out) {
    predict(densePredictorSource(x_test, numTestObservations,
                                 data().numPredictors),
            numTestObservations, out);
  }
  void predictVariance(const double* x_test, std::size_t numTestObservations,
                       double* out) {
    predictVariance(densePredictorSource(x_test, numTestObservations,
                                         data().numPredictors),
                    numTestObservations, out);
  }
  virtual void getState(SamplerStateData& state) = 0;
  /// currentPredictors supplies raw for a cross-grid restore's re-quantization
  /// (null for a same-spec continuation, which re-quantizes nothing).
  virtual bool setState(const SamplerStateData& state,
                        const double* currentPredictors) = 0;
  virtual WarmStartResult installForests(
      const SamplerStateData& donor,
      const std::vector<std::pair<std::size_t, int>>& sampleMap) = 0;
  virtual void sampleTreesFromPrior() = 0;
  virtual void sampleNodeParametersFromPrior() = 0;
  virtual void growFromRoot(std::size_t numSweeps) = 0;
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
  virtual const ColumnStore& data() const = 0;
  virtual const double* latents(std::size_t chainNum) const = 0;
  virtual double sigma(std::size_t chainNum) const = 0;
  /// BCF surface (docs/design/bcf.md); no-op/false off BCF. out receives
  /// {a, b0, b1}; forestTotalFits writes numObservations internal-scale fits.
  virtual void setTreatment(const double* z) = 0;
  /// Installs (or clears, at a null vector) a borrowed per-observation weight
  /// on forest forestIndex in every chain; false, installing nothing, when the
  /// coupling admits none or the index names no forest. Chain::setForestWeights
  /// states the semantics.
  virtual bool setForestWeights(std::size_t forestIndex,
                                const double* weights) = 0;
  /// Installs a borrowed replacement count matrix (n x K, category-major) and
  /// its per-observation trials in every chain; false, installing nothing, off
  /// a counts-owning coupling. n and K are fixed at creation, so the host
  /// validates the shape. Chain::setCounts states the semantics.
  virtual bool setCounts(const int* counts, const int* trials) = 0;
  /// Installs a borrowed n x K category offset (category-major) on the forest
  /// coupling's linear predictor in every chain, clearing it at a null pointer;
  /// false, installing nothing, off a counts-owning coupling. Distinct from
  /// setOffset, which the response model adds after the forests are combined.
  /// Chain::setCategoryOffset states the semantics.
  virtual bool setCategoryOffset(const double* offset) = 0;
  virtual bool bcfGlue(std::size_t chainNum, double* out) const = 0;
  virtual void forestTotalFits(std::size_t chainNum, std::size_t forestIndex,
                               double* out) const = 0;
  /// Forest forestIndex's per-predictor split usage into out (numPredictors
  /// entries); the per-forest analog of the recorded variable-count channel.
  virtual void forestVariableCounts(std::size_t chainNum,
                                    std::size_t forestIndex,
                                    std::uint32_t* out) const = 0;
  /// Tree count of forest forestIndex; equals SamplerShape::numTrees for
  /// forest 0.
  virtual std::size_t numTreesInForest(std::size_t forestIndex) const = 0;
};

template <IntegrableLeafModel L, typename ResidT = double>
class SamplerFacade final : public SamplerBase {
public:
  template <typename... Args>
  explicit SamplerFacade(Args&&... args) : impl_(std::forward<Args>(args)...) {}

  SamplerShape shape() const override {
    SamplerShape s;
    s.numObservations = impl_.numObservations();
    s.numPredictors = impl_.numPredictors();
    s.numTestObservations = impl_.numTestObservations();
    s.numChains = impl_.numChains();
    s.numThreads = impl_.numThreads();
    s.numTrees = impl_.numTrees();
    s.numForests = impl_.numForests();
    s.numGroups = impl_.numGroups();
    s.numLeafCovariates = impl_.numLeafCovariates();
    s.leafCovariateColumns = impl_.leafCovariateColumns();
    s.numReportedLocations = impl_.numReportedLocations();
    s.numVariableCountForests = impl_.numVariableCountForests();
    s.numCutpoints = impl_.numCutpoints();
    s.savedTreeCapacity = impl_.savedTreeCapacity();
    s.family = impl_.family();
    s.hasVarianceForest = impl_.hasVarianceForest();
    s.varianceLeafPrior = impl_.varianceLeafPrior();
    s.usesFunctionLeaves = Sampler<L, ResidT>::usesFunctionLeaves();
    s.kIsSampled = impl_.kIsSampled();
    s.usesDart = impl_.usesDart();
    s.supportsResponseMutation = impl_.supportsResponseMutation();
    s.supportsForestWeights = impl_.supportsForestWeights();
    s.supportsCountsMutation = impl_.supportsCountsMutation();
    s.testFitsAreDefined = impl_.testFitsAreDefined();
    s.forestReportingIsDefined = impl_.forestReportingIsDefined();
    return s;
  }

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
  bool setTestData(const PredictorSource& source) override {
    return impl_.setTestData(source);
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
  // the view-taking overrides hide the base's dense convenience spellings;
  // re-expose them so a dense caller still resolves
  using SamplerBase::setPredictor;
  using SamplerBase::updatePredictor;
  PredictorUpdateResult setPredictor(const PredictorSource& newX,
                                     bool forceUpdate,
                                     bool updateCutPoints) override {
    return impl_.setPredictor(newX, forceUpdate, updateCutPoints);
  }
  PredictorUpdateResult updatePredictor(
    const PredictorSource& newColumns, const std::size_t* columns,
    std::size_t numColumns, bool forceUpdate,
    bool updateCutPoints) override {
    return impl_.updatePredictor(newColumns, columns, numColumns, forceUpdate,
                                 updateCutPoints);
  }
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints,
                    const std::size_t* columns,
                    std::size_t numColumns,
                    const double* currentPredictors) override {
    impl_.setCutPoints(newCutPoints, numCutPoints, columns, numColumns,
                       currentPredictors);
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
  std::size_t currentSampleNum() const override {
    return impl_.currentSampleNum();
  }
  const std::vector<FlatNode>& savedTree(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex) const override {
    return impl_.savedTree(chainNum, slot, treeNum, forestIndex);
  }
  const std::vector<double>& savedTreeSlopes(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex) const override {
    return impl_.savedTreeSlopes(chainNum, slot, treeNum, forestIndex);
  }
  const std::vector<std::uint64_t>& savedTreeMasks(
    std::size_t chainNum, std::size_t slot, std::size_t treeNum,
    std::size_t forestIndex) const override {
    return impl_.savedTreeMasks(chainNum, slot, treeNum, forestIndex);
  }
  void flattenTree(std::size_t chainNum, std::size_t treeNum,
                   std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes,
                   std::vector<std::uint64_t>* masks,
                   std::size_t forestIndex) override {
    impl_.flattenTree(chainNum, treeNum, nodes, counts, slopes, masks,
                      forestIndex);
  }
  // as for setPredictor: the view-taking overrides hide the base's dense
  // convenience spellings, so re-expose them
  using SamplerBase::predict;
  using SamplerBase::predictVariance;
  void predict(const PredictorSource& source, std::size_t numTestObservations,
               double* out) override {
    impl_.predict(source, numTestObservations, out);
  }
  void predictVariance(const PredictorSource& source,
                       std::size_t numTestObservations, double* out) override {
    impl_.predictVariance(source, numTestObservations, out);
  }
  void getState(SamplerStateData& state) override { impl_.getState(state); }
  bool setState(const SamplerStateData& state,
                const double* currentPredictors) override {
    return impl_.setState(state, currentPredictors);
  }
  WarmStartResult installForests(
      const SamplerStateData& donor,
      const std::vector<std::pair<std::size_t, int>>& sampleMap) override {
    return impl_.installForests(donor, sampleMap);
  }
  void sampleTreesFromPrior() override { impl_.sampleTreesFromPrior(); }
  void sampleNodeParametersFromPrior() override {
    impl_.sampleNodeParametersFromPrior();
  }
  void growFromRoot(std::size_t numSweeps) override {
    impl_.growFromRoot(numSweeps);
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
  const ColumnStore& data() const override { return impl_.data(); }
  const double* latents(std::size_t chainNum) const override {
    return impl_.latents(chainNum);
  }
  double sigma(std::size_t chainNum) const override {
    return impl_.sigma(chainNum);
  }
  void setTreatment(const double* z) override { impl_.setTreatment(z); }
  bool setForestWeights(std::size_t forestIndex,
                        const double* weights) override {
    return impl_.setForestWeights(forestIndex, weights);
  }
  bool setCounts(const int* counts, const int* trials) override {
    return impl_.setCounts(counts, trials);
  }
  bool setCategoryOffset(const double* offset) override {
    return impl_.setCategoryOffset(offset);
  }
  bool bcfGlue(std::size_t chainNum, double* out) const override {
    return impl_.bcfGlue(chainNum, out);
  }
  void forestTotalFits(std::size_t chainNum, std::size_t forestIndex,
                       double* out) const override {
    impl_.forestTotalFits(chainNum, forestIndex, out);
  }
  void forestVariableCounts(std::size_t chainNum, std::size_t forestIndex,
                            std::uint32_t* out) const override {
    impl_.forestVariableCounts(chainNum, forestIndex, out);
  }
  std::size_t numTreesInForest(std::size_t forestIndex) const override {
    return impl_.numTreesInForest(forestIndex);
  }

  Sampler<L, ResidT>& impl() { return impl_; }
  const Sampler<L, ResidT>& impl() const { return impl_; }

private:
  Sampler<L, ResidT> impl_;
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

  std::size_t numObservations = samplers[0]->shape().numObservations;

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
///
/// The opt-in fp32 running residual (options.fp32Residual,
/// docs/design/reduced-precision-storage.md sec 3b) mints ONE extra
/// instantiation - the gaussian constant leaf with ResidT = float - and only
/// on that branch; every other family/leaf keeps the byte-identical fp64
/// engine. The host (R bridge) validates the gaussian-constant-leaf gate before
/// setting the flag; a non-gaussian request never reaches the float branch.
inline std::unique_ptr<SamplerBase> createConstantLeafSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  ResponseFamily family, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* const* rngs) {
  if (options.fp32Residual && family == ResponseFamily::gaussian)
    return std::make_unique<SamplerFacade<ConstantGaussianLeaf, float>>(
      x, y, numObservations, numPredictors, weights, offset, family,
      sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
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

/// Whether any predictor carries a nonzero monotone direction (the constrained
/// leaf's construction-time selector); a null spec or an all-zero vector keeps
/// the unchanged constant-leaf path. A direction on a categorical column is
/// refused (monotonicity is undefined on unordered codes).
inline bool monotoneConstraintIsActive(const SamplerOptions& options,
                                       std::size_t numPredictors) {
  if (options.monotoneDirections == nullptr) return false;
  bool active = false;
  for (std::size_t j = 0; j < numPredictors; ++j) {
    if (options.monotoneDirections[j] == 0) continue;
    if (options.predictors.typeOf(j) == ColumnType::categorical) return false;
    active = true;
  }
  return active;
}

/// Dispatch on the leaf model: designated leaf covariates select the
/// linear-leaf instantiation - or the GP one under options.gpLeaves -
/// anything else the constant leaf. Returns null on an invalid
/// designation - more than maxNumCovariates columns (both leaf models share
/// the bound), a column out of range, or a categorical column (category
/// codes are unordered; interact through splits instead) - which the host
/// turns into its own error.
inline std::unique_ptr<SamplerBase> createSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  ResponseFamily family, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* const* rngs) {
  // the heteroscedastic variance forest is gaussian + plain-constant-leaf only:
  // the latent families own the weight channel it routes through (a collision),
  // and v1 keeps the mean leaf constant. Refuse every other combination here,
  // before any Chain is built (docs/design/heteroscedastic.md section 5).
  if (options.numVarianceTrees > 0 &&
      (family != ResponseFamily::gaussian || options.numLeafCovariates != 0 ||
       monotoneConstraintIsActive(options, numPredictors)))
    return nullptr;
  if (options.numLeafCovariates == 0 &&
      monotoneConstraintIsActive(options, numPredictors))
    return std::make_unique<SamplerFacade<MonotoneConstantGaussianLeaf>>(
      x, y, numObservations, numPredictors, weights, offset, family,
      sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
  if (options.numLeafCovariates == 0)
    return createConstantLeafSampler(x, y, numObservations, numPredictors,
                                     weights, offset, family, sigmaEstimate,
                                     sigmaDf, sigmaRawScale, options, rngs);

  // CSC-backed columns hold no contiguous raw values for a leaf model: a
  // design whose every column is CSC-backed is refused outright, a mixed one
  // per designated column
  if (predictorSourceIsAllCsc(options.predictors, numPredictors))
    return nullptr;
  if (!leafCovariateDesignationIsValid(
        options, numPredictors,
        [&](std::size_t j) {
          return options.predictors.typeOf(j) == ColumnType::categorical;
        },
        [&](std::size_t j) { return options.predictors.sourceOf(j) >= 0; }))
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
  // only the single-forest chain builds a variance forest, so accepting the
  // option at these two factories would drop it silently; refuse it as
  // createSampler does
  if (options.numVarianceTrees > 0) return nullptr;
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, sigmaEstimate,
    sigmaDf, sigmaRawScale, options, spec, rngs);
}

/// A K-forest multinomial (softmax) sampler (docs/design/multinomial.md):
/// constant-leaf only, so the single instantiation. rngs supplies one
/// generator per chain.
inline std::unique_ptr<SamplerBase> createMultinomialSampler(
  const double* x, std::size_t numObservations, std::size_t numPredictors,
  const SamplerOptions& options, const MultinomialSpec& spec,
  ext_rng* const* rngs) {
  if (options.numVarianceTrees > 0) return nullptr;  // as createBCFSampler
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, numObservations, numPredictors, options, spec, rngs);
}

}  // namespace bartcore

#endif  // BARTCORE_FACADE_HPP
