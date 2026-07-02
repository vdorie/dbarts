#ifndef BARTCORE_FACADE_HPP
#define BARTCORE_FACADE_HPP

#include <cstddef>
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

  virtual void run(std::size_t numBurnIn, std::size_t numSamples,
                   Results& results) = 0;
  virtual void setOffset(const double* offset, bool updateScale) = 0;
  virtual void setResponse(const double* y) = 0;
  virtual void setSigma(double sigmaOriginalScale) = 0;
  virtual void setTestPredictors(const double* x_test,
                                 std::size_t numTestObservations) = 0;
  virtual void setData(const double* x, const double* y,
                       std::size_t numObservations, const double* weights,
                       const double* offset, const double* x_test,
                       std::size_t numTestObservations) = 0;
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
  virtual ext_rng* rng() const = 0;
  virtual const double* latents(std::size_t chainNum) const = 0;
  virtual double sigma(std::size_t chainNum) const = 0;
  virtual bool kIsSampled() const = 0;
  virtual std::size_t numChains() const = 0;
  virtual std::size_t numObservations() const = 0;
  virtual std::size_t numPredictors() const = 0;
  virtual std::size_t numTestObservations() const = 0;
};

template <IntegrableLeafModel L>
class SamplerFacade final : public SamplerBase {
public:
  template <typename... Args>
  explicit SamplerFacade(Args&&... args) : impl_(std::forward<Args>(args)...) {}

  void run(std::size_t numBurnIn, std::size_t numSamples,
           Results& results) override {
    impl_.run(numBurnIn, numSamples, results);
  }
  void setOffset(const double* offset, bool updateScale) override {
    impl_.setOffset(offset, updateScale);
  }
  void setResponse(const double* y) override { impl_.setResponse(y); }
  void setSigma(double sigmaOriginalScale) override {
    impl_.setSigma(sigmaOriginalScale);
  }
  void setTestPredictors(const double* x_test,
                         std::size_t numTestObservations) override {
    impl_.setTestPredictors(x_test, numTestObservations);
  }
  void setData(const double* x, const double* y, std::size_t numObservations,
               const double* weights, const double* offset,
               const double* x_test, std::size_t numTestObservations) override {
    impl_.setData(x, y, numObservations, weights, offset, x_test,
                  numTestObservations);
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
  ext_rng* rng() const override { return impl_.rng(); }
  const double* latents(std::size_t chainNum) const override {
    return impl_.latents(chainNum);
  }
  double sigma(std::size_t chainNum) const override {
    return impl_.sigma(chainNum);
  }
  bool kIsSampled() const override { return impl_.kIsSampled(); }
  std::size_t numChains() const override { return impl_.numChains(); }
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

/// The instantiation matrix, phase 2: one leaf model. rngs supplies one
/// generator per chain (options.numChains of them).
inline std::unique_ptr<SamplerBase> createClassicSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  bool responseIsBinary, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* const* rngs) {
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, responseIsBinary,
    sigmaEstimate, sigmaDf, sigmaRawScale, options, rngs);
}

}  // namespace bartcore

#endif  // BARTCORE_FACADE_HPP
