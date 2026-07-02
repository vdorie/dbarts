#ifndef BARTCORE_FACADE_HPP
#define BARTCORE_FACADE_HPP

#include <cstddef>
#include <memory>

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
  virtual const double* latents() const = 0;
  virtual double sigma() const = 0;
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
  const double* latents() const override { return impl_.latents(); }
  double sigma() const override { return impl_.sigma(); }
  std::size_t numObservations() const override { return impl_.numObservations(); }
  std::size_t numPredictors() const override { return impl_.numPredictors(); }
  std::size_t numTestObservations() const override {
    return impl_.numTestObservations();
  }

  Sampler<L>& impl() { return impl_; }

private:
  Sampler<L> impl_;
};

/// The instantiation matrix, phase 2: one leaf model.
inline std::unique_ptr<SamplerBase> createClassicSampler(
  const double* x, const double* y, std::size_t numObservations,
  std::size_t numPredictors, const double* weights, const double* offset,
  bool responseIsBinary, double sigmaEstimate, double sigmaDf,
  double sigmaRawScale, const SamplerOptions& options, ext_rng* rng) {
  return std::make_unique<SamplerFacade<ConstantGaussianLeaf>>(
    x, y, numObservations, numPredictors, weights, offset, responseIsBinary,
    sigmaEstimate, sigmaDf, sigmaRawScale, options, rng);
}

}  // namespace bartcore

#endif  // BARTCORE_FACADE_HPP
