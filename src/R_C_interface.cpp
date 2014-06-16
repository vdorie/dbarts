#include "config.hpp"
#include <bart/R_C_interface.hpp>

#include <new>
#include <cmath>

#include <bart/bartFit.hpp>
#include <bart/model.hpp>
#include <bart/results.hpp>

#include <external/stats.h>

using namespace bart;

extern "C" {
  BARTFit* bart_createFit(Control* control, Model* model, Data* data) {
    return new BARTFit(*control, *model, *data);
  }
  void bart_initializeFit(BARTFit* fit, Control* control, Model* model, Data* data) {
    new (fit) BARTFit(*control, *model, *data);
  }
  void bart_destroyFit(bart::BARTFit* fit) {
    delete fit;
  }
  void bart_invalidateFit(BARTFit* fit) {
    fit->~BARTFit();
  }
  
  
  Results* bart_runSampler(BARTFit* fit) {
    return fit->runSampler();
  }
  
  Results* bart_runSamplerForIterations(BARTFit* fit, size_t numBurnIn, size_t numSamples) {
    return fit->runSampler(numBurnIn, numSamples);
  }
  
  void bart_setResponse(BARTFit* fit, const double* newResponse) {
    fit->setResponse(newResponse);
  }
  
  CGMPrior* bart_createCGMPrior() {
    return new CGMPrior;
  }
  CGMPrior* bart_createCGMPriorFromOptions(double base, double power) {
    return new CGMPrior(base, power);
  }
  void bart_destroyCGMPrior(CGMPrior* prior) {
    delete prior;
  }
  void bart_initializeCGMPriorFromOptions(CGMPrior* prior, double base, double power)
  {
    new (prior) CGMPrior(base, power);
  }
  void bart_invalidateCGMPrior(CGMPrior* prior) {
    prior->~CGMPrior();
  }
  
  NormalPrior* bart_createNormalPrior() {
    return new NormalPrior;
  }
  NormalPrior* bart_createNormalPriorFromOptions(const Control* control, double k) {
    return new NormalPrior(*control, k);
  }
  void bart_destroyNormalPrior(NormalPrior* prior) {
    delete prior;
  }
  void bart_initializeNormalPriorFromOptions(NormalPrior* prior, const Control* control, double k)
  {
    new (prior) NormalPrior(*control, k);
  }
  void bart_invalidateNormalPrior(NormalPrior* prior) {
    prior->~NormalPrior();
  }
  
  ChiSquaredPrior* bart_createChiSquaredPrior() {
    return new ChiSquaredPrior;
  }
  ChiSquaredPrior* bart_createChiSquaredPriorFromOptions(double degreesOfFreedom, double quantile) {
    return new ChiSquaredPrior(degreesOfFreedom, quantile);
  }
  void bart_destroyChiSquaredPrior(ChiSquaredPrior* prior) {
    delete prior;
  }
  void bart_initializeChiSquaredPriorFromOptions(ChiSquaredPrior* prior, double degreesOfFreedom, double quantile)
  {
    new (prior) ChiSquaredPrior(degreesOfFreedom, quantile);
  }
  void bart_invalidateChiSquaredPrior(ChiSquaredPrior* prior) {
    prior->~ChiSquaredPrior();
  }
}
