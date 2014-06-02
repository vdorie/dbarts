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
  CGMPrior* bart_createCGMPriorFromControl(const Control* control) {
    return new CGMPrior(*control);
  }
  void bart_destroyCGMPrior(CGMPrior* prior) {
    delete prior;
  }
  void bart_initializeCGMPriorFromControl(CGMPrior* prior, const Control* control)
  {
    new (prior) CGMPrior(*control);
  }
  void bart_invalidateCGMPrior(CGMPrior* prior) {
    prior->~CGMPrior();
  }
  
  NormalPrior* bart_createNormalPrior() {
    return new NormalPrior;
  }
  NormalPrior* bart_createNormalPriorFromControl(const Control* control) {
    return new NormalPrior(*control);
  }
  void bart_destroyNormalPrior(NormalPrior* prior) {
    delete prior;
  }
  void bart_initializeNormalPriorFromControl(NormalPrior* prior, const Control* control)
  {
    new (prior) NormalPrior(*control);
  }
  void bart_invalidateNormalPrior(NormalPrior* prior) {
    prior->~NormalPrior();
  }
  
  ChiSquaredPrior* bart_createChiSquaredPrior() {
    return new ChiSquaredPrior;
  }
  ChiSquaredPrior* bart_createChiSquaredPriorFromControl(const Control* control) {
    return new ChiSquaredPrior(*control);
  }
  void bart_destroyChiSquaredPrior(ChiSquaredPrior* prior) {
    delete prior;
  }
  void bart_initializeChiSquaredPriorFromControl(ChiSquaredPrior* prior, const Control* control)
  {
    new (prior) ChiSquaredPrior(*control);
  }
  void bart_invalidateChiSquaredPrior(ChiSquaredPrior* prior) {
    prior->~ChiSquaredPrior();
  }
}
