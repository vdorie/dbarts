#include "config.hpp"
#include <dbarts/R_C_interface.hpp>

#include <new>
#include <cmath>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#include <external/stats.h>

using namespace dbarts;

extern "C" {
  BARTFit* dbarts_createFit(Control* control, Model* model, Data* data) {
    return new BARTFit(*control, *model, *data);
  }
  void dbarts_initializeFit(BARTFit* fit, Control* control, Model* model, Data* data) {
    new (fit) BARTFit(*control, *model, *data);
  }
  void dbarts_destroyFit(BARTFit* fit) {
    delete fit;
  }
  void dbarts_invalidateFit(BARTFit* fit) {
    fit->~BARTFit();
  }
  
  
  Results* dbarts_runSampler(BARTFit* fit) {
    return fit->runSampler();
  }
  
  Results* dbarts_runSamplerForIterations(BARTFit* fit, size_t numBurnIn, size_t numSamples) {
    return fit->runSampler(numBurnIn, numSamples);
  }
  
  void dbarts_setResponse(BARTFit* fit, const double* newResponse) {
    fit->setResponse(newResponse);
  }
  
  void dbarts_setOffset(BARTFit* fit, const double* newOffset) {
    fit->setOffset(newOffset);
  }
  
  void dbarts_setPredictor(BARTFit* fit, const double* newPredictor, size_t predictorColumn) {
    fit->setPredictor(newPredictor, predictorColumn);
  }
  
  void dbarts_setTestPredictor(BARTFit* fit, const double* newPredictor, size_t predictorColumn) {
    fit->setTestPredictor(newPredictor, predictorColumn);
  }
  
  void dbarts_setTestPredictors(BARTFit* fit, const double* newPredictors, size_t numTestObservations) {
    fit->setTestPredictors(newPredictors, numTestObservations);
  }
  
  void dbarts_setTestOffset(BARTFit* fit, const double* newTestOffset) {
    fit->setTestOffset(newTestOffset);
  }
  
  void dbarts_setTestPredictorsAndOffset(BARTFit* fit, const double* newPredictors, const double* newTestOffset, size_t numTestObservations) {
    fit->setTestPredictors(newPredictors, newTestOffset, numTestObservations);
  }
  
  CGMPrior* dbarts_createCGMPrior() {
    return new CGMPrior;
  }
  CGMPrior* dbarts_createCGMPriorFromOptions(double base, double power) {
    return new CGMPrior(base, power);
  }
  void dbarts_destroyCGMPrior(CGMPrior* prior) {
    delete prior;
  }
  void dbarts_initializeCGMPriorFromOptions(CGMPrior* prior, double base, double power)
  {
    new (prior) CGMPrior(base, power);
  }
  void dbarts_invalidateCGMPrior(CGMPrior* prior) {
    prior->~CGMPrior();
  }
  
  NormalPrior* dbarts_createNormalPrior() {
    return new NormalPrior;
  }
  NormalPrior* dbarts_createNormalPriorFromOptions(const Control* control, double k) {
    return new NormalPrior(*control, k);
  }
  void dbarts_destroyNormalPrior(NormalPrior* prior) {
    delete prior;
  }
  void dbarts_initializeNormalPriorFromOptions(NormalPrior* prior, const Control* control, double k)
  {
    new (prior) NormalPrior(*control, k);
  }
  void dbarts_invalidateNormalPrior(NormalPrior* prior) {
    prior->~NormalPrior();
  }
  
  ChiSquaredPrior* dbarts_createChiSquaredPrior() {
    return new ChiSquaredPrior;
  }
  ChiSquaredPrior* dbarts_createChiSquaredPriorFromOptions(double degreesOfFreedom, double quantile) {
    return new ChiSquaredPrior(degreesOfFreedom, quantile);
  }
  void dbarts_destroyChiSquaredPrior(ChiSquaredPrior* prior) {
    delete prior;
  }
  void dbarts_initializeChiSquaredPriorFromOptions(ChiSquaredPrior* prior, double degreesOfFreedom, double quantile)
  {
    new (prior) ChiSquaredPrior(degreesOfFreedom, quantile);
  }
  void dbarts_invalidateChiSquaredPrior(ChiSquaredPrior* prior) {
    prior->~ChiSquaredPrior();
  }
}
