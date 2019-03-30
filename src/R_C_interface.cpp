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
  
  void dbarts_setRNGState(BARTFit* fit, const void* const* uniformState, const void* const* normalState) {
    fit->setRNGState(uniformState, normalState);
  }
  
  Results* dbarts_runSampler(BARTFit* fit) {
    return fit->runSampler();
  }
  
  Results* dbarts_runSamplerForIterations(BARTFit* fit, size_t numBurnIn, size_t numSamples) {
    return fit->runSampler(numBurnIn, numSamples);
  }
  
  void dbarts_sampleTreesFromPrior(BARTFit* fit) {
    fit->sampleTreesFromPrior();
  }
  
  void dbarts_setResponse(BARTFit* fit, const double* newResponse) {
    fit->setResponse(newResponse);
  }
  
  void dbarts_setOffset(BARTFit* fit, const double* newOffset) {
    fit->setOffset(newOffset);
  }
  
  
  int dbarts_setPredictor(BARTFit* fit, const double* newPredictor, int forceUpdate, int updateCutPoints) {
    return fit->setPredictor(newPredictor, forceUpdate, updateCutPoints);
  }
  
  int dbarts_updatePredictor(BARTFit* fit, const double* newPredictor, size_t column, int forceUpdate, int updateCutPoints) {
    return fit->updatePredictor(newPredictor, &column, 1, forceUpdate, updateCutPoints);
  }
  
  int dbarts_updatePredictors(BARTFit* fit, const double* newPredictor, const size_t* columns, size_t numColumns, int forceUpdate, int updateCutPoints)
  {
    return fit->updatePredictor(newPredictor, columns, numColumns, forceUpdate, updateCutPoints);
  }
  
  void dbarts_setTestPredictor(BARTFit* fit, const double* newTestPredictor, size_t numTestObservations)
  {
    fit->setTestPredictor(newTestPredictor, numTestObservations);
  }
  
  void dbarts_setTestOffset(BARTFit* fit, const double* newTestOffset)
  {
    fit->setTestOffset(newTestOffset);
  }
  
  void dbarts_setTestPredictorAndOffset(BARTFit* fit, const double* newTestPredictor, const double* newTestOffset, size_t numTestObservations)
  {
    fit->setTestPredictorAndOffset(newTestPredictor, newTestOffset, numTestObservations);
  }
  
  void dbarts_updateTestPredictor(dbarts::BARTFit* fit, const double* newTestPredictor, size_t column)
  {
    fit->updateTestPredictor(newTestPredictor, column);
  }
  
  void dbarts_updateTestPredictors(dbarts::BARTFit* fit, const double* newTestPredictor, const size_t* columns, size_t numColumns)
  {
    fit->updateTestPredictors(newTestPredictor, columns, numColumns);
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
  NormalPrior* dbarts_createNormalPriorFromOptions(const Control* control, const Model* model, double k) {
    return new NormalPrior(*control, *model, k);
  }
  void dbarts_destroyNormalPrior(NormalPrior* prior) {
    delete prior;
  }
  void dbarts_initializeNormalPriorFromOptions(NormalPrior* prior, const Control* control, const Model* model, double k)
  {
    new (prior) NormalPrior(*control, *model, k);
  }
  void dbarts_invalidateNormalPrior(NormalPrior* prior) {
    prior->~NormalPrior();
  }
  
  NormalHyperprior* dbarts_createNormalHyperprior() {
    return new NormalHyperprior;
  }
  NormalHyperprior* dbarts_createNormalHyperpriorFromOptions(double scale) {
    return new NormalHyperprior(scale);
  }
  void dbarts_destroyNormalHyperprior(NormalHyperprior* prior) {
    delete prior;
  }
  void dbarts_initializeNormalHyperpriorFromOptions(NormalHyperprior* prior, double scale)
  {
    new (prior) NormalHyperprior(scale);
  }
  void dbarts_invalidateNormalHyperprior(NormalHyperprior* prior) {
    prior->~NormalHyperprior();
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
