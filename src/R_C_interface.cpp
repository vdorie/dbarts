#include "config.hpp"
#include <dbarts/R_C_interface.hpp>

#include <new>
#include <cmath>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#include <external/stats.h>

#include "R_interface_common.hpp"

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
  
  Control* dbarts_createControl(SEXP controlExpr) {
    Control* result = new Control;
    initializeControlFromExpression(*result, controlExpr);
    return result;
  }
  void dbarts_destroyControl(Control* control) {
    delete control;
  }
  void dbarts_initializeControl(Control* control, SEXP controlExpr) {
    initializeControlFromExpression(*control, controlExpr);
  }
  // void dbarts_invalidateControl(Control* control) { }
  
  void dbarts_setControl(dbarts::BARTFit* fit, const dbarts::Control* control) {
    fit->setControl(*control);
  }
  
  dbarts::Data* dbarts_createData(SEXP dataExpr) {
    Data* result = new Data;
    initializeDataFromExpression(*result, dataExpr);
    return result;
  }
  void dbarts_destroyData(dbarts::Data* data) {
    invalidateData(*data);
    delete data;
  }
  void dbarts_initializeData(dbarts::Data* data, SEXP dataExpr) {
    initializeDataFromExpression(*data, dataExpr);
  }
  void dbarts_invalidateData(dbarts::Data* data) {
    invalidateData(*data);
  }
  
  dbarts::Model* dbarts_createModel(SEXP modelExpr, const dbarts::Control* control, const dbarts::Data* data) {
    Model* result = new Model(control->responseIsBinary);
    initializeModelFromExpression(*result, modelExpr, *control, *data);
    return result;
  }
  void dbarts_destroyModel(dbarts::Model* model) {
    invalidateModel(*model);
    delete model;
  }
  void dbarts_initializeModel(dbarts::Model* model, SEXP modelExpr, const dbarts::Control* control, const dbarts::Data* data) {
    initializeModelFromExpression(*model, modelExpr, *control, *data);
  }
  void dbarts_invalidateModel(dbarts::Model* model) {
    invalidateModel(*model);
  }
  
  SEXP dbarts_createStateExpression(const dbarts::BARTFit* fit) {
    return createStateExpressionFromFit(*fit);
  }
  void dbarts_initializeState(dbarts::BARTFit* fit, SEXP stateExpr) {
    initializeStateFromExpression(*fit, stateExpr);
  }
  
  void dbarts_setRNGState(BARTFit* fit, const void* const* uniformState, const void* const* normalState) {
    fit->setRNGState(uniformState, normalState);
  }
  void dbarts_printInitialSummary(const dbarts::BARTFit* fit) {
    fit->printInitialSummary();
  }
  void dbarts_printTrees(const dbarts::BARTFit* fit,
                         const std::size_t* chainIndices,  std::size_t numChainIndices,
                         const std::size_t* sampleIndices, std::size_t numSampleIndices,
                         const std::size_t* treeIndices,   std::size_t numTreeIndices)
  {
    fit->printTrees(chainIndices,  numChainIndices,
                    sampleIndices, numSampleIndices,
                    treeIndices,   numTreeIndices);
  }

  FlattenedTrees* dbarts_getTrees(const dbarts::BARTFit* fit,
                                  const std::size_t* chainIndices,  std::size_t numChainIndices,
                                  const std::size_t* sampleIndices, std::size_t numSampleIndices,
                                  const std::size_t* treeIndices,   std::size_t numTreeIndices)
  {
    return fit->getFlattenedTrees(chainIndices,  numChainIndices,
                                  sampleIndices, numSampleIndices,
                                  treeIndices,   numTreeIndices);
  }

  Results* dbarts_runSampler(BARTFit* fit) {
    return fit->runSampler();
  }
  
  Results* dbarts_runSamplerForIterations(BARTFit* fit, size_t numBurnIn, size_t numSamples) {
    return fit->runSampler(numBurnIn, numSamples);
  }
  
  void dbarts_runSamplerWithResults(BARTFit* fit, size_t numBurnIn, Results* results) {
    fit->runSampler(numBurnIn, results);
  }
  
  void dbarts_sampleTreesFromPrior(BARTFit* fit) {
    fit->sampleTreesFromPrior();
  }
  
  void dbarts_sampleNodeParametersFromPrior(BARTFit* fit) {
    fit->sampleNodeParametersFromPrior();
  }
  
  void dbarts_predict(const BARTFit* fit, const double* x_test, std::size_t numTestObservations, const double* testOffset, double* result) {
    fit->predict(x_test, numTestObservations, testOffset, result);
  }
  
  void dbarts_setResponse(BARTFit* fit, const double* newResponse) {
    fit->setResponse(newResponse);
  }
  
  void dbarts_setOffset(BARTFit* fit, const double* newOffset, bool updateScale) {
    fit->setOffset(newOffset, updateScale);
  }
  
  void dbarts_setSigma(BARTFit* fit, const double* newSigma) {
    fit->setSigma(newSigma);
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
  
  void dbarts_storeLatents(const dbarts::BARTFit* fit, double* target)
  {
    fit->storeLatents(target);
  }
  
  CGMPrior* dbarts_createCGMPrior() {
    return new CGMPrior;
  }
  CGMPrior* dbarts_createCGMPriorFromOptions(double base, double power, const double* splitProbabilities) {
    return new CGMPrior(base, power, splitProbabilities);
  }
  void dbarts_destroyCGMPrior(CGMPrior* prior) {
    delete prior;
  }
  void dbarts_initializeCGMPriorFromOptions(CGMPrior* prior, double base, double power, const double* splitProbabilities)
  {
    new (prior) CGMPrior(base, power, splitProbabilities);
  }
  void dbarts_invalidateCGMPrior(CGMPrior* prior) {
    prior->~CGMPrior();
  }
  
  NormalPrior* dbarts_createNormalPrior() {
    return new NormalPrior;
  }
  NormalPrior* dbarts_createNormalPriorFromOptions(const Control* control, const Model* model) {
    return new NormalPrior(*control, *model);
  }
  void dbarts_destroyNormalPrior(NormalPrior* prior) {
    delete prior;
  }
  void dbarts_initializeNormalPriorFromOptions(NormalPrior* prior, const Control* control, const Model* model)
  {
    new (prior) NormalPrior(*control, *model);
  }
  void dbarts_invalidateNormalPrior(NormalPrior* prior) {
    prior->~NormalPrior();
  }
  
  ChiHyperprior* dbarts_createChiHyperprior() {
    return new ChiHyperprior;
  }
  ChiHyperprior* dbarts_createChiHyperpriorFromOptions(double degreesOfFreedom, double scale) {
    return new ChiHyperprior(degreesOfFreedom, scale);
  }
  void dbarts_destroyChiHyperprior(ChiHyperprior* prior) {
    delete prior;
  }
  void dbarts_initializeChiHyperpriorFromOptions(ChiHyperprior* prior, double degreesOfFreedom, double scale)
  {
    new (prior) ChiHyperprior(degreesOfFreedom, scale);
  }
  void dbarts_invalidateChiHyperprior(ChiHyperprior* prior) {
    prior->~ChiHyperprior();
  }
  
  FixedHyperprior* dbarts_createFixedHyperprior() {
    return new FixedHyperprior;
  }
  FixedHyperprior* dbarts_createFixedHyperpriorFromOptions(double k) {
    return new FixedHyperprior(k);
  }
  void dbarts_destroyFixedHyperprior(FixedHyperprior* prior) {
    delete prior;
  }
  void dbarts_initializeFixedHyperpriorFromOptions(FixedHyperprior* prior, double k)
  {
    new (prior) FixedHyperprior(k);
  }
  void dbarts_invalidateFixedHyperprior(FixedHyperprior* prior) {
    prior->~FixedHyperprior();
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
