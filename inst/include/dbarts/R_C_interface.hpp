#ifndef DBARTS_R_C_INTERFACE_HPP
#define DBARTS_R_C_INTERFACE_HPP

#include <cstddef> // size_t

namespace dbarts {
  struct Control;
  struct Model;
  struct Data;
  
  struct CGMPrior;
  struct NormalPrior;
  struct ChiSquaredPrior;
  
  struct BARTFit;
  struct Results;
}

extern "C" {
  // pair calls of create<->destroy, initialize<->invalidate
  dbarts::BARTFit* dbarts_createFit(dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_initializeFit(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_destroyFit(dbarts::BARTFit* fit);
  void dbarts_invalidateFit(dbarts::BARTFit* fit);
  
  dbarts::Results* dbarts_runSampler(dbarts::BARTFit* fit);
  dbarts::Results* dbarts_runSamplerForIterations(dbarts::BARTFit* fit, std::size_t numBurnIn, std::size_t numSamples);
  void dbarts_setResponse(dbarts::BARTFit* fit, const double* newResponse);
  void dbarts_setOffset(dbarts::BARTFit* fit, const double* newOffset);
  void dbarts_setPredictor(dbarts::BARTFit* fit, const double* newPredictor, std::size_t predictorColumn);
  void dbarts_setTestPredictor(dbarts::BARTFit* fit, const double* newPredictor, std::size_t predictorColumn);
  void dbarts_setTestPredictors(dbarts::BARTFit* fit, const double* newPredictors, std::size_t numTestObservations);
  
  dbarts::CGMPrior* dbarts_createCGMPrior();
  dbarts::CGMPrior* dbarts_createCGMPriorFromOptions(double base, double power);
  void dbarts_destroyCGMPrior(dbarts::CGMPrior* prior);
  void dbarts_initializeCGMPriorFromOptions(dbarts::CGMPrior* prior, double base, double power);
  void dbarts_invalidateCGMPrior(dbarts::CGMPrior* prior);
  
  dbarts::NormalPrior* dbarts_createNormalPrior();
  dbarts::NormalPrior* dbarts_createNormalPriorFromOptions(const dbarts::Control* control, double k);
  void dbarts_destroyNormalPrior(dbarts::NormalPrior* prior);
  void dbarts_initializeNormalPriorFromOptions(dbarts::NormalPrior* prior, const dbarts::Control* control, double k);
  void dbarts_invalidateNormalPrior(dbarts::NormalPrior* prior);
  
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPrior();
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPriorFromOptions(double degreesOfFreedom, double quantile);
  void dbarts_destroyChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
  void dbarts_initializeChiSquaredPriorFromOptions(dbarts::ChiSquaredPrior* prior, double degreesOfFreedom, double quantile);
  void dbarts_invalidateChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
}

#endif // DBARTS_R_C_INTERFACE_HPP
