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
  struct NormalHyperprior;
  
  struct BARTFit;
  struct Results;
}

extern "C" {
  // pair calls of create<->destroy, initialize<->invalidate
  dbarts::BARTFit* dbarts_createFit(dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_initializeFit(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_destroyFit(dbarts::BARTFit* fit);
  void dbarts_invalidateFit(dbarts::BARTFit* fit);
  
  void dbarts_setRNGState(dbarts::BARTFit* fit, const void* const* uniformState, const void* const* normalState);
  
  dbarts::Results* dbarts_runSampler(dbarts::BARTFit* fit);
  dbarts::Results* dbarts_runSamplerForIterations(dbarts::BARTFit* fit, std::size_t numBurnIn, std::size_t numSamples);
  void dbarts_sampleTreesFromPrior(dbarts::BARTFit* fit);
  
  // 'settors' simply replace local pointers to variables. dimensions much match
  // 'update' modifies the local copy (which may belong to someone else)
  void dbarts_setResponse(dbarts::BARTFit* fit, const double* newResponse);
  void dbarts_setOffset(dbarts::BARTFit* fit, const double* newOffset);
  
  // forceUpdate == true will cause the sampler to go through with the change even if it
  //   would leave the sampler in an invalid state, i.e. with a leaf having no observations;
  //   in that case, the empty leaf will be pruned;
  // forceUpdate == false rolls back invalid changes and can be used for rejection sampling
  // updateCutPoints == true takes the new predictor and uses the default rule to create
  //   a new set of cut points from it, and rebalances observations in nodes as necessary
  int dbarts_setPredictor(dbarts::BARTFit* fit, const double* newPredictor, int forceUpdate, int updateCutPoints);
  int dbarts_updatePredictor(dbarts::BARTFit* fit, const double* newPredictor, std::size_t column, int forceUpdate, int updateCutPoints);
  int dbarts_updatePredictors(dbarts::BARTFit* fit, const double* newPredictor, const std::size_t* columns, std::size_t numColumns, int forceUpdate, int updateCutPoints);
  
  void dbarts_setTestPredictor(dbarts::BARTFit* fit, const double* newTestPredictor, std::size_t numTestObservations);
  void dbarts_setTestOffset(dbarts::BARTFit* fit, const double* newTestOffset);
  void dbarts_setTestPredictorAndOffset(dbarts::BARTFit* fit, const double* newTestPredictor, const double* newTestOffset, std::size_t numTestObservations);
  
  void dbarts_updateTestPredictor(dbarts::BARTFit* fit, const double* newTestPredictor, std::size_t column);
  void dbarts_updateTestPredictors(dbarts::BARTFit* fit, const double* newTestPredictor, const std::size_t* columns, std::size_t numColumns);
    
  dbarts::CGMPrior* dbarts_createCGMPrior();
  dbarts::CGMPrior* dbarts_createCGMPriorFromOptions(double base, double power);
  void dbarts_destroyCGMPrior(dbarts::CGMPrior* prior);
  void dbarts_initializeCGMPriorFromOptions(dbarts::CGMPrior* prior, double base, double power);
  void dbarts_invalidateCGMPrior(dbarts::CGMPrior* prior);
  
  dbarts::NormalPrior* dbarts_createNormalPrior();
  dbarts::NormalPrior* dbarts_createNormalPriorFromOptions(const dbarts::Control* control, const dbarts::Model* model, double k);
  void dbarts_destroyNormalPrior(dbarts::NormalPrior* prior);
  void dbarts_initializeNormalPriorFromOptions(dbarts::NormalPrior* prior, const dbarts::Control* control, const dbarts::Model* model, double k);
  void dbarts_invalidateNormalPrior(dbarts::NormalPrior* prior);
  
  dbarts::NormalHyperprior* dbarts_createNormalHyperprior();
  dbarts::NormalHyperprior* dbarts_createNormalHyperpriorFromOptions(double scale);
  void dbarts_destroyNormalHyperprior(dbarts::NormalHyperprior* prior);
  void dbarts_initializeNormalHyperpriorFromOptions(dbarts::NormalHyperprior* prior, double scale);
  void dbarts_invalidateNormalHyperprior(dbarts::NormalHyperprior* prior);
  
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPrior();
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPriorFromOptions(double degreesOfFreedom, double quantile);
  void dbarts_destroyChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
  void dbarts_initializeChiSquaredPriorFromOptions(dbarts::ChiSquaredPrior* prior, double degreesOfFreedom, double quantile);
  void dbarts_invalidateChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
}

#endif // DBARTS_R_C_INTERFACE_HPP
