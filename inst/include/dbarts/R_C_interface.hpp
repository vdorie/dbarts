#ifndef DBARTS_R_C_INTERFACE_HPP
#define DBARTS_R_C_INTERFACE_HPP

#include <cstddef> // size_t

// imports Rinternals.h while doing the least to pollute namespaces

#include <Rversion.h>

#if R_VERSION <= R_Version(3,3,1)
// Rinternals.h includes R_ext/Memory.h and R_ext/Utils.h which reference size_t
// Rinternals.h also references FILE from stdio.h
#  define NO_C_HEADERS
#  include <climits>
// #  include <cstddef>
#  include <cstdio>
using std::size_t;
using std::FILE;
#endif

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning -> warning
#define R_NO_REMAP
#include <Rinternals.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP

namespace dbarts {
  struct Control;
  struct Model;
  struct Data;
  
  struct CGMPrior;
  struct NormalPrior;
  struct ChiSquaredPrior;
  struct ChiHyperprior;
  
  struct BARTFit;
  struct Results;
}

// pair calls of create<->destroy, initialize<->invalidate
extern "C" {
  // creates a control C++ object from a dbartsControl R structure
  dbarts::Control* dbarts_createControl(SEXP controlExpr);
  void dbarts_destroyControl(dbarts::Control* control);
  void dbarts_initializeControl(dbarts::Control* control, SEXP controlExpr);
  // void dbarts_invalidateControl(dbarts::Control* control); // invalidation not necessary, owns no memory
  
  dbarts::Data* dbarts_createData(SEXP dataExpr);
  void dbarts_destroyData(dbarts::Data* data);
  void dbarts_initializeData(dbarts::Data* data, SEXP dataExpr);
  void dbarts_invalidateData(dbarts::Data* data);
  
  dbarts::Model* dbarts_createModel(SEXP modelExpr, dbarts::Control* control);
  void dbarts_destroyModel(dbarts::Model* model);
  void dbarts_initializeModel(dbarts::Model* model, SEXP modelExpr, const dbarts::Control* control);
  void dbarts_invalidateModel(dbarts::Model* model);
  
  dbarts::BARTFit* dbarts_createFit(dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_initializeFit(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
  void dbarts_destroyFit(dbarts::BARTFit* fit);
  void dbarts_invalidateFit(dbarts::BARTFit* fit);
  
  void dbarts_setRNGState(dbarts::BARTFit* fit, const void* const* uniformState, const void* const* normalState);
  
  void dbarts_printInitialSummary(const dbarts::BARTFit* fit);
  dbarts::Results* dbarts_runSampler(dbarts::BARTFit* fit);
  dbarts::Results* dbarts_runSamplerForIterations(dbarts::BARTFit* fit, std::size_t numBurnIn, std::size_t numSamples);
  void dbarts_runSamplerWithResults(dbarts::BARTFit* fit, std::size_t numBurnIn, dbarts::Results* results);
  void dbarts_sampleTreesFromPrior(dbarts::BARTFit* fit);
  void dbarts_sampleNodeParametersFromPrior(dbarts::BARTFit* fit);
  
  // 'settors' simply replace local pointers to variables. dimensions much match
  // 'update' modifies the local copy (which may belong to someone else)
  void dbarts_setResponse(dbarts::BARTFit* fit, const double* newResponse);
  void dbarts_setOffset(dbarts::BARTFit* fit, const double* newOffset, bool updateScale);
  
  // one sigma for each chain
  void dbarts_setSigma(dbarts::BARTFit* fit, const double* newSigma);
  
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
  
  void dbarts_storeLatents(const dbarts::BARTFit* fit, double* target);
    
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
  
  dbarts::ChiHyperprior* dbarts_createChiHyperprior();
  dbarts::ChiHyperprior* dbarts_createChiHyperpriorFromOptions(double degreesOfFreedom, double scale);
  void dbarts_destroyChiHyperprior(dbarts::ChiHyperprior* prior);
  void dbarts_initializeChiHyperpriorFromOptions(dbarts::ChiHyperprior* prior, double degreesOfFreedom, double scale);
  void dbarts_invalidateChiHyperprior(dbarts::ChiHyperprior* prior);
  
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPrior();
  dbarts::ChiSquaredPrior* dbarts_createChiSquaredPriorFromOptions(double degreesOfFreedom, double quantile);
  void dbarts_destroyChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
  void dbarts_initializeChiSquaredPriorFromOptions(dbarts::ChiSquaredPrior* prior, double degreesOfFreedom, double quantile);
  void dbarts_invalidateChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
}

#endif // DBARTS_R_C_INTERFACE_HPP
