#ifndef DBARTS_R_C_INTERFACE_HPP
#define DBARTS_R_C_INTERFACE_HPP

#include <cstddef> // size_t

// imports Rinternals.h while doing the least to pollute namespaces

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#  define USE_FC_LEN_T
#endif

#if R_VERSION <= R_Version(3, 3, 1)
// Rinternals.h includes R_ext/Memory.h and R_ext/Utils.h which reference size_t
// Rinternals.h also references FILE from stdio.h
#  define NO_C_HEADERS
#  include <climits>
// #  include <cstddef>
#  include <cstdio>
using std::size_t;
using std::FILE;
#endif

// prevents R_ext/Error.h from mapping Rf_error -> error and Rf_warning ->
// warning
#ifndef R_NO_REMAP
#  define UNMAP_R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <Rinternals.h>

#undef NO_C_HEADERS
#ifdef UNMAP_R_NO_REMAP
#  undef R_NO_REMAP
#  undef UNMAP_R_NO_REMAP
#endif
#undef USE_FC_LEN_T

namespace dbarts {
  struct Control;
  struct Model;
  struct Data;
  struct State;

  struct CGMPrior;
  struct NormalPrior;
  struct ChiSquaredPrior;
  struct ChiHyperprior;
  struct FixedHyperprior;

  struct BARTFit;
  struct Results;

  struct FlattenedTrees;
} // namespace dbarts

// pair calls of create<->destroy, initialize<->invalidate
extern "C" {
/// Creates a control C++ object from a dbartsControl R structure.
dbarts::Control* dbarts_createControl(SEXP controlExpr);
void dbarts_destroyControl(dbarts::Control* control);
void dbarts_initializeControl(dbarts::Control* control, SEXP controlExpr);
// void dbarts_invalidateControl(dbarts::Control* control); // invalidation not
// necessary, allocates no memory
void dbarts_setControl(dbarts::BARTFit* fit, const dbarts::Control* control);

dbarts::Data* dbarts_createData(SEXP dataExpr);
void dbarts_destroyData(dbarts::Data* data);
void dbarts_initializeData(dbarts::Data* data, SEXP dataExpr);
void dbarts_invalidateData(dbarts::Data* data);

dbarts::Model* dbarts_createModel(
  SEXP modelExpr,
  const dbarts::Control* control,
  const dbarts::Data* data
);
void dbarts_destroyModel(dbarts::Model* model);
void dbarts_initializeModel(
  dbarts::Model* model,
  SEXP modelExpr,
  const dbarts::Control* control,
  const dbarts::Data* data
);
void dbarts_invalidateModel(dbarts::Model* model);

dbarts::BARTFit* dbarts_createFit(
  dbarts::Control* control,
  dbarts::Model* model,
  dbarts::Data* data
);
void dbarts_initializeFit(
  dbarts::BARTFit* fit,
  dbarts::Control* control,
  dbarts::Model* model,
  dbarts::Data* data
);
void dbarts_destroyFit(dbarts::BARTFit* fit);
void dbarts_invalidateFit(dbarts::BARTFit* fit);

SEXP dbarts_createStateExpression(const dbarts::BARTFit* fit);
void dbarts_initializeState(dbarts::BARTFit* fit, SEXP stateExpr);

void dbarts_setRNGState(
  dbarts::BARTFit* fit,
  const void* const* uniformState,
  const void* const* normalState
);

void dbarts_printInitialSummary(const dbarts::BARTFit* fit);
void dbarts_printTrees(
  const dbarts::BARTFit* fit,
  const std::size_t* chainIndices,
  std::size_t numChainIndices,
  const std::size_t* sampleIndices,
  std::size_t numSampleIndices,
  const std::size_t* treeIndices,
  std::size_t numTreeIndices
);
dbarts::FlattenedTrees* dbarts_getTrees(
  const dbarts::BARTFit* fit,
  const std::size_t* chainIndices,
  std::size_t numChainIndices,
  const std::size_t* sampleIndices,
  std::size_t numSampleIndices,
  const std::size_t* treeIndices,
  std::size_t numTreeIndices
);

dbarts::Results* dbarts_runSampler(dbarts::BARTFit* fit);
dbarts::Results* dbarts_runSamplerForIterations(
  dbarts::BARTFit* fit,
  std::size_t numBurnIn,
  std::size_t numSamples
);
void dbarts_runSamplerWithResults(
  dbarts::BARTFit* fit,
  std::size_t numBurnIn,
  dbarts::Results* results
);
dbarts::Results* dbarts_runSamplerMultithreadedForIterations(
  dbarts::BARTFit* fit,
  std::size_t numBurnIn,
  std::size_t numThreads,
  std::size_t numSamples
);
void dbarts_runSamplerMultithreadedWithResults(
  dbarts::BARTFit* fit,
  std::size_t numBurnIn,
  std::size_t numThreads,
  dbarts::Results* results
);

void dbarts_sampleTreesFromPrior(dbarts::BARTFit* fit);
void dbarts_sampleNodeParametersFromPrior(dbarts::BARTFit* fit);

void dbarts_predict(
  const dbarts::BARTFit* fit,
  const double* x_test,
  std::size_t numTestObservations,
  const double* testOffset,
  double* result
);
void dbarts_predictMultithreaded(
  const dbarts::BARTFit* fit,
  const double* x_test,
  std::size_t numTestObservations,
  const double* testOffset,
  std::size_t numThreads,
  double* result
);
// 'settors' simply replace local pointers to variables; dimensions must match.
// 'update' modifies the local copy (which may belong to someone else).
void dbarts_setResponse(dbarts::BARTFit* fit, const double* newResponse);
void dbarts_setOffset(
  dbarts::BARTFit* fit,
  const double* newOffset,
  bool updateScale
);

/// \param newSigma One sigma for each chain.
void dbarts_setSigma(dbarts::BARTFit* fit, const double* newSigma);

/// Replaces the predictor matrix, rolling back the change if it would leave a
/// leaf empty.
/// \param forceUpdate If \c true, go through with the change even when it
///   empties a leaf (the empty leaf is then pruned); if \c false, roll back
///   invalid changes, so this can be used for rejection sampling.
/// \param updateCutPoints If \c true, use the default rule to derive fresh cut
///   points from the new predictor and rebalance observations as needed.
/// \returns \c 1 if the change was kept, \c 0 if it was rolled back (only
///   possible when \p forceUpdate is \c false).
int dbarts_setPredictor(
  dbarts::BARTFit* fit,
  const double* newPredictor,
  int forceUpdate,
  int updateCutPoints
);
/// Replaces a single column of the predictor matrix; \p forceUpdate,
/// \p updateCutPoints, and the return value are as in \c dbarts_setPredictor.
int dbarts_updatePredictor(
  dbarts::BARTFit* fit,
  const double* newPredictor,
  std::size_t column,
  int forceUpdate,
  int updateCutPoints
);
/// Replaces the given columns of the predictor matrix; \p forceUpdate,
/// \p updateCutPoints, and the return value are as in \c dbarts_setPredictor.
int dbarts_updatePredictors(
  dbarts::BARTFit* fit,
  const double* newPredictor,
  const std::size_t* columns,
  std::size_t numColumns,
  int forceUpdate,
  int updateCutPoints
);
/// Per-observation rollback for a single column: installs each observation's
/// new value on its own, rejecting only those that would empty a leaf.
/// \param[out] installed Must point to space for numObservations ints; on
///   return each is \c 1 if the observation's new value was installed and \c 0
///   if it was rolled back (because installing it would empty a leaf).
void dbarts_updatePredictorPerObservation(
  dbarts::BARTFit* fit,
  const double* newColumn,
  std::size_t column,
  int* installed
);
/// Joint per-observation rollback across several index-aligned samplers: a
/// single sweep installs each observation in every sampler or in none, keeping
/// every sampler valid (no empty leaf) at once.
/// \param columns \p columns[j] is the column to replace in sampler \c j.
/// \param newColumn The shared length-numObservations proposal.
/// \param[out] installed Length numObservations; receives \c 1 / \c 0 per
///   observation, identical across all samplers.
void dbarts_updatePredictorPerObservationJointly(
  dbarts::BARTFit** fits,
  std::size_t numFits,
  const double* newColumn,
  const std::size_t* columns,
  int* installed
);

void dbarts_setTestPredictor(
  dbarts::BARTFit* fit,
  const double* newTestPredictor,
  std::size_t numTestObservations
);
void dbarts_setTestOffset(dbarts::BARTFit* fit, const double* newTestOffset);
void dbarts_setTestPredictorAndOffset(
  dbarts::BARTFit* fit,
  const double* newTestPredictor,
  const double* newTestOffset,
  std::size_t numTestObservations
);

void dbarts_updateTestPredictor(
  dbarts::BARTFit* fit,
  const double* newTestPredictor,
  std::size_t column
);
void dbarts_updateTestPredictors(
  dbarts::BARTFit* fit,
  const double* newTestPredictor,
  const std::size_t* columns,
  std::size_t numColumns
);

void dbarts_storeLatents(const dbarts::BARTFit* fit, double* target);

std::size_t dbarts_startThreads(dbarts::BARTFit* fit);
std::size_t
dbarts_startNumThreads(dbarts::BARTFit* fit, std::size_t numThreads);
void dbarts_stopThreads(dbarts::BARTFit* fit);

dbarts::CGMPrior* dbarts_createCGMPrior();
dbarts::CGMPrior* dbarts_createCGMPriorFromOptions(
  double base,
  double power,
  const double* splitProbabilities
);
void dbarts_destroyCGMPrior(dbarts::CGMPrior* prior);
void dbarts_initializeCGMPriorFromOptions(
  dbarts::CGMPrior* prior,
  double base,
  double power,
  const double* splitProbabilities
);
void dbarts_invalidateCGMPrior(dbarts::CGMPrior* prior);

dbarts::NormalPrior* dbarts_createNormalPrior();
dbarts::NormalPrior* dbarts_createNormalPriorFromOptions(
  const dbarts::Control* control,
  const dbarts::Model* model
);
void dbarts_destroyNormalPrior(dbarts::NormalPrior* prior);
void dbarts_initializeNormalPriorFromOptions(
  dbarts::NormalPrior* prior,
  const dbarts::Control* control,
  const dbarts::Model* model
);
void dbarts_invalidateNormalPrior(dbarts::NormalPrior* prior);

dbarts::ChiHyperprior* dbarts_createChiHyperprior();
dbarts::ChiHyperprior*
dbarts_createChiHyperpriorFromOptions(double degreesOfFreedom, double scale);
void dbarts_destroyChiHyperprior(dbarts::ChiHyperprior* prior);
void dbarts_initializeChiHyperpriorFromOptions(
  dbarts::ChiHyperprior* prior,
  double degreesOfFreedom,
  double scale
);
void dbarts_invalidateChiHyperprior(dbarts::ChiHyperprior* prior);

dbarts::FixedHyperprior* dbarts_createFixedHyperprior();
dbarts::FixedHyperprior* dbarts_createFixedHyperpriorFromOptions(double k);
void dbarts_destroyFixedHyperprior(dbarts::FixedHyperprior* prior);
void dbarts_initializeFixedHyperpriorFromOptions(
  dbarts::FixedHyperprior* prior,
  double k
);
void dbarts_invalidateFixedHyperprior(dbarts::FixedHyperprior* prior);

dbarts::ChiSquaredPrior* dbarts_createChiSquaredPrior();
dbarts::ChiSquaredPrior* dbarts_createChiSquaredPriorFromOptions(
  double degreesOfFreedom,
  double quantile
);
void dbarts_destroyChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
void dbarts_initializeChiSquaredPriorFromOptions(
  dbarts::ChiSquaredPrior* prior,
  double degreesOfFreedom,
  double quantile
);
void dbarts_invalidateChiSquaredPrior(dbarts::ChiSquaredPrior* prior);
}

#endif // DBARTS_R_C_INTERFACE_HPP
