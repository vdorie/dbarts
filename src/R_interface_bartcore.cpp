#include "config.hpp"
#include "R_interface_bartcore.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include <external/Rinternals.h>
#include <R_ext/Random.h> // GetRNGstate, PutRNGstate

#include <external/random.h>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/types.hpp>

#include "R_interface_common.hpp"

#include "bartcore/bartcore.hpp"

using std::size_t;
using namespace dbarts;

namespace {

struct BartcoreHolder {
  std::unique_ptr<bartcore::SamplerBase> sampler;
  std::vector<ext_rng*> rngs;  // one per chain
  bool keepTrainingFits;

  ~BartcoreHolder() {
    for (size_t c = rngs.size(); c > 0; --c)
      if (rngs[c - 1] != NULL) ext_rng_destroy(rngs[c - 1]);
  }
};

// The external pointer's protection slot pins the vectors the sampler
// borrows, one fixed slot per borrowable so replacements do not accumulate.
enum {
  PROT_DATA = 0,
  PROT_RESPONSE,
  PROT_OFFSET,
  PROT_PREDICTORS,
  PROT_TEST_PREDICTORS,
  PROT_COUNT
};

void retain(SEXP ptrExpr, int slot, SEXP value) {
  SET_VECTOR_ELT(R_ExternalPtrProtected(ptrExpr), slot, value);
}

void holderFinalizer(SEXP ptrExpr) {
  BartcoreHolder* holder =
    static_cast<BartcoreHolder*>(R_ExternalPtrAddr(ptrExpr));
  if (holder == NULL) return;
  delete holder;
  R_ClearExternalPtr(ptrExpr);
}

BartcoreHolder& holderFromExpression(SEXP ptrExpr) {
  BartcoreHolder* holder =
    static_cast<BartcoreHolder*>(R_ExternalPtrAddr(ptrExpr));
  if (holder == NULL)
    Rf_error("bartcore function called on NULL external pointer");
  return *holder;
}

} // namespace

extern "C" {

// The external pointer's protection slot pins everything the sampler
// borrows: the data expression at creation, and any replacement vectors the
// setters install later.
SEXP bartcore_create(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr) {
  Control control;
  Data data;
  Model model;

  initializeControlFromExpression(control, controlExpr);
  initializeDataFromExpression(data, dataExpr);
  initializeModelFromExpression(model, modelExpr, control, data);

  // Rf_error longjmps past destructors, so collect the reason, clean up,
  // and error at the end.
  const char* errorMessage = NULL;

  for (size_t j = 0; j < data.numPredictors && errorMessage == NULL; ++j) {
    if (data.variableTypes[j] != ORDINAL)
      errorMessage = "bartcore supports only ordinal predictors so far";
  }
  if (errorMessage == NULL && !control.responseIsBinary &&
      model.sigmaSqPrior->isFixed)
    errorMessage = "bartcore does not support fixed sigma for continuous responses";
  if (errorMessage == NULL && data.testOffset != NULL)
    errorMessage = "bartcore does not support test offsets";

  double k = 2.0, sigmaDf = 3.0, sigmaRawScale = 1.0;
  bool updateK = !model.kPrior->isFixed;
  double kDf = 1.25, kScale = HUGE_VAL;
  if (errorMessage == NULL) {
    if (updateK) {
      const ChiHyperprior& kPrior(*static_cast<ChiHyperprior*>(model.kPrior));
      kDf = kPrior.degreesOfFreedom;
      kScale = kPrior.scale;
    } else {
      k = static_cast<FixedHyperprior*>(model.kPrior)->getK();
    }
    if (!control.responseIsBinary) {
      const ChiSquaredPrior& sigmaSqPrior(
        *static_cast<ChiSquaredPrior*>(model.sigmaSqPrior));
      sigmaDf = sigmaSqPrior.degreesOfFreedom;
      sigmaRawScale = sigmaSqPrior.scale;
    }
  }

  if (errorMessage != NULL) {
    invalidateModel(model);
    invalidateData(data);
    Rf_error("%s", errorMessage);
  }

  const CGMPrior& treePrior(*static_cast<CGMPrior*>(model.treePrior));

  bartcore::SamplerOptions options;
  options.numTrees = control.numTrees;
  options.k = k;
  options.nodeScale = model.nodeScale;
  options.base = treePrior.base;
  options.power = treePrior.power;
  options.birthOrDeathProbability = model.birthOrDeathProbability;
  options.swapProbability = model.swapProbability;
  options.changeProbability = model.changeProbability;
  options.birthProbability = model.birthProbability;
  options.maxNumCutsPerVariable = data.maxNumCuts; // copied during build
  options.useQuantiles = control.useQuantiles;
  options.splitProbabilities = treePrior.splitProbabilities; // copied by ctor
  options.updateK = updateK;
  options.kHyperprior.degreesOfFreedom = kDf;
  options.kHyperprior.scale = kScale;
  options.numChains = control.numChains;
  options.numThreads = control.numThreads;
  options.numThin = control.treeThinningRate;

  // A single chain draws through R's generator; several chains each get a
  // Mersenne twister seeded from R's stream, so results do not depend on the
  // thread count and worker threads never touch the R API.
  std::vector<ext_rng*> rngs(options.numChains, static_cast<ext_rng*>(NULL));
  bool rngFailed = false;
  if (options.numChains == 1) {
    rngs[0] = ext_rng_createDefault(true);
    rngFailed = rngs[0] == NULL;
  } else {
    GetRNGstate();
    for (size_t c = 0; c < options.numChains && !rngFailed; ++c) {
      rngs[c] = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      rngFailed = rngs[c] == NULL ||
        ext_rng_setSeed(rngs[c], static_cast<std::uint_least32_t>(
                                   unif_rand() * 4294967295.0)) != 0;
    }
    PutRNGstate();
  }
  if (rngFailed) {
    for (size_t c = rngs.size(); c > 0; --c)
      if (rngs[c - 1] != NULL) ext_rng_destroy(rngs[c - 1]);
    invalidateModel(model);
    invalidateData(data);
    Rf_error("could not allocate rng");
  }

  std::unique_ptr<bartcore::SamplerBase> sampler = bartcore::createClassicSampler(
    data.x, data.y, data.numObservations, data.numPredictors, data.weights,
    data.offset, control.responseIsBinary, data.sigmaEstimate, sigmaDf,
    sigmaRawScale, options, rngs.data());

  if (data.numTestObservations > 0)
    sampler->setTestPredictors(data.x_test, data.numTestObservations);

  invalidateModel(model);
  invalidateData(data);

  BartcoreHolder* holder = new BartcoreHolder{std::move(sampler),
                                              std::move(rngs),
                                              control.keepTrainingFits};

  SEXP protExpr = PROTECT(Rf_allocVector(VECSXP, PROT_COUNT));
  SET_VECTOR_ELT(protExpr, PROT_DATA, dataExpr);
  SEXP result = PROTECT(R_MakeExternalPtr(holder, R_NilValue, protExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  UNPROTECT(2);
  return result;
}

SEXP bartcore_run(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  size_t numBurnIn = static_cast<size_t>(Rf_asInteger(numBurnInExpr));
  size_t numSamples = static_cast<size_t>(Rf_asInteger(numSamplesExpr));

  size_t numObservations = sampler.numObservations();
  size_t numPredictors = sampler.numPredictors();
  size_t numTestObservations = sampler.numTestObservations();
  size_t numChains = sampler.numChains();
  int numSamplesInt = static_cast<int>(numSamples);
  int numChainsInt = static_cast<int>(numChains);

  if (numSamples == 0) {
    bartcore::Results empty;
    GetRNGstate();
    sampler.run(numBurnIn, 0, empty);
    PutRNGstate();
    return R_NilValue;
  }

  // several chains add a trailing chain dimension, as the classic engine's
  // results do
  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP sigmaExpr = numChains == 1
    ? PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples)))
    : PROTECT(Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  SEXP trainExpr;
  if (!holder.keepTrainingFits) {
    trainExpr = PROTECT(R_NilValue);
  } else if (numChains == 1) {
    trainExpr = PROTECT(Rf_allocMatrix(
      REALSXP, static_cast<int>(numObservations), numSamplesInt));
  } else {
    trainExpr = PROTECT(Rf_alloc3DArray(
      REALSXP, static_cast<int>(numObservations), numSamplesInt,
      numChainsInt));
  }
  SEXP testExpr;
  if (numTestObservations > 0) {
    testExpr = numChains == 1
      ? PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numTestObservations),
                               numSamplesInt))
      : PROTECT(Rf_alloc3DArray(REALSXP,
                                static_cast<int>(numTestObservations),
                                numSamplesInt, numChainsInt));
  } else {
    testExpr = PROTECT(R_NilValue);
  }
  SEXP varcountExpr = numChains == 1
    ? PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(numPredictors),
                             numSamplesInt))
    : PROTECT(Rf_alloc3DArray(INTSXP, static_cast<int>(numPredictors),
                              numSamplesInt, numChainsInt));
  SEXP kExpr;
  if (sampler.kIsSampled()) {
    kExpr = numChains == 1
      ? PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples)))
      : PROTECT(Rf_allocMatrix(REALSXP, numSamplesInt, numChainsInt));
  } else {
    kExpr = PROTECT(R_NilValue);
  }

  std::vector<std::uint32_t> variableCounts(numPredictors * numSamples *
                                            numChains);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = holder.keepTrainingFits ? REAL(trainExpr) : NULL;
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();
  results.k = sampler.kIsSampled() ? REAL(kExpr) : NULL;

  GetRNGstate();
  sampler.run(numBurnIn, numSamples, results);
  PutRNGstate();

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numSamples * numChains; ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, testExpr);
  SET_VECTOR_ELT(resultExpr, 3, varcountExpr);
  SET_VECTOR_ELT(resultExpr, 4, kExpr);

  // named as the classic engine's run results are, so the engines are
  // drop-in replacements for each other
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(7);
  return resultExpr;
}

SEXP bartcore_setOffset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  const double* offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  holder.sampler->setOffset(offset, Rf_asLogical(updateScaleExpr) == TRUE);
  retain(ptrExpr, PROT_OFFSET, offsetExpr);
  return R_NilValue;
}

SEXP bartcore_setResponse(SEXP ptrExpr, SEXP yExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  GetRNGstate(); // probit latent redraw
  holder.sampler->setResponse(REAL(yExpr));
  PutRNGstate();
  retain(ptrExpr, PROT_RESPONSE, yExpr);
  return R_NilValue;
}

SEXP bartcore_setSigma(SEXP ptrExpr, SEXP sigmaExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  holder.sampler->setSigma(Rf_asReal(sigmaExpr));
  return R_NilValue;
}

SEXP bartcore_setData(SEXP ptrExpr, SEXP dataExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  bartcore::SamplerBase& sampler(*holder.sampler);

  if (!Rf_inherits(dataExpr, "dbartsData"))
    Rf_error("'data' argument to bartcore_setData not of class 'dbartsData'");

  Data data;
  initializeDataFromExpression(data, dataExpr);

  // Rf_error longjmps past destructors, so collect the reason, clean up,
  // and error at the end.
  const char* errorMessage = NULL;
  if (data.numPredictors != sampler.numPredictors())
    errorMessage = "bartcore setData requires the same predictors";
  for (size_t j = 0; j < data.numPredictors && errorMessage == NULL; ++j) {
    if (data.variableTypes[j] != ORDINAL)
      errorMessage = "bartcore supports only ordinal predictors so far";
  }
  if (errorMessage == NULL && data.testOffset != NULL)
    errorMessage = "bartcore does not support test offsets";

  if (errorMessage != NULL) {
    invalidateData(data);
    Rf_error("%s", errorMessage);
  }

  sampler.setData(data.x, data.y, data.numObservations, data.weights,
                  data.offset, data.x_test, data.numTestObservations);

  invalidateData(data);

  // everything the sampler borrows now comes from the new data object
  retain(ptrExpr, PROT_DATA, dataExpr);
  retain(ptrExpr, PROT_RESPONSE, R_NilValue);
  retain(ptrExpr, PROT_OFFSET, R_NilValue);
  retain(ptrExpr, PROT_PREDICTORS, R_NilValue);
  retain(ptrExpr, PROT_TEST_PREDICTORS, R_NilValue);

  return R_NilValue;
}

SEXP bartcore_setTestPredictor(SEXP ptrExpr, SEXP xTestExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  SEXP dims = Rf_getAttrib(xTestExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[1]) != holder.sampler->numPredictors())
    Rf_error("bartcore_setTestPredictor requires a matrix with matching columns");
  holder.sampler->setTestPredictors(
    REAL(xTestExpr), static_cast<size_t>(INTEGER(dims)[0]));
  retain(ptrExpr, PROT_TEST_PREDICTORS, xTestExpr);
  return R_NilValue;
}

SEXP bartcore_isValidPointer(SEXP ptrExpr) {
  return Rf_ScalarLogical(R_ExternalPtrAddr(ptrExpr) != NULL ? TRUE : FALSE);
}

SEXP bartcore_getSigmas(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numChains = holder.sampler->numChains();
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numChains)));
  for (size_t c = 0; c < numChains; ++c)
    REAL(result)[c] = holder.sampler->sigma(c);
  UNPROTECT(1);
  return result;
}

// Predictor mutation. The sampler borrows a full replacement matrix
// (R/bartcore.R retains it on success); column and per-observation updates
// write in place into the matrix the sampler currently borrows, aliasing the
// R-side data like the classic engine does.

SEXP bartcore_setPredictor(SEXP ptrExpr, SEXP xExpr, SEXP forceUpdateExpr,
                           SEXP updateCutPointsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2 ||
      static_cast<size_t>(INTEGER(dims)[0]) != holder.sampler->numObservations() ||
      static_cast<size_t>(INTEGER(dims)[1]) != holder.sampler->numPredictors())
    Rf_error("bartcore_setPredictor requires a matrix with matching dimensions");

  bartcore::PredictorUpdateResult result = holder.sampler->setPredictor(
    REAL(xExpr), Rf_asLogical(forceUpdateExpr) == TRUE,
    Rf_asLogical(updateCutPointsExpr) == TRUE);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  if (result == bartcore::PredictorUpdateResult::accepted)
    retain(ptrExpr, PROT_PREDICTORS, xExpr);  // installed by pointer swap
  return Rf_ScalarLogical(
    result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
}

SEXP bartcore_updatePredictor(SEXP ptrExpr, SEXP xExpr, SEXP columnsExpr,
                              SEXP forceUpdateExpr, SEXP updateCutPointsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numObservations = holder.sampler->numObservations();
  size_t numPredictors = holder.sampler->numPredictors();

  size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
  if (numColumns == 0 ||
      static_cast<size_t>(Rf_xlength(xExpr)) != numObservations * numColumns)
    Rf_error("bartcore_updatePredictor requires numObservations values per column");

  std::vector<size_t> columns(numColumns);
  for (size_t k = 0; k < numColumns; ++k) {
    int column = INTEGER(columnsExpr)[k];
    if (column < 1 || static_cast<size_t>(column) > numPredictors)
      Rf_error("bartcore_updatePredictor column out of range");
    columns[k] = static_cast<size_t>(column - 1);
  }

  bartcore::PredictorUpdateResult result = holder.sampler->updatePredictor(
    REAL(xExpr), columns.data(), numColumns,
    Rf_asLogical(forceUpdateExpr) == TRUE,
    Rf_asLogical(updateCutPointsExpr) == TRUE);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  return Rf_ScalarLogical(
    result == bartcore::PredictorUpdateResult::accepted ? TRUE : FALSE);
}

SEXP bartcore_setCutPoints(SEXP ptrExpr, SEXP cutPointsExpr, SEXP columnsExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numPredictors = holder.sampler->numPredictors();

  size_t numColumns = static_cast<size_t>(Rf_xlength(columnsExpr));
  if (numColumns == 0 ||
      static_cast<size_t>(Rf_xlength(cutPointsExpr)) != numColumns)
    Rf_error("bartcore_setCutPoints requires one cut point vector per column");

  std::vector<const double*> cutPoints(numColumns);
  std::vector<std::uint32_t> numCutPoints(numColumns);
  std::vector<size_t> columns(numColumns);
  for (size_t k = 0; k < numColumns; ++k) {
    int column = INTEGER(columnsExpr)[k];
    if (column < 1 || static_cast<size_t>(column) > numPredictors)
      Rf_error("bartcore_setCutPoints column out of range");
    columns[k] = static_cast<size_t>(column - 1);

    SEXP cutsExpr = VECTOR_ELT(cutPointsExpr, static_cast<R_xlen_t>(k));
    R_xlen_t numCuts = Rf_xlength(cutsExpr);
    if (numCuts > 65535)  // codes must fit xint_t, including numCuts itself
      Rf_error("bartcore_setCutPoints cut point vector too long");
    const double* cuts = REAL(cutsExpr);
    for (R_xlen_t i = 1; i < numCuts; ++i)
      if (cuts[i] <= cuts[i - 1])
        Rf_error("bartcore_setCutPoints requires strictly increasing cut "
                 "points");
    cutPoints[k] = cuts;
    numCutPoints[k] = static_cast<std::uint32_t>(numCuts);
  }

  holder.sampler->setCutPoints(cutPoints.data(), numCutPoints.data(),
                               columns.data(), numColumns);
  return R_NilValue;
}

SEXP bartcore_updatePredictorPerObservation(SEXP ptrExpr, SEXP xExpr,
                                            SEXP columnExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  size_t numObservations = holder.sampler->numObservations();

  if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
    Rf_error("bartcore_updatePredictorPerObservation requires one value per "
             "observation");
  int column = Rf_asInteger(columnExpr);
  if (column < 1 ||
      static_cast<size_t>(column) > holder.sampler->numPredictors())
    Rf_error("bartcore_updatePredictorPerObservation column out of range");

  std::unique_ptr<bool[]> installed(new bool[numObservations]);

  GetRNGstate();  // scan-order permutation
  bool treesAreValid = holder.sampler->updatePredictorPerObservation(
    REAL(xExpr), static_cast<size_t>(column - 1), installed.get());
  PutRNGstate();

  // The sequential guard admits no empty leaves, so an invalid rebuild is an
  // internal invariant violation; fail loudly rather than leave the sampler
  // in an invalid state.
  if (!treesAreValid)
    Rf_error("bartcore updatePredictorPerObservation produced a tree with an "
             "empty leaf");

  SEXP result = PROTECT(
    Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
  for (size_t i = 0; i < numObservations; ++i)
    LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
  UNPROTECT(1);
  return result;
}

SEXP bartcore_updatePredictorPerObservationJointly(SEXP ptrsExpr, SEXP xExpr,
                                                   SEXP columnsExpr) {
  size_t numSamplers = static_cast<size_t>(Rf_xlength(ptrsExpr));
  if (numSamplers == 0 ||
      static_cast<size_t>(Rf_xlength(columnsExpr)) != numSamplers)
    Rf_error("bartcore_updatePredictorPerObservationJointly requires one "
             "column per sampler");

  std::vector<bartcore::SamplerBase*> samplers(numSamplers);
  std::vector<size_t> columns(numSamplers);
  for (size_t k = 0; k < numSamplers; ++k) {
    BartcoreHolder& holder(
      holderFromExpression(VECTOR_ELT(ptrsExpr, static_cast<R_xlen_t>(k))));
    samplers[k] = holder.sampler.get();
    int column = INTEGER(columnsExpr)[k];
    if (column < 1 ||
        static_cast<size_t>(column) > samplers[k]->numPredictors())
      Rf_error("bartcore_updatePredictorPerObservationJointly column out of "
               "range");
    columns[k] = static_cast<size_t>(column - 1);
  }

  size_t numObservations = samplers[0]->numObservations();
  for (size_t k = 1; k < numSamplers; ++k)
    if (samplers[k]->numObservations() != numObservations)
      Rf_error("bartcore_updatePredictorPerObservationJointly requires "
               "index-aligned samplers");
  if (static_cast<size_t>(Rf_xlength(xExpr)) != numObservations)
    Rf_error("bartcore_updatePredictorPerObservationJointly requires one "
             "value per observation");

  std::unique_ptr<bool[]> installed(new bool[numObservations]);

  GetRNGstate();  // scan-order permutation
  bool treesAreValid = bartcore::updatePredictorPerObservationJointly(
    samplers.data(), numSamplers, REAL(xExpr), columns.data(), installed.get());
  PutRNGstate();

  if (!treesAreValid)
    Rf_error("bartcore updatePredictorPerObservationJointly produced a tree "
             "with an empty leaf");

  SEXP result = PROTECT(
    Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(numObservations)));
  for (size_t i = 0; i < numObservations; ++i)
    LOGICAL(result)[i] = installed[i] ? TRUE : FALSE;
  UNPROTECT(1);
  return result;
}

SEXP bartcore_getLatents(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  if (holder.sampler->latents(0) == NULL) return R_NilValue;

  size_t numObservations = holder.sampler->numObservations();
  size_t numChains = holder.sampler->numChains();

  SEXP result;
  if (numChains == 1) {
    result =
      PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numObservations)));
  } else {
    result = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numObservations),
                                    static_cast<int>(numChains)));
  }
  for (size_t c = 0; c < numChains; ++c)
    std::memcpy(REAL(result) + c * numObservations,
                holder.sampler->latents(c), numObservations * sizeof(double));
  UNPROTECT(1);
  return result;
}

} // extern "C"
