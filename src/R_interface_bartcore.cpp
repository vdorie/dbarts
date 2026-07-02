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
  ext_rng* rng;

  ~BartcoreHolder() {
    if (rng != NULL) ext_rng_destroy(rng);
  }
};

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

// The data expression is stashed in the external pointer's protection slot,
// as the sampler borrows its vectors. Vectors passed to setters later are
// the R caller's to keep alive (R/bartcore.R retains them).
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
    if (data.maxNumCuts[j] != data.maxNumCuts[0])
      errorMessage = "bartcore requires a single n.cuts for all predictors";
  }
  if (errorMessage == NULL && !model.kPrior->isFixed)
    errorMessage = "bartcore supports only fixed k so far";
  if (errorMessage == NULL && !control.responseIsBinary &&
      model.sigmaSqPrior->isFixed)
    errorMessage = "bartcore does not support fixed sigma for continuous responses";

  double k = 0.0, sigmaDf = 3.0, sigmaRawScale = 1.0;
  if (errorMessage == NULL) {
    k = static_cast<FixedHyperprior*>(model.kPrior)->getK();
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
  options.maxNumCuts = data.maxNumCuts[0];
  options.splitProbabilities = treePrior.splitProbabilities; // copied by ctor

  ext_rng* rng = ext_rng_createDefault(true);
  if (rng == NULL) {
    invalidateModel(model);
    invalidateData(data);
    Rf_error("could not allocate rng");
  }

  std::unique_ptr<bartcore::SamplerBase> sampler = bartcore::createClassicSampler(
    data.x, data.y, data.numObservations, data.numPredictors, data.weights,
    data.offset, control.responseIsBinary, data.sigmaEstimate, sigmaDf,
    sigmaRawScale, options, rng);

  if (data.numTestObservations > 0)
    sampler->setTestPredictors(data.x_test, data.numTestObservations);

  invalidateModel(model);
  invalidateData(data);

  BartcoreHolder* holder = new BartcoreHolder{std::move(sampler), rng};

  SEXP result = PROTECT(R_MakeExternalPtr(holder, R_NilValue, dataExpr));
  R_RegisterCFinalizerEx(result, holderFinalizer, static_cast<Rboolean>(FALSE));

  UNPROTECT(1);
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

  if (numSamples == 0) {
    bartcore::Results empty;
    GetRNGstate();
    sampler.run(numBurnIn, 0, empty);
    PutRNGstate();
    return R_NilValue;
  }

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP sigmaExpr =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numSamples)));
  SEXP trainExpr = PROTECT(Rf_allocMatrix(
    REALSXP, static_cast<int>(numObservations), static_cast<int>(numSamples)));
  SEXP testExpr = numTestObservations > 0
    ? PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(numTestObservations),
                             static_cast<int>(numSamples)))
    : PROTECT(R_NilValue);
  SEXP varcountExpr = PROTECT(Rf_allocMatrix(
    INTSXP, static_cast<int>(numPredictors), static_cast<int>(numSamples)));

  std::vector<std::uint32_t> variableCounts(numPredictors * numSamples);

  bartcore::Results results;
  results.sigma = REAL(sigmaExpr);
  results.trainingFits = REAL(trainExpr);
  results.testFits = numTestObservations > 0 ? REAL(testExpr) : NULL;
  results.variableCounts = variableCounts.data();

  GetRNGstate();
  sampler.run(numBurnIn, numSamples, results);
  PutRNGstate();

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < numPredictors * numSamples; ++i)
    varcountOut[i] = static_cast<int>(variableCounts[i]);

  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, testExpr);
  SET_VECTOR_ELT(resultExpr, 3, varcountExpr);

  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("yhat.train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("yhat.test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(6);
  return resultExpr;
}

SEXP bartcore_setOffset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  const double* offset = Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr);
  holder.sampler->setOffset(offset, Rf_asLogical(updateScaleExpr) == TRUE);
  return R_NilValue;
}

SEXP bartcore_setResponse(SEXP ptrExpr, SEXP yExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  GetRNGstate(); // probit latent redraw
  holder.sampler->setResponse(REAL(yExpr));
  PutRNGstate();
  return R_NilValue;
}

SEXP bartcore_setSigma(SEXP ptrExpr, SEXP sigmaExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  holder.sampler->setSigma(Rf_asReal(sigmaExpr));
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
  return R_NilValue;
}

SEXP bartcore_getLatents(SEXP ptrExpr) {
  BartcoreHolder& holder(holderFromExpression(ptrExpr));
  const double* latents = holder.sampler->latents();
  if (latents == NULL) return R_NilValue;

  size_t numObservations = holder.sampler->numObservations();
  SEXP result =
    PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(numObservations)));
  std::memcpy(REAL(result), latents, numObservations * sizeof(double));
  UNPROTECT(1);
  return result;
}

} // extern "C"
