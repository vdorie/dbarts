/* A minimal consumer of the flat C API (dbarts/dbarts.h), compiled by
 * test-capi.R with R CMD SHLIB against the installed package headers. It
 * resolves every entry point through R_GetCCallable, exactly as a
 * LinkingTo package would, and exposes .Call wrappers for the R-side
 * assertions. */

#include <dbarts/dbarts.h>

#include <string.h> /* memcpy */

#include <R_ext/Rdynload.h>

static int (*p_apiVersion)(void);
static dbarts_sampler* (*p_create)(SEXP, SEXP, SEXP, const char*);
static void (*p_destroy)(dbarts_sampler*);
static void (*p_run)(dbarts_sampler*, size_t, size_t, dbarts_results*);
static void (*p_sampleTreesFromPrior)(dbarts_sampler*);
static void (*p_setResponse)(dbarts_sampler*, const double*);
static void (*p_setOffset)(dbarts_sampler*, const double*, int);
static void (*p_setSigma)(dbarts_sampler*, double);
static int (*p_getLatents)(const dbarts_sampler*, double*);
static int (*p_setPredictor)(dbarts_sampler*, const double*, int, int);
static int (*p_updatePredictor)(dbarts_sampler*, const double*, const size_t*,
                                size_t, int, int);
static void (*p_setTestPredictors)(dbarts_sampler*, const double*, size_t);
static void (*p_predict)(dbarts_sampler*, const double*, size_t,
                         const double*, double*);
static void (*p_setTreeStorage)(dbarts_sampler*, int, size_t);
static SEXP (*p_getTrees)(dbarts_sampler*, const size_t*, size_t,
                          const size_t*, size_t, const size_t*, size_t, int);
static SEXP (*p_storeState)(dbarts_sampler*);
static void (*p_setState)(dbarts_sampler*, SEXP);
static size_t (*p_numObservations)(const dbarts_sampler*);
static size_t (*p_numPredictors)(const dbarts_sampler*);
static size_t (*p_numTestObservations)(const dbarts_sampler*);
static size_t (*p_numChains)(const dbarts_sampler*);
static size_t (*p_numTrees)(const dbarts_sampler*);
static size_t (*p_numSavedSamples)(const dbarts_sampler*);
static int (*p_kIsSampled)(const dbarts_sampler*);
static int (*p_usesDart)(const dbarts_sampler*);

static void initApi(void) {
  if (p_apiVersion != NULL) return;
#define LOOKUP(_T_, _P_, _N_) _P_ = (_T_) R_GetCCallable("dbarts", _N_)
  LOOKUP(int (*)(void), p_apiVersion, "dbarts_apiVersion");
  LOOKUP(dbarts_sampler* (*)(SEXP, SEXP, SEXP, const char*), p_create,
         "dbarts_sampler_create");
  LOOKUP(void (*)(dbarts_sampler*), p_destroy, "dbarts_sampler_destroy");
  LOOKUP(void (*)(dbarts_sampler*, size_t, size_t, dbarts_results*), p_run,
         "dbarts_sampler_run");
  LOOKUP(void (*)(dbarts_sampler*), p_sampleTreesFromPrior,
         "dbarts_sampler_sampleTreesFromPrior");
  LOOKUP(void (*)(dbarts_sampler*, const double*), p_setResponse,
         "dbarts_sampler_setResponse");
  LOOKUP(void (*)(dbarts_sampler*, const double*, int), p_setOffset,
         "dbarts_sampler_setOffset");
  LOOKUP(void (*)(dbarts_sampler*, double), p_setSigma,
         "dbarts_sampler_setSigma");
  LOOKUP(int (*)(const dbarts_sampler*, double*), p_getLatents,
         "dbarts_sampler_getLatents");
  LOOKUP(int (*)(dbarts_sampler*, const double*, int, int), p_setPredictor,
         "dbarts_sampler_setPredictor");
  LOOKUP(int (*)(dbarts_sampler*, const double*, const size_t*, size_t, int,
                 int),
         p_updatePredictor, "dbarts_sampler_updatePredictor");
  LOOKUP(void (*)(dbarts_sampler*, const double*, size_t),
         p_setTestPredictors, "dbarts_sampler_setTestPredictors");
  LOOKUP(void (*)(dbarts_sampler*, const double*, size_t, const double*,
                  double*),
         p_predict, "dbarts_sampler_predict");
  LOOKUP(void (*)(dbarts_sampler*, int, size_t), p_setTreeStorage,
         "dbarts_sampler_setTreeStorage");
  LOOKUP(SEXP (*)(dbarts_sampler*, const size_t*, size_t, const size_t*,
                  size_t, const size_t*, size_t, int),
         p_getTrees, "dbarts_sampler_getTrees");
  LOOKUP(SEXP (*)(dbarts_sampler*), p_storeState, "dbarts_sampler_storeState");
  LOOKUP(void (*)(dbarts_sampler*, SEXP), p_setState,
         "dbarts_sampler_setState");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numObservations,
         "dbarts_sampler_numObservations");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numPredictors,
         "dbarts_sampler_numPredictors");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numTestObservations,
         "dbarts_sampler_numTestObservations");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numChains,
         "dbarts_sampler_numChains");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numTrees,
         "dbarts_sampler_numTrees");
  LOOKUP(size_t (*)(const dbarts_sampler*), p_numSavedSamples,
         "dbarts_sampler_numSavedSamples");
  LOOKUP(int (*)(const dbarts_sampler*), p_kIsSampled,
         "dbarts_sampler_kIsSampled");
  LOOKUP(int (*)(const dbarts_sampler*), p_usesDart,
         "dbarts_sampler_usesDart");
#undef LOOKUP
}

static void samplerFinalizer(SEXP ptrExpr) {
  dbarts_sampler* sampler = (dbarts_sampler*) R_ExternalPtrAddr(ptrExpr);
  if (sampler == NULL) return;
  p_destroy(sampler);
  R_ClearExternalPtr(ptrExpr);
}

static dbarts_sampler* samplerFromExpr(SEXP ptrExpr) {
  dbarts_sampler* sampler = (dbarts_sampler*) R_ExternalPtrAddr(ptrExpr);
  if (sampler == NULL) Rf_error("consumer called on NULL sampler");
  return sampler;
}

SEXP capi_version(void) {
  initApi();
  return Rf_ScalarInteger(p_apiVersion());
}

SEXP capi_create(SEXP control, SEXP model, SEXP data, SEXP family) {
  initApi();
  const char* familyName =
    Rf_isNull(family) ? "" : CHAR(STRING_ELT(family, 0));
  dbarts_sampler* sampler = p_create(control, model, data, familyName);
  SEXP result = PROTECT(R_MakeExternalPtr(sampler, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, samplerFinalizer, FALSE);
  UNPROTECT(1);
  return result;
}

SEXP capi_dims(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 8));
  int* dims = INTEGER(result);
  dims[0] = (int) p_numObservations(sampler);
  dims[1] = (int) p_numPredictors(sampler);
  dims[2] = (int) p_numTestObservations(sampler);
  dims[3] = (int) p_numChains(sampler);
  dims[4] = (int) p_numTrees(sampler);
  dims[5] = (int) p_numSavedSamples(sampler);
  dims[6] = p_kIsSampled(sampler);
  dims[7] = p_usesDart(sampler);
  UNPROTECT(1);
  return result;
}

SEXP capi_sample_trees_from_prior(SEXP ptrExpr) {
  p_sampleTreesFromPrior(samplerFromExpr(ptrExpr));
  return R_NilValue;
}

/* run with caller-owned buffers; keepTrain/keepTest exercise the
 * null-means-skip contract */
SEXP capi_run(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr,
              SEXP keepTrainExpr, SEXP keepTestExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  int keepTrain = Rf_asLogical(keepTrainExpr) == TRUE;
  int keepTest = Rf_asLogical(keepTestExpr) == TRUE;

  size_t n = p_numObservations(sampler);
  size_t p = p_numPredictors(sampler);
  size_t nTest = p_numTestObservations(sampler);
  size_t chains = p_numChains(sampler);

  SEXP sigmaExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numSamples * chains)));
  SEXP trainExpr = keepTrain
    ? PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n * numSamples * chains)))
    : PROTECT(R_NilValue);
  SEXP testExpr = keepTest && nTest > 0
    ? PROTECT(
        Rf_allocVector(REALSXP, (R_xlen_t) (nTest * numSamples * chains)))
    : PROTECT(R_NilValue);
  SEXP varcountExpr = PROTECT(
    Rf_allocVector(INTSXP, (R_xlen_t) (p * numSamples * chains)));

  uint32_t* varcount =
    (uint32_t*) R_alloc(p * numSamples * chains, sizeof(uint32_t));

  dbarts_results results;
  results.sigma = REAL(sigmaExpr);
  results.train = keepTrain ? REAL(trainExpr) : NULL;
  results.test = keepTest && nTest > 0 ? REAL(testExpr) : NULL;
  results.varcount = varcount;
  results.k = NULL;
  results.varprobs = NULL;

  p_run(sampler, numBurnIn, numSamples, &results);

  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < p * numSamples * chains; ++i)
    varcountOut[i] = (int) varcount[i];

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, testExpr);
  SET_VECTOR_ELT(resultExpr, 3, varcountExpr);
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(6);
  return resultExpr;
}

SEXP capi_set_response(SEXP ptrExpr, SEXP yExpr) {
  p_setResponse(samplerFromExpr(ptrExpr), REAL(yExpr));
  return R_NilValue;
}

SEXP capi_set_offset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  p_setOffset(samplerFromExpr(ptrExpr),
              Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr),
              Rf_asLogical(updateScaleExpr) == TRUE);
  return R_NilValue;
}

SEXP capi_set_sigma(SEXP ptrExpr, SEXP sigmaExpr) {
  p_setSigma(samplerFromExpr(ptrExpr), Rf_asReal(sigmaExpr));
  return R_NilValue;
}

SEXP capi_get_latents(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t length = p_numObservations(sampler) * p_numChains(sampler);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  int haveLatents = p_getLatents(sampler, REAL(result));
  UNPROTECT(1);
  return haveLatents ? result : R_NilValue;
}

SEXP capi_set_predictor(SEXP ptrExpr, SEXP xExpr) {
  return Rf_ScalarLogical(
    p_setPredictor(samplerFromExpr(ptrExpr), REAL(xExpr), FALSE, TRUE));
}

SEXP capi_update_predictor(SEXP ptrExpr, SEXP xExpr, SEXP columnExpr) {
  size_t column = (size_t) Rf_asInteger(columnExpr); /* already 0-based */
  return Rf_ScalarLogical(p_updatePredictor(
    samplerFromExpr(ptrExpr), REAL(xExpr), &column, 1, FALSE, TRUE));
}

SEXP capi_set_test_predictors(SEXP ptrExpr, SEXP xTestExpr) {
  if (Rf_isNull(xTestExpr)) {
    p_setTestPredictors(samplerFromExpr(ptrExpr), NULL, 0);
  } else {
    size_t numRows = (size_t) INTEGER(Rf_getAttrib(xTestExpr, R_DimSymbol))[0];
    p_setTestPredictors(samplerFromExpr(ptrExpr), REAL(xTestExpr), numRows);
  }
  return R_NilValue;
}

SEXP capi_set_tree_storage(SEXP ptrExpr, SEXP keepTreesExpr,
                           SEXP numSamplesExpr) {
  p_setTreeStorage(samplerFromExpr(ptrExpr),
                   Rf_asLogical(keepTreesExpr) == TRUE,
                   (size_t) Rf_asInteger(numSamplesExpr));
  return R_NilValue;
}

SEXP capi_predict(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numRows = (size_t) INTEGER(Rf_getAttrib(xTestExpr, R_DimSymbol))[0];
  size_t saved = p_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length = numRows * numSamples * p_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  p_predict(sampler, REAL(xTestExpr), numRows,
            Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr), REAL(result));
  UNPROTECT(1);
  return result;
}

SEXP capi_get_trees(SEXP ptrExpr, SEXP useLiveTreesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  int useLiveTrees = Rf_asLogical(useLiveTreesExpr) == TRUE;

  size_t numChains = p_numChains(sampler);
  size_t numSaved = useLiveTrees ? 0 : p_numSavedSamples(sampler);
  size_t numTrees = p_numTrees(sampler);

  size_t* chainIndices = (size_t*) R_alloc(numChains, sizeof(size_t));
  for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
  size_t* sampleIndices =
    numSaved > 0 ? (size_t*) R_alloc(numSaved, sizeof(size_t)) : NULL;
  for (size_t i = 0; i < numSaved; ++i) sampleIndices[i] = i;
  size_t* treeIndices = (size_t*) R_alloc(numTrees, sizeof(size_t));
  for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;

  return p_getTrees(sampler, chainIndices, numChains, sampleIndices, numSaved,
                    treeIndices, numTrees, useLiveTrees);
}

SEXP capi_store_state(SEXP ptrExpr) {
  return p_storeState(samplerFromExpr(ptrExpr));
}

SEXP capi_set_state(SEXP ptrExpr, SEXP stateExpr) {
  p_setState(samplerFromExpr(ptrExpr), stateExpr);
  return R_NilValue;
}
