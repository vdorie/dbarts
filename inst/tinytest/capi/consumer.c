/* A minimal consumer of the flat C API (dbarts/dbarts.h), compiled by
 * test-capi.R with R CMD SHLIB against the installed package headers. It drives
 * the entry points through the header's DBARTS_USE_STUBS stubs - the supported
 * LinkingTo path, where each dbarts_sampler_* call binds to a cached
 * R_GetCCallable pointer generated inside dbarts.h - and exposes .Call wrappers
 * for the R-side assertions. One entry point is still resolved by hand as a
 * deliberate canary (see p_apiVersion_raw). */

#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

#include <string.h> /* memcpy */

#include <R_ext/Rdynload.h> /* R_GetCCallable, for the raw canary below */

/* Deliberate canary: dbarts_apiVersion is ALSO resolved the old way, by hand
 * through R_GetCCallable with a hand-written cast - the un-stubbed per-symbol
 * path a consumer that declines DBARTS_USE_STUBS (or a diagnostic tool) still
 * relies on. Everything else goes through the stubs, so this one raw path
 * guards that plain R_RegisterCCallable registration keeps working on its own. */
static int (*p_apiVersion_raw)(void);

static void initCanary(void) {
  if (p_apiVersion_raw == NULL)
    p_apiVersion_raw =
      (int (*)(void)) R_GetCCallable("dbarts", "dbarts_apiVersion");
}

static void samplerFinalizer(SEXP ptrExpr) {
  dbarts_sampler* sampler = (dbarts_sampler*) R_ExternalPtrAddr(ptrExpr);
  if (sampler == NULL) return;
  dbarts_sampler_destroy(sampler);
  R_ClearExternalPtr(ptrExpr);
}

static dbarts_sampler* samplerFromExpr(SEXP ptrExpr) {
  dbarts_sampler* sampler = (dbarts_sampler*) R_ExternalPtrAddr(ptrExpr);
  if (sampler == NULL) Rf_error("consumer called on NULL sampler");
  return sampler;
}

/* the packed version, through the raw canary path */
SEXP capi_version(void) {
  initCanary();
  return Rf_ScalarInteger(p_apiVersion_raw());
}

/* the packed integer plus the two components, through the stubs; the R side
 * checks all three agree with the header macros */
SEXP capi_versions(void) {
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 3));
  int* v = INTEGER(result);
  v[0] = dbarts_apiVersion();
  v[1] = dbarts_apiMajorVersion();
  v[2] = dbarts_apiMinorVersion();
  UNPROTECT(1);
  return result;
}

SEXP capi_create(SEXP control, SEXP model, SEXP data, SEXP family) {
  const char* familyName =
    Rf_isNull(family) ? "" : CHAR(STRING_ELT(family, 0));
  dbarts_sampler* sampler =
    dbarts_sampler_create(control, model, data, familyName);
  SEXP result = PROTECT(R_MakeExternalPtr(sampler, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, samplerFinalizer, FALSE);
  UNPROTECT(1);
  return result;
}

SEXP capi_dims(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 8));
  int* dims = INTEGER(result);
  dims[0] = (int) dbarts_sampler_numObservations(sampler);
  dims[1] = (int) dbarts_sampler_numPredictors(sampler);
  dims[2] = (int) dbarts_sampler_numTestObservations(sampler);
  dims[3] = (int) dbarts_sampler_numChains(sampler);
  dims[4] = (int) dbarts_sampler_numTrees(sampler);
  dims[5] = (int) dbarts_sampler_numSavedSamples(sampler);
  dims[6] = dbarts_sampler_kIsSampled(sampler);
  dims[7] = dbarts_sampler_usesDart(sampler);
  UNPROTECT(1);
  return result;
}

SEXP capi_sample_trees_from_prior(SEXP ptrExpr) {
  dbarts_sampler_sampleTreesFromPrior(samplerFromExpr(ptrExpr));
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

  size_t n = dbarts_sampler_numObservations(sampler);
  size_t p = dbarts_sampler_numPredictors(sampler);
  size_t nTest = dbarts_sampler_numTestObservations(sampler);
  size_t chains = dbarts_sampler_numChains(sampler);

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

  dbarts_results results = {0};
  results.structSize = sizeof results;
  results.sigma = REAL(sigmaExpr);
  results.train = keepTrain ? REAL(trainExpr) : NULL;
  results.test = keepTest && nTest > 0 ? REAL(testExpr) : NULL;
  results.varcount = varcount;
  results.k = NULL;
  results.varprobs = NULL;

  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);

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

/* the write-guard canary: simulate an OLD, smaller caller by pinning
 * structSize to the offset of `test`, so only sigma and train are present-by-
 * size. Every field past that boundary is set to a poisoned pointer; the
 * gaussian sampler would produce varcount (and test), so a guard that wrote
 * through those slots by size-blindly trusting the pointers would dereference
 * the poison and crash. Returns TRUE iff the run completed (guard skipped the
 * out-of-bounds fields) with sigma still filled. */
SEXP capi_run_guard(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  size_t chains = dbarts_sampler_numChains(sampler);

  double* sigma = (double*) R_alloc(numSamples * chains, sizeof(double));
  double* train = (double*) R_alloc(n * numSamples * chains, sizeof(double));
  double* poison = (double*) (uintptr_t) 0x1;

  dbarts_results results = {0};
  results.structSize = offsetof(dbarts_results, test);
  results.sigma = sigma;
  results.train = train;
  results.test = poison;
  results.varcount = (uint32_t*) poison;
  results.k = poison;
  results.varprobs = poison;
  results.tau = poison;
  results.groupEffects = poison;

  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);

  int ok = 1;
  for (size_t i = 0; i < numSamples * chains; ++i)
    if (!(sigma[i] > 0.0) || sigma[i] != sigma[i]) ok = 0;
  return Rf_ScalarLogical(ok);
}

/* per-sweep callback state: sets sigma[sweepIndex] when sigmas is non-null,
 * counts invocations, and returns 0 (stop) once the count reaches stopAt */
typedef struct {
  const double* sigmas;
  size_t numSigmas;
  int stopAt; /* < 0 disables early stop */
  int count;
} callbackState;

static int sweepCallback(void* userData, dbarts_sampler* sampler,
                         size_t chainIndex, size_t sweepIndex, int isBurnIn) {
  callbackState* state = (callbackState*) userData;
  (void) chainIndex;
  (void) isBurnIn;
  ++state->count;
  if (state->sigmas != NULL && sweepIndex < state->numSigmas)
    dbarts_sampler_setSigma(sampler, state->sigmas[sweepIndex]);
  if (state->stopAt >= 0 && state->count >= state->stopAt) return 0;
  return 1;
}

/* registers sweepCallback for one run then clears it, so the borrowed sigmas
 * stay live throughout; returns sigma/train/varcount plus the invocation
 * count. sigmasExpr may be NULL; stopAt < 0 disables early stop. */
SEXP capi_run_with_callback(SEXP ptrExpr, SEXP numBurnInExpr,
                            SEXP numSamplesExpr, SEXP sigmasExpr,
                            SEXP stopAtExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  size_t p = dbarts_sampler_numPredictors(sampler);
  size_t chains = dbarts_sampler_numChains(sampler);

  callbackState state;
  state.sigmas = Rf_isNull(sigmasExpr) ? NULL : REAL(sigmasExpr);
  state.numSigmas =
    Rf_isNull(sigmasExpr) ? 0 : (size_t) Rf_xlength(sigmasExpr);
  state.stopAt = Rf_asInteger(stopAtExpr);
  state.count = 0;

  SEXP sigmaExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numSamples * chains)));
  SEXP trainExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (n * numSamples * chains)));
  uint32_t* varcount =
    (uint32_t*) R_alloc(p * numSamples * chains, sizeof(uint32_t));

  dbarts_results results = {0};
  results.structSize = sizeof results;
  results.sigma = REAL(sigmaExpr);
  results.train = REAL(trainExpr);
  results.test = NULL;
  results.varcount = varcount;
  results.k = NULL;
  results.varprobs = NULL;
  results.tau = NULL;
  results.groupEffects = NULL;

  dbarts_sampler_setCallback(sampler, sweepCallback, &state);
  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);
  dbarts_sampler_setCallback(sampler, NULL, NULL);

  SEXP varcountExpr = PROTECT(
    Rf_allocVector(INTSXP, (R_xlen_t) (p * numSamples * chains)));
  int* varcountOut = INTEGER(varcountExpr);
  for (size_t i = 0; i < p * numSamples * chains; ++i)
    varcountOut[i] = (int) varcount[i];

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, varcountExpr);
  SET_VECTOR_ELT(resultExpr, 3, Rf_ScalarInteger(state.count));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("varcount"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("count"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(5);
  return resultExpr;
}

/* a grouped run: fills and returns tau (numSamples x chains) and groupEffects
 * (numGroups x numSamples x chains), the buffers dbarts_results now carries */
SEXP capi_run_grouped(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr,
                      SEXP numGroupsExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t numGroups = (size_t) Rf_asInteger(numGroupsExpr);
  size_t chains = dbarts_sampler_numChains(sampler);

  SEXP sigmaExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numSamples * chains)));
  SEXP tauExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numSamples * chains)));
  SEXP ranefExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numGroups * numSamples * chains)));

  dbarts_results results = {0};
  results.structSize = sizeof results;
  results.sigma = REAL(sigmaExpr);
  results.train = NULL;
  results.test = NULL;
  results.varcount = NULL;
  results.k = NULL;
  results.varprobs = NULL;
  results.tau = REAL(tauExpr);
  results.groupEffects = REAL(ranefExpr);

  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, tauExpr);
  SET_VECTOR_ELT(resultExpr, 2, ranefExpr);
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("tau"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("ranef"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(5);
  return resultExpr;
}

/* runs with the logLikelihood channel set alongside sigma and train, and
 * returns all three, so the R side can check the per-draw log-likelihood
 * against a density recomputed on the same sigma/train draws */
SEXP capi_run_loglik(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  size_t chains = dbarts_sampler_numChains(sampler);

  SEXP sigmaExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numSamples * chains)));
  SEXP trainExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (n * numSamples * chains)));
  SEXP loglikExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (n * numSamples * chains)));

  dbarts_results results = {0};
  results.structSize = sizeof results;
  results.sigma = REAL(sigmaExpr);
  results.train = REAL(trainExpr);
  results.logLikelihood = REAL(loglikExpr);

  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);

  SEXP resultExpr = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(resultExpr, 0, sigmaExpr);
  SET_VECTOR_ELT(resultExpr, 1, trainExpr);
  SET_VECTOR_ELT(resultExpr, 2, loglikExpr);
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("loglik"));
  Rf_setAttrib(resultExpr, R_NamesSymbol, namesExpr);

  UNPROTECT(5);
  return resultExpr;
}

SEXP capi_sample_node_parameters_from_prior(SEXP ptrExpr) {
  dbarts_sampler_sampleNodeParametersFromPrior(samplerFromExpr(ptrExpr));
  return R_NilValue;
}

SEXP capi_set_response(SEXP ptrExpr, SEXP yExpr) {
  dbarts_sampler_setResponse(samplerFromExpr(ptrExpr), REAL(yExpr));
  return R_NilValue;
}

SEXP capi_set_weights(SEXP ptrExpr, SEXP weightsExpr) {
  dbarts_sampler_setWeights(samplerFromExpr(ptrExpr), REAL(weightsExpr));
  return R_NilValue;
}

SEXP capi_set_test_offset(SEXP ptrExpr, SEXP offsetExpr) {
  dbarts_sampler_setTestOffset(samplerFromExpr(ptrExpr),
                               Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr));
  return R_NilValue;
}

/* prints the first tree of the first chain, exercising the entry point;
 * the R side captures the console output */
SEXP capi_print_trees(SEXP ptrExpr) {
  size_t chainIndex = 0, treeIndex = 0;
  dbarts_sampler_printTrees(samplerFromExpr(ptrExpr), &chainIndex, 1, NULL, 0,
                            &treeIndex, 1);
  return R_NilValue;
}

SEXP capi_set_run_controls(SEXP ptrExpr, SEXP numThreadsExpr,
                           SEXP numThinExpr, SEXP verboseExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  dbarts_sampler_setNumThreads(sampler, (size_t) Rf_asInteger(numThreadsExpr));
  dbarts_sampler_setNumThin(sampler, (size_t) Rf_asInteger(numThinExpr));
  dbarts_sampler_setVerbose(sampler, Rf_asLogical(verboseExpr) == TRUE, 100);
  return R_NilValue;
}

SEXP capi_set_offset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  dbarts_sampler_setOffset(samplerFromExpr(ptrExpr),
                           Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr),
                           Rf_asLogical(updateScaleExpr) == TRUE);
  return R_NilValue;
}

SEXP capi_set_sigma(SEXP ptrExpr, SEXP sigmaExpr) {
  dbarts_sampler_setSigma(samplerFromExpr(ptrExpr), Rf_asReal(sigmaExpr));
  return R_NilValue;
}

SEXP capi_get_latents(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t length =
    dbarts_sampler_numObservations(sampler) * dbarts_sampler_numChains(sampler);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  int haveLatents = dbarts_sampler_getLatents(sampler, REAL(result));
  UNPROTECT(1);
  return haveLatents ? result : R_NilValue;
}

SEXP capi_set_predictor(SEXP ptrExpr, SEXP xExpr) {
  return Rf_ScalarLogical(
    dbarts_sampler_setPredictor(samplerFromExpr(ptrExpr), REAL(xExpr), FALSE,
                                TRUE));
}

SEXP capi_update_predictor(SEXP ptrExpr, SEXP xExpr, SEXP columnExpr) {
  size_t column = (size_t) Rf_asInteger(columnExpr); /* already 0-based */
  return Rf_ScalarLogical(dbarts_sampler_updatePredictor(
    samplerFromExpr(ptrExpr), REAL(xExpr), &column, 1, FALSE, TRUE));
}

SEXP capi_set_test_predictors(SEXP ptrExpr, SEXP xTestExpr) {
  if (Rf_isNull(xTestExpr)) {
    dbarts_sampler_setTestPredictors(samplerFromExpr(ptrExpr), NULL, 0);
  } else {
    size_t numRows = (size_t) INTEGER(Rf_getAttrib(xTestExpr, R_DimSymbol))[0];
    dbarts_sampler_setTestPredictors(samplerFromExpr(ptrExpr), REAL(xTestExpr),
                                     numRows);
  }
  return R_NilValue;
}

SEXP capi_set_tree_storage(SEXP ptrExpr, SEXP keepTreesExpr,
                           SEXP numSamplesExpr) {
  dbarts_sampler_setTreeStorage(samplerFromExpr(ptrExpr),
                                Rf_asLogical(keepTreesExpr) == TRUE,
                                (size_t) Rf_asInteger(numSamplesExpr));
  return R_NilValue;
}

SEXP capi_predict(SEXP ptrExpr, SEXP xTestExpr, SEXP offsetExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numRows = (size_t) INTEGER(Rf_getAttrib(xTestExpr, R_DimSymbol))[0];
  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length = numRows * numSamples * dbarts_sampler_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  dbarts_sampler_predict(sampler, REAL(xTestExpr), numRows,
                         Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr),
                         REAL(result));
  UNPROTECT(1);
  return result;
}

/* sampleNumsExpr null reads every saved sample; otherwise its 1-based indices
 * select the saved samples, so the caller can compare the all-samples table
 * against a per-sample gather (the consistency stan4bart's extract relies on). */
SEXP capi_get_trees(SEXP ptrExpr, SEXP useLiveTreesExpr, SEXP sampleNumsExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  int useLiveTrees = Rf_asLogical(useLiveTreesExpr) == TRUE;

  size_t numChains = dbarts_sampler_numChains(sampler);
  size_t numSaved = useLiveTrees ? 0 : dbarts_sampler_numSavedSamples(sampler);
  size_t numTrees = dbarts_sampler_numTrees(sampler);

  size_t* chainIndices = (size_t*) R_alloc(numChains, sizeof(size_t));
  for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
  size_t* treeIndices = (size_t*) R_alloc(numTrees, sizeof(size_t));
  for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;

  size_t numSampleIndices;
  size_t* sampleIndices;
  if (useLiveTrees || Rf_isNull(sampleNumsExpr)) {
    numSampleIndices = numSaved;
    sampleIndices =
      numSaved > 0 ? (size_t*) R_alloc(numSaved, sizeof(size_t)) : NULL;
    for (size_t i = 0; i < numSaved; ++i) sampleIndices[i] = i;
  } else {
    numSampleIndices = (size_t) Rf_xlength(sampleNumsExpr);
    sampleIndices = (size_t*) R_alloc(numSampleIndices, sizeof(size_t));
    for (size_t i = 0; i < numSampleIndices; ++i)
      sampleIndices[i] = (size_t) INTEGER(sampleNumsExpr)[i] - 1;
  }

  return dbarts_sampler_getTrees(sampler, chainIndices, numChains,
                                 sampleIndices, numSampleIndices, treeIndices,
                                 numTrees, useLiveTrees);
}

SEXP capi_store_state(SEXP ptrExpr) {
  return dbarts_sampler_storeState(samplerFromExpr(ptrExpr));
}

SEXP capi_set_state(SEXP ptrExpr, SEXP stateExpr) {
  dbarts_sampler_setState(samplerFromExpr(ptrExpr), stateExpr);
  return R_NilValue;
}
