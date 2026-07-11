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
static void (*p_sampleNodeParametersFromPrior)(dbarts_sampler*);
static void (*p_setResponse)(dbarts_sampler*, const double*);
static void (*p_setOffset)(dbarts_sampler*, const double*, int);
static void (*p_setWeights)(dbarts_sampler*, const double*);
static void (*p_setSigma)(dbarts_sampler*, double);
static void (*p_setCallback)(dbarts_sampler*, dbarts_sampler_callback, void*);
static void (*p_setTestOffset)(dbarts_sampler*, const double*);
static void (*p_printTrees)(dbarts_sampler*, const size_t*, size_t,
                            const size_t*, size_t, const size_t*, size_t);
static void (*p_setNumThreads)(dbarts_sampler*, size_t);
static void (*p_setNumThin)(dbarts_sampler*, size_t);
static void (*p_setVerbose)(dbarts_sampler*, int, uint32_t);
static int (*p_getLatents)(const dbarts_sampler*, double*);
static int (*p_setPredictor)(dbarts_sampler*, const double*, int, int);
static int (*p_updatePredictor)(dbarts_sampler*, const double*, const size_t*,
                                size_t, int, int);
static void (*p_setTestPredictors)(dbarts_sampler*, const double*, size_t);
static void (*p_predict)(dbarts_sampler*, const double*, size_t,
                         const double*, double*);
static void (*p_setTreeStorage)(dbarts_sampler*, int, size_t);
static SEXP (*p_getTrees)(dbarts_sampler*, const size_t*, size_t,
                          const size_t*, size_t, const size_t*, size_t, int,
                          const double*);
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
  LOOKUP(void (*)(dbarts_sampler*), p_sampleNodeParametersFromPrior,
         "dbarts_sampler_sampleNodeParametersFromPrior");
  LOOKUP(void (*)(dbarts_sampler*, const double*), p_setResponse,
         "dbarts_sampler_setResponse");
  LOOKUP(void (*)(dbarts_sampler*, const double*, int), p_setOffset,
         "dbarts_sampler_setOffset");
  LOOKUP(void (*)(dbarts_sampler*, const double*), p_setWeights,
         "dbarts_sampler_setWeights");
  LOOKUP(void (*)(dbarts_sampler*, double), p_setSigma,
         "dbarts_sampler_setSigma");
  LOOKUP(void (*)(dbarts_sampler*, dbarts_sampler_callback, void*),
         p_setCallback, "dbarts_sampler_setCallback");
  LOOKUP(void (*)(dbarts_sampler*, const double*), p_setTestOffset,
         "dbarts_sampler_setTestOffset");
  LOOKUP(void (*)(dbarts_sampler*, const size_t*, size_t, const size_t*,
                  size_t, const size_t*, size_t),
         p_printTrees, "dbarts_sampler_printTrees");
  LOOKUP(void (*)(dbarts_sampler*, size_t), p_setNumThreads,
         "dbarts_sampler_setNumThreads");
  LOOKUP(void (*)(dbarts_sampler*, size_t), p_setNumThin,
         "dbarts_sampler_setNumThin");
  LOOKUP(void (*)(dbarts_sampler*, int, uint32_t), p_setVerbose,
         "dbarts_sampler_setVerbose");
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

  dbarts_results results = {0};
  results.structSize = sizeof results;
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
  size_t n = p_numObservations(sampler);
  size_t chains = p_numChains(sampler);

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

  p_run(sampler, numBurnIn, numSamples, &results);

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
    p_setSigma(sampler, state->sigmas[sweepIndex]);
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
  size_t n = p_numObservations(sampler);
  size_t p = p_numPredictors(sampler);
  size_t chains = p_numChains(sampler);

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

  p_setCallback(sampler, sweepCallback, &state);
  p_run(sampler, numBurnIn, numSamples, &results);
  p_setCallback(sampler, NULL, NULL);

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
  size_t chains = p_numChains(sampler);

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

  p_run(sampler, numBurnIn, numSamples, &results);

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
  size_t n = p_numObservations(sampler);
  size_t chains = p_numChains(sampler);

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

  p_run(sampler, numBurnIn, numSamples, &results);

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
  p_sampleNodeParametersFromPrior(samplerFromExpr(ptrExpr));
  return R_NilValue;
}

SEXP capi_set_response(SEXP ptrExpr, SEXP yExpr) {
  p_setResponse(samplerFromExpr(ptrExpr), REAL(yExpr));
  return R_NilValue;
}

SEXP capi_set_weights(SEXP ptrExpr, SEXP weightsExpr) {
  p_setWeights(samplerFromExpr(ptrExpr), REAL(weightsExpr));
  return R_NilValue;
}

SEXP capi_set_test_offset(SEXP ptrExpr, SEXP offsetExpr) {
  p_setTestOffset(samplerFromExpr(ptrExpr),
                  Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr));
  return R_NilValue;
}

/* prints the first tree of the first chain, exercising the entry point;
 * the R side captures the console output */
SEXP capi_print_trees(SEXP ptrExpr) {
  size_t chainIndex = 0, treeIndex = 0;
  p_printTrees(samplerFromExpr(ptrExpr), &chainIndex, 1, NULL, 0,
               &treeIndex, 1);
  return R_NilValue;
}

SEXP capi_set_run_controls(SEXP ptrExpr, SEXP numThreadsExpr,
                           SEXP numThinExpr, SEXP verboseExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  p_setNumThreads(sampler, (size_t) Rf_asInteger(numThreadsExpr));
  p_setNumThin(sampler, (size_t) Rf_asInteger(numThinExpr));
  p_setVerbose(sampler, Rf_asLogical(verboseExpr) == TRUE, 100);
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
                    treeIndices, numTrees, useLiveTrees, NULL);
}

SEXP capi_store_state(SEXP ptrExpr) {
  return p_storeState(samplerFromExpr(ptrExpr));
}

SEXP capi_set_state(SEXP ptrExpr, SEXP stateExpr) {
  p_setState(samplerFromExpr(ptrExpr), stateExpr);
  return R_NilValue;
}
