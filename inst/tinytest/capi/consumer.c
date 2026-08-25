/* A minimal consumer of the flat C API (dbarts/dbarts.h), compiled by
 * test-capi.R with R CMD SHLIB against the installed package headers. It drives
 * the entry points through the header's DBARTS_USE_STUBS stubs - the supported
 * LinkingTo path, where each dbarts_sampler_* call binds to a cached
 * R_GetCCallable pointer generated inside dbarts.h - and exposes .Call wrappers
 * for the R-side assertions. One entry point is still resolved by hand as a
 * deliberate canary (see p_apiHash_raw). */

#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

#include <math.h>   /* fabs */
#include <stdio.h>  /* snprintf */
#include <string.h> /* memcpy, strcmp */

#include <R_ext/Rdynload.h> /* R_GetCCallable, for the raw canary below */

/* Deliberate canary: dbarts_apiHash is ALSO resolved the old way, by hand
 * through R_GetCCallable with a hand-written cast - the un-stubbed per-symbol
 * path a consumer that declines DBARTS_USE_STUBS (or a diagnostic tool) still
 * relies on. Everything else goes through the stubs, so this one raw path
 * guards that plain R_RegisterCCallable registration keeps working on its own.
 * It guards the signature token specifically because that token is the only
 * runtime signal a stale consumer binary trips while the version constants
 * stay put. */
static uint64_t (*p_apiHash_raw)(void);

static void initCanary(void) {
  if (p_apiHash_raw == NULL)
    p_apiHash_raw =
      (uint64_t (*)(void)) R_GetCCallable("dbarts", "dbarts_apiHash");
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

/* the signature token as text (it does not fit an R integer), plus whether the
 * raw canary path and the stubs agree on it, plus whether the installed
 * library's token equals the one this consumer compiled against */
SEXP capi_hash(void) {
  char text[32];
  uint64_t stubbed = dbarts_apiHash();
  initCanary();
  snprintf(text, sizeof(text), "0x%016llx", (unsigned long long) stubbed);

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, Rf_mkString(text));
  SET_VECTOR_ELT(result, 1, Rf_ScalarLogical(p_apiHash_raw() == stubbed));
  SET_VECTOR_ELT(result, 2, Rf_ScalarLogical(stubbed == DBARTS_C_API_HASH));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("text"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("raw.agrees"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("matches.header"));
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(2);
  return result;
}

/* the two version components, through the stubs; the R side checks they agree
 * with the header macros */
SEXP capi_versions(void) {
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 2));
  int* v = INTEGER(result);
  v[0] = dbarts_apiMajorVersion();
  v[1] = dbarts_apiMinorVersion();
  UNPROTECT(1);
  return result;
}

/* An R-built dbarts_predictor_source: the list's elements map onto the
 * struct's members one for one, so the R side can hand the entries a dense
 * block, a CSC triple, a mixed map, or a deliberately malformed argument. The
 * source borrows every array from the list, which stays protected by the
 * caller for the duration of the entry-point call. */
static SEXP getElement(SEXP list, const char* name) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_isNull(names)) return R_NilValue;
  for (R_xlen_t i = 0; i < Rf_xlength(list); ++i)
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
      return VECTOR_ELT(list, i);
  return R_NilValue;
}

static dbarts_predictor_source sourceFromList(SEXP spec) {
  dbarts_predictor_source source = DBARTS_PREDICTOR_SOURCE_INIT;
  SEXP element;

  source.numRows = (size_t) Rf_asInteger(getElement(spec, "numRows"));
  source.numColumns = (size_t) Rf_asInteger(getElement(spec, "numColumns"));

  element = getElement(spec, "dense");
  if (!Rf_isNull(element)) source.denseValues = REAL(element);

  element = getElement(spec, "cscColumnPointers");
  if (!Rf_isNull(element)) {
    source.numCscColumns = (size_t) (Rf_xlength(element) - 1);
    source.cscColumnPointers = INTEGER(element);
    source.cscRowIndices = INTEGER(getElement(spec, "cscRowIndices"));
    source.cscValues = REAL(getElement(spec, "cscValues"));
  }
  /* a caller may declare fewer CSC columns than its pointer array carries,
   * which is how the out-of-range decode is driven */
  element = getElement(spec, "numCscColumns");
  if (!Rf_isNull(element))
    source.numCscColumns = (size_t) Rf_asInteger(element);

  element = getElement(spec, "columnSources");
  if (!Rf_isNull(element)) source.columnSources = (const int32_t*) INTEGER(element);
  element = getElement(spec, "columnTypes");
  if (!Rf_isNull(element)) source.columnTypes = (const int32_t*) INTEGER(element);
  element = getElement(spec, "categoryCounts");
  if (!Rf_isNull(element))
    source.categoryCounts = (const uint32_t*) INTEGER(element);
  element = getElement(spec, "referenceCodes");
  if (!Rf_isNull(element))
    source.referenceCodes = (const int32_t*) INTEGER(element);
  return source;
}

/* the string -> dbarts_family table every wrapper below that used to hand the
 * flat entries a bare family string now goes through, since dbarts_sampler_
 * create/drawLatents/workingResponse take the enum by value */
static int familyFromString(const char* name) {
  if (name[0] == '\0') return DBARTS_FAMILY_AUTO;
  if (strcmp(name, "gaussian") == 0) return DBARTS_FAMILY_GAUSSIAN;
  if (strcmp(name, "probit") == 0) return DBARTS_FAMILY_PROBIT;
  if (strcmp(name, "logistic") == 0) return DBARTS_FAMILY_LOGISTIC;
  if (strcmp(name, "aft") == 0) return DBARTS_FAMILY_AFT;
  if (strcmp(name, "ordinal") == 0) return DBARTS_FAMILY_ORDINAL;
  if (strcmp(name, "nbinom") == 0) return DBARTS_FAMILY_NBINOM;
  if (strcmp(name, "student") == 0) return DBARTS_FAMILY_STUDENT;
  if (strcmp(name, "multinomial") == 0) return DBARTS_FAMILY_MULTINOMIAL;
  Rf_error("familyFromString: unrecognized family \"%s\"", name);
  return DBARTS_FAMILY_AUTO; /* unreached */
}

SEXP capi_create(SEXP control, SEXP model, SEXP data, SEXP family) {
  const char* familyName =
    Rf_isNull(family) ? "" : CHAR(STRING_ELT(family, 0));
  dbarts_sampler* sampler = dbarts_sampler_create(
    control, model, data, familyFromString(familyName));
  SEXP result = PROTECT(R_MakeExternalPtr(sampler, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, samplerFinalizer, FALSE);
  UNPROTECT(1);
  return result;
}

/* the admission probes: an unmapped int driven straight through, exactly as a
 * miscompiled or hand-rolled caller (never going through familyFromString)
 * would send it */
SEXP capi_create_raw_family(SEXP control, SEXP model, SEXP data,
                            SEXP familyInt) {
  dbarts_sampler* sampler = dbarts_sampler_create(
    control, model, data, Rf_asInteger(familyInt));
  SEXP result = PROTECT(R_MakeExternalPtr(sampler, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, samplerFinalizer, FALSE);
  UNPROTECT(1);
  return result;
}

/* every dbarts_family value, in header order, so the R side can check that
 * this consumer's compiled-in numbering agrees with the installed header's */
SEXP capi_family_constants(void) {
  static const char* const names[9] = {
    "auto", "gaussian", "probit", "logistic", "aft",
    "ordinal", "nbinom", "student", "multinomial"
  };
  int values[9] = {
    DBARTS_FAMILY_AUTO, DBARTS_FAMILY_GAUSSIAN, DBARTS_FAMILY_PROBIT,
    DBARTS_FAMILY_LOGISTIC, DBARTS_FAMILY_AFT, DBARTS_FAMILY_ORDINAL,
    DBARTS_FAMILY_NBINOM, DBARTS_FAMILY_STUDENT, DBARTS_FAMILY_MULTINOMIAL
  };
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 9));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 9));
  for (int i = 0; i < 9; ++i) {
    INTEGER(result)[i] = values[i];
    SET_STRING_ELT(namesExpr, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(2);
  return result;
}

SEXP capi_sampler_family(SEXP ptrExpr) {
  return Rf_ScalarInteger(dbarts_sampler_family(samplerFromExpr(ptrExpr)));
}

SEXP capi_dims(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  SEXP result = PROTECT(Rf_allocVector(INTSXP, 8));
  int* dims = INTEGER(result);
  dims[0] = (int) dbarts_sampler_numObservations(sampler);
  dims[1] = (int) dbarts_sampler_numPredictors(sampler);
  dims[2] = (int) dbarts_sampler_numTestObservations(sampler);
  dims[3] = (int) dbarts_sampler_numChains(sampler);
  dims[4] = (int) dbarts_sampler_numTrees(sampler, 0);
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

  dbarts_results results = DBARTS_RESULTS_INIT;
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

/* the zero-structSize guard: a caller that value-inits dbarts_results but
 * forgets to set structSize (rather than using DBARTS_RESULTS_INIT) must be
 * rejected outright - never silently handed an all-skip no-op with an
 * uninitialized buffer. This run must Rf_error, so from R it surfaces as an
 * error; reaching the return would mean the guard failed to fire. */
SEXP capi_run_zero_structsize(SEXP ptrExpr, SEXP numBurnInExpr, SEXP numSamplesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  size_t chains = dbarts_sampler_numChains(sampler);

  double* sigma = (double*) R_alloc(numSamples * chains, sizeof(double));
  double* train = (double*) R_alloc(n * numSamples * chains, sizeof(double));

  dbarts_results results = {0}; /* structSize deliberately left 0 */
  results.sigma = sigma;
  results.train = train;

  dbarts_sampler_run(sampler, numBurnIn, numSamples, &results);
  return Rf_ScalarLogical(1); /* unreachable: the guard must have errored */
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

  dbarts_results results = DBARTS_RESULTS_INIT;
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

  dbarts_results results = DBARTS_RESULTS_INIT;
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

  dbarts_results results = DBARTS_RESULTS_INIT;
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

SEXP capi_set_response(SEXP ptrExpr, SEXP yExpr, SEXP updateScaleExpr) {
  dbarts_sampler_setResponse(samplerFromExpr(ptrExpr), REAL(yExpr),
                             Rf_asLogical(updateScaleExpr) == TRUE);
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

/* prints the first tree of the first chain of the named forest, exercising the
 * entry point; the R side captures the console output. useLiveTrees forwards
 * unchanged: the saved-sample count only matters when it is FALSE. */
SEXP capi_print_trees(SEXP ptrExpr, SEXP useLiveTreesExpr, SEXP forestExpr) {
  size_t chainIndex = 0, treeIndex = 0, sampleIndex = 0;
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  int useLiveTrees = Rf_asLogical(useLiveTreesExpr) == TRUE;
  size_t numSamples =
    !useLiveTrees && dbarts_sampler_numSavedSamples(sampler) > 0 ? 1 : 0;
  dbarts_sampler_printTrees(sampler, &chainIndex, 1, &sampleIndex, numSamples,
                            &treeIndex, 1, useLiveTrees,
                            (size_t) Rf_asInteger(forestExpr));
  return R_NilValue;
}

SEXP capi_num_trees(SEXP ptrExpr, SEXP forestExpr) {
  return Rf_ScalarInteger((int) dbarts_sampler_numTrees(
    samplerFromExpr(ptrExpr), (size_t) Rf_asInteger(forestExpr)));
}

SEXP capi_num_forests(SEXP ptrExpr) {
  return Rf_ScalarInteger(
    (int) dbarts_sampler_numForests(samplerFromExpr(ptrExpr)));
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

/* the dispersion channel a count host reads, both spellings. The recorded slot
 * is NA-poisoned before the run, so a library that never fills it reads back as
 * NA rather than as a plausible number, and the second run pins structSize
 * below the appended field over a poisoned pointer: a size-blind write would
 * dereference it on the one family that HAS a dispersion to write. */
SEXP capi_run_dispersion(SEXP ptrExpr, SEXP numBurnInExpr,
                         SEXP numSamplesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t length = numSamples * dbarts_sampler_numChains(sampler);

  SEXP recorded = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  for (size_t i = 0; i < length; ++i) REAL(recorded)[i] = NA_REAL;
  double* sigma = (double*) R_alloc(length, sizeof(double));

  dbarts_results older = DBARTS_RESULTS_INIT;
  older.structSize = offsetof(dbarts_results, dispersion);
  older.sigma = sigma;
  older.dispersion = (double*) (uintptr_t) 0x1;
  dbarts_sampler_run(sampler, numBurnIn, numSamples, &older);

  /* second, so the state the getter reads afterwards is this run's last draw */
  dbarts_results results = DBARTS_RESULTS_INIT;
  results.sigma = sigma;
  results.dispersion = REAL(recorded);
  dbarts_sampler_run(sampler, 0, numSamples, &results);

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, recorded);
  SET_VECTOR_ELT(
    result, 1, Rf_ScalarLogical(DBARTS_RESULTS_HAS(&results, dispersion)));
  SET_VECTOR_ELT(result, 2, Rf_ScalarLogical(1));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("recorded"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("present"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("guarded"));
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(3);
  return result;
}

/* the Student-t df channel a robust host reads: the results slot appended to
 * dbarts_results after the dispersion one. Same discipline as the dispersion
 * shim above - the recorded slot is NA-poisoned before the run, so an error law
 * that never fills it reads back as NA rather than as a plausible number, and
 * the first run pins structSize below the appended field over a poisoned
 * pointer, which a size-blind write would dereference. */
SEXP capi_run_residual_df(SEXP ptrExpr, SEXP numBurnInExpr,
                          SEXP numSamplesExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numBurnIn = (size_t) Rf_asInteger(numBurnInExpr);
  size_t numSamples = (size_t) Rf_asInteger(numSamplesExpr);
  size_t length = numSamples * dbarts_sampler_numChains(sampler);

  SEXP recorded = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  for (size_t i = 0; i < length; ++i) REAL(recorded)[i] = NA_REAL;
  double* sigma = (double*) R_alloc(length, sizeof(double));

  dbarts_results older = DBARTS_RESULTS_INIT;
  older.structSize = offsetof(dbarts_results, residualDf);
  older.sigma = sigma;
  older.residualDf = (double*) (uintptr_t) 0x1;
  dbarts_sampler_run(sampler, numBurnIn, numSamples, &older);

  dbarts_results results = DBARTS_RESULTS_INIT;
  results.sigma = sigma;
  results.residualDf = REAL(recorded);
  dbarts_sampler_run(sampler, 0, numSamples, &results);

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, recorded);
  SET_VECTOR_ELT(
    result, 1, Rf_ScalarLogical(DBARTS_RESULTS_HAS(&results, residualDf)));
  SET_VECTOR_ELT(result, 2, Rf_ScalarLogical(1));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("recorded"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("present"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("guarded"));
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(3);
  return result;
}

/* the mid-sweep getter, on the capability contract every reader here follows:
 * NULL stands for the 0 return a family carrying no dispersion answers with */
SEXP capi_dispersion(SEXP ptrExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  SEXP result = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) dbarts_sampler_numChains(sampler)));
  int carries = dbarts_sampler_getDispersion(sampler, REAL(result));
  UNPROTECT(1);
  return carries ? result : R_NilValue;
}

/* the wrapped augmentation entries, over caller-supplied arrays: an R NULL is
 * the absent argument, and a scalar's NA is the one no law reads */
static const double* optionalReal(SEXP x) {
  return Rf_isNull(x) ? NULL : REAL(x);
}

SEXP capi_draw_latents(SEXP familyExpr, SEXP fitExpr, SEXP yExpr,
                       SEXP weightsExpr, SEXP offsetExpr, SEXP sigmaExpr,
                       SEXP dispersionExpr, SEXP cutpointsExpr, SEXP dfExpr) {
  size_t n = (size_t) Rf_xlength(fitExpr);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n));
  dbarts_drawLatents(familyFromString(CHAR(STRING_ELT(familyExpr, 0))), n,
                     REAL(fitExpr), REAL(yExpr), optionalReal(weightsExpr),
                     optionalReal(offsetExpr), Rf_asReal(sigmaExpr),
                     Rf_asReal(dispersionExpr), optionalReal(cutpointsExpr),
                     (size_t) Rf_xlength(cutpointsExpr), Rf_asReal(dfExpr),
                     REAL(result));
  UNPROTECT(1);
  return result;
}

SEXP capi_working_response(SEXP familyExpr, SEXP latentExpr, SEXP yExpr,
                           SEXP weightsExpr, SEXP offsetExpr,
                           SEXP dispersionExpr) {
  size_t n = (size_t) Rf_xlength(latentExpr);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n));
  dbarts_workingResponse(familyFromString(CHAR(STRING_ELT(familyExpr, 0))), n,
                         REAL(latentExpr), REAL(yExpr),
                         optionalReal(weightsExpr),
                         optionalReal(offsetExpr), Rf_asReal(dispersionExpr),
                         REAL(result));
  UNPROTECT(1);
  return result;
}

/* the dense spelling every wrapper below hands the entries: the header's own
 * constructor over an R matrix */
static dbarts_predictor_source denseSource(SEXP xExpr) {
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  return dbarts_dense_predictor_source(
    REAL(xExpr), (size_t) INTEGER(dims)[0], (size_t) INTEGER(dims)[1]);
}

SEXP capi_set_predictor(SEXP ptrExpr, SEXP xExpr) {
  dbarts_predictor_source source = denseSource(xExpr);
  return Rf_ScalarLogical(
    dbarts_sampler_setPredictor(samplerFromExpr(ptrExpr), &source, FALSE,
                                TRUE));
}

SEXP capi_update_predictor(SEXP ptrExpr, SEXP xExpr, SEXP columnExpr) {
  size_t column = (size_t) Rf_asInteger(columnExpr); /* already 0-based */
  dbarts_predictor_source source = denseSource(xExpr);
  return Rf_ScalarLogical(dbarts_sampler_updatePredictor(
    samplerFromExpr(ptrExpr), &source, &column, 1, FALSE, TRUE));
}

/* a transactional column update against the EXISTING cut grid: the flavor a
 * replacement that would empty leaves can reach at all, since refreshing the
 * cuts from a collapsed column is refused for its cut count first */
SEXP capi_update_predictor_fixed_cuts(SEXP ptrExpr, SEXP xExpr,
                                      SEXP columnExpr) {
  size_t column = (size_t) Rf_asInteger(columnExpr); /* already 0-based */
  dbarts_predictor_source source = denseSource(xExpr);
  return Rf_ScalarLogical(dbarts_sampler_updatePredictor(
    samplerFromExpr(ptrExpr), &source, &column, 1, FALSE, FALSE));
}

/* the forced flavors of the two above, exercising the transactional guard's
 * accept arm: forceUpdate = TRUE bypasses the empty-leaf veto and always
 * installs */
SEXP capi_set_predictor_forced(SEXP ptrExpr, SEXP xExpr) {
  dbarts_predictor_source source = denseSource(xExpr);
  return Rf_ScalarLogical(
    dbarts_sampler_setPredictor(samplerFromExpr(ptrExpr), &source, TRUE,
                                TRUE));
}

SEXP capi_update_predictor_forced(SEXP ptrExpr, SEXP xExpr, SEXP columnExpr) {
  size_t column = (size_t) Rf_asInteger(columnExpr); /* already 0-based */
  dbarts_predictor_source source = denseSource(xExpr);
  return Rf_ScalarLogical(dbarts_sampler_updatePredictor(
    samplerFromExpr(ptrExpr), &source, &column, 1, TRUE, TRUE));
}

SEXP capi_set_test_predictors(SEXP ptrExpr, SEXP xTestExpr) {
  if (Rf_isNull(xTestExpr)) {
    dbarts_sampler_setTestPredictors(samplerFromExpr(ptrExpr), NULL);
  } else {
    dbarts_predictor_source source = denseSource(xTestExpr);
    dbarts_sampler_setTestPredictors(samplerFromExpr(ptrExpr), &source);
  }
  return R_NilValue;
}

/* the source-shaped flavors, over an R-built spec: dense, CSC, mixed, or
 * malformed. Every entry that takes a source has one, since the refusals are
 * stated once and must fire at every funnel. */
SEXP capi_set_predictor_source(SEXP ptrExpr, SEXP specExpr) {
  dbarts_predictor_source source = sourceFromList(specExpr);
  return Rf_ScalarLogical(dbarts_sampler_setPredictor(
    samplerFromExpr(ptrExpr), &source, FALSE, TRUE));
}

SEXP capi_update_predictor_source(SEXP ptrExpr, SEXP specExpr,
                                  SEXP columnsExpr) {
  dbarts_predictor_source source = sourceFromList(specExpr);
  size_t numColumns = (size_t) Rf_xlength(columnsExpr);
  size_t* columns = (size_t*) R_alloc(numColumns, sizeof(size_t));
  for (size_t k = 0; k < numColumns; ++k)
    columns[k] = (size_t) INTEGER(columnsExpr)[k]; /* already 0-based */
  return Rf_ScalarLogical(dbarts_sampler_updatePredictor(
    samplerFromExpr(ptrExpr), &source, columns, numColumns, FALSE, TRUE));
}

SEXP capi_set_test_predictors_source(SEXP ptrExpr, SEXP specExpr) {
  dbarts_predictor_source source = sourceFromList(specExpr);
  dbarts_sampler_setTestPredictors(samplerFromExpr(ptrExpr), &source);
  return R_NilValue;
}

SEXP capi_predict_source(SEXP ptrExpr, SEXP specExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  dbarts_predictor_source source = sourceFromList(specExpr);
  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length =
    source.numRows * numSamples * dbarts_sampler_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  dbarts_sampler_predict(sampler, &source, NULL, 0, REAL(result));
  UNPROTECT(1);
  return result;
}

/* the input-side write guard, inverted for a READ: an old, smaller caller
 * pins structSize below columnTypes and poisons every field past that
 * boundary, so an entry that read a member it was not handed would fault on
 * the unmapped page. Returns the fits, which must equal the dense answer. */
SEXP capi_predict_truncated(SEXP ptrExpr, SEXP xTestExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  SEXP dims = Rf_getAttrib(xTestExpr, R_DimSymbol);
  void* poison = (void*) (uintptr_t) 0x1;

  dbarts_predictor_source source;
  memset(&source, 0, sizeof(source));
  source.structSize = offsetof(dbarts_predictor_source, columnTypes);
  source.numRows = (size_t) INTEGER(dims)[0];
  source.numColumns = (size_t) INTEGER(dims)[1];
  source.denseValues = REAL(xTestExpr);
  source.columnTypes = (const int32_t*) poison;
  source.categoryCounts = (const uint32_t*) poison;
  source.referenceCodes = (const int32_t*) poison;

  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length =
    source.numRows * numSamples * dbarts_sampler_numChains(sampler);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  dbarts_sampler_predict(sampler, &source, NULL, 0, REAL(result));
  UNPROTECT(1);
  return result;
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
  dbarts_predictor_source source = denseSource(xTestExpr);
  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length =
    source.numRows * numSamples * dbarts_sampler_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  /* 0 is the header's "the sampler's own count", so every assertion that
   * runs through here covers that resolution as well */
  dbarts_sampler_predict(sampler, &source,
                         Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr), 0,
                         REAL(result));
  UNPROTECT(1);
  return result;
}

/* the same replay at an explicit per-call count. The count does not persist
 * and cannot move a value, so the answer must equal capi_predict's at every
 * one of them. */
SEXP capi_predict_threads(SEXP ptrExpr, SEXP xTestExpr, SEXP nThreadsExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  dbarts_predictor_source source = denseSource(xTestExpr);
  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length =
    source.numRows * numSamples * dbarts_sampler_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  dbarts_sampler_predict(sampler, &source, NULL,
                         (size_t) Rf_asInteger(nThreadsExpr), REAL(result));
  UNPROTECT(1);
  return result;
}

/* sampleNumsExpr null reads every saved sample; otherwise its 1-based indices
 * select the saved samples, so the caller can compare the all-samples table
 * against a per-sample gather (the consistency stan4bart's extract relies on). */
SEXP capi_get_trees(SEXP ptrExpr, SEXP useLiveTreesExpr, SEXP sampleNumsExpr,
                    SEXP forestExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  int useLiveTrees = Rf_asLogical(useLiveTreesExpr) == TRUE;
  size_t forest = (size_t) Rf_asInteger(forestExpr);

  size_t numChains = dbarts_sampler_numChains(sampler);
  size_t numSaved = useLiveTrees ? 0 : dbarts_sampler_numSavedSamples(sampler);
  size_t numTrees = dbarts_sampler_numTrees(sampler, forest);

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
                                 numTrees, useLiveTrees, forest);
}

SEXP capi_store_state(SEXP ptrExpr) {
  return dbarts_sampler_storeState(samplerFromExpr(ptrExpr));
}

SEXP capi_set_state(SEXP ptrExpr, SEXP stateExpr) {
  dbarts_sampler_setState(samplerFromExpr(ptrExpr), stateExpr);
  return R_NilValue;
}

/* the per-observation 0/1 active-row mask; a NULL clears it. The values are
 * consumed during the call, so nothing here has to outlive it. */
SEXP capi_set_active_rows(SEXP ptrExpr, SEXP activeExpr) {
  return Rf_ScalarInteger(dbarts_sampler_setActiveRows(
    samplerFromExpr(ptrExpr),
    Rf_isNull(activeExpr) ? NULL : REAL(activeExpr)));
}

/* the per-forest precision weight. BORROWED until replaced, unlike the mask
 * above, so the R side keeps its vector alive for the sampler's life. */
SEXP capi_set_forest_weights(SEXP ptrExpr, SEXP forestExpr,
                             SEXP weightsExpr) {
  return Rf_ScalarInteger(dbarts_sampler_setForestWeights(
    samplerFromExpr(ptrExpr), (size_t) Rf_asInteger(forestExpr),
    Rf_isNull(weightsExpr) ? NULL : REAL(weightsExpr)));
}

/* the mean channel: the basis a forest's amplitudes multiply, column-major
 * numObservations x numColumns and copied by the entry */
SEXP capi_set_forest_basis(SEXP ptrExpr, SEXP forestExpr, SEXP basisExpr) {
  SEXP dims = Rf_getAttrib(basisExpr, R_DimSymbol);
  size_t numColumns =
    Rf_isNull(dims) ? 1 : (size_t) INTEGER(dims)[1];
  return Rf_ScalarInteger(dbarts_sampler_setForestBasis(
    samplerFromExpr(ptrExpr), (size_t) Rf_asInteger(forestExpr),
    REAL(basisExpr), numColumns));
}

/* the ragged amplitude read: the count first, so the caller sizes its own
 * buffer, then the values as numForestAmplitudes x numChains */
SEXP capi_forest_amplitudes(SEXP ptrExpr, SEXP forestExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t forest = (size_t) Rf_asInteger(forestExpr);
  size_t numAmplitudes = dbarts_sampler_numForestAmplitudes(sampler, forest);
  size_t numChains = dbarts_sampler_numChains(sampler);

  SEXP valuesExpr = PROTECT(
    Rf_allocVector(REALSXP, (R_xlen_t) (numAmplitudes * numChains)));
  int accepted =
    dbarts_sampler_getForestAmplitudes(sampler, forest, REAL(valuesExpr));

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, Rf_ScalarInteger((int) numAmplitudes));
  SET_VECTOR_ELT(result, 1, valuesExpr);
  SET_VECTOR_ELT(result, 2, Rf_ScalarInteger(accepted));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("count"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("values"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("accepted"));
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(3);
  return result;
}

/* one forest's calibration through the size-first output struct. mode is the
 * caller's structSize: 0 carries the whole struct, 1 stops below leafModel
 * (the omitting caller the presence test exists for) and 2 stops at the
 * calibration-map append boundary - a PRE-APPEND caller, compiled against the
 * struct as it shipped before the five map fields, which must still read the
 * eight original members and leave the five it does not carry alone. Every
 * member past a mode's boundary is poisoned, so a fill that ignored structSize
 * would write through 0x1 and crash rather than quietly pass. skipK leaves the
 * k pointer null, the null-member-skips half of the same contract. */
SEXP capi_forest_calibration(SEXP ptrExpr, SEXP forestExpr, SEXP modeExpr,
                             SEXP skipKExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t numChains = dbarts_sampler_numChains(sampler);
  int mode = Rf_asInteger(modeExpr);
  int skipK = Rf_asLogical(skipKExpr) == TRUE;

  SEXP priorScaleExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP priorSdExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP priorMeanExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP kExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP scaleExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP shiftExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP hyperExpr = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) numChains));
  SEXP leafExpr = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) numChains));
  SEXP aVarExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP aScaleExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP factorExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP divisorExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  SEXP rowNormExpr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) numChains));
  for (size_t c = 0; c < numChains; ++c) {
    REAL(kExpr)[c] = -1.0;
    INTEGER(leafExpr)[c] = -1;
    REAL(aVarExpr)[c] = -1.0;
    REAL(aScaleExpr)[c] = -1.0;
    REAL(factorExpr)[c] = -1.0;
    REAL(divisorExpr)[c] = -1.0;
    REAL(rowNormExpr)[c] = -1.0;
  }

  double* poison = (double*) (uintptr_t) 0x1; /* never read, never written */
  dbarts_forest_calibration calibration = DBARTS_FOREST_CALIBRATION_INIT;
  calibration.priorScale = REAL(priorScaleExpr);
  calibration.priorSd = REAL(priorSdExpr);
  calibration.priorMean = REAL(priorMeanExpr);
  calibration.k = skipK ? NULL : REAL(kExpr);
  calibration.responseScale = REAL(scaleExpr);
  calibration.responseShift = REAL(shiftExpr);
  calibration.kHasHyperprior = INTEGER(hyperExpr);
  calibration.leafModel = mode == 1 ? (int32_t*) (uintptr_t) 0x1
                                    : INTEGER(leafExpr);
  if (mode == 0) {
    calibration.amplitudePriorVariance = REAL(aVarExpr);
    calibration.amplitudePriorScale = REAL(aScaleExpr);
    calibration.nodeScaleFactor = REAL(factorExpr);
    calibration.nodeScaleDivisor = REAL(divisorExpr);
    calibration.basisRowNorm = REAL(rowNormExpr);
  } else {
    calibration.structSize =
      mode == 1
        ? offsetof(dbarts_forest_calibration, leafModel)
        : offsetof(dbarts_forest_calibration, amplitudePriorVariance);
    calibration.amplitudePriorVariance = poison;
    calibration.amplitudePriorScale = poison;
    calibration.nodeScaleFactor = poison;
    calibration.nodeScaleDivisor = poison;
    calibration.basisRowNorm = poison;
  }

  int accepted =
    dbarts_sampler_getForestCalibration(sampler, (size_t) Rf_asInteger(forestExpr),
                                        &calibration);

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 14));
  SET_VECTOR_ELT(result, 0, priorScaleExpr);
  SET_VECTOR_ELT(result, 1, priorSdExpr);
  SET_VECTOR_ELT(result, 2, priorMeanExpr);
  SET_VECTOR_ELT(result, 3, kExpr);
  SET_VECTOR_ELT(result, 4, scaleExpr);
  SET_VECTOR_ELT(result, 5, shiftExpr);
  SET_VECTOR_ELT(result, 6, hyperExpr);
  SET_VECTOR_ELT(result, 7, leafExpr);
  SET_VECTOR_ELT(result, 8, aVarExpr);
  SET_VECTOR_ELT(result, 9, aScaleExpr);
  SET_VECTOR_ELT(result, 10, factorExpr);
  SET_VECTOR_ELT(result, 11, divisorExpr);
  SET_VECTOR_ELT(result, 12, rowNormExpr);
  SET_VECTOR_ELT(result, 13, Rf_ScalarInteger(accepted));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, 14));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("prior.scale"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("prior.sd"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("prior.mean"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("k"));
  SET_STRING_ELT(namesExpr, 4, Rf_mkChar("response.scale"));
  SET_STRING_ELT(namesExpr, 5, Rf_mkChar("response.shift"));
  SET_STRING_ELT(namesExpr, 6, Rf_mkChar("k.has.hyperprior"));
  SET_STRING_ELT(namesExpr, 7, Rf_mkChar("leaf.model"));
  SET_STRING_ELT(namesExpr, 8, Rf_mkChar("amplitude.prior.variance"));
  SET_STRING_ELT(namesExpr, 9, Rf_mkChar("amplitude.prior.scale"));
  SET_STRING_ELT(namesExpr, 10, Rf_mkChar("node.scale.factor"));
  SET_STRING_ELT(namesExpr, 11, Rf_mkChar("node.scale.divisor"));
  SET_STRING_ELT(namesExpr, 12, Rf_mkChar("basis.row.norm"));
  SET_STRING_ELT(namesExpr, 13, Rf_mkChar("accepted"));
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(15);
  return result;
}

/* the zero-structSize guard on the calibration buffers, the read-side twin of
 * capi_run_zero_structsize: this call must error rather than fill nothing */
SEXP capi_forest_calibration_zero_structsize(SEXP ptrExpr) {
  dbarts_forest_calibration calibration;
  memset(&calibration, 0, sizeof(calibration)); /* structSize left 0 */
  return Rf_ScalarInteger(dbarts_sampler_getForestCalibration(
    samplerFromExpr(ptrExpr), 0, &calibration));
}

SEXP capi_set_forest_prior_scale(SEXP ptrExpr, SEXP forestExpr,
                                 SEXP priorScaleExpr) {
  return Rf_ScalarInteger(dbarts_sampler_setForestPriorScale(
    samplerFromExpr(ptrExpr), (size_t) Rf_asInteger(forestExpr),
    Rf_asReal(priorScaleExpr)));
}

/* The two-forest (BCF) surface. Every verdict below is
 * reached HERE rather than in R: what is under test is that the flat API and
 * the R bridge apply one rule, so each acceptance, each refusal, and the reason
 * a refusal names are checked in the consumer, and the R side only reads the
 * per-leg results off the returned vector. A refusal arrives as an R error,
 * which longjmps out of the entry point, so every leg runs under
 * R_tryCatchError. */
enum {
  LEG_NUM_FORESTS,
  LEG_RESPONSE_PINNED,
  LEG_RESPONSE_RESCALED,
  LEG_OFFSET_PINNED,
  LEG_OFFSET_RESCALED,
  LEG_WEIGHTS,
  LEG_TEST_OFFSET,
  LEG_TEST_PREDICTORS,
  LEG_PREDICT,
  LEG_BASIS,
  LEG_BASIS_FOREST_0,
  LEG_BASIS_WIDTH,
  LEG_BASIS_CONTINUOUS,
  LEG_BASIS_RANGE,
  LEG_FOREST_FITS,
  LEG_AMPLITUDES,
  LEG_FOREST_WEIGHTS,
  LEG_FOREST_WEIGHTS_RANGE,
  LEG_CALIBRATION,
  LEG_COUNT
};

static const char* const legNames[LEG_COUNT] = {
  "numForests",
  "response.pinned",
  "response.rescaled",
  "offset.pinned",
  "offset.rescaled",
  "weights",
  "testOffset",
  "setTestPredictors",
  "predict",
  "setForestBasis",
  "setForestBasis.forest0",
  "setForestBasis.width",
  "setForestBasis.continuous",
  "setForestBasis.range",
  "forestFits",
  "forestAmplitudes",
  "setForestWeights",
  "setForestWeights.range",
  "forestCalibration"
};

/* The refusal each leg must draw, or NULL where it must be accepted.
 *
 * The response, offset and weight refusals come from the guard the R bridge
 * shares with this surface, which has a second branch naming a coupling whose
 * response is its own count matrix. BCF is not that coupling - it opts into the
 * response conduit - so these legs are also the pin that the branch stays
 * conditioned on the capability rather than on the forest count: were it to
 * fire here, the three accepting legs below would refuse. The branch's own
 * message is unreachable from this surface, which has no multinomial creation
 * entry (dbarts_sampler_create builds single-forest and BCF samplers only);
 * inst/tinytest/test-multinomial-counts-mutation.R pins the text. */
static const char* const legRefusals[LEG_COUNT] = {
  NULL,
  NULL,
  "a response swap only with updateScale = FALSE",
  NULL,
  "an offset swap only with updateScale = FALSE",
  NULL,
  "multi-forest sampler fixes its data at creation",
  "forest amplitudes have no off-sample basis",
  "forest amplitudes have no off-sample basis",
  NULL,
  NULL,
  "a basis needs at least one column",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

typedef struct {
  dbarts_sampler* sampler;
  const double* y;
  const double* offset;
  const double* weights;
  const double* z;
  const double* xTest;
  size_t numTestObservations;
  double* out;
  int leg;
  int accepted;
  int errored;
  char message[256];
} bcfLegs;

/* Reinstalls the constructed (1 - z, z) pair on forest 1, so a leg that moved
 * the basis leaves the layout the legs after it read. Reinstalling is the same
 * one operation an install is - there is no second, restoring route - and it is
 * the bitwise identity on the amplitudes while the width does not move. */
static void restoreIndicatorBasis(bcfLegs* legs, size_t n) {
  double* basis = (double*) R_alloc(2 * n, sizeof(double));
  for (size_t i = 0; i < n; ++i) {
    basis[2 * i] = 1.0 - legs->z[i];
    basis[2 * i + 1] = legs->z[i];
  }
  if (dbarts_sampler_setForestBasis(legs->sampler, 1, basis, 2) != 1)
    legs->accepted = 0;
}

static SEXP bcfLegBody(void* data) {
  bcfLegs* legs = (bcfLegs*) data;
  size_t n = dbarts_sampler_numObservations(legs->sampler);
  size_t chains = dbarts_sampler_numChains(legs->sampler);
  switch (legs->leg) {
    case LEG_NUM_FORESTS:
      legs->accepted = dbarts_sampler_numForests(legs->sampler) == 2;
      break;
    case LEG_RESPONSE_PINNED:
      dbarts_sampler_setResponse(legs->sampler, legs->y, 0);
      break;
    case LEG_RESPONSE_RESCALED:
      dbarts_sampler_setResponse(legs->sampler, legs->y, 1);
      break;
    case LEG_OFFSET_PINNED:
      dbarts_sampler_setOffset(legs->sampler, legs->offset, 0);
      break;
    case LEG_OFFSET_RESCALED:
      dbarts_sampler_setOffset(legs->sampler, legs->offset, 1);
      break;
    case LEG_WEIGHTS:
      dbarts_sampler_setWeights(legs->sampler, legs->weights);
      break;
    case LEG_TEST_OFFSET:
      dbarts_sampler_setTestOffset(legs->sampler, legs->offset);
      break;
    case LEG_TEST_PREDICTORS: {
      dbarts_predictor_source source = dbarts_dense_predictor_source(
        legs->xTest, legs->numTestObservations,
        dbarts_sampler_numPredictors(legs->sampler));
      dbarts_sampler_setTestPredictors(legs->sampler, &source);
      break;
    }
    case LEG_PREDICT: {
      dbarts_predictor_source source = dbarts_dense_predictor_source(
        legs->xTest, legs->numTestObservations,
        dbarts_sampler_numPredictors(legs->sampler));
      dbarts_sampler_predict(legs->sampler, &source, NULL, 0, legs->out);
      break;
    }
    case LEG_BASIS: {
      /* the two-column complementary indicator (1 - z, z) the amplitudes
       * contrast on, laid ROW-major (row i at basis + i * numColumns) as the
       * header states and the engine's own contraction reads it */
      double* basis = (double*) R_alloc(2 * n, sizeof(double));
      for (size_t i = 0; i < n; ++i) {
        basis[2 * i] = 1.0 - legs->z[i];
        basis[2 * i + 1] = legs->z[i];
      }
      legs->accepted = dbarts_sampler_setForestBasis(legs->sampler, 1, basis, 2);
      break;
    }
    case LEG_BASIS_FOREST_0: {
      /* forest 0 takes a basis like any other forest - an ACCEPTANCE, where
       * this leg once pinned a capability answer - while an index past the
       * last forest stays a capability answer rather than a raise */
      double* basis = (double*) R_alloc(2 * n, sizeof(double));
      for (size_t i = 0; i < n; ++i) {
        basis[2 * i] = 1.0 - legs->z[i];
        basis[2 * i + 1] = legs->z[i];
      }
      double* ones = (double*) R_alloc(n, sizeof(double));
      for (size_t i = 0; i < n; ++i) ones[i] = 1.0;
      legs->accepted =
        dbarts_sampler_setForestBasis(legs->sampler, 0, basis, 2) == 1 &&
        dbarts_sampler_numForestAmplitudes(legs->sampler, 0) == 2 &&
        dbarts_sampler_setForestBasis(legs->sampler, 2, basis, 2) == 0 &&
        /* narrowing back is the same one operation, so the legs that follow
         * see the constructed layout again */
        dbarts_sampler_setForestBasis(legs->sampler, 0, ones, 1) == 1 &&
        dbarts_sampler_numForestAmplitudes(legs->sampler, 0) == 1;
      break;
    }
    case LEG_BASIS_WIDTH:
      /* a zero-width basis is malformed and RAISES */
      dbarts_sampler_setForestBasis(legs->sampler, 1, legs->z, 0);
      break;
    case LEG_BASIS_CONTINUOUS: {
      /* two columns, not complementary 0/1: LEGAL now, where this leg once
       * pinned a raise. The forest keeps its two amplitudes and moves onto
       * the general conditional, which no return value reports */
      double* basis = (double*) R_alloc(2 * n, sizeof(double));
      for (size_t i = 0; i < n; ++i) {
        basis[2 * i] = 0.25;
        basis[2 * i + 1] = 0.75;
      }
      legs->accepted =
        dbarts_sampler_setForestBasis(legs->sampler, 1, basis, 2) == 1 &&
        dbarts_sampler_numForestAmplitudes(legs->sampler, 1) == 2;
      restoreIndicatorBasis(legs, n);
      break;
    }
    case LEG_BASIS_RANGE: {
      /* a THREE-column basis widens the block, and the ragged read follows */
      double* basis = (double*) R_alloc(3 * n, sizeof(double));
      for (size_t i = 0; i < n; ++i) {
        basis[3 * i] = 1.0;
        basis[3 * i + 1] = legs->z[i];
        basis[3 * i + 2] = 1.0 - legs->z[i];
      }
      legs->accepted =
        dbarts_sampler_setForestBasis(legs->sampler, 1, basis, 3) == 1 &&
        dbarts_sampler_numForestAmplitudes(legs->sampler, 1) == 3 &&
        dbarts_sampler_getForestAmplitudes(legs->sampler, 1, legs->out) == 1;
      for (size_t i = 0; i < 3 * chains; ++i)
        if (!R_FINITE(legs->out[i])) legs->accepted = 0;
      restoreIndicatorBasis(legs, n);
      break;
    }
    case LEG_FOREST_FITS:
      /* both forests read, an index past the last one refuses, and every
       * value written is finite */
      legs->accepted =
        dbarts_sampler_getForestFits(legs->sampler, 0, legs->out) &&
        dbarts_sampler_getForestFits(legs->sampler, 1, legs->out + n * chains) &&
        !dbarts_sampler_getForestFits(legs->sampler, 2, legs->out);
      for (size_t i = 0; i < 2 * n * chains; ++i)
        if (!R_FINITE(legs->out[i])) legs->accepted = 0;
      break;
    case LEG_AMPLITUDES:
      /* the ragged pair: one amplitude for the intercept forest, two for the
       * indicator basis, every value finite, and an index past the last
       * forest refuses the read while the count entry raises */
      legs->accepted =
        dbarts_sampler_numForestAmplitudes(legs->sampler, 0) == 1 &&
        dbarts_sampler_numForestAmplitudes(legs->sampler, 1) == 2 &&
        dbarts_sampler_getForestAmplitudes(legs->sampler, 0, legs->out) &&
        dbarts_sampler_getForestAmplitudes(legs->sampler, 1,
                                        legs->out + chains) &&
        !dbarts_sampler_getForestAmplitudes(legs->sampler, 2, legs->out);
      for (size_t i = 0; i < 3 * chains; ++i)
        if (!R_FINITE(legs->out[i])) legs->accepted = 0;
      break;
    case LEG_FOREST_WEIGHTS:
      /* 1 = accepted, and a null clears; the weights are BORROWED, so the
       * caller's vector outlives the call (the R side holds it) */
      legs->accepted =
        dbarts_sampler_setForestWeights(legs->sampler, 1, legs->weights) == 1 &&
        dbarts_sampler_setForestWeights(legs->sampler, 1, NULL) == 1;
      break;
    case LEG_FOREST_WEIGHTS_RANGE:
      /* 0 = refused, without touching the sampler: a forest past the last one
       * is a capability answer, not a raise */
      legs->accepted =
        dbarts_sampler_setForestWeights(legs->sampler, 2, legs->weights) == 0;
      break;
    case LEG_CALIBRATION: {
      /* both forests read, an index past the last one refuses, and the write
       * is refused on every forest because the two-forest map owns it */
      dbarts_forest_calibration calibration = DBARTS_FOREST_CALIBRATION_INIT;
      calibration.priorScale = legs->out;
      calibration.k = legs->out + chains;
      /* the calibration map's own three, which only a mapped sampler reports:
       * the prognostic forest leaves all three at 1 and the treatment forest
       * carries the half-normal median as its divisor over the synthesized
       * (1 - z, z) pair, whose rows are unit norm */
      calibration.nodeScaleFactor = legs->out + 2 * chains;
      calibration.nodeScaleDivisor = legs->out + 3 * chains;
      calibration.basisRowNorm = legs->out + 4 * chains;
      legs->accepted =
        dbarts_sampler_getForestCalibration(legs->sampler, 0, &calibration);
      double muAnchor = legs->out[0] * legs->out[3 * chains] *
                        legs->out[4 * chains] / legs->out[2 * chains];
      for (size_t c = 0; c < chains; ++c)
        if (legs->out[2 * chains + c] != 1.0 ||
            legs->out[3 * chains + c] != 1.0 ||
            legs->out[4 * chains + c] != 1.0)
          legs->accepted = 0;
      if (!dbarts_sampler_getForestCalibration(legs->sampler, 1, &calibration))
        legs->accepted = 0;
      for (size_t c = 0; c < chains; ++c)
        if (legs->out[3 * chains + c] != 0.674 ||
            legs->out[4 * chains + c] != 1.0)
          legs->accepted = 0;
      /* and the anchor the map states both node scales against is the same one
       * either forest's decomposition recovers */
      double tauAnchor = legs->out[0] * legs->out[3 * chains] *
                         legs->out[4 * chains] / legs->out[2 * chains];
      if (!(fabs(tauAnchor - muAnchor) <= 1.0e-12 * fabs(muAnchor)))
        legs->accepted = 0;
      if (dbarts_sampler_getForestCalibration(legs->sampler, 2, &calibration) ||
          dbarts_sampler_setForestPriorScale(legs->sampler, 0, 2.5) ||
          dbarts_sampler_setForestPriorScale(legs->sampler, 1, 2.5))
        legs->accepted = 0;
      for (size_t i = 0; i < 2 * chains; ++i)
        if (!R_FINITE(legs->out[i]) || legs->out[i] <= 0.0)
          legs->accepted = 0;
      break;
    }
  }
  return R_NilValue;
}

static SEXP bcfLegHandler(SEXP condExpr, void* data) {
  bcfLegs* legs = (bcfLegs*) data;
  legs->errored = 1;
  if (Rf_isVectorList(condExpr) && Rf_xlength(condExpr) > 0) {
    SEXP messageExpr = VECTOR_ELT(condExpr, 0); /* conditionMessage */
    if (Rf_isString(messageExpr) && Rf_xlength(messageExpr) > 0) {
      strncpy(legs->message, CHAR(STRING_ELT(messageExpr, 0)),
              sizeof(legs->message) - 1);
      legs->message[sizeof(legs->message) - 1] = '\0';
    }
  }
  return R_NilValue;
}

/* an accepting leg must return 1 without erroring; a refusing one must error
 * with a message naming its reason */
static int runBCFLeg(bcfLegs* legs, int leg) {
  legs->leg = leg;
  legs->accepted = 1;
  legs->errored = 0;
  legs->message[0] = '\0';
  R_tryCatchError(bcfLegBody, legs, bcfLegHandler, legs);
  if (legRefusals[leg] == NULL) return !legs->errored && legs->accepted;
  return legs->errored && strstr(legs->message, legRefusals[leg]) != NULL;
}

SEXP capi_bcf_surface(SEXP ptrExpr, SEXP yExpr, SEXP offsetExpr,
                      SEXP weightsExpr, SEXP zExpr, SEXP xTestExpr) {
  bcfLegs legs;
  legs.sampler = samplerFromExpr(ptrExpr);
  legs.y = REAL(yExpr);
  legs.offset = REAL(offsetExpr);
  legs.weights = REAL(weightsExpr);
  legs.z = REAL(zExpr);
  legs.xTest = REAL(xTestExpr);
  legs.numTestObservations =
    (size_t) INTEGER(Rf_getAttrib(xTestExpr, R_DimSymbol))[0];
  /* the widest leg is the two forests' fits; the glue is 3 x chains and
   * predict, refused before it reads anything, would want fewer rows still */
  legs.out = (double*) R_alloc(
    2 * dbarts_sampler_numObservations(legs.sampler) *
      dbarts_sampler_numChains(legs.sampler),
    sizeof(double));

  SEXP result = PROTECT(Rf_allocVector(LGLSXP, LEG_COUNT));
  SEXP namesExpr = PROTECT(Rf_allocVector(STRSXP, LEG_COUNT));
  for (int leg = 0; leg < LEG_COUNT; ++leg) {
    LOGICAL(result)[leg] = runBCFLeg(&legs, leg);
    SET_STRING_ELT(namesExpr, leg, Rf_mkChar(legNames[leg]));
  }
  Rf_setAttrib(result, R_NamesSymbol, namesExpr);
  UNPROTECT(2);
  return result;
}
