/* A minimal consumer of the flat C API (dbarts/dbarts.h), compiled by
 * test-capi.R with R CMD SHLIB against the installed package headers. It drives
 * the entry points through the header's DBARTS_USE_STUBS stubs - the supported
 * LinkingTo path, where each dbarts_sampler_* call binds to a cached
 * R_GetCCallable pointer generated inside dbarts.h - and exposes .Call wrappers
 * for the R-side assertions. One entry point is still resolved by hand as a
 * deliberate canary (see p_apiHash_raw). */

/* dbarts.h's prototype view is plain C and brings no R header with it, so a
 * consumer that speaks SEXP - as every .Call wrapper below does - includes R's
 * own headers itself, in whatever order it likes. R_NO_REMAP is this file's
 * choice, not the header's. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

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

/* THE HANDLE: ptrExpr is the external pointer an R dbartsSampler object hands
 * out (its getPointer method), and the address inside it IS the
 * dbarts_sampler* the flat entries take - the whole creation route this
 * consumer has. Nothing here owns it: the R object does, so this file
 * registers no finalizer and frees nothing. */
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
  /* the code channel and its declared width, which a caller may understate to
   * drive the out-of-range refusal exactly as numCscColumns does */
  element = getElement(spec, "denseCodes");
  if (!Rf_isNull(element)) {
    source.denseCodes = (const int32_t*) INTEGER(element);
    source.numDenseCodeColumns =
      (size_t) (Rf_xlength(element) / (R_xlen_t) source.numRows);
  }
  element = getElement(spec, "numDenseCodeColumns");
  if (!Rf_isNull(element))
    source.numDenseCodeColumns = (size_t) Rf_asInteger(element);
  return source;
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

/* the conditioning setters answer a capability status, which the R side reads
 * as an integer: 1 on a mutation, 0 where the sampler carries no such channel
 * at all and nothing was touched */
SEXP capi_set_response(SEXP ptrExpr, SEXP yExpr, SEXP updateScaleExpr) {
  return Rf_ScalarInteger(dbarts_sampler_setResponse(
    samplerFromExpr(ptrExpr), REAL(yExpr),
    Rf_asLogical(updateScaleExpr) == TRUE));
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
  dbarts_sampler_printTrees(sampler, (size_t) Rf_asInteger(forestExpr),
                            &chainIndex, 1, &sampleIndex, numSamples,
                            &treeIndex, 1, useLiveTrees);
  return R_NilValue;
}

SEXP capi_num_trees(SEXP ptrExpr, SEXP forestExpr) {
  return Rf_ScalarInteger((int) dbarts_sampler_numTrees(
    samplerFromExpr(ptrExpr), (size_t) Rf_asInteger(forestExpr)));
}

/* printEvery gets its own entrance: the run-control setter below always hands
 * over a legal one, so the entry point's refusal at 0 is otherwise unreachable
 * from here */
SEXP capi_set_verbose(SEXP ptrExpr, SEXP verboseExpr, SEXP printEveryExpr) {
  dbarts_sampler_setVerbose(samplerFromExpr(ptrExpr),
                            Rf_asLogical(verboseExpr) == TRUE,
                            (size_t) Rf_asInteger(printEveryExpr));
  return R_NilValue;
}

SEXP capi_set_run_controls(SEXP ptrExpr, SEXP numThreadsExpr,
                           SEXP verboseExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  dbarts_sampler_setNumThreads(sampler, (size_t) Rf_asInteger(numThreadsExpr));
  dbarts_sampler_setVerbose(sampler, Rf_asLogical(verboseExpr) == TRUE, 100);
  return R_NilValue;
}

SEXP capi_set_offset(SEXP ptrExpr, SEXP offsetExpr, SEXP updateScaleExpr) {
  return Rf_ScalarInteger(dbarts_sampler_setOffset(
    samplerFromExpr(ptrExpr),
    Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr),
    Rf_asLogical(updateScaleExpr) == TRUE));
}

SEXP capi_set_sigma(SEXP ptrExpr, SEXP sigmaExpr) {
  return Rf_ScalarInteger(
    dbarts_sampler_setSigma(samplerFromExpr(ptrExpr), Rf_asReal(sigmaExpr)));
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

/* the dense spelling every wrapper below hands the entries: the header's own
 * constructor over an R matrix */
static dbarts_predictor_source denseSource(SEXP xExpr) {
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  return dbarts_dense_predictor_source(
    REAL(xExpr), (size_t) INTEGER(dims)[0], (size_t) INTEGER(dims)[1]);
}

/* the predict wrappers hand back NULL on a capability 0, the shape
 * capi_get_latents already has: the buffer is untouched there, so returning it
 * would report uninitialized memory as fits */
SEXP capi_predict_source(SEXP ptrExpr, SEXP specExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  dbarts_predictor_source source = sourceFromList(specExpr);
  size_t saved = dbarts_sampler_numSavedSamples(sampler);
  size_t numSamples = saved > 0 ? saved : 1;
  size_t length =
    source.numRows * numSamples * dbarts_sampler_numChains(sampler);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) length));
  int predicted = dbarts_sampler_predict(sampler, &source, NULL, 0,
                                         REAL(result));
  UNPROTECT(1);
  return predicted ? result : R_NilValue;
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
  int predicted = dbarts_sampler_predict(sampler, &source, NULL, 0,
                                         REAL(result));
  UNPROTECT(1);
  return predicted ? result : R_NilValue;
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
  int predicted =
    dbarts_sampler_predict(sampler, &source,
                           Rf_isNull(offsetExpr) ? NULL : REAL(offsetExpr), 0,
                           REAL(result));
  UNPROTECT(1);
  return predicted ? result : R_NilValue;
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
  int predicted =
    dbarts_sampler_predict(sampler, &source, NULL,
                           (size_t) Rf_asInteger(nThreadsExpr), REAL(result));
  UNPROTECT(1);
  return predicted ? result : R_NilValue;
}


/* COPY-ON-SET: the entry is handed a buffer THIS consumer owns, and that
 * buffer is overwritten with clobber before returning. A setter that retained
 * the pointer would leave the sampler conditioned on the clobber values; one
 * that copied leaves it conditioned on y, which is what the R side's two runs
 * compare. R_alloc is the right storage: it outlives the entry call and is
 * released when this .Call returns, by which point the sampler holds its own
 * copy and nothing points here. */
SEXP capi_set_response_clobber(SEXP ptrExpr, SEXP yExpr, SEXP clobberExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  double* buffer = (double*) R_alloc(n, sizeof(double));
  int status;
  memcpy(buffer, REAL(yExpr), n * sizeof(double));
  status = dbarts_sampler_setResponse(sampler, buffer, FALSE);
  memcpy(buffer, REAL(clobberExpr), n * sizeof(double));
  return Rf_ScalarInteger(status);
}

/* the offset twin, on the same contract */
SEXP capi_set_offset_clobber(SEXP ptrExpr, SEXP offsetExpr, SEXP clobberExpr) {
  dbarts_sampler* sampler = samplerFromExpr(ptrExpr);
  size_t n = dbarts_sampler_numObservations(sampler);
  double* buffer = (double*) R_alloc(n, sizeof(double));
  int status;
  memcpy(buffer, REAL(offsetExpr), n * sizeof(double));
  status = dbarts_sampler_setOffset(sampler, buffer, FALSE);
  memcpy(buffer, REAL(clobberExpr), n * sizeof(double));
  return Rf_ScalarInteger(status);
}

/* the early release: the engine goes, the R object keeps the pointer and
 * reads dead through its own methods. Called twice by the R side, since a
 * second destroy is the one call a destroyed handle still takes. */
SEXP capi_destroy(SEXP ptrExpr) {
  dbarts_sampler_destroy(samplerFromExpr(ptrExpr));
  return R_NilValue;
}
