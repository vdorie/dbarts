#ifndef DBARTS_DBARTS_H
#define DBARTS_DBARTS_H

/// \file dbarts.h
/// Flat C interface to the dbarts sampler, for packages that drive BART as
/// a conditional model inside a larger sampler (LinkingTo: dbarts).
///
/// Every function is registered with R_RegisterCCallable under its own name.
/// A consumer has two ways to reach an entry point. It can look each one up by
/// hand with R_GetCCallable("dbarts", "<name>") and cast the DL_FUNC to the
/// matching signature, or - the supported path - it can define DBARTS_USE_STUBS
/// before including this header, which replaces the prototypes below with
/// same-name static inline stubs that resolve and cache the pointer on first
/// call. The stubs are generated from DBARTS_C_API_LIST, the single source of
/// truth for the surface, so a rebuild always re-derives the consumer's call
/// types from this header and a stale hand-rolled signature cannot drift.
///
/// The version is two components: check dbarts_apiMajorVersion() ==
/// DBARTS_C_API_MAJOR && dbarts_apiMinorVersion() >= DBARTS_C_API_MINOR at load
/// time (major = incompatible change, minor = additive). dbarts_apiHash() adds
/// an exact signature token beside that handshake, for a consumer that wants
/// lockstep rather than compatibility. The interface only ever grows: names and
/// signatures below are stable and function additions arrive under new names (a
/// minor bump), while dbarts_results grows in place by appending fields - its
/// leading structSize keeps callers compiled against an older layout safe, and
/// the caller-filled structs below (dbarts_predictor_source,
/// dbarts_forest_calibration) carry the same leading member for the same
/// reason, in the other direction.
///
/// Contracts common to all entry points:
/// - Invalid arguments raise R errors (Rf_error), which longjmp through the
///   caller's frames; call from contexts that are safe to unwind.
/// - Functions that draw (creation, run, the prior samplers, predictor
///   updates, dbarts_drawLatents) manage R's RNG state internally and must be
///   called from the main R thread. Do not wrap them in a GetRNGstate/
///   PutRNGstate bracket that spans your own draws through R's API.
///   dbarts_sampler_predict and dbarts_sampler_setTestPredictors are
///   main-R-thread-only for a separate reason: both are R_alloc-backed
///   internally, and R_alloc is unsafe off that thread.
/// - Creation preserves the data specification object against garbage
///   collection for the sampler's lifetime, and the engine borrows its
///   predictors only to quantize them into owned codes at construction; the
///   preserved data object is thereafter the sampler's own predictor GC anchor
///   and the call-time raw source for saved-tree replay and state restore.
///   Raw-array setters borrow for the call's duration only: setResponse,
///   setOffset, setWeights and setForestWeights keep the array alive until
///   replaced; the predictor setters re-quantize the borrow and retain no
///   pointer, as setActiveRows consumes its mask, and setForestBasis copies
///   the basis columns it is handed.
/// - Matrices are column-major. Result and prediction layouts put samples
///   and then chains in trailing dimensions.

#include <stddef.h> // size_t
#include <stdint.h> // uint32_t

// imports Rinternals.h for SEXP while doing the least to pollute namespaces
#include <Rversion.h>
#if R_VERSION >= R_Version(3, 6, 2)
#  define USE_FC_LEN_T
#endif
#ifndef R_NO_REMAP
#  define DBARTS_UNMAP_R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <Rinternals.h>
#ifdef DBARTS_UNMAP_R_NO_REMAP
#  undef R_NO_REMAP
#  undef DBARTS_UNMAP_R_NO_REMAP
#endif
#undef USE_FC_LEN_T

/// The C ABI version, two components. major changes are
/// incompatible; minor changes are additive-only. The safe consumer handshake
/// is major-equality with a minor floor:
///   dbarts_apiMajorVersion() == DBARTS_C_API_MAJOR &&
///   dbarts_apiMinorVersion() >= DBARTS_C_API_MINOR
/// The constants become a compatibility contract at the first release: no
/// version of this API has shipped yet, so whatever they read then simply IS
/// the initial contract, and they do not move before it.
#define DBARTS_C_API_MAJOR 1
#define DBARTS_C_API_MINOR 0

/// FNV-1a hash of the stringized DBARTS_C_API_LIST signatures, baked here and
/// static_assert'd against a recomputation in dbarts's own C++ build. Any
/// signature change in the list fails dbarts's compile until this literal is
/// re-baked; that re-bake is the mechanical acknowledgment of an ABI change.
/// Plain hex literal, usable from C.
///
/// dbarts_apiHash() == DBARTS_C_API_HASH is the exact-signature check: while
/// the constants above do not move, it is the only runtime signal that
/// distinguishes a consumer binary built against a different header from one
/// built against this one, so a lockstep consumer checks it alongside the
/// major/minor handshake.
///
/// The token covers the entry-point SIGNATURES - return types, names and
/// parameter lists - and NOT the layout of the structs those signatures name:
/// two headers differing only in a dbarts_results field's type hash the same.
/// The layout contract is each struct's leading structSize plus the exact
/// offset locks dbarts's own build asserts, so a struct layout change is not
/// self-detecting and is announced to consumers by hand.
#define DBARTS_C_API_HASH 0x85bd1ef04beb3848ULL

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// ABI types, shared by both the prototype view and the stub view below. Every
// type that crosses the ABI is defined here (dbarts_results and the callback
// are inline; control/model/data/state cross as SEXP) so that the single-source
// list and the compile-time token see the whole surface.
// ---------------------------------------------------------------------------

/// Opaque sampler handle.
typedef struct dbarts_sampler_t dbarts_sampler;

/// Caller-owned, growable output buffers for dbarts_sampler_run. The caller
/// MUST set structSize to sizeof(dbarts_results) as the caller compiled it;
/// the library fills only fields whose end offset falls within structSize,
/// so a caller built against an older (smaller) header is never written
/// past. Fields append monotonically below the marked boundary and never
/// reorder across releases; an append after 1.0-0 bumps DBARTS_C_API_MINOR,
/// while a pre-1.0-0 append extends the initial field set and moves no version
/// constant. A field is
/// filled only when both present-by-size and non-null: a null member skips
/// that quantity, and a zero or unset structSize makes dbarts_sampler_run error
/// rather than silently produce no output. k requires a k
/// hyperprior (dbarts_sampler_kIsSampled), varprobs a DART tree prior
/// (dbarts_sampler_usesDart), tau/groupEffects a grouped
/// random-intercept sampler, dispersion a count (nbinom) response, and
/// residualDf a Student-t residual law; each is
/// left untouched otherwise. logLikelihood
/// carries the per-draw training-data log-likelihood for the gaussian,
/// binary, and aft families; aft reports the log density for events and the
/// log survival tail for right-censored observations, and a grouped
/// (random-intercept) sampler composes it over the per-group fits. It is
/// NaN-filled wherever the combined per-observation location is not visible to
/// the response model to score - any sampler whose forests combine through
/// amplitudes, at any forest count, and the multinomial softmax - and skipping
/// it (null or absent-by-size) elides all of its computation.
/// varcount reports the sampler's REPORTED forest, which on a multi-forest
/// model is the prognostic forest (the first): this struct declares no forest
/// count, so the engine writes exactly the numPredictors x numSamples x
/// numChains slab documented below whatever the sampler's forest count is. A
/// caller wanting every forest's split counts drives the sampler from R, whose
/// run channel carries a forest axis.
/// Value-initialize with DBARTS_RESULTS_INIT (sets structSize, zeroes the rest):
///   dbarts_results results = DBARTS_RESULTS_INIT;
typedef struct dbarts_results_t {
  size_t structSize;  ///< caller sets to sizeof(dbarts_results)
  double* sigma;      ///< numSamples x numChains
  double* train;      ///< numObservations x numSamples x numChains
  double* test;       ///< numTestObservations x numSamples x numChains
  uint32_t* varcount; ///< numPredictors x numSamples x numChains
  double* k;          ///< numSamples x numChains
  double* varprobs;   ///< numPredictors x numSamples x numChains
  double* tau;        ///< numSamples x numChains
  double* groupEffects; ///< numGroups x numSamples x numChains
  double* logLikelihood; ///< numObservations x numSamples x numChains
  double* dispersion;    ///< numSamples x numChains, the nbinom r per draw
  double* residualDf;    ///< numSamples x numChains, the Student-t nu per draw
  /* 1.0-0 field boundary: every future append goes below this line, never
     above. An append after 1.0-0 bumps DBARTS_C_API_MINOR; a pre-1.0-0 one
     extends the initial field set above it and moves no version constant. */
} dbarts_results;

/// True when the caller's struct (per structSize) actually carries `field`.
/// The sizeof operand is unevaluated, so this never dereferences past the
/// caller's buffer. Every size-first struct here reads its optional members
/// through it, whether the library writes them (dbarts_results,
/// dbarts_forest_calibration) or reads them (dbarts_predictor_source).
#define DBARTS_HAS_FIELD(type, ptr, field) \
  ((ptr)->structSize >= offsetof(type, field) + sizeof((ptr)->field))

/// The dbarts_results spelling of DBARTS_HAS_FIELD.
#define DBARTS_RESULTS_HAS(r, field) DBARTS_HAS_FIELD(dbarts_results, r, field)

/// Value-initializer: sets structSize (the leading member, offset 0) and zeroes
/// the field pointers. Prefer it to hand-setting structSize - dbarts_sampler_run
/// rejects a zero structSize rather than silently producing no output.
///   dbarts_results results = DBARTS_RESULTS_INIT;
#define DBARTS_RESULTS_INIT { sizeof(dbarts_results) }

/// A predictor column's type. Ordinal columns are cut on their values;
/// categorical ones carry 0-based category codes.
typedef enum {
  DBARTS_COLUMN_ORDINAL = 0,
  DBARTS_COLUMN_CATEGORICAL = 1
} dbarts_column_type;

/// A borrowed, self-describing view of predictor values: what the four
/// predictor entries take instead of a bare pointer. The caller MUST set
/// structSize to sizeof(dbarts_predictor_source) as it compiled it; the library
/// reads only fields whose end offset falls within structSize, so a caller
/// built against an older (smaller) header is never read past, and a zero
/// structSize is an error rather than a source of silent nulls. Fields append
/// monotonically below the marked boundary and never reorder. Nothing here is
/// retained: the entries quantize or replay the values during the call.
///
/// numRows x numColumns is the shape the argument declares of ITSELF, which is
/// what a caller's own width must agree with - the entries refuse a numColumns
/// that disagrees with the sampler rather than reading whatever lies past the
/// caller's matrix.
///
/// Storage is dense, compressed-column (CSC), or a mix. Without a map
/// (columnSources null) the view is the plain column-major block denseValues
/// points at. With one, column j reads dense column columnSources[j] when that
/// is >= 0 and CSC column ~columnSources[j] when it is < 0; the CSC triple is
/// the usual (column pointers of length numCscColumns + 1, row indices,
/// values), and a CSC column's absent rows read its declared reference code
/// when the sampler holds that column categorical, 0 otherwise.
///
/// Value-initialize with DBARTS_PREDICTOR_SOURCE_INIT (sets structSize, zeroes
/// the rest), or build the dense case with dbarts_dense_predictor_source():
///   dbarts_predictor_source x = DBARTS_PREDICTOR_SOURCE_INIT;
typedef struct dbarts_predictor_source_t {
  size_t structSize;              ///< caller sets to sizeof(dbarts_predictor_source)
  size_t numRows;
  size_t numColumns;
  const double* denseValues;      ///< column-major, numRows x numColumns
  size_t numCscColumns;           ///< CSC columns; 0 when there is no CSC part
  const int* cscColumnPointers;   ///< length numCscColumns + 1
  const int* cscRowIndices;
  const double* cscValues;
  const int32_t* columnSources;   ///< NULL = identity; >= 0 dense col; < 0 CSC col ~v
  const int32_t* columnTypes;     ///< dbarts_column_type per column; NULL = all ordinal
  const uint32_t* categoryCounts; ///< declared K per column; NULL/0 = infer
  const int32_t* referenceCodes;  ///< per column; < 0 = declared none
  /* 1.0-0 field boundary: appends go below, never above. */
} dbarts_predictor_source;

/// Value-initializer: sets structSize (the leading member, offset 0) and zeroes
/// everything else, which reads as "dense block, no map, nothing declared".
#define DBARTS_PREDICTOR_SOURCE_INIT \
  { sizeof(dbarts_predictor_source), 0, 0, NULL, 0, NULL, NULL, NULL, NULL, \
    NULL, NULL, NULL }

/// The dense spelling: a plain column-major numRows x numColumns block, which
/// is the shape the predictor entries took before this struct existed.
static inline dbarts_predictor_source
dbarts_dense_predictor_source(const double* values, size_t numRows,
                              size_t numColumns) {
  dbarts_predictor_source source = DBARTS_PREDICTOR_SOURCE_INIT;
  source.numRows = numRows;
  source.numColumns = numColumns;
  source.denseValues = values;
  return source;
}

/// The leaf model a sampler's forests carry, which qualifies what a reported
/// prior sd means (see dbarts_sampler_forestCalibration).
typedef enum {
  DBARTS_LEAF_CONSTANT = 0,
  DBARTS_LEAF_MONOTONE = 1,
  DBARTS_LEAF_LINEAR   = 2,
  DBARTS_LEAF_GP       = 3
} dbarts_leaf_model;

/// Caller-owned output buffers for dbarts_sampler_forestCalibration, on the
/// dbarts_results contract: set structSize; EVERY member is a pointer and is
/// filled only when both present-by-size and non-null, each over numChains; a
/// zero structSize errors. Fields append below the marked boundary and never
/// reorder. All quantities are in RESPONSE units (the family's latent units
/// where the response is not rescaled).
///   dbarts_forest_calibration calibration = DBARTS_FOREST_CALIBRATION_INIT;
typedef struct dbarts_forest_calibration_t {
  size_t  structSize;      ///< caller sets to sizeof(dbarts_forest_calibration)
  double* priorScale;      ///< numChains; forest-total prior sd at k = 1
  double* priorSd;         ///< numChains; priorScale / k at the current k
  double* priorMean;       ///< numChains; prior mean of the forest total
  double* k;               ///< numChains
  double* responseScale;   ///< numChains; internal-to-response multiplier
  double* responseShift;   ///< numChains; internal-to-response offset
  int*    kHasHyperprior;  ///< numChains; THIS FOREST's own k law (not the
                           ///< sampler-wide dbarts_sampler_kIsSampled,
                           ///< which reads the sampler option and
                           ///< disagrees on BCF and multinomial)
  int*    leafModel;       ///< numChains; dbarts_leaf_model, qualifying
                           ///< priorSd and priorMean (see below)
  /* The multi-forest CALIBRATION MAP's decomposition of priorScale, which is
   * factor * s / (divisor * rowNorm) at the family's latent anchor s: NaN on
   * every forest with no map entry (any single-forest or multinomial sampler),
   * so a caller reads "not map-derived" rather than a plausible 1.0. */
  double* amplitudePriorVariance;  ///< numChains; the forest's amplitude prior
                                   ///< variance, or NaN where the forest
                                   ///< carries the half-Cauchy scale mixture
                                   ///< instead - the two are EXCLUSIVE
  double* amplitudePriorScale;     ///< numChains; the half-Cauchy median, or
                                   ///< NaN on a fixed-variance forest
  double* nodeScaleFactor;         ///< numChains; NaN once a state install has
                                   ///< brought a leaf scale this map did not
                                   ///< derive, until setForestBasis re-imposes
                                   ///< it
  double* nodeScaleDivisor;        ///< numChains; NaN on the same rule
  double* basisRowNorm;            ///< numChains; median nonzero row norm of
                                   ///< the forest's basis IN FORCE
  /* 1.0-0 field boundary: every future append goes below this line, never
     above. An append after 1.0-0 bumps DBARTS_C_API_MINOR; a pre-1.0-0 one
     extends the initial field set above it and moves no version constant. */
} dbarts_forest_calibration;

/// Value-initializer: sets structSize (the leading member, offset 0) and
/// zeroes the field pointers, so a caller fills in only what it wants.
#define DBARTS_FOREST_CALIBRATION_INIT { sizeof(dbarts_forest_calibration) }

/// Per-sweep conditioning callback. dbarts_sampler_run invokes it on the
/// calling thread before every sweep - each of the (numBurnIn + numSamples) x
/// numThin iterations - passing the chain index, the 0-based sweep counter,
/// and 1 while the sweep is discarded burn-in. Mutate conditioning state
/// (dbarts_sampler_setSigma, dbarts_sampler_setOffset, ...) from inside it to
/// reproduce a setState-then-run(0, 1) loop exactly, at no per-sweep R round
/// trip; return 0 to stop the run early (the results filled so far are then
/// undefined). It fires before dbarts_sampler_setSigma's held sigma or the
/// gaussian sigma draw enters the sweep, so a value set here conditions it.
typedef int (*dbarts_sampler_callback)(void* userData, dbarts_sampler* sampler,
                                       size_t chainIndex, size_t sweepIndex,
                                       int isBurnIn);

// ---------------------------------------------------------------------------
// The single source of truth for the entry-point surface. Each
// entry is X(returnType, name, (parameterList), (argumentList)): the parameter
// list carries names so it also spells the forwarding stub's signature, and the
// argument list forwards those names. Registration (R_interface.cpp), the
// consumer stubs below, the provider-side binding asserts, and the compile-time
// token (C_interface.cpp) are all expansions of this one list, so a signature
// stated here is the only place it is stated. The readable Doxygen prototypes
// kept in the #else branch below are compile-time bound to this list in
// dbarts's own build, so any drift between them fails dbarts's compile.
// ---------------------------------------------------------------------------
#define DBARTS_C_API_LIST(X) \
  X(int, dbarts_apiMajorVersion, (void), ()) \
  X(int, dbarts_apiMinorVersion, (void), ()) \
  X(dbarts_sampler*, dbarts_sampler_create, \
    (SEXP control, SEXP model, SEXP data, const char* family), \
    (control, model, data, family)) \
  X(void, dbarts_sampler_destroy, (dbarts_sampler* sampler), (sampler)) \
  X(void, dbarts_sampler_run, \
    (dbarts_sampler* sampler, size_t numBurnIn, size_t numSamples, \
     dbarts_results* results), \
    (sampler, numBurnIn, numSamples, results)) \
  X(void, dbarts_sampler_sampleTreesFromPrior, (dbarts_sampler* sampler), \
    (sampler)) \
  X(void, dbarts_sampler_sampleNodeParametersFromPrior, \
    (dbarts_sampler* sampler), (sampler)) \
  X(void, dbarts_sampler_setResponse, \
    (dbarts_sampler* sampler, const double* y, int updateScale), \
    (sampler, y, updateScale)) \
  X(void, dbarts_sampler_setOffset, \
    (dbarts_sampler* sampler, const double* offset, int updateScale), \
    (sampler, offset, updateScale)) \
  X(void, dbarts_sampler_setWeights, \
    (dbarts_sampler* sampler, const double* weights), (sampler, weights)) \
  X(void, dbarts_sampler_setSigma, \
    (dbarts_sampler* sampler, double sigma), (sampler, sigma)) \
  X(void, dbarts_sampler_setCallback, \
    (dbarts_sampler* sampler, dbarts_sampler_callback callback, \
     void* userData), \
    (sampler, callback, userData)) \
  X(int, dbarts_sampler_getLatents, \
    (const dbarts_sampler* sampler, double* out), (sampler, out)) \
  X(int, dbarts_sampler_setPredictor, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* x, \
     int forceUpdate, int updateCutPoints), \
    (sampler, x, forceUpdate, updateCutPoints)) \
  X(int, dbarts_sampler_updatePredictor, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* x, \
     const size_t* columns, size_t numColumns, int forceUpdate, \
     int updateCutPoints), \
    (sampler, x, columns, numColumns, forceUpdate, updateCutPoints)) \
  X(void, dbarts_sampler_setTestPredictors, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest), \
    (sampler, xTest)) \
  X(void, dbarts_sampler_setTestOffset, \
    (dbarts_sampler* sampler, const double* offsetTest), \
    (sampler, offsetTest)) \
  X(void, dbarts_sampler_predict, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest, \
     const double* offsetTest, double* out), \
    (sampler, xTest, offsetTest, out)) \
  X(void, dbarts_sampler_setTreeStorage, \
    (dbarts_sampler* sampler, int keepTrees, size_t numSamplesToStore), \
    (sampler, keepTrees, numSamplesToStore)) \
  X(SEXP, dbarts_sampler_getTrees, \
    (dbarts_sampler* sampler, const size_t* chainIndices, \
     size_t numChainIndices, const size_t* sampleIndices, \
     size_t numSampleIndices, const size_t* treeIndices, \
     size_t numTreeIndices, int useLiveTrees, size_t forest), \
    (sampler, chainIndices, numChainIndices, sampleIndices, numSampleIndices, \
     treeIndices, numTreeIndices, useLiveTrees, forest)) \
  X(void, dbarts_sampler_printTrees, \
    (dbarts_sampler* sampler, const size_t* chainIndices, \
     size_t numChainIndices, const size_t* sampleIndices, \
     size_t numSampleIndices, const size_t* treeIndices, \
     size_t numTreeIndices, size_t forest), \
    (sampler, chainIndices, numChainIndices, sampleIndices, numSampleIndices, \
     treeIndices, numTreeIndices, forest)) \
  X(SEXP, dbarts_sampler_storeState, (dbarts_sampler* sampler), (sampler)) \
  X(void, dbarts_sampler_setState, \
    (dbarts_sampler* sampler, SEXP state), (sampler, state)) \
  X(void, dbarts_sampler_setNumThreads, \
    (dbarts_sampler* sampler, size_t numThreads), (sampler, numThreads)) \
  X(void, dbarts_sampler_setNumThin, \
    (dbarts_sampler* sampler, size_t numThin), (sampler, numThin)) \
  X(void, dbarts_sampler_setVerbose, \
    (dbarts_sampler* sampler, int verbose, uint32_t printEvery), \
    (sampler, verbose, printEvery)) \
  X(size_t, dbarts_sampler_numObservations, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numPredictors, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numTestObservations, \
    (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numChains, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numTrees, \
    (const dbarts_sampler* sampler, size_t forest), (sampler, forest)) \
  X(size_t, dbarts_sampler_numSavedSamples, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(int, dbarts_sampler_kIsSampled, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(int, dbarts_sampler_usesDart, (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numForests, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(int, dbarts_sampler_setForestBasis, \
    (dbarts_sampler* sampler, size_t forest, const double* basis, \
     size_t numColumns), \
    (sampler, forest, basis, numColumns)) \
  X(int, dbarts_sampler_forestFits, \
    (const dbarts_sampler* sampler, size_t forest, double* out), \
    (sampler, forest, out)) \
  X(size_t, dbarts_sampler_numForestAmplitudes, \
    (const dbarts_sampler* sampler, size_t forest), (sampler, forest)) \
  X(int, dbarts_sampler_forestAmplitudes, \
    (const dbarts_sampler* sampler, size_t forest, double* out), \
    (sampler, forest, out)) \
  X(uint64_t, dbarts_apiHash, (void), ()) \
  X(int, dbarts_sampler_setForestWeights, \
    (dbarts_sampler* sampler, size_t forest, const double* weights), \
    (sampler, forest, weights)) \
  X(int, dbarts_sampler_forestCalibration, \
    (const dbarts_sampler* sampler, size_t forest, \
     dbarts_forest_calibration* out), \
    (sampler, forest, out)) \
  X(int, dbarts_sampler_setForestPriorScale, \
    (dbarts_sampler* sampler, size_t forest, double priorScale), \
    (sampler, forest, priorScale)) \
  X(int, dbarts_sampler_setActiveRows, \
    (dbarts_sampler* sampler, const double* active), (sampler, active)) \
  X(int, dbarts_sampler_dispersion, \
    (const dbarts_sampler* sampler, double* out), (sampler, out)) \
  X(void, dbarts_drawLatents, \
    (const char* family, size_t numObservations, const double* fit, \
     const double* y, const double* weights, const double* offset, \
     double sigma, double dispersion, const double* cutpoints, \
     size_t numCutpoints, double df, double* out), \
    (family, numObservations, fit, y, weights, offset, sigma, dispersion, \
     cutpoints, numCutpoints, df, out)) \
  X(void, dbarts_workingResponse, \
    (const char* family, size_t numObservations, const double* latent, \
     const double* y, const double* weights, const double* offset, \
     double dispersion, double* out), \
    (family, numObservations, latent, y, weights, offset, dispersion, out))

/// One stringized "returnType name(parameterList);" per list entry, adjacent
/// string literals that concatenate into the full declaration text the
/// compile-time token hashes.
#define DBARTS_API_STRINGIZE(ret, name, params, args) #ret " " #name #params ";"
#define DBARTS_C_API_DECLS DBARTS_C_API_LIST(DBARTS_API_STRINGIZE)

#ifdef DBARTS_USE_STUBS

// Same-name cached-pointer forwarders generated from the list, one per entry,
// in place of the extern prototypes (the xts inst/include/xtsAPI.h idiom). The
// first call resolves the symbol through R_GetCCallable and caches it in a
// block-scope static; later calls forward directly. A consumer defines
// DBARTS_USE_STUBS to opt in, and then never restates a signature.
#include <R_ext/Rdynload.h> // R_GetCCallable, DL_FUNC

// void-return detection: ISO C forbids `return <expr>;` in a void function, so
// a void stub must forward without `return`. These helpers pick the right body
// from the entry's return type; all are #undef'd right after the expansion so
// none leak into the consumer's macro namespace.
#define DBARTS_CAT(a, b) DBARTS_CAT_(a, b)
#define DBARTS_CAT_(a, b) a##b
#define DBARTS_CHECK_N(x, n, ...) n
#define DBARTS_CHECK(...) DBARTS_CHECK_N(__VA_ARGS__, 0,)
#define DBARTS_PROBE(x) x, 1,
#define DBARTS_VOID_void DBARTS_PROBE(~)
#define DBARTS_IS_VOID(ret) DBARTS_CHECK(DBARTS_CAT(DBARTS_VOID_, ret))
#define DBARTS_IIF(c) DBARTS_CAT(DBARTS_IIF_, c)
#define DBARTS_IIF_0(t, f) f
#define DBARTS_IIF_1(t, f) t

// The resolution retypes DL_FUNC (void *(*)(void)) to the entry's own
// signature through void (*)(void), the universal function-pointer type the
// -Wextra cast diagnostics exempt (gcc -Wcast-function-type, clang
// -Wcast-function-type-mismatch): cast directly, and every consumer
// translation unit warns once per list entry. The round trip is a retype
// only and generates no code.
#define DBARTS_API_STUB(ret, name, params, args) \
  static inline ret name params { \
    static ret (*dbarts_stub_fn) params = NULL; \
    if (dbarts_stub_fn == NULL) \
      dbarts_stub_fn = \
        (ret (*) params) (void (*)(void)) R_GetCCallable("dbarts", #name); \
    DBARTS_IIF(DBARTS_IS_VOID(ret))(dbarts_stub_fn args;, \
                                    return dbarts_stub_fn args;) \
  }
DBARTS_C_API_LIST(DBARTS_API_STUB)

#undef DBARTS_API_STUB
#undef DBARTS_IIF_1
#undef DBARTS_IIF_0
#undef DBARTS_IIF
#undef DBARTS_IS_VOID
#undef DBARTS_VOID_void
#undef DBARTS_PROBE
#undef DBARTS_CHECK
#undef DBARTS_CHECK_N
#undef DBARTS_CAT_
#undef DBARTS_CAT

#else // !DBARTS_USE_STUBS: the readable, Doxygen-documented prototypes.

/// Returns DBARTS_C_API_MAJOR of the installed package (incompatible-change
/// component of the version; equal to the caller's for a usable library).
int dbarts_apiMajorVersion(void);
/// Returns DBARTS_C_API_MINOR of the installed package (additive component;
/// at least the caller's for a usable library).
int dbarts_apiMinorVersion(void);
/// Returns the installed package's DBARTS_C_API_HASH: the FNV-1a token over
/// the entry-point signatures. Equality with the caller's DBARTS_C_API_HASH
/// says the two were built from the same declared surface - the lockstep
/// check, and the only runtime signal that moves while the version constants
/// do not. It sees signatures, never struct layout (see DBARTS_C_API_HASH).
uint64_t dbarts_apiHash(void);

/// Creates a sampler from the R specification objects (dbartsControl,
/// dbartsModel, dbartsData). family selects the response model for binary
/// responses, on a single forest or on K of them: "" or "probit" give probit
/// latents, "logistic" the
/// Polya-Gamma sampler; continuous responses are gaussian and accept "" or
/// "gaussian", or "aft" for an accelerated failure time (log-normal) survival
/// fit, in which case the response holds log survival/censoring times and a
/// numeric status vector (1 = event, 0 = right-censored) is read from the
/// control's "bartcore.survival" attribute. A verbose control prints the
/// initial summary here.
///
/// The K-forest amplitude family (docs/design/multiplier-combiner.md) is
/// created through this same entry point. Each forest carries its own basis,
/// whose row contracts with that forest's own amplitude vector into the scalar
/// the forest's fit is multiplied by, so the location is
/// sum_f dot(a_f, B_f(i, .)) * forestFits(f)[i]; a Bayesian causal forest
/// (docs/design/bcf.md) is the K = 2 instance, forest 0 prognostic over an
/// implicit intercept and forest 1 over a treatment indicator. The R
/// specification dbartsSpec(data, control, forests = list(forest(), forest(
/// basis = ~ factor(z)))) puts each declared forest's basis columns on the data
/// object and those forests' configuration on the control, and creation reads
/// both halves (either alone is an error); a K-length forests list builds K
/// forests, and dbartsData(bases = ) carries a numeric basis directly. Such a
/// sampler reports numForests == K, swaps any forest's basis with
/// dbarts_sampler_setForestBasis, reads each forest's surface with
/// dbarts_sampler_forestFits and its ragged amplitude block with
/// dbarts_sampler_numForestAmplitudes and dbarts_sampler_forestAmplitudes,
/// takes a scale-pinned response, offset or weight swap, and refuses the whole
/// test surface (setTestPredictors, setTestOffset, predict), whose blend is
/// undefined without a test basis. Gaussian, probit and logistic responses;
/// the other families are refused at creation, naming what each is missing.
/// Under a latent family the combined location IS the index, on the link's own
/// fixed scale - so every forest's prior scale is stated in latent sd units,
/// sigma is pinned, and there is no response transform to rescale.
dbarts_sampler* dbarts_sampler_create(SEXP control, SEXP model, SEXP data,
                                      const char* family);
void dbarts_sampler_destroy(dbarts_sampler* sampler);

/// Runs numBurnIn discarded then numSamples recorded iterations per chain,
/// recording into results (which may be null when numSamples is 0). Set
/// results->structSize before calling; fields whose end offset exceeds it
/// are skipped (never read, never written). Thinning applies within
/// recorded iterations at the rate the control set.
void dbarts_sampler_run(dbarts_sampler* sampler, size_t numBurnIn,
                        size_t numSamples, dbarts_results* results);
void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler);
void dbarts_sampler_sampleNodeParametersFromPrior(dbarts_sampler* sampler);

/// Registers callback for subsequent runs, or clears it when null; userData is
/// passed back unchanged and its lifetime is the caller's. Raises an error
/// while chains would run on worker threads (numThreads > 1 and numChains > 1),
/// which must never call into R. Inline multi-chain runs execute chains
/// sequentially - chain 0 completes all its sweeps before chain 1 starts - so
/// the callback conditions each chain to completion in turn and cannot see one
/// chain's progress while advancing another. (See dbarts_sampler_callback.)
void dbarts_sampler_setCallback(dbarts_sampler* sampler,
                                dbarts_sampler_callback callback,
                                void* userData);

/// y has numObservations values, which must lie in the family's support: 0/1
/// for probit and logistic, an integer category index in [1, K] for ordinal, a
/// finite non-negative integer count no larger than 1e6 for nbinom (the
/// dispersion grid's count histogram is sized from the largest count, so a
/// larger one allocates without bound). Out-of-support values are an
/// error, as they are at creation; gaussian and aft (log survival times)
/// constrain nothing. updateScale re-derives the internal response transform
/// from the new response, as dbarts_sampler_setOffset's argument does (gaussian
/// only); pass false once burnt in so fits stay comparable. true is refused on
/// any multi-forest sampler, at any forest count, whose per-forest leaf
/// calibrations are stated against the transform it was built with, and on a
/// heteroscedastic one, whose variance forest is calibrated the same way. The
/// swap itself is refused outright on a coupling that caches per-forest state
/// across sweeps rather than re-deriving it.
void dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y,
                                int updateScale);
/// offset has numObservations values or is null to remove. updateScale
/// rescales the internal response transform to the offset-adjusted range
/// (gaussian only); pass false once burnt in so fits stay comparable. A
/// multi-forest sampler, at any forest count, or a heteroscedastic one refuses
/// true (see setResponse).
void dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                              int updateScale);
/// weights has numObservations values. Post-creation weight changes are for
/// gaussian responses only - the binary families read a weight as a copy count
/// fixed at creation - but a logistic sampler accepts positive integer counts
/// AT creation, at any forest count. There is no scale to pin, so a
/// multi-forest sampler takes this as it stands, at any forest count, provided
/// its coupling admits the conduit at all.
void dbarts_sampler_setWeights(dbarts_sampler* sampler,
                               const double* weights);
/// Holds the residual standard deviation at sigma (original response scale)
/// until the next call or gaussian draw; the Gibbs conditioning hook.
void dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma);
/// Copies the current draw of the augmentation variable (numObservations x
/// numChains) into out. Returns 1, or 0 without touching out when the family
/// augments nothing - a plain gaussian response, and a multinomial one.
///
/// WHAT the variable is depends on the family and is not uniform. A LOCATION,
/// on the sampler's own latent scale, for probit (the truncated normal z),
/// ordinal (the same z under the cut points) and aft (the imputed log survival
/// time): a host regresses on these directly. A PRECISION, one per
/// observation, for logistic and nbinom (the Polya-Gamma omega) and for a
/// Student-t residual distribution (the scale-mixing lambda): these WEIGHT a
/// working response and are not on the response scale at all. Note the last
/// case - a sampler whose family is gaussian but whose residual distribution
/// is Student-t DOES report latents, and they are precisions.
int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out);

/// Copies the dispersion r in force on each chain (numChains values) into out:
/// the same scalar dbarts_results::dispersion records once per kept draw, read
/// mid-sweep and without serializing state, so a host composing a count block
/// reads it between sweeps instead of through dbarts_sampler_storeState.
/// Returns 1, or 0 without touching out when the response family carries no
/// dispersion - every family but nbinom - so a caller tests the channel rather
/// than a filler value. Under a fixed dispersion it repeats the value the
/// sampler was created with; otherwise it is that sweep's grid draw.
int dbarts_sampler_dispersion(const dbarts_sampler* sampler, double* out);

/// Replaces the full predictor matrix from a borrowed source: x declares
/// numRows == numObservations and numColumns == numPredictors, and any source
/// disagreeing with either is refused rather than read to the sampler's own
/// width. Transactional: returns 1 when every tree still has non-empty leaves
/// under the new values (or forceUpdate), 0 after rolling back. Errors when cut
/// points would invalidate existing splits and updateCutPoints is false.
///
/// A MUTATION reads its source through one dense materialization, because
/// every mutation kernel indexes values column-major: CSC storage buys the
/// caller a uniform argument, an explicit shape and the validation below, not
/// resident sparse mutation. The test-side entries below are where compressed
/// storage stays compressed.
int dbarts_sampler_setPredictor(dbarts_sampler* sampler,
                                const dbarts_predictor_source* x,
                                int forceUpdate, int updateCutPoints);
/// Replaces numColumns columns, 0-indexed by columns; x declares those same
/// numColumns columns, in ARGUMENT order, over numObservations rows. Same
/// transaction contract as setPredictor.
int dbarts_sampler_updatePredictor(dbarts_sampler* sampler,
                                   const dbarts_predictor_source* x,
                                   const size_t* columns, size_t numColumns,
                                   int forceUpdate, int updateCutPoints);

/// Installs test data from a borrowed source declaring numPredictors columns
/// over its own numRows, or null to remove test data (clearing any test
/// offset). A CSC-backed source is consumed as it stands - no dense
/// materialization - except that a designated leaf covariate column must be
/// dense, which is an error naming the repair. Refused on any sampler whose
/// test blend is undefined - the predicate is the blend, not the forest count,
/// so a multi-forest model whose test location IS defined passes through (see
/// dbarts_sampler_create).
void dbarts_sampler_setTestPredictors(dbarts_sampler* sampler,
                                      const dbarts_predictor_source* xTest);
/// offsetTest has numTestObservations values or is null to remove.
void dbarts_sampler_setTestOffset(dbarts_sampler* sampler,
                                  const double* offsetTest);

/// Fits for new data on the original response scale (binary families give
/// the latent scale), from a borrowed source declaring numPredictors columns
/// over the rows to predict. With tree storage out is xTest->numRows x
/// numSavedSamples x numChains from the saved trees; without, one set per
/// chain from the live trees. offsetTest, when non-null, is added to
/// every sample's fits. A CSC-backed source routes its rows resident, without
/// a dense materialization. Refused on any sampler whose blend is undefined -
/// the predicate is the blend, not the forest count (see
/// dbarts_sampler_create): read the forests separately with
/// dbarts_sampler_forestFits and combine them with the amplitudes.
void dbarts_sampler_predict(dbarts_sampler* sampler,
                            const dbarts_predictor_source* xTest,
                            const double* offsetTest, double* out);

/// Turns saved-tree storage on or off; numSamplesToStore sizes the buffer
/// when on. Turn on for recorded iterations to predict from them later.
void dbarts_sampler_setTreeStorage(dbarts_sampler* sampler, int keepTrees,
                                   size_t numSamplesToStore);
/// A data.frame of tree structure: pre-order rows of ([chain,] [sample,]
/// tree, n, var, value), var 1-based with -1 marking leaves, leaf values on
/// the engine's internal response scale. Indices are 0-based here. Saved
/// trees are read unless useLiveTrees; sampleIndices is ignored when
/// reading live trees. The caller must protect the result.
///
/// forest is a 0-based index in [0, dbarts_sampler_numForests(sampler)), so 0
/// is the only legal value on a single-forest sampler (on a two-forest bcf one,
/// forest 1 is the one carrying the basis - an example, not the range).
/// treeIndices are read against THAT forest's tree count
/// (dbarts_sampler_numTrees(sampler, forest)). An index at or above
/// dbarts_sampler_numForests is an error.
///
/// The n column replays the creation specification's predictors through each
/// saved tree (the engine keeps no predictor matrix); it is left unpopulated
/// for a sparse creation spec, and live trees carry their own counts.
SEXP dbarts_sampler_getTrees(dbarts_sampler* sampler,
                             const size_t* chainIndices,
                             size_t numChainIndices,
                             const size_t* sampleIndices,
                             size_t numSampleIndices,
                             const size_t* treeIndices,
                             size_t numTreeIndices, int useLiveTrees,
                             size_t forest);
/// Prints forest number forest's trees to R's console, on the same index
/// contract as dbarts_sampler_getTrees.
void dbarts_sampler_printTrees(dbarts_sampler* sampler,
                               const size_t* chainIndices,
                               size_t numChainIndices,
                               const size_t* sampleIndices,
                               size_t numSampleIndices,
                               const size_t* treeIndices,
                               size_t numTreeIndices, size_t forest);

/// A serializable R object holding the complete sampler state (trees,
/// parameters, latents, rng); the caller must protect it. Restoring into a
/// sampler created from the same specification continues bitwise
/// identically; use for save/load and predict-after-reload. States are one
/// R list element per chain, so single-chain states may be gathered into a
/// multi-chain restore; a stored generator of a different kind than the
/// destination chain's is then left unrestored (the destination keeps its
/// own stream), which only forfeits cross-kind bitwise continuation.
///
/// dbarts_sampler_setState enforces an encoding floor: a state whose format
/// version predates the running library's minimum readable version is refused
/// (naming both versions) rather than silently misread, and there is no
/// cross-version migration. Blocks it does not recognize by name are skipped
/// rather than errored, so a state written by a newer library restores its
/// shared blocks and leaves the unknown ones to the destination's own draw.
SEXP dbarts_sampler_storeState(dbarts_sampler* sampler);
void dbarts_sampler_setState(dbarts_sampler* sampler, SEXP state);

void dbarts_sampler_setNumThreads(dbarts_sampler* sampler, size_t numThreads);
void dbarts_sampler_setNumThin(dbarts_sampler* sampler, size_t numThin);
void dbarts_sampler_setVerbose(dbarts_sampler* sampler, int verbose,
                               uint32_t printEvery);

size_t dbarts_sampler_numObservations(const dbarts_sampler* sampler);
size_t dbarts_sampler_numPredictors(const dbarts_sampler* sampler);
size_t dbarts_sampler_numTestObservations(const dbarts_sampler* sampler);
size_t dbarts_sampler_numChains(const dbarts_sampler* sampler);
/// Forest number forest's tree count; an index at or above
/// dbarts_sampler_numForests is an error, since a size_t probe carries no
/// refusal a caller could tell from a legitimate answer.
size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler, size_t forest);
/// The saved-tree buffer size, 0 without tree storage; the sample count
/// predict produces.
size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler);
int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler);
int dbarts_sampler_usesDart(const dbarts_sampler* sampler);
/// The number of forests the model fits, always at least 1: K for a sampler
/// built from a K-length forests = specification (a Bayesian causal forest,
/// whose forest 0 is prognostic and forest 1 the treatment effect, is the K = 2
/// instance - see dbarts_sampler_create), K for a multinomial softmax sampler's
/// category forests, and 1 otherwise. A heteroscedastic sampler's variance
/// forest is NOT counted here - it is a separate member, not one of the
/// combined forests - so such a sampler reports 1.
size_t dbarts_sampler_numForests(const dbarts_sampler* sampler);

/// Replaces the basis forest number forest's amplitudes multiply; basis is
/// ROW-major numObservations x numColumns - row i at basis + i * numColumns -
/// and is COPIED, so the caller's array need not outlive the call (a continuous
/// basis cannot be coerced-and-copied incidentally the way a 0/1 indicator can,
/// which is why the contract is stated). Row-major because the contraction with
/// the forest's amplitude vector is the only read the engine makes of it. The
/// one whole-data swap a multi-forest sampler supports - it routes through the
/// combiner rather than rebuilding a forest - and the SOLE route by which any
/// basis changes after creation.
///
/// Two channels, as everywhere here: a CAPABILITY answer is the return value -
/// 0, touching nothing, when the model carries no amplitudes at all or forest
/// names no forest - while a MALFORMED BASIS raises. Any forest takes a basis
/// of any width from 1 up; the values must be finite. A basis that is not one
/// of the canonical shapes (a dense all-ones column, or a complementary 0/1
/// pair) moves that forest onto the general amplitude conditional, which is a
/// model fact and not revertible by a later data swap.
int dbarts_sampler_setForestBasis(dbarts_sampler* sampler, size_t forest,
                                  const double* basis, size_t numColumns);
/// Copies forest number forest's current function values, numObservations x
/// numChains, into out. The values are on the engine's internal response scale
/// (the scale saved-tree leaf values carry), which the stored state's fit.scale
/// maps back to the response scale. Returns 1, or 0 without touching out when
/// forest names no forest.
int dbarts_sampler_forestFits(const dbarts_sampler* sampler, size_t forest,
                              double* out);
/// The number of amplitudes forest number forest carries, so a reader can size
/// its own buffer for a ragged vector: it is that forest's basis column count,
/// whatever the forest count (on a two-forest bcf sampler, 1 for forest 0, the
/// free amplitude a, and 2 for forest 1, the pair b0 and b1); 0 where the model
/// composes its forests without amplitudes. An index at or above
/// dbarts_sampler_numForests is an error (a size_t probe carries no refusal).
size_t dbarts_sampler_numForestAmplitudes(const dbarts_sampler* sampler,
                                          size_t forest);
/// Copies forest number forest's amplitudes,
/// numForestAmplitudes(forest) x numChains, into out. Observation i's
/// internal-scale location is sum_f dot(a_f, B_f(i, .)) * forestFits(f)[i],
/// contracting each forest's amplitude block with the row of the basis last
/// installed on that forest. On a two-forest bcf sampler that is
/// a * forestFits(0)[i] + b_{z_i} * forestFits(1)[i], where a is forest 0's
/// single amplitude and (b0, b1) forest 1's pair against its 0/1 basis.
/// Returns 1, or 0 without touching out when forest names no forest or the
/// model carries no amplitudes.
int dbarts_sampler_forestAmplitudes(const dbarts_sampler* sampler,
                                    size_t forest, double* out);

/// Installs (or clears, at a null vector) a per-forest, per-observation weight
/// on forest number forest: a multiplicative precision factor on that forest's
/// own leaf conditionals, composing with the case weights rather than widening
/// either channel. weights has numObservations finite non-negative values and
/// is BORROWED until replaced - the one place a flat entry's contract differs
/// from the R bridge's, which copies - so the caller owns the array for the
/// sampler's life. Returns 1, or 0 without touching the sampler when the
/// forest coupling admits no such weight or forest names no forest; a
/// non-finite or negative element raises.
///
/// The precision channel and the mean channel stay two entries: this one
/// scales forest f's own leaf conditionals and never enters the combined fit,
/// while dbarts_sampler_setForestBasis scales forest f's contribution to it.
int dbarts_sampler_setForestWeights(dbarts_sampler* sampler, size_t forest,
                                    const double* weights);

/// Fills out with forest number forest's leaf-prior calibration in RESPONSE
/// units, each member over numChains: priorScale (the forest total's prior sd
/// at k = 1, the quantity dbarts_sampler_setForestPriorScale writes), priorSd
/// (priorScale / k at the chain's current k, which moves every sweep under
/// kHasHyperprior while priorScale does not), priorMean, k, the response
/// transform's multiplier and offset, this FOREST's own k law, the leaf
/// model tag, and the multi-forest calibration map's own five: the forest's
/// amplitude prior in whichever of its two EXCLUSIVE spellings it carries
/// (amplitudePriorVariance or amplitudePriorScale, the other NaN), the node
/// scale's factor and divisor, and the basis row norm the map divides out.
/// Those five are NaN on every forest with no map entry - any single-forest
/// or multinomial sampler - and the two node-scale members are NaN once a
/// state install has brought a leaf scale the map did not derive, until
/// dbarts_sampler_setForestBasis re-imposes it. They decompose priorScale as
/// factor * s / (divisor * rowNorm), so the family's latent anchor s is
/// recovered as priorScale * divisor * rowNorm / factor whenever
/// nodeScaleFactor is not NaN.
///
/// prior.scale and prior.sd describe the LEAF-PARAMETER scale of the forest
/// total, which equals the prior sd of f(x) at every x for the constant leaf
/// (DBARTS_LEAF_CONSTANT) and bounds it otherwise. DBARTS_LEAF_LINEAR: a LOWER
/// bound, attained at the standardized covariate origin, larger elsewhere by
/// sqrt(1 + ||z(x)||^2). DBARTS_LEAF_GP: an UPPER bound, attained at rows
/// reproducing a leaf member and on over-cap leaves, elsewhere
/// priorSd^2 c(x)' C^-1 c(x), decaying to 0 as x leaves the leaf's data cloud,
/// where every draw equals priorMean. DBARTS_LEAF_MONOTONE: a LOWER bound in
/// the interior (the realized sd runs a few per cent to about 20% above it),
/// and its priorMean is NOT the prior mean of f(x) under an active constraint,
/// that marginal being skew with an x-dependent mean spanning several priorSd
/// along the constrained axis. priorMean is exact for the constant, linear and
/// gp leaves.
///
/// Returns 1, or 0 without touching out when forest names no forest; errors on
/// a zero structSize.
int dbarts_sampler_forestCalibration(const dbarts_sampler* sampler,
                                     size_t forest,
                                     dbarts_forest_calibration* out);
/// Restates forest number forest's leaf prior on EVERY chain so that the
/// forest total's prior sd at k = 1 is priorScale, in response units. k, the
/// response transform, sigma and the tree prior are untouched; it takes effect
/// on the NEXT sweep and reinterprets no leaf value already drawn; a write
/// reproducing the current internal scale bitwise is skipped, so a
/// read-then-write is inert. To move the prior MEAN, shift the reported fit
/// with dbarts_sampler_setOffset instead. The leaf model qualifies the write
/// exactly as it qualifies the read (see dbarts_sampler_forestCalibration).
///
/// Two channels: a CAPABILITY answer is the return value - 0, touching
/// nothing, when forest names no forest, or a combiner owns this forest's
/// calibration - any sampler whose forests are combined, which is every
/// amplitude-family sampler at any forest count and every multinomial one -
/// while a MALFORMED VALUE, a non-finite or non-positive priorScale, raises.
/// Returns 1 on a write.
int dbarts_sampler_setForestPriorScale(dbarts_sampler* sampler, size_t forest,
                                       double priorScale);

/// Installs (or clears, at a null pointer) a per-observation 0/1 mask saying
/// which rows are in the data set each sweep: an inactive row leaves every
/// sufficient statistic, every family-level parameter update and its own
/// latent draw, while keeping its leaf occupancy and its fitted value. active
/// has numObservations values, each exactly 0.0 or 1.0 (the entry reads the
/// row count from the sampler, so there is no length argument); an all-ones
/// vector is accepted and installs nothing.
///
/// The values are CONSUMED during the call - the entry retains no pointer, so
/// the caller's array is free the moment it returns. Returns 1, or 0 without
/// touching the sampler when the response family implements no mask or an
/// element is neither 0 nor 1 (a fractional value is a weighted likelihood,
/// which the latent families have no coherent form for; use
/// dbarts_sampler_setWeights for that).
int dbarts_sampler_setActiveRows(dbarts_sampler* sampler,
                                 const double* active);

/// Draws the augmentation variable a family's sweep draws, at a caller's fit
/// rather than inside a sampler, and writes numObservations values into out.
/// family is one of "probit", "logistic", "ordinal", "aft", "nbinom" or
/// "student"; WHAT the draw is per family is dbarts_sampler_getLatents's own
/// contract, a location for the first three plus aft and a precision for the
/// rest. fit is the location WITHOUT the offset (the convention
/// dbarts_results::train reports), so the linear predictor is fit + offset and
/// a null offset is zero.
///
/// Each remaining argument belongs to one family's law, which has no default to
/// fall back on: weights are the logistic observation counts (null is unit),
/// sigma the aft and student scale, dispersion the nbinom r (a positive WHOLE
/// number - a fractional one is rounded into a different-shape draw, not
/// refused), cutpoints the ordinal's numCutpoints strictly increasing interior
/// cut points (unordered ones corrupt the draw silently, the rejection loop
/// degrading to its fallback rather than erroring), and df the Student-t
/// degrees of freedom. Omitting what a family's law requires is an error; what
/// it does not read is ignored. Out-of-support y raises, as it does at
/// dbarts_sampler_setResponse - and, by the same rule, y is checked only where
/// the family's support constrains it, so under "aft" and "student" a
/// non-finite y propagates NaN into out.
///
/// The draw comes from R's own stream through a per-call generator, so set.seed
/// reproduces it and it advances the stream for whatever R draws next, while a
/// sampler's chain generators are untouched.
void dbarts_drawLatents(const char* family, size_t numObservations,
                        const double* fit, const double* y,
                        const double* weights, const double* offset,
                        double sigma, double dispersion,
                        const double* cutpoints, size_t numCutpoints, double df,
                        double* out);
/// Turns a drawn latent into the quantity a host regresses on, numObservations
/// values into out: for "logistic" and "nbinom" the Polya-Gamma kappa divided
/// by the drawn precision (which must be positive), for "student" the response
/// itself (its latent weights the row instead), and for the location families
/// the latent - each less the offset. Draws nothing; the arguments carry the
/// meanings dbarts_drawLatents states.
void dbarts_workingResponse(const char* family, size_t numObservations,
                            const double* latent, const double* y,
                            const double* weights, const double* offset,
                            double dispersion, double* out);

#endif // DBARTS_USE_STUBS

#ifdef __cplusplus
}
#endif

#endif // DBARTS_DBARTS_H
