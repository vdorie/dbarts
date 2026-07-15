#ifndef DBARTS_DBARTS_H
#define DBARTS_DBARTS_H

/// \file dbarts.h
/// Flat C interface to the dbarts sampler, for packages that drive BART as
/// a conditional model inside a larger sampler (LinkingTo: dbarts).
///
/// Every function is registered with R_RegisterCCallable under its own
/// name; look entry points up at load time with
/// R_GetCCallable("dbarts", "<name>") and cast the DL_FUNC to the matching
/// signature. Check dbarts_apiVersion() against DBARTS_C_API_VERSION before
/// using the rest. The interface only ever grows: names and signatures
/// below are stable and function additions arrive under new names, while
/// dbarts_results grows in place by appending fields - its leading
/// structSize keeps callers compiled against an older layout safe.
///
/// Contracts common to all entry points:
/// - Invalid arguments raise R errors (Rf_error), which longjmp through the
///   caller's frames; call from contexts that are safe to unwind.
/// - Functions that draw (creation, run, the prior samplers, predictor
///   updates) manage R's RNG state internally and must be called from the
///   main R thread. Do not wrap them in a GetRNGstate/PutRNGstate bracket
///   that spans your own draws through R's API.
/// - Creation preserves the data specification object against garbage
///   collection for the sampler's lifetime, and the engine borrows its
///   predictors only to quantize them into owned codes at construction; the
///   preserved data object is thereafter the sampler's own predictor GC anchor
///   and the call-time raw source for saved-tree replay and state restore.
///   Raw-array setters borrow for the call's duration only: setResponse,
///   setOffset, and setWeights keep the array alive until replaced; the
///   predictor setters re-quantize the borrow and retain no pointer.
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

#define DBARTS_C_API_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

/// Opaque sampler handle.
typedef struct dbarts_sampler_t dbarts_sampler;

/// Caller-owned, growable output buffers for dbarts_sampler_run. The caller
/// MUST set structSize to sizeof(dbarts_results) as the caller compiled it;
/// the library fills only fields whose end offset falls within structSize,
/// so a caller built against an older (smaller) header is never written
/// past. Fields append monotonically below the marked boundary and never
/// reorder across releases; an append bumps DBARTS_C_API_VERSION. A field is
/// filled only when both present-by-size and non-null: a null member skips
/// that quantity, and a zero structSize skips everything. k requires a k
/// hyperprior (dbarts_sampler_kIsSampled), varprobs a DART tree prior
/// (dbarts_sampler_usesDart), and tau/groupEffects a grouped
/// random-intercept sampler; each is left untouched otherwise. logLikelihood
/// carries the per-draw training-data log-likelihood for the gaussian,
/// binary, and aft families; aft reports the log density for events and the
/// log survival tail for right-censored observations, and a grouped
/// (random-intercept) sampler composes it over the per-group fits. It is
/// NaN-filled where the per-observation location is not fully recorded (a BCF
/// two-forest fit), and skipping it (null or absent-by-size) elides all of
/// its computation.
/// Value-initialize and set the size:
///   dbarts_results results = {0}; results.structSize = sizeof results;
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
  /* 1.0-0 field boundary: every future append goes below this line, never
     above, and bumps DBARTS_C_API_VERSION. */
} dbarts_results;

/// True when the caller's struct (per structSize) actually carries `field`.
/// The sizeof operand is unevaluated, so this never dereferences past the
/// caller's buffer.
#define DBARTS_RESULTS_HAS(r, field) \
  ((r)->structSize >= offsetof(dbarts_results, field) + sizeof((r)->field))

/// Returns DBARTS_C_API_VERSION of the installed package.
int dbarts_apiVersion(void);

/// Creates a sampler from the R specification objects (dbartsControl,
/// dbartsModel, dbartsData). family selects the response model for binary
/// responses: "" or "probit" give probit latents, "logistic" the
/// Polya-Gamma sampler; continuous responses are gaussian and accept "" or
/// "gaussian", or "aft" for an accelerated failure time (log-normal) survival
/// fit, in which case the response holds log survival/censoring times and a
/// numeric status vector (1 = event, 0 = right-censored) is read from the
/// control's "bartcore.survival" attribute. A verbose control prints the
/// initial summary here.
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
/// Registers callback for subsequent runs, or clears it when null; userData is
/// passed back unchanged and its lifetime is the caller's. Raises an error
/// while chains would run on worker threads (numThreads > 1 and numChains > 1),
/// which must never call into R. Inline multi-chain runs execute chains
/// sequentially - chain 0 completes all its sweeps before chain 1 starts - so
/// the callback conditions each chain to completion in turn and cannot see one
/// chain's progress while advancing another.
void dbarts_sampler_setCallback(dbarts_sampler* sampler,
                                dbarts_sampler_callback callback,
                                void* userData);

/// y has numObservations values; for binary families they must be 0/1.
void dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y);
/// offset has numObservations values or is null to remove. updateScale
/// rescales the internal response transform to the offset-adjusted range
/// (gaussian only); pass false once burnt in so fits stay comparable.
void dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                              int updateScale);
/// weights has numObservations values; gaussian responses only.
void dbarts_sampler_setWeights(dbarts_sampler* sampler,
                               const double* weights);
/// Holds the residual standard deviation at sigma (original response scale)
/// until the next call or gaussian draw; the Gibbs conditioning hook.
void dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma);
/// Copies the latent response values (numObservations x numChains) into
/// out. Returns 1, or 0 without touching out when the family has no
/// latents (gaussian).
int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out);

/// Replaces the full numObservations x numPredictors matrix. Transactional:
/// returns 1 when every tree still has non-empty leaves under the new
/// values (or forceUpdate), 0 after rolling back. Errors when cut points
/// would invalidate existing splits and updateCutPoints is false.
int dbarts_sampler_setPredictor(dbarts_sampler* sampler, const double* x,
                                int forceUpdate, int updateCutPoints);
/// Replaces numColumns columns, 0-indexed by columns; x holds them
/// contiguously. Same transaction contract as setPredictor.
int dbarts_sampler_updatePredictor(dbarts_sampler* sampler, const double* x,
                                   const size_t* columns, size_t numColumns,
                                   int forceUpdate, int updateCutPoints);

/// x_test is numTestObservations x numPredictors, or null with 0 rows to
/// remove test data (clearing any test offset).
void dbarts_sampler_setTestPredictors(dbarts_sampler* sampler,
                                      const double* x_test,
                                      size_t numTestObservations);
/// offset_test has numTestObservations values or is null to remove.
void dbarts_sampler_setTestOffset(dbarts_sampler* sampler,
                                  const double* offset_test);

/// Fits for new data on the original response scale (binary families give
/// the latent scale). With tree storage out is numTestObservations x
/// numSavedSamples x numChains from the saved trees; without, one set per
/// chain from the live trees. offset_test, when non-null, is added to
/// every sample's fits.
void dbarts_sampler_predict(dbarts_sampler* sampler, const double* x_test,
                            size_t numTestObservations,
                            const double* offset_test, double* out);

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
/// The n column replays the creation specification's predictors through each
/// saved tree (the engine keeps no predictor matrix); it is left unpopulated
/// for a sparse creation spec, and live trees carry their own counts.
SEXP dbarts_sampler_getTrees(dbarts_sampler* sampler,
                             const size_t* chainIndices,
                             size_t numChainIndices,
                             const size_t* sampleIndices,
                             size_t numSampleIndices,
                             const size_t* treeIndices,
                             size_t numTreeIndices, int useLiveTrees);
void dbarts_sampler_printTrees(dbarts_sampler* sampler,
                               const size_t* chainIndices,
                               size_t numChainIndices,
                               const size_t* sampleIndices,
                               size_t numSampleIndices,
                               const size_t* treeIndices,
                               size_t numTreeIndices);

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
size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler);
/// The saved-tree buffer size, 0 without tree storage; the sample count
/// predict produces.
size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler);
int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler);
int dbarts_sampler_usesDart(const dbarts_sampler* sampler);

#ifdef __cplusplus
}
#endif

#endif // DBARTS_DBARTS_H
