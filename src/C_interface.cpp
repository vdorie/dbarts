// Implementation of the flat C API (inst/include/dbarts/dbarts.h) over the
// bartcore engine; entry points are registered with R_RegisterCCallable in
// R_interface.cpp. Argument validation is minimal - consumers are compiled
// packages - except where an engine invariant is at stake (categorical
// category codes, column ranges).

#include <dbarts/dbarts.h>

#include <cstddef> // size_t
#include <cstdint> // uint64_t
#include <cstring> // memcpy
#include <type_traits> // is_same

#include <external/Rinternals.h>

#include <misc/linearAlgebra.h> // misc_addVectorsInPlace

#include "R_interface_bartcore_common.hpp"

using std::size_t;
using bartcore_bridge::BartcoreHolder;
using bartcore_bridge::refuseMultiForestMutation;
using bartcore_bridge::refuseMultiForestTransactionalUpdate;
using bartcore_bridge::refusePinnedSigmaChange;
using bartcore_bridge::validateColumnValues;

struct dbarts_sampler_t {
  BartcoreHolder* holder;
  SEXP data; // preserved against collection for the columns it lends
  dbarts_sampler_callback callback = NULL; // per-sweep conditioning hook
  void* callbackData = NULL;
};

namespace {

inline bartcore::SamplerBase& samplerOf(dbarts_sampler* sampler) {
  return *sampler->holder->sampler;
}
inline const bartcore::SamplerBase& samplerOf(const dbarts_sampler* sampler) {
  return *sampler->holder->sampler;
}

// The dense predictor matrix of a retained dbartsData spec (@x), or null when
// it is absent or non-dense (a sparse creation spec): the call-time raw source
// for saved-tree replay and cross-grid state restore. The spec is preserved
// against collection for the sampler's lifetime, so the borrow is valid.
const double* predictorsFromDataExpr(SEXP dataExpr) {
  if (dataExpr == NULL || Rf_isNull(dataExpr)) return NULL;
  SEXP xExpr = Rf_getAttrib(dataExpr, Rf_install("x"));
  if (!Rf_isReal(xExpr)) return NULL;
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2) return NULL;
  return REAL(xExpr);
}

} // namespace

// Layout lock for dbarts_results (dbarts.h): the growable ABI. Fields append
// monotonically and never reorder, so the library's offsetof matches every
// caller's. Every field's exact offset is pinned (all trailing members are
// pointer-sized), so a mid-struct insertion shifts a downstream offset and
// fails here, a reorder fails at the swapped pair, and the size assert forces
// an author who appends a field to update it (and bump DBARTS_C_API_MINOR).
static_assert(offsetof(dbarts_results, structSize) == 0);
static_assert(offsetof(dbarts_results, sigma) == sizeof(size_t) + 0 * sizeof(double*));
static_assert(offsetof(dbarts_results, train) == sizeof(size_t) + 1 * sizeof(double*));
static_assert(offsetof(dbarts_results, test) == sizeof(size_t) + 2 * sizeof(double*));
static_assert(offsetof(dbarts_results, varcount) == sizeof(size_t) + 3 * sizeof(double*));
static_assert(offsetof(dbarts_results, k) == sizeof(size_t) + 4 * sizeof(double*));
static_assert(offsetof(dbarts_results, varprobs) == sizeof(size_t) + 5 * sizeof(double*));
static_assert(offsetof(dbarts_results, tau) == sizeof(size_t) + 6 * sizeof(double*));
static_assert(offsetof(dbarts_results, groupEffects) == sizeof(size_t) + 7 * sizeof(double*));
static_assert(offsetof(dbarts_results, logLikelihood) == sizeof(size_t) + 8 * sizeof(double*));
static_assert(sizeof(dbarts_results) == sizeof(size_t) + 9 * sizeof(double*),
              "dbarts_results layout changed; update these offsets, and bump "
              "DBARTS_C_API_MINOR if a field was appended after 1.0-0");

// Compile-time signature token: FNV-1a over the stringized
// DBARTS_C_API_LIST, checked against the baked DBARTS_C_API_HASH. A changed
// signature moves the hash and fails this assert until DBARTS_C_API_HASH is
// re-baked - the mechanical acknowledgment that the ABI surface changed.
namespace {
constexpr std::uint64_t dbarts_fnv1a(const char* text) {
  std::uint64_t hash = 0xcbf29ce484222325ULL; // FNV-1a 64-bit offset basis
  while (*text != '\0') {
    hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(*text));
    hash *= 0x100000001b3ULL; // FNV-1a 64-bit prime
    ++text;
  }
  return hash;
}
} // namespace
static_assert(dbarts_fnv1a(DBARTS_C_API_DECLS) == DBARTS_C_API_HASH,
              "dbarts.h C API signatures changed; re-bake DBARTS_C_API_HASH in "
              "inst/include/dbarts/dbarts.h (and bump DBARTS_C_API_MAJOR or "
              "DBARTS_C_API_MINOR as the change warrants)");

extern "C" {

int dbarts_apiVersion(void) { return DBARTS_C_API_VERSION; }
int dbarts_apiMajorVersion(void) { return DBARTS_C_API_MAJOR; }
int dbarts_apiMinorVersion(void) { return DBARTS_C_API_MINOR; }

dbarts_sampler* dbarts_sampler_create(SEXP control, SEXP model, SEXP data,
                                      const char* family) {
  BartcoreHolder* holder = bartcore_bridge::createHolder(
    control, model, data, family == NULL ? "" : family);
  R_PreserveObject(data);
  return new dbarts_sampler_t{holder, data};
}

void dbarts_sampler_destroy(dbarts_sampler* sampler) {
  if (sampler == NULL) return;
  delete sampler->holder;
  R_ReleaseObject(sampler->data);
  delete sampler;
}

void dbarts_sampler_run(dbarts_sampler* sampler, size_t numBurnIn,
                        size_t numSamples, dbarts_results* results) {
  bartcore::SamplerShape shape = samplerOf(sampler).shape();
  bartcore::Results engineResults;
  // the internal location stride (invisible to the frozen dbarts_results ABI):
  // 1 for every dbarts.h-created sampler, since the flat C API builds no
  // multi-location model, so the caller's n x numSamples train/test hold
  engineResults.numReportedLocations = shape.numReportedLocations;
  // A zero structSize means the caller forgot to set it (see DBARTS_RESULTS_INIT):
  // reject loudly instead of silently skipping every field and handing back an
  // uninitialized buffer - the flat-API footgun that fed garbage draws to a
  // consumer's Gibbs loop. A nonzero older/smaller structSize stays valid.
  if (results != NULL && results->structSize == 0)
    Rf_error("dbarts_sampler_run: results.structSize is 0 - set it to "
             "sizeof(dbarts_results) (e.g. dbarts_results r = DBARTS_RESULTS_INIT)");

  if (results != NULL && numSamples > 0) {
    // A field is filled only when present-by-size AND non-null. offsetof is
    // against the library's (newest) layout; fields only append, so it
    // equals the caller's offset and structSize bounds the buffer.
#define FILL(field, member) \
  engineResults.member = DBARTS_RESULTS_HAS(results, field) ? results->field : NULL
    FILL(sigma, sigma);
    FILL(train, trainingFits);
    FILL(test, testFits);
    FILL(varcount, variableCounts);
    FILL(k, k);
    FILL(varprobs, splitProbabilities);
    FILL(tau, tau);
    FILL(groupEffects, groupEffects);
    FILL(logLikelihood, logLikelihood);
#undef FILL
  }

  bartcore::SweepCallback onSweep;
  if (sampler->callback != NULL) {
    if (shape.numThreads > 1 && shape.numChains > 1)
      Rf_error("dbarts_sampler_run: a per-sweep callback cannot run while "
               "chains execute on worker threads");
    dbarts_sampler_callback fn = sampler->callback;
    void* userData = sampler->callbackData;
    onSweep = [sampler, fn, userData](size_t chainIndex, size_t sweepIndex,
                                      bool isBurnIn) -> bool {
      return fn(userData, sampler, chainIndex, sweepIndex,
                isBurnIn ? 1 : 0) == 0;
    };
  }

  // The engine samples only from each chain's own Mersenne Twister (seeded
  // from R's stream once at creation), never from R's stream during a run, so
  // no GetRNGstate/PutRNGstate bracket is needed here - and none is left
  // unbalanced by a longjmp out of the engine.
  samplerOf(sampler).run(numBurnIn, numSamples, engineResults, {}, onSweep);
}

void dbarts_sampler_setCallback(dbarts_sampler* sampler,
                                dbarts_sampler_callback callback,
                                void* userData) {
  if (callback != NULL) {
    bartcore::SamplerShape shape = samplerOf(sampler).shape();
    if (shape.numThreads > 1 && shape.numChains > 1)
      Rf_error("dbarts_sampler_setCallback: a per-sweep callback requires "
               "chains to run inline (numThreads == 1 or numChains == 1)");
  }
  sampler->callback = callback;
  sampler->callbackData = userData;
}

void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler) {
  // draws from the chain RNG only, not R's stream (see dbarts_sampler_run)
  samplerOf(sampler).sampleTreesFromPrior();
}

void dbarts_sampler_sampleNodeParametersFromPrior(dbarts_sampler* sampler) {
  // draws from the chain RNG only, not R's stream (see dbarts_sampler_run)
  samplerOf(sampler).sampleNodeParametersFromPrior();
}

void dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y) {
  // unreachable today (ResponseFamily has no multi-forest member on the flat
  // C API's creation path) but guarded defensively, matching the bridge entry
  refuseMultiForestMutation(samplerOf(sampler), "dbarts_sampler_setResponse");
  // the probit latent redraw draws from the chain RNG, not R's stream
  samplerOf(sampler).setResponse(y, true); // flat ABI keeps re-anchoring
}

void dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                              int updateScale) {
  // the offset is the response-side swap under a different pointer; unreachable
  // today, guarded defensively, see dbarts_sampler_setResponse
  refuseMultiForestMutation(samplerOf(sampler), "dbarts_sampler_setOffset");
  samplerOf(sampler).setOffset(offset, updateScale != 0);
}

void dbarts_sampler_setWeights(dbarts_sampler* sampler,
                               const double* weights) {
  // unreachable today, guarded defensively; see dbarts_sampler_setResponse
  refuseMultiForestMutation(samplerOf(sampler), "dbarts_sampler_setWeights");
  samplerOf(sampler).setWeights(weights);
}

void dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma) {
  // reachable here: this entry creates probit/logistic/ordinal/nbinom samplers
  // by family name, and dbartsSpec(variance = ) hands a consumer a
  // heteroscedastic control (multinomial has no flat creation path)
  refusePinnedSigmaChange(samplerOf(sampler), "dbarts_sampler_setSigma");
  samplerOf(sampler).setSigma(sigma);
}

int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  if (engine.latents(0) == NULL) return 0;

  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  for (size_t c = 0; c < shape.numChains; ++c)
    std::memcpy(out + c * numObservations, engine.latents(c),
                numObservations * sizeof(double));
  return 1;
}

int dbarts_sampler_setPredictor(dbarts_sampler* sampler, const double* x,
                                int forceUpdate, int updateCutPoints) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // the variance clause is reachable today: dbartsSpec(variance = ) hands a
  // consumer a flat-creatable heteroscedastic sampler, and an unforced
  // transactional call here used to accept and silently misroute s^2(x); the
  // numForests >= 2 clause is unreachable until BCF gains a flat creation
  // path, and becomes load-bearing from that tip on
  refuseMultiForestTransactionalUpdate(engine, "dbarts_sampler_setPredictor",
                                       forceUpdate != 0);
  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  for (size_t j = 0; j < shape.numPredictors; ++j)
    validateColumnValues(engine.data(), j, x + j * numObservations,
                         numObservations);

  bartcore::PredictorUpdateResult result =
    engine.setPredictor(x, forceUpdate != 0, updateCutPoints != 0);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  return result == bartcore::PredictorUpdateResult::accepted ? 1 : 0;
}

int dbarts_sampler_updatePredictor(dbarts_sampler* sampler, const double* x,
                                   const size_t* columns, size_t numColumns,
                                   int forceUpdate, int updateCutPoints) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // same guard as dbarts_sampler_setPredictor above
  refuseMultiForestTransactionalUpdate(
    engine, "dbarts_sampler_updatePredictor", forceUpdate != 0);
  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  for (size_t k = 0; k < numColumns; ++k) {
    if (columns[k] >= shape.numPredictors)
      Rf_error("dbarts_sampler_updatePredictor column out of range");
    validateColumnValues(engine.data(), columns[k], x + k * numObservations,
                         numObservations);
  }

  bartcore::PredictorUpdateResult result = engine.updatePredictor(
    x, columns, numColumns, forceUpdate != 0, updateCutPoints != 0);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  return result == bartcore::PredictorUpdateResult::accepted ? 1 : 0;
}

void dbarts_sampler_setTestPredictors(dbarts_sampler* sampler,
                                      const double* xTest,
                                      size_t numTestObservations) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  if (xTest == NULL) {
    engine.setTestPredictors(NULL, 0);
    return;
  }
  for (size_t j = 0; j < engine.shape().numPredictors; ++j)
    validateColumnValues(engine.data(), j, xTest + j * numTestObservations,
                         numTestObservations);
  engine.setTestPredictors(xTest, numTestObservations);
}

void dbarts_sampler_setTestOffset(dbarts_sampler* sampler,
                                  const double* offsetTest) {
  // a multi-forest test offset lands after the forests are blended; unreachable
  // today, guarded defensively, see dbarts_sampler_setResponse
  refuseMultiForestMutation(samplerOf(sampler), "dbarts_sampler_setTestOffset");
  samplerOf(sampler).setTestOffset(offsetTest);
}

void dbarts_sampler_predict(dbarts_sampler* sampler, const double* xTest,
                            size_t numTestObservations,
                            const double* offsetTest, double* out) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  for (size_t j = 0; j < shape.numPredictors; ++j)
    validateColumnValues(engine.data(), j, xTest + j * numTestObservations,
                         numTestObservations);

  engine.predict(xTest, numTestObservations, out);

  if (offsetTest != NULL) {
    size_t capacity = shape.savedTreeCapacity;
    size_t numSamples = capacity > 0 ? capacity : 1;
    for (size_t slab = 0; slab < numSamples * shape.numChains; ++slab)
      misc_addVectorsInPlace(offsetTest, numTestObservations,
                             out + slab * numTestObservations);
  }
}

void dbarts_sampler_setTreeStorage(dbarts_sampler* sampler, int keepTrees,
                                   size_t numSamplesToStore) {
  samplerOf(sampler).setTreeStorage(keepTrees != 0, numSamplesToStore);
}

SEXP dbarts_sampler_getTrees(dbarts_sampler* sampler,
                             const size_t* chainIndices,
                             size_t numChainIndices,
                             const size_t* sampleIndices,
                             size_t numSampleIndices,
                             const size_t* treeIndices, size_t numTreeIndices,
                             int useLiveTrees) {
  // the n column replays the retained creation spec's predictors through each
  // saved tree; the engine keeps no matrix, and a caller that mutated
  // predictors since creation sees the pre-mutation spec
  const double* replay = predictorsFromDataExpr(sampler->data);
  return bartcore_bridge::getTrees(
    samplerOf(sampler), chainIndices, numChainIndices, sampleIndices,
    numSampleIndices, treeIndices, numTreeIndices, useLiveTrees != 0, NULL, 0,
    replay, samplerOf(sampler).shape().numObservations, 0,
    "dbarts_sampler_getTrees");
}

void dbarts_sampler_printTrees(dbarts_sampler* sampler,
                               const size_t* chainIndices,
                               size_t numChainIndices,
                               const size_t* sampleIndices,
                               size_t numSampleIndices,
                               const size_t* treeIndices,
                               size_t numTreeIndices) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  for (size_t i = 0; i < numChainIndices; ++i) {
    if (chainIndices[i] >= shape.numChains)
      Rf_error("dbarts_sampler_printTrees chain number out of range");
  }
  for (size_t i = 0; i < numSampleIndices; ++i) {
    if (sampleIndices[i] >= shape.savedTreeCapacity)
      Rf_error("dbarts_sampler_printTrees sample number out of range");
  }
  for (size_t i = 0; i < numTreeIndices; ++i) {
    if (treeIndices[i] >= shape.numTrees)
      Rf_error("dbarts_sampler_printTrees tree number out of range");
  }
  engine.printTrees(chainIndices, numChainIndices, sampleIndices,
                    numSampleIndices, treeIndices, numTreeIndices);
}

SEXP dbarts_sampler_storeState(dbarts_sampler* sampler) {
  return bartcore_bridge::storeState(samplerOf(sampler));
}

void dbarts_sampler_setState(dbarts_sampler* sampler, SEXP state) {
  // a cross-grid restore re-quantizes dense columns from the retained creation
  // spec; a same-spec continuation (the contract) skips per column and reads no
  // raw. A flat-C caller that mutated predictors since sees the creation spec.
  bartcore_bridge::setState(samplerOf(sampler), state,
                            predictorsFromDataExpr(sampler->data));
}

void dbarts_sampler_setNumThreads(dbarts_sampler* sampler,
                                  size_t numThreads) {
  samplerOf(sampler).setNumThreads(numThreads);
}

void dbarts_sampler_setNumThin(dbarts_sampler* sampler, size_t numThin) {
  samplerOf(sampler).setNumThin(numThin);
}

void dbarts_sampler_setVerbose(dbarts_sampler* sampler, int verbose,
                               uint32_t printEvery) {
  samplerOf(sampler).setVerbose(verbose != 0, printEvery);
}

size_t dbarts_sampler_numObservations(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numObservations;
}

size_t dbarts_sampler_numPredictors(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numPredictors;
}

size_t dbarts_sampler_numTestObservations(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numTestObservations;
}

size_t dbarts_sampler_numChains(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numChains;
}

size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numTrees;
}

size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().savedTreeCapacity;
}

int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().kIsSampled ? 1 : 0;
}

int dbarts_sampler_usesDart(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().usesDart ? 1 : 0;
}

// Provider-side binding: each real function's address must
// have exactly the type the single-source DBARTS_C_API_LIST ascribes to it, so
// any drift between the readable prototypes above and the list fails dbarts's
// own compile. Placed inside extern "C" so the list-formed pointer type carries
// the same C language linkage as the resolved &function.
#define DBARTS_BIND_ASSERT(ret, name, params, args) \
  static_assert(std::is_same<decltype(&name), ret (*) params>::value, \
                #name " signature drifted from DBARTS_C_API_LIST");
DBARTS_C_API_LIST(DBARTS_BIND_ASSERT)
#undef DBARTS_BIND_ASSERT

} // extern "C"
