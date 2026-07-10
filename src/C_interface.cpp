// Implementation of the flat C API (inst/include/dbarts/dbarts.h) over the
// bartcore engine; entry points are registered with R_RegisterCCallable in
// R_interface.cpp. Argument validation is minimal - consumers are compiled
// packages - except where an engine invariant is at stake (categorical
// category codes, column ranges).

#include <dbarts/dbarts.h>

#include <cstddef> // size_t
#include <cstring> // memcpy

#include <external/Rinternals.h>
#include <R_ext/Random.h> // GetRNGstate, PutRNGstate

#include <misc/linearAlgebra.h> // misc_addVectorsInPlace

#include "R_interface_bartcore_common.hpp"

using std::size_t;
using bartcore_bridge::BartcoreHolder;
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

} // namespace

// Layout lock for dbarts_results (dbarts.h): the growable ABI. Fields append
// monotonically and never reorder, so the library's offsetof matches every
// caller's; the exact-size assert forces an author who appends a field to
// update it (and, by convention, bump DBARTS_C_API_VERSION).
static_assert(offsetof(dbarts_results, structSize) == 0);
static_assert(offsetof(dbarts_results, sigma) == sizeof(size_t));
static_assert(offsetof(dbarts_results, sigma) < offsetof(dbarts_results, train));
static_assert(offsetof(dbarts_results, train) < offsetof(dbarts_results, test));
static_assert(offsetof(dbarts_results, test) < offsetof(dbarts_results, varcount));
static_assert(offsetof(dbarts_results, varcount) < offsetof(dbarts_results, k));
static_assert(offsetof(dbarts_results, k) < offsetof(dbarts_results, varprobs));
static_assert(offsetof(dbarts_results, varprobs) < offsetof(dbarts_results, tau));
static_assert(offsetof(dbarts_results, tau) <
              offsetof(dbarts_results, groupEffects));
static_assert(sizeof(dbarts_results) == sizeof(size_t) + 8 * sizeof(double*),
              "dbarts_results grew; update this size and DBARTS_C_API_VERSION");

extern "C" {

int dbarts_apiVersion(void) { return DBARTS_C_API_VERSION; }

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
  bartcore::Results engineResults;
  if (results != NULL && numSamples > 0) {
    // A field is filled only when present-by-size AND non-null; a zero
    // structSize (caller forgot to set it) is an all-skip no-op. offsetof is
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
#undef FILL
  }

  bartcore::SweepCallback onSweep;
  if (sampler->callback != NULL) {
    bartcore::SamplerBase& engine(samplerOf(sampler));
    if (engine.numThreads() > 1 && engine.numChains() > 1)
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

  GetRNGstate();
  samplerOf(sampler).run(numBurnIn, numSamples, engineResults, {}, onSweep);
  PutRNGstate();
}

void dbarts_sampler_setCallback(dbarts_sampler* sampler,
                                dbarts_sampler_callback callback,
                                void* userData) {
  if (callback != NULL) {
    const bartcore::SamplerBase& engine(samplerOf(sampler));
    if (engine.numThreads() > 1 && engine.numChains() > 1)
      Rf_error("dbarts_sampler_setCallback: a per-sweep callback requires "
               "chains to run inline (numThreads == 1 or numChains == 1)");
  }
  sampler->callback = callback;
  sampler->callbackData = userData;
}

void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler) {
  GetRNGstate();
  samplerOf(sampler).sampleTreesFromPrior();
  PutRNGstate();
}

void dbarts_sampler_sampleNodeParametersFromPrior(dbarts_sampler* sampler) {
  GetRNGstate();
  samplerOf(sampler).sampleNodeParametersFromPrior();
  PutRNGstate();
}

void dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y) {
  GetRNGstate(); // probit latent redraw
  samplerOf(sampler).setResponse(y, true); // flat ABI keeps re-anchoring
  PutRNGstate();
}

void dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                              int updateScale) {
  samplerOf(sampler).setOffset(offset, updateScale != 0);
}

void dbarts_sampler_setWeights(dbarts_sampler* sampler,
                               const double* weights) {
  samplerOf(sampler).setWeights(weights);
}

void dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma) {
  samplerOf(sampler).setSigma(sigma);
}

int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  if (engine.latents(0) == NULL) return 0;

  size_t numObservations = engine.numObservations();
  for (size_t c = 0; c < engine.numChains(); ++c)
    std::memcpy(out + c * numObservations, engine.latents(c),
                numObservations * sizeof(double));
  return 1;
}

int dbarts_sampler_setPredictor(dbarts_sampler* sampler, const double* x,
                                int forceUpdate, int updateCutPoints) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  size_t numObservations = engine.numObservations();
  for (size_t j = 0; j < engine.numPredictors(); ++j)
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
  size_t numObservations = engine.numObservations();
  for (size_t k = 0; k < numColumns; ++k) {
    if (columns[k] >= engine.numPredictors())
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
                                      const double* x_test,
                                      size_t numTestObservations) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  if (x_test == NULL) {
    engine.setTestPredictors(NULL, 0);
    return;
  }
  for (size_t j = 0; j < engine.numPredictors(); ++j)
    validateColumnValues(engine.data(), j, x_test + j * numTestObservations,
                         numTestObservations);
  engine.setTestPredictors(x_test, numTestObservations);
}

void dbarts_sampler_setTestOffset(dbarts_sampler* sampler,
                                  const double* offset_test) {
  samplerOf(sampler).setTestOffset(offset_test);
}

void dbarts_sampler_predict(dbarts_sampler* sampler, const double* x_test,
                            size_t numTestObservations,
                            const double* offset_test, double* out) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  for (size_t j = 0; j < engine.numPredictors(); ++j)
    validateColumnValues(engine.data(), j, x_test + j * numTestObservations,
                         numTestObservations);

  engine.predict(x_test, numTestObservations, out);

  if (offset_test != NULL) {
    size_t capacity = engine.savedTreeCapacity();
    size_t numSamples = capacity > 0 ? capacity : 1;
    for (size_t slab = 0; slab < numSamples * engine.numChains(); ++slab)
      misc_addVectorsInPlace(offset_test, numTestObservations,
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
  return bartcore_bridge::getTrees(
    samplerOf(sampler), chainIndices, numChainIndices, sampleIndices,
    numSampleIndices, treeIndices, numTreeIndices, useLiveTrees != 0, NULL, 0,
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
  for (size_t i = 0; i < numChainIndices; ++i) {
    if (chainIndices[i] >= engine.numChains())
      Rf_error("dbarts_sampler_printTrees chain number out of range");
  }
  for (size_t i = 0; i < numSampleIndices; ++i) {
    if (sampleIndices[i] >= engine.savedTreeCapacity())
      Rf_error("dbarts_sampler_printTrees sample number out of range");
  }
  for (size_t i = 0; i < numTreeIndices; ++i) {
    if (treeIndices[i] >= engine.numTrees())
      Rf_error("dbarts_sampler_printTrees tree number out of range");
  }
  engine.printTrees(chainIndices, numChainIndices, sampleIndices,
                    numSampleIndices, treeIndices, numTreeIndices);
}

SEXP dbarts_sampler_storeState(dbarts_sampler* sampler) {
  return bartcore_bridge::storeState(samplerOf(sampler));
}

void dbarts_sampler_setState(dbarts_sampler* sampler, SEXP state) {
  bartcore_bridge::setState(samplerOf(sampler), state);
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
  return samplerOf(sampler).numObservations();
}

size_t dbarts_sampler_numPredictors(const dbarts_sampler* sampler) {
  return samplerOf(sampler).numPredictors();
}

size_t dbarts_sampler_numTestObservations(const dbarts_sampler* sampler) {
  return samplerOf(sampler).numTestObservations();
}

size_t dbarts_sampler_numChains(const dbarts_sampler* sampler) {
  return samplerOf(sampler).numChains();
}

size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler) {
  return samplerOf(sampler).numTrees();
}

size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler) {
  return samplerOf(sampler).savedTreeCapacity();
}

int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler) {
  return samplerOf(sampler).kIsSampled() ? 1 : 0;
}

int dbarts_sampler_usesDart(const dbarts_sampler* sampler) {
  return samplerOf(sampler).usesDart() ? 1 : 0;
}

} // extern "C"
