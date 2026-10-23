#ifndef BARTCORE_SAMPLER_HPP
#define BARTCORE_SAMPLER_HPP

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#ifndef _WIN32
#include <csignal>  // block SIGINT in worker threads; see run()
#endif
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <external/random.h>

#include "chain.hpp"
#include "data.hpp"
#include "model.hpp"

namespace bartcore {

/// The test-only channel for an out-of-sample replay's fan-out. Everything it
/// reports is written on the calling thread before any worker starts, so a
/// test can prove that a thread argument reached the engine and that the
/// partition covers every slab exactly once - deterministically, without
/// measuring time. cutoffOverride is the one field a test WRITES: it lets a
/// fixture small enough to run in a unit test still fan out, which is what
/// makes an identity-across-thread-counts assertion non-vacuous there.
///
/// One object serves every sampler because the replay entry points are
/// contractually main-thread-only (they are R_alloc-backed at the boundary),
/// so nothing writes it concurrently. None of it can change an answer.
struct PredictPartitionChannel {
  /// The count resolved from the argument and the sampler, floored at 1 and
  /// BEFORE the slab and cutoff clamps below.
  std::size_t resolvedThreads = 0;
  /// The workers actually run, min(resolvedThreads, slabs), or 1 below cutoff.
  std::size_t numWorkers = 0;
  /// workerForSlab[s] is the worker index that replayed slab s.
  std::vector<std::size_t> workerForSlab;
  /// When non-zero, the traversal cutoff in force instead of the derived one.
  std::size_t cutoffOverride = 0;
};
inline PredictPartitionChannel predictPartition;

/// Prints progress lines as they arrive; only safe on the main thread.
struct DirectProgressSink final : ProgressSink {
  void report(const char* line) override { ext_printf("%s", line); }
};

/// Buffers progress lines from worker threads for the main thread to flush;
/// workers must never call into R.
struct QueuedProgressSink final : ProgressSink {
  void report(const char* line) override {
    std::lock_guard<std::mutex> lock(mutex_);
    lines_.push_back(line);
  }
  void flush() {
    std::vector<std::string> pending;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      pending.swap(lines_);
    }
    for (const std::string& line : pending) ext_printf("%s", line.c_str());
  }

private:
  std::mutex mutex_;
  std::vector<std::string> lines_;
};

/// Outcome of a transactional predictor change. invalidCutPoints reports a
/// quantile-mode cut refresh whose new column would induce fewer cuts than
/// existing splits require; unlike the pre-1.0 engine, which errored midway
/// through installation, nothing has been modified.
enum class PredictorUpdateResult { accepted, rolledBack, invalidCutPoints };

/// Outcome of a response or offset swap. scaleKept: the swap asked for the
/// response transform to be re-derived while a row was out of the likelihood
/// and the rows in held fewer than two distinct values of what it is derived
/// from, so the transform in force was kept and the swap made as with
/// updateScale false. Nothing is raised: the caller decides what to say.
enum class ResponseUpdateResult { accepted, scaleKept };

/// A whole sampler's serializable state: per-chain states, the cut grid they
/// were read on, the saved-tree write position, and the draw count that
/// position is read against.
struct SamplerStateData {
  std::vector<ChainStateData> chains;
  /// The store's cut points when the state was read, an empty vector per
  /// column that splits by subset: the grid every live threshold of the
  /// chains is a point of. A record of where the trees were stored, which an
  /// install compares with the store's grid or, asked to adopt it, makes the
  /// store's (Sampler::setState).
  std::vector<std::vector<double>> cutPoints;
  /// That grid's weights (ColumnStore::cutMass), empty per unweighted column;
  /// the whole vector empty where none is weighted, as a state stored without
  /// weights reads. Written to the store only by an install that adopts the
  /// grid.
  std::vector<std::vector<double>> cutMass;
  size_t currentSampleNum = 0;
  /// Recorded draws the store retains, capped at its capacity. Required: it
  /// cannot be inferred from the buffer, whose unwritten slots hold zero-leaf
  /// trees that read as legitimate draws.
  size_t recordedDraws = 0;
  /// A byte per predictor: whether the column could hold a missing value when
  /// the state was read (ColumnStore::hasMissing), whatever it then held.
  /// Empty where the state carries no record, which does not say that none
  /// could.
  std::vector<std::uint8_t> missingColumns;
};

/// What a warm start (installForests) did: ok for a donor installed, declined
/// for one an unforced call did not install because a tree or a leaf value
/// would have had to change, and otherwise why it was refused. A single donor
/// forest can seed several chains, so the donor's chain count need not match.
/// gridMismatch: another number of columns, a column that splits by subset on
/// one side only, a donor grid that repeats a point, or, on a column the
/// store holds as an ordered factor, a donor grid other than the store's.
/// rebuildFailed: a donor tree does not build, on the grid the donor records
/// or on the store's; found on scratch before anything is written, and
/// undone where an install fails all the same. leafKindMismatch: gp leaves
/// on one side only; leafRowsMismatch: gp leaves fit to another number of
/// rows; leafBlockMismatch: a tree's leaf block not of the size the donor's
/// own leaves give.
enum class WarmStartResult {
  ok, declined, shapeMismatch, gridMismatch, varianceMismatch,
  varianceSlotMismatch, rebuildFailed, leafKindMismatch, leafRowsMismatch,
  leafBlockMismatch
};

/// A warm-start donor's leaves as its own blocks give them: constant with
/// none, linear with numCovariates slopes a leaf, gp with numObservations
/// fits a tree.
struct DonorLeafShape {
  enum class Kind { constant, linear, gp } kind = Kind::constant;
  std::size_t numCovariates = 0, numObservations = 0;
};

/// What a warm start installed, filled by installForests for a donor that
/// went in: a chain's counts for its mean forests and for its variance
/// forest, and the three ways the leaves it holds are not the donor's -
/// function fits started at 0 for a donor recorded on other cut points, a
/// linear donor's slopes left behind, and linear coefficients kept as built
/// across a covariate of one value on each side. An install is altered
/// exactly when a count names a change or one of those is up. The rest say
/// what the donor is read as and change nothing of it: slopes of 0 under a
/// constant donor, and DART or a variance forest on one side only, the
/// sampler's own started as a new fit's or the donor's left behind.
struct WarmStartReport {
  std::vector<WarmStartTreeCounts> mean, variance;
  bool fitsCold = false, slopesLeftBehind = false, coefficientsKept = false;
  bool slopesStartedAtZero = false, dartStartedFresh = false,
       dartLeftBehind = false, varianceStartedFresh = false,
       varianceLeftBehind = false;
  bool changed() const {
    bool any = fitsCold || slopesLeftBehind || coefficientsKept;
    for (const WarmStartTreeCounts& counts : mean) any = any || counts.changed();
    for (const WarmStartTreeCounts& counts : variance)
      any = any || counts.changed();
    return any;
  }
};

/// A sequential per-observation predictor update: stage one observation's
/// leaf moves, test that no leaf empties, then commit or skip, with a single
/// fits rebuild at the end. Type-erased so several samplers sharing an
/// index-aligned column can be swept jointly, committing all-or-none.
/// Relays a host's per-draw observer and records that one returned nonzero, so
/// a true return from Sampler::run can be told apart from an interrupt. One
/// relay serves every chain: it holds no per-chain state and the flag is a
/// relaxed set-once, read only after every worker has joined.
struct DrawRelay {
  DrawHook hook;
  std::atomic<bool> stopped{false};
  static int call(void* context, const DrawInfo* draw) {
    DrawRelay& relay = *static_cast<DrawRelay*>(context);
    int result = relay.hook.fn(relay.hook.context, draw);
    if (result != 0) relay.stopped.store(true, std::memory_order_relaxed);
    return result;
  }
};

class PredictorUpdateSession {
public:
  virtual ~PredictorUpdateSession() = default;
  /// Stages observation i's leaf moves against the running occupancy counts;
  /// true unless some tree's leaf would be left empty or, the value being the
  /// column's first missing one, some monotone tree's leaf values would be
  /// left out of order. A first missing value is judged with the column's
  /// rules holding the directions its arrival draws; they stay on, with the
  /// column's flag, only if the observation is then committed, and are taken
  /// back, generators included, by the next call here or by finalize
  /// otherwise. A false stages nothing the caller has to undo.
  virtual bool observationWouldRemainValid(size_t i) = 0;
  /// Installs observation i's staged value and occupancy moves.
  virtual void commitObservation(size_t i) = 0;
  /// Re-routes all observations and rebuilds fits. The sequential guard
  /// admits no empty leaf and no tree out of order, so false is an internal
  /// invariant violation.
  virtual bool finalize() = 0;
};

/// Starts numWorkers threads, worker w running body(w), with SIGINT blocked
/// across the spawn so the workers inherit the block and a Ctrl-C reaches only
/// the calling thread, whose mask is restored before return. body is copied
/// into each thread, so it must capture only what outlives the caller's join.
/// A spawn that fails (std::system_error at the process thread limit) must not
/// unwind past a joinable std::thread, which is std::terminate: stopStarted
/// tells the workers already running to finish early, they are joined, and the
/// exception is rethrown on the caller's thread with workers left empty.
template <typename Body, typename StopStarted>
void spawnWorkers(std::vector<std::thread>& workers, std::size_t numWorkers,
                  const Body& body, StopStarted&& stopStarted) {
  workers.reserve(numWorkers);
#ifndef _WIN32
  sigset_t interruptSet, previousSet;
  sigemptyset(&interruptSet);
  sigaddset(&interruptSet, SIGINT);
  pthread_sigmask(SIG_BLOCK, &interruptSet, &previousSet);
#endif
  try {
    for (std::size_t w = 0; w < numWorkers; ++w) workers.emplace_back(body, w);
  } catch (...) {
#ifndef _WIN32
    pthread_sigmask(SIG_SETMASK, &previousSet, nullptr);
#endif
    stopStarted();
    for (std::thread& worker : workers) worker.join();
    workers.clear();
    throw;
  }
#ifndef _WIN32
  pthread_sigmask(SIG_SETMASK, &previousSet, nullptr);
#endif
}

/// Sampler::setState's adoptCapacity default: judge a state's saved trees
/// against the live store and leave its capacity alone.
inline constexpr std::size_t keepStoreCapacity = static_cast<std::size_t>(-1);

/// The sampler proper: a shared column store and one or more chains over it.
/// Chains run independently (optionally on worker threads) and hold all
/// per-chain state; the sampler owns the data and orchestrates transactions,
/// which are all-or-none across every tree of every chain.
template <IntegrableLeafModel L, typename ResidT = double>
class Sampler {
public:
  /// rngs supplies one generator per chain (options.numChains of them). With
  /// numThreads > 1 no rng may call into R (or otherwise share state).
  /// weights apply only to the gaussian family.
  Sampler(const double* x, const double* y, size_t numObservations,
          size_t numPredictors, const double* weights, const double* offset,
          ResponseFamily family, double sigmaEstimate, double sigmaDf,
          double sigmaRawScale, const SamplerOptions& options,
          ext_rng* const* rngs)
    : options_(options), family_(family) {
    data_.sparseDensityCutoff = options.sparseDensityThreshold;
    if (!data_.build(creationPredictorSource(options.predictors, x,
                                             numObservations, numPredictors),
                     options.maxNumCutsPerVariable, options.maxNumCuts,
                     options.useQuantiles, options.leafCovariateColumns,
                     options.numLeafCovariates)) {
      ingestionRefused_ = true;
      return;
    }
    options_.maxNumCutsPerVariable = nullptr;  // borrowed; consumed by build
    // borrowed; consumed by build, which retains what the store needs (a
    // mapped build's CSC slices, its own copy of the dense block)
    options_.predictors = {};
    adoptMissingSeen();

    initializeChains(y, weights, offset, sigmaEstimate, sigmaDf,
                     sigmaRawScale, rngs);
  }

  /// A sampler over a pre-built store (ColumnStore::buildFromParent): a
  /// row-subset view carrying the parent's cut grid, so folds bin
  /// identically to the full data. y/weights/offset are the subset's
  /// vectors, kept alive by the caller. Views hold no raw predictor values,
  /// so the raw-x mutation surface (setPredictor, updatePredictor and the
  /// per-observation sessions, setData, setCutPoints, setState) must not be
  /// called on one; the bridge refuses beforehand.
  Sampler(ColumnStore&& store, const double* y, const double* weights,
          const double* offset, ResponseFamily family, double sigmaEstimate,
          double sigmaDf, double sigmaRawScale, const SamplerOptions& options,
          ext_rng* const* rngs)
    : options_(options), family_(family) {
    data_ = std::move(store);
    options_.maxNumCutsPerVariable = nullptr;
    // the view carries the parent's grid, so its types and counts are fixed
    options_.predictors = {};
    options_.useQuantiles = data_.useQuantiles;
    options_.missingSeen = nullptr;

    initializeChains(y, weights, offset, sigmaEstimate, sigmaDf,
                     sigmaRawScale, rngs);
  }

  /// A K-forest combining sampler over dense predictors, bcf's
  /// prognostic-plus-treatment pair being its K = 2 instance. The family
  /// rides the spec (gaussian, probit or logistic); family_ takes it rather
  /// than a pin, since every family-keyed predicate a host reads - the pinned
  /// sigma and binary weight refusals, the response-support test - answers
  /// through it. The CSC/mixed and view ingestion paths are not offered here.
  Sampler(const double* x, const double* y, size_t numObservations,
          size_t numPredictors, const double* weights, const double* offset,
          double sigmaEstimate, double sigmaDf, double sigmaRawScale,
          const SamplerOptions& options, const AmplitudeSpec& spec,
          ext_rng* const* rngs)
    : options_(options), family_(spec.family) {
    data_.sparseDensityCutoff = options.sparseDensityThreshold;
    if (!data_.build(denseCreationPredictorSource(options.predictors, x,
                                                 numObservations,
                                                 numPredictors),
                     options.maxNumCutsPerVariable, options.maxNumCuts,
                     options.useQuantiles)) {
      ingestionRefused_ = true;
      return;
    }
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
    adoptMissingSeen();
    // single-forest queries (numTrees, savedTree, printTrees) address forest 0,
    // whose count the K-length spec carries wherever it came from
    options_.numTrees = expandForestSpecs(spec)[0].forest.numTrees;

    chains_.reserve(options_.numChains);
    for (size_t c = 0; c < options_.numChains; ++c)
      chains_.push_back(std::make_unique<Chain<L, ResidT>>(
        data_, y, weights, offset, sigmaEstimate, sigmaDf, sigmaRawScale,
        options_, spec, rngs[c]));
    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
    recordAnchor();
  }

  /// A K-forest multinomial (softmax) sampler over dense predictors:
  /// constant leaf only, no response/weights/
  /// offset (the category labels ride the spec; only single-trial is
  /// supported here). The
  /// CSC/mixed and view ingestion paths are not offered here.
  Sampler(const double* x, size_t numObservations, size_t numPredictors,
          const SamplerOptions& options, const MultinomialSpec& spec,
          ext_rng* const* rngs)
    : options_(options), family_(ResponseFamily::logistic) {
    data_.sparseDensityCutoff = options.sparseDensityThreshold;
    if (!data_.build(denseCreationPredictorSource(options.predictors, x,
                                                 numObservations,
                                                 numPredictors),
                     options.maxNumCutsPerVariable, options.maxNumCuts,
                     options.useQuantiles)) {
      ingestionRefused_ = true;
      return;
    }
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
    adoptMissingSeen();
    // single-forest queries (numTrees, savedTree, printTrees) address forest 0
    options_.numTrees = spec.forest.numTrees;

    chains_.reserve(options_.numChains);
    for (size_t c = 0; c < options_.numChains; ++c)
      chains_.push_back(
        std::make_unique<Chain<L, ResidT>>(data_, options_, spec, rngs[c]));

    // saved-tree storage under keepTrees, as initializeChains does for the
    // other constructors (this multinomial path builds its chains directly and
    // so must init it here); each chain allocates all K forests' slots. Skipped
    // off keepTrees, so a non-keepTrees multinomial layout is unchanged.
    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
    recordAnchor();
  }

  // chains reference the store member, so the sampler's address is pinned
  Sampler(const Sampler&) = delete;
  Sampler& operator=(const Sampler&) = delete;

  /// Whether the creation build refused a predictor cell the store cannot
  /// represent. True leaves the sampler unusable - no chain was built - and
  /// the factory answers its host with a null sampler rather than one over an
  /// unchecked store.
  bool ingestionRefused() const { return ingestionRefused_; }

  /// The GP fallback census of the last run, summed over chains.
  GPFallbackTally gpFallbackTally() const {
    GPFallbackTally total;
    for (const auto& chain : chains_) total.add(chain->gpFallbackTally());
    return total;
  }

  /// The slow-count tally of the last run, summed over chains; empty off a
  /// normalized leaf.
  SlowCountTally slowCountTally() const {
    SlowCountTally total;
    for (const auto& chain : chains_) total.add(chain->slowCountTally());
    return total;
  }

  /// The monotone directions (numPredictors, borrowed from the leaf for the
  /// sampler's lifetime) and prior the sampler was built with; null and 0
  /// off the monotone leaf.
  const std::int8_t* monotoneDirections() const {
    if constexpr (requires(const L& leaf) { leaf.directions; })
      return chains_[0]->leaf().directions.data();
    else
      return nullptr;
  }
  std::uint8_t monotonePrior() const {
    if constexpr (requires(const L& leaf) { leaf.directions; })
      return static_cast<std::uint8_t>(chains_[0]->leaf().prior);
    else
      return 0;
  }

  /// Replace the test predictors from a borrowed view, keeping any test offset:
  /// the caller guarantees the row count still matches it (the bridge refuses
  /// otherwise). Passing a new offset too goes through setTestOffset. The test
  /// store shares the training cut grid and owns its raw, so the view need not
  /// outlive the call. Refuses (returns false, test store untouched) a
  /// designated leaf covariate that would be CSC-backed, since leaf models
  /// gather dense raw test covariates that sparse storage does not serve, and
  /// a factor cell that is not a level code of the training table.
  bool setTestData(const PredictorSource& source) {
    size_t numCovariates = numLeafCovariates();
    const size_t* covariateColumns = leafCovariateColumns();
    for (size_t k = 0; k < numCovariates; ++k)
      if (source.sourceOf(covariateColumns[k]) < 0) return false;
    const double* testOffset = data_.testOffset;
    if (!data_.buildTest(source)) return false;
    data_.testOffset = testOffset;
    for (auto& chain : chains_) chain->resizeTestStorage();
    return true;
  }

  /// Dense convenience spelling: a plain column-major test matrix. Every
  /// column is dense, so the CSC-backed leaf-covariate refusal cannot fire.
  void setTestPredictors(const double* x_test, size_t numTestObservations) {
    setTestData(densePredictorSource(x_test, numTestObservations,
                                     data_.numPredictors));
  }

  /// Borrowed, length numTestObservations (the caller validates); null
  /// clears. Added to recorded test fits; predictions take their offset as
  /// an argument instead.
  void setTestOffset(const double* testOffset) {
    data_.testOffset = testOffset;
  }

  /// Runs every chain numBurnIn + numSamples iterations, filling per-chain
  /// slabs of the (chain-major) results arrays; see Results. Chains execute
  /// on up to min(numThreads, numChains) worker threads.
  ///
  /// pollInterrupt, if set, is called only on this (the caller's) thread to
  /// ask whether the user requested cancellation; run returns true if it
  /// stopped early, having joined every worker first so nothing outlives the
  /// call. It is called at most every ~100ms, so a normal run pays only a
  /// clock read per sweep and stays bitwise identical. With worker chains the
  /// call returns as soon as the last chain finishes, not at the next poll,
  /// so a one-sweep run in a host's outer loop costs no added latency.
  ///
  /// onSweep, if set, is the host's per-sweep conditioning hook; it runs only
  /// when chains run inline (min(numThreads, numChains) <= 1), so the caller
  /// must not set it alongside worker-thread chains. Chains then run
  /// sequentially, so onSweep sees chain c completed before chain c + 1 begins.
  ///
  /// onDraw, if set, is the host's per-draw observer; unlike onSweep it is
  /// forwarded UNCHANGED on both paths - no lock, no queue, no hop to this
  /// thread - so it must be callable from a worker thread. A nonzero return
  /// aborts the run exactly as an interrupt does, and sets *stoppedByCallback
  /// (when non-null) so the caller can tell the two apart on a true return.
  /// An abort leaves the sampler INCONSISTENT with the results: this call
  /// returns before advancing the sample cursors, while storeSavedTreeRecord
  /// has already written saved trees into the slots those cursors count, so
  /// the caller must discard both the results and any saved trees.
  ///
  /// An exception a chain throws (an allocation failure in a monotone order
  /// count) stops every chain as a cancel does and leaves this call after
  /// the join, on the caller's thread.
  bool run(size_t numBurnIn, size_t numSamples, Results& results,
           const std::function<bool()>& pollInterrupt = {},
           const SweepCallback& onSweep = {},
           const DrawHook& onDraw = DrawHook(),
           bool* stoppedByCallback = nullptr) {
    size_t numChains = chains_.size();
    if (stoppedByCallback != nullptr) *stoppedByCallback = false;
    // the census covers THIS run, so a caller reading it after run() sees the
    // draws it just took and not everything the sampler has ever done
    for (auto& chain : chains_) chain->resetGPFallbackTally();
    for (auto& chain : chains_) chain->resetSlowCountTally();
    for (auto& chain : chains_) chain->setSavedSlotBase(currentSampleNum_);
    // the per-observation fits carry numReportedLocations channels per sample
    // (one everywhere but a multi-location combiner), so the per-chain slab
    // stride folds it in; L = 1 leaves the exact current chain-major stride
    size_t numLocations = results.numReportedLocations;
    // the per-sample varcount slab carries numVariableCountForests forests (1
    // for a single-forest model, up to K for multinomial and for a multi-forest
    // amplitude model), so the per-chain varcount stride folds it in; count 1
    // leaves the exact current chain-major stride. The caller declares it, and
    // the clamp to what the coupling can actually report happens HERE, once:
    // the stride below and storeSample's writes both read this one value, so
    // they agree by construction and storeSample needs no clamp of its own.
    size_t numVarCountForests = results.numVariableCountForests;
    if (numVarCountForests > numVariableCountForests())
      numVarCountForests = numVariableCountForests();
    // the per-sample threshold slab carries numOrdinalThresholds entries (0 off
    // ordinal), so the per-chain threshold stride folds it in
    size_t numOrdinalThresholds = results.numOrdinalThresholds;
    // per-draw strides for the per-observation channels, resolved ONCE from the
    // caller's fields: a chain's slab is numSamples of them, except at a stride
    // of zero, where the chain owns a single draw's buffer and the C buffers
    // are laid out contiguously. storeSample resolves the same field the same
    // way, which is why the stride is a field rather than a repeated formula.
    size_t n = data_.numObservations, nTest = data_.numTestObservations;
    size_t trainDraw = n * numLocations, testDraw = nTest * numLocations;
    size_t forestDraw = n * numForests();
    size_t trainStride = Results::strideOr(results.trainingFitsStride, trainDraw);
    size_t testStride = Results::strideOr(results.testFitsStride, testDraw);
    size_t varianceStride = Results::strideOr(results.varianceFitsStride, n);
    size_t varianceTestStride =
      Results::strideOr(results.varianceTestFitsStride, nTest);
    size_t forestStride =
      Results::strideOr(results.forestFitsStride, forestDraw);
    auto chainSlab = [numSamples](size_t stride, size_t drawExtent) {
      return stride == 0 ? drawExtent : numSamples * stride;
    };
    std::vector<Results> chainResults(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      Results& r(chainResults[c]);
      r.numReportedLocations = numLocations;
      r.numVariableCountForests = numVarCountForests;
      r.numOrdinalThresholds = numOrdinalThresholds;
      r.trainingFitsStride = trainStride;
      r.testFitsStride = testStride;
      r.varianceFitsStride = varianceStride;
      r.varianceTestFitsStride = varianceTestStride;
      r.forestFitsStride = forestStride;
      if (results.sigma != nullptr) r.sigma = results.sigma + c * numSamples;
      if (results.k != nullptr) r.k = results.k + c * numSamples;
      if (results.trainingFits != nullptr)
        r.trainingFits =
          results.trainingFits + c * chainSlab(trainStride, trainDraw);
      if (results.testFits != nullptr)
        r.testFits = results.testFits + c * chainSlab(testStride, testDraw);
      if (results.variableCounts != nullptr)
        r.variableCounts = results.variableCounts +
          c * numSamples * data_.numPredictors * numVarCountForests;
      if (results.splitProbabilities != nullptr)
        r.splitProbabilities =
          results.splitProbabilities + c * numSamples * data_.numPredictors;
      if (results.logLikelihood != nullptr)
        r.logLikelihood =
          results.logLikelihood + c * numSamples * data_.numObservations;
      if (results.ordinalThresholds != nullptr)
        r.ordinalThresholds =
          results.ordinalThresholds + c * numSamples * numOrdinalThresholds;
      if (results.shape != nullptr)
        r.shape = results.shape + c * numSamples;
      if (results.residualDf != nullptr)
        r.residualDf = results.residualDf + c * numSamples;
      if (results.varianceFits != nullptr)
        r.varianceFits =
          results.varianceFits + c * chainSlab(varianceStride, n);
      if (results.varianceTestFits != nullptr)
        r.varianceTestFits =
          results.varianceTestFits + c * chainSlab(varianceTestStride, nTest);
      // the per-forest slab carries every forest's own fits per sample, and the
      // glue slab the ragged amplitude vector; both are filled only by a
      // coupling that defines them, so the strides are reached on that path
      // alone
      if (results.forestFits != nullptr)
        r.forestFits =
          results.forestFits + c * chainSlab(forestStride, forestDraw);
      if (results.glue != nullptr)
        r.glue = results.glue + c * numSamples * totalAmplitudes();
    }

    if (options_.verbose) ext_printf("Running mcmc loop:\n");
    std::chrono::steady_clock::time_point startTime =
      std::chrono::steady_clock::now();

    size_t numWorkers = options_.numThreads < numChains ? options_.numThreads
                                                        : numChains;
    // the observer is forwarded to every chain on BOTH paths, unwrapped except
    // for the relay that records its stop; it is a plain function pointer over
    // plain arrays, so a worker thread can call it where onSweep cannot
    DrawRelay drawRelay;
    drawRelay.hook = onDraw;
    DrawHook relayHook;
    if (onDraw.fn != nullptr) {
      relayHook.fn = &DrawRelay::call;
      relayHook.context = &drawRelay;
    }
    const DrawHook* onDrawPtr = onDraw.fn != nullptr ? &relayHook : nullptr;
    bool cancelled = false;
    if (numWorkers <= 1) {
      // chains run on the main thread, so progress prints directly and the
      // interrupt poll runs inline between sweeps, throttled to ~100ms so a
      // fast problem is not dominated by the poll
      DirectProgressSink progress;
      // start one interval in the past so the first sweep polls immediately
      auto lastPoll =
        std::chrono::steady_clock::now() - std::chrono::milliseconds(100);
      std::function<bool()> shouldCancel;
      if (pollInterrupt) {
        shouldCancel = [&pollInterrupt, &lastPoll]() -> bool {
          auto now = std::chrono::steady_clock::now();
          if (now - lastPoll < std::chrono::milliseconds(100)) return false;
          lastPoll = now;
          return pollInterrupt();
        };
      }
      const std::function<bool()>* shouldCancelPtr =
        shouldCancel ? &shouldCancel : nullptr;
      const SweepCallback* onSweepPtr = onSweep ? &onSweep : nullptr;
      for (size_t c = 0; c < numChains && !cancelled; ++c)
        cancelled = chains_[c]->run(numBurnIn, numSamples, chainResults[c],
                                    &progress, c, shouldCancelPtr, onSweepPtr,
                                    onDrawPtr);
    } else {
      // workers never call into R: progress lines queue and the main thread
      // flushes them and polls for interrupts at least every 0.1 seconds,
      // setting the cancel flag the workers read (a relaxed atomic load per
      // sweep)
      QueuedProgressSink progress;
      // the caller blocks on chainsDone rather than on a timer: a worker takes
      // chainsMutex for its decrement, and wait_for re-tests the count under
      // that same lock before it sleeps, so a chain that finishes while the
      // lock is free - it is released for the flush and the poll - is seen
      // rather than slept through, and the last chain's notify cannot be
      // missed. The timeout that remains bounds only the progress flush and
      // the interrupt poll, both main-thread-only work.
      std::mutex chainsMutex;
      std::condition_variable chainsDone;
      size_t numChainsRunning = numChains;
      std::atomic<bool> cancelFlag(false);
      // the first exception a chain throws, rethrown on this thread after
      // the join; guarded by chainsMutex
      std::exception_ptr firstError;
      std::function<bool()> workerCancel = [&cancelFlag]() {
        return cancelFlag.load(std::memory_order_relaxed);
      };
      // the workers are spawned with SIGINT blocked so a Ctrl-C is delivered
      // only to this (the main) thread, whose poll turns it into a
      // cooperative cancel; a worker running R's interrupt handler could
      // longjmp across threads. (On Windows R's console Ctrl-C already runs
      // on the main thread, so no masking is needed.) A failed spawn cancels
      // the chains already started, which stop at their next sweep boundary.
      std::vector<std::thread> workers;
      auto worker = [this, numWorkers, numChains, numBurnIn, numSamples,
                     &chainResults, &progress, &chainsMutex, &chainsDone,
                     &numChainsRunning, &workerCancel, &cancelFlag,
                     &firstError, onDrawPtr](size_t w) {
        for (size_t c = w; c < numChains; c += numWorkers) {
          // a chain that stopped itself - its observer returned nonzero -
          // publishes the stop, so every other chain, this worker's own
          // remaining ones included, sees it at its next sweep boundary.
          // A chain that throws is stopped the same way, its exception
          // kept for the caller. The catch must stay inside this body:
          // the decrement below has to run for every chain, or the wait
          // for numChainsRunning never ends.
          try {
            if (chains_[c]->run(numBurnIn, numSamples, chainResults[c],
                                &progress, c, &workerCancel, nullptr,
                                onDrawPtr))
              cancelFlag.store(true, std::memory_order_relaxed);
          } catch (...) {
            {
              std::lock_guard<std::mutex> lock(chainsMutex);
              if (!firstError) firstError = std::current_exception();
            }
            cancelFlag.store(true, std::memory_order_relaxed);
          }
          bool last;
          {
            std::lock_guard<std::mutex> lock(chainsMutex);
            last = --numChainsRunning == 0;
          }
          // notified with the mutex released: the one waiter wakes to an
          // unheld lock, and the count it re-tests is already zero
          if (last) chainsDone.notify_one();
        }
      };
      spawnWorkers(workers, numWorkers, worker, [&cancelFlag]() {
        cancelFlag.store(true, std::memory_order_relaxed);
      });
      {
        std::unique_lock<std::mutex> lock(chainsMutex);
        while (numChainsRunning > 0) {
          // the flush prints and the poll may re-enter the host (R's interrupt
          // check), so both run off the lock, exactly as the timer tick did;
          // the wait re-tests the count on reacquiring it
          lock.unlock();
          if (options_.verbose) progress.flush();
          if (pollInterrupt && !cancelFlag.load(std::memory_order_relaxed) &&
              pollInterrupt())
            cancelFlag.store(true, std::memory_order_relaxed);
          lock.lock();
          chainsDone.wait_for(lock, std::chrono::milliseconds(100),
                              [&numChainsRunning]() {
                                return numChainsRunning == 0;
                              });
        }
      }
      for (std::thread& worker : workers) worker.join();
      if (options_.verbose) progress.flush();
      if (firstError) std::rethrow_exception(firstError);
      cancelled = cancelFlag.load(std::memory_order_relaxed);
    }

    if (stoppedByCallback != nullptr)
      *stoppedByCallback = drawRelay.stopped.load(std::memory_order_relaxed);

    // a cancelled run is aborted by the caller, so skip the completion cursor
    // advance and terminal summary that assume every requested sample landed
    if (cancelled) return true;

    runningTime_ += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - startTime).count();

    size_t capacity = savedTreeCapacity();
    if (capacity > 0 && numSamples > 0) {
      currentSampleNum_ = (currentSampleNum_ + numSamples) % capacity;
      recordedDraws_ = recordedDraws_ + numSamples < capacity
        ? recordedDraws_ + numSamples : capacity;
    }

    if (options_.verbose) printTerminalSummary();
    return false;
  }

  /// End-of-run report: accumulated loop time, leaf counts per tree, and
  /// split-variable usage, all from the current state.
  void printTerminalSummary() const {
    ext_printf("total seconds in loop: %f\n", runningTime_);

    ext_printf("\nTree sizes, last iteration:\n");
    std::vector<int32_t> bottoms;
    for (size_t c = 0; c < chains_.size(); ++c) {
      size_t linePrintCount = 0;
      ext_printf("[%lu] ", static_cast<unsigned long>(c + 1));
      // the "[c] " prefix counts as two slots toward the 20-per-line wrap
      linePrintCount += 2;
      for (size_t t = 0; t < options_.numTrees; ++t) {
        bottoms.clear();
        chains_[c]->tree(t).fillBottom(0, bottoms);
        ext_printf("%lu ", static_cast<unsigned long>(bottoms.size()));
        if ((linePrintCount++ + 1) % 20 == 0) ext_printf("\n");
      }
      if ((linePrintCount % 20) != 0) ext_printf("\n");
    }
    ext_printf("\n");

    ext_printf("Variable Usage, last iteration (var:count):\n");
    std::vector<std::uint32_t> variableCounts(data_.numPredictors, 0);
    for (size_t c = 0; c < chains_.size(); ++c)
      for (size_t t = 0; t < options_.numTrees; ++t)
        chains_[c]->tree(t).countVariableUses(variableCounts.data());
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      ext_printf("(%lu: %u) ", static_cast<unsigned long>(j + 1),
                 variableCounts[j]);
      if ((j + 1) % 5 == 0) ext_printf("\n");
    }

    ext_printf("\nDONE BART\n\n");
  }

  // Saved trees (keepTrees) and prediction.

  size_t savedTreeCapacity() const {
    return chains_[0]->savedTreeCapacity();
  }
  size_t currentSampleNum() const { return currentSampleNum_; }
  /// Recorded draws the store retains: min(draws recorded since the last
  /// reset, capacity). Every saved-tree read reports exactly this many draws,
  /// so a store still filling reports what it holds instead of padding with
  /// slots nothing has written - which read as zero-leaf trees, not as an
  /// absence.
  size_t filledSavedDraws() const { return recordedDraws_; }
  /// The slot holding output draw drawIndex, in [0, filledSavedDraws()). The
  /// retained draws are read OLDEST FIRST, ending at the slot before the write
  /// cursor, so the order is chronological however many run() calls filled the
  /// ring; on a full store this is (currentSampleNum + drawIndex) % capacity.
  size_t savedSlotForDraw(size_t drawIndex) const {
    size_t capacity = savedTreeCapacity();
    if (capacity == 0) return 0;
    return (currentSampleNum_ + capacity - recordedDraws_ + drawIndex) %
      capacity;
  }
  /// The retained draws as the ring slots holding them, for a chain call that
  /// rewrites them.
  SavedDrawSlots savedDrawSlots() const {
    return SavedDrawSlots{savedSlotForDraw(0), recordedDraws_};
  }

  const std::vector<FlatNode>& savedTree(size_t chainNum, size_t slot,
                                         size_t treeNum,
                                         size_t forestIndex = 0) const {
    return chains_[chainNum]->savedTree(slot, treeNum, forestIndex);
  }
  const std::vector<double>& savedTreeSlopes(size_t chainNum, size_t slot,
                                             size_t treeNum,
                                             size_t forestIndex = 0) const {
    return chains_[chainNum]->savedTreeSlopes(slot, treeNum, forestIndex);
  }
  const std::vector<std::uint64_t>& savedTreeMasks(
    size_t chainNum, size_t slot, size_t treeNum,
    size_t forestIndex = 0) const {
    return chains_[chainNum]->savedTreeMasks(slot, treeNum, forestIndex);
  }
  /// slopes, when non-null, receives a vector-parameter tree's slopes
  /// (numParams - 1 per leaf, pre-order); cleared for scalar leaf models.
  /// masks, when non-null, receives the wide categorical side channel.
  /// forestIndex selects the forest to read (0 for every non-BCF sampler).
  void flattenTree(size_t chainNum, size_t treeNum,
                   std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes = nullptr,
                   std::vector<std::uint64_t>* masks = nullptr,
                   size_t forestIndex = 0) {
    chains_[chainNum]->flattenTree(treeNum, nodes, counts, slopes, masks,
                                   forestIndex);
  }

  /// The traversal count below which an out-of-sample replay runs inline on
  /// the caller's thread. A traversal is one (row, tree, slab) descent. Each
  /// entry point counts ITS OWN traversals, not forest 0's.
  ///
  /// Nothing breaks either side: the partition gives each slab its own output
  /// range and reduces nothing across workers, so a replay is bitwise
  /// identical at any worker count and the cutoff buys only time. The value
  /// is the measured crossover: spawn and join cost tens of microseconds, so
  /// the fan-out pays once a replay's serial work is a small multiple of that.
  /// A caller who measures otherwise on their own machine moves it through
  /// SamplerOptions::predictParallelCutoff; this is the default that field's
  /// 0 selects.
  static constexpr size_t predictParallelCutoff = 50000;

  /// Every forest's tree count summed, the traversal weight of a replay that
  /// walks all of them (multi-location and per-forest).
  size_t totalTreesAcrossForests() const {
    size_t total = 0;
    for (size_t f = 0; f < numForests(); ++f) total += numTreesInForest(f);
    return total;
  }

  /// Runs body(slab, scratch) for every slab in [0, numSlabs), block-
  /// partitioned contiguously over min(resolved threads, numSlabs) workers.
  ///
  /// numThreads == 0 means the sampler's own count; a resolved count below 1 -
  /// which a sampler whose count was set to 0 produces, that setter storing
  /// what it is given - is floored at 1, because zero workers would leave the
  /// caller's output buffer entirely unwritten. traversals is this entry's own
  /// (row x tree x slab) descent count, and below predictParallelCutoff the
  /// replay stays inline rather than paying for a fan-out it cannot amortize.
  ///
  /// Bitwise identity at every thread count is a property of the partition,
  /// not a tolerance: each slab owns its output range whole, the per-row tree
  /// sum runs in tree order inside one slab, and nothing is reduced across
  /// workers, so no addend order depends on how the slabs were dealt out.
  ///
  /// Scratch is one struct per WORKER, not one per slab: the replay's buffers
  /// grow to their working size on a worker's first slab and are reused for
  /// the rest of its block, so nothing is shared between workers and the
  /// allocation count is bounded by the worker count rather than by the slab
  /// count. Worker bodies must not call into R: an
  /// exception escaping a std::thread body is std::terminate, and a large
  /// replay is exactly where bad_alloc is plausible, so each body catches and
  /// the first message is rethrown after the join, on the caller's thread.
  ///
  /// That rethrow is a C++ exception and not an R error on purpose. Raising
  /// through R here would longjmp from the bottom of the replay, past this
  /// frame's worker buffers and past the rank bitmaps a sparse view built two
  /// frames up; throwing unwinds both. Every entrance that can reach this
  /// converts the exception back into an R error at its own boundary, where
  /// nothing is left to free.
  ///
  /// There is no interrupt poll inside the fan-out. R's unix SIGINT handler
  /// only sets a pending flag; the longjmp happens in R_CheckUserInterrupt,
  /// which no replay path calls, so a Ctrl-C during a join cannot unwind out
  /// of a worker. Blocking SIGINT across the spawn keeps the signal on this
  /// thread regardless, as run() does.
  template <typename Body>
  void fanOutPredictSlabs(size_t numSlabs, size_t numThreads,
                          size_t traversals, Body&& body) {
    size_t resolved = numThreads != 0 ? numThreads : options_.numThreads;
    if (resolved < 1) resolved = 1;
    predictPartition.resolvedThreads = resolved;
    if (numSlabs == 0) {
      predictPartition.numWorkers = 0;
      predictPartition.workerForSlab.clear();
      return;
    }
    // the test seam wins over the caller's setting, which wins over the
    // calibrated default
    size_t cutoff = predictPartition.cutoffOverride != 0
      ? predictPartition.cutoffOverride
      : (options_.predictParallelCutoff != 0 ? options_.predictParallelCutoff
                                             : predictParallelCutoff);
    size_t numWorkers = resolved < numSlabs ? resolved : numSlabs;
    if (traversals < cutoff) numWorkers = 1;

    size_t base = numSlabs / numWorkers, remainder = numSlabs % numWorkers;
    predictPartition.numWorkers = numWorkers;
    predictPartition.workerForSlab.assign(numSlabs, 0);
    for (size_t w = 0, begin = 0; w < numWorkers; ++w) {
      size_t end = begin + base + (w < remainder ? 1 : 0);
      for (size_t slab = begin; slab < end; ++slab)
        predictPartition.workerForSlab[slab] = w;
      begin = end;
    }

    if (numWorkers == 1) {
      PredictScratch scratch;
      for (size_t slab = 0; slab < numSlabs; ++slab) body(slab, scratch);
      return;
    }

    std::vector<PredictScratch> scratch(numWorkers);
    std::vector<char> failed(numWorkers, 0);
    std::vector<std::string> firstError(numWorkers);
    std::vector<std::thread> workers;
    // worker w's block starts after every earlier worker's, the first
    // `remainder` of which carry one extra slab
    auto worker = [&body, &scratch, &failed, &firstError, base,
                   remainder](size_t w) {
      size_t begin = w * base + (w < remainder ? w : remainder);
      size_t end = begin + base + (w < remainder ? 1 : 0);
      try {
        for (size_t slab = begin; slab < end; ++slab) body(slab, scratch[w]);
      } catch (const std::exception& error) {
        failed[w] = 1;
        firstError[w] = error.what();
      } catch (...) {
        failed[w] = 1;
        firstError[w] = "unknown exception";
      }
    };
    // a failed spawn lets the started workers finish their own blocks
    spawnWorkers(workers, numWorkers, worker, []() {});
    for (std::thread& worker : workers) worker.join();
    for (size_t w = 0; w < numWorkers; ++w)
      if (failed[w] != 0)
        throw std::runtime_error("predict worker " + std::to_string(w) +
                                 " failed: " + firstError[w]);
  }

  /// Fits for raw column-major test rows, on the original response scale
  /// (offsets are the caller's problem). With saved trees, out is
  /// numTestObservations (x numReportedLocations) x filledSavedDraws() x
  /// numChains, chain-major like Results, the retained draws oldest first;
  /// without, one slab per chain from the live trees. A multi-location combiner (multinomial: K softmax channels)
  /// inserts the K location dimension between the observations and the slots;
  /// L = 1 keeps the exact numTestObservations-per-slot byte layout.
  void predict(const double* x_test, size_t numTestObservations,
               size_t numThreads, double* out) {
    predictColumns(DenseColumns{x_test, numTestObservations},
                   numTestObservations, nullptr, numThreads, out);
  }

  /// The same over a borrowed view: a dense block reads through the reader the
  /// raw entry builds, so its replay is unchanged, while a CSC-backed source
  /// routes rows off rank bitmaps built here, for this call only. Nothing is
  /// retained - the view is the caller's, and the sampler's own test store is
  /// untouched.
  ///
  /// categoryOffset, when non-null, is the caller's per-call nNew x K offset
  /// for a multi-location surface, added to the raw fits before the softmax.
  /// It is never read from the sampler: predict's rows are the caller's, so
  /// its offset must be too.
  void predict(const PredictorSource& source, size_t numTestObservations,
               const double* categoryOffset, size_t numThreads, double* out) {
    if (source.isDenseBlock())
      predictColumns(DenseColumns{source.denseValues, numTestObservations},
                     numTestObservations, categoryOffset, numThreads, out);
    else
      predictColumns(PredictorSourceColumns(source, data_.types.data()),
                     numTestObservations, categoryOffset, numThreads, out);
  }

  template <typename Columns>
  void predictColumns(const Columns& columns, size_t numTestObservations,
                      const double* categoryOffset, size_t numThreads,
                      double* out) {
    size_t capacity = savedTreeCapacity();
    size_t numDraws = recordedDraws_;
    size_t numLocations = numReportedLocations();
    size_t slab = numTestObservations * numLocations;
    size_t numChains = chains_.size();
    // a multi-location replay walks every forest's trees, a single-location
    // one only forest 0's
    size_t treesPerSlab =
      numLocations > 1 ? totalTreesAcrossForests() : numTreesInForest(0);
    if (capacity > 0) {
      // slab s is chain s / numDraws at draw s % numDraws, which is exactly
      // the (c * numDraws + i) offset the output is laid out by
      fanOutPredictSlabs(
        numChains * numDraws, numThreads,
        numChains * numDraws * treesPerSlab * numTestObservations,
        [&](size_t s, PredictScratch& scratch) {
          size_t c = s / numDraws;
          size_t slot = savedSlotForDraw(s - c * numDraws);
          double* dst = out + s * slab;
          if (numLocations > 1)
            chains_[c]->predictFromSavedSampleMulti(
              slot, columns, numTestObservations, categoryOffset, scratch, dst);
          else
            chains_[c]->predictFromSavedSample(slot, columns,
                                               numTestObservations, scratch,
                                               dst);
        });
    } else {
      // without a store there is exactly ONE slab per chain, so no two workers
      // reach the same Chain: the live-tree replay flattens through the
      // chain's own buffers, and a finer (row-axis) partition would have to
      // privatize them first
      assert(numChains > 0);
      fanOutPredictSlabs(
        numChains, numThreads,
        numChains * treesPerSlab * numTestObservations,
        [&](size_t c, PredictScratch& scratch) {
          if (numLocations > 1)
            chains_[c]->predictFromCurrentTreesMulti(
              columns, numTestObservations, categoryOffset, scratch,
              out + c * slab);
          else
            chains_[c]->predictFromCurrentTrees(columns, numTestObservations,
                                                scratch, out + c * slab);
        });
    }
  }

  /// Per-forest RAW fits for new rows of a borrowed predictor view: forest f's
  /// own internal-scale total, with no amplitude glue, no response transform
  /// and no offset (Chain::predictPerForestFromSavedSample states why). out is
  /// numTestObservations x numForests x filledSavedDraws() x numChains,
  /// chain-major like predict's and in the same oldest-first draw order;
  /// without saved trees, one slab per chain from the live trees. Nothing resident is read or retained: the rows and the
  /// recombination both belong to the caller.
  void predictPerForest(const PredictorSource& source,
                        size_t numTestObservations, size_t numThreads,
                        double* out) {
    if (source.isDenseBlock())
      predictPerForestColumns(
        DenseColumns{source.denseValues, numTestObservations},
        numTestObservations, numThreads, out);
    else
      predictPerForestColumns(PredictorSourceColumns(source, data_.types.data()),
                              numTestObservations, numThreads, out);
  }

  /// Dense convenience spelling: a plain column-major block of new rows.
  void predictPerForest(const double* x_test, size_t numTestObservations,
                        size_t numThreads, double* out) {
    predictPerForestColumns(DenseColumns{x_test, numTestObservations},
                            numTestObservations, numThreads, out);
  }

  template <typename Columns>
  void predictPerForestColumns(const Columns& columns,
                               size_t numTestObservations, size_t numThreads,
                               double* out) {
    size_t capacity = savedTreeCapacity();
    size_t numDraws = recordedDraws_;
    size_t numChains = chains_.size();
    size_t slab = numTestObservations * numForests();
    size_t treesPerSlab = totalTreesAcrossForests();
    if (capacity > 0) {
      fanOutPredictSlabs(
        numChains * numDraws, numThreads,
        numChains * numDraws * treesPerSlab * numTestObservations,
        [&](size_t s, PredictScratch& scratch) {
          size_t c = s / numDraws;
          chains_[c]->predictPerForestFromSavedSample(
            savedSlotForDraw(s - c * numDraws), columns, numTestObservations,
            scratch, out + s * slab);
        });
    } else {
      // one slab per chain, as in predictColumns and for the same reason
      assert(numChains > 0);
      fanOutPredictSlabs(
        numChains, numThreads,
        numChains * treesPerSlab * numTestObservations,
        [&](size_t c, PredictScratch& scratch) {
          chains_[c]->predictPerForestFromCurrentTrees(
            columns, numTestObservations, scratch, out + c * slab);
        });
    }
  }

  /// Whether this sampler carries a heteroscedastic variance forest, and its
  /// tree count; the run/predict bridges gate the s(x) channels on this.
  bool hasVarianceForest() const { return chains_[0]->hasVarianceForest(); }
  size_t numVarianceTrees() const { return chains_[0]->numVarianceTrees(); }

  /// Chain chainNum's CURRENT variance surface on the original scale, training
  /// or test rows; Chain::currentVarianceFits states the contract and owns the
  /// two refusals.
  bool currentVarianceFits(size_t chainNum, bool test, double* out) {
    return chains_[chainNum]->currentVarianceFits(test, out);
  }

  /// The variance surface s^2(x) for raw column-major new rows, original scale,
  /// mirroring predict: out is numTestObservations x filledSavedDraws() x
  /// numChains, chain-major, same draw order. Requires saved trees
  /// (keepTrees).
  void predictVariance(const double* x_test, size_t numTestObservations,
                       size_t numThreads, double* out) {
    predictVarianceColumns(DenseColumns{x_test, numTestObservations},
                           numTestObservations, numThreads, out);
  }

  /// The variance twin of the view-taking predict.
  void predictVariance(const PredictorSource& source,
                       size_t numTestObservations, size_t numThreads,
                       double* out) {
    if (source.isDenseBlock())
      predictVarianceColumns(
        DenseColumns{source.denseValues, numTestObservations},
        numTestObservations, numThreads, out);
    else
      predictVarianceColumns(PredictorSourceColumns(source, data_.types.data()),
                             numTestObservations, numThreads, out);
  }

  template <typename Columns>
  void predictVarianceColumns(const Columns& columns,
                              size_t numTestObservations, size_t numThreads,
                              double* out) {
    size_t numDraws = recordedDraws_;
    size_t numChains = chains_.size();
    fanOutPredictSlabs(
      numChains * numDraws, numThreads,
      numChains * numDraws * numVarianceTrees() * numTestObservations,
      [&](size_t s, PredictScratch& scratch) {
        size_t c = s / numDraws;
        chains_[c]->predictVarianceFromSavedSample(
          savedSlotForDraw(s - c * numDraws), columns, numTestObservations,
          scratch, out + s * numTestObservations);
      });
  }

  // State serialization: getState captures everything needed to reconstruct
  // the sampler's posterior state in a fresh instance over the same data;
  // setState validates every chain before touching anything, so failure
  // leaves the sampler unchanged. The cut grid and its weights are the
  // store's: a state records the ones it was read on, and only an install
  // told to adopt them writes them.

  void getState(SamplerStateData& state) {
    state.chains.resize(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->getState(state.chains[c]);
    state.cutPoints = data_.cutPoints;
    state.cutMass.clear();
    for (size_t j = 0; j < data_.numPredictors && state.cutMass.empty(); ++j)
      if (data_.cutsWeighted(j)) state.cutMass = data_.cutMass;
    state.currentSampleNum = currentSampleNum_;
    state.recordedDraws = recordedDraws_;
    state.missingColumns = data_.hasMissing;
  }

  /// Installs \p state's chains. The cut grid and its weights are the
  /// store's, and \p adoptGrid says which grid the state goes onto.
  ///
  /// With adoptGrid the store first takes the points and weights the state
  /// records on every column that is not a factor, those columns
  /// re-quantized against them, and the state is judged and installed there:
  /// the call of a sampler just created, which holds no grid of its own to
  /// keep. A factor column keeps the grid of its level table, which no
  /// install writes, whatever the state records on it. A state that is
  /// refused or declined puts the store's grid, weights and codes back.
  ///
  /// Without it no point and no weight of the store is written. A state
  /// whose recorded points are the store's is judged and installed as it is,
  /// whatever weights it records. So is one recorded on other points every
  /// live threshold of which, mean and variance, is a point of the store's
  /// grid (Chain::liveThresholdsAreOnGrid): its trees resolve there as
  /// stored. A state holding a live threshold the store's grid lacks cannot
  /// stand as stored and is declined without force; with force it is
  /// installed across the grids (installAcrossGrids): each tree is built and
  /// its rows routed on the recorded grid, then every split on a column whose
  /// grids differ is moved to the nearest point of the store's grid inside
  /// the interval the splits above it leave, a split left no point merging
  /// with everything beneath it into one leaf at the mean of the leaves
  /// beneath weighted by the rows each held, and a leaf the moved splits
  /// leave empty is merged as after a forced predictor update. That is the
  /// state installed on its own grid followed by setCutPoints back to the
  /// store's by value, with the store's weights kept. The kept draws are
  /// copied as stored either way and never moved.
  ///
  /// currentPredictors is the call-time predictor matrix a column whose grid
  /// is replaced re-quantizes from (data@x, or the retained creation spec's
  /// @x); null for CSC/mixed stores. Only an adopted grid other than the
  /// store's and an install across grids read it.
  ///
  /// A refusal - false with notClean down - is the same with and without
  /// force and leaves the sampler exactly as it was: a state of another
  /// shape, a malformed block or tree, a grid that repeats a point, latents
  /// the family cannot hold, saved-tree blocks that name different draw
  /// counts, a record of missable columns of another length, a live
  /// threshold that is no point of the grid the state itself records, and,
  /// without adoptGrid, on a column the store holds as an ordered factor,
  /// whose grid its level count fixes, a recorded grid other than the
  /// store's. With adoptGrid that record is not held against the state: no
  /// install writes a factor column's grid
  /// (ColumnStore::setCutPointsForColumn), so the store of a sampler just
  /// created keeps the grid of the levels its data declares, which a
  /// whole-data replacement can leave other than the one the sampler held
  /// when it stored the state, and the state is judged on the store's.
  /// lengthscaleRefused names a state whose saved gp draws were made
  /// under other lengthscales than the ones this sampler was given
  /// (Chain::lengthscaleStateFeasible), asked only where a saved draw would
  /// be kept. Everything is judged before any chain is touched, against the
  /// grid the state is installed on or, for an install across grids, the one
  /// it records.
  ///
  /// A state that records the columns that could hold a missing value where
  /// it was stored (missingColumns) is first completed for the columns that
  /// can here and could not there: every rule on one, in its live trees and
  /// its recorded saved draws, holds no drawn direction and takes one with
  /// probability one half, each chain's coins from the generator the state
  /// carries for it (Chain::completeMissingDirections). The draw is made on a
  /// copy before anything of the sampler is written, and everything below
  /// judges and installs the copy, its generators standing past the coins. A
  /// state without the record is installed as it is, as is one whose record
  /// names every such column; a column the record names and this sampler
  /// has never seen missing plays no part. No flag is raised.
  ///
  /// A state past those is installed as stored unless one of its live trees,
  /// mean or variance, cannot stand in this sampler as its flat form holds
  /// it:
  /// - a bottom node no row of the current data reaches, as a state stored
  ///   before setPredictor(forceUpdate) leaves;
  /// - a split outside the interval its ancestors leave;
  /// - a missing direction on a column, not pooled, that has never held a
  ///   missing value here;
  /// - a split on a column its forest's mask forbids (a moderator subset, a
  ///   blocks() row, a restricted variance forest), which
  ///   splitVariableLogProbability would price against a menu that drops the
  ///   column for as long as the split lives;
  /// - a split that breaks its forest's interaction constraint;
  /// - a monotone tree whose leaf values, as stored, leave the cone.
  /// With force such a tree is repaired on the way in: the bottom, the split
  /// outside its interval and the first forbidden split on each path from
  /// the root are merged into one leaf at the weighted mean of the leaves
  /// beneath (Chain::rebuildLiveForest, Chain::rebuildVarianceForest), the
  /// direction is dropped, the monotone tree is reseeded at 0. Without force
  /// the state is not installed: an adopted grid is put back, notClean, when
  /// non-null, is set, and false is returned with the chains, the store, the
  /// generators and the latents untouched. Only LIVE trees are held to any of
  /// this, a saved slot being a replay target routed over new rows rather
  /// than over this partition.
  ///
  /// No value is converted: trees, leaf values, k, sigma and the saved draws
  /// go in as stored, numbers on the internal scale read against this
  /// sampler's response transform, whatever transform the state records. The
  /// transform is model and is not moved.
  ///
  /// The kept draws are output, and a state installs whatever store it was
  /// taken from. A state whose store is the size of this sampler's has its
  /// ring copied slot for slot, write position included. From a store of
  /// another size the newest of the state's recorded draws that fit go to
  /// slots 0 on, oldest first, and the write position follows them; the
  /// capacity does not change, and a sampler keeping none takes none.
  ///
  /// adoptCapacity, when not keepStoreCapacity, resizes the store to that
  /// many samples before the kept draws go in, only once the state is
  /// accepted: a refused or declined state leaves the store and its draws as
  /// they were. An allocation failure in the resize throws with the cut grid
  /// restored and nothing else changed.
  ///
  /// altered, when non-null, reports whether a state was installed other
  /// than as stored: true exactly when a live tree was repaired as above or
  /// moved across grids, so a forced install sets it exactly when the same
  /// call without force declines. Only the installing build reports: the
  /// scratch builds that validate a state never do, and a refusal leaves it
  /// false.
  /// Not reported, since the chain is still the stored one: a value the
  /// sampler holds fixed or the state lacks, a generator of another kind left
  /// in place, a pooled categorical rule keeping a missing bit its column no
  /// longer routes, and the saved draws, which are copied.
  ///
  /// gridDiffers, when non-null, reports whether an install that does not
  /// adopt them put the state on cut points or weights other than the ones
  /// it records: the trees, or the tree prior they run under, are then not
  /// the stored chain's. An adopting install leaves it false, also where an
  /// ordered factor's column kept other points than the state records, and
  /// so does a state not installed.
  bool setState(const SamplerStateData& state,
                const double* currentPredictors,
                bool* lengthscaleRefused = nullptr,
                bool* altered = nullptr,
                size_t adoptCapacity = keepStoreCapacity,
                bool force = true,
                bool* notClean = nullptr,
                bool adoptGrid = true,
                bool* gridDiffers = nullptr) {
    if (lengthscaleRefused != nullptr) *lengthscaleRefused = false;
    if (altered != nullptr) *altered = false;
    if (notClean != nullptr) *notClean = false;
    if (gridDiffers != nullptr) *gridDiffers = false;
    bool installAltered = false;
    if (state.chains.size() != chains_.size()) return false;
    if (state.cutPoints.size() != data_.numPredictors) return false;
    bool stateWeighted = !state.cutMass.empty();
    if (stateWeighted && state.cutMass.size() != data_.numPredictors)
      return false;
    if (!state.missingColumns.empty() &&
        state.missingColumns.size() != data_.numPredictors)
      return false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.splitsBySubset(j)) {
        if (!state.cutPoints[j].empty()) return false;
        if (stateWeighted && !state.cutMass[j].empty()) return false;
      } else if (!cutGridIsValid(state.cutPoints[j].data(),
                                 state.cutPoints[j].size()) ||
                 (stateWeighted &&
                  !cutMassIsValid(state.cutMass[j],
                                  state.cutPoints[j].size()))) {
        return false;
      } else if (!adoptGrid && data_.types[j] == ColumnKind::orderedFactor &&
                 state.cutPoints[j] != data_.cutPoints[j]) {
        // the level count fixes the grid the store keeps, so another is
        // another column's
        return false;
      }
    }
    // the store the state was taken from, and the newest of its recorded
    // draws that the store it goes into holds
    size_t stateCapacity = 0;
    bool blocksAgree = true;
    for (size_t c = 0; c < chains_.size() && blocksAgree; ++c) {
      size_t chainCapacity;
      blocksAgree =
        chains_[c]->savedStateCapacity(state.chains[c], chainCapacity) &&
        (c == 0 || chainCapacity == stateCapacity);
      stateCapacity = chainCapacity;
    }
    size_t liveCapacity = savedTreeCapacity();
    bool resize =
      adoptCapacity != keepStoreCapacity && adoptCapacity != liveCapacity;
    size_t capacity = resize ? adoptCapacity : liveCapacity;
    size_t keptDraws = std::min({state.recordedDraws, stateCapacity, capacity});
    // blocks that disagree are refused below, after the lengthscales
    bool keepsDraws = blocksAgree ? keptDraws > 0 : state.recordedDraws > 0;
    for (size_t c = 0; c < chains_.size(); ++c) {
      if (!chains_[c]->lengthscaleStateFeasible(state.chains[c], keepsDraws)) {
        if (lengthscaleRefused != nullptr) *lengthscaleRefused = true;
        return false;
      }
    }
    // the directions the state was stored without, drawn on a copy of its
    // chains that stands in for them from here on. Nothing of the sampler
    // has been written and each chain's generator is back where it was, so
    // a refusal, a decline or a throw below has no draw to undo. A tree the
    // draw cannot walk leaves the state as given for the refusal below.
    size_t stateCursor =
      stateCapacity > 0 ? state.currentSampleNum % stateCapacity : 0;
    std::vector<ChainStateData> completed;
    if (blocksAgree && !state.missingColumns.empty()) {
      std::vector<std::uint8_t> raised(data_.numPredictors);
      bool draws = false;
      for (size_t j = 0; j < data_.numPredictors; ++j) {
        raised[j] =
          data_.hasMissing[j] != 0 && state.missingColumns[j] == 0 ? 1 : 0;
        draws = draws || raised[j] != 0;
      }
      for (size_t c = 0; c < chains_.size() && draws; ++c)
        draws = chains_[c]->flatTreesAreWalkable(state.chains[c]);
      if (draws) {
        completed = state.chains;
        size_t recorded = std::min(state.recordedDraws, stateCapacity);
        SavedDrawSlots stateKept;
        if (recorded > 0)
          stateKept = SavedDrawSlots{
            (stateCursor + stateCapacity - recorded) % stateCapacity, recorded};
        for (size_t c = 0; c < chains_.size(); ++c)
          chains_[c]->completeMissingDirections(completed[c], raised.data(),
                                                stateKept);
      }
    }
    const std::vector<ChainStateData>& chainStates =
      completed.empty() ? state.chains : completed;
    // where the recorded grid and the store's part: the points of a column,
    // or the weights of equal points
    bool pointsDiffer = false, massDiffers = false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.splitsBySubset(j)) continue;
      if (state.cutPoints[j] != data_.cutPoints[j])
        pointsDiffer = true;
      else if (stateWeighted ? state.cutMass[j] != data_.cutMass[j]
                             : data_.cutsWeighted(j))
        massDiffers = true;
    }
    // a state kept off its recorded points: a live threshold the record
    // lacks is the refusal, and one the store's grid lacks takes the install
    // across the grids
    bool acrossGrids = false;
    if (!adoptGrid && pointsDiffer) {
      for (size_t c = 0; c < chains_.size(); ++c) {
        if (!chains_[c]->liveThresholdsAreOnGrid(chainStates[c],
                                                 state.cutPoints))
          return false;
        acrossGrids = acrossGrids ||
          !chains_[c]->liveThresholdsAreOnGrid(chainStates[c],
                                               data_.cutPoints);
      }
    }
    // adopting: install the state's cuts, snapshotting for rollback, tree
    // validity being defined against them
    std::optional<GridSnapshot> adopted;
    if (adoptGrid) {
      adopted.emplace(data_);
      for (size_t j = 0; j < data_.numPredictors; ++j) {
        if (data_.splitsBySubset(j)) continue;
        const double* mass = stateWeighted && !state.cutMass[j].empty()
          ? state.cutMass[j].data() : nullptr;
        // a restored grid equal to the live one leaves the codes already
        // correct (the continuation contract), so skip its re-quantization
        // and its raw; its weights are the state's either way
        if (state.cutPoints[j] == data_.cutPoints[j]) {
          if (mass != nullptr)
            data_.cutMass[j] = state.cutMass[j];
          else
            data_.cutMass[j].clear();
          continue;
        }
        data_.setCutPointsForColumn(j, state.cutPoints[j].data(),
                                    static_cast<std::uint32_t>(
                                      state.cutPoints[j].size()),
                                    currentPredictors, mass);
      }
    }

    // validity for every chain, then the verdict where the call is not
    // forced: both judge the state against the grid it is installed on, and a
    // state refused for one chain is refused whatever another's trees need.
    // Across grids the trees resolve only on the recorded points, swapped in
    // for the builds, and no verdict is asked: every such state is moved
    bool allValid = blocksAgree;
    bool needsRepair = acrossGrids;
    if (acrossGrids) {
      ScopedCutGrid recorded(data_, state.cutPoints);
      for (size_t c = 0; c < chains_.size() && allValid; ++c)
        allValid = chains_[c]->stateIsValid(chainStates[c]);
    } else {
      for (size_t c = 0; c < chains_.size() && allValid; ++c)
        allValid = chains_[c]->stateIsValid(chainStates[c],
                                            force ? nullptr : &needsRepair);
    }

    if (!allValid || (needsRepair && !force)) {
      if (adopted) adopted->restore(data_);
      if (notClean != nullptr) *notClean = allValid;
      return false;
    }
    if (resize) {
      try {
        resizeSavedTrees(adoptCapacity);
      } catch (...) {
        if (adopted) adopted->restore(data_);
        throw;
      }
      options_.keepTrees = adoptCapacity > 0;
      options_.numSamplesToStore = adoptCapacity;
    }

    // a store of the state's own size continues its ring; any other takes
    // the newest draws at its start (Chain::setState)
    bool repack = stateCapacity != capacity;
    SavedDrawSlots newest;
    if (repack && keptDraws > 0)
      newest = SavedDrawSlots{
        (stateCursor + stateCapacity - keptDraws) % stateCapacity, keptDraws};
    auto installChains = [&]() {
      for (size_t c = 0; c < chains_.size(); ++c)
        if (!chains_[c]->setState(chainStates[c], &installAltered, newest))
          return false;
      return true;
    };
    if (acrossGrids) {
      if (!installAcrossGrids(state.cutPoints, currentPredictors,
                              installChains))
        return false;
      installAltered = true;
    } else if (!installChains()) {
      return false;
    }
    currentSampleNum_ =
      capacity > 0 ? (repack ? keptDraws : state.currentSampleNum) % capacity
                   : 0;
    recordedDraws_ = keptDraws;
    if (altered != nullptr) *altered = installAltered;
    if (gridDiffers != nullptr)
      *gridDiffers = !adoptGrid && (pointsDiffer || massDiffers);
    return true;
  }

  /// Warm start: seed each chain's live forest from a donor sample. sampleMap
  /// selects, for destination chain c, a donor (chain, slot); slot < 0 takes
  /// the donor chain's live trees, else its saved slot (slot-major, one forest
  /// per slot). Only trees, DART, and the sigma, k, concentration and
  /// amplitudes this sampler draws transfer - rng and auxiliary state stay
  /// fresh, so each chain evolves independently from its own stream. The cut
  /// grid, its weights and the codes are the store's, the same before and
  /// after, byte for byte. The donor must share this sampler's forest count;
  /// its own tree counts set the stride of a saved slot.
  ///
  /// Parts and leaves. DART weights go only between two DART fits and a
  /// variance forest only where both have one: a part the donor alone has
  /// is not read, and one this sampler alone has starts as a new fit's once
  /// the donor is in (Chain::startDartFresh, Chain::startVarianceForestFresh),
  /// neither judged. \p donorLeaves is the donor's leaf kind as its own blocks
  /// give it, null for this sampler's own, and each tree's block is held to
  /// it. gp leaves go only to gp leaves on as many rows. A constant donor
  /// into linear leaves takes slopes of 0, which is clean; a linear one into
  /// constant leaves, or into linear leaves on another number of covariates,
  /// leaves its slopes behind, each leaf at the intercept it stores and not
  /// restated, the donor's standardization being of other columns.
  ///
  /// Another tree count. A forest of this sampler holding more trees than
  /// the donor's takes all of them, in the donor's order, and then single
  /// leaves at 0, with slopes and fits of 0, or at 1 in the variance forest,
  /// as a new fit starts them; that is no change. Under blocks(), whatever
  /// the count, each block in turn takes in that order the trees its slots
  /// hold as stored (Chain::placeDonorTrees). A forest holding fewer, or
  /// under blocks one that leaves a tree out, is not clean: forced, it takes
  /// the trees whose fits vary most over this sampler's rows, chosen by
  /// trial installs (chooseDonorTrees), in the donor's order and as
  /// installed. Each chain chooses among the trees of the donor sample it
  /// takes, the mean and the variance forest apart. The directions below
  /// are drawn once, over every tree of the sample, and every round reads
  /// them.
  ///
  /// How a donor is read, none of it a change to it. A donor that records
  /// the columns that could hold a missing value where its trees were grown
  /// (missingColumns) holds no direction on the others: for each of those
  /// this sampler can hold one in, every rule on it in the trees a chain
  /// takes draws its direction with probability one half from that chain's
  /// own generator, before the trees are judged
  /// (Chain::drawDonorMissingDirections); a donor with no record draws
  /// nothing. Leaf values, k and sigma go in as the donor stored them,
  /// numbers on the internal scale read against this sampler's response
  /// transform, never converted from the donor's; a linear leaf's
  /// coefficients drawn under another covariate standardization are restated
  /// in this sampler's, whose standardization and lengthscales stay its own.
  /// Under a variance forest the scale surface rides along, taken from the
  /// same slot as the mean forest.
  ///
  /// A donor is clean when every tree a chain takes stands in this sampler
  /// as the donor holds it and every leaf is the donor's. It is not where
  /// - a live threshold, mean or variance, is no point of the store's grid;
  /// - a tree would be installed other than as its flat form holds it, by
  ///   the staging and for the reasons a state's verdict gives
  ///   (Chain::donorTreesStage, Sampler::setState), judged after the
  ///   directions are drawn;
  /// - the leaves are function leaves and the donor's recorded points are
  ///   not the store's, equal points being the one sign that the fits, values
  ///   by row, are of these rows: the fits then start at 0;
  /// - a linear donor's slopes are left behind;
  /// - a linear leaf's covariate holds one value in the donor's data and
  ///   another single value here under a slope that is not 0, the
  ///   coefficients being kept as built (Chain::convertDonorStandardization);
  /// - a forest holds fewer trees than the donor's or, under blocks, leaves
  ///   one out.
  /// Without \p force such a donor is not installed: declined is returned
  /// with nothing written and the generators back where they were. With it
  /// the donor goes in with the changes. Where every live threshold is a
  /// point of the store's grid the trees are built there and repaired as a
  /// forced state install repairs them (Chain::rebuildLiveForest,
  /// Chain::rebuildVarianceForest). Otherwise they go in across the grids
  /// (installAcrossGrids): built and routed on the donor's recorded points,
  /// then every split on a column whose grids differ moved to the nearest
  /// point of the store's grid inside the interval the splits above it
  /// leave, a split left no point and a leaf left no row merging at the mean
  /// of the leaves beneath weighted by the rows each held. On either arm a
  /// monotone tree out of order is collapsed where the order breaks
  /// (collapseMonotoneTree), where a state's is reset. A clean donor is
  /// installed the same with and without force.
  ///
  /// \p currentPredictors is the call-time predictor matrix the columns
  /// re-quantize from on the way across grids (null for CSC/mixed stores),
  /// read by no other install.
  ///
  /// A refusal is the same in both forms and leaves the sampler as it was,
  /// generators included: every chain's trees are built on scratch before
  /// any answer is given, on the donor's recorded points where they are not
  /// the store's and on the store's unless the install goes across grids. An
  /// install that fails all the same, or throws, restores every chain
  /// touched from a snapshot taken on the store's grid, after that grid is
  /// back in place.
  ///
  /// An install drops the kept draws, which belong to another fit's
  /// position, and leaves each chain's generator past its coins. \p altered,
  /// when non-null, is false on every exit but an install and then true
  /// exactly where the donor was not clean, so a forced call sets it exactly
  /// when the same call without force declines; on the store's grid the
  /// installing builds report it, never a scratch build. \p report, when
  /// non-null, is filled for a donor that went in and left as it was on every
  /// other exit: the counts are read off the installed trees
  /// (Chain::countDonorTrees), so they and altered are two readings of one
  /// install.
  WarmStartResult installForests(
      const SamplerStateData& donor,
      const std::vector<std::pair<size_t, int>>& sampleMap,
      const double* currentPredictors, bool force = true,
      bool* altered = nullptr, WarmStartReport* report = nullptr,
      const DonorLeafShape* donorLeaves = nullptr) {
    if (altered != nullptr) *altered = false;
    if (sampleMap.size() != chains_.size())
      return WarmStartResult::shapeMismatch;
    using Kind = DonorLeafShape::Kind;
    DonorLeafShape own;
    if constexpr (L::hasFunctionParams) {
      own.kind = Kind::gp;
      own.numObservations = data_.numObservations;
    } else if constexpr (L::hasVectorParams) {
      own.kind = Kind::linear;
      own.numCovariates = chains_[0]->singleLeafParams(0);
    }
    const DonorLeafShape& given = donorLeaves != nullptr ? *donorLeaves : own;
    if ((given.kind == Kind::gp) != (own.kind == Kind::gp))
      return WarmStartResult::leafKindMismatch;
    if (own.kind == Kind::gp && given.numObservations != own.numObservations)
      return WarmStartResult::leafRowsMismatch;
    bool slopesLeft = given.kind == Kind::linear &&
                      (own.kind != Kind::linear ||
                       given.numCovariates != own.numCovariates);
    bool slopesZero = given.kind == Kind::constant && own.kind == Kind::linear;
    if (!donor.missingColumns.empty() &&
        donor.missingColumns.size() != data_.numPredictors)
      return WarmStartResult::shapeMismatch;

    // a donor recorded on other points is still about these columns unless
    // the two disagree on what a column is
    bool pointsDiffer = donor.cutPoints != data_.cutPoints;
    if (pointsDiffer) {
      if (donor.cutPoints.size() != data_.numPredictors)
        return WarmStartResult::gridMismatch;
      for (size_t j = 0; j < data_.numPredictors; ++j) {
        bool categorical = data_.splitsBySubset(j);
        if (categorical != donor.cutPoints[j].empty())
          return WarmStartResult::gridMismatch;
        if (!categorical && !cutGridIsValid(donor.cutPoints[j].data(),
                                            donor.cutPoints[j].size()))
          return WarmStartResult::gridMismatch;
        // the level count fixes the grid the store keeps, and no column
        // re-quantizes its levels against another
        if (data_.types[j] == ColumnKind::orderedFactor &&
            donor.cutPoints[j] != data_.cutPoints[j])
          return WarmStartResult::gridMismatch;
      }
    }

    // the scale-leaf rule state validation holds a variance block to, every
    // leaf positive and finite, applied at BOTH install arms below for the
    // reason it applies there: a donor state is hand-buildable whichever arm
    // reaches it, and a rebuild scatters the leaf straight into a divisor, so
    // a factor that is not positive and finite is a broken model rather than
    // a bad number
    auto scaleLeavesArePositive =
      [](const std::vector<std::vector<FlatNode>>& trees) {
        for (const std::vector<FlatNode>& tree : trees)
          for (const FlatNode& node : tree)
            if (flatKindOf(node) == FlatKind::leaf &&
                !(node.value > 0.0 && std::isfinite(node.value)))
              return false;
        return true;
      };

    std::vector<ChainStateData> install(chains_.size());
    bool coefficientsKept = false;
    // by chain, a part this sampler has and the donor sample does not, which
    // starts as a new fit's once the donor is in: 1 DART, 2 a variance forest
    std::vector<int> fresh(chains_.size(), 0);
    bool dartLeft = false, varianceLeft = false;
    for (size_t c = 0; c < chains_.size(); ++c) {
      size_t dc = sampleMap[c].first;
      int slot = sampleMap[c].second;
      if (dc >= donor.chains.size()) return WarmStartResult::shapeMismatch;
      const ChainStateData& src = donor.chains[dc];
      if (src.forests.size() != chains_[c]->numForests())
        return WarmStartResult::shapeMismatch;
      bool donorDart = !src.dartProbabilities.empty();
      bool donorVariance = !src.varianceTrees.empty();
      bool ownVariance = chains_[c]->hasVarianceForest();
      fresh[c] = (chains_[c]->usesDart() && !donorDart ? 1 : 0) |
                 (ownVariance && !donorVariance ? 2 : 0);
      dartLeft = dartLeft || (!chains_[c]->usesDart() && donorDart);
      varianceLeft = varianceLeft || (!ownVariance && donorVariance);

      ChainStateData& dst = install[c];
      dst.forests.resize(src.forests.size());
      for (size_t f = 0; f < src.forests.size(); ++f) {
        const ForestStateData& sfs = src.forests[f];
        ForestStateData& dfs = dst.forests[f];
        // the donor forest's own tree count sets the saved slot's stride
        size_t nt = sfs.trees.size();
        if (nt == 0) return WarmStartResult::shapeMismatch;
        if (slot < 0) {
          dfs.trees = sfs.trees;
          dfs.treeParams = sfs.treeParams;
          dfs.treeMasks = sfs.treeMasks;
        } else {
          size_t base = static_cast<size_t>(slot) * nt;
          if (sfs.savedTrees.size() % nt != 0 ||
              sfs.savedTrees.size() < base + nt)
            return WarmStartResult::shapeMismatch;
          dfs.trees.assign(sfs.savedTrees.begin() + base,
                           sfs.savedTrees.begin() + base + nt);
          if (!sfs.savedTreeParams.empty())
            dfs.treeParams.assign(sfs.savedTreeParams.begin() + base,
                                  sfs.savedTreeParams.begin() + base + nt);
          if (!sfs.savedTreeMasks.empty())
            dfs.treeMasks.assign(sfs.savedTreeMasks.begin() + base,
                                 sfs.savedTreeMasks.begin() + base + nt);
        }
        dfs.k = sfs.k;
        // the leaf-covariate calibration is not copied: this sampler keeps
        // its own, and the donor's linear coefficients are restated in it
        // below. Each tree's block is held to the donor's own leaves; slopes
        // left behind, or never given, are 0 here
        for (size_t t = 0; t < nt && given.kind != Kind::constant; ++t) {
          size_t size = given.kind == Kind::gp
            ? given.numObservations
            : ((dfs.trees[t].size() + 1) / 2) * given.numCovariates;
          if (dfs.treeParams.size() != nt || dfs.treeParams[t].size() != size)
            return WarmStartResult::leafBlockMismatch;
        }
        if (own.kind == Kind::constant) {
          dfs.treeParams.clear();
        } else if (slopesLeft || slopesZero) {
          dfs.treeParams.resize(nt);
          for (size_t t = 0; t < nt; ++t)
            dfs.treeParams[t].assign(
              ((dfs.trees[t].size() + 1) / 2) * own.numCovariates, 0.0);
        }
      }
      dst.sigma = src.sigma;
      dst.fitMin = src.fitMin;
      dst.fitMax = src.fitMax;
      if (chains_[c]->usesDart()) {
        // one weight a predictor: installDrawnScalars copies the block whole
        if (donorDart && src.dartProbabilities.size() != data_.numPredictors)
          return WarmStartResult::shapeMismatch;
        dst.dartProbabilities = src.dartProbabilities;
        dst.dartAlpha = src.dartAlpha;
        dst.dartNumUpdatesSkipped = src.dartNumUpdatesSkipped;
      }
      dst.hasAmplitudes = src.hasAmplitudes;
      dst.amplitudeWidths = src.amplitudeWidths;
      dst.amplitudes = src.amplitudes;
      dst.amplitudeVariances = src.amplitudeVariances;
      // The scale surface travels with the mean forest it was drawn beside: a
      // slot-sourced install takes that slot of the donor's saved variance
      // buffer, index-aligned with the mean slice above because one slot base
      // drives both buffers and both index by the sample number. A live-sourced
      // install takes the donor's live pair.
      if (ownVariance && slot < 0) {
        dst.varianceTrees = src.varianceTrees;
        dst.varianceTreeMasks = src.varianceTreeMasks;
        if (!scaleLeavesArePositive(dst.varianceTrees))
          return WarmStartResult::varianceMismatch;
      } else if (ownVariance && donorVariance) {
        // The saved buffer's STRIDE is the donor's own variance tree count, so
        // a size-only bound would let a state whose live block is shorter than
        // the stride slice ACROSS slot boundaries and install a mixture of two
        // sweeps' trees - which the count check downstream cannot see. Require
        // the buffer to be exactly capacity x stride, capacity being the mean
        // side's, the quantity that bounds slot.
        size_t nvt = src.varianceTrees.size();
        size_t capacity =
          src.forests[0].savedTrees.size() / src.forests[0].trees.size();
        size_t base = static_cast<size_t>(slot) * nvt;
        if (src.savedVarianceTrees.size() != capacity * nvt ||
            base + nvt > src.savedVarianceTrees.size())
          return WarmStartResult::varianceSlotMismatch;
        dst.varianceTrees.assign(src.savedVarianceTrees.begin() + base,
                                 src.savedVarianceTrees.begin() + base + nvt);
        if (!src.savedVarianceTreeMasks.empty()) {
          if (src.savedVarianceTreeMasks.size() != src.savedVarianceTrees.size())
            return WarmStartResult::varianceSlotMismatch;
          dst.varianceTreeMasks.assign(
            src.savedVarianceTreeMasks.begin() + base,
            src.savedVarianceTreeMasks.begin() + base + nvt);
        }
        if (!scaleLeavesArePositive(dst.varianceTrees))
          return WarmStartResult::varianceMismatch;
      }
      // the donor's linear coefficients in this sampler's standardization;
      // a record that does not fit the leaf is the donor's shape
      if (!slopesLeft && !slopesZero &&
          !chains_[c]->convertDonorStandardization(dst, src,
                                                   coefficientsKept))
        return WarmStartResult::shapeMismatch;
    }

    // the directions the donor's trees were grown without, each chain's from
    // its own generator. They are drawn on the scratch, before the builds
    // below, which then judge the trees as they will be installed, and before
    // the install's collapse of what the rows leave empty. A tree the draw
    // cannot walk is left for those builds to refuse. Every exit but the
    // last puts the generators back.
    struct GeneratorGuard {
      std::vector<std::unique_ptr<Chain<L, ResidT>>>& chains;
      std::vector<std::vector<unsigned char>> bytes;
      ~GeneratorGuard() {
        for (size_t c = 0; c < bytes.size(); ++c)
          if (!bytes[c].empty())
            ext_rng_readSerializedState(chains[c]->rng(), bytes[c].data());
      }
    } generators{chains_, {}};
    if (!donor.missingColumns.empty()) {
      std::vector<std::uint8_t> raised(data_.numPredictors);
      bool draws = false;
      for (size_t j = 0; j < data_.numPredictors; ++j) {
        raised[j] =
          data_.hasMissing[j] != 0 && donor.missingColumns[j] == 0 ? 1 : 0;
        draws = draws || raised[j] != 0;
      }
      for (size_t c = 0; c < chains_.size() && draws; ++c)
        draws = chains_[c]->flatTreesAreWalkable(install[c]);
      if (draws) {
        generators.bytes.resize(chains_.size());
        for (size_t c = 0; c < chains_.size(); ++c) {
          ext_rng* rng = chains_[c]->rng();
          generators.bytes[c].resize(ext_rng_getSerializedStateLength(rng));
          if (!generators.bytes[c].empty())
            ext_rng_writeSerializedState(rng, generators.bytes[c].data());
          chains_[c]->drawDonorMissingDirections(install[c], raised.data());
        }
      }
    }

    bool fitsCold = false;
    if constexpr (L::hasFunctionParams) {
      if (pointsDiffer) {
        for (size_t c = 0; c < chains_.size(); ++c)
          chains_[c]->coldStartFunctionFits(install[c]);
        fitsCold = true;
      }
    }

    // the donor's own tree counts, which the report gives as offered. A
    // forest of this sampler holding as many or more takes the donor's trees
    // in its order and single leaves after them; one holding fewer is not
    // clean, and chooses among them below
    std::vector<std::vector<size_t>> offered(chains_.size());
    std::vector<size_t> varianceOffered(chains_.size(), 0);
    bool fewer = false;
    for (size_t c = 0; c < chains_.size(); ++c) {
      for (size_t f = 0; f < install[c].forests.size(); ++f) {
        offered[c].push_back(install[c].forests[f].trees.size());
        fewer = fewer || offered[c][f] > chains_[c]->numTreesInForest(f);
        chains_[c]->appendSingleLeaves(install[c].forests[f], f);
      }
      if (chains_[c]->hasVarianceForest()) {
        varianceOffered[c] = install[c].varianceTrees.size();
        fewer = fewer || varianceOffered[c] > chains_[c]->numVarianceTrees();
        chains_[c]->appendSingleVarianceLeaves(install[c]);
      }
    }

    // every chain's trees are built before any answer is given, with and
    // without force, so a donor that does not build is refused and never
    // declined: on the points it records where those are not the store's,
    // swapped in for the builds, and then on the store's, with the verdict
    // where the call is not forced. Across grids the trees resolve on the
    // recorded points alone and no verdict is asked: every such donor moves
    bool acrossGrids = false;
    for (size_t c = 0; c < chains_.size() && pointsDiffer; ++c)
      acrossGrids = acrossGrids || !chains_[c]->liveThresholdsAreOnGrid(
                                     install[c], data_.cutPoints);
    bool builds = true, varianceFailed = false, needsRepair = false;
    if (pointsDiffer) {
      ScopedCutGrid recorded(data_, donor.cutPoints);
      for (size_t c = 0; c < chains_.size() && builds; ++c)
        builds = chains_[c]->donorTreesStage(install[c], nullptr,
                                             varianceFailed);
    }
    // under blocks a forest the clean placement places goes in as placed
    // (Chain::placeDonorTrees); any other, and every one across grids, is
    // not clean and is chosen by rank. taken is for the report
    std::vector<std::vector<char>> byRank(chains_.size());
    std::vector<std::vector<size_t>> taken(chains_.size());
    bool unplaced = false;
    for (size_t c = 0; c < chains_.size() && builds; ++c) {
      size_t numForests = install[c].forests.size();
      std::vector<std::vector<size_t>> placed(numForests);
      bool anyPlaced = false;
      byRank[c].assign(numForests, 0);
      for (size_t f = 0; f < numForests; ++f) {
        taken[c].push_back(std::min(offered[c][f], numTreesInForest(f)));
        if (!chains_[c]->hasBlocks(f)) continue;
        byRank[c][f] = acrossGrids || !chains_[c]->placeDonorTrees(
                                        f, install[c].forests[f],
                                        offered[c][f], placed[f]);
        if (byRank[c][f]) placed[f].clear();
        unplaced = unplaced || byRank[c][f];
        anyPlaced = anyPlaced || !byRank[c][f];
      }
      if (anyPlaced)
        install[c] = chains_[c]->donorTreesOf(install[c], placed, {});
    }
    for (size_t c = 0; c < chains_.size() && builds && !acrossGrids; ++c)
      builds = chains_[c]->donorTreesStage(
        install[c], force ? nullptr : &needsRepair, varianceFailed);
    if (!builds)
      return varianceFailed ? WarmStartResult::varianceMismatch
                            : WarmStartResult::rebuildFailed;
    bool notClean = acrossGrids || needsRepair || fitsCold || slopesLeft ||
                    coefficientsKept || fewer || unplaced;
    if (notClean && !force) return WarmStartResult::declined;

    // Every check above ran before this point, so the installs below are not
    // expected to fail. They are still undoable: each chain's live state is
    // snapshotted first, on the store's grid and with the store excluded
    // since nothing here writes it, and a rebuild that fails anyway, or an
    // exception (an allocation), restores every chain touched - from the
    // flattened trees, as a save and load does, so its fits match the
    // snapshot to rounding - before the refusal or the rethrow. A warm start
    // installs no model value, so the state is all there is to put back.
    // Across grids every chain counts as touched from the start, and the
    // restore runs only once installAcrossGrids has put the store's grid
    // back: the snapshots' thresholds resolve on that grid alone.
    std::vector<ChainStateData> snapshot(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->getState(snapshot[c], false);
    size_t touched = acrossGrids ? chains_.size() : 0;
    auto restoreTouched = [&]() {
      for (size_t c = 0; c < touched; ++c) chains_[c]->setState(snapshot[c]);
    };

    // one install of every chain's \p states on the arm the donor takes,
    // false after a failure, which \p failure names
    bool installAltered = false;
    WarmStartResult failure = WarmStartResult::ok;
    auto installAll = [&](const std::vector<ChainStateData>& states) {
      auto installChains = [&]() {
        for (size_t c = 0; c < chains_.size() && failure ==
               WarmStartResult::ok; ++c) {
          touched = std::max(touched, c + 1);
          if (!chains_[c]->installForest(states[c], &installAltered))
            failure = WarmStartResult::rebuildFailed;
        }
        // the variance half, separately so a refusal names the variance
        // forest
        for (size_t c = 0; c < chains_.size() && failure ==
               WarmStartResult::ok; ++c) {
          if (chains_[c]->hasVarianceForest() &&
              !chains_[c]->installVarianceForest(states[c].varianceTrees,
                                                 states[c].varianceTreeMasks,
                                                 &installAltered))
            failure = WarmStartResult::varianceMismatch;
        }
        return failure == WarmStartResult::ok;
      };
      if (acrossGrids)
        installAcrossGrids(donor.cutPoints, currentPredictors, installChains,
                           MonotoneRepair::collapse);
      else
        installChains();
      return failure == WarmStartResult::ok;
    };
    std::vector<ChainStateData> chosen;
    try {
      if (!fewer && !unplaced)
        installAll(install);
      else if (chooseDonorTrees(install, offered, varianceOffered, byRank,
                                taken, installAll, chosen))
        installAll(chosen);
    } catch (...) {
      restoreTouched();
      throw;
    }
    if (failure != WarmStartResult::ok) {
      restoreTouched();
      return failure;
    }
    const std::vector<ChainStateData>& installed =
      fewer || unplaced ? chosen : install;
    for (size_t c = 0; c < chains_.size(); ++c) {
      if ((fresh[c] & 1) != 0) chains_[c]->startDartFresh();
      if ((fresh[c] & 2) != 0) chains_[c]->startVarianceForestFresh();
    }
    // the store itself is left alone, so the draws it holds belong to the
    // donor's fit, not this one's: drop them rather than let a read replay
    // another sampler's posterior
    currentSampleNum_ = 0;
    recordedDraws_ = 0;
    generators.bytes.clear();
    if (altered != nullptr) *altered = notClean || installAltered;
    if (report != nullptr) {
      report->mean.assign(chains_.size(), WarmStartTreeCounts());
      report->variance.assign(chains_.size(), WarmStartTreeCounts());
      for (size_t c = 0; c < chains_.size(); ++c)
        chains_[c]->countDonorTrees(installed[c], offered[c], taken[c],
                                    varianceOffered[c], report->mean[c],
                                    report->variance[c]);
      report->fitsCold = fitsCold;
      report->slopesLeftBehind = slopesLeft;
      report->coefficientsKept = coefficientsKept;
      report->slopesStartedAtZero = slopesZero;
      report->dartStartedFresh = report->varianceStartedFresh = false;
      for (int parts : fresh) {
        report->dartStartedFresh = report->dartStartedFresh || (parts & 1);
        report->varianceStartedFresh =
          report->varianceStartedFresh || (parts & 2);
      }
      report->dartLeftBehind = dartLeft;
      report->varianceLeftBehind = varianceLeft;
    }
    return WarmStartResult::ok;
  }

  // Between-sample mutation, fanned out to every chain; new-vector lifetimes
  // are the caller's problem. A re-anchor restates the retained draws in the
  // transform it leaves, so a replay returns what it did before.
  // Whether a re-derivation has too little to read is decided here, once,
  // from chain 0's response model against the incoming vector: rows are in
  // by membership, which every chain shares. Where it has, every chain takes
  // the swap with updateScale false, so none of them branches on the keep.
  ResponseUpdateResult setOffset(const double* offset, bool updateScale) {
    bool keep = updateScale && chains_[0]->lacksDistinctRowsIn(
                                 IncomingVector::offset, offset);
    updateScale = updateScale && !keep;
    SavedDrawSlots kept = savedDrawSlots();
    for (auto& chain : chains_) chain->setOffset(offset, updateScale, kept);
    if (updateScale) recordAnchor();
    return keep ? ResponseUpdateResult::scaleKept
                : ResponseUpdateResult::accepted;
  }
  ResponseUpdateResult setResponse(const double* y, bool updateScale) {
    bool keep = updateScale && chains_[0]->lacksDistinctRowsIn(
                                 IncomingVector::response, y);
    updateScale = updateScale && !keep;
    SavedDrawSlots kept = savedDrawSlots();
    for (auto& chain : chains_) chain->setResponse(y, updateScale, kept);
    if (updateScale) recordAnchor();
    return keep ? ResponseUpdateResult::scaleKept
                : ResponseUpdateResult::accepted;
  }

  /// Re-derives the transform from the response and the offset in force over
  /// the rows now in the likelihood: setOffset at updateScale true handed the
  /// offset every chain already holds, which is what a caller that has just
  /// installed weights or a mask asks for. Nothing is drawn and no latent
  /// moves. No chain is called where the family carries no units, which has
  /// nothing to derive, or where the rows in have too little to read
  /// (scaleKept): the transform in force stays, whatever rows it was read
  /// from.
  ResponseUpdateResult rederiveScale() {
    double min, max;
    chains_[0]->getScale(min, max);
    if (!chains_[0]->carriesUnits(min, max))
      return ResponseUpdateResult::accepted;
    const double* offset = chains_[0]->offset();
    if (chains_[0]->lacksDistinctRowsIn(IncomingVector::offset, offset))
      return ResponseUpdateResult::scaleKept;
    SavedDrawSlots kept = savedDrawSlots();
    for (auto& chain : chains_) chain->setOffset(offset, true, kept);
    recordAnchor();
    return ResponseUpdateResult::accepted;
  }

  /// Makes (min, max) the sampler's transform - the units its chains hold
  /// their numbers in, and the range the leaf prior's k is relative to - as
  /// a sampler re-created from its host's record takes the one recorded.
  /// Every chain is moved there now (Chain::moveScale), which requires trees
  /// that are still creation's: no install moves a chain, so a state that
  /// follows goes into the recorded transform only because this put the
  /// chains there. A chain already there is left untouched. Ignored where
  /// the family carries no such units.
  void setAnchor(double min, double max) {
    if (!chains_[0]->carriesUnits(min, max)) return;
    anchorMin_ = min;
    anchorMax_ = max;
    for (auto& chain : chains_) {
      double chainMin, chainMax;
      chain->getScale(chainMin, chainMax);
      if (chainMin != min || chainMax != max) chain->moveScale(min, max);
    }
  }
  void getAnchor(double& min, double& max) const {
    min = anchorMin_;
    max = anchorMax_;
  }
  /// Case weights: a pointer swap entering the next iteration's node
  /// statistics and sigma draw with nothing rescaled, plus, for logistic -
  /// whose weights are the Polya-Gamma counts - an immediate latent refresh
  /// off each chain's OWN generator, so the draws stay independent of the
  /// thread count. The host rejects the families that carry no weights.
  void setWeights(const double* weights) {
    for (auto& chain : chains_) chain->setWeights(weights);
  }
  /// A digest of the weights in force. Weights are chain-invariant - setWeights
  /// fans the same pointer to every chain - so chain 0 answers for all.
  std::uint64_t weightsDigest() const { return chains_[0]->weightsDigest(); }
  /// Re-derive every chain's weight-dependent latents against the weights
  /// already in force, each off its OWN generator, so the draws stay
  /// independent of the thread count exactly as setWeights's do.
  /// storedZeroWeightRows is the installed state's zeroWeightRows record, or
  /// null when it carries none; every chain reads the same one.
  void reapplyWeights(const unsigned char* storedZeroWeightRows = nullptr) {
    for (auto& chain : chains_) chain->reapplyWeights(storedZeroWeightRows);
  }
  /// The rows at a weight that is not positive, a byte each, for a family
  /// whose state names them; false for every other. Weights are
  /// chain-invariant, so chain 0 answers for all.
  bool zeroWeightRows(unsigned char* flags) const {
    return chains_[0]->zeroWeightRows(flags);
  }
  /// A new per-observation censoring status, fanned to every chain: each
  /// rebuilds its own censoring structure against its own observed times.
  /// Nothing is drawn here; the setResponse that follows redraws each chain's
  /// censored latents off its OWN generator, the setWeights rule.
  void setSurvivalStatus(const double* status) {
    for (auto& chain : chains_) chain->setSurvivalStatus(status);
  }
  /// A digest of the censoring structure in force. The status is
  /// chain-invariant - setSurvivalStatus fans the same vector to every chain -
  /// so chain 0 answers for all, as the weights digest does.
  std::uint64_t survivalDigest() const { return chains_[0]->survivalDigest(); }
  /// Re-derive every chain's censored latents against the structure already in
  /// force, each off its OWN generator, exactly as reapplyWeights does.
  void reapplySurvivalStatus() {
    for (auto& chain : chains_) chain->reapplySurvivalStatus();
  }
  void setSigma(double sigmaOriginalScale) {
    for (auto& chain : chains_) chain->setSigma(sigmaOriginalScale);
  }
  const double* latents(size_t chainNum = 0) const {
    return chains_[chainNum]->latents();
  }

  /// Between-run reconfiguration. Chain count, tree count, generators, and the
  /// cut grid are fixed at creation;
  /// the host refuses changes to those.
  void setNumThreads(size_t numThreads) {
    options_.numThreads = numThreads;
    for (auto& chain : chains_) chain->setNumThreads(numThreads);
  }
  void setNumThin(size_t numThin) {
    options_.numThin = numThin;
    for (auto& chain : chains_) chain->setNumThin(numThin);
  }
  void setVerbose(bool verbose, size_t printEvery) {
    options_.verbose = verbose;
    options_.printEvery = printEvery;
    for (auto& chain : chains_) chain->setVerbose(verbose, printEvery);
  }
  /// The multiplier taking internal-scale fits to the original response
  /// scale: the response range for gaussian, 1 for the binary families.
  double fitScale() const { return chains_[0]->fitScale(); }

  /// Reconfigure saved-tree storage: toggling keepTrees or changing the
  /// capacity reallocates every chain's slots and resets the write position;
  /// a no-op when nothing changes, preserving stored samples.
  ///
  /// Every chain's new store is built aside before any is swapped in, so an
  /// allocation that fails (a capacity too large to hold) throws with the
  /// store, its capacity and its recorded draws exactly as they were. The old
  /// and new stores coexist for the length of the swap.
  void setTreeStorage(bool keepTrees, size_t numSamplesToStore) {
    size_t capacity =
      keepTrees ? (numSamplesToStore > 0 ? numSamplesToStore : 1) : 0;
    if (keepTrees == options_.keepTrees && capacity == savedTreeCapacity())
      return;
    resizeSavedTrees(capacity);
    options_.keepTrees = keepTrees;
    options_.numSamplesToStore = numSamplesToStore;
  }

  /// Builds every chain's store at capacity, then swaps them all in and
  /// empties the store's draw count; throws, changing nothing, if any build
  /// fails.
  void resizeSavedTrees(size_t capacity) {
    using Store = typename Chain<L, ResidT>::SavedTreeStore;
    std::vector<Store> stores;
    stores.reserve(chains_.size());
    for (auto& chain : chains_)
      stores.push_back(chain->prepareSavedTrees(capacity));
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->installSavedTrees(stores[c]);
    currentSampleNum_ = 0;
    recordedDraws_ = 0;
  }

  /// Install a replacement prior on every chain; see ModelParameters.
  void setModel(const ModelParameters& model) {
    for (auto& chain : chains_) chain->setModel(model);
  }

  double sumOfSquaredResiduals(size_t chainNum) {
    return chains_[chainNum]->sumOfSquaredResiduals();
  }

  /// Prior draws, every chain from its own generator: tree structures alone
  /// (fits left stale), or leaf parameters with the fits rebuilt to match.
  void sampleTreesFromPrior() {
    for (auto& chain : chains_) chain->sampleTreesFromPrior();
  }
  void sampleNodeParametersFromPrior() {
    for (auto& chain : chains_) chain->sampleNodeParametersFromPrior();
  }
  /// The variance forest's own prior draw, structures and leaf factors; a
  /// no-op on every homoscedastic sampler.
  void sampleVarianceForestFromPrior() {
    for (auto& chain : chains_) chain->sampleVarianceForestFromPrior();
  }

  /// Warm-start initializer: run numSweeps of grow-from-root in place on every
  /// chain, fanning across chains on up to min(numThreads, numChains) workers
  /// exactly as run() does, so the grown forest is thread-count-independent
  /// (each chain draws only on its own generator). A single chain runs inline
  /// on the caller's thread; multi-chain workers run on spawned threads.
  /// Every chain, inline or threaded, draws from its own Mersenne Twister,
  /// seeded from R's stream once at sampler creation, and never touches R's
  /// stream during sampling. Constant leaf only (Chain no-ops for
  /// vector/function leaves; the R surface refuses them).
  void growFromRoot(size_t numSweeps) {
    size_t numChains = chains_.size();
    size_t numWorkers = options_.numThreads < numChains ? options_.numThreads
                                                        : numChains;
    if (numWorkers <= 1) {
      for (auto& chain : chains_) chain->growForestFromRoot(numSweeps);
      return;
    }

    // spawned with SIGINT blocked, as run() spawns. An exception escaping a
    // thread body is std::terminate, so each worker keeps the first one a
    // chain throws and it is rethrown here after the join; a failed spawn
    // lets the started workers finish their chains.
    std::mutex errorMutex;
    std::exception_ptr firstError;
    std::vector<std::thread> workers;
    auto worker = [this, numWorkers, numChains, numSweeps, &errorMutex,
                   &firstError](size_t w) {
      for (size_t c = w; c < numChains; c += numWorkers) {
        try {
          chains_[c]->growForestFromRoot(numSweeps);
        } catch (...) {
          std::lock_guard<std::mutex> lock(errorMutex);
          if (!firstError) firstError = std::current_exception();
        }
      }
    };
    spawnWorkers(workers, numWorkers, worker, []() {});
    for (std::thread& worker : workers) worker.join();
    if (firstError) std::rethrow_exception(firstError);
  }

  /// Info dump of forest forestIndex; the per-node output format is R-visible
  /// and pinned by tests. Without keepTrees, or with useLiveTrees, the live
  /// trees print and sample indices are ignored; otherwise sampleIndices are
  /// DRAW numbers on the same oldest-first axis predict and getTrees report,
  /// and the draws they name print in the saved-tree format. treeIndices are
  /// read against the NAMED forest's tree count, which a multi-forest sampler
  /// states per forest; range-checking is the bridge's, as everywhere here.
  void printTrees(const size_t* chainIndices, size_t numChainIndices,
                  const size_t* sampleIndices, size_t numSampleIndices,
                  const size_t* treeIndices, size_t numTreeIndices,
                  size_t forestIndex, bool useLiveTrees) {
    int indent = 0;
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      if (numChainIndices > 1) {
        ext_printf("chain %lu\n", static_cast<unsigned long>(chainNum + 1));
        indent += 2;
      }
      if (useLiveTrees || !options_.keepTrees) {
        for (size_t k = 0; k < numTreeIndices; ++k)
          chains_[chainNum]->printTree(treeIndices[k], indent, forestIndex);
      } else {
        for (size_t j = 0; j < numSampleIndices; ++j) {
          size_t sampleNum = sampleIndices[j];
          if (numSampleIndices > 1) {
            ext_printf("%*ssample %lu\n", indent, "",
                       static_cast<unsigned long>(sampleNum + 1));
            indent += 2;
          }
          size_t slot = savedSlotForDraw(sampleNum);
          for (size_t k = 0; k < numTreeIndices; ++k)
            chains_[chainNum]->printSavedTree(slot, treeIndices[k], indent,
                                              forestIndex);
          if (numSampleIndices > 1)
            indent -= 2;
        }
      }
      if (numChainIndices > 1)
        indent -= 2;
    }
  }

  /// Replace the entire data set: predictors, response, and optionally
  /// weights, offset, and test predictors, with a possibly different number
  /// of observations (all borrowed; the predictor count is fixed). Not
  /// transactional: cut points are rebuilt from scratch, existing splits are
  /// remapped onto the value-nearest new cuts, and any subtree left invalid
  /// or empty collapses. A column's missing-value flag carries over, and one
  /// the replacement first raises draws its rules' directions
  /// (Chain::applyNewData). Gaussian chains keep sigma and the variance prior
  /// fixed on the original scale. The retained draws stay the functions they
  /// were (Chain::applyNewData), a gp leaf's excepted.
  ///
  /// False REFUSES the replacement and the sampler is UNTOUCHED: a factor
  /// cell on one side or the other is not a level code of its column's table.
  /// Everything either build can refuse is decided here, before either is
  /// installed - the level table is fixed at creation and a replacement does
  /// not move it, so the answer does not depend on the order they run in.
  /// That matters because they cannot both be rolled back: the test store is
  /// quantized against the grid the training replacement rebuilds, so it must
  /// be built second, and a refusal there would leave the store holding new
  /// codes and a new observation count while the chains were still sized for
  /// the old ones. The two status tests below are backstops that the checks
  /// here make unreachable.
  [[nodiscard]] bool setData(const double* x, const double* y,
                             size_t numObservations, const double* weights,
                             const double* offset, const double* x_test,
                             size_t numTestObservations,
                             const double* testOffset = nullptr) {
    bool hasTest = x_test != nullptr && numTestObservations > 0;
    if (!data_.replacementLevelCodesAreValid(x, numObservations)) return false;
    if (hasTest &&
        !data_.testSourceIsIngestible(x_test, numTestObservations))
      return false;
    // recover parameters against the old fits and partitions before anything
    // moves; the old cut values drive the split remap. The variance forest's
    // factors ride the same phase: their per-observation slab is strided by
    // the store's observation count, which the replacement moves.
    std::vector<typename Chain<L, ResidT>::TreeParameters> params(chains_.size());
    std::vector<typename Chain<L, ResidT>::TreeParameters>
      varianceParams(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c) {
      chains_[c]->recoverTreeParameters(params[c]);
      chains_[c]->recoverVarianceParameters(varianceParams[c]);
    }

    std::vector<std::vector<double>> oldCutPoints(data_.cutPoints);
    std::vector<std::uint8_t> raisedMissing(data_.hasMissing);

    if (!data_.setData(x, numObservations)) return false;
    if (hasTest) {
      if (!data_.buildTest(x_test, numTestObservations)) return false;
      data_.testOffset = testOffset;
    } else {
      data_.resetTestStorage();
    }

    // the columns the replacement first gave a missing value; the flags
    // themselves carry over
    bool anyRaised = false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      raisedMissing[j] =
        data_.hasMissing[j] != 0 && raisedMissing[j] == 0 ? 1 : 0;
      anyRaised = anyRaised || raisedMissing[j] != 0;
    }

    SavedDrawSlots kept = savedDrawSlots();
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->applyNewData(y, weights, offset, oldCutPoints, params[c],
                               varianceParams[c], kept,
                               anyRaised ? raisedMissing.data() : nullptr);
    recordAnchor();
    return true;
  }

  /// Replace the predictor matrix from a borrowed view of numObservations
  /// rows, read for the call only (the store keeps its old values on
  /// failure). Unless forceUpdate, a leaf that would empty in any tree of any
  /// chain, or a monotone tree whose leaf values would be out of order, rolls
  /// the whole change back; forceUpdate instead collapses emptied leaves into
  /// their parents and sets every leaf of a monotone tree left out of order
  /// to zero. Any view shape is consumed column by column
  /// (ColumnStore::mutateColumnFromSource): a CSC column onto a CSC-backed
  /// one as its entries, never densified, and a CSC or coded column onto a
  /// dense-backed one through one reused column of scratch.
  ///
  /// updateCutPoints re-derives each numeric column's grid for its new
  /// values, as creation would; placement then says where the splits on a
  /// column whose grid changed go (runPredictorTransaction).
  PredictorUpdateResult setPredictor(
      const PredictorSource& newX, bool forceUpdate, bool updateCutPoints,
      SplitPlacement placement = SplitPlacement::byPosition) {
    WholeMatrixUpdate strategy{data_, newX};
    return runPredictorTransaction(strategy, forceUpdate, updateCutPoints,
                                   placement);
  }

  /// Overwrite a subset of columns in place; view column k, numObservations
  /// rows, fills store column columns[k]. Same transaction and view semantics
  /// as setPredictor; a column named twice is applied in order and, on a
  /// reject, restored exactly.
  PredictorUpdateResult updatePredictor(
      const PredictorSource& newColumns, const size_t* columns,
      size_t numColumns, bool forceUpdate, bool updateCutPoints,
      SplitPlacement placement = SplitPlacement::byPosition) {
    SubsetUpdate strategy{data_, newColumns, columns, numColumns};
    return runPredictorTransaction(strategy, forceUpdate, updateCutPoints,
                                   placement);
  }

  /// Dense convenience spellings: plain column-major blocks over the store's
  /// own row count.
  PredictorUpdateResult setPredictor(
      const double* newX, bool forceUpdate, bool updateCutPoints,
      SplitPlacement placement = SplitPlacement::byPosition) {
    return setPredictor(densePredictorSource(newX, data_.numObservations,
                                             data_.numPredictors),
                        forceUpdate, updateCutPoints, placement);
  }
  PredictorUpdateResult updatePredictor(
      const double* newColumns, const size_t* columns, size_t numColumns,
      bool forceUpdate, bool updateCutPoints,
      SplitPlacement placement = SplitPlacement::byPosition) {
    return updatePredictor(
      densePredictorSource(newColumns, data_.numObservations, numColumns),
      columns, numColumns, forceUpdate, updateCutPoints, placement);
  }

  /// Install externally chosen cut points (strictly increasing) for a subset
  /// of columns, unweighted unless a column is handed back the grid it holds,
  /// and unconditionally refresh the trees: each split on a named
  /// column goes where placement says (splitsToMove), and one whose interval
  /// empties or that loses its observations collapses into its parent,
  /// exactly as in a forced predictor update.
  /// currentPredictors is the call-time predictor matrix the columns
  /// re-quantize from (data@x); null for CSC/mixed stores, whose columns read
  /// their retained slices instead.
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints, const size_t* columns,
                    size_t numColumns, const double* currentPredictors,
                    SplitPlacement placement = SplitPlacement::byPosition) {
    std::vector<std::vector<double>> oldCutPoints(data_.numPredictors);
    for (size_t k = 0; k < numColumns; ++k)
      if (oldCutPoints[columns[k]].empty())
        oldCutPoints[columns[k]] = data_.cutPoints[columns[k]];
    // a grid other than the one the column holds goes in unweighted; the held
    // grid, handed back bit for bit, keeps its weights, so the call is a no-op
    for (size_t k = 0; k < numColumns; ++k) {
      size_t j = columns[k];
      bool held = numCutPoints[k] == data_.numCuts[j] &&
        std::memcmp(newCutPoints[k], data_.cutPoints[j].data(),
                    numCutPoints[k] * sizeof(double)) == 0;
      std::vector<double> mass;
      if (held) mass = data_.cutMass[j];
      data_.setCutPointsForColumn(j, newCutPoints[k], numCutPoints[k],
                                  currentPredictors,
                                  mass.empty() ? nullptr : mass.data());
    }
    const std::vector<std::vector<double>>* moved =
      splitsToMove(oldCutPoints, placement);
    for (auto& chain : chains_) chain->forceRefreshTrees(moved, placement);
  }

  /// Narrows \p oldCutPoints, the grids some columns held before a change
  /// (empty for every other column), to the columns whose splits must move
  /// onto the grid now held: by position one whose count changed, by value
  /// one whose grid changed at all. Null when no column is left, so a change
  /// that moves no split takes the path of one that replaced no grid.
  const std::vector<std::vector<double>>* splitsToMove(
      std::vector<std::vector<double>>& oldCutPoints,
      SplitPlacement placement) const {
    bool any = false;
    for (size_t j = 0; j < oldCutPoints.size(); ++j) {
      std::vector<double>& old = oldCutPoints[j];
      if (old.empty()) continue;
      const std::vector<double>& held = data_.cutPoints[j];
      bool same = old.size() == held.size() &&
        (placement == SplitPlacement::byPosition ||
         std::memcmp(old.data(), held.data(),
                     old.size() * sizeof(double)) == 0);
      if (same) old.clear();
      else any = true;
    }
    return any ? &oldCutPoints : nullptr;
  }

  /// Install one column's new values observation-by-observation in random
  /// scan order, declining exactly those whose move would empty a leaf in any
  /// tree of any forest of any chain and, where the column holds no missing
  /// value and one would leave a monotone tree's leaf values out of order,
  /// every missing one; installed must have room for a flag per observation.
  /// Returns finalize() validity, which the guard makes true by construction.
  bool updatePredictorPerObservation(const double* newColumn, size_t column,
                                     bool* installed) {
    size_t n = data_.numObservations;
    UpdateSessionImpl session(*this, newColumn, column);

    std::vector<size_t> scanOrder(n);
    ext_rng_drawPermutation(chains_[0]->rng(), scanOrder.data(), n);

    for (size_t i = 0; i < n; ++i) {
      size_t j = scanOrder[i];
      bool valid = session.observationWouldRemainValid(j);
      installed[j] = valid;
      if (valid) session.commitObservation(j);
    }
    return session.finalize();
  }

  std::unique_ptr<PredictorUpdateSession> beginPredictorUpdate(
    const double* newColumn, size_t column) {
    return std::make_unique<UpdateSessionImpl>(*this, newColumn, column);
  }

  ext_rng* rng() const { return chains_[0]->rng(); }

  /// The leaf covariate designation (linear and GP leaves); 0/null for
  /// scalar leaf models.
  size_t numLeafCovariates() const {
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      return chains_[0]->leaf().numCovariates();
    else
      return 0;
  }
  const size_t* leafCovariateColumns() const {
    if constexpr (L::hasVectorParams || L::hasFunctionParams)
      return chains_[0]->leaf().covariateColumns().data();
    else
      return nullptr;
  }
  /// True for function-valued (GP) leaf models, whose state and reporting
  /// layouts differ from the vector-parameter ones. A property of the leaf
  /// model alone, so constexpr.
  static constexpr bool usesFunctionLeaves() { return L::hasFunctionParams; }

  ResponseFamily family() const { return family_; }
  double sigma(size_t chainNum = 0) const { return chains_[chainNum]->sigma(); }
  double k(size_t chainNum = 0) const { return chains_[chainNum]->k(); }
  bool kIsSampled() const { return options_.updateK; }
  bool usesDart() const { return options_.useDart; }
  const ColumnStore& data() const { return data_; }
  const Chain<L, ResidT>& chain(size_t chainNum) const { return *chains_[chainNum]; }
  Chain<L, ResidT>& chain(size_t chainNum) { return *chains_[chainNum]; }
  size_t numChains() const { return chains_.size(); }
  size_t numThreads() const { return options_.numThreads; }

  // BCF surface, fanned to every chain; benign on
  // single-forest samplers, where numForests() is 1 and no amplitudes exist.
  size_t numForests() const { return chains_[0]->numForests(); }
  /// Whether the forest coupling permits a whole-response swap; chain 0
  /// answers for all, as every chain carries the same combiner and family.
  bool supportsResponseMutation() const {
    return chains_[0]->supportsResponseMutation();
  }
  size_t numReportedLocations() const {
    return chains_[0]->numReportedLocations();
  }
  size_t numVariableCountForests() const {
    return chains_[0]->numVariableCountForests();
  }
  size_t numOrdinalThresholds() const {
    return chains_[0]->numOrdinalThresholds();
  }
  /// Whether the response family carries a shape r; chain 0 answers for
  /// all, as every chain carries the same family.
  bool carriesShape() const { return chains_[0]->carriesShape(); }
  /// Chain chainNum's shape r in force, 0 off a family carrying one.
  double shape(size_t chainNum) const {
    return chains_[chainNum]->shape();
  }
  /// Whether the response family carries a residual df nu; chain 0 answers for
  /// all, as every chain carries the same error law.
  bool carriesResidualDf() const { return chains_[0]->carriesResidualDf(); }
  bool testFitsAreDefined() const {
    return chains_[0]->testFitsAreDefined();
  }
  /// Whether the recorded per-forest fits and glue channels carry defined
  /// values; chain 0 answers for all, as every chain carries the same combiner.
  bool forestReportingIsDefined() const {
    return chains_[0]->forestReportingIsDefined();
  }
  /// Installs forest forestIndex's n x numColumns ROW-major amplitude basis in
  /// every chain; false, installing nothing, on a refusal. Every chain refuses
  /// on the same conditions, so the fan-out cannot land half applied.
  bool setForestBasis(size_t forestIndex, const double* values,
                      size_t numColumns) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setForestBasis(forestIndex, values, numColumns) &&
                  installed;
    return installed;
  }
  /// Whether the forest coupling admits a caller-supplied per-forest weight;
  /// chain 0 answers for all, as every chain carries the same combiner.
  bool supportsForestWeights() const {
    return chains_[0]->supportsForestWeights();
  }
  /// Installs (or clears, at null) a borrowed per-observation weight on forest
  /// forestIndex in every chain, each of which composes it into its own
  /// scratch; false, installing nothing, on a refusal. Every chain refuses on
  /// the same two conditions, so the fan-out cannot land half applied.
  bool setForestWeights(size_t forestIndex, const double* weights) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setForestWeights(forestIndex, weights) && installed;
    return installed;
  }
  /// The leaf model in force; a property of the type, so every chain agrees.
  static constexpr LeafModelKind leafModel() { return leafModelKindOf<L>(); }
  /// Forest forestIndex's calibration on ONE chain: the chains carry their own
  /// transforms and their own drawn k, so they can and do disagree, and the
  /// caller reads each rather than a flattened summary.
  /// Chain::forestCalibration states the semantics.
  ForestCalibration forestCalibration(size_t chainNum,
                                      size_t forestIndex) const {
    return chains_[chainNum]->forestCalibration(forestIndex);
  }
  /// Names forest forestIndex's leaf-prior sd on EVERY chain; false, writing
  /// nothing, on a refusal. Every chain refuses on the same two conditions, so
  /// the fan-out cannot land half applied. Chain::setForestNamedSd states the
  /// semantics.
  bool setForestNamedSd(size_t forestIndex, double namedSd) {
    bool written = true;
    for (auto& chain : chains_)
      written = chain->setForestNamedSd(forestIndex, namedSd) && written;
    return written;
  }
  /// Restates forest forestIndex's fixed k on EVERY chain; false, writing
  /// nothing, on a refusal, which every chain makes on the same predicate.
  /// Chain::setForestFixedK states the semantics.
  bool setForestFixedK(size_t forestIndex, double k) {
    bool written = true;
    for (auto& chain : chains_)
      written = chain->setForestFixedK(forestIndex, k) && written;
    return written;
  }
  /// Restates map forest forestIndex's spread on EVERY chain; false, writing
  /// nothing, on a refusal, which every chain makes on the same predicate.
  /// Chain::setForestMapSd states the semantics.
  bool setForestMapSd(size_t forestIndex, double sd) {
    bool written = true;
    for (auto& chain : chains_)
      written = chain->setForestMapSd(forestIndex, sd) && written;
    return written;
  }
  /// Whether the response family implements the active-row channel; chain 0
  /// answers for all, as every chain carries the same family.
  bool supportsActiveRows() const { return chains_[0]->supportsActiveRows(); }
  /// Installs (or clears, at null) the 0/1 active-row mask in every chain,
  /// each of which copies it into its own family buffer; false, installing
  /// nothing, on a refusal. Every chain refuses on the same family predicate
  /// and the same value scan, so the fan-out cannot land half applied.
  bool setActiveRows(const double* active) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setActiveRows(active) && installed;
    return installed;
  }
  /// Whether the forest coupling owns a replaceable count-matrix response;
  /// chain 0 answers for all, as every chain carries the same combiner.
  bool supportsCountsMutation() const {
    return chains_[0]->supportsCountsMutation();
  }
  /// Installs the borrowed replacement counts and trials in every chain, so a
  /// multi-chain sampler cannot be left half swapped; false, installing
  /// nothing, on a refusal. Every chain refuses on the same predicate.
  bool setCounts(const int* counts, const int* trials) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setCounts(counts, trials) && installed;
    return installed;
  }
  /// Installs (or clears, at null) the borrowed n x K category offset in every
  /// chain, so a multi-chain sampler cannot be left half offset; false,
  /// installing nothing, on a refusal. Every chain refuses on the same
  /// predicate.
  bool setCategoryOffset(const double* offset) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setCategoryOffset(offset) && installed;
    return installed;
  }
  /// The same for the borrowed nTest x K TEST offset, which every chain's
  /// reported test blend reads; false, installing nothing, on a refusal.
  bool setCategoryTestOffset(const double* offset) {
    bool installed = true;
    for (auto& chain : chains_)
      installed = chain->setCategoryTestOffset(offset) && installed;
    return installed;
  }
  /// The amplitude channel; chain 0 answers the two shape questions for all,
  /// as every chain carries the same combiner and the same bases.
  size_t totalAmplitudes() const { return chains_[0]->totalAmplitudes(); }
  size_t numForestAmplitudes(size_t forestIndex) const {
    return chains_[0]->numForestAmplitudes(forestIndex);
  }
  bool amplitudes(size_t chainNum, double* out) const {
    return chains_[chainNum]->amplitudes(out);
  }
  void forestTotalFits(size_t chainNum, size_t forestIndex, double* out) const {
    chains_[chainNum]->forestTotalFits(forestIndex, out);
  }
  /// Chain chainNum's combined per-observation location on the response scale
  /// and without the offset; false, writing nothing, on a multi-location
  /// coupling. Non-const: the read refills the combiner's scratch buffer
  /// (Chain::fitsWithoutOffset states the choice).
  bool fitsWithoutOffset(size_t chainNum, double* out) {
    return chains_[chainNum]->fitsWithoutOffset(out);
  }
  void forestVariableCounts(size_t chainNum, size_t forestIndex,
                            std::uint32_t* out) const {
    chains_[chainNum]->forestVariableCounts(forestIndex, out);
  }
  size_t numTrees() const { return options_.numTrees; }
  size_t numTreesInForest(size_t forestIndex) const {
    return chains_[0]->numTreesInForest(forestIndex);
  }
  size_t numObservations() const { return data_.numObservations; }
  size_t numPredictors() const { return data_.numPredictors; }
  size_t numTestObservations() const { return data_.numTestObservations; }

private:
  friend struct TestPeer;

  /// The store's cut grid, its weights and everything quantized against it,
  /// copied so that a grid installed for a while can be taken out again byte
  /// for byte. Rank columns re-quantize into their own storage and not into
  /// codes, and a rank-backed test column into test.sparseColumns (empty,
  /// hence a no-op restore, on every dense-test path), so both are held.
  struct GridSnapshot {
    std::vector<std::vector<double>> cutPoints, cutMass;
    std::vector<std::uint32_t> numCuts;
    std::vector<xint_t> codes, testCodes;
    std::vector<SparseColumnData> sparseColumns, testSparseColumns;
    explicit GridSnapshot(const ColumnStore& data)
      : cutPoints(data.cutPoints), cutMass(data.cutMass),
        numCuts(data.numCuts), codes(data.train.codes),
        testCodes(data.test.codes), sparseColumns(data.train.sparseColumns),
        testSparseColumns(data.test.sparseColumns) {}
    /// Moves the copy back; the snapshot is spent.
    void restore(ColumnStore& data) {
      data.cutPoints = std::move(cutPoints);
      data.numCuts = std::move(numCuts);
      data.cutMass = std::move(cutMass);
      data.train.codes = std::move(codes);
      data.test.codes = std::move(testCodes);
      data.train.sparseColumns = std::move(sparseColumns);
      data.test.sparseColumns = std::move(testSparseColumns);
    }
  };

  /// A warm start's choice by rank (installForests): \p install holds each
  /// chain's donor sample, every forest filled to at least this sampler's
  /// count, and \p offered and \p varianceOffered the donor's counts. A
  /// forest given more trees than it holds, or one \p byRank names, tries
  /// every tree in the slots of each of its blocks() blocks, a block's count
  /// at a time in the donor's order, by \p installAll, the install the donor
  /// takes, single leaves after the last; without blocks, and in the
  /// variance forest, the slots are one block. After each round each tree
  /// tried is ranked by the variance over the rows of its fits as installed
  /// in that block (Chain::liveTreeVariation), of the log of its factors in
  /// the variance forest. Each block in turn keeps up to its count of the
  /// trees that vary most there that no earlier block kept, a tie to the
  /// earlier, in the donor's order, then single leaves; under \p byRank a
  /// tree that block holds as a single leaf and the donor as more is not
  /// kept there. Each chain chooses from its own sample; \p chosen receives
  /// the trees and \p taken, by forest, how many of the donor's each holds.
  /// Any other forest goes through the rounds as given. False where a
  /// round's install fails.
  template <typename Install>
  bool chooseDonorTrees(const std::vector<ChainStateData>& install,
                        const std::vector<std::vector<size_t>>& offered,
                        const std::vector<size_t>& varianceOffered,
                        const std::vector<std::vector<char>>& byRank,
                        std::vector<std::vector<size_t>>& taken,
                        Install&& installAll,
                        std::vector<ChainStateData>& chosen) {
    size_t numChains = chains_.size(), numForests = chains_[0]->numForests();
    constexpr size_t leaf = Chain<L, ResidT>::singleLeaf;
    // forest f of a chain, the variance forest at f numForests: its slots by
    // block and each tried tree's variation in each block, by donor index,
    // below 0 where it is not kept there; no blocks for a forest that takes
    // all it is given
    struct Choice {
      std::vector<std::vector<size_t>> blocks;
      std::vector<std::vector<double>> variation;
      size_t slots = 0;
    };
    std::vector<std::vector<Choice>> choices(numChains);
    size_t numRounds = 1;
    for (size_t c = 0; c < numChains; ++c) {
      choices[c].resize(numForests + 1);
      for (size_t f = 0; f <= numForests; ++f) {
        Choice& choice(choices[c][f]);
        size_t given = f < numForests ? offered[c][f] : varianceOffered[c];
        choice.slots = f < numForests ? chains_[c]->numTreesInForest(f)
                                      : chains_[c]->numVarianceTrees();
        if (given <= choice.slots && !(f < numForests && byRank[c][f]))
          continue;
        if (f < numForests) {
          choice.blocks = chains_[c]->slotsByBlock(f);
        } else {
          choice.blocks.resize(1);
          for (size_t t = 0; t < choice.slots; ++t)
            choice.blocks[0].push_back(t);
        }
        choice.variation.assign(choice.blocks.size(),
                                std::vector<double>(given, 0.0));
        for (const std::vector<size_t>& block : choice.blocks)
          numRounds =
            std::max(numRounds, (given + block.size() - 1) / block.size());
      }
    }
    // the trees of round r by slot, or of the choice where r is numRounds; a
    // block with no round r left tries its first trees again, and is not read
    auto treesOf = [&](size_t c, size_t r) {
      std::vector<std::vector<size_t>> mean(numForests);
      std::vector<size_t> variance;
      for (size_t f = 0; f <= numForests; ++f) {
        const Choice& choice(choices[c][f]);
        if (choice.blocks.empty()) continue;
        std::vector<size_t>& trees(f < numForests ? mean[f] : variance);
        trees.assign(choice.slots, leaf);
        size_t given = choice.variation[0].size();
        std::vector<bool> kept(given, false);
        if (r == numRounds && f < numForests) taken[c][f] = 0;
        bool byBlock = f < numForests && byRank[c][f];
        for (size_t b = 0; b < choice.blocks.size(); ++b) {
          const std::vector<size_t>& block(choice.blocks[b]);
          std::vector<size_t> picked;
          if (r < numRounds) {
            size_t first = r * block.size() < given ? r * block.size() : 0;
            for (size_t i = first; i < given && i < first + block.size(); ++i)
              picked.push_back(i);
          } else {
            const std::vector<double>& ranked(choice.variation[b]);
            for (size_t i = 0; i < given; ++i)
              if (!kept[i] && (!byBlock || ranked[i] >= 0.0))
                picked.push_back(i);
            std::stable_sort(picked.begin(), picked.end(),
                             [&](size_t i, size_t j) {
                               return ranked[i] > ranked[j];
                             });
            picked.resize(std::min(picked.size(), block.size()));
            std::sort(picked.begin(), picked.end());
            for (size_t i : picked) kept[i] = true;
            if (f < numForests) taken[c][f] += picked.size();
          }
          for (size_t k = 0; k < picked.size(); ++k)
            trees[block[k]] = picked[k];
        }
      }
      return chains_[c]->donorTreesOf(install[c], mean, variance);
    };
    std::vector<ChainStateData> round(numChains);
    std::vector<double> ranks;
    for (size_t r = 0; r < numRounds; ++r) {
      for (size_t c = 0; c < numChains; ++c) round[c] = treesOf(c, r);
      if (!installAll(round)) return false;
      for (size_t c = 0; c < numChains; ++c) {
        for (size_t f = 0; f <= numForests; ++f) {
          Choice& choice(choices[c][f]);
          if (choice.blocks.empty()) continue;
          if (f < numForests)
            chains_[c]->liveTreeVariation(f, ranks);
          else
            chains_[c]->liveVarianceTreeVariation(ranks);
          for (size_t b = 0; b < choice.blocks.size(); ++b) {
            const std::vector<size_t>& block(choice.blocks[b]);
            std::vector<double>& ranked(choice.variation[b]);
            size_t first = r * block.size();
            for (size_t k = 0; first + k < ranked.size() && k < block.size();
                 ++k)
              ranked[first + k] =
                f < numForests && byRank[c][f] &&
                    chains_[c]->liveTreeIsLeaf(f, block[k]) &&
                    install[c].forests[f].trees[first + k].size() > 1
                  ? -1.0
                  : ranks[block[k]];
          }
        }
      }
    }
    chosen.resize(numChains);
    for (size_t c = 0; c < numChains; ++c) chosen[c] = treesOf(c, numRounds);
    return true;
  }

  /// Installs trees stored on \p storedCutPoints, a grid other than the
  /// store's, and moves them onto the store's, which with its weights and
  /// codes is the same before and after, byte for byte. \p storedCutPoints
  /// holds one entry per predictor, empty exactly where the column splits by
  /// subset and a valid grid elsewhere.
  ///
  /// The stored grid goes in on every column where the two differ,
  /// unweighted, its codes re-quantized from \p currentPredictors (null for
  /// CSC/mixed stores); \p install, called as install() and returning
  /// whether it succeeded, then builds and partitions every chain's live
  /// trees there, so each stored leaf holds the rows the stored tree sends
  /// it; the store's own grid is put back; and each chain moves its splits
  /// on the differing columns by value and re-routes
  /// (Chain::forceRefreshTrees), a merged leaf weighing the leaves beneath it
  /// by the rows install gave them. The trees \p install leaves must be
  /// whole trees with fits on the stored grid, as Chain::setState leaves
  /// them. A false from install puts the store's grid back, moves nothing
  /// and is returned. A throw, from install or from the re-quantization
  /// before it, leaves with the store's grid, weights and codes back and
  /// nothing moved. The chains are then as install left them: a chain it
  /// had finished holds trees built on the stored grid over a store that is
  /// back on its own, so after a throw that follows an install the sampler
  /// must not be read or run until a state is installed. \p repair is what
  /// a monotone tree the move leaves out of order gets: the reset for a
  /// state, the collapse for a warm start.
  template <typename Install>
  bool installAcrossGrids(
      const std::vector<std::vector<double>>& storedCutPoints,
      const double* currentPredictors, Install&& install,
      MonotoneRepair repair = MonotoneRepair::reset) {
    std::vector<std::vector<double>> moved(data_.numPredictors);
    for (size_t j = 0; j < data_.numPredictors; ++j)
      if (data_.splitsByThreshold(j) &&
          storedCutPoints[j] != data_.cutPoints[j])
        moved[j] = storedCutPoints[j];
    GridSnapshot own(data_);
    bool installed = false;
    try {
      for (size_t j = 0; j < data_.numPredictors; ++j)
        if (!moved[j].empty())
          data_.setCutPointsForColumn(
            j, moved[j].data(), static_cast<std::uint32_t>(moved[j].size()),
            currentPredictors);
      installed = install();
    } catch (...) {
      own.restore(data_);
      throw;
    }
    own.restore(data_);
    if (!installed) return false;
    for (auto& chain : chains_)
      chain->forceRefreshTrees(&moved, SplitPlacement::byValue, repair);
    return true;
  }

  /// Shared constructor tail: one chain per rng over the built store, plus
  /// saved-tree storage under keepTrees.
  void initializeChains(const double* y, const double* weights,
                        const double* offset, double sigmaEstimate,
                        double sigmaDf, double sigmaRawScale,
                        ext_rng* const* rngs) {
    chains_.reserve(options_.numChains);
    for (size_t c = 0; c < options_.numChains; ++c)
      chains_.push_back(std::make_unique<Chain<L, ResidT>>(
        data_, y, weights, offset, family_, sigmaEstimate, sigmaDf,
        sigmaRawScale, options_, rngs[c]));
    options_.survivalStatus = nullptr;  // borrowed; consumed by the chains

    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
    recordAnchor();
  }

  /// The sampler's transform becomes the one its chains were just put in;
  /// every chain shares it after a creation or a re-anchor.
  void recordAnchor() { chains_[0]->getScale(anchorMin_, anchorMax_); }

  /// Raises the flags options_.missingSeen names over the ones the build read
  /// from the data. Runs before any chain is built: no tree holds a rule yet,
  /// so no direction is owed a draw.
  void adoptMissingSeen() {
    if (options_.missingSeen != nullptr)
      for (size_t j = 0; j < data_.numPredictors; ++j)
        data_.hasMissing[j] |= options_.missingSeen[j] != 0 ? 1 : 0;
    options_.missingSeen = nullptr;  // borrowed
  }

  /// The first-sight draw across the chains, each from its own generator
  /// (Chain::drawMissingDirections): for the columns whose flag is up now and
  /// was down in \p flagsBefore, the store's flags from before the change
  /// that has just been applied. A change that raised none draws nothing and
  /// touches nothing. \p records, when non-null, receives one record per
  /// chain for takeBackMissingDirections, and is left empty where nothing was
  /// drawn.
  void drawFirstSightDirections(
      const std::vector<std::uint8_t>& flagsBefore,
      std::vector<MissingDirectionDraw>* records = nullptr) {
    // flags only rise, so the two differ exactly where one rose
    if (data_.hasMissing == flagsBefore) return;
    std::vector<std::uint8_t> raised(data_.numPredictors);
    for (size_t j = 0; j < data_.numPredictors; ++j)
      raised[j] = data_.hasMissing[j] != flagsBefore[j] ? 1 : 0;
    SavedDrawSlots kept = savedDrawSlots();
    if (records != nullptr) records->resize(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->drawMissingDirections(
        raised.data(), kept, records != nullptr ? &(*records)[c] : nullptr);
  }

  /// Undoes drawFirstSightDirections for a change that is refused: every
  /// direction, live and saved, and every chain's generator as they were.
  void takeBackMissingDirections(std::vector<MissingDirectionDraw>& records) {
    for (size_t c = 0; c < records.size(); ++c) {
      chains_[c]->swapMissingDirections(records[c]);
      chains_[c]->releaseMissingDirections(records[c]);
    }
    records.clear();
  }

  /// Two-phase transaction over every chain: validate all trees of all
  /// forests of all chains first - every leaf occupied and every monotone
  /// tree's leaf values in order - then rebuild fits only if everything holds,
  /// so a failure in a late chain never leaves an early chain's fits
  /// overwritten. columns names the touched columns (null for a whole-matrix
  /// swap), which prunes the trees of forests past the first that cannot be
  /// affected; the survivor lists travel from phase one to phase two inside
  /// the handoff, so both phases quantify over the same set.
  ///
  /// movedCutPoints, when non-null, holds the old grid of each column whose
  /// splits move onto the grid now held (splitsToMove); phase one moves them
  /// and fails where a split's interval is empty. A failure puts every moved
  /// position back before it returns, so the caller's restore of the store
  /// and a repartition leave the trees as they were.
  bool revalidateAllChains(
      const size_t* columns, size_t numColumns,
      const std::vector<std::vector<double>>* movedCutPoints = nullptr,
      SplitPlacement placement = SplitPlacement::byPosition) {
    size_t numChains = chains_.size();
    std::vector<typename Chain<L, ResidT>::ForestRevalidation> state(numChains);

    bool allValid = true;
    for (size_t c = 0; c < numChains && allValid; ++c)
      allValid = chains_[c]->revalidateTrees(state[c], columns, numColumns,
                                             movedCutPoints, placement);
    if (!allValid) {
      for (size_t c = 0; c < numChains; ++c)
        Chain<L, ResidT>::undoSplitMoves(state[c]);
      return false;
    }

    for (size_t c = 0; c < numChains; ++c)
      chains_[c]->rebuildFitsFromParameters(state[c]);
    return true;
  }

  /// The two ways a predictor-replacement transaction moves its columns: the
  /// column list each spans (a null columns pointer means every predictor in
  /// order, view column j filling store column j) and the
  /// snapshot/apply/restore triple runPredictorTransaction drives.
  /// WholeMatrixUpdate swaps the whole code storage aside; SubsetUpdate
  /// journals a named subset column by column.
  struct WholeMatrixUpdate {
    ColumnStore& data;
    PredictorSource source;
    const size_t* columns = nullptr;
    std::vector<xint_t> oldCodes;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
    std::vector<std::uint32_t> oldNumCuts;
    std::vector<std::vector<double>> oldCutMass;
    // CSC/mixed rollback: a sparse column's storage lives outside train.codes
    // (rank bitmaps in sparseColumns, the borrowed/owned slice in sources, the
    // mutation-owned nonzero buffers), so snapshot it alongside the codes when
    // the store carries CSC-backed columns
    bool csc = false;
    std::vector<SparseColumnData> oldSparseColumns;
    std::vector<ColumnSource> oldSources;
    std::vector<std::vector<int>> oldOwnedRows;
    std::vector<std::vector<double>> oldOwnedValues;
    std::vector<std::uint8_t> oldCscOwned;

    size_t numColumns() const { return data.numPredictors; }
    void applyForced(bool updateCuts) { data.setPredictors(source, updateCuts); }
    void snapshotApply(bool updateCuts) {
      // move the live codes aside and rebuild into fresh storage: a reject swaps
      // them back and an accept drops them, so no whole-matrix copy survives
      oldCodes = std::move(data.train.codes);
      data.train.codes.assign(oldCodes.size(), 0);
      oldHasMissing = data.hasMissing;
      if (updateCuts) {
        // a refreshed grid may hold another count than the one it replaces
        oldCuts = data.cutPoints;
        oldNumCuts = data.numCuts;
        oldCutMass = data.cutMass;
      }
      csc = data.builtFromCsc;
      if (csc) {
        oldSparseColumns = data.train.sparseColumns;
        oldSources = data.train.sources;
        oldOwnedRows = data.ownedCscRows;
        oldOwnedValues = data.ownedCscValues;
        oldCscOwned = data.cscColumnOwned;
      }
      data.setPredictors(source, updateCuts);
    }
    void restore(bool updateCuts) {
      data.train.codes = std::move(oldCodes);
      data.hasMissing = std::move(oldHasMissing);
      if (updateCuts) {
        data.cutPoints = std::move(oldCuts);
        data.numCuts = std::move(oldNumCuts);
        data.cutMass = std::move(oldCutMass);
      }
      if (csc) {
        data.train.sparseColumns = std::move(oldSparseColumns);
        data.train.sources = std::move(oldSources);
        data.ownedCscRows = std::move(oldOwnedRows);
        data.ownedCscValues = std::move(oldOwnedValues);
        data.cscColumnOwned = std::move(oldCscOwned);
        // the restored source slices still point at the pre-change buffers,
        // which the moves relocated; repoint the owned ones (borrowed slices
        // keep their valid R pointer)
        for (size_t j = 0; j < data.numPredictors; ++j)
          if (data.cscColumnOwned[j]) data.repointOwnedSlice(j);
      }
    }
  };

  struct SubsetUpdate {
    ColumnStore& data;
    PredictorSource source;
    const size_t* columns;
    size_t count;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
    std::vector<std::uint32_t> oldNumCuts;
    std::vector<std::vector<double>> oldCutMass;
    std::vector<ColumnStore::ColumnCodeRollback> records;
    // CSC-backed columns of the subset snapshot their sparse storage here
    // instead of the dense per-cell journal (which has no codes[] to journal
    // for a rank column)
    std::vector<ColumnStore::CscColumnRollback> cscRecords;

    size_t numColumns() const { return count; }
    void applyForced(bool updateCuts) {
      data.setColumns(source, columns, count, updateCuts);
    }
    void snapshotApply(bool updateCuts) {
      // journal each touched column's changed cells (past a quarter changed the
      // record falls back to a whole pre-change copy); the small hasMissing and
      // cut pieces snapshot per column. A CSC-backed column instead snapshots
      // its rank/densified storage and owned slice whole.
      std::vector<double> scratch;
      oldHasMissing.resize(count);
      oldCuts.resize(updateCuts ? count : 0);
      oldNumCuts.resize(updateCuts ? count : 0);
      oldCutMass.resize(updateCuts ? count : 0);
      records.resize(count);
      cscRecords.resize(count);
      for (size_t k = 0; k < count; ++k) {
        size_t j = columns[k];
        oldHasMissing[k] = data.hasMissing[j];
        if (updateCuts) {
          oldCuts[k] = data.cutPoints[j];
          oldNumCuts[k] = data.numCuts[j];
          oldCutMass[k] = data.cutMass[j];
        }
        bool csc = data.columnIsCscBacked(j);
        if (csc) data.snapshotCscColumn(j, cscRecords[k]);
        data.mutateColumnFromSource(j, source, k, updateCuts, scratch,
                                    csc ? nullptr : &records[k]);
      }
    }
    /// Walks the records in reverse: record k was taken after records 0 to
    /// k - 1 were applied, so a column named twice unwinds to its pre-call
    /// state only in that order.
    void restore(bool updateCuts) {
      for (size_t k = count; k-- > 0;) {
        size_t j = columns[k];
        if (data.columnIsCscBacked(j))
          data.restoreCscColumn(j, cscRecords[k]);
        else
          data.restoreColumn(j, records[k]);
        data.hasMissing[j] = oldHasMissing[k];
        if (updateCuts) {
          data.cutPoints[j] = std::move(oldCuts[k]);
          data.numCuts[j] = oldNumCuts[k];
          data.cutMass[j] = std::move(oldCutMass[k]);
        }
      }
    }
  };

  /// One predictor-replacement transaction over a strategy's column list:
  /// precheck that the cuts stay representable, then either collapse emptied
  /// leaves under forceUpdate or snapshot-apply and keep the change when
  /// revalidateAllChains holds, else restore the snapshot and repartition;
  /// rolledBack reports an emptied leaf and a monotone tree out of order
  /// alike. The engine keeps no predictor matrix, so the raw a reject must put
  /// back is snapshotted here: the gathered leaf-covariate copies and, for a
  /// mixed store, the owned dense block of the columns the strategy touches
  /// (the strategy's own records carry no raw). The strategy owns the codes,
  /// missing flags, and cut grids it moves.
  ///
  /// A column the change first gives a missing value draws the direction of
  /// every rule already on it, live and saved (drawFirstSightDirections),
  /// before the trees are refreshed or revalidated, so a collapse merges, and
  /// occupancy and the monotone order are judged, under the directions the
  /// rows are routed by. A rollback takes the draw back with the rest.
  ///
  /// A refreshed grid is derived from the new values alone, so it may differ
  /// from the one it replaces in its points and in their count. Where it
  /// does, placement moves the column's splits onto it (splitsToMove): by
  /// position each keeps its place, rescaled when the count changed and
  /// untouched otherwise; by value each takes the point nearest its old
  /// threshold. Either stays inside the interval the split's ancestors leave.
  /// A split whose interval is empty after the move is merged under
  /// forceUpdate and rolls the transaction back otherwise, the moved
  /// positions put back with the grids.
  template <typename Strategy>
  PredictorUpdateResult runPredictorTransaction(
      Strategy& strategy, bool forceUpdate, bool updateCutPoints,
      SplitPlacement placement = SplitPlacement::byPosition) {
    std::vector<double> scratch;
    for (size_t k = 0; k < strategy.numColumns(); ++k) {
      size_t j = strategy.columns ? strategy.columns[k] : k;
      // a factor column of either kind must hold representable level codes
      // whether or not cut points refresh: its grid is its level table, so a
      // value outside the table has no position on it
      if ((updateCutPoints || data_.isFactor(j)) &&
          !data_.cutsWouldRemainValid(j, strategy.source, k, scratch))
        return PredictorUpdateResult::invalidCutPoints;
    }

    // the grid each refreshed numeric column holds before the change
    std::vector<std::vector<double>> oldCutPoints;
    if (updateCutPoints) {
      oldCutPoints.resize(data_.numPredictors);
      for (size_t k = 0; k < strategy.numColumns(); ++k) {
        size_t j = strategy.columns ? strategy.columns[k] : k;
        if (!data_.isFactor(j) && oldCutPoints[j].empty())
          oldCutPoints[j] = data_.cutPoints[j];
      }
    }

    const std::vector<std::uint8_t> flagsBefore(data_.hasMissing);

    if (forceUpdate) {
      strategy.applyForced(updateCutPoints);
      drawFirstSightDirections(flagsBefore);
      const std::vector<std::vector<double>>* moved =
        updateCutPoints ? splitsToMove(oldCutPoints, placement) : nullptr;
      for (auto& chain : chains_) chain->forceRefreshTrees(moved, placement);
      requantizeTestColumns(strategy, updateCutPoints);
      return PredictorUpdateResult::accepted;
    }

    std::vector<double> oldGatheredRaw = data_.gatheredRawValues;
    ColumnStore::OwnedDenseRollback oldOwnedDense;
    data_.snapshotOwnedDenseColumns(strategy.columns, strategy.numColumns(),
                                    oldOwnedDense);
    strategy.snapshotApply(updateCutPoints);
    std::vector<MissingDirectionDraw> firstSight;
    drawFirstSightDirections(flagsBefore, &firstSight);
    const std::vector<std::vector<double>>* moved =
      updateCutPoints ? splitsToMove(oldCutPoints, placement) : nullptr;
    if (!revalidateAllChains(strategy.columns, strategy.numColumns(), moved,
                             placement)) {
      takeBackMissingDirections(firstSight);
      data_.gatheredRawValues = std::move(oldGatheredRaw);
      data_.restoreOwnedDenseColumns(oldOwnedDense);
      strategy.restore(updateCutPoints);
      for (auto& chain : chains_) chain->repartitionTrees();
      return PredictorUpdateResult::rolledBack;
    }
    requantizeTestColumns(strategy, updateCutPoints);
    return PredictorUpdateResult::accepted;
  }

  /// Re-quantize the transaction's columns into the test block after a
  /// cut-refreshing accept, matching the train-side requantize.
  template <typename Strategy>
  void requantizeTestColumns(Strategy& strategy, bool updateCutPoints) {
    if (updateCutPoints && data_.numTestObservations > 0)
      for (size_t k = 0; k < strategy.numColumns(); ++k)
        data_.quantizeTestColumn(strategy.columns ? strategy.columns[k] : k);
  }

  /// Caches each observation's leaf and per-leaf occupancy for every tree of
  /// every forest - and of the variance forest - of every chain that splits on
  /// the session's column, by routing codes through the split structure, so
  /// staging a move is a descent and a count check per cached tree.
  class UpdateSessionImpl final : public PredictorUpdateSession {
  public:
    UpdateSessionImpl(Sampler& sampler, const double* newColumn, size_t column)
      : sampler_(sampler), column_(column), newColumn_(newColumn) {
      size_t n = sampler_.data_.numObservations;

      newCodes_.resize(n);
      for (size_t i = 0; i < n; ++i)
        newCodes_[i] = sampler_.data_.codeFor(column_, newColumn_[i]);

      // The cached set spans every forest and the variance forest, and is
      // pruned to the trees the column can move, so it is sparse: an explicit
      // (chain, forest, tree) table rather than offset arithmetic over a
      // rectangular tree count.
      std::vector<std::uint32_t> census;
      std::vector<size_t> splitting;
      for (size_t c = 0; c < sampler_.chains_.size(); ++c) {
        const Chain<L, ResidT>& chain(*sampler_.chains_[c]);
        for (size_t f = 0; f < chain.numForests(); ++f) {
          chain.treesSplittingOnColumn(f, column_, census, splitting);
          for (size_t t : splitting) cached_.push_back(CachedTree{c, f, t});
        }
        if (!chain.hasVarianceForest()) continue;
        chain.varianceTreesSplittingOnColumn(column_, census, splitting);
        for (size_t j : splitting)
          cached_.push_back(CachedTree{c, varianceForest, j});
      }

      size_t numCachedTrees = cached_.size();
      observationLeaf_.resize(numCachedTrees * n);
      leafCounts_.resize(numCachedTrees);
      pendingNewLeaf_.assign(numCachedTrees, invalidNode);
      pendingOldLeaf_.resize(numCachedTrees);

      for (size_t t = 0; t < numCachedTrees; ++t) {
        const Tree& tree(treeAt(t));
        leafCounts_[t].assign(tree.nodes.size(), 0);
        for (size_t i = 0; i < n; ++i) {
          int32_t leaf = tree.findBottomNodeForObservation(sampler_.data_, i);
          observationLeaf_[t * n + i] = leaf;
          ++leafCounts_[t][static_cast<size_t>(leaf)];
        }
      }
    }

    ~UpdateSessionImpl() override { settleFirstSight(); }

    bool observationWouldRemainValid(size_t i) override {
      size_t n = sampler_.data_.numObservations;
      takeBackUncommitted();
      // the column's first missing value, by the test setCell marks the
      // column with: the flag goes up and the rules on the column draw their
      // directions now, so that this row's leaves and the monotone order are
      // judged under them. Nothing is written that a decline does not take
      // back.
      bool firstMissing =
        isNA(newColumn_[i]) && !sampler_.data_.hasMissing[column_];
      if (firstMissing) {
        if (orderWithMissing_ == OrderWithMissing::breaks) return false;
        putOnFirstSight();
        if (!orderHoldsWithMissing()) {
          takeBackFirstSight();
          return false;
        }
      }
      bool valid = true;
      for (size_t t = 0; t < leafCounts_.size() && valid; ++t) {
        const Tree& tree(treeAt(t));
        int32_t oldLeaf = observationLeaf_[t * n + i];
        int32_t newLeaf = tree.findBottomNodeForObservation(
          sampler_.data_, i, static_cast<int32_t>(column_), newCodes_[i]);

        if (newLeaf == oldLeaf) {
          pendingNewLeaf_[t] = invalidNode;
        } else if (leafCounts_[t][static_cast<size_t>(oldLeaf)] == 1) {
          valid = false;  // last occupant; the move would empty the leaf
        } else {
          pendingOldLeaf_[t] = oldLeaf;
          pendingNewLeaf_[t] = newLeaf;
        }
      }
      if (firstMissing) {
        if (valid) firstSight_ = FirstSight::uncommitted;
        else takeBackFirstSight();
      }
      return valid;
    }

    void commitObservation(size_t i) override {
      // a first missing value's draw stands with its row
      if (firstSight_ == FirstSight::uncommitted) {
        firstSight_ = FirstSight::none;
        firstSightDraws_.clear();
      }
      // the session's one cell-write: setCell re-quantizes newColumn_[i] to the
      // code the valid-move check descended on (cuts fixed) and marks hasMissing
      sampler_.data_.setCell(i, column_, newColumn_[i]);
      for (size_t t = 0; t < leafCounts_.size(); ++t) {
        if (pendingNewLeaf_[t] != invalidNode) {
          --leafCounts_[t][static_cast<size_t>(pendingOldLeaf_[t])];
          ++leafCounts_[t][static_cast<size_t>(pendingNewLeaf_[t])];
        }
      }
    }

    /// The session moves exactly one column, so the revalidation prunes
    /// against that column alone - the same predicate that built the cache.
    /// The cache is the subset of the revalidated set that drops forest 0's
    /// non-splitting trees, whose partitions the revalidation reproduces
    /// unchanged, so no leaf it did not guard can empty and this returns true
    /// by construction. The variance arm is guarded and revalidated on the same
    /// predicate with no exemption, so the two sets coincide there. The
    /// monotone order the revalidation also judges can move only with the
    /// column's has-missing flag, which the guard lets rise only where the
    /// order holds under the directions drawn with it.
    bool finalize() override {
      settleFirstSight();
      return sampler_.revalidateAllChains(&column_, 1);
    }

  private:
    /// One cached tree, addressed as the engine holds it. The variance forest
    /// has no index in forests_, so it takes a sentinel in the forest slot.
    static constexpr size_t varianceForest = static_cast<size_t>(-1);
    struct CachedTree {
      size_t chain, forest, tree;
    };

    /// Whether every tree of every chain stays in the monotone order with
    /// the column's flag up and its first-sight directions on, which is how
    /// it is asked. Leaf values do not move within a session, the order reads
    /// no partition and the directions are the same each time they are put
    /// on, so the answer is the same for whichever row asks and is found
    /// once.
    bool orderHoldsWithMissing() {
      if constexpr (TreeDrawLeafModel<L>) {
        if (orderWithMissing_ == OrderWithMissing::unknown) {
          bool holds = true;
          for (size_t c = 0; c < sampler_.chains_.size() && holds; ++c)
            holds = sampler_.chains_[c]->monotoneLeavesInOrder();
          orderWithMissing_ =
            holds ? OrderWithMissing::holds : OrderWithMissing::breaks;
        }
        return orderWithMissing_ == OrderWithMissing::holds;
      } else {
        return true;
      }
    }

    /// Raises the column's flag and puts its first-sight directions on: drawn
    /// the first time, each chain from its own generator, and the same draw
    /// restored after that, so every row that asks is judged under the
    /// directions the first one drew and the generators end where one draw
    /// leaves them.
    void putOnFirstSight() {
      sampler_.data_.hasMissing[column_] = 1;
      size_t numChains = sampler_.chains_.size();
      if (firstSightDraws_.empty()) {
        std::vector<std::uint8_t> raised(sampler_.data_.numPredictors, 0);
        raised[column_] = 1;
        SavedDrawSlots kept = sampler_.savedDrawSlots();
        firstSightDraws_.resize(numChains);
        for (size_t c = 0; c < numChains; ++c)
          sampler_.chains_[c]->drawMissingDirections(raised.data(), kept,
                                                     &firstSightDraws_[c]);
      } else {
        for (size_t c = 0; c < numChains; ++c)
          sampler_.chains_[c]->swapMissingDirections(firstSightDraws_[c]);
      }
    }

    /// Lowers the flag and takes the directions and the generators back to
    /// where they were, keeping the draw for the next row that needs it.
    void takeBackFirstSight() {
      for (size_t c = 0; c < firstSightDraws_.size(); ++c)
        sampler_.chains_[c]->swapMissingDirections(firstSightDraws_[c]);
      sampler_.data_.hasMissing[column_] = 0;
      firstSight_ = FirstSight::none;
    }

    /// A row found valid whose commit never came - another session of a
    /// joint update declined it - leaves its draw on; it comes off before
    /// anything else is judged.
    void takeBackUncommitted() {
      if (firstSight_ == FirstSight::uncommitted) takeBackFirstSight();
    }

    /// Ends the session's hold on a draw no row kept: the directions are
    /// already back, and the pool words the draw took go with it.
    void settleFirstSight() {
      takeBackUncommitted();
      for (size_t c = 0; c < firstSightDraws_.size(); ++c)
        sampler_.chains_[c]->releaseMissingDirections(firstSightDraws_[c]);
      firstSightDraws_.clear();
    }

    const Tree& treeAt(size_t t) const {
      const CachedTree& entry(cached_[t]);
      const Chain<L, ResidT>& chain(*sampler_.chains_[entry.chain]);
      return entry.forest == varianceForest
        ? chain.varianceTree(entry.tree)
        : chain.treeInForest(entry.forest, entry.tree);
    }

    Sampler& sampler_;
    size_t column_;
    const double* newColumn_;
    std::vector<xint_t> newCodes_;
    std::vector<CachedTree> cached_;
    std::vector<int32_t> observationLeaf_;  // cached_.size() x numObservations
    std::vector<std::vector<std::uint32_t>> leafCounts_;  // arena-id indexed
    std::vector<int32_t> pendingNewLeaf_;  // invalidNode when no move staged
    std::vector<int32_t> pendingOldLeaf_;
    enum class OrderWithMissing { unknown, holds, breaks };
    OrderWithMissing orderWithMissing_ = OrderWithMissing::unknown;
    // the column's first-sight draw, one record per chain: empty until a
    // first missing row draws it and again once a row commits it, and
    // otherwise held, on or off, for the session's later missing rows
    std::vector<MissingDirectionDraw> firstSightDraws_;
    // uncommitted: the draw is on for a row found valid and not yet committed
    enum class FirstSight { none, uncommitted };
    FirstSight firstSight_ = FirstSight::none;
  };

  SamplerOptions options_;
  ResponseFamily family_;
  ColumnStore data_;
  std::vector<std::unique_ptr<Chain<L, ResidT>>> chains_;
  size_t currentSampleNum_ = 0;  // next saved-tree slot, wrapping circularly
  // the sampler's response transform, held apart from each chain's: what the
  // leaf prior's k is relative to and what an installed state is read
  // against. Set by creation, a re-anchor and setAnchor; an install never
  // moves it. (0, 0) on a scale-free family.
  double anchorMin_ = 0.0, anchorMax_ = 0.0;
  size_t recordedDraws_ = 0;     // slots written since the last reset, capped
                                 // at capacity; the extent of every read
  double runningTime_ = 0.0;     // seconds accumulated across runs
  // set when the creation build refused a predictor cell; the sampler is then
  // an empty shell its factory discards, so no chain is built over the store
  bool ingestionRefused_ = false;
};

using ConstantLeafSampler = Sampler<ConstantGaussianLeaf>;

}  // namespace bartcore

#endif  // BARTCORE_SAMPLER_HPP
