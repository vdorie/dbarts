#ifndef BARTCORE_SAMPLER_HPP
#define BARTCORE_SAMPLER_HPP

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <mutex>
#ifndef _WIN32
#include <csignal>  // block SIGINT in worker threads; see run()
#endif
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

/// A whole sampler's serializable state: per-chain states, the store's cut
/// points (setCutPoints may have replaced the ones creation induces), the
/// saved-tree write position, and the draw count that position is read
/// against.
struct SamplerStateData {
  std::vector<ChainStateData> chains;
  std::vector<std::vector<double>> cutPoints;  // empty vector per categorical column
  size_t currentSampleNum = 0;
  /// Recorded draws the store retains, capped at its capacity. Required: it
  /// cannot be inferred from the buffer, whose unwritten slots hold zero-leaf
  /// trees that read as legitimate draws.
  size_t recordedDraws = 0;
};

/// Why a warm start (installForests) refused; ok on success. A single donor
/// forest can seed several chains, so the donor's chain count need not match.
enum class WarmStartResult {
  ok, shapeMismatch, gridMismatch, dartMismatch, interactionMismatch,
  columnMaskMismatch, varianceMismatch, varianceSlotMismatch
};

/// A sequential per-observation predictor update: stage one observation's
/// leaf moves, test that no leaf empties, then commit or skip, with a single
/// fits rebuild at the end. Type-erased so several samplers sharing an
/// index-aligned column can be swept jointly, committing all-or-none.
class PredictorUpdateSession {
public:
  virtual ~PredictorUpdateSession() = default;
  /// Stages observation i's leaf moves against the running occupancy counts;
  /// true unless some tree's leaf would be left empty.
  virtual bool observationWouldRemainValid(size_t i) = 0;
  /// Installs observation i's staged value and occupancy moves.
  virtual void commitObservation(size_t i) = 0;
  /// Re-routes all observations and rebuilds fits. The sequential guard
  /// admits no empty leaves, so false is an internal invariant violation.
  virtual bool finalize() = 0;
};

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
    data_.build(creationPredictorSource(options.predictors, x,
                                        numObservations, numPredictors),
                options.maxNumCutsPerVariable, options.maxNumCuts,
                options.useQuantiles, options.leafCovariateColumns,
                options.numLeafCovariates);
    options_.maxNumCutsPerVariable = nullptr;  // borrowed; consumed by build
    // borrowed; consumed by build, which retains what the store needs (a
    // mapped build's CSC slices, its own copy of the dense block)
    options_.predictors = {};

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
    data_.build(denseCreationPredictorSource(options.predictors, x,
                                             numObservations, numPredictors),
                options.maxNumCutsPerVariable, options.maxNumCuts,
                options.useQuantiles);
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
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
    data_.build(denseCreationPredictorSource(options.predictors, x,
                                             numObservations, numPredictors),
                options.maxNumCutsPerVariable, options.maxNumCuts,
                options.useQuantiles);
    options_.maxNumCutsPerVariable = nullptr;
    options_.predictors = {};
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
  }

  // chains reference the store member, so the sampler's address is pinned
  Sampler(const Sampler&) = delete;
  Sampler& operator=(const Sampler&) = delete;

  /// Replace the test predictors from a borrowed view, keeping any test offset:
  /// the caller guarantees the row count still matches it (the bridge refuses
  /// otherwise). Passing a new offset too goes through setTestOffset. The test
  /// store shares the training cut grid and owns its raw, so the view need not
  /// outlive the call. Refuses (returns false, test store untouched) a
  /// designated leaf covariate that would be CSC-backed, since leaf models
  /// gather dense raw test covariates that sparse storage does not serve.
  bool setTestData(const PredictorSource& source) {
    size_t numCovariates = numLeafCovariates();
    const size_t* covariateColumns = leafCovariateColumns();
    for (size_t k = 0; k < numCovariates; ++k)
      if (source.sourceOf(covariateColumns[k]) < 0) return false;
    const double* testOffset = data_.testOffset;
    data_.buildTest(source);
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
  /// clock read per sweep and stays bitwise identical.
  ///
  /// onSweep, if set, is the host's per-sweep conditioning hook; it runs only
  /// when chains run inline (min(numThreads, numChains) <= 1), so the caller
  /// must not set it alongside worker-thread chains. Chains then run
  /// sequentially, so onSweep sees chain c completed before chain c + 1 begins.
  bool run(size_t numBurnIn, size_t numSamples, Results& results,
           const std::function<bool()>& pollInterrupt = {},
           const SweepCallback& onSweep = {}) {
    size_t numChains = chains_.size();
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
    std::vector<Results> chainResults(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      Results& r(chainResults[c]);
      r.numReportedLocations = numLocations;
      r.numVariableCountForests = numVarCountForests;
      r.numOrdinalThresholds = numOrdinalThresholds;
      if (results.sigma != nullptr) r.sigma = results.sigma + c * numSamples;
      if (results.k != nullptr) r.k = results.k + c * numSamples;
      if (results.trainingFits != nullptr)
        r.trainingFits = results.trainingFits +
          c * numSamples * data_.numObservations * numLocations;
      if (results.testFits != nullptr)
        r.testFits = results.testFits +
          c * numSamples * data_.numTestObservations * numLocations;
      if (results.variableCounts != nullptr)
        r.variableCounts = results.variableCounts +
          c * numSamples * data_.numPredictors * numVarCountForests;
      if (results.splitProbabilities != nullptr)
        r.splitProbabilities =
          results.splitProbabilities + c * numSamples * data_.numPredictors;
      if (results.tau != nullptr) r.tau = results.tau + c * numSamples;
      if (results.groupEffects != nullptr)
        r.groupEffects =
          results.groupEffects + c * numSamples * options_.numGroups;
      if (results.logLikelihood != nullptr)
        r.logLikelihood =
          results.logLikelihood + c * numSamples * data_.numObservations;
      if (results.ordinalThresholds != nullptr)
        r.ordinalThresholds =
          results.ordinalThresholds + c * numSamples * numOrdinalThresholds;
      if (results.dispersion != nullptr)
        r.dispersion = results.dispersion + c * numSamples;
      if (results.residualDf != nullptr)
        r.residualDf = results.residualDf + c * numSamples;
      if (results.varianceFits != nullptr)
        r.varianceFits =
          results.varianceFits + c * numSamples * data_.numObservations;
      if (results.varianceTestFits != nullptr)
        r.varianceTestFits = results.varianceTestFits +
          c * numSamples * data_.numTestObservations;
      // the per-forest slab carries every forest's own fits per sample, and the
      // glue slab the ragged amplitude vector; both are filled only by a
      // coupling that defines them, so the strides are reached on that path
      // alone
      if (results.forestFits != nullptr)
        r.forestFits = results.forestFits +
          c * numSamples * data_.numObservations * numForests();
      if (results.glue != nullptr)
        r.glue = results.glue + c * numSamples * totalAmplitudes();
    }

    if (options_.verbose) ext_printf("Running mcmc loop:\n");
    std::chrono::steady_clock::time_point startTime =
      std::chrono::steady_clock::now();

    size_t numWorkers = options_.numThreads < numChains ? options_.numThreads
                                                        : numChains;
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
                                    &progress, c, shouldCancelPtr, onSweepPtr);
    } else {
      // workers never call into R: progress lines queue and the main thread
      // flushes them every 0.1 seconds and polls for interrupts, setting the
      // cancel flag the workers read (a relaxed atomic load per sweep)
      QueuedProgressSink progress;
      std::atomic<size_t> numChainsRunning(numChains);
      std::atomic<bool> cancelFlag(false);
      std::function<bool()> workerCancel = [&cancelFlag]() {
        return cancelFlag.load(std::memory_order_relaxed);
      };
      std::vector<std::thread> workers;
      workers.reserve(numWorkers);
#ifndef _WIN32
      // spawn the workers with SIGINT blocked so they inherit the block and a
      // Ctrl-C is delivered only to this (the main) thread, whose poll turns
      // it into a cooperative cancel. A worker running R's interrupt handler
      // could longjmp across threads. The main thread's mask is restored right
      // after, before it polls. (On Windows R's console Ctrl-C already runs on
      // the main thread, so no masking is needed.)
      sigset_t interruptSet, previousSet;
      sigemptyset(&interruptSet);
      sigaddset(&interruptSet, SIGINT);
      pthread_sigmask(SIG_BLOCK, &interruptSet, &previousSet);
#endif
      for (size_t w = 0; w < numWorkers; ++w) {
        workers.emplace_back([this, w, numWorkers, numChains, numBurnIn,
                              numSamples, &chainResults, &progress,
                              &numChainsRunning, &workerCancel]() {
          for (size_t c = w; c < numChains; c += numWorkers) {
            chains_[c]->run(numBurnIn, numSamples, chainResults[c], &progress,
                            c, &workerCancel);
            numChainsRunning.fetch_sub(1);
          }
        });
      }
#ifndef _WIN32
      pthread_sigmask(SIG_SETMASK, &previousSet, nullptr);
#endif
      while (numChainsRunning.load() > 0) {
        if (pollInterrupt && !cancelFlag.load(std::memory_order_relaxed) &&
            pollInterrupt())
          cancelFlag.store(true, std::memory_order_relaxed);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        if (options_.verbose) progress.flush();
      }
      for (std::thread& worker : workers) worker.join();
      if (options_.verbose) progress.flush();
      cancelled = cancelFlag.load(std::memory_order_relaxed);
    }

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
  /// the caller's thread. A traversal is one (row, tree, slab) descent, and
  /// the measured cost of one is ~2.6 ns, so 1e7 of them is ~26 ms of work -
  /// comfortably above the few tens of microseconds a spawn and a join cost,
  /// and small enough that no interactive call waits on the threshold. Each
  /// entry point counts ITS OWN traversals, not forest 0's.
  static constexpr size_t predictParallelCutoff = 10000000;

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
  /// the first message is raised after the join, on the caller's thread.
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
    size_t cutoff = predictPartition.cutoffOverride != 0
      ? predictPartition.cutoffOverride : predictParallelCutoff;
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
    workers.reserve(numWorkers);
#ifndef _WIN32
    sigset_t interruptSet, previousSet;
    sigemptyset(&interruptSet);
    sigaddset(&interruptSet, SIGINT);
    pthread_sigmask(SIG_BLOCK, &interruptSet, &previousSet);
#endif
    for (size_t w = 0, begin = 0; w < numWorkers; ++w) {
      size_t end = begin + base + (w < remainder ? 1 : 0);
      workers.emplace_back([&body, &scratch, &failed, &firstError, w, begin,
                            end]() {
        try {
          for (size_t slab = begin; slab < end; ++slab)
            body(slab, scratch[w]);
        } catch (const std::exception& error) {
          failed[w] = 1;
          firstError[w] = error.what();
        } catch (...) {
          failed[w] = 1;
          firstError[w] = "unknown exception";
        }
      });
      begin = end;
    }
#ifndef _WIN32
    pthread_sigmask(SIG_SETMASK, &previousSet, nullptr);
#endif
    for (std::thread& worker : workers) worker.join();
    for (size_t w = 0; w < numWorkers; ++w)
      if (failed[w] != 0)
        ext_throwError("predict worker %lu failed: %s",
                       static_cast<unsigned long>(w), firstError[w].c_str());
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
  /// Test hook passthrough: chain c's variance tree j.
  const Tree& varianceTreeForTesting(size_t chainNum, size_t j) const {
    return chains_[chainNum]->varianceTree(j);
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
  // setState validates every chain against the state's cut points before
  // touching anything, so failure leaves the sampler unchanged.

  void getState(SamplerStateData& state) {
    state.chains.resize(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->getState(state.chains[c]);
    state.cutPoints = data_.cutPoints;
    state.currentSampleNum = currentSampleNum_;
    state.recordedDraws = recordedDraws_;
  }

  /// currentPredictors is the call-time predictor matrix a cross-grid column
  /// re-quantizes from (data@x, or the retained creation spec's @x); null for
  /// CSC/mixed stores. A same-spec restore (the continuation contract) matches
  /// the live grid column for column, so its per-column skip guard re-quantizes
  /// nothing and needs no raw.
  ///
  /// The containment rule, one law for both tree-install entries: a state may
  /// install a tree only if every split it carries lies inside the recipient
  /// forest's own column mask - each forest of forests_ (a moderator subset, a
  /// blocks() row) and the variance forest alike - because
  /// splitVariableLogProbability prices a rule against collectAvailableVariables,
  /// which drops the masked columns, so a forbidden split is mis-scored for as
  /// long as it lives. setState and installForests therefore share the one
  /// predicate, columnMaskStateFeasible, and report the one refusal; only LIVE
  /// trees are held to it, a saved slot being a replay target routed over new
  /// rows rather than over this partition. It runs against the state's own cut
  /// grid, where the state's splits resolve, and before any chain is touched,
  /// so a refusal leaves the sampler exactly as it was. columnMaskRefused, when
  /// non-null, separates that refusal from every other invalid state so the
  /// host can name it.
  bool setState(const SamplerStateData& state,
                const double* currentPredictors,
                bool* columnMaskRefused = nullptr) {
    if (columnMaskRefused != nullptr) *columnMaskRefused = false;
    if (state.chains.size() != chains_.size()) return false;
    if (state.cutPoints.size() != data_.numPredictors) return false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.splitsBySubset(j)) {
        if (!state.cutPoints[j].empty()) return false;
      } else if (state.cutPoints[j].empty() ||
                 state.cutPoints[j].size() > 65535) {
        return false;
      } else {
        for (size_t k = 1; k < state.cutPoints[j].size(); ++k)
          if (state.cutPoints[j][k] <= state.cutPoints[j][k - 1]) return false;
      }
    }

    // install the state's cuts, snapshotting for rollback: tree validity is
    // defined against them
    std::vector<std::vector<double>> oldCutPoints(data_.cutPoints);
    std::vector<std::uint32_t> oldNumCuts(data_.numCuts);
    std::vector<std::uint32_t> oldMaxNumCuts(data_.maxNumCuts);
    std::vector<xint_t> oldCodes(data_.train.codes);
    std::vector<xint_t> oldTestCodes(data_.test.codes);
    // rank columns re-quantize into their own storage, not codes; a rank-backed
    // test column likewise re-quantizes into test.sparseColumns (empty, hence a
    // no-op restore, on every dense-test path)
    std::vector<SparseColumnData> oldSparseColumns(data_.train.sparseColumns);
    std::vector<SparseColumnData> oldTestSparseColumns(data_.test.sparseColumns);
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.splitsBySubset(j)) continue;
      // a restored grid equal to the live one leaves the codes already correct
      // (the continuation contract), so skip its re-quantization and its raw
      if (state.cutPoints[j] == data_.cutPoints[j]) continue;
      data_.setCutPointsForColumn(j, state.cutPoints[j].data(),
                                  static_cast<std::uint32_t>(
                                    state.cutPoints[j].size()),
                                  currentPredictors);
    }

    // containment first, validity second: both judge the state against the grid
    // just installed, and the split-variable refusal is separated only so the
    // host can name it as installForests does
    bool columnMaskOk = true;
    for (size_t c = 0; c < chains_.size() && columnMaskOk; ++c)
      columnMaskOk = chains_[c]->columnMaskStateFeasible(state.chains[c]);
    bool allValid = columnMaskOk;
    for (size_t c = 0; c < chains_.size() && allValid; ++c)
      allValid = chains_[c]->stateIsValid(state.chains[c]);

    if (!allValid) {
      data_.cutPoints = std::move(oldCutPoints);
      data_.numCuts = std::move(oldNumCuts);
      data_.maxNumCuts = std::move(oldMaxNumCuts);
      data_.train.codes = std::move(oldCodes);
      data_.test.codes = std::move(oldTestCodes);
      data_.train.sparseColumns = std::move(oldSparseColumns);
      data_.test.sparseColumns = std::move(oldTestSparseColumns);
      if (columnMaskRefused != nullptr) *columnMaskRefused = !columnMaskOk;
      return false;
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      if (!chains_[c]->setState(state.chains[c])) return false;
    size_t capacity = savedTreeCapacity();
    currentSampleNum_ = capacity > 0 ? state.currentSampleNum % capacity : 0;
    recordedDraws_ =
      state.recordedDraws < capacity ? state.recordedDraws : capacity;
    return true;
  }

  /// Warm start: seed each chain's live forest from a donor sample. sampleMap
  /// selects, for destination chain c, a donor (chain, slot); slot < 0 takes
  /// the donor chain's live trees, else its saved slot (slot-major, one forest
  /// per slot). Only trees, sigma, k, and DART transfer - rng and auxiliary
  /// state stay fresh, so each chain evolves independently from its own
  /// stream. Under a variance forest the scale surface rides along, taken from
  /// the same slot as the mean forest. A donor on a different cut grid has its
  /// splits remapped onto this sampler's grid (starved splits collapse), as
  /// setData remaps a data replacement; the donor must still share this
  /// sampler's per-forest tree counts and DART mode. On any mismatch nothing is
  /// touched.
  WarmStartResult installForests(
      const SamplerStateData& donor,
      const std::vector<std::pair<size_t, int>>& sampleMap) {
    if (sampleMap.size() != chains_.size())
      return WarmStartResult::shapeMismatch;

    // A donor grown on a different cut grid is remapped onto this sampler's
    // grid at install (as setData remaps a data replacement's old splits onto
    // the freshly rebuilt cuts), collapsing splits the new grid starves; a
    // same-grid donor installs verbatim. Only a structural predictor mismatch
    // the remap cannot bridge - a categorical/continuous disagreement, or a
    // malformed donor grid - is refused here as gridMismatch.
    bool crossGrid = donor.cutPoints != data_.cutPoints;
    if (crossGrid) {
      if (donor.cutPoints.size() != data_.numPredictors)
        return WarmStartResult::gridMismatch;
      for (size_t j = 0; j < data_.numPredictors; ++j) {
        bool categorical = data_.splitsBySubset(j);
        if (categorical != donor.cutPoints[j].empty())
          return WarmStartResult::gridMismatch;
        for (size_t k = 1; k < donor.cutPoints[j].size(); ++k)
          if (donor.cutPoints[j][k] <= donor.cutPoints[j][k - 1])
            return WarmStartResult::gridMismatch;
      }
    }

    // the scale-leaf positivity law state validation holds a variance block to,
    // applied at BOTH install arms below for the reason it applies there: a
    // donor state is hand-buildable whichever arm reaches it, and a rebuild
    // scatters the leaf straight into a divisor, so a non-positive factor is a
    // broken model rather than a bad number
    auto scaleLeavesArePositive =
      [](const std::vector<std::vector<FlatNode>>& trees) {
        for (const std::vector<FlatNode>& tree : trees)
          for (const FlatNode& node : tree)
            if (flatKindOf(node) == FlatKind::leaf && !(node.value > 0.0))
              return false;
        return true;
      };

    std::vector<ChainStateData> install(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c) {
      size_t dc = sampleMap[c].first;
      int slot = sampleMap[c].second;
      if (dc >= donor.chains.size()) return WarmStartResult::shapeMismatch;
      const ChainStateData& src = donor.chains[dc];
      if (src.forests.size() != chains_[c]->numForests())
        return WarmStartResult::shapeMismatch;
      if (chains_[c]->usesDart() != !src.dartProbabilities.empty())
        return WarmStartResult::dartMismatch;
      // a donor and destination that disagree here would drop a scale surface
      // or leave one cold-started under a mean forest fitted with it; the
      // trees themselves ride the reassembled state below
      if (src.varianceTrees.empty() == chains_[c]->hasVarianceForest())
        return WarmStartResult::shapeMismatch;

      ChainStateData& dst = install[c];
      dst.forests.resize(src.forests.size());
      for (size_t f = 0; f < src.forests.size(); ++f) {
        const ForestStateData& sfs = src.forests[f];
        ForestStateData& dfs = dst.forests[f];
        size_t nt = chains_[c]->numTreesInForest(f);
        // the donor's tree count sets the saved slot's stride, so it must
        // match this forest's before either path proceeds
        if (sfs.trees.size() != nt) return WarmStartResult::shapeMismatch;
        if (slot < 0) {
          dfs.trees = sfs.trees;
          dfs.treeParams = sfs.treeParams;
          dfs.treeMasks = sfs.treeMasks;
        } else {
          size_t base = static_cast<size_t>(slot) * nt;
          if (sfs.savedTrees.size() < base + nt)
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
        // the per-forest leaf scale rides the warm start with k; this
        // reassembled ForestStateData is the ONLY thing installForest sees, so
        // dropping it here would leave that path's install arm dead
        dfs.leafScale = sfs.leafScale;
      }
      dst.sigma = src.sigma;
      dst.fitMin = src.fitMin;
      dst.fitMax = src.fitMax;
      dst.dartProbabilities = src.dartProbabilities;
      dst.dartAlpha = src.dartAlpha;
      dst.dartNumUpdatesSkipped = src.dartNumUpdatesSkipped;
      dst.hasAmplitudes = src.hasAmplitudes;
      dst.amplitudeWidths = src.amplitudeWidths;
      dst.amplitudes = src.amplitudes;
      dst.amplitudeVariances = src.amplitudeVariances;
      // The scale surface travels with the mean forest it was drawn beside: a
      // slot-sourced install takes that slot of the donor's saved variance
      // buffer, index-aligned with the mean slice above because one slot base
      // drives both buffers and both index by the sample number. A live-sourced
      // install takes the donor's live pair.
      if (slot < 0) {
        dst.varianceTrees = src.varianceTrees;
        dst.varianceTreeMasks = src.varianceTreeMasks;
        if (!scaleLeavesArePositive(dst.varianceTrees))
          return WarmStartResult::varianceMismatch;
      } else if (!src.varianceTrees.empty()) {
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
    }

    // containment (design "Containment"): a donor grown under a different (or
    // no) interaction constraint may hold a tree this sampler's constraint
    // forbids, and one under a split-variable restriction (BCF moderators, a
    // column-restricted variance forest) may split on a forbidden column;
    // either would be mis-scored, so refuse before touching live state. A no-op
    // when no chain is constrained. On the cross-grid path the donor's flat
    // trees resolve only against the donor grid, so install it over the store
    // for the scratch builds these feasibility checks make (ScopedCutGrid
    // restores the live grid on scope exit); the remap only ever collapses
    // splits, so a donor feasible pre-remap stays feasible after.
    auto checkContainment = [&]() -> WarmStartResult {
      for (size_t c = 0; c < chains_.size(); ++c)
        if (!chains_[c]->interactionStateFeasible(install[c]))
          return WarmStartResult::interactionMismatch;
      for (size_t c = 0; c < chains_.size(); ++c)
        if (!chains_[c]->columnMaskStateFeasible(install[c]))
          return WarmStartResult::columnMaskMismatch;
      return WarmStartResult::ok;
    };
    WarmStartResult containment;
    if (crossGrid) {
      ScopedCutGrid donorGrid(data_, donor.cutPoints);
      containment = checkContainment();
    } else {
      containment = checkContainment();
    }
    if (containment != WarmStartResult::ok) return containment;

    const std::vector<std::vector<double>>* donorGridPtr =
      crossGrid ? &donor.cutPoints : nullptr;
    for (size_t c = 0; c < chains_.size(); ++c)
      if (!chains_[c]->installForest(install[c], donorGridPtr, &data_))
        return WarmStartResult::shapeMismatch;
    // the variance half, separately so a refusal names the variance forest;
    // the shape gate above pairs hasVarianceForest with a non-empty block
    for (size_t c = 0; c < chains_.size(); ++c)
      if (chains_[c]->hasVarianceForest() &&
          !chains_[c]->installVarianceForest(install[c].varianceTrees,
                                             install[c].varianceTreeMasks,
                                             donorGridPtr, &data_))
        return WarmStartResult::varianceMismatch;
    // the store itself is left alone, so the draws it holds belong to the
    // donor's fit, not this one's: drop them rather than let a read replay
    // another sampler's posterior
    currentSampleNum_ = 0;
    recordedDraws_ = 0;
    return WarmStartResult::ok;
  }

  // Between-sample mutation, fanned out to every chain; new-vector lifetimes
  // are the caller's problem.
  void setOffset(const double* offset, bool updateScale) {
    for (auto& chain : chains_) chain->setOffset(offset, updateScale);
  }
  void setResponse(const double* y, bool updateScale) {
    for (auto& chain : chains_) chain->setResponse(y, updateScale);
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
  void reapplyWeights() {
    for (auto& chain : chains_) chain->reapplyWeights();
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
  void setTreeStorage(bool keepTrees, size_t numSamplesToStore) {
    size_t capacity =
      keepTrees ? (numSamplesToStore > 0 ? numSamplesToStore : 1) : 0;
    if (keepTrees == options_.keepTrees && capacity == savedTreeCapacity())
      return;
    options_.keepTrees = keepTrees;
    options_.numSamplesToStore = numSamplesToStore;
    for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
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

    std::vector<std::thread> workers;
    workers.reserve(numWorkers);
#ifndef _WIN32
    // spawn with SIGINT blocked so a Ctrl-C during the grow phase reaches only
    // the main thread, never a worker that has no R interrupt handler; the
    // main thread's mask is restored right after the spawn (mirrors run())
    sigset_t interruptSet, previousSet;
    sigemptyset(&interruptSet);
    sigaddset(&interruptSet, SIGINT);
    pthread_sigmask(SIG_BLOCK, &interruptSet, &previousSet);
#endif
    for (size_t w = 0; w < numWorkers; ++w) {
      workers.emplace_back([this, w, numWorkers, numChains, numSweeps]() {
        for (size_t c = w; c < numChains; c += numWorkers)
          chains_[c]->growForestFromRoot(numSweeps);
      });
    }
#ifndef _WIN32
    pthread_sigmask(SIG_SETMASK, &previousSet, nullptr);
#endif
    for (std::thread& worker : workers) worker.join();
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
  /// or empty collapses. Gaussian chains keep sigma and the variance prior
  /// fixed on the original scale.
  void setData(const double* x, const double* y, size_t numObservations,
               const double* weights, const double* offset,
               const double* x_test, size_t numTestObservations,
               const double* testOffset = nullptr) {
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

    data_.setData(x, numObservations);
    if (x_test != nullptr && numTestObservations > 0) {
      data_.buildTest(x_test, numTestObservations);
      data_.testOffset = testOffset;
    } else {
      data_.resetTestStorage();
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->applyNewData(y, weights, offset, oldCutPoints, params[c],
                               varianceParams[c]);
  }

  /// Replace the predictor matrix from a borrowed view (the store keeps its
  /// old values on failure). Unless forceUpdate, a leaf that would empty in
  /// any tree of any chain rolls the whole change back; forceUpdate instead
  /// collapses emptied leaves into their parents. The mutation kernels index
  /// the values as one column-major block, so the view must be a dense one
  /// (PredictorSource::isDenseBlock): a mapped or CSC-valued source has no
  /// kernel, and its holder materializes it before the transaction begins.
  PredictorUpdateResult setPredictor(const PredictorSource& newX,
                                     bool forceUpdate, bool updateCutPoints) {
    WholeMatrixUpdate strategy{data_, newX.denseValues};
    return runPredictorTransaction(strategy, forceUpdate, updateCutPoints);
  }

  /// Overwrite a subset of columns in place; the view holds the replacement
  /// block, column-major numObservations x numColumns, and columns names the
  /// store columns it fills. Same transaction and dense-view semantics as
  /// setPredictor.
  PredictorUpdateResult updatePredictor(const PredictorSource& newColumns,
                                        const size_t* columns,
                                        size_t numColumns, bool forceUpdate,
                                        bool updateCutPoints) {
    SubsetUpdate strategy{data_, newColumns.denseValues, columns, numColumns};
    return runPredictorTransaction(strategy, forceUpdate, updateCutPoints);
  }

  /// Dense convenience spellings: plain column-major blocks over the store's
  /// own row count.
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) {
    return setPredictor(densePredictorSource(newX, data_.numObservations,
                                             data_.numPredictors),
                        forceUpdate, updateCutPoints);
  }
  PredictorUpdateResult updatePredictor(const double* newColumns,
                                        const size_t* columns,
                                        size_t numColumns, bool forceUpdate,
                                        bool updateCutPoints) {
    return updatePredictor(
      densePredictorSource(newColumns, data_.numObservations, numColumns),
      columns, numColumns, forceUpdate, updateCutPoints);
  }

  /// Install externally chosen cut points (ascending) for a subset of
  /// columns and unconditionally refresh the trees: splits that fall out of
  /// range or lose their observations collapse into their parents, exactly
  /// as a forced predictor update does.
  /// currentPredictors is the call-time predictor matrix the columns
  /// re-quantize from (data@x); null for CSC/mixed stores, whose columns read
  /// their retained slices instead.
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints, const size_t* columns,
                    size_t numColumns, const double* currentPredictors) {
    for (size_t k = 0; k < numColumns; ++k)
      data_.setCutPointsForColumn(columns[k], newCutPoints[k], numCutPoints[k],
                                  currentPredictors);
    for (auto& chain : chains_) chain->forceRefreshTrees();
  }

  /// Install one column's new values observation-by-observation in random
  /// scan order, declining exactly those whose move would empty a leaf in any
  /// tree of any forest of any chain; installed must have room for a flag per
  /// observation. Returns finalize() validity, which the guard makes true by
  /// construction.
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
  /// Grouped random intercepts: the group count, 0 when ungrouped.
  size_t numGroups() const { return options_.numGroups; }
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
  // single-forest samplers, where numForests() is 1 and bcfGlue reports none.
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
  /// Whether the response family carries a dispersion r; chain 0 answers for
  /// all, as every chain carries the same family.
  bool carriesDispersion() const { return chains_[0]->carriesDispersion(); }
  /// Chain chainNum's dispersion r in force, 0 off a family carrying one.
  double dispersion(size_t chainNum) const {
    return chains_[chainNum]->dispersion();
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
  /// Restates forest forestIndex's leaf prior on EVERY chain; false, writing
  /// nothing, on a refusal. Every chain refuses on the same two conditions, so
  /// the fan-out cannot land half applied - and every chain skips a write
  /// reproducing its own current value independently, so a read-then-write on
  /// diverged chains is inert on each.
  bool setForestPriorScale(size_t forestIndex, double priorScale) {
    bool written = true;
    for (auto& chain : chains_)
      written = chain->setForestPriorScale(forestIndex, priorScale) && written;
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
    options_.groupIndices = nullptr;    // borrowed; consumed by the chains
    options_.survivalStatus = nullptr;  // borrowed; consumed by the chains

    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
  }

  /// Two-phase transaction over every chain: validate all trees of all
  /// forests of all chains first, then rebuild fits only if everything holds,
  /// so a failure in a late chain never leaves an early chain's fits
  /// overwritten. columns names the touched columns (null for a whole-matrix
  /// swap), which prunes the trees of forests past the first that cannot be
  /// affected; the survivor lists travel from phase one to phase two inside
  /// the handoff, so both phases quantify over the same set.
  bool revalidateAllChains(const size_t* columns, size_t numColumns) {
    size_t numChains = chains_.size();
    std::vector<typename Chain<L, ResidT>::ForestRevalidation> state(numChains);

    bool allValid = true;
    for (size_t c = 0; c < numChains && allValid; ++c)
      allValid = chains_[c]->revalidateTrees(state[c], columns, numColumns);
    if (!allValid) return false;

    for (size_t c = 0; c < numChains; ++c)
      chains_[c]->rebuildFitsFromParameters(state[c]);
    return true;
  }

  /// The two ways a predictor-replacement transaction moves its columns: the
  /// column list each spans (a null columns pointer means every predictor in
  /// order, over column-major values) and the snapshot/apply/restore triple
  /// runPredictorTransaction drives. WholeMatrixUpdate swaps the whole code
  /// storage aside; SubsetUpdate journals a named subset column by column.
  struct WholeMatrixUpdate {
    ColumnStore& data;
    const double* values;
    const size_t* columns = nullptr;
    std::vector<xint_t> oldCodes;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
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
    void applyForced(bool updateCuts) { data.setPredictors(values, updateCuts); }
    void snapshotApply(bool updateCuts) {
      // move the live codes aside and rebuild into fresh storage: a reject swaps
      // them back and an accept drops them, so no whole-matrix copy survives
      oldCodes = std::move(data.train.codes);
      data.train.codes.assign(oldCodes.size(), 0);
      oldHasMissing = data.hasMissing;
      if (updateCuts) oldCuts = data.cutPoints;
      csc = data.builtFromCsc;
      if (csc) {
        oldSparseColumns = data.train.sparseColumns;
        oldSources = data.train.sources;
        oldOwnedRows = data.ownedCscRows;
        oldOwnedValues = data.ownedCscValues;
        oldCscOwned = data.cscColumnOwned;
      }
      data.setPredictors(values, updateCuts);
    }
    void restore(bool updateCuts) {
      data.train.codes = std::move(oldCodes);
      data.hasMissing = std::move(oldHasMissing);
      if (updateCuts) data.cutPoints = std::move(oldCuts);
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
    const double* values;
    const size_t* columns;
    size_t count;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
    std::vector<ColumnStore::ColumnCodeRollback> records;
    // CSC-backed columns of the subset snapshot their sparse storage here
    // instead of the dense per-cell journal (which has no codes[] to journal
    // for a rank column)
    std::vector<ColumnStore::CscColumnRollback> cscRecords;

    size_t numColumns() const { return count; }
    void applyForced(bool updateCuts) {
      data.setColumns(values, columns, count, updateCuts);
    }
    void snapshotApply(bool updateCuts) {
      // journal each touched column's changed cells (past a quarter changed the
      // record falls back to a whole pre-change copy); the small hasMissing and
      // cut pieces snapshot per column. A CSC-backed column instead snapshots
      // its rank/densified storage and owned slice and rebuilds from the dense
      // column.
      size_t n = data.numObservations;
      oldHasMissing.resize(count);
      oldCuts.resize(updateCuts ? count : 0);
      records.resize(count);
      cscRecords.resize(count);
      for (size_t k = 0; k < count; ++k) {
        size_t j = columns[k];
        oldHasMissing[k] = data.hasMissing[j];
        if (updateCuts) oldCuts[k] = data.cutPoints[j];
        if (data.columnIsCscBacked(j)) {
          data.snapshotCscColumn(j, cscRecords[k]);
          data.mutateCscColumnFromDense(j, values + k * n, updateCuts);
        } else {
          data.setColumnJournaled(j, values + k * n, updateCuts, n / 4,
                                  records[k]);
        }
      }
    }
    void restore(bool updateCuts) {
      for (size_t k = 0; k < count; ++k) {
        size_t j = columns[k];
        if (data.columnIsCscBacked(j))
          data.restoreCscColumn(j, cscRecords[k]);
        else
          data.restoreColumn(j, records[k]);
        data.hasMissing[j] = oldHasMissing[k];
        if (updateCuts) data.cutPoints[j] = std::move(oldCuts[k]);
      }
    }
  };

  /// One predictor-replacement transaction over a strategy's column list:
  /// precheck that the cuts stay representable, then either collapse emptied
  /// leaves under forceUpdate or snapshot-apply and keep the change when
  /// revalidateAllChains holds, else restore the snapshot and repartition. The
  /// engine keeps no predictor matrix, so the raw a reject must put back is
  /// snapshotted here: the gathered leaf-covariate copies and, for a mixed
  /// store, the owned dense block of the columns the strategy touches (the
  /// strategy's own records carry no raw). The strategy owns the codes, missing
  /// flags, and cut grids it moves.
  template <typename Strategy>
  PredictorUpdateResult runPredictorTransaction(Strategy& strategy,
                                                bool forceUpdate,
                                                bool updateCutPoints) {
    for (size_t k = 0; k < strategy.numColumns(); ++k) {
      size_t j = strategy.columns ? strategy.columns[k] : k;
      // a factor column of either kind must hold representable level codes
      // whether or not cut points refresh: its grid is its level table, so a
      // value outside the table has no position on it
      if ((updateCutPoints || data_.isFactor(j)) &&
          !data_.cutsWouldRemainValid(
            j, strategy.values + k * data_.numObservations))
        return PredictorUpdateResult::invalidCutPoints;
    }

    if (forceUpdate) {
      strategy.applyForced(updateCutPoints);
      for (auto& chain : chains_) chain->forceRefreshTrees();
      requantizeTestColumns(strategy, updateCutPoints);
      return PredictorUpdateResult::accepted;
    }

    std::vector<double> oldGatheredRaw = data_.gatheredRawValues;
    ColumnStore::OwnedDenseRollback oldOwnedDense;
    data_.snapshotOwnedDenseColumns(strategy.columns, strategy.numColumns(),
                                    oldOwnedDense);
    strategy.snapshotApply(updateCutPoints);
    if (!revalidateAllChains(strategy.columns, strategy.numColumns())) {
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

    bool observationWouldRemainValid(size_t i) override {
      size_t n = sampler_.data_.numObservations;
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
      return valid;
    }

    void commitObservation(size_t i) override {
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
    /// predicate with no exemption, so the two sets coincide there.
    bool finalize() override {
      return sampler_.revalidateAllChains(&column_, 1);
    }

  private:
    /// One cached tree, addressed as the engine holds it. The variance forest
    /// has no index in forests_, so it takes a sentinel in the forest slot.
    static constexpr size_t varianceForest = static_cast<size_t>(-1);
    struct CachedTree {
      size_t chain, forest, tree;
    };

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
  };

  SamplerOptions options_;
  ResponseFamily family_;
  ColumnStore data_;
  std::vector<std::unique_ptr<Chain<L, ResidT>>> chains_;
  size_t currentSampleNum_ = 0;  // next saved-tree slot, wrapping circularly
  size_t recordedDraws_ = 0;     // slots written since the last reset, capped
                                 // at capacity; the extent of every read
  double runningTime_ = 0.0;     // seconds accumulated across runs
};

using ConstantLeafSampler = Sampler<ConstantGaussianLeaf>;

}  // namespace bartcore

#endif  // BARTCORE_SAMPLER_HPP
