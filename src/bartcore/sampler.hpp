#ifndef BARTCORE_SAMPLER_HPP
#define BARTCORE_SAMPLER_HPP

#include <atomic>
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
/// points (setCutPoints may have replaced the ones creation induces), and
/// the saved-tree write position.
struct SamplerStateData {
  std::vector<ChainStateData> chains;
  std::vector<std::vector<double>> cutPoints;  // empty vector per categorical column
  size_t currentSampleNum = 0;
};

/// Why a warm start (installForests) refused; ok on success. A single donor
/// forest can seed several chains, so the donor's chain count need not match.
enum class WarmStartResult { ok, shapeMismatch, gridMismatch, dartMismatch };

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
template <IntegrableLeafModel L>
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
    if (options.columnSources != nullptr) {
      data_.buildMixed(options.mixedDenseValues, options.cscColumnPointers,
                       options.cscRowIndices, options.cscValues,
                       options.columnSources, numObservations, numPredictors,
                       options.maxNumCutsPerVariable, options.maxNumCuts,
                       options.useQuantiles, options.columnTypes);
    } else if (options.cscColumnPointers != nullptr) {
      data_.buildFromCsc(options.cscColumnPointers, options.cscRowIndices,
                         options.cscValues, numObservations, numPredictors,
                         options.maxNumCutsPerVariable, options.maxNumCuts,
                         options.useQuantiles);
    } else if (options.maxNumCutsPerVariable != nullptr) {
      data_.build(x, numObservations, numPredictors,
                  options.maxNumCutsPerVariable, options.useQuantiles,
                  options.columnTypes);
    } else {
      data_.build(x, numObservations, numPredictors, options.maxNumCuts,
                  options.useQuantiles, options.columnTypes);
    }
    options_.maxNumCutsPerVariable = nullptr;  // borrowed; consumed by build
    options_.columnTypes = nullptr;
    // the CSC slices and dense-source pointers live on in the store
    options_.cscColumnPointers = nullptr;
    options_.cscRowIndices = nullptr;
    options_.cscValues = nullptr;
    options_.mixedDenseValues = nullptr;
    options_.columnSources = nullptr;

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
    options_.columnTypes = nullptr;
    options_.useQuantiles = data_.useQuantiles;

    initializeChains(y, weights, offset, sigmaEstimate, sigmaDf,
                     sigmaRawScale, rngs);
  }

  /// A BCF two-forest sampler over dense predictors (docs/design/bcf.md):
  /// gaussian only, one prognostic and one treatment forest per chain. The
  /// CSC/mixed and view ingestion paths are not offered here.
  Sampler(const double* x, const double* y, size_t numObservations,
          size_t numPredictors, const double* weights, const double* offset,
          double sigmaEstimate, double sigmaDf, double sigmaRawScale,
          const SamplerOptions& options, const BCFSpec& spec,
          ext_rng* const* rngs)
    : options_(options), family_(ResponseFamily::gaussian) {
    if (options.maxNumCutsPerVariable != nullptr)
      data_.build(x, numObservations, numPredictors,
                  options.maxNumCutsPerVariable, options.useQuantiles,
                  options.columnTypes);
    else
      data_.build(x, numObservations, numPredictors, options.maxNumCuts,
                  options.useQuantiles, options.columnTypes);
    options_.maxNumCutsPerVariable = nullptr;
    options_.columnTypes = nullptr;
    // single-forest queries (numTrees, savedTree, printTrees) address the
    // prognostic forest
    options_.numTrees = spec.mu.numTrees;

    chains_.reserve(options_.numChains);
    for (size_t c = 0; c < options_.numChains; ++c)
      chains_.push_back(std::make_unique<Chain<L>>(
        data_, y, weights, offset, sigmaEstimate, sigmaDf, sigmaRawScale,
        options_, spec, rngs[c]));
    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
  }

  // chains reference the store member, so the sampler's address is pinned
  Sampler(const Sampler&) = delete;
  Sampler& operator=(const Sampler&) = delete;

  /// Replace the test predictors, keeping any test offset: the caller
  /// guarantees the row count still matches it (the bridge refuses
  /// otherwise). Passing a new offset too goes through setTestOffset.
  void setTestPredictors(const double* x_test, size_t numTestObservations) {
    const double* testOffset = data_.testOffset;
    data_.buildTest(x_test, numTestObservations);
    data_.testOffset = testOffset;
    for (auto& chain : chains_) chain->resizeTestStorage();
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
    std::vector<Results> chainResults(numChains);
    for (size_t c = 0; c < numChains; ++c) {
      Results& r(chainResults[c]);
      if (results.sigma != nullptr) r.sigma = results.sigma + c * numSamples;
      if (results.k != nullptr) r.k = results.k + c * numSamples;
      if (results.trainingFits != nullptr)
        r.trainingFits =
          results.trainingFits + c * numSamples * data_.numObservations;
      if (results.testFits != nullptr)
        r.testFits =
          results.testFits + c * numSamples * data_.numTestObservations;
      if (results.variableCounts != nullptr)
        r.variableCounts =
          results.variableCounts + c * numSamples * data_.numPredictors;
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
    if (capacity > 0 && numSamples > 0)
      currentSampleNum_ = (currentSampleNum_ + numSamples) % capacity;

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
      // prefix every chain's line unconditionally
      ext_printf("[%lu] ", static_cast<unsigned long>(c + 1));
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
  void setCurrentSampleNum(size_t sampleNum) { currentSampleNum_ = sampleNum; }

  const std::vector<FlatNode>& savedTree(size_t chainNum, size_t slot,
                                         size_t treeNum) const {
    return chains_[chainNum]->savedTree(slot, treeNum);
  }
  const std::vector<double>& savedTreeSlopes(size_t chainNum, size_t slot,
                                             size_t treeNum) const {
    return chains_[chainNum]->savedTreeSlopes(slot, treeNum);
  }
  const std::vector<std::uint64_t>& savedTreeMasks(size_t chainNum,
                                                   size_t slot,
                                                   size_t treeNum) const {
    return chains_[chainNum]->savedTreeMasks(slot, treeNum);
  }
  /// slopes, when non-null, receives a vector-parameter tree's slopes
  /// (numParams - 1 per leaf, pre-order); cleared for scalar leaf models.
  /// masks, when non-null, receives the wide categorical side channel.
  void flattenTree(size_t chainNum, size_t treeNum,
                   std::vector<FlatNode>& nodes,
                   std::vector<std::uint32_t>& counts,
                   std::vector<double>* slopes = nullptr,
                   std::vector<std::uint64_t>* masks = nullptr) {
    chains_[chainNum]->flattenTree(treeNum, nodes, counts, slopes, masks);
  }

  /// Fits for raw column-major test rows, on the original response scale
  /// (offsets are the caller's problem). With saved trees, out is
  /// numTestObservations x savedTreeCapacity x numChains, chain-major like
  /// Results; without, one slab per chain from the live trees.
  void predict(const double* x_test, size_t numTestObservations, double* out) {
    size_t capacity = savedTreeCapacity();
    for (size_t c = 0; c < chains_.size(); ++c) {
      if (capacity > 0) {
        for (size_t slot = 0; slot < capacity; ++slot)
          chains_[c]->predictFromSavedSample(
            slot, x_test, numTestObservations,
            out + (c * capacity + slot) * numTestObservations);
      } else {
        chains_[c]->predictFromCurrentTrees(x_test, numTestObservations,
                                            out + c * numTestObservations);
      }
    }
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
  }

  bool setState(const SamplerStateData& state) {
    if (state.chains.size() != chains_.size()) return false;
    if (state.cutPoints.size() != data_.numPredictors) return false;
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.types[j] == ColumnType::categorical) {
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
    std::vector<xint_t> oldCodes(data_.codes);
    std::vector<xint_t> oldTestCodes(data_.testCodes);
    // rank columns re-quantize into their own storage, not codes
    std::vector<SparseColumnData> oldSparseColumns(data_.sparseColumns);
    for (size_t j = 0; j < data_.numPredictors; ++j) {
      if (data_.types[j] == ColumnType::categorical) continue;
      data_.setCutPointsForColumn(j, state.cutPoints[j].data(),
                                  static_cast<std::uint32_t>(
                                    state.cutPoints[j].size()));
    }

    bool allValid = true;
    for (size_t c = 0; c < chains_.size() && allValid; ++c)
      allValid = chains_[c]->stateIsValid(state.chains[c]);

    if (!allValid) {
      data_.cutPoints = std::move(oldCutPoints);
      data_.numCuts = std::move(oldNumCuts);
      data_.maxNumCuts = std::move(oldMaxNumCuts);
      data_.codes = std::move(oldCodes);
      data_.testCodes = std::move(oldTestCodes);
      data_.sparseColumns = std::move(oldSparseColumns);
      return false;
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      if (!chains_[c]->setState(state.chains[c])) return false;
    size_t capacity = savedTreeCapacity();
    currentSampleNum_ = capacity > 0 ? state.currentSampleNum % capacity : 0;
    return true;
  }

  /// Warm start: seed each chain's live forest from a donor sample. sampleMap
  /// selects, for destination chain c, a donor (chain, slot); slot < 0 takes
  /// the donor chain's live trees, else its saved slot (slot-major, one forest
  /// per slot). Only trees, sigma, k, and DART transfer - rng and auxiliary
  /// state stay fresh, so each chain evolves independently from its own
  /// stream. The donor must share this sampler's cut grid, per-forest tree
  /// counts, and DART mode; on any mismatch nothing is touched.
  WarmStartResult installForests(
      const SamplerStateData& donor,
      const std::vector<std::pair<size_t, int>>& sampleMap) {
    if (sampleMap.size() != chains_.size())
      return WarmStartResult::shapeMismatch;
    if (donor.cutPoints != data_.cutPoints)
      return WarmStartResult::gridMismatch;

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
      }
      dst.sigma = src.sigma;
      dst.fitMin = src.fitMin;
      dst.fitMax = src.fitMax;
      dst.dartProbabilities = src.dartProbabilities;
      dst.dartAlpha = src.dartAlpha;
      dst.dartNumUpdatesSkipped = src.dartNumUpdatesSkipped;
      dst.hasBCF = src.hasBCF;
      dst.a = src.a;
      dst.aVariance = src.aVariance;
      dst.b0 = src.b0;
      dst.b1 = src.b1;
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      if (!chains_[c]->installForest(install[c]))
        return WarmStartResult::shapeMismatch;
    currentSampleNum_ = 0;
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
  /// Case weights, gaussian only (the host rejects binary families): a bare
  /// pointer swap, entering the next iteration's node statistics and sigma
  /// draw with nothing rescaled.
  void setWeights(const double* weights) {
    for (auto& chain : chains_) chain->setWeights(weights);
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
  void setVerbose(bool verbose, std::uint32_t printEvery) {
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

  /// Info dump; the per-node output format is R-visible and pinned by tests.
  /// Without keepTrees the live trees print and sample indices are ignored;
  /// with it, the requested saved slots print in the saved-tree format.
  void printTrees(const size_t* chainIndices, size_t numChainIndices,
                  const size_t* sampleIndices, size_t numSampleIndices,
                  const size_t* treeIndices, size_t numTreeIndices) {
    int indent = 0;
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      if (numChainIndices > 1) {
        ext_printf("chain %lu\n", static_cast<unsigned long>(chainNum + 1));
        indent += 2;
      }
      if (!options_.keepTrees) {
        for (size_t k = 0; k < numTreeIndices; ++k)
          chains_[chainNum]->printTree(treeIndices[k], indent);
      } else {
        for (size_t j = 0; j < numSampleIndices; ++j) {
          size_t sampleNum = sampleIndices[j];
          if (numSampleIndices > 1) {
            ext_printf("%*ssample %lu\n", indent, "",
                       static_cast<unsigned long>(sampleNum + 1));
            indent += 2;
          }
          for (size_t k = 0; k < numTreeIndices; ++k)
            chains_[chainNum]->printSavedTree(sampleNum, treeIndices[k], indent);
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
    // moves; the old cut values drive the split remap
    std::vector<typename Chain<L>::TreeParameters> params(chains_.size());
    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->recoverTreeParameters(params[c]);

    std::vector<std::vector<double>> oldCutPoints(data_.cutPoints);

    data_.setData(x, numObservations);
    if (x_test != nullptr && numTestObservations > 0) {
      data_.buildTest(x_test, numTestObservations);
      data_.testOffset = testOffset;
    } else {
      data_.clearTest();
    }

    for (size_t c = 0; c < chains_.size(); ++c)
      chains_[c]->applyNewData(y, weights, offset, oldCutPoints, params[c]);
  }

  /// Replace the predictor matrix (borrowed, column-major; the old pointer is
  /// kept on failure). Unless forceUpdate, a leaf that would empty in any
  /// tree of any chain rolls the whole change back; forceUpdate instead
  /// collapses emptied leaves into their parents.
  PredictorUpdateResult setPredictor(const double* newX, bool forceUpdate,
                                     bool updateCutPoints) {
    size_t n = data_.numObservations;
    // categorical columns must hold representable codes whether or not cut
    // points refresh
    for (size_t j = 0; j < data_.numPredictors; ++j)
      if ((updateCutPoints ||
           data_.types[j] == ColumnType::categorical) &&
          !data_.cutsWouldRemainValid(j, newX + j * n))
        return PredictorUpdateResult::invalidCutPoints;

    const double* oldX = data_.x;
    std::vector<xint_t> oldCodes;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
    if (!forceUpdate) {
      oldCodes = data_.codes;
      // re-quantizing recomputes hasMissing; a rollback must restore it too so
      // rules stay consistent with the reachable gauge
      oldHasMissing = data_.hasMissing;
      if (updateCutPoints) oldCuts = data_.cutPoints;
    }

    data_.setPredictors(newX, updateCutPoints);

    if (forceUpdate) {
      for (auto& chain : chains_) chain->forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t j = 0; j < data_.numPredictors; ++j)
          data_.quantizeTestColumn(j);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateAllChains()) {
      data_.x = oldX;
      data_.codes = std::move(oldCodes);
      data_.hasMissing = std::move(oldHasMissing);
      if (updateCutPoints) data_.cutPoints = std::move(oldCuts);
      for (auto& chain : chains_) chain->repartitionTrees();
      return PredictorUpdateResult::rolledBack;
    }

    if (updateCutPoints && data_.numTestObservations > 0)
      for (size_t j = 0; j < data_.numPredictors; ++j)
        data_.quantizeTestColumn(j);
    return PredictorUpdateResult::accepted;
  }

  /// Overwrite a subset of columns in place; newColumns is column-major,
  /// numObservations x numColumns. Same transaction semantics as
  /// setPredictor.
  PredictorUpdateResult updatePredictor(const double* newColumns,
                                        const size_t* columns,
                                        size_t numColumns, bool forceUpdate,
                                        bool updateCutPoints) {
    size_t n = data_.numObservations;
    for (size_t k = 0; k < numColumns; ++k)
      if ((updateCutPoints ||
           data_.types[columns[k]] == ColumnType::categorical) &&
          !data_.cutsWouldRemainValid(columns[k], newColumns + k * n))
        return PredictorUpdateResult::invalidCutPoints;

    std::vector<double> oldValues;
    std::vector<xint_t> oldCodes;
    std::vector<std::uint8_t> oldHasMissing;
    std::vector<std::vector<double>> oldCuts;
    if (!forceUpdate) {
      oldValues.resize(n * numColumns);
      oldCodes.resize(n * numColumns);
      oldHasMissing.resize(numColumns);
      if (updateCutPoints) oldCuts.resize(numColumns);
      for (size_t k = 0; k < numColumns; ++k) {
        std::memcpy(oldValues.data() + k * n, data_.x + columns[k] * n,
                    n * sizeof(double));
        std::memcpy(oldCodes.data() + k * n,
                    data_.codes.data() + data_.codeOffsets[columns[k]],
                    n * sizeof(xint_t));
        // re-quantizing rebuilds hasMissing; a rollback must restore it too so
        // rules stay consistent with the reachable gauge
        oldHasMissing[k] = data_.hasMissing[columns[k]];
        if (updateCutPoints) oldCuts[k] = data_.cutPoints[columns[k]];
      }
    }

    data_.setColumns(newColumns, columns, numColumns, updateCutPoints);

    if (forceUpdate) {
      for (auto& chain : chains_) chain->forceRefreshTrees();
      if (updateCutPoints && data_.numTestObservations > 0)
        for (size_t k = 0; k < numColumns; ++k)
          data_.quantizeTestColumn(columns[k]);
      return PredictorUpdateResult::accepted;
    }

    if (!revalidateAllChains()) {
      double* x_mutable = const_cast<double*>(data_.x);
      for (size_t k = 0; k < numColumns; ++k) {
        size_t j = columns[k];
        std::memcpy(x_mutable + j * n, oldValues.data() + k * n,
                    n * sizeof(double));
        std::memcpy(data_.codes.data() + data_.codeOffsets[j],
                    oldCodes.data() + k * n, n * sizeof(xint_t));
        data_.hasMissing[j] = oldHasMissing[k];
        if (updateCutPoints) data_.cutPoints[j] = std::move(oldCuts[k]);
      }
      for (auto& chain : chains_) chain->repartitionTrees();
      return PredictorUpdateResult::rolledBack;
    }

    if (updateCutPoints && data_.numTestObservations > 0)
      for (size_t k = 0; k < numColumns; ++k)
        data_.quantizeTestColumn(columns[k]);
    return PredictorUpdateResult::accepted;
  }

  /// Install externally chosen cut points (ascending) for a subset of
  /// columns and unconditionally refresh the trees: splits that fall out of
  /// range or lose their observations collapse into their parents, exactly
  /// as a forced predictor update does.
  void setCutPoints(const double* const* newCutPoints,
                    const std::uint32_t* numCutPoints, const size_t* columns,
                    size_t numColumns) {
    for (size_t k = 0; k < numColumns; ++k)
      data_.setCutPointsForColumn(columns[k], newCutPoints[k],
                                  numCutPoints[k]);
    for (auto& chain : chains_) chain->forceRefreshTrees();
  }

  /// Install one column's new values observation-by-observation in random
  /// scan order, rolling back exactly those whose move would empty a leaf in
  /// any tree of any chain; installed must have room for a flag per
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

  ResponseFamily family() const { return family_; }
  /// Grouped random intercepts: the group count, 0 when ungrouped.
  size_t numGroups() const { return options_.numGroups; }
  double sigma(size_t chainNum = 0) const { return chains_[chainNum]->sigma(); }
  double k(size_t chainNum = 0) const { return chains_[chainNum]->k(); }
  bool kIsSampled() const { return options_.updateK; }
  bool usesDart() const { return options_.useDart; }
  const ColumnStore& data() const { return data_; }
  const Chain<L>& chain(size_t chainNum) const { return *chains_[chainNum]; }
  Chain<L>& chain(size_t chainNum) { return *chains_[chainNum]; }
  size_t numChains() const { return chains_.size(); }
  size_t numThreads() const { return options_.numThreads; }

  // BCF surface, fanned to every chain (docs/design/bcf.md); benign on
  // single-forest samplers, where numForests() is 1 and bcfGlue reports none.
  size_t numForests() const { return chains_[0]->numForests(); }
  void setTreatment(const double* z) {
    for (auto& chain : chains_) chain->setTreatment(z);
  }
  bool bcfGlue(size_t chainNum, double* out) const {
    return chains_[chainNum]->bcfGlue(out[0], out[1], out[2]);
  }
  void forestTotalFits(size_t chainNum, size_t forestIndex, double* out) const {
    chains_[chainNum]->forestTotalFits(forestIndex, out);
  }
  size_t numTrees() const { return options_.numTrees; }
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
      chains_.push_back(std::make_unique<Chain<L>>(
        data_, y, weights, offset, family_, sigmaEstimate, sigmaDf,
        sigmaRawScale, options_, rngs[c]));
    options_.groupIndices = nullptr;  // borrowed; consumed by the chains

    if (options_.keepTrees) {
      size_t capacity =
        options_.numSamplesToStore > 0 ? options_.numSamplesToStore : 1;
      for (auto& chain : chains_) chain->initializeSavedTrees(capacity);
    }
  }

  /// Two-phase transaction over every chain: validate all trees of all
  /// chains first, then rebuild fits only if everything holds, so a failure
  /// in a late chain never leaves an early chain's fits overwritten.
  bool revalidateAllChains() {
    size_t numChains = chains_.size();
    std::vector<typename Chain<L>::TreeParameters> params(numChains);

    bool allValid = true;
    for (size_t c = 0; c < numChains && allValid; ++c)
      allValid = chains_[c]->revalidateTrees(params[c]);
    if (!allValid) return false;

    for (size_t c = 0; c < numChains; ++c)
      chains_[c]->rebuildFitsFromParameters(params[c]);
    return true;
  }

  /// Caches each observation's leaf and per-leaf occupancy for every tree of
  /// every chain by routing codes through the split structure, so staging a
  /// move is a descent and a count check per tree.
  class UpdateSessionImpl final : public PredictorUpdateSession {
  public:
    UpdateSessionImpl(Sampler& sampler, const double* newColumn, size_t column)
      : sampler_(sampler), column_(column), newColumn_(newColumn) {
      size_t n = sampler_.data_.numObservations;
      size_t totalNumTrees =
        sampler_.options_.numTrees * sampler_.chains_.size();

      newCodes_.resize(n);
      for (size_t i = 0; i < n; ++i)
        newCodes_[i] = sampler_.data_.codeFor(column_, newColumn_[i]);

      observationLeaf_.resize(totalNumTrees * n);
      leafCounts_.resize(totalNumTrees);
      pendingNewLeaf_.assign(totalNumTrees, invalidNode);
      pendingOldLeaf_.resize(totalNumTrees);

      for (size_t t = 0; t < totalNumTrees; ++t) {
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
      ColumnStore& data(sampler_.data_);
      size_t n = data.numObservations;
      const_cast<double*>(data.x)[i + column_ * n] = newColumn_[i];
      data.codes[data.codeOffsets[column_] + i] = newCodes_[i];
      for (size_t t = 0; t < leafCounts_.size(); ++t) {
        if (pendingNewLeaf_[t] != invalidNode) {
          --leafCounts_[t][static_cast<size_t>(pendingOldLeaf_[t])];
          ++leafCounts_[t][static_cast<size_t>(pendingNewLeaf_[t])];
        }
      }
    }

    bool finalize() override { return sampler_.revalidateAllChains(); }

  private:
    const Tree& treeAt(size_t t) const {
      size_t numTrees = sampler_.options_.numTrees;
      return sampler_.chains_[t / numTrees]->tree(t % numTrees);
    }

    Sampler& sampler_;
    size_t column_;
    const double* newColumn_;
    std::vector<xint_t> newCodes_;
    std::vector<int32_t> observationLeaf_;  // totalNumTrees x numObservations
    std::vector<std::vector<std::uint32_t>> leafCounts_;  // arena-id indexed
    std::vector<int32_t> pendingNewLeaf_;  // invalidNode when no move staged
    std::vector<int32_t> pendingOldLeaf_;
  };

  SamplerOptions options_;
  ResponseFamily family_;
  ColumnStore data_;
  std::vector<std::unique_ptr<Chain<L>>> chains_;
  size_t currentSampleNum_ = 0;  // next saved-tree slot, wrapping circularly
  double runningTime_ = 0.0;     // seconds accumulated across runs
};

using ConstantLeafSampler = Sampler<ConstantGaussianLeaf>;

}  // namespace bartcore

#endif  // BARTCORE_SAMPLER_HPP
