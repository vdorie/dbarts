// Allocation failures inside a sweep. The global allocation functions are
// replaced for the whole binary: unarmed they count, armed they fail the N-th
// allocation once. bartcore is header-only, so every engine vector allocates
// through them. Each arm builds a fresh sampler for every N - scratch
// capacity survives setState, so a reused one would no longer make the
// allocations it was counted making - installs one chain state, fails the
// N-th allocation of the step for every N that step makes, and checks that
// the chain is left whole: the totals those of the trees, a state that reads
// and installs, a next sweep that is finite and consistent, and after it the
// fits the trees route the rows to.

#include "common.hpp"

#include <cstdlib>
#include <initializer_list>
#include <new>

namespace {

std::atomic<bool> allocationCounting{false};
std::atomic<long> allocationCount{0};
std::atomic<bool> allocationArmed{false};
std::atomic<long> allocationCountdown{0};

/// Null where the armed countdown reaches this allocation.
void* allocateOrFail(std::size_t size, std::size_t alignment) {
  if (allocationCounting.load(std::memory_order_relaxed))
    allocationCount.fetch_add(1, std::memory_order_relaxed);
  if (allocationArmed.load(std::memory_order_relaxed) &&
      allocationCountdown.fetch_sub(1, std::memory_order_relaxed) == 1) {
    allocationArmed.store(false, std::memory_order_relaxed);
    return nullptr;
  }
  if (size == 0) size = 1;
  if (alignment <= alignof(std::max_align_t)) return std::malloc(size);
  void* p = nullptr;
  return posix_memalign(&p, alignment, size) == 0 ? p : nullptr;
}

void* allocateOrThrow(std::size_t size, std::size_t alignment) {
  void* p = allocateOrFail(size, alignment);
  if (p == nullptr) throw std::bad_alloc();
  return p;
}

}  // namespace

void* operator new(std::size_t size) {
  return allocateOrThrow(size, alignof(std::max_align_t));
}
void* operator new[](std::size_t size) {
  return allocateOrThrow(size, alignof(std::max_align_t));
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocateOrThrow(size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocateOrThrow(size, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  return allocateOrFail(size, alignof(std::max_align_t));
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  return allocateOrFail(size, alignof(std::max_align_t));
}
void* operator new(std::size_t size, std::align_val_t alignment,
                   const std::nothrow_t&) noexcept {
  return allocateOrFail(size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
  return allocateOrFail(size, static_cast<std::size_t>(alignment));
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
void operator delete(void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept {
  std::free(p);
}
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept {
  std::free(p);
}
void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept {
  std::free(p);
}
void operator delete(void* p, std::align_val_t,
                     const std::nothrow_t&) noexcept {
  std::free(p);
}
void operator delete[](void* p, std::align_val_t,
                       const std::nothrow_t&) noexcept {
  std::free(p);
}

namespace {

/// The allocations \p body makes.
template <typename F>
long countAllocations(F&& body) {
  allocationCount.store(0, std::memory_order_relaxed);
  allocationCounting.store(true, std::memory_order_relaxed);
  body();
  allocationCounting.store(false, std::memory_order_relaxed);
  return allocationCount.load(std::memory_order_relaxed);
}

/// Whether \p body threw an allocation failure with its N-th allocation
/// failed; reached says whether it made an N-th. Any other exception is a
/// failure of the test.
template <typename F>
bool failsAt(long n, F&& body, bool& otherException, bool& reached) {
  allocationCountdown.store(n, std::memory_order_relaxed);
  allocationArmed.store(true, std::memory_order_relaxed);
  bool threw = false;
  try {
    body();
  } catch (const std::bad_alloc&) {
    threw = true;
  } catch (...) {
    otherException = true;
  }
  reached = !allocationArmed.load(std::memory_order_relaxed);
  allocationArmed.store(false, std::memory_order_relaxed);
  return threw;
}

/// Every chain's every forest's totalFits against its trees' fits summed in
/// tree order, and every value finite.
template <typename S>
bool totalsAreTheTrees(S& sampler, double tolerance) {
  std::size_t n = sampler.data().numObservations;
  std::size_t numForests = sampler.numForests();
  for (std::size_t c = 0; c < sampler.numChains(); ++c) {
    auto& chain = sampler.chain(c);
    for (std::size_t f = 0; f < numForests; ++f) {
      std::size_t numTrees = chain.numTreesInForest(f);
      std::vector<double> fits(n * numTrees);
      TestPeer::forestTreeFits(chain, f, fits.data());
      const std::vector<double>& total = TestPeer::totalFitsInForest(chain, f);
      for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t t = 0; t < numTrees; ++t) sum += fits[t * n + i];
        if (!std::isfinite(total[i]) || !std::isfinite(sum) ||
            std::fabs(total[i] - sum) > tolerance * (1.0 + std::fabs(sum)))
          return false;
      }
    }
  }
  return true;
}

/// Whether every column that splits by threshold holds as many points as its
/// count says, increasing, with no weights or one per interval.
bool gridIsWhole(const ColumnStore& data) {
  for (std::size_t j = 0; j < data.numPredictors; ++j) {
    if (!data.splitsByThreshold(j)) continue;
    const std::vector<double>& points = data.cutPoints[j];
    if (points.size() != data.numCuts[j]) return false;
    for (std::size_t k = 1; k < points.size(); ++k)
      if (!(points[k - 1] < points[k])) return false;
    if (!data.cutMass[j].empty() && data.cutMass[j].size() != points.size() + 1)
      return false;
  }
  return true;
}

/// Whether every column that splits by threshold holds the codes of one of
/// \p sources' columns (column-major, the store's row count each, null
/// entries skipped) quantized on the grid the store holds: a step that moves
/// the store leaves each column moved or not, never part way.
bool codesAreTheGrids(const ColumnStore& data,
                      std::initializer_list<const double*> sources) {
  std::size_t n = data.numObservations;
  for (std::size_t j = 0; j < data.numPredictors; ++j) {
    if (!data.splitsByThreshold(j)) continue;
    bool any = false;
    for (const double* source : sources) {
      if (source == nullptr || any) continue;
      bool all = true;
      for (std::size_t i = 0; i < n && all; ++i)
        all = data.codeAt(j, i) == data.codeFor(j, source[j * n + i]);
      any = all;
    }
    if (!any) return false;
  }
  return true;
}

struct HarnessTally {
  long allocations = 0, failed = 0, absorbed = 0, consistent = 0,
       described = 0, readable = 0, continued = 0, routed = 0, kept = 0;
  bool otherException = false;
  bool checksKept = false;
};

/// Whether two flattened trees split alike: the same nodes, variables,
/// directions and cut points or level masks, whatever their leaf values.
bool sameStructure(const std::vector<FlatNode>& a,
                   const std::vector<FlatNode>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i].variable != b[i].variable ||
        a[i].numMaskWords != b[i].numMaskWords || a[i].flags != b[i].flags)
      return false;
    if (a[i].variable != invalidVariable && a[i].mask != b[i].mask)
      return false;
  }
  return true;
}

/// Whether every tree of \p state splits as the same tree does in one of
/// \p sweeps, the states an unfailed step passes through at each sweep's
/// end: a move that throws puts its tree back as it was before the
/// proposal, so a failed sweep leaves each tree as the sweep before it or
/// the sweep it was in left it.
bool everyTreeIsOneSwept(const SamplerStateData& state,
                         const std::vector<SamplerStateData>& sweeps) {
  for (std::size_t c = 0; c < state.chains.size(); ++c)
    for (std::size_t f = 0; f < state.chains[c].forests.size(); ++f) {
      const auto& trees = state.chains[c].forests[f].trees;
      for (std::size_t t = 0; t < trees.size(); ++t) {
        bool found = false;
        for (const SamplerStateData& sweep : sweeps)
          found = found ||
            sameStructure(trees[t], sweep.chains[c].forests[f].trees[t]);
        if (!found) return false;
      }
    }
  return true;
}

/// Every forest's per-tree fits, and the variance surface, against those of
/// a new sampler \p make builds and hands \p sampler's state: that one
/// routes the rows through each tree's structure and its leaf values as the
/// state holds them, so a leaf map, table or fit left behind its tree shows
/// here, where the totals, re-summed from the map, agree with it.
template <typename S, typename Make>
bool fitsAreTheTreesRouted(S& sampler, Make&& make) {
  SamplerStateData state;
  sampler.getState(state);
  auto fresh = make();
  if (!fresh->setState(state, nullptr)) return false;
  auto near = [](double a, double b) {
    return std::fabs(a - b) <= 1e-9 * (1.0 + std::fabs(b));
  };
  std::size_t n = sampler.data().numObservations;
  for (std::size_t f = 0; f < sampler.numForests(); ++f) {
    std::size_t numTrees = sampler.chain(0).numTreesInForest(f);
    std::vector<double> live(n * numTrees), routed(n * numTrees);
    TestPeer::forestTreeFits(sampler.chain(0), f, live.data());
    TestPeer::forestTreeFits(fresh->chain(0), f, routed.data());
    for (std::size_t i = 0; i < live.size(); ++i)
      if (!near(live[i], routed[i])) return false;
  }
  if (sampler.chain(0).hasVarianceForest()) {
    const auto& live = TestPeer::varianceForest(sampler.chain(0));
    const auto& routed = TestPeer::varianceForest(fresh->chain(0));
    for (std::size_t i = 0; i < n; ++i)
      if (!near(live.combinedVariance[i], routed.combinedVariance[i]))
        return false;
  }
  return true;
}

/// For every allocation \p step makes from \p start on a sampler \p make
/// builds, a run of it with that allocation failed. Each run is on a new
/// sampler, since a sampler keeps the room a step grew: a second step from
/// the same state allocates less than the first.
/// Where a step leaves the rows its checks route: in the store a new
/// sampler holds, or, for a step that moves the store, in the failed
/// chain's own store alone, a state read from it then installed on a new
/// sampler over \p predictors, the predictors it was built on, since its
/// grid may be another.
enum class RowsIn { newSampler, ownStore };

///
/// \p sweeps, where given, holds the states an unfailed step passes
/// through at each sweep's end, and each failure must leave every tree as
/// one of them has it (everyTreeIsOneSwept). Where \p predictors is given,
/// each column's codes must be those of it or of \p moved, the predictors
/// the step installs, on the store's grid (codesAreTheGrids).
template <typename Make, typename Step>
HarnessTally failEachAllocation(Make&& make, const SamplerStateData& start,
                                Step&& step,
                                RowsIn rowsIn = RowsIn::newSampler,
                                const double* predictors = nullptr,
                                const std::vector<SamplerStateData>* sweeps =
                                  nullptr,
                                const double* moved = nullptr) {
  // the fits the trees route the rows to, through each chain's own store
  // and, where it has not moved, through a new sampler given its state
  auto routed = [&](auto& sampler) {
    for (std::size_t c = 0; c < sampler.numChains(); ++c)
      if (!TestPeer::treesRouteTheirRows(sampler.chain(c))) return false;
    return rowsIn == RowsIn::ownStore || fitsAreTheTreesRouted(sampler, make);
  };
  auto described = [&](auto& sampler) {
    if (!gridIsWhole(sampler.data())) return false;
    if (predictors != nullptr &&
        !codesAreTheGrids(sampler.data(), {predictors, moved}))
      return false;
    for (std::size_t c = 0; c < sampler.numChains(); ++c)
      if (!TestPeer::leafTablesCoverTrees(sampler.chain(c)) ||
          !TestPeer::monotoneLeavesInOrder(sampler.chain(c)) ||
          !TestPeer::everyLeafHoldsARow(sampler.chain(c)) ||
          !TestPeer::leafCovariatesAreStandardized(sampler.chain(c)) ||
          !TestPeer::testStorageCoversTheStore(sampler.chain(c)))
        return false;
    return routed(sampler);
  };
  HarnessTally tally;
  tally.checksKept = sweeps != nullptr;
  Results none;
  {
    auto sampler = make();
    check(sampler->setState(start, nullptr),
          "allocation harness: the start state installs");
    tally.allocations = countAllocations([&]() { step(*sampler); });
  }
  for (long n = 1; n <= tally.allocations; ++n) {
    auto sampler = make();
    sampler->setState(start, nullptr);
    bool reached = false;
    if (!failsAt(n, [&]() { step(*sampler); }, tally.otherException,
                 reached)) {
      if (reached) ++tally.absorbed;
      continue;
    }
    ++tally.failed;
    if (totalsAreTheTrees(*sampler, 1e-12)) ++tally.consistent;
    // the state read before any further sweep is each chain: its tables
    // cover its trees, its monotone leaves are in order, every leaf holds a
    // row, and its fits are those its trees route the rows to
    if (described(*sampler)) ++tally.described;
    // the state installs on a new sampler, so the failed one runs on from
    // what the failure left, not from maps an install rebuilds
    SamplerStateData after;
    sampler->getState(after);
    if (sweeps != nullptr && everyTreeIsOneSwept(after, *sweeps)) ++tally.kept;
    if (make()->setState(after, predictors)) ++tally.readable;
    sampler->run(0, 1, none);
    if (totalsAreTheTrees(*sampler, 1e-9)) ++tally.continued;
    if (routed(*sampler)) ++tally.routed;
  }
  return tally;
}

void report(const HarnessTally& tally, const std::string& name) {
  // a step from a grown state may allocate nothing, as the linear leaf's
  // prior draw does, and then has no failure to test
  check(tally.allocations == 0 || tally.failed > 0,
        (name + ": a failure reaches the caller").c_str());
  check(!tally.otherException,
        (name + ": every failure leaves as an allocation failure").c_str());
  check(tally.consistent == tally.failed,
        (name + ": after each failure the totals are the trees'").c_str());
  check(tally.described == tally.failed,
        (name + ": after each failure the state it reads describes it, the "
                "store's grid whole and its codes its values', its leaf tables "
                "covering its trees, its leaf covariates standardized, its "
                "test storage covering the store's test rows, its "
                "monotone leaves in order, "
                "every leaf holding a row and its fits those its trees route "
                "the rows to").c_str());
  check(!tally.checksKept || tally.kept == tally.failed,
        (name + ": after each failure every tree is as a sweep's end left "
                "it, the proposal in flight taken back").c_str());
  check(tally.readable == tally.failed,
        (name + ": after each failure the state reads and installs").c_str());
  check(tally.continued == tally.failed,
        (name + ": after each failure the next sweep is finite and "
                "consistent").c_str());
  check(tally.routed == tally.failed,
        (name + ": after each failure and a sweep the fits are those the "
                "trees route the rows to").c_str());
  printf("ok: %s (%ld allocations, %ld failures reached the caller, %ld "
         "absorbed)\n", name.c_str(), tally.allocations, tally.failed,
         tally.absorbed);
}

struct HarnessData {
  static constexpr std::size_t n = 80, p = 2;
  std::vector<double> x, y;
  HarnessData() : x(n * p), y(n) {
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = static_cast<double>((i * 37) % n) / static_cast<double>(n);
      x[n + i] = static_cast<double>((i * 11) % n) / static_cast<double>(n);
      y[i] = 2.0 * x[i] + std::sin(6.0 * x[n + i]) +
             0.1 * std::cos(static_cast<double>(i));
    }
  }
};

/// The steps that move the store under the trees: grids of 7 and 5 cuts on
/// columns 0 and 1 (setCutPoints); column 0 shifted, so rows cross a cut,
/// for setPredictor and the per-row update; column 0 reversed, with a new
/// response, for the data replacement; and every column squeezed, for an
/// install that repairs.
struct MovedData {
  std::vector<double> cuts0, cuts1, y, x, column, replaced, squeezed;
  const double* cutColumns[2];
  std::uint32_t numCuts[2];
  std::size_t columns[2] = {0, 1};
  explicit MovedData(const HarnessData& data)
    : y(data.y), x(data.x), column(HarnessData::n), replaced(data.x),
      squeezed(data.x) {
    for (int k = 1; k < 8; ++k) cuts0.push_back(k / 8.0 + 0.003);
    for (int k = 1; k < 6; ++k) cuts1.push_back(k / 6.0 + 0.002);
    cutColumns[0] = cuts0.data();
    cutColumns[1] = cuts1.data();
    numCuts[0] = static_cast<std::uint32_t>(cuts0.size());
    numCuts[1] = static_cast<std::uint32_t>(cuts1.size());
    // the column's values are 1/80 apart, its cuts between them, so a shift
    // of 0.009 takes a row past one
    for (std::size_t i = 0; i < HarnessData::n; ++i) {
      y[i] += 0.3 * std::cos(static_cast<double>(3 * i));
      x[i] += 0.009 * (static_cast<double>(i % 3) - 1.0);
      column[i] = x[i];
      // a replacement reverses the column, which a grid derived again
      // cannot absorb
      replaced[i] = 1.0 - data.x[i];
    }
    // and columns squeezed into their middle tenth empty the leaves past it
    for (std::size_t k = 0; k < HarnessData::n * HarnessData::p; ++k)
      squeezed[k] = 0.45 + 0.1 * data.x[k];
  }
  /// A sampler of \p make's with its columns squeezed under their grids,
  /// which a state grown on their own values cannot stand in unrepaired.
  template <typename Make>
  auto makeSqueezed(Make& make) const {
    auto sampler = make();
    (void) sampler->setPredictor(squeezed.data(), true, false);
    return sampler;
  }
  /// Whether \p state goes into a squeezed sampler only repaired.
  template <typename Make>
  bool installRepairs(Make& make, const SamplerStateData& state) const {
    auto sampler = makeSqueezed(make);
    bool altered = false;
    return sampler->setState(state, squeezed.data(), nullptr, &altered) &&
           altered;
  }
  template <typename S>
  void setCutPoints(S& s, const HarnessData& data) const {
    s.setCutPoints(cutColumns, numCuts, columns, 2, data.x.data());
  }
  template <typename S>
  void setData(S& s) const {
    (void) s.setData(replaced.data(), y.data(), HarnessData::n, nullptr,
                     nullptr, nullptr, 0);
  }
  /// The data replacement with its rows as test rows too.
  template <typename S>
  void setDataWithTest(S& s) const {
    (void) s.setData(replaced.data(), y.data(), HarnessData::n, nullptr,
                     nullptr, replaced.data(), HarnessData::n);
  }
};

template <typename L>
std::unique_ptr<Sampler<L>> makeHarnessSampler(const HarnessData& data,
                                               SamplerOptions options,
                                               ext_rng** rng) {
  options.numTrees = 4;
  return std::make_unique<Sampler<L>>(
    data.x.data(), data.y.data(), HarnessData::n, HarnessData::p, nullptr,
    nullptr, ResponseFamily::gaussian, 1.0, 3.0, 0.37804942330213542, options,
    rng);
}

/// Births and deaths only, from stumps, so the failures land around trees
/// that grow (a birth probability of 1 would leave a birth no reverse move
/// and reject it); then the prior draw of every leaf.
template <typename L>
void runLeafArms(const HarnessData& data, SamplerOptions options,
                 const std::string& leaf, std::uint32_t seed) {
  SamplerOptions everyMove = options;
  options.birthOrDeathProbability = 1.0;
  options.birthProbability = 0.5;
  options.changeProbability = 0.0;
  options.swapProbability = 0.0;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    rngs.push_back(rng);
    return makeHarnessSampler<L>(data, options, &rngs.back());
  };
  auto sampler = make();
  SamplerStateData stumps;
  sampler->getState(stumps);
  report(failEachAllocation(make, stumps,
                            [](Sampler<L>& s) {
                              Results none;
                              s.run(0, 2, none);
                            }),
         "allocation failure in a sweep, " + leaf + " leaf, births");

  Results none;
  sampler->run(5, 0, none);
  SamplerStateData grown;
  sampler->getState(grown);
  report(failEachAllocation(make, grown,
                            [](Sampler<L>& s) {
                              s.sampleNodeParametersFromPrior();
                            }),
         "allocation failure in a prior leaf draw, " + leaf + " leaf");
  report(failEachAllocation(make, grown,
                            [](Sampler<L>& s) { s.sampleTreesFromPrior(); }),
         "allocation failure in a prior tree draw, " + leaf + " leaf");
  // DART's joint prior draw of the weights and every tree, at budgets of 1000
  // and of one (where a refused draw falls to the per-tree draw)
  SamplerOptions dartOptions = options;
  dartOptions.useDart = true;
  auto makeDart = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    rngs.push_back(rng);
    return makeHarnessSampler<L>(data, dartOptions, &rngs.back());
  };
  SamplerStateData dartGrown;
  {
    auto dart = makeDart();
    dart->run(5, 0, none);
    dart->getState(dartGrown);
  }
  for (std::size_t budget : {std::size_t(1000), std::size_t(1)})
    report(failEachAllocation(makeDart, dartGrown,
                              [&](Sampler<L>& s) {
                                (void) s.sampleHyperparametersFromPrior(budget);
                              }),
           "allocation failure in a DART joint prior draw, " + leaf +
             " leaf, budget " + std::to_string(budget));
  report(failEachAllocation(make, stumps,
                            [&](Sampler<L>& s) { s.setState(grown, nullptr); }),
         "allocation failure in a state install, " + leaf + " leaf");
  const std::vector<std::pair<std::size_t, int>> live = {{0, -1}};
  report(failEachAllocation(make, stumps,
                            [&](Sampler<L>& s) {
                              (void) s.installForests(grown, live,
                                                      data.x.data());
                            }),
         "allocation failure in a warm start, " + leaf + " leaf");
  if constexpr (!L::hasVectorParams && !L::hasFunctionParams)
    report(failEachAllocation(make, grown,
                              [](Sampler<L>& s) { s.growFromRoot(1); }),
           "allocation failure in a grow from the root, " + leaf + " leaf");

  // every move the sampler makes by default, from the grown trees, so the
  // failures land in change and swap proposals too
  auto makeEveryMove = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, seed);
    rngs.push_back(rng);
    return makeHarnessSampler<L>(data, everyMove, &rngs.back());
  };
  report(failEachAllocation(makeEveryMove, grown,
                            [](Sampler<L>& s) {
                              Results none;
                              s.run(0, 2, none);
                            }),
         "allocation failure in a sweep, " + leaf + " leaf, every move");

  // one move kind at a time, three sweeps from deeper trees, so swaps find
  // a split below a split, each tree checked against the trees the unfailed
  // sweeps start and end with
  SamplerStateData deep;
  {
    auto deeper = makeEveryMove();
    deeper->setState(grown, nullptr);
    deeper->run(20, 0, none);
    deeper->getState(deep);
  }
  struct MoveKind {
    const char* name;
    double birthOrDeath, swap, perturb, ruleGibbs, change;
  };
  const MoveKind kinds[] = {{"births and deaths", 1, 0, 0, 0, 0},
                            {"changes", 0, 0, 0, 0, 1},
                            {"swaps", 0, 1, 0, 0, 0},
                            {"perturbs", 0, 0, 1, 0, 0},
                            {"rule Gibbs steps", 0, 0, 0, 1, 0}};
  for (const MoveKind& kind : kinds) {
    SamplerOptions only = everyMove;
    only.birthOrDeathProbability = kind.birthOrDeath;
    only.swapProbability = kind.swap;
    only.perturbProbability = kind.perturb;
    only.ruleGibbsProbability = kind.ruleGibbs;
    only.changeProbability = kind.change;
    auto makeOnly = [&]() {
      ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rng, seed);
      rngs.push_back(rng);
      return makeHarnessSampler<L>(data, only, &rngs.back());
    };
    std::vector<SamplerStateData> sweeps(4);
    {
      auto unfailed = makeOnly();
      unfailed->setState(deep, nullptr);
      unfailed->getState(sweeps[0]);
      for (std::size_t k = 1; k < sweeps.size(); ++k) {
        Results none;
        unfailed->run(0, 1, none);
        unfailed->getState(sweeps[k]);
      }
    }
    report(failEachAllocation(makeOnly, deep,
                              [](Sampler<L>& s) {
                                Results none;
                                for (int k = 0; k < 3; ++k) s.run(0, 1, none);
                              },
                              RowsIn::newSampler, nullptr, &sweeps),
           "allocation failure in a move, " + leaf + " leaf, " + kind.name);
  }

  MovedData moved(data);
  report(failEachAllocation(
           make, grown, [&](Sampler<L>& s) { moved.setCutPoints(s, data); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced tree refresh, " + leaf + " leaf");
  report(failEachAllocation(
           make, grown, [&](Sampler<L>& s) { moved.setData(s); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.replaced.data()),
         "allocation failure in a data replacement, " + leaf + " leaf");
  report(failEachAllocation(
           make, grown, [&](Sampler<L>& s) { moved.setDataWithTest(s); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.replaced.data()),
         "allocation failure in a data replacement with test rows, " + leaf +
           " leaf");
  report(failEachAllocation(make, grown,
                            [&](Sampler<L>& s) {
                              s.setTestPredictors(moved.replaced.data(),
                                                  HarnessData::n);
                            }),
         "allocation failure in a test data replacement, " + leaf + " leaf");
  report(failEachAllocation(
           make, grown,
           [&](Sampler<L>& s) { s.setPredictor(moved.x.data(), false, false); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a setPredictor, " + leaf + " leaf");
  report(failEachAllocation(
           make, grown,
           [&](Sampler<L>& s) { s.setPredictor(moved.x.data(), true, true); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced setPredictor with new grids, " + leaf +
           " leaf");
  auto squeezedMake = [&]() { return moved.makeSqueezed(make); };
  check(moved.installRepairs(make, grown),
        ("allocation failure in a repairing state install, " + leaf +
         " leaf: the install repairs a tree").c_str());
  report(failEachAllocation(
           squeezedMake, stumps,
           [&](Sampler<L>& s) { s.setState(grown, moved.squeezed.data()); },
           RowsIn::ownStore, moved.squeezed.data()),
         "allocation failure in a repairing state install, " + leaf + " leaf");
  if constexpr (L::hasVectorParams || L::hasFunctionParams)
    // over the squeezed values, so the constants move from the ones the
    // leaf keeps
    report(failEachAllocation(
             squeezedMake, grown,
             [](Sampler<L>& s) { (void) s.rederiveStandardization(); },
             RowsIn::ownStore, moved.squeezed.data()),
           "allocation failure in a standardization rederive, " + leaf +
             " leaf");
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
}

/// The predictor paths over a constant leaf, and over a variance forest
/// when \p options asks for one: setPredictor forced and not, the per-row
/// update, and a state recorded on another grid installed with its splits
/// moved onto this one; with more than one chain, the forced tree refresh
/// and the data replacement too, so a throw in one chain finds another the
/// store has moved under.
void runPredictorArms(const HarnessData& data, SamplerOptions options,
                      const std::string& name, std::uint32_t seed) {
  using S = Sampler<ConstantGaussianLeaf>;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    std::vector<ext_rng*> own;
    for (std::size_t c = 0; c < options.numChains; ++c) {
      ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
      ext_rng_setSeed(rng, seed + static_cast<std::uint32_t>(c));
      rngs.push_back(rng);
      own.push_back(rng);
    }
    return makeHarnessSampler<ConstantGaussianLeaf>(data, options, own.data());
  };
  MovedData moved(data);
  SamplerStateData grown, otherGrid;
  {
    auto sampler = make();
    Results none;
    sampler->run(5, 0, none);
    sampler->getState(grown);
    moved.setCutPoints(*sampler, data);
    sampler->run(2, 0, none);
    sampler->getState(otherGrid);
  }
  report(failEachAllocation(
           make, grown,
           [&](S& s) { s.setPredictor(moved.x.data(), true, false); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced setPredictor, " + name);
  report(failEachAllocation(
           make, grown,
           [&](S& s) { s.setPredictor(moved.x.data(), false, false); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a setPredictor, " + name);
  report(failEachAllocation(
           make, grown,
           [&](S& s) { s.setPredictor(moved.x.data(), true, true); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced setPredictor with new grids, " + name);
  // from the 7 and 5 cuts setCutPoints gave, so the grids derived again are
  // longer than the ones they replace
  auto makeFewCuts = [&]() {
    auto sampler = make();
    moved.setCutPoints(*sampler, data);
    return sampler;
  };
  report(failEachAllocation(
           makeFewCuts, otherGrid,
           [&](S& s) { s.setPredictor(moved.x.data(), true, true); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced setPredictor growing the grids, " +
           name);
  report(failEachAllocation(
           make, grown,
           [&](S& s) { s.setPredictor(moved.x.data(), false, true); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a setPredictor with new grids, " + name);
  report(failEachAllocation(
           make, grown,
           [&](S& s) {
             bool installed[HarnessData::n] = {};
             (void) s.updatePredictorPerObservation(moved.column.data(), 0,
                                                    installed);
           },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a per-row predictor update, " + name);
  report(failEachAllocation(
           make, grown,
           [&](S& s) { (void) s.setState(otherGrid, data.x.data()); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a state install across grids, " + name);
  std::vector<std::pair<std::size_t, int>> live;
  for (std::size_t c = 0; c < options.numChains; ++c)
    live.push_back({c, -1});
  report(failEachAllocation(
           make, grown,
           [&](S& s) {
             (void) s.installForests(otherGrid, live, data.x.data());
           },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a warm start across grids, " + name);
  if (options.numChains > 1) {
    report(failEachAllocation(
             make, grown, [&](S& s) { moved.setCutPoints(s, data); },
             RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
           "allocation failure in a forced tree refresh, " + name);
    report(failEachAllocation(
             make, grown, [&](S& s) { moved.setData(s); },
             RowsIn::ownStore, data.x.data(), nullptr,
             moved.replaced.data()),
           "allocation failure in a data replacement, " + name);
  }
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
}

/// A constant leaf with a variance forest: a sweep of births and deaths from
/// stumps, the variance forest's prior draw, a state's install and a grow
/// from the root, the variance surface checked with the fits.
void runVarianceArms(const HarnessData& data) {
  SamplerOptions options;
  options.numVarianceTrees = 3;
  options.birthOrDeathProbability = 1.0;
  options.birthProbability = 0.5;
  options.changeProbability = 0.0;
  options.swapProbability = 0.0;
  using S = Sampler<ConstantGaussianLeaf>;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 20261121u);
    rngs.push_back(rng);
    return makeHarnessSampler<ConstantGaussianLeaf>(data, options,
                                                    &rngs.back());
  };
  auto sampler = make();
  check(sampler->chain(0).hasVarianceForest(),
        "allocation failure, variance forest: the sampler has one");
  SamplerStateData stumps, grown;
  sampler->getState(stumps);
  Results none;
  sampler->run(5, 0, none);
  sampler->getState(grown);
  report(failEachAllocation(make, stumps,
                            [](S& s) {
                              Results none;
                              s.run(0, 2, none);
                            }),
         "allocation failure in a sweep, variance forest");
  report(failEachAllocation(make, grown,
                            [](S& s) { s.sampleVarianceForestFromPrior(); }),
         "allocation failure in a variance prior draw");
  report(failEachAllocation(make, stumps,
                            [&](S& s) { s.setState(grown, nullptr); }),
         "allocation failure in a state install, variance forest");
  report(failEachAllocation(make, grown, [](S& s) { s.growFromRoot(1); }),
         "allocation failure in a grow from the root, variance forest");
  MovedData moved(data);
  report(failEachAllocation(
           make, grown, [&](S& s) { moved.setCutPoints(s, data); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.x.data()),
         "allocation failure in a forced tree refresh, variance forest");
  report(failEachAllocation(
           make, grown, [&](S& s) { moved.setData(s); },
           RowsIn::ownStore, data.x.data(), nullptr, moved.replaced.data()),
         "allocation failure in a data replacement, variance forest");
  report(failEachAllocation(make, grown,
                            [&](S& s) {
                              s.setTestPredictors(moved.replaced.data(),
                                                  HarnessData::n);
                            }),
         "allocation failure in a test data replacement, variance forest");
  auto squeezedMake = [&]() { return moved.makeSqueezed(make); };
  check(moved.installRepairs(make, grown),
        "allocation failure in a repairing state install, variance forest: "
        "the install repairs a tree");
  report(failEachAllocation(
           squeezedMake, stumps,
           [&](S& s) { s.setState(grown, moved.squeezed.data()); },
           RowsIn::ownStore, moved.squeezed.data()),
         "allocation failure in a repairing state install, variance forest");
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
}

/// The monotone leaf's settle: a grown tree's leaf values put out of the
/// order, as a move a failure stops part way can leave a new leaf's slot,
/// then settled with that tree in flight. Its leaves are back in order, the
/// totals are the trees' and the next sweep keeps both.
void runMonotoneSettleArm(const HarnessData& data) {
  const std::int8_t directions[] = {1, 0};
  SamplerOptions options;
  options.monotoneDirections = directions;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 20261122u);
  auto sampler =
    makeHarnessSampler<MonotoneConstantGaussianLeaf>(data, options, &rng);
  Results none;
  sampler->run(5, 0, none);
  auto& chain = sampler->chain(0);
  std::size_t numTrees = chain.numTreesInForest(0), found = numTrees;
  for (std::size_t t = 0; t < numTrees && found == numTrees; ++t)
    for (double sign : {1.0, -1.0}) {
      std::vector<double>& mu = TestPeer::muByTree(chain, t);
      std::vector<double> kept = mu;
      for (std::size_t i = 0; i < mu.size(); ++i)
        mu[i] = (i % 2 == 0 ? -sign : sign) * static_cast<double>(i + 1);
      if (!TestPeer::monotoneLeavesInOrder(chain)) {
        found = t;
        break;
      }
      mu = kept;
    }
  check(found < numTrees,
        "allocation failure, monotone settle: a tree's leaves leave the order");
  TestPeer::settleTreeInFlight(chain, 0, found);
  check(TestPeer::monotoneLeavesInOrder(chain) &&
          totalsAreTheTrees(*sampler, 1e-12),
        "allocation failure, monotone settle: the tree in flight is back in "
        "order and the totals are the trees'");
  sampler->run(0, 1, none);
  check(TestPeer::monotoneLeavesInOrder(chain) &&
          totalsAreTheTrees(*sampler, 1e-9),
        "allocation failure, monotone settle: the next sweep keeps both");
  printf("ok: allocation failure, monotone settle\n");
  ext_rng_destroy(rng);
}

/// A tree's mask-pool compaction under each allocation failure it meets: two
/// live pooled rules stored past garbage, so the copy moves both, and a
/// failure part way through it would leave a rule pointing into a pool never
/// swapped in. After each failure every rule reads the words it held.
void runMaskPoolCompactionArm() {
  const std::size_t n = 60, K = 70;  // past 63 categories: two pooled words
  std::vector<double> x(n), y(n);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i * 37) % K);
    y[i] = static_cast<double>(i % 5);
  }
  const ColumnKind types[] = {ColumnKind::categorical};
  ColumnStore store;
  check(store.build(x.data(), n, 1, 10, false, types) &&
          store.columnIsPooled(0),
        "allocation failure, mask pool: the factor's masks are pooled");
  std::vector<index_t> indices(n);
  // the root sends categories 1 and 40 right, its left child 2 and 3
  auto grow = [&](Tree& tree) {
    tree.initialize(indices.data(), n);
    for (int g = 0; g < 200; ++g) tree.allocateMask(2);
    std::int32_t node = 0;
    for (std::uint32_t first : {1u, 2u}) {
      std::size_t offset = tree.allocateMask(2);
      std::uint64_t* words = tree.mutableMaskWordsFor(offset);
      maskSetBit(words, first);
      maskSetBit(words, first == 1u ? 40u : 3u);
      Rule rule;
      rule.variableIndex = 0;
      rule.setMaskOffset(offset);
      tree.birth(store, node, rule, y.data(), nullptr);
      node = tree.at(node).leftChild;
    }
  };
  auto ruleWords = [](const Tree& tree) {
    std::vector<std::uint64_t> words;
    for (std::int32_t i : {0, tree.at(0).leftChild}) {
      const std::uint64_t* w = tree.maskWordsFor(tree.at(i).rule);
      words.insert(words.end(), w, w + 2);
    }
    return words;
  };
  long allocations = 0;
  {
    Tree tree;
    grow(tree);
    std::vector<std::uint64_t> before = ruleWords(tree);
    allocations =
      countAllocations([&]() { tree.compactMaskPoolIfNeeded(store); });
    check(tree.maskPool.size() == 4 && ruleWords(tree) == before,
          "allocation failure, mask pool: the compaction keeps the two live "
          "masks");
  }
  long failed = 0, intact = 0;
  bool otherException = false;
  for (long k = 1; k <= allocations; ++k) {
    Tree tree;
    grow(tree);
    std::vector<std::uint64_t> before = ruleWords(tree);
    bool reached = false;
    if (!failsAt(k, [&]() { tree.compactMaskPoolIfNeeded(store); },
                 otherException, reached))
      continue;
    ++failed;
    if (ruleWords(tree) == before) ++intact;
  }
  check(allocations > 0 && failed == allocations && !otherException,
        "allocation failure, mask pool: each failure reaches the caller");
  check(intact == failed, "allocation failure, mask pool: after each failure "
                          "every rule reads the words it held");
  printf("ok: allocation failure in a mask-pool compaction (%ld allocations, "
         "%ld failures reached the caller)\n", allocations, failed);
}

/// The harness's rows with column 1 a factor of 60000 declared categories,
/// 40 of them observed, which the response follows. Past 63 categories a
/// rule's mask lives in its tree's pool, 938 words a rule here, so a tree's
/// first accepted split on the factor crosses the pool's first high-water
/// mark and the sweep compacts the pool between that move and the tree's
/// leaf map catching up with it.
HarnessData pooledHarnessData() {
  HarnessData data;
  const std::size_t n = HarnessData::n;
  for (std::size_t i = 0; i < n; ++i) {
    double level = static_cast<double>((i * 7) % 40);
    data.x[n + i] = level;
    data.y[i] = data.x[i] + 2.0 * std::sin(level) +
                0.1 * std::cos(static_cast<double>(i));
  }
  return data;
}

/// Births and deaths from stumps on the pooled factor's rows: a failure in the
/// compaction leaves an accepted move's tree ahead of its leaf map, which the
/// settle marks for a rebuild at the tree's next draw.
void runPooledFactorArm() {
  HarnessData data = pooledHarnessData();
  const ColumnKind types[HarnessData::p] = {ColumnKind::numeric,
                                            ColumnKind::categorical};
  const std::uint32_t counts[HarnessData::p] = {0, 60000};
  SamplerOptions options;
  options.predictors.columnTypes = types;
  options.predictors.categoryCounts = counts;
  options.birthOrDeathProbability = 1.0;
  options.birthProbability = 0.5;
  options.changeProbability = 0.0;
  options.swapProbability = 0.0;
  using S = Sampler<ConstantGaussianLeaf>;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 20261113u);
    rngs.push_back(rng);
    return makeHarnessSampler<ConstantGaussianLeaf>(data, options,
                                                    &rngs.back());
  };
  SamplerStateData stumps;
  long pooledSplits = 0;
  {
    auto sampler = make();
    check(sampler->data().columnIsPooled(1),
          "allocation failure, pooled factor: the factor's masks are pooled");
    sampler->getState(stumps);
    Results none;
    sampler->run(0, 3, none);
    for (std::size_t t = 0; t < sampler->chain(0).numTrees(); ++t)
      for (const Node& node : sampler->chain(0).tree(t).nodes)
        if (!node.isBottom() && node.rule.variableIndex == 1) ++pooledSplits;
  }
  check(pooledSplits > 0,
        "allocation failure, pooled factor: the step splits on the factor");
  report(failEachAllocation(make, stumps,
                            [](S& s) {
                              Results none;
                              s.run(0, 3, none);
                            }),
         "allocation failure in a sweep, constant leaf, a pooled factor");

  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
}

/// A full ring of kept draws: a sweep that fails part way drops the slot it
/// began, and the store serves the others, oldest first, as run 1 wrote
/// them.
void runSavedRingArm(const HarnessData& data) {
  const std::size_t capacity = 6;
  SamplerOptions options;
  options.keepTrees = true;
  options.numSamplesToStore = capacity;
  using S = Sampler<ConstantGaussianLeaf>;
  std::vector<ext_rng*> rngs;
  auto make = [&]() {
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    ext_rng_setSeed(rng, 20261106u);
    rngs.push_back(rng);
    return makeHarnessSampler<ConstantGaussianLeaf>(data, options,
                                                    &rngs.back());
  };
  auto served = [](S& sampler) {
    std::vector<std::vector<std::vector<FlatNode>>> draws;
    for (std::size_t d = 0; d < sampler.filledSavedDraws(); ++d) {
      std::size_t slot = sampler.savedSlotForDraw(d);
      std::vector<std::vector<FlatNode>> trees;
      for (std::size_t t = 0; t < sampler.chain(0).numTrees(); ++t)
        trees.push_back(sampler.savedTree(0, slot, t));
      draws.push_back(trees);
    }
    return draws;
  };
  auto sameNode = [](const FlatNode& a, const FlatNode& b) {
    return a.variable == b.variable && a.numMaskWords == b.numMaskWords &&
           a.mask == b.mask && a.flags == b.flags;
  };
  auto sameDraw = [&](const std::vector<std::vector<FlatNode>>& a,
                      const std::vector<std::vector<FlatNode>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t t = 0; t < a.size(); ++t) {
      if (a[t].size() != b[t].size()) return false;
      for (std::size_t i = 0; i < a[t].size(); ++i)
        if (!sameNode(a[t][i], b[t][i])) return false;
    }
    return true;
  };
  Results none;
  SamplerStateData full;
  std::vector<std::vector<std::vector<FlatNode>>> before;
  {
    auto sampler = make();
    sampler->run(20, capacity, none);
    sampler->getState(full);
    before = served(*sampler);
  }
  check(before.size() == capacity,
        "allocation failure, kept draws: the ring is full");
  long allocations;
  {
    auto sampler = make();
    sampler->setState(full, nullptr);
    allocations = countAllocations([&]() { sampler->run(0, 1, none); });
  }
  long failed = 0, kept = 0, dropped = 0;
  bool otherException = false;
  for (long n = 1; n <= allocations; ++n) {
    auto sampler = make();
    sampler->setState(full, nullptr);
    bool reached = false;
    if (!failsAt(n, [&]() { sampler->run(0, 1, none); }, otherException,
                 reached))
      continue;
    ++failed;
    auto after = served(*sampler);
    bool tail = after.size() <= capacity;
    for (std::size_t d = 0; d < after.size() && tail; ++d)
      tail = sameDraw(after[d], before[before.size() - after.size() + d]);
    if (tail) ++kept;
    if (after.size() < capacity) ++dropped;
  }
  check(failed > 0 && dropped > 0 && !otherException,
        "allocation failure, kept draws: failures reach the caller, some "
        "inside the kept sweep");
  check(kept == failed,
        "allocation failure, kept draws: the store serves run 1's newest "
        "draws the stopped sweep did not touch, oldest first");
  printf("ok: allocation failure with kept draws (%ld allocations, %ld "
         "failures, %ld inside the kept sweep)\n", allocations, failed,
         dropped);
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
}

}  // namespace

void runAllocTests() {
  HarnessData data;
  runLeafArms<ConstantGaussianLeaf>(data, SamplerOptions(), "constant",
                                    20261106u);
  const std::size_t covariates[] = {0};
  SamplerOptions linear;
  linear.leafCovariateColumns = covariates;
  linear.numLeafCovariates = 1;
  runLeafArms<LinearGaussianLeaf>(data, linear, "linear", 20261107u);
  SamplerOptions gp = linear;
  gp.gpMaxLeafSize = 100;
  runLeafArms<GPGaussianLeaf>(data, gp, "gp", 20261108u);
  const std::int8_t directions[] = {1, 0};
  SamplerOptions monotone;
  monotone.monotoneDirections = directions;
  runLeafArms<MonotoneConstantGaussianLeaf>(data, monotone, "monotone",
                                            20261109u);
  runVarianceArms(data);
  runPredictorArms(data, SamplerOptions(), "constant leaf", 20261124u);
  SamplerOptions variance;
  variance.numVarianceTrees = 3;
  runPredictorArms(data, variance, "variance forest", 20261125u);
  // a throw in one chain's rebuild, and in the store before any chain's,
  // with another chain the store has moved under
  SamplerOptions twoChains = variance;
  twoChains.numChains = 2;
  runPredictorArms(data, twoChains, "two chains", 20261126u);
  runMonotoneSettleArm(data);
  runMaskPoolCompactionArm();
  runPooledFactorArm();
  runSavedRingArm(data);
}
