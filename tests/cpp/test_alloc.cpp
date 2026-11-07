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

/// Every forest's totalFits against its trees' fits summed in tree order,
/// and every value finite.
template <typename S>
bool totalsAreTheTrees(S& sampler, double tolerance) {
  auto& chain = sampler.chain(0);
  std::size_t n = sampler.data().numObservations;
  std::size_t numForests = sampler.numForests();
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
  return true;
}

struct HarnessTally {
  long allocations = 0, failed = 0, absorbed = 0, consistent = 0,
       readable = 0, continued = 0, routed = 0;
  bool otherException = false;
};

/// Every forest's per-tree fits against those of a new sampler \p make
/// builds and hands \p sampler's state: that one routes the rows through
/// each tree's structure, so a leaf map left behind its tree shows here,
/// where the totals, re-summed from the map, agree with it.
template <typename S, typename Make>
bool fitsAreTheTreesRouted(S& sampler, Make&& make) {
  SamplerStateData state;
  sampler.getState(state);
  auto fresh = make();
  if (!fresh->setState(state, nullptr)) return false;
  std::size_t n = sampler.data().numObservations;
  for (std::size_t f = 0; f < sampler.numForests(); ++f) {
    std::size_t numTrees = sampler.chain(0).numTreesInForest(f);
    std::vector<double> live(n * numTrees), routed(n * numTrees);
    TestPeer::forestTreeFits(sampler.chain(0), f, live.data());
    TestPeer::forestTreeFits(fresh->chain(0), f, routed.data());
    for (std::size_t i = 0; i < live.size(); ++i)
      if (!(std::fabs(live[i] - routed[i]) <=
            1e-9 * (1.0 + std::fabs(routed[i]))))
        return false;
  }
  return true;
}

/// For every allocation \p step makes from \p start on a sampler \p make
/// builds, a run of it with that allocation failed. Each run is on a new
/// sampler, since a sampler keeps the room a step grew: a second step from
/// the same state allocates less than the first.
template <typename Make, typename Step>
HarnessTally failEachAllocation(Make&& make, const SamplerStateData& start,
                                Step&& step) {
  HarnessTally tally;
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
    // the state installs on a new sampler, so the failed one runs on from
    // what the failure left, not from maps an install rebuilds
    SamplerStateData after;
    sampler->getState(after);
    if (make()->setState(after, nullptr)) ++tally.readable;
    sampler->run(0, 1, none);
    if (totalsAreTheTrees(*sampler, 1e-9)) ++tally.continued;
    if (fitsAreTheTreesRouted(*sampler, make)) ++tally.routed;
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
  for (ext_rng* rng : rngs) ext_rng_destroy(rng);
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
  runMaskPoolCompactionArm();
  runPooledFactorArm();
  runSavedRingArm(data);
}
