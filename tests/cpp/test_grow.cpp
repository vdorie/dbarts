// Component tests for grow-from-root (grow.hpp + Chain::growForestFromRoot):
// seeded determinism and the documented per-node draw count, the well-formed /
// legal-chain-state invariants of a grown tree, the categorical rules it now
// places and their gauge, the prior mass an enumerated cut candidate carries on
// a column with missing values, the law the categorical root draw realizes in
// both enumeration branches and the closed-form mass invariant behind it, the
// 1 + A coins a categorical rule assembly spends, and grow-then-MH
// continuation through a live sampler.
#include "common.hpp"

#include <ctime>

namespace {

// deterministic grows allocate nodes in the same order, so the arenas compare
// index by index
bool treeStructureEqual(const Tree& a, const Tree& b) {
  if (a.nodes.size() != b.nodes.size()) return false;
  for (size_t i = 0; i < a.nodes.size(); ++i) {
    const Node& x(a.nodes[i]);
    const Node& y(b.nodes[i]);
    if (x.leftChild != y.leftChild || x.parent != y.parent) return false;
    if (!x.isBottom() &&
        (x.rule.variableIndex != y.rule.variableIndex || x.rule.bits != y.rule.bits))
      return false;
  }
  return true;
}

// every node the recursion draws at has positive prior growth probability; on a
// no-missing, all-ordinal build there are no missing coins and no categorical
// orientation or absent-position coins, so this is the exact uniform count
size_t positiveGrowthNodeCount(const Tree& tree, const ColumnStore& store,
                               const CGMTreePrior& prior) {
  size_t count = 0;
  for (size_t i = 0; i < tree.nodes.size(); ++i)
    if (prior.growthProbability(tree, store,
                                static_cast<int32_t>(i)) > 0.0)
      ++count;
  return count;
}

void growOnce(const ColumnStore& store, const CGMTreePrior& prior,
              const ConstantGaussianLeaf& leaf, ext_rng* rng, Tree& tree,
              std::vector<index_t>& indexBuffer, const double* y, size_t n) {
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y, nullptr);
  GrowScratch scratch;
  growTreeFromRoot(store, prior, leaf, rng, tree, 0, y, nullptr, 2.0, 0.9,
                   scratch);
}

void makeStepData(std::vector<double>& x, std::vector<double>& y, size_t n) {
  x.resize(n);
  y.resize(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    y[i] = (x[i] < 0.4 ? -2.0 : (x[i] < 0.7 ? 0.5 : 3.0)) + 0.1 * (runif01() - 0.5);
  }
}

void testDeterminismAndDrawCount() {
  const size_t n = 300;
  std::vector<double> x, y;
  makeStepData(x, y, n);
  ColumnStore store;
  store.build(x.data(), n, 1, 50);
  check(store.hasMissing[0] == 0 && store.types[0] == ColumnType::ordinal,
        "no-missing, all-ordinal build for the draw-count census");

  CGMTreePrior prior;  // base 0.95, power 2
  ConstantGaussianLeaf leaf{0.5};

  std::vector<index_t> bufA(n), bufB(n), bufC(n);
  Tree treeA, treeB, treeC;

  ext_rng* rngA = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng* rngB = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng* rngC = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngA, 20260710u);
  ext_rng_setSeed(rngB, 20260710u);
  ext_rng_setSeed(rngC, 99u);

  growOnce(store, prior, leaf, rngA, treeA, bufA, y.data(), n);
  growOnce(store, prior, leaf, rngB, treeB, bufB, y.data(), n);
  growOnce(store, prior, leaf, rngC, treeC, bufC, y.data(), n);

  check(treeStructureEqual(treeA, treeB),
        "same seed grows the same tree (determinism)");
  check(treeA.nodes.size() > 1, "the strong signal grows past the root");

  // exact draw count: replay the census many continuous uniforms on a fresh
  // generator and confirm it reaches the same serialized state as the grow.
  // One uniform per positive-growth node is the WHOLE count only because this
  // fixture spends no post-draw coin of either kind
  size_t expectedDraws = positiveGrowthNodeCount(treeA, store, prior);
  ext_rng* replay = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(replay, 20260710u);
  for (size_t d = 0; d < expectedDraws; ++d)
    (void) ext_rng_simulateContinuousUniform(replay);

  size_t len = ext_rng_getSerializedStateLength(rngA);
  std::vector<unsigned char> stateGrow(len), stateReplay(len);
  ext_rng_writeSerializedState(rngA, stateGrow.data());
  ext_rng_writeSerializedState(replay, stateReplay.data());
  check(stateGrow == stateReplay,
        "grow spends exactly one uniform per positive-growth node");

  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngC);
  ext_rng_destroy(replay);
  printf("ok: grow determinism and draw count\n");
}

void testVetoDrawsNothing() {
  const size_t n = 40;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) { x[i] = runif01(); y[i] = runif01(); }
  ColumnStore store;
  store.build(x.data(), n, 1, 20);

  // base 0 -> growthProbability is zero at every node, so the root draws
  // nothing and stays a single leaf
  CGMTreePrior prior;
  prior.base = 0.0;
  ConstantGaussianLeaf leaf{0.5};

  std::vector<index_t> buf(n);
  Tree tree;
  tree.initialize(buf.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 7u);
  size_t len = ext_rng_getSerializedStateLength(rng);
  std::vector<unsigned char> before(len), after(len);
  ext_rng_writeSerializedState(rng, before.data());

  GrowScratch scratch;
  growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0, 0.9,
                   scratch);
  ext_rng_writeSerializedState(rng, after.data());

  check(tree.hasSingleNode(), "a growth-vetoed root stays a leaf");
  check(before == after, "a growth-vetoed node draws nothing");

  ext_rng_destroy(rng);
  printf("ok: grow veto draws nothing\n");
}

void testGrownTreeWellFormed() {
  const size_t n = 400;
  std::vector<double> x, y;
  makeStepData(x, y, n);
  ColumnStore store;
  store.build(x.data(), n, 1, 40);

  CGMTreePrior prior;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> buf(n);
  Tree tree;

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 123u);
  growOnce(store, prior, leaf, rng, tree, buf, y.data(), n);
  ext_rng_destroy(rng);

  check(tree.nodes.size() > 1, "the tree grew");
  check(tree.bottomNodesAreOccupied(), "no grown leaf is empty (occupancy)");

  // partitions cover all observations exactly once and every ordinal rule is in
  // gauge (its split index within the ancestor-constrained interval)
  std::vector<int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  size_t covered = 0;
  bool occupied = true;
  for (int32_t b : bottoms) {
    covered += tree.at(b).numObservations();
    occupied &= tree.at(b).numObservations() > 0;
  }
  check(covered == n && occupied, "leaves partition every observation");

  std::vector<int32_t> internal;
  tree.fillNotBottom(0, internal);
  bool inGauge = true;
  for (int32_t node : internal) {
    int32_t var = tree.at(node).rule.variableIndex;
    int32_t left, right;
    tree.splitInterval(store, node, var, &left, &right);
    int32_t split = tree.at(node).rule.splitIndex();
    inGauge &= split >= left && split <= right;
  }
  check(inGauge, "every grown ordinal rule is in gauge");

  printf("ok: grown tree well-formed\n");
}

// The v1 "categorical predictors are never split here" contract, inverted: a
// categorical column is scanned and split like any other. The fixture keeps the
// signal ORDINAL, so categorical rules appear only because they carry their
// family's prior mass - the weaker and more informative direction. Every one
// must be in gauge (nonzero, a strict subset of the node's reachable set,
// unequal to it) so buildFromFlat would accept the grown tree.
void testCategoricalSplits() {
  const size_t n = 300;
  const uint32_t numLevels = 5, numObserved = 4;
  // column 0 categorical, column 1 ordinal with the signal
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % numObserved);
    x[i + n] = runif01();
    y[i] = (x[i + n] < 0.5 ? -2.0 : 2.0) + 0.1 * (runif01() - 0.5);
  }
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  // one declared level no row observes: reachable everywhere, present nowhere,
  // so A = 1 at the root and the absent-position coin runs
  uint32_t categoryCounts[] = {numLevels, 0};
  ColumnStore store;
  store.build(x.data(), n, 2, 20, false, types, nullptr, 0, categoryCounts);
  check(store.numCuts[0] == numLevels && !store.columnIsPooled(0),
        "the categorical fixture declares five inline levels");

  CGMTreePrior prior;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> buf(n);
  Tree tree;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 555u);
  GrowScratch scratch;

  size_t numCategorical = 0, numGrown = 0, absentRight = 0, absentSeen = 0;
  bool inGauge = true, occupied = true, everOrdinal = false;
  for (int iter = 0; iter < 200; ++iter) {
    tree.initialize(buf.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0,
                     0.9, scratch);
    occupied &= tree.bottomNodesAreOccupied();
    std::vector<int32_t> internal;
    tree.fillNotBottom(0, internal);
    if (!internal.empty()) ++numGrown;
    for (int32_t node : internal) {
      if (tree.at(node).rule.variableIndex != 0) { everOrdinal = true; continue; }
      ++numCategorical;
      uint64_t mask = tree.at(node).rule.categoryDirections();
      uint64_t reachable = tree.reachableCategories(store, node, 0);
      inGauge &= mask != 0 && (mask & ~reachable) == 0 && mask != reachable;
      if ((reachable >> (numLevels - 1)) & 1ull) {
        ++absentSeen;  // the never-observed level, drawn by its own coin
        if ((mask >> (numLevels - 1)) & 1ull) ++absentRight;
      }
    }
  }
  ext_rng_destroy(rng);

  check(numGrown > 0 && everOrdinal, "the ordinal signal grew the tree");
  check(numCategorical > 0,
        "grow-from-root places categorical rules (the v1 contract inverted)");
  check(inGauge, "every grown categorical mask is in gauge");
  check(occupied, "categorical grow keeps every leaf occupied");
  check(absentSeen > 0 && absentRight > 0 && absentRight < absentSeen,
        "an absent reachable category is drawn, not pinned to one side");

  printf("ok: grow categorical splits (%zu rules over %zu grown trees)\n",
         numCategorical, numGrown);
}

// The pooled tier: past 63 categories a mask no longer fits the rule word, so
// grow allocates it in the tree's own pool and writes it through
// mutableMaskWordsFor. Same gauge, read through reachableCategoriesWide. Its
// missing pseudo-category sits at position K rather than 63, and the fixture
// leaves categories reachable but ABSENT at every node, so the absent-position
// coins run.
void testPooledCategoricalGrow() {
  // a private generator and a saved global-rng snapshot keep the shared stream
  // (and the seed-pinned suites that follow) bitwise intact
  uint64_t savedRngState = rngState;
  rngState = 31415u;
  const size_t n = 420;
  const uint32_t numLevels = 70, numObserved = 60;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % numObserved);
    y[i] = (i % numObserved < 30 ? -1.5 : 1.5) + 0.1 * (runif01() - 0.5);
  }
  x[5] = std::nan("");
  ColumnType types[] = {ColumnType::categorical};
  // the declared level count exceeds what any row observes, so ten categories
  // are reachable but absent at every node and their coins really run
  uint32_t categoryCounts[] = {numLevels};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types, nullptr, 0, categoryCounts);
  check(store.columnIsPooled(0) && store.hasMissing[0] == 1 &&
          store.numCuts[0] == numLevels,
        "the wide fixture pools its mask and carries a missing value");

  CGMTreePrior prior;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> buf(n);
  Tree tree;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 8080u);
  GrowScratch scratch;
  size_t numWords = maskWordsForCount(numLevels);
  std::vector<uint64_t> reachable(numWords);

  size_t numRules = 0, absentRight = 0, absentSeen = 0;
  bool inGauge = true, occupied = true;
  for (int iter = 0; iter < 40; ++iter) {
    tree.initialize(buf.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0,
                     0.9, scratch);
    occupied &= tree.bottomNodesAreOccupied();
    std::vector<int32_t> internal;
    tree.fillNotBottom(0, internal);
    for (int32_t node : internal) {
      ++numRules;
      const uint64_t* mask = tree.maskWordsFor(tree.at(node).rule);
      tree.reachableCategoriesWide(store, node, 0, reachable.data());
      inGauge &= !maskIsZero(mask, numWords) &&
                 maskIsSubsetOf(mask, reachable.data(), numWords) &&
                 !maskEquals(mask, reachable.data(), numWords);
      // the never-observed levels: a coin each, so both directions appear
      for (uint32_t c = numObserved; c < numLevels; ++c)
        if (maskTestBit(reachable.data(), c)) {
          ++absentSeen;
          if (maskTestBit(mask, c)) ++absentRight;
        }
    }
  }
  ext_rng_destroy(rng);

  check(numRules > 0, "the pooled column grows past the root");
  check(inGauge, "every grown pooled mask is in gauge");
  check(occupied, "pooled grow keeps every leaf occupied");
  check(absentSeen > 0 && absentRight > 0 && absentRight < absentSeen,
        "absent reachable categories are drawn, not pinned to one side");

  rngState = savedRngState;
  printf("ok: grow pooled categorical (%zu rules)\n", numRules);
}

// grow in place through a live sampler, then confirm the forest is a legal
// chain state MH continues from with coherent fits
void testGrowThenContinue(ext_rng* rng) {
  const size_t n = 250;
  std::vector<double> x, f;
  makeMutationData(x, f, n);

  SamplerOptions options;
  options.numTrees = 20;
  options.nodeScale = 3.0;
  ConstantLeafSampler sampler(x.data(), f.data(), n, 2, nullptr, nullptr,
                              ResponseFamily::gaussian, 1.0, 3.0,
                              0.37804942330213542, options, &rng);

  sampler.chain(0).growForestFromRoot(3);

  bool occupied = true, grew = false;
  for (size_t t = 0; t < options.numTrees; ++t) {
    occupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
    grew |= sampler.chain(0).tree(t).nodes.size() > 1;
  }
  check(occupied, "grown forest has no empty leaves");
  check(grew, "grow produced structure");

  const std::vector<double>& total = sampler.chain(0).totalFits();
  bool finiteTotal = true;
  for (double v : total) finiteTotal &= std::isfinite(v);
  check(finiteTotal, "grown forest total fits are finite");
  check(std::isfinite(sampler.chain(0).sumOfSquaredResiduals()),
        "grown forest residual is finite");

  // the exact MH sweeps own the forest from here
  std::vector<double> trainingFits(n * 5);
  Results results;
  results.trainingFits = trainingFits.data();
  sampler.run(0, 5, results);
  bool finite = true;
  for (double v : trainingFits) finite &= std::isfinite(v);
  check(finite, "MH continues from the grown forest with finite fits");

  bool stillOccupied = true;
  for (size_t t = 0; t < options.numTrees; ++t)
    stillOccupied &= sampler.chain(0).tree(t).bottomNodesAreOccupied();
  check(stillOccupied, "post-continuation forest stays occupied");

  printf("ok: grow then continue\n");
}

// grow-from-root is a FIFTH split generator (design "Containment"); it routes
// availability through collectAvailableVariables, so it is correct-by-
// construction ONCE the predicate carries the constraint - but that must be
// named and tested. Under a max-order cap + a forbidden pair, every grown tree
// must satisfy the whole-subtree feasibility walk and honor the order cap.
void testGrowHonorsInteraction() {
  // a private generator and a saved global-rng snapshot keep the shared stream
  // (and the seed-pinned suites that follow) bitwise intact
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, 4242);
  std::uint64_t savedRngState = rngState;
  rngState = 271828u;
  const size_t n = 400, p = 4;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j)
      x[i + j * n] = static_cast<double>((i * (2 * j + 1) + j) % 13) / 13.0;
    // a multi-way signal so the free grow would reach for several variables
    y[i] = (x[i] > 0.5 ? 1.0 : -1.0) + (x[i + n] > 0.5 ? 0.7 : -0.7) +
           (x[i + 2 * n] > 0.5 ? 0.5 : -0.5) + 0.1 * (runif01() - 0.5);
  }
  ColumnStore store;
  store.build(x.data(), n, p, 12);

  size_t forbiddenPair[] = {0, 1};       // x0 and x1 may not co-occur
  InteractionConstraint constraint;
  constraint.build(p, 2, forbiddenPair, 1);  // max.order 2 AND forbid(x0, x1)

  CGMTreePrior prior;
  prior.base = 0.95;
  prior.power = 1.5;  // encourage deeper trees so the cap actually bites
  ConstantGaussianLeaf leaf{0.5};

  std::vector<index_t> indexBuffer(n);
  Tree tree;
  bool allFeasible = true, everSplit = false, everDeep = false;
  for (int iter = 0; iter < 400; ++iter) {
    tree.initialize(indexBuffer.data(), n);
    tree.setInteractionConstraint(&constraint);
    tree.computeLeafStats(0, y.data(), nullptr);
    GrowScratch scratch;
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0,
                     0.5, scratch);
    allFeasible &= tree.interactionSubtreeIsValid(0);
    everSplit |= !tree.at(0).isBottom();
    if (tree.nodes.size() > 3) everDeep = true;
  }
  tree.setInteractionConstraint(nullptr);
  check(allFeasible, "grow-from-root: every grown tree is interaction-feasible");
  check(everSplit && everDeep,
        "grow-from-root: the constrained grow still builds non-trivial trees");
  ext_rng_destroy(rng);
  rngState = savedRngState;
  printf("ok: grow-from-root honors the interaction predicate\n");
}

// ---------------------------------------------------------------------------
// What an enumerated cut candidate stands for on a column with missing values.
//
// growTreeFromRoot enumerates one candidate per cut and gives it the prior mass
// -log(numCuts), which is the mass CGM puts on the PAIR of rules (missing left,
// missing right) of a column that routes missing values. After the draw it
// spends a fair coin to decide which of the two the winning candidate became.
// The candidate therefore stands for a two-rule group and carries that group's
// whole mass, the coin picking uniformly within it. This measures that the
// realized law is really so, by comparing three laws over the same outcome
// space {no-split} U {(cut, missing side)}:
//
//   group    the weights grow.hpp assembles, halved by the fair coin: each cut
//            candidate carrying its two-rule group's mass
//   halved   the same with each candidate instead carrying the mass of ONE of
//            its two rules, half the group's - the convention this kernel
//            implemented until 2026-08-12, kept as the arm the realized draws
//            must now REJECT
//   exact    all 2 * numCuts rules enumerated separately, each carrying its own
//            likelihood with the missing rows PLACED on the side its rule
//            names - no scan approximation. This is the law the CGM prior and
//            the leaf marginal define on the full rule set.
//
// Decision rule, fixed before the first run and re-pointed onto the group law
// when the halving was deleted: chi-square goodness of fit at alpha = 1e-3,
// df = cells - 1, over 2e5 grows. Realized root rules matching the group law
// while REJECTING the halved law confirms that a cut candidate carries its
// group's mass; matching the halved law instead would mean the halving is back
// in the weight. The exact law's own draws are the calibration control and must
// not reject. The realized draws must still reject the EXACT law: the scan
// omits the missing rows from a split's likelihood while the no-split term
// counts them, a separate inconsistency carried on its own ticket.
//
// The fixture is deliberately signal-free and gives its missing rows the node's
// mean response, so the two rules of a cut's group have near-equal exact
// likelihood: that is the regime where the group's total mass is unambiguous
// and the prior-mass question is not confounded with the scan's separate
// omission of the missing rows from the split likelihood. The residual between
// the group law and the exact law measures that second effect.

// deterministic fixture noise from a local counter generator, so the fixture
// disturbs neither the shared ext_rng stream nor the global rngState the
// seed-pinned suites depend on
double fixtureUniform(std::uint64_t& state) {
  state = state * 6364136223846793005ull + 1442695040888963407ull;
  return static_cast<double>((state >> 11) & ((1ull << 53) - 1)) *
         (1.0 / 9007199254740992.0);
}

std::vector<double> normalizedFromLogWeights(
    const std::vector<double>& logWeights) {
  double maxLogWeight = *std::max_element(logWeights.begin(), logWeights.end());
  std::vector<double> probabilities(logWeights.size());
  double sum = 0.0;
  for (size_t i = 0; i < logWeights.size(); ++i) {
    probabilities[i] = std::exp(logWeights[i] - maxLogWeight);
    sum += probabilities[i];
  }
  for (double& probability : probabilities) probability /= sum;
  return probabilities;
}

// Pearson goodness of fit against a fully specified law: df = cells - 1
double chiSquareStatistic(const std::vector<double>& counts,
                          const std::vector<double>& probabilities,
                          double numDraws) {
  double statistic = 0.0;
  for (size_t i = 0; i < counts.size(); ++i) {
    double expected = numDraws * probabilities[i];
    double deviation = counts[i] - expected;
    statistic += deviation * deviation / expected;
  }
  return statistic;
}

// Regularized upper incomplete gamma Q(a, x): the series for P below the
// crossover, the Lentz continued fraction for Q above it. Coded here because
// libR's own pchisq silently returns zero without an initialized R runtime,
// which this standalone host is not; agrees with R's to 6 figures.
double upperIncompleteGamma(double a, double x) {
  double logGammaA = std::lgamma(a);
  if (x < a + 1.0) {
    double term = 1.0 / a, sum = term;
    for (int i = 1; i < 1000; ++i) {
      term *= x / (a + i);
      sum += term;
      if (std::fabs(term) < std::fabs(sum) * 1e-16) break;
    }
    return 1.0 - sum * std::exp(-x + a * std::log(x) - logGammaA);
  }
  const double tiny = 1e-300;
  double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
  for (int i = 1; i < 1000; ++i) {
    double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    double delta = d * c;
    h *= delta;
    if (std::fabs(delta - 1.0) < 1e-16) break;
  }
  return h * std::exp(-x + a * std::log(x) - logGammaA);
}

double chiSquareUpperTail(double statistic, double df) {
  return upperIncompleteGamma(0.5 * df, 0.5 * statistic);
}

double totalVariation(const std::vector<double>& p,
                      const std::vector<double>& q) {
  double distance = 0.0;
  for (size_t i = 0; i < p.size(); ++i) distance += std::fabs(p[i] - q[i]);
  return 0.5 * distance;
}

void testOrdinalMissingRuleGroupWeight() {
  const size_t n = 64, numMissing = 8, numDraws = 200000;
  const double alpha = 1e-3, k = 2.0, sigma = 0.9;
  double residualVariance = sigma * sigma;

  std::vector<double> x(n), y(n);
  std::uint64_t generator = 20260812u;
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 5);
    y[i] = 0.6 * (2.0 * fixtureUniform(generator) - 1.0);  // no signal
  }
  for (size_t m = 0; m < numMissing; ++m) {
    x[m * (n / numMissing)] = std::nan("");
    y[m * (n / numMissing)] = 0.0;  // the node's mean response
  }

  ColumnStore store;
  store.build(x.data(), n, 1, 4);
  size_t numCuts = store.numCuts[0];
  check(store.hasMissing[0] == 1,
        "the measured fixture's only splittable column carries missing values");
  check(numCuts == 4, "the measured fixture has four ordinal cuts");

  // the node's suffstats, per code and over the missing rows
  const xint_t* codes = store.column(0);
  std::vector<ConstantLeafScanBin> bins(numCuts + 1);
  ConstantLeafScanBin missing, nodeTotal, observedTotal;
  for (size_t i = 0; i < n; ++i) {
    nodeTotal.addObservation(1.0, y[i]);
    if (codes[i] == naCode) missing.addObservation(1.0, y[i]);
    else bins[codes[i]].addObservation(1.0, y[i]);
  }
  for (const ConstantLeafScanBin& bin : bins) observedTotal.addBin(bin);
  check(missing.count == static_cast<double>(numMissing),
        "every missing row reaches the node as naCode");

  CGMTreePrior prior;
  prior.power = 8.0;  // depth-1 growth ~ 0.004: the root draw is the measurement
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  double growth = prior.growthProbability(tree, store, 0);
  double logGrowth = std::log(growth);
  // CGM's uniform over the 2 * numCuts rules, i.e. ONE rule's mass; grow.hpp's
  // logCut is twice this (the group's) and the coin halves it back
  double logRule = -std::log(static_cast<double>(2 * numCuts));

  auto marginal = [&](const ConstantLeafScanBin& bin) {
    return leaf.logIntegratedLikelihood(k, residualVariance, bin.sumWeights,
                                        bin.sumWeightedResponse);
  };

  size_t numCells = 1 + 2 * numCuts;
  std::vector<double> logGroup(numCells), logHalved(numCells);
  std::vector<double> logExact(numCells);
  double noSplit = std::log(1.0 - growth) + marginal(nodeTotal);
  logGroup[0] = logHalved[0] = logExact[0] = noSplit;

  ConstantLeafScanBin left;
  for (size_t cut = 0; cut < numCuts; ++cut) {
    left.addBin(bins[cut]);
    ConstantLeafScanBin right;
    right.count = observedTotal.count - left.count;
    right.sumWeights = observedTotal.sumWeights - left.sumWeights;
    right.sumWeightedResponse =
      observedTotal.sumWeightedResponse - left.sumWeightedResponse;
    check(left.count > 0.0 && right.count > 0.0,
          "every cut of the measured fixture is occupancy-nonempty");
    double scanned = marginal(left) + marginal(right);
    for (size_t side = 0; side < 2; ++side) {
      ConstantLeafScanBin exactLeft(left), exactRight(right);
      (side == 0 ? exactLeft : exactRight).addBin(missing);
      size_t cell = 1 + 2 * cut + side;
      logGroup[cell] = logGrowth + logRule + scanned;
      logHalved[cell] = logGrowth + logRule + scanned - std::log(2.0);
      logExact[cell] =
        logGrowth + logRule + marginal(exactLeft) + marginal(exactRight);
    }
  }

  std::vector<double> group(normalizedFromLogWeights(logGroup));
  std::vector<double> halved(normalizedFromLogWeights(logHalved));
  std::vector<double> exact(normalizedFromLogWeights(logExact));

  // realized root rules: the shipped kernel, then the exact law's own draws
  std::vector<double> shippedCounts(numCells, 0.0), exactCounts(numCells, 0.0);
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20260812u);
  GrowScratch scratch;
  for (size_t draw = 0; draw < numDraws; ++draw) {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, k,
                     sigma, scratch);
    size_t cell = 0;
    if (!tree.at(0).isBottom())
      cell = 1 + 2 * static_cast<size_t>(tree.at(0).rule.splitIndex()) +
             (tree.at(0).rule.missingGoesRight() ? 1u : 0u);
    shippedCounts[cell] += 1.0;
  }
  for (size_t draw = 0; draw < numDraws; ++draw)
    exactCounts[ext_rng_drawFromDiscreteDistribution(rng, exact.data(),
                                                     numCells)] += 1.0;
  ext_rng_destroy(rng);

  double total = static_cast<double>(numDraws);
  double df = static_cast<double>(numCells - 1);
  double x2ShippedVsGroup = chiSquareStatistic(shippedCounts, group, total);
  double x2ShippedVsHalved = chiSquareStatistic(shippedCounts, halved, total);
  double x2ShippedVsExact = chiSquareStatistic(shippedCounts, exact, total);
  double x2ExactVsExact = chiSquareStatistic(exactCounts, exact, total);
  auto pValue = [df](double statistic) {
    return chiSquareUpperTail(statistic, df);
  };

  printf("  missing-rule-group weight, %zu cells, df %.0f, %zu grows\n",
         numCells, df, numDraws);
  printf("    no-split probability: group %.5f halved %.5f exact %.5f, "
         "realized %.5f\n", group[0], halved[0], exact[0],
         shippedCounts[0] / total);
  printf("    chi2 shipped draws vs group law   %.2f (p %.3g)\n",
         x2ShippedVsGroup, pValue(x2ShippedVsGroup));
  printf("    chi2 shipped draws vs halved law  %.2f (p %.3g)\n",
         x2ShippedVsHalved, pValue(x2ShippedVsHalved));
  printf("    chi2 shipped draws vs exact law   %.2f (p %.3g)\n",
         x2ShippedVsExact, pValue(x2ShippedVsExact));
  printf("    chi2 exact draws vs exact law     %.2f (p %.3g)\n",
         x2ExactVsExact, pValue(x2ExactVsExact));
  printf("    total variation: group-exact %.5f halved-exact %.5f "
         "group-halved %.5f\n", totalVariation(group, exact),
         totalVariation(halved, exact), totalVariation(group, halved));

  // the calibration control: the chi-square is calibrated at this cell count
  // and draw count
  check(pValue(x2ExactVsExact) >= alpha,
        "the exact law's own draws match it (chi-square calibration)");

  // the measurement, pinned as measured: a cut candidate on a missing-bearing
  // column carries its two-rule group's whole prior mass and the post-draw coin
  // picks uniformly within the group, so the realized law is the group law and
  // neither the halved law nor the exact law
  check(pValue(x2ShippedVsGroup) >= alpha,
        "realized root rules match a cut candidate carrying its group's mass");
  check(pValue(x2ShippedVsHalved) < alpha,
        "the realized law is not the one that halves the group to one rule");
  check(pValue(x2ShippedVsExact) < alpha,
        "the realized law is not the exact law on the full rule set");
  // exactly a factor of two, read off the split-to-no-split odds of the two
  // reconstructed laws rather than inferred from the chi-square
  double groupOdds = (1.0 - group[0]) / group[0];
  double halvedOdds = (1.0 - halved[0]) / halved[0];
  checkNear(groupOdds / halvedOdds, 2.0, 1e-9,
            "the group's mass is exactly twice one of its two rules'");

  printf("ok: grow ordinal missing rule-group weight\n");
}

// ---------------------------------------------------------------------------
// What a categorical candidate carries, and where the orientation coin puts it.
//
// At P = R - every declared category present at the node - with no missing
// values the enumerated candidate set covers the WHOLE categorical rule space:
// the counter emits each unordered partition of the present set exactly once,
// the orientation coin picks between that partition's two direction masks, and
// no reachable position is left absent to need a coin of its own. The realized
// root rule must then follow the ideal law
//
//   no-split   (1 - growth) L(node)
//   mask m     growth / (2^R - 2) * L(left of m) L(right of m)
//
// whose rule factor is the CGM prior's OWN uniform over the non-degenerate
// direction assignments, written here from that definition rather than from
// categoricalGroupLogProbability, and whose likelihood is exact - unlike the
// ordinal path a categorical scan drops no member. One statistic falsifies the
// group mass (it sets the split-to-no-split odds), the orientation coin (a
// partition's two masks must appear equally often) and any disagreement
// between the candidate that was SCORED and the mask that was ASSEMBLED.

// per-category (count, sum w, sum wz) over the raw rows, plus the node total
void categoryBinsFromRaw(const ColumnStore& store, size_t variable,
                         const double* y, size_t n,
                         std::vector<ConstantLeafScanBin>& bins,
                         ConstantLeafScanBin& total) {
  bins.assign(static_cast<size_t>(store.numCuts[variable]), ConstantLeafScanBin{});
  for (size_t i = 0; i < n; ++i) {
    bins[store.codeAt(variable, i)].addObservation(1.0, y[i]);
    total.addObservation(1.0, y[i]);
  }
}

void testCategoricalExactDrawLaw() {
  const size_t n = 64, numDraws = 200000;
  const std::uint32_t numLevels = 4;
  const double alpha = 1e-3, k = 2.0, sigma = 0.9;
  double residualVariance = sigma * sigma;

  std::vector<double> x(n), y(n);
  std::uint64_t generator = 20260813u;
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % numLevels);
    // a mild per-category shift: enough that the seven partitions carry
    // visibly different weight, so a scrambled decode would show; not so much
    // that a cell's expected count stops supporting the chi-square
    y[i] = 0.3 * (x[i] - 1.5) + 1.2 * (fixtureUniform(generator) - 0.5);
  }
  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);
  check(store.numCuts[0] == numLevels && store.hasMissing[0] == 0 &&
          !store.columnIsPooled(0),
        "the exact-branch fixture is one inline four-level column, no missing");

  CGMTreePrior prior;
  prior.power = 8.0;  // depth-1 growth ~ 0.004: the root draw is the measurement
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  check(std::popcount(tree.reachableCategories(store, 0, 0)) ==
          static_cast<int>(numLevels),
        "every declared category reaches the root, so P = R and A = 0");

  std::vector<ConstantLeafScanBin> bins;
  ConstantLeafScanBin nodeTotal;
  categoryBinsFromRaw(store, 0, y.data(), n, bins, nodeTotal);
  auto marginal = [&](const ConstantLeafScanBin& bin) {
    return leaf.logIntegratedLikelihood(k, residualVariance, bin.sumWeights,
                                        bin.sumWeightedResponse);
  };

  // cell 0 is no-split, cell m the direction mask m over 1 .. 2^R - 2. The
  // one predictor is the only available one, so the split-variable factor is
  // log 1 and does not appear.
  double growth = prior.growthProbability(tree, store, 0);
  size_t numCells = (static_cast<size_t>(1) << numLevels) - 1;
  std::vector<double> logWeight(numCells, 0.0);
  logWeight[0] = std::log(1.0 - growth) + marginal(nodeTotal);
  double logRule =
    -std::log(std::pow(2.0, static_cast<double>(numLevels)) - 2.0);
  for (size_t mask = 1; mask < numCells; ++mask) {
    ConstantLeafScanBin left, right;
    for (std::uint32_t c = 0; c < numLevels; ++c)
      (((mask >> c) & 1u) != 0 ? right : left).addBin(bins[c]);
    logWeight[mask] =
      std::log(growth) + logRule + marginal(left) + marginal(right);
  }
  std::vector<double> law(normalizedFromLogWeights(logWeight));

  std::vector<double> counts(numCells, 0.0), control(numCells, 0.0);
  size_t offGauge = 0;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20260813u);
  GrowScratch scratch;
  for (size_t draw = 0; draw < numDraws; ++draw) {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, k,
                     sigma, scratch);
    size_t cell = 0;
    if (!tree.at(0).isBottom()) {
      std::uint64_t mask = tree.at(0).rule.categoryDirections();
      if (mask == 0 || mask >= numCells) ++offGauge;
      else cell = static_cast<size_t>(mask);
    }
    counts[cell] += 1.0;
  }
  // the calibration control: the law's own draws, at this cell and draw count
  for (size_t draw = 0; draw < numDraws; ++draw)
    control[ext_rng_drawFromDiscreteDistribution(rng, law.data(), numCells)] +=
      1.0;
  ext_rng_destroy(rng);

  double total = static_cast<double>(numDraws);
  double df = static_cast<double>(numCells - 1);
  double x2Grown = chiSquareStatistic(counts, law, total);
  double x2Control = chiSquareStatistic(control, law, total);
  double minExpected = total;
  for (double probability : law)
    if (total * probability < minExpected) minExpected = total * probability;

  printf("  categorical exact branch, P = R = %u, %zu cells, df %.0f, %zu "
         "grows\n", numLevels, numCells, df, numDraws);
  printf("    no-split probability: law %.5f realized %.5f; smallest expected "
         "cell %.0f\n", law[0], counts[0] / total, minExpected);
  printf("    chi2 grown root rules vs the exact CGM law %.2f (p %.3g)\n",
         x2Grown, chiSquareUpperTail(x2Grown, df));
  printf("    chi2 the law's own draws vs itself       %.2f (p %.3g)\n",
         x2Control, chiSquareUpperTail(x2Control, df));

  check(offGauge == 0,
        "every grown root mask is a non-degenerate assignment of the four "
        "categories");
  check(minExpected >= 20.0,
        "every cell of the exact law is measurable at this draw count");
  check(chiSquareUpperTail(x2Control, df) >= alpha,
        "the exact law's own draws match it (chi-square calibration)");
  check(chiSquareUpperTail(x2Grown, df) >= alpha,
        "grown categorical root rules follow the exact CGM law at P = R");

  printf("ok: grow categorical exact-branch draw law\n");
}

// ---------------------------------------------------------------------------
// The same measurement above the cap, where the enumeration is a sorted-prefix
// family of P - 1 candidates rather than all 2^(P-1) - 1 partitions. The
// family's TOTAL prior mass is unchanged - that is the mass repair - so each
// retained prefix carries (2^P - 2)/(2^R - 2) / (P - 1), and the orientation
// coin halves it onto the prefix's two masks. The emitted set now depends on
// the sort order, so the test rebuilds that order from the raw rows: a wrong
// key, a wrong tie-break or a prefix decode that disagrees with the assembly
// all move the realized mask frequencies.
void testCategoricalPrefixDrawLaw() {
  const size_t n = 120, numDraws = 200000;
  const std::uint32_t numLevels = 12;
  const double alpha = 1e-3, k = 2.0, sigma = 0.9;
  double residualVariance = sigma * sigma;

  std::vector<double> x(n), y(n);
  std::uint64_t generator = 20260814u;
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % numLevels);
    y[i] = 0.04 * (x[i] - 5.5) + 1.2 * (fixtureUniform(generator) - 0.5);
  }
  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 20, false, types);
  check(store.numCuts[0] == numLevels && store.hasMissing[0] == 0 &&
          numLevels > categoricalExhaustiveCap,
        "the prefix-branch fixture declares twelve levels, past the cap");

  CGMTreePrior prior;
  prior.power = 8.0;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);

  std::vector<ConstantLeafScanBin> bins;
  ConstantLeafScanBin nodeTotal;
  categoryBinsFromRaw(store, 0, y.data(), n, bins, nodeTotal);
  auto marginal = [&](const ConstantLeafScanBin& bin) {
    return leaf.logIntegratedLikelihood(k, residualVariance, bin.sumWeights,
                                        bin.sumWeightedResponse);
  };

  // the enumeration's order, rebuilt here: the singleton-category leaf
  // posterior mean, ties broken on the code
  double priorVariance = (k / leaf.scale) * (k / leaf.scale) * residualVariance;
  std::vector<std::uint32_t> order(numLevels);
  for (std::uint32_t c = 0; c < numLevels; ++c) order[c] = c;
  auto sortKey = [&](std::uint32_t c) {
    return bins[c].sumWeightedResponse / (priorVariance + bins[c].sumWeights);
  };
  std::sort(order.begin(), order.end(),
            [&](std::uint32_t a, std::uint32_t b) {
              return sortKey(a) != sortKey(b) ? sortKey(a) < sortKey(b) : a < b;
            });
  bool strictOrder = true;
  for (size_t c = 0; c + 1 < numLevels; ++c)
    strictOrder &= sortKey(order[c]) < sortKey(order[c + 1]);
  check(strictOrder, "the prefix fixture's category keys are strictly ordered");

  // cell 0 is no-split; cells 1 + 2c and 2 + 2c are prefix c's two masks
  double growth = prior.growthProbability(tree, store, 0);
  size_t numEmitted = numLevels - 1;
  size_t numCells = 1 + 2 * numEmitted;
  // the family's whole mass - the 2^P - 2 masks that induce a non-degenerate
  // partition of the present set, out of the 2^R - 2 the prior counts - spread
  // over what the enumeration keeps, then halved by the orientation coin
  double numPresent = static_cast<double>(numLevels);
  double numReachable = static_cast<double>(numLevels);  // A = 0 at the root
  double family = (std::pow(2.0, numPresent) - 2.0) /
                  (std::pow(2.0, numReachable) - 2.0);
  double logRule = std::log(family / static_cast<double>(numEmitted)) -
                   std::log(2.0);
  std::vector<double> logWeight(numCells, 0.0);
  std::vector<std::uint64_t> cellMask(numCells, 0);
  logWeight[0] = std::log(1.0 - growth) + marginal(nodeTotal);
  for (size_t c = 0; c + 1 < numLevels; ++c) {
    ConstantLeafScanBin left, right;
    std::uint64_t prefix = 0, complement = 0;
    for (size_t position = 0; position < numLevels; ++position) {
      bool held = position <= c;
      (held ? left : right).addBin(bins[order[position]]);
      (held ? prefix : complement) |= 1ull << order[position];
    }
    logWeight[1 + 2 * c] = logWeight[2 + 2 * c] =
      std::log(growth) + logRule + marginal(left) + marginal(right);
    cellMask[1 + 2 * c] = prefix;
    cellMask[2 + 2 * c] = complement;
  }
  std::vector<double> law(normalizedFromLogWeights(logWeight));

  std::vector<double> counts(numCells, 0.0), control(numCells, 0.0);
  size_t unexpectedMask = 0;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20260814u);
  GrowScratch scratch;
  for (size_t draw = 0; draw < numDraws; ++draw) {
    tree.initialize(indexBuffer.data(), n);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, k,
                     sigma, scratch);
    size_t cell = 0;
    if (!tree.at(0).isBottom()) {
      std::uint64_t mask = tree.at(0).rule.categoryDirections();
      cell = numCells;
      for (size_t candidate = 1; candidate < numCells; ++candidate)
        if (cellMask[candidate] == mask) cell = candidate;
      if (cell == numCells) { ++unexpectedMask; cell = 0; }
    }
    counts[cell] += 1.0;
  }
  for (size_t draw = 0; draw < numDraws; ++draw)
    control[ext_rng_drawFromDiscreteDistribution(rng, law.data(), numCells)] +=
      1.0;
  ext_rng_destroy(rng);

  double total = static_cast<double>(numDraws);
  double df = static_cast<double>(numCells - 1);
  double x2Grown = chiSquareStatistic(counts, law, total);
  double x2Control = chiSquareStatistic(control, law, total);
  double minExpected = total;
  for (double probability : law)
    if (total * probability < minExpected) minExpected = total * probability;

  printf("  categorical prefix branch, P = R = %u over a cap of %zu, %zu "
         "cells, df %.0f, %zu grows\n", numLevels, categoricalExhaustiveCap,
         numCells, df, numDraws);
  printf("    no-split probability: law %.5f realized %.5f; smallest expected "
         "cell %.0f\n", law[0], counts[0] / total, minExpected);
  printf("    chi2 grown root rules vs the prefix law %.2f (p %.3g)\n",
         x2Grown, chiSquareUpperTail(x2Grown, df));
  printf("    chi2 the law's own draws vs itself      %.2f (p %.3g)\n",
         x2Control, chiSquareUpperTail(x2Control, df));

  check(unexpectedMask == 0,
        "every grown root mask is a sorted prefix or its complement");
  check(minExpected >= 20.0,
        "every cell of the prefix law is measurable at this draw count");
  check(chiSquareUpperTail(x2Control, df) >= alpha,
        "the prefix law's own draws match it (chi-square calibration)");
  check(chiSquareUpperTail(x2Grown, df) >= alpha,
        "grown categorical root rules follow the mass-spread prefix law");

  printf("ok: grow categorical prefix-branch draw law\n");
}

// ---------------------------------------------------------------------------
// The same weight in closed form, with no sampling: the invariant that keeps a
// variable's realized split mass continuous across the cap.
//
// Grouping the 2^R - 2 direction masks by the partition they induce on the P
// present categories leaves 2^(A + 1) masks per unordered partition, so the
// enumerable family carries M(R, P) = (2^P - 2) 2^A / (2^R - 2) whatever the
// enumeration keeps of it. Below the cap that spreads over 2^(P-1) - 1
// candidates and must reproduce ONE partition's group mass exactly; above it,
// over the P - 1 prefixes, where only the total is pinned. The reference
// spellings are deliberately not the kernel's: expm1 rather than
// log1p/exp2, and the elementary 2^R - 2 ratio wherever that is still exact.
void testCategoricalGroupMassClosedForm() {
  auto oneMinusHalfPower = [](double count) {  // 1 - 2^(1 - count)
    return -std::expm1(-(count - 1.0) * std::log(2.0));
  };
  const std::size_t presentGrid[] = {2, 3, 4, 9, 10, 11, 12, 16, 24};
  const std::size_t absentGrid[] = {0, 1, 2, 7};
  const std::size_t wideGrid[] = {64, 1024, 65537};
  double worstTotal = 0.0, worstGroup = 0.0;
  std::size_t numGridCells = 0;
  for (std::size_t numPresent : presentGrid) {
    std::vector<std::size_t> reaches;
    for (std::size_t absent : absentGrid) reaches.push_back(numPresent + absent);
    for (std::size_t wide : wideGrid)
      if (wide >= numPresent) reaches.push_back(wide);
    for (std::size_t numReachable : reaches) {
      double emitted =
        static_cast<double>(categoricalNumEmitted(numPresent));
      double perCandidate = std::exp(
        CGMTreePrior::categoricalGroupLogProbability(numReachable, numPresent,
                                                     emitted));
      double family = oneMinusHalfPower(static_cast<double>(numPresent)) /
                      oneMinusHalfPower(static_cast<double>(numReachable));
      double relative =
        std::fabs(emitted * perCandidate - family) / family;
      if (relative > worstTotal) worstTotal = relative;
      ++numGridCells;
      // below the cap a candidate IS one partition, so it must carry the
      // 2^(A + 1) masks it groups over the 2^R - 2 the prior counts
      if (numPresent <= categoricalExhaustiveCap && numReachable <= 40) {
        double group =
          std::pow(2.0, static_cast<double>(numReachable - numPresent + 1)) /
          (std::pow(2.0, static_cast<double>(numReachable)) - 2.0);
        double relativeGroup = std::fabs(perCandidate - group) / group;
        if (relativeGroup > worstGroup) worstGroup = relativeGroup;
      }
    }
  }

  // A = 0: the family is every non-degenerate mask, so a partition's group
  // mass is exactly TWICE one mask's - the quantity
  // ruleForVariableLogProbability returns. Relative rather than bitwise
  // because one path goes through pow(2, R) - 2 and the other through
  // log1p(-exp2(1 - R)); R = 63 and 70 cross the R > 54 branch, and 70 is
  // pooled, so the cross-check meets both spellings and both mask tiers.
  const std::uint32_t crossCheckGrid[] = {3, 4, 10, 12, 20, 63, 70};
  CGMTreePrior prior;
  double worstCrossCheck = 0.0;
  for (std::uint32_t numReachable : crossCheckGrid) {
    std::vector<double> x(numReachable);
    for (std::uint32_t code = 0; code < numReachable; ++code)
      x[code] = static_cast<double>(code);
    ColumnType types[] = {ColumnType::categorical};
    ColumnStore store;
    store.build(x.data(), numReachable, 1, 10, false, types);
    check(store.numCuts[0] == numReachable &&
            store.columnIsPooled(0) == (numReachable > 63),
          "the cross-check store declares every category at its own tier");
    std::vector<index_t> indexBuffer(numReachable);
    Tree tree;
    tree.initialize(indexBuffer.data(), numReachable);
    tree.at(0).rule.variableIndex = 0;
    double oneMask = prior.ruleForVariableLogProbability(tree, store, 0);
    double group = CGMTreePrior::categoricalGroupLogProbability(
      numReachable, numReachable,
      std::ldexp(1.0, static_cast<int>(numReachable) - 1) - 1.0);
    double relative =
      std::fabs((group - oneMask) / std::log(2.0) - 1.0);
    if (relative > worstCrossCheck) worstCrossCheck = relative;
  }

  printf("  categorical group mass over %zu (R, P) cells spanning both "
         "branches and the cap\n", numGridCells);
  printf("    worst relative error: family total %.2e, per-partition group "
         "%.2e, A = 0 cross-check %.2e\n", worstTotal, worstGroup,
         worstCrossCheck);

  check(worstTotal <= 1e-12,
        "the emitted candidates' masses sum to the family's total, both "
        "branches");
  check(worstGroup <= 1e-12,
        "below the cap a candidate carries exactly its partition's group mass");
  check(worstCrossCheck <= 1e-12,
        "at A = 0 a group is exactly two of the rule prior's own masks");

  printf("ok: grow categorical group mass closed form\n");
}

// ---------------------------------------------------------------------------
// Gauge, legality and the draw count of a categorical grow, on an inline and a
// pooled fixture in a no-missing and an NA-bearing variant. Deliberately NOT
// testGrownTreeWellFormed's loop, which reads splitInterval/splitIndex
// unguarded by column type - meaningless on a rule whose bits are a mask or a
// pool offset.
//
// Per internal node: the mask is nonzero, a strict subset of the node's
// reachable set and unequal to it, and both children are occupied. Per grow:
// the generator advances by exactly one uniform per positive-growth node plus
// 1 + A per categorical rule placed (one orientation coin, then one per
// category reachable at the node but absent from its members), which is the
// only oracle that contract has. Both mask tiers route their missing
// pseudo-category as a real category - present at a node it is a histogram bin
// the partition places, absent from one it is another absent-position coin.
struct CategoricalGrowFixture {
  const char* name;
  std::uint32_t numLevels;    // declared
  std::uint32_t numObserved;  // fewer, so categories stay reachable but absent
  bool missing;
  size_t numIterations;
};

// Gauge and coin bookkeeping for one grown internal node, returning the coins
// its rule assembly must have spent.
struct GaugeTally {
  size_t nodes = 0, absentSeen = 0, absentRight = 0;
  size_t missingPresent = 0, missingAbsent = 0;
  bool inGauge = true, occupied = true, roundTrip = true;
};

size_t tallyCategoricalNode(const Tree& tree, const ColumnStore& store,
                            size_t variable, std::int32_t nodeIndex,
                            GaugeTally& tally,
                            std::vector<std::uint64_t>& reachable,
                            std::vector<std::uint8_t>& present) {
  std::int32_t column = static_cast<std::int32_t>(variable);
  size_t numWords = maskWordsForCount(store.numCuts[variable]);
  reachable.assign(numWords, 0);
  std::uint64_t inlineMask = 0;
  const std::uint64_t* mask = &inlineMask;
  if (store.columnIsPooled(variable)) {
    tree.reachableCategoriesWide(store, nodeIndex, column, reachable.data());
    mask = tree.maskWordsFor(tree.at(nodeIndex).rule);
  } else {
    reachable[0] = tree.reachableCategories(store, nodeIndex, column);
    inlineMask = tree.at(nodeIndex).rule.categoryDirections();
  }
  tally.inGauge &= !maskIsZero(mask, numWords) &&
                   maskIsSubsetOf(mask, reachable.data(), numWords) &&
                   !maskEquals(mask, reachable.data(), numWords);

  std::uint32_t missingCode = missingCategoryCode(store.numCuts[variable]);
  present.assign(static_cast<size_t>(missingCode) + 1, 0);
  size_t numPresent = 0;
  const Node& node = tree.at(nodeIndex);
  for (size_t i = node.begin; i < node.end; ++i) {
    size_t code = store.codeAt(variable, tree.indices[i]);
    if (present[code] == 0) { present[code] = 1; ++numPresent; }
  }
  size_t numReachable = maskPopcount(reachable.data(), numWords);
  for (std::uint32_t category = 0; category <= missingCode; ++category) {
    if (!maskTestBit(reachable.data(), category)) continue;
    if (present[category] != 0) {
      if (category == missingCode) ++tally.missingPresent;
      continue;
    }
    ++tally.absentSeen;
    if (maskTestBit(mask, category)) ++tally.absentRight;
    if (category == missingCode) ++tally.missingAbsent;
  }
  ++tally.nodes;
  return 1 + (numReachable - numPresent);  // orientation, then one per absent
}

void testCategoricalGrowGaugeAndCoins() {
  const CategoricalGrowFixture fixtures[] = {
    {"inline K = 6", 6, 5, false, 400},
    {"inline K = 6 with NA", 6, 5, true, 400},
    {"pooled K = 70", 70, 62, false, 200},
    {"pooled K = 70 with NA", 70, 60, true, 200},
  };
  const size_t n = 420;

  for (const CategoricalGrowFixture& fixture : fixtures) {
    std::uint64_t generator = 20260815u;
    std::vector<double> x(n * 2), y(n);
    for (size_t i = 0; i < n; ++i) {
      x[i] = static_cast<double>(i % fixture.numObserved);
      x[i + n] = static_cast<double>(i % 7) / 7.0;  // a companion ordinal
      // a per-category mean, so the recursion keeps finding structure and the
      // deep nodes - where A grows as P falls - are really reached
      y[i] = 2.0 * x[i] / static_cast<double>(fixture.numObserved - 1) - 1.0 +
             (x[i + n] < 0.5 ? -0.6 : 0.6) +
             0.4 * (fixtureUniform(generator) - 0.5);
      // the missing rows all sit high in the ordinal column, so a node under a
      // low ordinal cut holds none of them while the missing pseudo-category
      // is still reachable there: the position is then drawn by its own coin
      if (fixture.missing && i % 7 >= 5 && i % 3 == 0) x[i] = std::nan("");
    }

    ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
    std::uint32_t categoryCounts[] = {fixture.numLevels, 0};
    ColumnStore store;
    store.build(x.data(), n, 2, 10, false, types, nullptr, 0, categoryCounts);
    check(store.numCuts[0] == fixture.numLevels &&
            (store.hasMissing[0] != 0) == fixture.missing &&
            store.hasMissing[1] == 0 &&
            store.columnIsPooled(0) == (fixture.numLevels > 63),
          "the gauge fixture builds at its declared tier");

    CGMTreePrior prior;
    prior.power = 1.5;  // deeper trees, so A varies down the recursion
    ConstantGaussianLeaf leaf{0.5};
    std::vector<index_t> indexBuffer(n), rebuiltBuffer(n);
    Tree tree, rebuilt;
    ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
    ext_rng_setSeed(rng, 606u);
    GrowScratch scratch;
    GaugeTally tally;
    std::vector<std::uint64_t> reachable, masks;
    std::vector<std::uint8_t> presentScratch;
    std::vector<double> params, rebuiltParams;
    std::vector<FlatNode> flat;
    std::vector<std::int32_t> internal;
    size_t expectedUniforms = 0;
    std::clock_t started = std::clock();

    for (size_t iteration = 0; iteration < fixture.numIterations; ++iteration) {
      tree.initialize(indexBuffer.data(), n);
      tree.computeLeafStats(0, y.data(), nullptr);
      growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0,
                       0.9, scratch);
      expectedUniforms += positiveGrowthNodeCount(tree, store, prior);
      tally.occupied &= tree.bottomNodesAreOccupied();
      internal.clear();  // fillNotBottom appends
      tree.fillNotBottom(0, internal);
      for (std::int32_t node : internal) {
        std::int32_t left = tree.at(node).leftChild;
        tally.occupied &= tree.at(left).numObservations() > 0 &&
                          tree.at(left + 1).numObservations() > 0;
        // the companion ordinal carries no missing values, so an ordinal rule
        // here spends nothing beyond the node's own discrete draw
        if (tree.at(node).rule.variableIndex != 0) continue;
        expectedUniforms += tallyCategoricalNode(tree, store, 0, node, tally,
                                                 reachable, presentScratch);
      }
      // the round trip re-checks every mask against the reachable set it was
      // drawn from and refuses a rule out of gauge, pooled words included
      params.assign(tree.nodes.size(), 0.0);
      tree.flatten(store, params.data(), flat, nullptr, 1, nullptr, &masks);
      rebuilt.initialize(rebuiltBuffer.data(), n);
      tally.roundTrip &=
        rebuilt.buildFromFlat(store, flat.data(), flat.size(), rebuiltParams, 1,
                              nullptr, masks.data(), masks.size());
    }
    double elapsed = static_cast<double>(std::clock() - started) /
                     static_cast<double>(CLOCKS_PER_SEC);

    // every consumer on this path takes exactly one continuous uniform, so a
    // census of the discrete draws and the coins replays the whole stream
    size_t len = ext_rng_getSerializedStateLength(rng);
    std::vector<unsigned char> grownState(len), replayState(len);
    ext_rng_writeSerializedState(rng, grownState.data());
    ext_rng* replay = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
                                     nullptr);
    ext_rng_setSeed(replay, 606u);
    for (size_t draw = 0; draw < expectedUniforms; ++draw)
      (void) ext_rng_simulateContinuousUniform(replay);
    ext_rng_writeSerializedState(replay, replayState.data());
    ext_rng_destroy(replay);
    ext_rng_destroy(rng);

    printf("  %s: %zu rules over %zu grows, %zu absent positions drawn (%zu "
           "right), missing bin present/absent %zu/%zu, %.3f ms per grow\n",
           fixture.name, tally.nodes, fixture.numIterations, tally.absentSeen,
           tally.absentRight, tally.missingPresent, tally.missingAbsent,
           1000.0 * elapsed / static_cast<double>(fixture.numIterations));

    check(tally.nodes > 0, "the categorical fixture grows past the root");
    check(tally.inGauge, "every grown categorical mask is in gauge");
    check(tally.occupied, "every child of a grown categorical rule is occupied");
    check(tally.absentSeen > 0 && tally.absentRight > 0 &&
            tally.absentRight < tally.absentSeen,
          "reachable-but-absent categories are drawn, not pinned to one side");
    check(tally.roundTrip, "every grown tree rebuilds from its flattened form");
    check(grownState == replayState,
          "a categorical grow spends one uniform per positive-growth node plus "
          "1 + A per rule");
    if (fixture.missing)
      check(tally.missingPresent > 0 && tally.missingAbsent > 0,
            "the missing pseudo-category is placed as a bin at some nodes and "
            "drawn as an absent position at others");
  }

  printf("ok: grow categorical gauge, legality and coin count\n");
}

// The interaction constraint is the FIFTH split generator's veto (design
// "Containment"), and grow reaches it only through
// collectAvailableVariables - which is exactly why a categorical column under
// a constraint needs its own assertion: the availability predicate is the only
// thing standing between a categorical scan and a forbidden co-occurrence.
void testCategoricalGrowHonorsInteraction() {
  const size_t n = 360, p = 3;
  std::uint64_t generator = 20260816u;
  std::vector<double> x(n * p), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 5);            // categorical
    x[i + n] = static_cast<double>((i / 5) % 4);  // categorical
    x[i + 2 * n] = fixtureUniform(generator);     // ordinal
    y[i] = (x[i] < 2.0 ? -1.0 : 1.0) + (x[i + n] < 2.0 ? 0.8 : -0.8) +
           (x[i + 2 * n] < 0.5 ? 0.5 : -0.5) +
           0.2 * (fixtureUniform(generator) - 0.5);
  }
  ColumnType types[] = {ColumnType::categorical, ColumnType::categorical,
                        ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, p, 12, false, types);

  size_t forbiddenPair[] = {0, 1};  // the two categorical columns
  InteractionConstraint constraint;
  constraint.build(p, 2, forbiddenPair, 1);  // max.order 2 AND forbid(x0, x1)

  CGMTreePrior prior;
  prior.power = 1.5;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> indexBuffer(n);
  Tree tree;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 909u);
  GrowScratch scratch;
  bool allFeasible = true;
  size_t numCategorical = 0, numDeep = 0;
  std::vector<std::int32_t> internal;
  for (int iteration = 0; iteration < 300; ++iteration) {
    tree.initialize(indexBuffer.data(), n);
    tree.setInteractionConstraint(&constraint);
    tree.computeLeafStats(0, y.data(), nullptr);
    growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0,
                     0.5, scratch);
    allFeasible &= tree.interactionSubtreeIsValid(0);
    if (tree.nodes.size() > 3) ++numDeep;
    internal.clear();  // fillNotBottom appends
    tree.fillNotBottom(0, internal);
    for (std::int32_t node : internal)
      if (store.types[static_cast<size_t>(tree.at(node).rule.variableIndex)] ==
          ColumnType::categorical)
        ++numCategorical;
  }
  tree.setInteractionConstraint(nullptr);
  ext_rng_destroy(rng);

  check(allFeasible && numCategorical > 0 && numDeep > 0,
        "a constrained grow places categorical rules and every grown tree "
        "stays interaction-feasible");
  printf("ok: grow categorical under an interaction constraint (%zu rules)\n",
         numCategorical);
}

}  // namespace

void runGrowTests(ext_rng* rng) {
  testDeterminismAndDrawCount();
  testVetoDrawsNothing();
  testGrownTreeWellFormed();
  testCategoricalSplits();
  testPooledCategoricalGrow();
  testGrowThenContinue(rng);
  testGrowHonorsInteraction();
  testOrdinalMissingRuleGroupWeight();
  testCategoricalExactDrawLaw();
  testCategoricalPrefixDrawLaw();
  testCategoricalGroupMassClosedForm();
  testCategoricalGrowGaugeAndCoins();
  testCategoricalGrowHonorsInteraction();
}
