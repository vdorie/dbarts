// Component tests for grow-from-root (grow.hpp + Chain::growForestFromRoot):
// seeded determinism and the documented per-node draw count, the well-formed /
// legal-chain-state invariants of a grown tree, the ordinal-only categorical
// contract, the prior mass an enumerated cut candidate carries on a column
// with missing values, and grow-then-MH continuation through a live sampler.
#include "common.hpp"

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

// every node the recursion draws at has positive prior growth probability; a
// no-missing build spends no missing coins, so this is the exact uniform count
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
  check(store.hasMissing[0] == 0, "no-missing build for the draw-count census");

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
  // generator and confirm it reaches the same serialized state as the grow
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

void testCategoricalNeverSplit() {
  const size_t n = 300;
  // column 0 categorical (4 levels), column 1 ordinal with the signal
  std::vector<double> x(n * 2), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    x[i + n] = runif01();
    y[i] = (x[i + n] < 0.5 ? -2.0 : 2.0) + 0.1 * (runif01() - 0.5);
  }
  ColumnType types[] = {ColumnType::categorical, ColumnType::ordinal};
  ColumnStore store;
  store.build(x.data(), n, 2, 20, false, types);

  CGMTreePrior prior;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<index_t> buf(n);
  Tree tree;
  tree.initialize(buf.data(), n);
  tree.computeLeafStats(0, y.data(), nullptr);
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 555u);
  GrowScratch scratch;
  growTreeFromRoot(store, prior, leaf, rng, tree, 0, y.data(), nullptr, 2.0, 0.9,
                   scratch);
  ext_rng_destroy(rng);

  std::vector<int32_t> internal;
  tree.fillNotBottom(0, internal);
  bool ordinalOnly = true;
  for (int32_t node : internal)
    ordinalOnly &= tree.at(node).rule.variableIndex == 1;
  check(!internal.empty(), "the ordinal signal grew the tree");
  check(ordinalOnly, "grow-from-root never splits a categorical column (v1)");
  check(tree.bottomNodesAreOccupied(), "categorical-present grow stays occupied");

  printf("ok: grow categorical never split\n");
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
// -log(numCuts) - log(2), which is the mass CGM puts on ONE rule of a column
// that routes missing values. After the draw it spends a fair coin to decide
// which of the two rules (missing left, missing right) the winning candidate
// became. The candidate therefore stands for a two-rule group while carrying
// one rule's mass. This measures whether that is really so before any weight
// changes, by comparing three laws over the same outcome space
// {no-split} U {(cut, missing side)}:
//
//   shipped  the weights grow.hpp assembles, halved by the fair coin
//   group    the same, each cut candidate instead carrying its two-rule
//            group's mass (the shipped weight times two)
//   exact    all 2 * numCuts rules enumerated separately, each carrying its own
//            likelihood with the missing rows PLACED on the side its rule
//            names - no scan approximation. This is the law the CGM prior and
//            the leaf marginal define on the full rule set.
//
// Decision rule, fixed before the first run: chi-square goodness of fit at
// alpha = 1e-3, df = cells - 1, over 2e5 grows. Realized root rules matching
// the shipped law but REJECTING the group law confirms that the shipped weight
// is one rule's, not the group's; not rejecting the group law refutes it and
// the shipped weight stands. The exact law's own draws are the calibration
// control and must not reject.
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
  // CGM's uniform over the 2 * numCuts rules; grow.hpp's logCut is this exactly
  double logRule = -std::log(static_cast<double>(2 * numCuts));

  auto marginal = [&](const ConstantLeafScanBin& bin) {
    return leaf.logIntegratedLikelihood(k, residualVariance, bin.sumWeights,
                                        bin.sumWeightedResponse);
  };

  size_t numCells = 1 + 2 * numCuts;
  std::vector<double> logShipped(numCells), logGroup(numCells);
  std::vector<double> logExact(numCells);
  double noSplit = std::log(1.0 - growth) + marginal(nodeTotal);
  logShipped[0] = logGroup[0] = logExact[0] = noSplit;

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
      logShipped[cell] = logGrowth + logRule + scanned - std::log(2.0);
      logGroup[cell] = logGrowth + logRule + scanned;
      logExact[cell] =
        logGrowth + logRule + marginal(exactLeft) + marginal(exactRight);
    }
  }

  std::vector<double> shipped(normalizedFromLogWeights(logShipped));
  std::vector<double> group(normalizedFromLogWeights(logGroup));
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
  double x2ShippedVsShipped = chiSquareStatistic(shippedCounts, shipped, total);
  double x2ShippedVsGroup = chiSquareStatistic(shippedCounts, group, total);
  double x2ShippedVsExact = chiSquareStatistic(shippedCounts, exact, total);
  double x2ExactVsExact = chiSquareStatistic(exactCounts, exact, total);
  auto pValue = [df](double statistic) {
    return chiSquareUpperTail(statistic, df);
  };

  printf("  missing-rule-group weight, %zu cells, df %.0f, %zu grows\n",
         numCells, df, numDraws);
  printf("    no-split probability: shipped %.5f group %.5f exact %.5f, "
         "realized %.5f\n", shipped[0], group[0], exact[0],
         shippedCounts[0] / total);
  printf("    chi2 shipped draws vs shipped law %.2f (p %.3g)\n",
         x2ShippedVsShipped, pValue(x2ShippedVsShipped));
  printf("    chi2 shipped draws vs group law   %.2f (p %.3g)\n",
         x2ShippedVsGroup, pValue(x2ShippedVsGroup));
  printf("    chi2 shipped draws vs exact law   %.2f (p %.3g)\n",
         x2ShippedVsExact, pValue(x2ShippedVsExact));
  printf("    chi2 exact draws vs exact law     %.2f (p %.3g)\n",
         x2ExactVsExact, pValue(x2ExactVsExact));
  printf("    total variation: shipped-exact %.5f group-exact %.5f "
         "shipped-group %.5f\n", totalVariation(shipped, exact),
         totalVariation(group, exact), totalVariation(shipped, group));

  // the calibration controls: the reconstructed shipped law is the kernel's
  // own, and the chi-square is calibrated at this cell count and draw count
  check(pValue(x2ShippedVsShipped) >= alpha,
        "realized root rules match the shipped weights as reconstructed");
  check(pValue(x2ExactVsExact) >= alpha,
        "the exact law's own draws match it (chi-square calibration)");

  // the measurement, pinned as measured: a cut candidate on a missing-bearing
  // column carries one rule's prior mass while standing for a two-rule group,
  // so the realized law is not the group law and not the exact law
  check(pValue(x2ShippedVsGroup) < alpha,
        "a cut candidate carries one rule's prior mass, not its group's");
  check(pValue(x2ShippedVsExact) < alpha,
        "the realized law is not the exact law on the full rule set");
  // exactly a factor of two, read off the split-to-no-split odds of the two
  // reconstructed laws rather than inferred from the chi-square
  double shippedOdds = (1.0 - shipped[0]) / shipped[0];
  double groupOdds = (1.0 - group[0]) / group[0];
  checkNear(groupOdds / shippedOdds, 2.0, 1e-9,
            "the group's mass is exactly twice the mass the candidate carries");

  printf("ok: grow ordinal missing rule-group weight\n");
}

}  // namespace

void runGrowTests(ext_rng* rng) {
  testDeterminismAndDrawCount();
  testVetoDrawsNothing();
  testGrownTreeWellFormed();
  testCategoricalNeverSplit();
  testGrowThenContinue(rng);
  testGrowHonorsInteraction();
  testOrdinalMissingRuleGroupWeight();
}
