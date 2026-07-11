// Component tests for grow-from-root (grow.hpp + Chain::growForestFromRoot):
// seeded determinism and the documented per-node draw count, the well-formed /
// legal-chain-state invariants of a grown tree, the ordinal-only categorical
// contract, and grow-then-MH continuation through a live sampler.
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
              std::vector<size_t>& indexBuffer, const double* y, size_t n) {
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

  std::vector<size_t> bufA(n), bufB(n), bufC(n);
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

  std::vector<size_t> buf(n);
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
  std::vector<size_t> buf(n);
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
  std::vector<size_t> buf(n);
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

}  // namespace

void runGrowTests(ext_rng* rng) {
  testDeterminismAndDrawCount();
  testVetoDrawsNothing();
  testGrownTreeWellFormed();
  testCategoricalNeverSplit();
  testGrowThenContinue(rng);
}
