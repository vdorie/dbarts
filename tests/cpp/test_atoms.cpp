#include "common.hpp"

#include <cstring>
#include <vector>

// ---------------------------------------------------------------------------
// Block-fusion atom-map component tests (docs/plans/block-fusion-stage-a.md).
//
// The oracle is the LIVE tree grown by the real move code: the atom map is
// built from a tree's current leaf partition and checked against that same
// partition (member slices, topology, obs -> atom round-trip). Aggregation
// (commit ii) then checks the atom (A, G, Q) against the node cache the live
// computeLeafStats writes, BITWISE.

// Raw-bit equality of two doubles: the tightest suffstat oracle. -0.0 vs +0.0
// and NaN payloads are treated as distinct, so a differently-ordered but
// mathematically-equal sum still fails.
static bool bitEqual(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// The four build-oracle invariants against a live tree:
//  - members aliases tree.indices (DESIGN A);
//  - exactly one atom per non-empty leaf, in fillBottom order, no atom for an
//    empty leaf, atom count == non-empty-leaf count;
//  - each atom's member slice == tree.indices[leaf.begin..leaf.end) (memcmp);
//  - atomOf round-trips: obs -> atom -> the atom's slice contains obs.
static void checkBuildInvariants(const AtomMap& map, const Tree& tree, size_t n,
                                 const char* label) {
  char message[128];
  std::snprintf(message, sizeof(message), "%s: members aliases tree.indices",
                label);
  check(map.members == tree.indices, message);

  std::vector<int32_t> leaves;
  tree.fillBottom(0, leaves);
  size_t nonEmpty = 0;
  for (int32_t leaf : leaves)
    if (tree.at(leaf).numObservations() > 0) ++nonEmpty;
  std::snprintf(message, sizeof(message),
                "%s: atom count == non-empty leaf count", label);
  check(map.numAtoms == nonEmpty, message);

  uint32_t atomId = 0;
  bool topologyOk = true, sliceOk = true;
  for (int32_t leaf : leaves) {
    const Node& node(tree.at(leaf));
    if (node.numObservations() == 0) continue;
    topologyOk &= map.leafTuple[atomId] == leaf;
    topologyOk &= map.atomBegin[atomId] == node.begin;
    topologyOk &= map.atomEnd[atomId] == node.end;
    // an independent copy of the leaf's index slice, so the memcmp proves the
    // alias points at the right buffer AND that atomBegin/atomEnd are the
    // leaf's own range, not a trivially-true self comparison
    std::vector<size_t> reference(tree.indices + node.begin,
                                  tree.indices + node.end);
    sliceOk &= std::memcmp(map.members + map.atomBegin[atomId],
                           reference.data(),
                           reference.size() * sizeof(size_t)) == 0;
    ++atomId;
  }
  std::snprintf(message, sizeof(message),
                "%s: atom topology matches the leaf partition", label);
  check(topologyOk, message);
  std::snprintf(message, sizeof(message),
                "%s: atom member slice matches the leaf index slice", label);
  check(sliceOk, message);

  bool noEmptyAtom = true;
  for (size_t c = 0; c < map.numAtoms; ++c)
    noEmptyAtom &= tree.at(map.leafTuple[c]).numObservations() > 0;
  std::snprintf(message, sizeof(message), "%s: no atom maps to an empty leaf",
                label);
  check(noEmptyAtom, message);

  bool roundTrip = true;
  for (size_t i = 0; i < n && roundTrip; ++i) {
    uint32_t c = map.atomOf[i];
    if (c >= map.numAtoms) { roundTrip = false; break; }
    bool found = false;
    for (size_t k = map.atomBegin[c]; k < map.atomEnd[c]; ++k)
      if (map.members[k] == i) { found = true; break; }
    roundTrip &= found;
  }
  std::snprintf(message, sizeof(message),
                "%s: atomOf round-trips obs -> atom -> slice", label);
  check(roundTrip, message);
}

// Real grown forest: burn a small sampler in, then build the map from each
// tree's live leaf partition and check the invariants.
static void testAtomBuildBurnedIn(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();

  size_t multiLeaf = 0;
  for (size_t t = 0; t < sampler->numTrees(); ++t) {
    Tree tree = sampler->chain(0).tree(t);  // copy; indices alias the live buffer
    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    checkBuildInvariants(map, tree, nObs, "burned-in tree");
    std::vector<int32_t> leaves;
    tree.fillBottom(0, leaves);
    if (leaves.size() > 1) ++multiLeaf;
  }
  check(multiLeaf > 0, "burned-in forest grew multi-leaf trees");
  printf("ok: atom build (burned-in forest)\n");
}

// Root-only stump: one atom spanning all observations, exercising the build
// path with no partition permutation.
static void testAtomBuildStump() {
  const size_t n = 64;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) { x[i] = runif01(); y[i] = runif01(); }
  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);

  AtomMap map;
  map.initialize(n);
  map.buildForTree(tree, store);
  check(map.numAtoms == 1, "stump has one atom");
  check(map.leafTuple[0] == 0, "stump atom is the root leaf");
  check(map.atomBegin[0] == 0 && map.atomEnd[0] == n,
        "stump atom spans all observations");
  checkBuildInvariants(map, tree, n, "stump");
  printf("ok: atom build (root stump)\n");
}

// A tree with a genuinely empty leaf: a categorical mask routes a category with
// no members in the left child down one side, emptying a grandchild. The build
// must skip it - one atom per non-empty leaf only.
static void testAtomBuildEmptyLeaf() {
  const size_t n = 160;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 4);
    y[i] = runif01();
  }
  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  // root: categories {2, 3} right, {0, 1} left
  Rule root;
  root.variableIndex = 0;
  root.setCategoryDirections((1u << 2) | (1u << 3));
  tree.birth(store, 0, root, y.data(), nullptr);
  // left child ({0, 1}): send category 3 right - none present -> right empty
  Rule child;
  child.variableIndex = 0;
  child.setCategoryDirections(1u << 3);
  tree.birth(store, tree.at(0).leftChild, child, y.data(), nullptr);

  int32_t leftChild = tree.at(0).leftChild;
  int32_t emptyLeaf = tree.at(leftChild).leftChild + 1;
  check(tree.at(emptyLeaf).numObservations() == 0, "fixture built an empty leaf");

  AtomMap map;
  map.initialize(n);
  map.buildForTree(tree, store);

  std::vector<int32_t> leaves;
  tree.fillBottom(0, leaves);
  size_t emptyCount = 0, nonEmptyCount = 0;
  for (int32_t leaf : leaves)
    (tree.at(leaf).numObservations() == 0 ? emptyCount : nonEmptyCount) += 1;
  check(emptyCount >= 1, "fixture exposes the empty-leaf path");
  check(map.numAtoms == nonEmptyCount, "empty leaves excluded from atoms");
  checkBuildInvariants(map, tree, n, "empty-leaf tree");
  printf("ok: atom build (empty leaf excluded)\n");
}

// The aggregation oracle: build + aggregate the atom SoA, then run the LIVE
// computeLeafStats over the same tree and residual and assert the atom
// (A, G, Q) equals the node cache (sumWeights, sumWeightedResponse,
// sumWeightedResponseSq) for every non-empty leaf, BITWISE. This is the
// tightest oracle and the crux of Stage A. `tree` must be mutable because
// computeLeafStats writes the node caches.
static void checkAggregation(Tree& tree, const ColumnStore& store, size_t n,
                             const double* g, const double* weights,
                             const char* label) {
  AtomMap map;
  map.initialize(n);
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g, weights);  // atom path fills the SoA

  std::vector<int32_t> leaves;
  tree.fillBottom(0, leaves);
  uint32_t atomId = 0;
  bool ok = true;
  for (int32_t leaf : leaves) {
    if (tree.at(leaf).numObservations() == 0) continue;
    tree.computeLeafStats(leaf, g, weights);  // live writer -> node cache
    const Node& node(tree.at(leaf));
    ok &= bitEqual(map.A[atomId], node.sumWeights);
    ok &= bitEqual(map.G[atomId], node.sumWeightedResponse);
    ok &= bitEqual(map.Q[atomId], node.sumWeightedResponseSq);
    ++atomId;
  }
  check(ok, label);
}

static void fillRandomResidual(std::vector<double>& g) {
  for (double& v : g) v = 2.0 * runif01() - 1.0;
}

static void fillRandomWeights(std::vector<double>& w) {
  for (double& v : w) v = 0.1 + runif01();
}

// Root/stump aggregation: exercises the NON-indexed kernel dispatch (isRoot),
// weighted and unweighted.
static void testAtomAggregationStump() {
  const size_t n = 96;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = runif01();
  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);

  std::vector<double> g(n), w(n);
  fillRandomResidual(g);
  fillRandomWeights(w);
  Tree unweighted = tree, weighted = tree;
  checkAggregation(unweighted, store, n, g.data(), nullptr,
                   "stump aggregation unweighted is bitwise");
  checkAggregation(weighted, store, n, g.data(), w.data(),
                   "stump aggregation weighted is bitwise");
  printf("ok: atom aggregation (root stump, weighted + unweighted)\n");
}

// Aggregation over a tree with an empty grandchild: the empty leaf is skipped,
// the non-empty leaves aggregate through the INDEXED kernel.
static void testAtomAggregationEmptyLeaf() {
  const size_t n = 160;
  std::vector<double> x(n), y(n, 0.0);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i % 4);
  ColumnType types[] = {ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 1, 10, false, types);

  std::vector<size_t> indexBuffer(n);
  Tree tree;
  tree.initialize(indexBuffer.data(), n);
  Rule root;
  root.variableIndex = 0;
  root.setCategoryDirections((1u << 2) | (1u << 3));
  tree.birth(store, 0, root, y.data(), nullptr);
  Rule child;
  child.variableIndex = 0;
  child.setCategoryDirections(1u << 3);
  tree.birth(store, tree.at(0).leftChild, child, y.data(), nullptr);

  std::vector<double> g(n), w(n);
  fillRandomResidual(g);
  fillRandomWeights(w);
  Tree unweighted = tree, weighted = tree;
  checkAggregation(unweighted, store, n, g.data(), nullptr,
                   "empty-leaf aggregation unweighted is bitwise");
  checkAggregation(weighted, store, n, g.data(), w.data(),
                   "empty-leaf aggregation weighted is bitwise");
  printf("ok: atom aggregation (empty leaf skipped, weighted + unweighted)\n");
}

// Fuzzed aggregation: burn a forest in, then for each grown tree draw several
// fresh random residuals (and a weighted case) and assert every leaf's atom
// (A, G, Q) is bitwise the live node cache. The permuted index buffers make the
// indexed kernel path the real test.
static void testAtomAggregationGrown(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();

  std::vector<double> g(nObs), w(nObs);
  for (size_t t = 0; t < sampler->numTrees(); ++t) {
    for (int trial = 0; trial < 4; ++trial) {
      fillRandomResidual(g);
      Tree tree = sampler->chain(0).tree(t);
      checkAggregation(tree, store, nObs, g.data(), nullptr,
                       "grown-tree aggregation unweighted is bitwise");
    }
    fillRandomResidual(g);
    fillRandomWeights(w);
    Tree tree = sampler->chain(0).tree(t);
    checkAggregation(tree, store, nObs, g.data(), w.data(),
                     "grown-tree aggregation weighted is bitwise");
  }
  printf("ok: atom aggregation (fuzzed grown forest, weighted + unweighted)\n");
}

// Replicate sampleParametersAndSetFits' constant-leaf draw + scatter over a tree
// whose node caches are already filled: one normal per NON-EMPTY leaf in
// fillBottom order (empties forced to 0.0, no draw), scattered across the leaf's
// members, recording params keyed by node arena index. Reseeding the rng before
// each call fixes the RNG stream so two cache sources can be compared draw-wise.
static void drawConstantLeaf(const ConstantGaussianLeaf& leaf, ext_rng* rng,
                             Tree& tree, double k, double sigmaSq,
                             std::vector<double>& params, double* fits,
                             size_t n) {
  std::fill(fits, fits + n, 0.0);
  params.assign(tree.nodes.size(), 0.0);
  std::vector<int32_t> bottoms;
  tree.fillBottom(0, bottoms);
  for (int32_t i : bottoms) {
    const Node& node(tree.at(i));
    double param = node.numObservations() == 0
      ? 0.0
      : leaf.drawFromPosteriorForNode(rng, tree, k, sigmaSq, i);
    params[static_cast<size_t>(i)] = param;
    if (node.parent == invalidNode)
      misc_setVectorToConstant(fits, node.numObservations(), param);
    else
      misc_setIndexedVectorToConstant(fits, tree.indices + node.begin,
                                      node.numObservations(), param);
  }
}

// Commit (iii) end-to-end oracle. For each grown tree fill a fresh residual,
// then run the draw twice from the SAME reseeded RNG stream: once with the node
// caches filled by the live setNodeAverages (flag OFF), once by the atom path
// aggregateTree + writeNodeCaches (flag ON). Assert bitwise-equal per-tree fits,
// the summed total across trees, and the drawn params - the differential from
// the task. Then assert S-consistency: after setInBlockFits, S(c) equals the
// scattered fit of every member of atom c (bitwise).
static void testAtomDrawDifferential(ext_rng* rng, bool weighted) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();

  ConstantGaussianLeaf leaf;
  leaf.scale = 0.5;
  const double k = 2.0, sigmaSq = 0.7 * 0.7;

  std::vector<double> g(nObs), w(nObs);
  std::vector<double> fitsOff(nObs), fitsOn(nObs), paramsOff, paramsOn;
  std::vector<double> totalOff(nObs, 0.0), totalOn(nObs, 0.0);

  bool fitsOk = true, paramsOk = true, sOk = true;
  for (size_t t = 0; t < sampler->numTrees(); ++t) {
    fillRandomResidual(g);
    const double* weights = nullptr;
    if (weighted) {
      fillRandomWeights(w);
      weights = w.data();
    }
    Tree tree = sampler->chain(0).tree(t);
    uint_least32_t seed = 918273u + static_cast<uint_least32_t>(t);

    // flag OFF: the live per-leaf writer fills the node caches
    tree.setNodeAverages(g.data(), weights);
    ext_rng_setSeed(rng, seed);
    drawConstantLeaf(leaf, rng, tree, k, sigmaSq, paramsOff, fitsOff.data(),
                     nObs);

    // flag ON: the atom path re-sources the SAME caches; same seed => same draw
    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    map.aggregateTree(tree, g.data(), weights);
    map.writeNodeCaches(tree);
    ext_rng_setSeed(rng, seed);
    drawConstantLeaf(leaf, rng, tree, k, sigmaSq, paramsOn, fitsOn.data(), nObs);

    for (size_t i = 0; i < nObs; ++i) {
      fitsOk &= bitEqual(fitsOff[i], fitsOn[i]);
      totalOff[i] += fitsOff[i];
      totalOn[i] += fitsOn[i];
    }
    for (size_t i = 0; i < paramsOff.size(); ++i)
      paramsOk &= bitEqual(paramsOff[i], paramsOn[i]);

    map.setInBlockFits(paramsOn);
    for (size_t c = 0; c < map.numAtoms; ++c)
      for (size_t m = map.atomBegin[c]; m < map.atomEnd[c]; ++m)
        sOk &= bitEqual(map.S[c], fitsOn[map.members[m]]);
  }
  bool totalOk = true;
  for (size_t i = 0; i < nObs; ++i) totalOk &= bitEqual(totalOff[i], totalOn[i]);

  const char* tag = weighted ? "weighted" : "unweighted";
  char message[128];
  std::snprintf(message, sizeof(message),
                "atom draw differential per-tree fits bitwise (%s)", tag);
  check(fitsOk, message);
  std::snprintf(message, sizeof(message),
                "atom draw differential total fits bitwise (%s)", tag);
  check(totalOk, message);
  std::snprintf(message, sizeof(message),
                "atom draw differential params bitwise (%s)", tag);
  check(paramsOk, message);
  std::snprintf(message, sizeof(message),
                "atom S(c) == member fit bitwise (%s)", tag);
  check(sOk, message);
  printf("ok: atom draw differential + S-consistency (%s)\n", tag);
}

void runAtomsTests(ext_rng* rng) {
  testAtomBuildStump();
  testAtomBuildEmptyLeaf();
  testAtomBuildBurnedIn(rng);
  testAtomAggregationStump();
  testAtomAggregationEmptyLeaf();
  testAtomAggregationGrown(rng);
  testAtomDrawDifferential(rng, false);
  testAtomDrawDifferential(rng, true);
}
