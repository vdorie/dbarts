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

// ---------------------------------------------------------------------------
// Commit (iv) birth-split fuzzer (block-fusion-stage-a.md 3(iv);
// docs/plans/mutation-fuzzing.md testMutationFuzzer). Drives seeded random
// birth sequences through the ATOM path (Tree::birthStructure +
// AtomMap::splitAtom) - exactly what the wired sweep runs under the flag - and
// after every birth re-checks the invariant against an independent oracle:
//  - after each birth the two child node caches splitAtom wrote are bitwise the
//    live computeLeafStats(left) then (right) (the equivalence-critical seam);
//  - an ACCEPTED birth leaves the map matching a from-scratch buildForTree +
//    aggregateTree rebuild - numAtoms, per-leaf (A, G, Q) bitwise, member
//    slices, and atomOf all consistent with the live tree.indices;
//  - a REJECTED birth (undoBirth + undoSplit) leaves the map bitwise-identical
//    to before (the fingerprint pattern).
// Only birth is atom-routed in this commit; after a rejected birth the map is
// resynced by a fresh rebuild (the sweep does this per tree), so the sequence
// never mixes a stale post-reject map into the next accept oracle.

// A deep bitwise fingerprint of the live map state a rejected birth must leave
// untouched: numAtoms + the [0, numAtoms) SoA + atomOf.
struct MapFingerprint {
  size_t numAtoms = 0;
  std::vector<size_t> begin, end;
  std::vector<int32_t> leaf;
  std::vector<double> A, G, Q, S;
  std::vector<uint32_t> atomOf;
};

static MapFingerprint captureFingerprint(const AtomMap& m, size_t n) {
  MapFingerprint f;
  f.numAtoms = m.numAtoms;
  f.begin.assign(m.atomBegin.begin(), m.atomBegin.begin() + m.numAtoms);
  f.end.assign(m.atomEnd.begin(), m.atomEnd.begin() + m.numAtoms);
  f.leaf.assign(m.leafTuple.begin(), m.leafTuple.begin() + m.numAtoms);
  f.A.assign(m.A.begin(), m.A.begin() + m.numAtoms);
  f.G.assign(m.G.begin(), m.G.begin() + m.numAtoms);
  f.Q.assign(m.Q.begin(), m.Q.begin() + m.numAtoms);
  f.S.assign(m.S.begin(), m.S.begin() + m.numAtoms);
  f.atomOf.assign(m.atomOf.begin(), m.atomOf.begin() + n);
  return f;
}

static bool fingerprintsEqual(const MapFingerprint& a, const MapFingerprint& b) {
  if (a.numAtoms != b.numAtoms) return false;
  if (a.begin != b.begin || a.end != b.end || a.leaf != b.leaf ||
      a.atomOf != b.atomOf)
    return false;
  for (size_t i = 0; i < a.numAtoms; ++i)
    if (!bitEqual(a.A[i], b.A[i]) || !bitEqual(a.G[i], b.G[i]) ||
        !bitEqual(a.Q[i], b.Q[i]) || !bitEqual(a.S[i], b.S[i]))
      return false;
  return true;
}

// The accepted-birth oracle: the incremental map M must equal a from-scratch
// rebuild R of the same tree/residual, matched by leaf (atom ids differ under
// recycling, so match on leafTuple, not on index). Also checks atomOf ties each
// obs to the atom of its actual leaf and that obs lies in that atom's slice.
static bool mapMatchesRebuild(const AtomMap& M, AtomMap& R, size_t n) {
  if (M.numAtoms != R.numAtoms) return false;
  for (size_t c = 0; c < M.numAtoms; ++c) {
    uint32_t r = R.atomForLeaf(M.leafTuple[c]);
    if (r == AtomMap::invalidAtom) return false;
    if (!bitEqual(M.A[c], R.A[r]) || !bitEqual(M.G[c], R.G[r]) ||
        !bitEqual(M.Q[c], R.Q[r]))
      return false;
    if (M.atomBegin[c] != R.atomBegin[r] || M.atomEnd[c] != R.atomEnd[r])
      return false;
  }
  for (size_t i = 0; i < n; ++i) {
    uint32_t c = M.atomOf[i], rc = R.atomOf[i];
    if (c >= M.numAtoms || rc >= R.numAtoms) return false;
    if (M.leafTuple[c] != R.leafTuple[rc]) return false;
    bool inSlice = false;
    for (size_t k = M.atomBegin[c]; k < M.atomEnd[c]; ++k)
      if (M.members[k] == i) { inSlice = true; break; }
    if (!inSlice) return false;
  }
  return true;
}

static void testAtomBirthSplitFuzz(ext_rng* rng, int numSeeds) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  size_t numTrees = sampler->numTrees();

  CGMTreePrior prior;  // base 0.95, power 2, as the sampler default
  std::vector<double> g(nObs), w(nObs);

  bool cachesOk = true, acceptOk = true, rejectOk = true;
  int births = 0, accepts = 0, rejects = 0;

  ext_rng* op = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  for (int sd = 0; sd < numSeeds; ++sd) {
    ext_rng_setSeed(op, 1315423u + static_cast<uint_least32_t>(sd) * 2654435761u);
    bool weighted = (sd % 3) == 0;
    for (double& v : g) v = 2.0 * ext_rng_simulateContinuousUniform(op) - 1.0;
    if (weighted)
      for (double& v : w) v = 0.1 + ext_rng_simulateContinuousUniform(op);
    const double* weights = weighted ? w.data() : nullptr;

    Tree tree = sampler->chain(0).tree(sd % numTrees);
    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    map.aggregateTree(tree, g.data(), weights);
    map.writeNodeCaches(tree);

    for (int step = 0; step < 8; ++step) {
      std::vector<int32_t> leaves;
      tree.fillBottom(0, leaves);
      std::vector<int32_t> birthable;
      for (int32_t leaf : leaves)
        if (tree.hasAnyAvailableVariable(store, leaf)) birthable.push_back(leaf);
      if (birthable.empty()) break;
      int32_t L = birthable[ext_rng_simulateUnsignedIntegerUniformInRange(
        op, 0, birthable.size())];

      MapFingerprint before = captureFingerprint(map, nObs);
      Node oldNode = tree.at(L);
      size_t maskMark = tree.maskPoolMark();
      Rule rule = prior.drawRuleAndVariable(tree, store, op, L);

      // atom path: attach the pair, then source both child suffstats from the
      // split rather than the live computeLeafStats
      tree.birthStructure(L, rule);
      map.splitAtom(tree, store, g.data(), weights, L);
      ++births;

      int32_t leftNode = tree.at(L).leftChild;
      int32_t rightNode = leftNode + 1;
      bool leftEmpty = tree.at(leftNode).numObservations() == 0;
      bool rightEmpty = tree.at(rightNode).numObservations() == 0;

      // the two child caches splitAtom wrote must be bitwise the live writer's;
      // recomputing over the same slice is bitwise-idempotent, so it is safe to
      // overwrite with computeLeafStats and compare against the saved values
      double la = tree.at(leftNode).sumWeights,
             lg = tree.at(leftNode).sumWeightedResponse,
             lq = tree.at(leftNode).sumWeightedResponseSq;
      double ra = tree.at(rightNode).sumWeights,
             rg = tree.at(rightNode).sumWeightedResponse,
             rq = tree.at(rightNode).sumWeightedResponseSq;
      tree.computeLeafStats(leftNode, g.data(), weights);
      tree.computeLeafStats(rightNode, g.data(), weights);
      cachesOk &= bitEqual(la, tree.at(leftNode).sumWeights) &&
                  bitEqual(lg, tree.at(leftNode).sumWeightedResponse) &&
                  bitEqual(lq, tree.at(leftNode).sumWeightedResponseSq) &&
                  bitEqual(ra, tree.at(rightNode).sumWeights) &&
                  bitEqual(rg, tree.at(rightNode).sumWeightedResponse) &&
                  bitEqual(rq, tree.at(rightNode).sumWeightedResponseSq);
      if (!cachesOk) {
        printf("FAIL: birth-split seed %d step %d leaf %d: child cache diverged\n",
               sd, step, L);
        break;
      }

      // an empty child would be vetoed in the real sweep; force reject so the
      // tree never carries an empty leaf. Otherwise flip a coin.
      bool reject = leftEmpty || rightEmpty ||
                    ext_rng_simulateContinuousUniform(op) < 0.4;

      if (!reject) {
        AtomMap rebuild;
        rebuild.initialize(nObs);
        rebuild.buildForTree(tree, store);
        rebuild.aggregateTree(tree, g.data(), weights);
        bool ok = mapMatchesRebuild(map, rebuild, nObs);
        acceptOk &= ok;
        ++accepts;
        if (!ok) {
          printf("FAIL: birth-split seed %d step %d leaf %d: accepted map != "
                 "rebuild\n", sd, step, L);
          break;
        }
      } else {
        tree.undoBirth(L);
        tree.truncateMaskPool(maskMark);
        tree.at(L).sumWeights = oldNode.sumWeights;
        tree.at(L).sumWeightedResponse = oldNode.sumWeightedResponse;
        tree.at(L).sumWeightedResponseSq = oldNode.sumWeightedResponseSq;
        map.undoSplit();
        ++rejects;
        MapFingerprint after = captureFingerprint(map, nObs);
        bool ok = fingerprintsEqual(before, after);
        rejectOk &= ok;
        if (!ok) {
          printf("FAIL: birth-split seed %d step %d leaf %d: rejected map "
                 "mutated\n", sd, step, L);
          break;
        }
        // resync the map over the now-permuted indices before the next birth,
        // exactly as the sweep rebuilds it fresh for the next tree
        map.buildForTree(tree, store);
        map.aggregateTree(tree, g.data(), weights);
      }
      if (!cachesOk || !acceptOk || !rejectOk) break;
    }
    if (!cachesOk || !acceptOk || !rejectOk) break;
  }
  ext_rng_destroy(op);

  check(cachesOk, "birth-split child caches bitwise the live writer");
  check(acceptOk, "accepted birth map matches from-scratch rebuild");
  check(rejectOk, "rejected birth leaves map bitwise-identical");
  printf("ok: atom birth-split fuzzer (%d seeds, %d births: %d accept, %d "
         "reject)\n", numSeeds, births, accepts, rejects);
}

// ---------------------------------------------------------------------------
// Commit (v) death / change / swap oracles (block-fusion-stage-a.md 3(v), 4.1).
//
// Two facts shape the oracles below, both from the design (block-fusion.md 4.1):
//  1. DEATH is ADDITIVE. mergeAtoms sets the merged parent atom's (A, G, Q) to
//     the children's sums left-then-right (pin 5), bitwise the parent NODE cache
//     Tree::orphanChildren forms. A full-slice kernel rescan of the parent slice
//     regroups the mod-5 sum and rounds differently, so a merged atom is
//     intentionally NON-canonical vs a from-scratch aggregate. The oracle for a
//     death therefore compares TOPOLOGY against a rebuild and (A, G, Q) against
//     the maintained NODE caches (soaMatchesCaches) + the captured child sums.
//  2. A REJECTED birth leaves the tree's index segment permuted (undoBirth by
//     design, moves.hpp), so a touched atom's stored (A, G, Q) can be a valid
//     sum in a stale member ORDER. Both effects mean per-leaf (A, G, Q) are
//     verified against the node caches the atom path maintains (the ground truth
//     the equivalence anchor proves), never re-derived by a full-slice rebuild.
// TOPOLOGY (numAtoms, leaf<->atom bijection, slices, atomOf) is always canonical
// and is checked against a from-scratch buildForTree.

// Topology-only rebuild match (no A/G/Q): numAtoms equal, each live atom's leaf
// has a rebuild atom over the same slice, and atomOf ties each obs to the atom
// of its actual leaf with the obs inside that atom's slice. Matched by leaf, not
// index (recycling permutes ids).
static bool mapTopologyMatchesRebuild(const AtomMap& M, AtomMap& R, size_t n) {
  if (M.numAtoms != R.numAtoms) return false;
  for (size_t c = 0; c < M.numAtoms; ++c) {
    uint32_t r = R.atomForLeaf(M.leafTuple[c]);
    if (r == AtomMap::invalidAtom) return false;
    if (M.atomBegin[c] != R.atomBegin[r] || M.atomEnd[c] != R.atomEnd[r])
      return false;
  }
  for (size_t i = 0; i < n; ++i) {
    uint32_t c = M.atomOf[i], rc = R.atomOf[i];
    if (c >= M.numAtoms || rc >= R.numAtoms) return false;
    if (M.leafTuple[c] != R.leafTuple[rc]) return false;
    bool inSlice = false;
    for (size_t k = M.atomBegin[c]; k < M.atomEnd[c]; ++k)
      if (M.members[k] == i) { inSlice = true; break; }
    if (!inSlice) return false;
  }
  return true;
}

// Every live atom's (A, G, Q) equals its leaf's node suffstat cache bitwise -
// the invariant every structural move maintains (splitAtom/refreshSubtree write
// both; death's mergeAtoms mirrors orphanChildren's additive cache). Holds
// regardless of a death's non-canonical additive value or a stale member order.
static bool soaMatchesCaches(const AtomMap& M, const Tree& tree) {
  for (size_t c = 0; c < M.numAtoms; ++c) {
    const Node& nd(tree.at(M.leafTuple[c]));
    if (!bitEqual(M.A[c], nd.sumWeights) ||
        !bitEqual(M.G[c], nd.sumWeightedResponse) ||
        !bitEqual(M.Q[c], nd.sumWeightedResponseSq))
      return false;
  }
  return true;
}

// Focused mergeAtoms oracle: on each burned-in tree with a death-able node,
// aggregate, run the live orphanChildren + the atom mergeAtoms, and assert the
// merged atom is the children summed left-then-right (== the orphaned cache),
// the topology matches a rebuild, and the SoA tracks the caches.
static void testAtomMergeDeath(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  std::vector<double> g(nObs);
  for (double& v : g) v = 2.0 * runif01() - 1.0;

  int tested = 0;
  bool ok = true;
  for (size_t t = 0; t < sampler->numTrees() && tested < 6; ++t) {
    Tree tree = sampler->chain(0).tree(t);
    std::vector<int32_t> noGrand;
    tree.fillNoGrand(0, noGrand);
    if (noGrand.empty()) continue;

    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    map.aggregateTree(tree, g.data(), nullptr);
    map.writeNodeCaches(tree);

    int32_t P = noGrand[0];
    int32_t L = tree.at(P).leftChild, R = L + 1;
    double expA = tree.at(L).sumWeights + tree.at(R).sumWeights;
    double expG = tree.at(L).sumWeightedResponse + tree.at(R).sumWeightedResponse;
    double expQ =
      tree.at(L).sumWeightedResponseSq + tree.at(R).sumWeightedResponseSq;

    Node oldP = tree.at(P);
    tree.orphanChildren(P);          // live additive parent cache
    map.mergeAtoms(tree, P, L);      // atom SoA merge
    tree.releasePair(oldP.leftChild);

    uint32_t mc = map.atomForLeaf(P);
    ok &= mc != AtomMap::invalidAtom;
    ok &= mc != AtomMap::invalidAtom && bitEqual(map.A[mc], expA) &&
          bitEqual(map.G[mc], expG) && bitEqual(map.Q[mc], expQ);
    ok &= bitEqual(tree.at(P).sumWeights, expA);  // orphanChildren == left+right

    AtomMap rebuild;
    rebuild.initialize(nObs);
    rebuild.buildForTree(tree, store);
    rebuild.aggregateTree(tree, g.data(), nullptr);
    ok &= mapTopologyMatchesRebuild(map, rebuild, nObs);
    ok &= soaMatchesCaches(map, tree);
    ++tested;
  }
  check(tested > 0, "mergeAtoms: found death-able burned-in nodes");
  check(ok, "mergeAtoms merges children additively, topology matches rebuild");
  printf("ok: atom merge (death, %d nodes)\n", tested);
}

// Focused refreshSubtree differential: run the atom refreshSubtree and the live
// Tree::refreshSubtree from the SAME pre-move state and rule, and assert every
// subtree leaf cache is BITWISE-equal - the tightest proof the atom path's DFS
// partition + aggregation order mirror the live path (pin 6). Also asserts
// snapshot/restore returns the map to its fingerprint.
static void testAtomRefreshSubtreeDifferential(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  CGMTreePrior prior;
  ext_rng* op = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(op, 0xC0FFEEu);
  std::vector<double> g(nObs);
  for (double& v : g) v = 2.0 * ext_rng_simulateContinuousUniform(op) - 1.0;

  int tested = 0;
  bool diffOk = true, restoreOk = true;
  for (size_t t = 0; t < sampler->numTrees() && tested < 8; ++t) {
    Tree tree = sampler->chain(0).tree(t);
    std::vector<int32_t> notBottom;
    tree.fillNotBottom(0, notBottom);
    int32_t P = invalidNode;
    for (int32_t node : notBottom)
      if (tree.hasAnyAvailableVariable(store, node)) { P = node; break; }
    if (P == invalidNode) continue;

    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    map.aggregateTree(tree, g.data(), nullptr);
    map.writeNodeCaches(tree);

    Rule newRule = prior.drawRuleAndVariable(tree, store, op, P);
    MapFingerprint before = captureFingerprint(map, nObs);

    // atom path over a snapshot
    Tree::SubtreeSnapshot snap;
    tree.snapshotSubtree(P, snap);
    map.snapshotSubtree(tree, P);
    tree.at(P).rule = newRule;
    map.refreshSubtree(tree, store, g.data(), nullptr, P);

    std::vector<int32_t> leaves;
    tree.fillBottom(P, leaves);
    std::vector<double> atomA, atomG, atomQ;
    for (int32_t leaf : leaves) {
      atomA.push_back(tree.at(leaf).sumWeights);
      atomG.push_back(tree.at(leaf).sumWeightedResponse);
      atomQ.push_back(tree.at(leaf).sumWeightedResponseSq);
    }

    // restore both, assert the map is bitwise back, then run the live path
    tree.restoreSubtree(snap);
    map.restoreSubtree();
    restoreOk &= fingerprintsEqual(before, captureFingerprint(map, nObs));

    tree.at(P).rule = newRule;
    tree.refreshSubtree(store, P, g.data(), nullptr);
    std::vector<int32_t> leaves2;
    tree.fillBottom(P, leaves2);
    diffOk &= leaves2.size() == leaves.size();
    for (size_t k = 0; k < leaves2.size() && diffOk; ++k) {
      diffOk &= bitEqual(atomA[k], tree.at(leaves2[k]).sumWeights) &&
                bitEqual(atomG[k], tree.at(leaves2[k]).sumWeightedResponse) &&
                bitEqual(atomQ[k], tree.at(leaves2[k]).sumWeightedResponseSq);
    }
    ++tested;
  }
  ext_rng_destroy(op);
  check(tested > 0, "refreshSubtree: found changeable burned-in nodes");
  check(diffOk, "atom refreshSubtree caches are bitwise the live refreshSubtree");
  check(restoreOk, "refreshSubtree snapshot/restore returns the map bitwise");
  printf("ok: atom refresh-subtree differential (%d nodes)\n", tested);
}

// The oracle must not be vacuous: perturb a restored map and confirm the
// fingerprint localizes the change, then confirm the revert clears it.
static void testAtomRestoreSabotage(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  std::vector<double> g(nObs);
  for (double& v : g) v = 2.0 * runif01() - 1.0;

  Tree tree = sampler->chain(0).tree(0);
  AtomMap map;
  map.initialize(nObs);
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  map.writeNodeCaches(tree);
  check(map.numAtoms > 0, "sabotage: map built");

  MapFingerprint before = captureFingerprint(map, nObs);
  double saved = map.G[0];
  map.G[0] = saved + 1.0;                       // perturb one atom's G
  check(!fingerprintsEqual(before, captureFingerprint(map, nObs)),
        "sabotage: fingerprint localizes a perturbed G");
  map.G[0] = saved;                             // revert
  check(fingerprintsEqual(before, captureFingerprint(map, nObs)),
        "sabotage: revert clears the perturbation");

  uint32_t savedAtom = map.atomOf[0];           // perturb one atomOf entry
  map.atomOf[0] = (savedAtom + 1u) % static_cast<uint32_t>(map.numAtoms);
  check(!fingerprintsEqual(before, captureFingerprint(map, nObs)),
        "sabotage: fingerprint localizes a perturbed atomOf");
  map.atomOf[0] = savedAtom;
  check(fingerprintsEqual(before, captureFingerprint(map, nObs)),
        "sabotage: atomOf revert clears the perturbation");
  printf("ok: atom restore sabotage guard\n");
}

// ---------------------------------------------------------------------------
// Full-vocabulary move fuzzer (block-fusion-stage-a.md 3(v);
// docs/plans/mutation-fuzzing.md). Drives long seeded random sequences over the
// WHOLE move set - birth, death, change, swap, weighted - through the atom path
// (Tree op + the matching AtomMap op, exactly what the wired sweep runs), on ONE
// live map that is NEVER rebuilt between moves. After every move:
//  - accepted -> topology matches a from-scratch buildForTree rebuild, the SoA
//    tracks the node caches, and each touched (canonical) leaf's cache is
//    bitwise the live computeLeafStats;
//  - rejected -> the map is bitwise-identical to before (deep fingerprint),
//    including the specific rejected-change and rejected-swap cases.
// Only births/deaths change atom count; the map staying live across all four
// move types (no inter-move rebuild) is the property under test.

enum FullOp { FOP_BIRTH, FOP_DEATH, FOP_CHANGE, FOP_SWAP, FOP_NUM };
static const int fullOpWeight[FOP_NUM] = {5, 4, 3, 3};

static int pickFullOp(ext_rng* r) {
  int total = 0;
  for (int o = 0; o < FOP_NUM; ++o) total += fullOpWeight[o];
  int draw =
    static_cast<int>(ext_rng_simulateUnsignedIntegerUniformInRange(r, 0, total));
  for (int o = 0; o < FOP_NUM; ++o) {
    draw -= fullOpWeight[o];
    if (draw < 0) return o;
  }
  return FOP_BIRTH;
}

static bool subtreeOccupied(const Tree& tree, int32_t root) {
  std::vector<int32_t> leaves;
  tree.fillBottom(root, leaves);
  for (int32_t leaf : leaves)
    if (tree.at(leaf).numObservations() == 0) return false;
  return true;
}

// Touched-leaf independent check: each canonical leaf's cache equals the live
// computeLeafStats bitwise (recompute is idempotent over the same slice).
static bool touchedLeavesMatchKernel(Tree& tree, int32_t root, const double* g,
                                     const double* w) {
  std::vector<int32_t> leaves;
  tree.fillBottom(root, leaves);
  for (int32_t leaf : leaves) {
    double a = tree.at(leaf).sumWeights, gg = tree.at(leaf).sumWeightedResponse,
           q = tree.at(leaf).sumWeightedResponseSq;
    tree.computeLeafStats(leaf, g, w);
    if (!bitEqual(a, tree.at(leaf).sumWeights) ||
        !bitEqual(gg, tree.at(leaf).sumWeightedResponse) ||
        !bitEqual(q, tree.at(leaf).sumWeightedResponseSq))
      return false;
  }
  return true;
}

// Apply one random full-vocabulary move (birth/death/change/swap) to the live
// (tree, map) pair through the atom path and check the after-move invariant:
//  - accepted -> topology matches a from-scratch rebuild, the SoA tracks the
//    node caches, and the touched subtree's caches are bitwise the kernel;
//  - rejected -> the map is bitwise-identical to before (deep fingerprint).
// Returns false on the first invariant failure. Shared by the within-sweep
// fuzzer and the cross-sweep persistence fuzzer so both drive identical moves.
static bool applyAndCheckRandomMove(Tree& tree, AtomMap& map,
                                    const ColumnStore& store, size_t nObs,
                                    const double* g, const double* weights,
                                    CGMTreePrior& prior, ext_rng* op,
                                    int opCount[FOP_NUM], int accept[FOP_NUM],
                                    int reject[FOP_NUM], int sd, int step) {
  bool ok = true;
  MapFingerprint before = captureFingerprint(map, nObs);
  int chosen = pickFullOp(op);

  // resolve a target; fall back to birth (the always-available move) so a
  // step is rarely wasted
  std::vector<int32_t> bottoms, noGrand, notBottom, swappable;
  tree.fillBottom(0, bottoms);
  tree.fillNoGrand(0, noGrand);
  tree.fillNotBottom(0, notBottom);
  tree.fillSwappable(0, swappable);
  std::vector<int32_t> birthable;
  for (int32_t leaf : bottoms)
    if (tree.hasAnyAvailableVariable(store, leaf)) birthable.push_back(leaf);

  if (chosen == FOP_DEATH && noGrand.empty()) chosen = FOP_BIRTH;
  if (chosen == FOP_SWAP && swappable.empty()) chosen = FOP_BIRTH;
  if (chosen == FOP_CHANGE && notBottom.empty()) chosen = FOP_BIRTH;
  if (chosen == FOP_BIRTH && birthable.empty()) {
    if (!noGrand.empty()) chosen = FOP_DEATH;
    else return true;  // stump with no available split: nothing to do
  }

  if (chosen == FOP_BIRTH) {
    int32_t L = birthable[ext_rng_simulateUnsignedIntegerUniformInRange(
      op, 0, birthable.size())];
    Node oldNode = tree.at(L);
    size_t maskMark = tree.maskPoolMark();
    Rule rule = prior.drawRuleAndVariable(tree, store, op, L);
    tree.birthStructure(L, rule);
    map.splitAtom(tree, store, g, weights, L);
    ++opCount[FOP_BIRTH];

    int32_t left = tree.at(L).leftChild, right = left + 1;
    bool empty = tree.at(left).numObservations() == 0 ||
                 tree.at(right).numObservations() == 0;
    bool acc = !empty && ext_rng_simulateContinuousUniform(op) >= 0.4;
    if (acc) {
      ++accept[FOP_BIRTH];
      AtomMap rb;
      rb.initialize(nObs);
      rb.buildForTree(tree, store);
      rb.aggregateTree(tree, g, weights);
      ok &= mapTopologyMatchesRebuild(map, rb, nObs);
      ok &= soaMatchesCaches(map, tree);
      ok &= touchedLeavesMatchKernel(tree, L, g, weights);
      if (!ok) printf("FAIL: full-fuzz seed %d step %d BIRTH-accept L=%d\n",
                      sd, step, L);
    } else {
      tree.undoBirth(L);
      tree.truncateMaskPool(maskMark);
      tree.at(L).sumWeights = oldNode.sumWeights;
      tree.at(L).sumWeightedResponse = oldNode.sumWeightedResponse;
      tree.at(L).sumWeightedResponseSq = oldNode.sumWeightedResponseSq;
      map.undoSplit();
      ++reject[FOP_BIRTH];
      ok &= fingerprintsEqual(before, captureFingerprint(map, nObs));
      if (!ok) printf("FAIL: full-fuzz seed %d step %d BIRTH-reject L=%d\n",
                      sd, step, L);
    }
  } else if (chosen == FOP_DEATH) {
    int32_t P = noGrand[ext_rng_simulateUnsignedIntegerUniformInRange(
      op, 0, noGrand.size())];
    int32_t L = tree.at(P).leftChild, R = L + 1;
    double expA = tree.at(L).sumWeights + tree.at(R).sumWeights;
    double expG =
      tree.at(L).sumWeightedResponse + tree.at(R).sumWeightedResponse;
    double expQ =
      tree.at(L).sumWeightedResponseSq + tree.at(R).sumWeightedResponseSq;
    Node oldNode = tree.at(P);
    tree.orphanChildren(P);
    ++opCount[FOP_DEATH];
    // death is (almost) always acceptable (parent absorbs both children);
    // flip a coin so both the merge and its rollback are exercised
    bool acc = ext_rng_simulateContinuousUniform(op) >= 0.4;
    if (acc) {
      tree.releasePair(oldNode.leftChild);
      map.mergeAtoms(tree, P, L);
      ++accept[FOP_DEATH];
      uint32_t mc = map.atomForLeaf(P);
      ok &= mc != AtomMap::invalidAtom && bitEqual(map.A[mc], expA) &&
            bitEqual(map.G[mc], expG) && bitEqual(map.Q[mc], expQ);
      AtomMap rb;
      rb.initialize(nObs);
      rb.buildForTree(tree, store);
      rb.aggregateTree(tree, g, weights);
      ok &= mapTopologyMatchesRebuild(map, rb, nObs);
      ok &= soaMatchesCaches(map, tree);
      if (!ok) printf("FAIL: full-fuzz seed %d step %d DEATH-accept P=%d\n",
                      sd, step, P);
    } else {
      tree.at(P) = oldNode;  // reattach children, restore parent cache
      ++reject[FOP_DEATH];
      ok &= fingerprintsEqual(before, captureFingerprint(map, nObs));
      if (!ok) printf("FAIL: full-fuzz seed %d step %d DEATH-reject P=%d\n",
                      sd, step, P);
    }
  } else if (chosen == FOP_CHANGE || chosen == FOP_SWAP) {
    int32_t P;
    Rule savedRule, savedChildRule;
    int32_t swapChild = invalidNode;
    size_t maskMark = tree.maskPoolMark();
    if (chosen == FOP_CHANGE) {
      P = notBottom[ext_rng_simulateUnsignedIntegerUniformInRange(
        op, 0, notBottom.size())];
      if (!tree.hasAnyAvailableVariable(store, P)) return true;  // no proposal
    } else {
      P = swappable[ext_rng_simulateUnsignedIntegerUniformInRange(
        op, 0, swappable.size())];
      int32_t l = tree.at(P).leftChild;
      swapChild = !tree.at(l).isBottom() ? l : l + 1;  // an internal child
      savedRule = tree.at(P).rule;
      savedChildRule = tree.at(swapChild).rule;
    }

    Tree::SubtreeSnapshot snap;
    tree.snapshotSubtree(P, snap);
    map.snapshotSubtree(tree, P);
    if (chosen == FOP_CHANGE) {
      tree.at(P).rule = prior.drawRuleAndVariable(tree, store, op, P);
      ++opCount[FOP_CHANGE];
    } else {
      tree.at(P).rule = savedChildRule;
      tree.at(swapChild).rule = savedRule;
      ++opCount[FOP_SWAP];
    }
    map.refreshSubtree(tree, store, g, weights, P);

    // an emptied leaf would be vetoed in the real sweep (-1e7); force reject
    bool acc =
      subtreeOccupied(tree, P) && ext_rng_simulateContinuousUniform(op) >= 0.4;
    if (acc) {
      ++accept[chosen];
      AtomMap rb;
      rb.initialize(nObs);
      rb.buildForTree(tree, store);
      rb.aggregateTree(tree, g, weights);
      ok &= mapTopologyMatchesRebuild(map, rb, nObs);
      ok &= soaMatchesCaches(map, tree);
      ok &= touchedLeavesMatchKernel(tree, P, g, weights);
      if (!ok)
        printf("FAIL: full-fuzz seed %d step %d %s-accept P=%d\n", sd, step,
               chosen == FOP_CHANGE ? "CHANGE" : "SWAP", P);
    } else {
      tree.restoreSubtree(snap);
      map.restoreSubtree();
      tree.truncateMaskPool(maskMark);
      ++reject[chosen];
      ok &= fingerprintsEqual(before, captureFingerprint(map, nObs));
      if (!ok)
        printf("FAIL: full-fuzz seed %d step %d %s-reject P=%d\n", sd, step,
               chosen == FOP_CHANGE ? "CHANGE" : "SWAP", P);
    }
  }
  return ok;
}

static void testAtomFullMoveFuzz(ext_rng* rng, int numSeeds) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  size_t numTrees = sampler->numTrees();
  CGMTreePrior prior;
  std::vector<double> g(nObs), w(nObs);

  bool ok = true;
  int opCount[FOP_NUM] = {0, 0, 0, 0};
  int accept[FOP_NUM] = {0, 0, 0, 0};
  int reject[FOP_NUM] = {0, 0, 0, 0};

  ext_rng* op = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  for (int sd = 0; sd < numSeeds && ok; ++sd) {
    ext_rng_setSeed(op, 2246822519u + static_cast<uint_least32_t>(sd) * 40503u);
    bool weighted = (sd % 3) == 0;
    for (double& v : g) v = 2.0 * ext_rng_simulateContinuousUniform(op) - 1.0;
    if (weighted)
      for (double& v : w) v = 0.1 + ext_rng_simulateContinuousUniform(op);
    const double* weights = weighted ? w.data() : nullptr;

    Tree tree = sampler->chain(0).tree(sd % numTrees);
    AtomMap map;
    map.initialize(nObs);
    map.buildForTree(tree, store);
    map.aggregateTree(tree, g.data(), weights);
    map.writeNodeCaches(tree);

    for (int step = 0; step < 40 && ok; ++step)
      ok &= applyAndCheckRandomMove(tree, map, store, nObs, g.data(), weights,
                                    prior, op, opCount, accept, reject, sd, step);
  }
  ext_rng_destroy(op);

  check(ok, "full-vocabulary fuzzer: map stays live and correct over all moves");
  check(accept[FOP_DEATH] > 0 && accept[FOP_CHANGE] > 0 && accept[FOP_SWAP] > 0,
        "full-vocabulary fuzzer: every move type accepted at least once");
  check(reject[FOP_CHANGE] > 0 && reject[FOP_SWAP] > 0,
        "full-vocabulary fuzzer: change and swap rejected at least once");
  printf("ok: atom full-vocabulary fuzzer (%d seeds; birth %d/%d, death %d/%d, "
         "change %d/%d, swap %d/%d accept/reject)\n",
         numSeeds, accept[FOP_BIRTH], reject[FOP_BIRTH], accept[FOP_DEATH],
         reject[FOP_DEATH], accept[FOP_CHANGE], reject[FOP_CHANGE],
         accept[FOP_SWAP], reject[FOP_SWAP]);
}

// ---------------------------------------------------------------------------
// Commit (vi) cross-sweep persistence of the residual-INDEPENDENT A
// (block-fusion-stage-a.md 3(vi), 2.5; the LinearGaussianLeaf crossproduct-cache
// template). A(c) is cached across sweeps keyed by the leaf's ordered member
// list and re-validated by the same memcmp aggregateAtom already runs; G/Q/S
// stay residual-dependent and rescan every sweep. The oracles below prove: the
// persisted-A path lands node caches BITWISE-identical to a force-rebuild
// control (aCacheBypass); a membership change invalidates EXACTLY the changed
// leaves; the served A is load-bearing (a corrupted cache diverges from the
// kernel); and the cache stays correct across sweeps interleaved with the full
// move vocabulary. Membership is the memcmp's; only the weight axis (BCF/latent,
// held static in these component tests) needs the sweep's clearACache hook.

// Every live atom's SoA (its served A + this sweep's fresh G/Q) equals the leaf's
// live computeLeafStats bitwise - the ground-truth control the whole anchor
// rides on. Overwrites the node caches with the kernel values (== the SoA when
// the check passes), so the caller re-primes them from the SoA afterward.
static bool aggregateSoAMatchesKernel(const AtomMap& map, Tree& tree,
                                      const double* g, const double* w) {
  std::vector<int32_t> leaves;
  tree.fillBottom(0, leaves);
  uint32_t atomId = 0;
  for (int32_t leaf : leaves) {
    if (tree.at(leaf).numObservations() == 0) continue;
    tree.computeLeafStats(leaf, g, w);
    const Node& nd(tree.at(leaf));
    if (!bitEqual(map.A[atomId], nd.sumWeights) ||
        !bitEqual(map.G[atomId], nd.sumWeightedResponse) ||
        !bitEqual(map.Q[atomId], nd.sumWeightedResponseSq))
      return false;
    ++atomId;
  }
  return true;
}

// Pick a multi-leaf burned-in tree so the cache holds several persisting atoms.
static size_t pickMultiLeafTree(const ConstantLeafSampler& sampler,
                                size_t minLeaves) {
  for (size_t t = 0; t < sampler.numTrees(); ++t) {
    Tree probe = sampler.chain(0).tree(t);
    std::vector<int32_t> lv;
    probe.fillBottom(0, lv);
    if (lv.size() >= minLeaves) return t;
  }
  return 0;
}

// (a) Persistence-vs-control: N sweeps of aggregation over a FIXED tree (static
// weights, fresh residual each sweep) through a persistent map and through a
// force-rebuild control (aCacheBypass, recomputes A every sweep). Their per-atom
// (A, G, Q) - hence node caches, hence draws - must be bitwise-equal every sweep,
// and the persistent map must actually SERVE the cache after the warm-up sweep.
static void testAtomACachePersistence(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  Tree tree = sampler->chain(0).tree(pickMultiLeafTree(*sampler, 3));

  const int sweeps = 6;
  bool ok = true;
  size_t servedTotal = 0;
  for (int weighted = 0; weighted < 2; ++weighted) {
    std::vector<double> w(nObs);
    for (double& v : w) v = 0.1 + runif01();
    const double* weights = weighted ? w.data() : nullptr;

    AtomMap persistent, control;
    persistent.initialize(nObs);
    control.initialize(nObs);
    control.aCacheBypass = true;  // force-rebuild every sweep

    std::vector<double> g(nObs);
    for (int s = 0; s < sweeps; ++s) {
      for (double& v : g) v = 2.0 * runif01() - 1.0;  // fresh residual, weights fixed
      persistent.buildForTree(tree, store);
      persistent.aggregateTree(tree, g.data(), weights);
      control.buildForTree(tree, store);
      control.aggregateTree(tree, g.data(), weights);
      ok &= persistent.numAtoms == control.numAtoms;
      for (size_t c = 0; ok && c < persistent.numAtoms; ++c)
        ok &= bitEqual(persistent.A[c], control.A[c]) &&
              bitEqual(persistent.G[c], control.G[c]) &&
              bitEqual(persistent.Q[c], control.Q[c]);
    }
    // the control never touches the cache; the persistent map misses only the
    // warm-up sweep and serves every leaf on each of the remaining sweeps
    ok &= control.aCacheHits == 0 && control.aCacheMisses == 0;
    ok &= persistent.aCacheMisses == persistent.numAtoms;
    ok &= persistent.aCacheHits ==
          persistent.numAtoms * static_cast<size_t>(sweeps - 1);
    servedTotal += persistent.aCacheHits;
  }
  check(ok, "persisted-A path matches the force-rebuild control bitwise");
  check(servedTotal > 0, "cross-sweep cache serves A after the warm-up sweep");
  printf("ok: atom A-cache persistence vs force-rebuild control (%zu served)\n",
         servedTotal);
}

// (b) Invalidation precision: over a fixed tree, sweep 1 misses every leaf,
// sweep 2 (unchanged) serves every leaf, then a within-leaf reorder of ONE leaf
// changes its ordered member list; sweep 3 must miss EXACTLY that leaf and serve
// the K-1 others. The served caches stay bitwise the kernel throughout.
static void testAtomACacheInvalidationPrecision(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  Tree tree = sampler->chain(0).tree(pickMultiLeafTree(*sampler, 3));
  std::vector<double> g(nObs);
  for (double& v : g) v = 2.0 * runif01() - 1.0;

  AtomMap map;
  map.initialize(nObs);
  bool ok = true;

  // sweep 1: cold cache misses every leaf
  size_t m = map.aCacheMisses, h = map.aCacheHits;
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  size_t K = map.numAtoms;
  ok &= (map.aCacheMisses - m) == K && (map.aCacheHits - h) == 0;
  ok &= aggregateSoAMatchesKernel(map, tree, g.data(), nullptr);
  map.writeNodeCaches(tree);

  // sweep 2: unchanged membership serves every leaf
  m = map.aCacheMisses;
  h = map.aCacheHits;
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  ok &= (map.aCacheHits - h) == K && (map.aCacheMisses - m) == 0;
  ok &= aggregateSoAMatchesKernel(map, tree, g.data(), nullptr);
  map.writeNodeCaches(tree);

  // reorder one leaf's members: same set, different ordered list -> that leaf's
  // memcmp must fail while the untouched leaves still match
  int32_t target = invalidNode;
  {
    std::vector<int32_t> leaves;
    tree.fillBottom(0, leaves);
    for (int32_t leaf : leaves)
      if (tree.at(leaf).numObservations() >= 2) { target = leaf; break; }
  }
  check(target != invalidNode, "invalidation: found a leaf with >= 2 members");
  std::swap(tree.indices[tree.at(target).begin],
            tree.indices[tree.at(target).begin + 1]);

  // sweep 3: exactly the reordered leaf misses; the others persist
  m = map.aCacheMisses;
  h = map.aCacheHits;
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  ok &= (map.aCacheMisses - m) == 1 && (map.aCacheHits - h) == (K - 1);
  ok &= aggregateSoAMatchesKernel(map, tree, g.data(), nullptr);

  check(ok, "a membership change invalidates exactly the changed leaf's A");
  printf("ok: atom A-cache invalidation precision (%zu leaves, 1 invalidated)\n",
         K);
}

// (b') Sabotage: a corrupted cached A is SERVED on an unchanged-membership hit,
// so the served value must diverge from the live kernel - proving the memcmp-
// gated serve is load-bearing (a stale A after a real membership change would
// diverge the same way; the memcmp gate is the only guard). Restoring the entry
// returns the served value to the kernel's.
static void testAtomACacheSabotage(ext_rng* rng) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  Tree tree = sampler->chain(0).tree(pickMultiLeafTree(*sampler, 2));
  std::vector<double> g(nObs);
  for (double& v : g) v = 2.0 * runif01() - 1.0;

  AtomMap map;
  map.initialize(nObs);
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);  // warm the cache
  check(map.numAtoms > 0, "sabotage: cache warmed");

  std::vector<int32_t> leaves;
  tree.fillBottom(0, leaves);
  int32_t target = invalidNode;
  for (int32_t leaf : leaves)
    if (tree.at(leaf).numObservations() > 0) { target = leaf; break; }

  AtomMap::TreeACache& tc = map.aCacheForTree(tree);
  double good = tc.nodes[static_cast<size_t>(target)].a;
  tc.nodes[static_cast<size_t>(target)].a = good + 13.0;  // corrupt the stored A

  // unchanged membership -> hit -> serves the corrupted A into the SoA
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  uint32_t c = map.atomForLeaf(target);
  tree.computeLeafStats(target, g.data(), nullptr);
  check(c != AtomMap::invalidAtom &&
        !bitEqual(map.A[c], tree.at(target).sumWeights),
        "a corrupted cached A is served on a hit (the serve is load-bearing)");

  tc.nodes[static_cast<size_t>(target)].a = good;  // restore
  map.buildForTree(tree, store);
  map.aggregateTree(tree, g.data(), nullptr);
  c = map.atomForLeaf(target);
  tree.computeLeafStats(target, g.data(), nullptr);
  check(c != AtomMap::invalidAtom &&
        bitEqual(map.A[c], tree.at(target).sumWeights),
        "restoring the cached A returns the served value to the kernel's");
  printf("ok: atom A-cache sabotage guard\n");
}

// (c) Cross-sweep persistence stress: run the full move vocabulary on ONE
// persistent map across SWEEP boundaries. Each sweep re-aggregates against a
// fresh residual (weights fixed per seed, so the A cache legitimately persists),
// asserts every served leaf equals the kernel bitwise, then drives a batch of
// random moves with the same per-move invariants. The cache is never rebuilt or
// cleared across sweeps or moves - the property under test is that it stays
// correct while both the residual and the topology evolve.
static void testAtomACacheCrossSweepFuzz(ext_rng* rng, int numSeeds) {
  const size_t n = 200;
  std::vector<double> x, y;
  makeMutationData(x, y, n);
  auto sampler = makeBurnedInSampler(x, y, n, rng);
  const ColumnStore& store = sampler->data();
  size_t nObs = sampler->numObservations();
  size_t numTrees = sampler->numTrees();
  CGMTreePrior prior;
  std::vector<double> g(nObs), w(nObs);

  bool ok = true;
  size_t servedTotal = 0;
  int sweepsRun = 0;
  int opCount[FOP_NUM] = {0, 0, 0, 0};
  int accept[FOP_NUM] = {0, 0, 0, 0};
  int reject[FOP_NUM] = {0, 0, 0, 0};

  ext_rng* op = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  for (int sd = 0; sd < numSeeds && ok; ++sd) {
    ext_rng_setSeed(op, 3221225473u + static_cast<uint_least32_t>(sd) * 2246822519u);
    bool weighted = (sd % 3) == 0;
    if (weighted)
      for (double& v : w) v = 0.1 + ext_rng_simulateContinuousUniform(op);
    const double* weights = weighted ? w.data() : nullptr;

    Tree tree = sampler->chain(0).tree(sd % numTrees);
    AtomMap map;
    map.initialize(nObs);

    for (int sweep = 0; sweep < 6 && ok; ++sweep) {
      // fresh residual each sweep; the weight-fixed A cache persists across it
      for (double& v : g) v = 2.0 * ext_rng_simulateContinuousUniform(op) - 1.0;
      map.buildForTree(tree, store);
      map.aggregateTree(tree, g.data(), weights);
      ok &= aggregateSoAMatchesKernel(map, tree, g.data(), weights);
      if (!ok) {
        printf("FAIL: cross-sweep seed %d sweep %d aggregate != kernel\n", sd,
               sweep);
        break;
      }
      map.writeNodeCaches(tree);
      ++sweepsRun;
      for (int step = 0; step < 6 && ok; ++step)
        ok &= applyAndCheckRandomMove(tree, map, store, nObs, g.data(), weights,
                                      prior, op, opCount, accept, reject, sd,
                                      step);
    }
    servedTotal += map.aCacheHits;
  }
  ext_rng_destroy(op);

  check(ok, "A cache stays correct across sweeps interleaved with all moves");
  check(servedTotal > 0, "cross-sweep fuzzer serves persisted A across sweeps");
  check(accept[FOP_BIRTH] > 0 && accept[FOP_DEATH] > 0 &&
        accept[FOP_CHANGE] > 0 && accept[FOP_SWAP] > 0,
        "cross-sweep fuzzer accepts every move type at least once");
  printf("ok: atom A-cache cross-sweep fuzzer (%d seeds, %d sweeps, %zu served; "
         "birth %d/%d death %d/%d change %d/%d swap %d/%d)\n",
         numSeeds, sweepsRun, servedTotal, accept[FOP_BIRTH], reject[FOP_BIRTH],
         accept[FOP_DEATH], reject[FOP_DEATH], accept[FOP_CHANGE],
         reject[FOP_CHANGE], accept[FOP_SWAP], reject[FOP_SWAP]);
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
  testAtomBirthSplitFuzz(rng, 24);
  testAtomMergeDeath(rng);
  testAtomRefreshSubtreeDifferential(rng);
  testAtomRestoreSabotage(rng);
  testAtomFullMoveFuzz(rng, 24);
  testAtomACachePersistence(rng);
  testAtomACacheInvalidationPrecision(rng);
  testAtomACacheSabotage(rng);
  testAtomACacheCrossSweepFuzz(rng, 24);
}
