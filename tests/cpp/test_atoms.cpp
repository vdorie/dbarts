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
}
