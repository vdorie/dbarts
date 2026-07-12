#ifndef BARTCORE_ATOMS_HPP
#define BARTCORE_ATOMS_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <vector>

#include <misc/stats.h>

#include "data.hpp"
#include "tree.hpp"

/// Block-fusion Stage A build knob (block-fusion-stage-a.md 5.1). A compile-time
/// switch, DEFAULT ON as of the Stage-A milestone: the shipped constant-leaf
/// steady-state sweep sources each sweep's per-leaf sufficient statistics
/// through the b=1 atom path (the fused buildAggregateWrite pass) instead of the
/// live setNodeAverages/computeLeafStats writer. The atom aggregation is bitwise
/// the same kernel over the same member order and lands the same values in the
/// same node cache, so every downstream reader and the draws are unchanged
/// (equivalence 22/22 exact). The path is constant-leaf + steady-state-run()
/// only: vector-param (linear) and function-param (GP) leaves and the
/// grow-from-root warm start stay on the legacy writer (see Chain::leafIsConstant
/// and the warm-start loop).
///
/// The knob stays so a build can define BARTCORE_BLOCK_FUSION=0 (e.g.
/// -DBARTCORE_BLOCK_FUSION=0) for a one-line abort back to the legacy writer,
/// AND it is Stage B's b-dispatch seam. Both writers stay compiled regardless
/// of the switch, so tests/cpp can drive them side by side.
#ifndef BARTCORE_BLOCK_FUSION
#define BARTCORE_BLOCK_FUSION 1
#endif

namespace bartcore {

/// Block-fusion atom map for the constant-leaf Gaussian sweep
/// (docs/design/block-fusion.md; docs/plans/block-fusion-stage-a.md).
///
/// An ATOM is a maximal set of observations that every tree of a block routes
/// to the same leaf tuple. At the b=1 special case wired here an atom is
/// exactly one non-empty leaf of one tree, its member slice is that leaf's
/// [begin, end) range of the tree's index buffer, and its per-atom sufficient
/// statistics are the leaf's node-cache values. This b=1 map is the bitwise
/// de-risking anchor: it reproduces the current node suffstat cache exactly by
/// aggregating through the atom vocabulary, before Stage B regroups the sum at
/// b > 1.
///
/// This struct is NOT templated on the leaf model (constant-leaf only): it
/// operates on plain doubles, a Tree, and the ColumnStore. One AtomMap is owned
/// per Forest and reused across sweeps.
///
/// DESIGN A (block-fusion-stage-a.md 2.4): at b=1 `members` ALIASES the block
/// tree's `tree.indices` rather than owning a second n-sized buffer, so the map
/// adds only atomOf (4n bytes) plus the tiny per-atom SoA and no member-copy
/// pass. The owned per-block buffer arrives with Stage B.
struct AtomMap {
  // --- topology (persists across sweeps; patched by moves in later commits) ---
  /// obs id -> atom id. Size n. Every observation lands in exactly one
  /// non-empty leaf, so every entry is written by buildForTree.
  std::vector<std::uint32_t> atomOf;
  /// per-atom slice [atomBegin, atomEnd) into `members`.
  std::vector<std::size_t> atomBegin;
  std::vector<std::size_t> atomEnd;
  /// per-atom leaf id in the block tree (b=1: one arena index per atom).
  std::vector<std::int32_t> leafTuple;

  // --- per-atom sufficient-statistic SoA ---
  /// A: weight mass (sum w; count when unweighted). G: sum w*g. Q: sum w*g^2
  /// (carried at b=1 for the bitwise anchor; dropped at Stage B). S: in-block
  /// fit (b=1: the atom's leaf mean; filled at the draw in a later commit).
  std::vector<double> A, G, Q, S;

  // --- member buffer ---
  /// At b=1 (DESIGN A) `members` aliases the block tree's index buffer. At
  /// blockSize > 1 it points at the OWNED buffer atomMembersOwned (below): the
  /// joint map is a permutation of 0..n-1 grouped by atom that no single tree
  /// owns (block-fusion.md 2.3). buildForBlock and the b>1 move kernels keep it
  /// grouped by atom; the b=1 path never touches atomMembersOwned.
  std::size_t* members = nullptr;
  std::vector<std::size_t> atomMembersOwned;

  /// Block width b. Stays 1 for the shipped path (leafTuple then holds one id
  /// per atom, indexed as leafTuple[c]); at b>1 leafTuple holds b ids per atom
  /// in row-major order (leafOf(c, j) = leafTuple[c * blockWidth + j], the
  /// atom's leaf in block tree j).
  std::size_t blockWidth = 1;

  /// Live atom count == non-empty-leaf count at b=1.
  std::size_t numAtoms = 0;

  /// Maintain the obs->atom map `atomOf`. This is Stage-B move bookkeeping: at
  /// the b=1 shipped path NO consumer reads atomOf (every node-cache value is
  /// produced from leafTuple + the A/G/Q SoA, and the leaf<->atom lookup scans
  /// leafTuple), so the shipped sweep leaves it OFF and buildForTree becomes
  /// O(#leaves) rather than paying an O(n) obs->atom scatter every tree, and the
  /// moves skip their O(n_leaf) atomOf writes. DEFAULT ON so the component tests
  /// and the fuzzers - which assert atomOf round-trips and rebuild-equality -
  /// and Stage B get the full map with no per-test opt-in; Chain::run turns it
  /// off for the shipped b=1 path. atomOf is inert either way (it feeds no draw),
  /// so toggling it moves no byte of the equivalence anchor.
  bool trackAtomOf = true;

  /// Sentinel for "no atom" (an empty leaf holds none).
  static constexpr std::uint32_t invalidAtom = 0xFFFFFFFFu;

  /// Reused fillBottom scratch so aggregation does not allocate per call.
  std::vector<std::int32_t> bottomScratch;

  // --- cross-sweep residual-INDEPENDENT persistence: the A cache ---
  //
  // A(c), a leaf's weight mass, is residual-INDEPENDENT within a fixed
  // weighting, so it persists across sweeps EXACTLY like LinearGaussianLeaf's
  // crossproduct cache (block-fusion-stage-a.md 2.5; model.hpp:458-539).
  // aggregateAtom re-validates on lookup by comparing the cached ordered member
  // list against the leaf's live members[begin..end) (std::memcmp): a structural
  // move that alters membership - or a rejected move whose rollback restores the
  // order - falls out of the compare with no per-move hook (model.hpp:459-462).
  // A is WEIGHT-dependent, so the sweep drops the cache wholesale on any working
  // weight change (clearACache, the atom analog of invalidateStatistics); the
  // member-list memcmp handles the membership axis. G/Q/S stay residual-DEPENDENT
  // and are rescanned from the current treeY every block entry, exactly as the
  // U'WU cache always rescans U'Wz / z'Wz.
  //
  /// One leaf's cached weight mass, tagged with the ordered member list that
  /// produced it; the memcmp against members[begin..end) is the sole coherence
  /// gate (rollback-stable, no per-move hook).
  struct CachedLeafA {
    std::vector<std::size_t> members;
    double a = 0.0;
  };
  /// Per-tree A cache, arena-indexed by node id (like TreeStatisticsCache).
  struct TreeACache {
    const Tree* tree = nullptr;
    std::vector<CachedLeafA> nodes;
  };
  std::vector<TreeACache> aCaches_;
  /// The b>1 A cache, keyed by ATOM id rather than (tree, node): a joint atom
  /// is not one tree's leaf, and several atoms share a block tree's leaf, so
  /// the b=1 (tree, node) key is ambiguous at b>1. The memcmp is the coherence
  /// gate exactly as above, but it validates against the atom's OWNED member
  /// slice (members + atomBegin[c], the atomMembersOwned buffer) rather than
  /// tree.indices -- at b>1 the members live there, so comparing tree.indices
  /// would validate the wrong buffer entirely.
  std::vector<CachedLeafA> aCacheByAtom_;
  /// Byte ceiling over cached member lists; when spent, further leaves rescan
  /// (still correct, just uncached), so the cache can never bloat unbounded.
  std::size_t aCacheUsedBytes_ = 0;
  std::size_t aCacheBudgetBytes_ = static_cast<std::size_t>(256) << 20;
  /// Skip the cross-sweep A cache entirely: aggregation recomputes A every sweep
  /// and stores nothing, so lookupOrStoreA returns the freshly kernelled mass
  /// with NO member-list memcmp and NO store. This is both the persistence-vs-
  /// control test's force-rebuild baseline AND the shipped Stage-A setting: at
  /// b=1 the monolithic suffstat kernel rescans G/Q every sweep and drags A
  /// along, so the cache saves no work and only adds an O(n) per-leaf memcmp
  /// pass - pure overhead. Chain::run sets this on the shipped b=1 path; the
  /// A-cache tests leave it off to exercise the Stage-B persistence machinery.
  /// Served A is byte-identical to the fresh mass, so the toggle moves no draw.
  bool aCacheBypass = false;
  /// Serve/rebuild counters (test instrumentation, negligible in the sweep): a
  /// hit served a persisted A, a miss recomputed and stored.
  std::size_t aCacheHits = 0, aCacheMisses = 0;

  /// Rollback record for one pending atom split (mirrors Tree::SubtreeSnapshot
  /// for the SoA). splitAtom fills it before mutating; undoSplit restores the
  /// parent atom bitwise and drops the child atom so a rejected birth leaves
  /// the map byte-for-byte as it was.
  struct SplitUndo {
    bool active = false;        // a split is pending rollback
    bool hasParent = false;     // the split leaf held an atom (always true for
                                // a chain-state birth target; empty leaves are
                                // vetoed out of the state)
    bool createdAtom = false;   // numAtoms was incremented (both children full)
    std::uint32_t parentAtom = invalidAtom;
    std::uint32_t rightAtom = invalidAtom;
    double A = 0.0, G = 0.0, Q = 0.0, S = 0.0;  // parent atom's pre-split stats
    std::size_t begin = 0, end = 0;             // parent atom's pre-split slice
    std::int32_t leaf = 0;                      // parent atom's pre-split leaf
    std::size_t numAtomsBefore = 0;
    std::vector<std::size_t> movedMembers;      // obs whose atomOf changed
  } undo_;

  /// Size the n-indexed state once at forest init. The per-atom SoA is grown
  /// lazily by buildForTree, so a fresh map starts empty there.
  void initialize(std::size_t numObservations) {
    atomOf.assign(numObservations, 0u);
    atomBegin.clear();
    atomEnd.clear();
    leafTuple.clear();
    A.clear();
    G.clear();
    Q.clear();
    S.clear();
    members = nullptr;
    atomMembersOwned.clear();
    blockWidth = 1;
    numAtoms = 0;
    clearACache();
    aCacheHits = 0;
    aCacheMisses = 0;
  }

  /// Grow every per-atom vector to hold k atoms. leafTuple carries blockWidth
  /// ids per atom (1 at b=1, so this is unchanged for the shipped path).
  void ensureAtomCapacity(std::size_t k) {
    atomBegin.resize(k);
    atomEnd.resize(k);
    leafTuple.resize(k * blockWidth);
    A.resize(k);
    G.resize(k);
    Q.resize(k);
    S.resize(k);
  }

  /// Build the b=1 map from a tree's leaf partition: one atom per NON-EMPTY
  /// leaf, taken in fillBottom (DFS, left child before right) order so atom ids
  /// follow the same leaf order the draw consumes RNG in. Each atom's
  /// [atomBegin, atomEnd) is the leaf's [begin, end), leafTuple is the leaf id,
  /// and atomOf[i] is the atom of every member obs i. `members` aliases
  /// tree.indices (DESIGN A). A/G/Q/S are left for aggregation.
  void buildForTree(const Tree& tree, const ColumnStore& data) {
    members = tree.indices;
    std::size_t n = data.numObservations;
    if (atomOf.size() != n) atomOf.assign(n, 0u);

    bottomScratch.clear();
    tree.fillBottom(0, bottomScratch);

    numAtoms = 0;
    for (std::int32_t leaf : bottomScratch)
      if (tree.at(leaf).numObservations() > 0) ++numAtoms;
    ensureAtomCapacity(numAtoms);

    std::uint32_t atomId = 0;
    for (std::int32_t leaf : bottomScratch) {
      const Node& node(tree.at(leaf));
      if (node.numObservations() == 0) continue;
      atomBegin[atomId] = node.begin;
      atomEnd[atomId] = node.end;
      leafTuple[atomId] = leaf;
      // The obs->atom scatter is the map's only O(n) build cost and no b=1
      // consumer reads it (Stage-B bookkeeping): the shipped sweep leaves it off
      // so the build is O(#leaves).
      if (trackAtomOf)
        for (std::size_t k = node.begin; k < node.end; ++k)
          atomOf[members[k]] = atomId;
      ++atomId;
    }
  }

  /// Aggregate one atom's (A, G, Q) = (sum w, sum w*g, sum w*g^2) over its
  /// member slice. This dispatches to EXACTLY the misc kernel
  /// Tree::computeLeafStats would pick for the atom's leaf - weighted vs
  /// unweighted, and the root/stump NON-indexed kernel vs the non-root INDEXED
  /// kernel - over the SAME member order, so the result is bitwise-identical to
  /// the node cache by construction (block-fusion-stage-a.md 1.3, pins 2-3).
  /// Never hand-roll the reduction: a differently-ordered but mathematically
  /// equal sum rounds differently and breaks the b=1 anchor.
  void aggregateAtom(const Tree& tree, const double* g, const double* weights,
                     std::uint32_t atomId) {
    std::int32_t leafId = leafTuple[atomId];
    std::size_t begin = atomBegin[atomId];
    std::size_t length = atomEnd[atomId] - begin;
    bool isRoot = tree.at(leafId).parent == invalidNode;

    double a, gMass, q;
    if (weights == nullptr) {
      if (isRoot)
        misc_computeSufficientStatisticsFast(g, length, &a, &gMass, &q);
      else
        misc_computeIndexedSufficientStatisticsFast(g, members + begin, length,
                                                    &a, &gMass, &q);
    } else {
      if (isRoot)
        misc_computeWeightedSufficientStatisticsFast(g, length, weights, &a,
                                                     &gMass, &q);
      else
        misc_computeIndexedWeightedSufficientStatisticsFast(
          g, members + begin, length, weights, &a, &gMass, &q);
    }
    // A (weight mass) is residual-independent: serve it from the cross-sweep
    // cache when the leaf's ordered member list is unchanged, else record the
    // fresh mass. The served A is byte-identical to `a` (same members, same
    // weighting, same kernel), so it moves no draw byte; G/Q are residual-
    // dependent and always come from this sweep's scan above.
    A[atomId] = lookupOrStoreA(tree, leafId, begin, length, a);
    G[atomId] = gMass;
    Q[atomId] = q;
  }

  /// Aggregate every atom, visiting leaves in fillBottom (DFS, left-first)
  /// order (pin 1) so the aggregation order matches the draw's leaf/RNG order
  /// at b=1. Atom ids were assigned in this same order by buildForTree, so the
  /// running counter selects the atom that owns the current leaf. Fills the
  /// SoA only; the node cache stays the live writer's until a later commit
  /// wires the source over.
  void aggregateTree(const Tree& tree, const double* g, const double* weights) {
    bottomScratch.clear();
    tree.fillBottom(0, bottomScratch);
    std::uint32_t atomId = 0;
    for (std::int32_t leaf : bottomScratch) {
      if (tree.at(leaf).numObservations() == 0) continue;
      aggregateAtom(tree, g, weights, atomId);
      ++atomId;
    }
  }

  /// The per-tree A cache slot (created on first use), keyed by tree identity
  /// like LinearGaussianLeaf::statisticsCacheForTree. The forest owns a stable
  /// tree per slot, so the linear scan is over the tiny forest-tree count.
  TreeACache& aCacheForTree(const Tree& tree) {
    for (TreeACache& c : aCaches_)
      if (c.tree == &tree) return c;
    aCaches_.emplace_back();
    aCaches_.back().tree = &tree;
    return aCaches_.back();
  }

  /// Record a freshly kernelled weight mass for one leaf, subject to the byte
  /// budget (mirrors storeCrossproduct). Over budget: free the entry so the next
  /// lookup misses and rescans - never a stale serve, even if accounting drifts.
  void storeLeafA(CachedLeafA& entry, std::size_t begin, std::size_t length,
                  double a) {
    std::size_t oldBytes = entry.members.size() * sizeof(std::size_t);
    std::size_t newBytes = length * sizeof(std::size_t);
    if (aCacheUsedBytes_ - oldBytes + newBytes > aCacheBudgetBytes_) {
      if (!entry.members.empty()) {
        aCacheUsedBytes_ -= oldBytes;
        entry.members.clear();
      }
      return;
    }
    aCacheUsedBytes_ += newBytes - oldBytes;
    entry.members.assign(members + begin, members + begin + length);
    entry.a = a;
  }

  /// Serve leaf `leafId`'s weight mass A from the cross-sweep cache when its
  /// ordered member list is unchanged - the same memcmp gate
  /// LinearGaussianLeaf::lookupCrossproduct uses against tree.indices[begin..end)
  /// - else record `freshA` and return it. On a hit the returned value is
  /// byte-identical to `freshA`, since the entry was stored by the identical
  /// kernel over the same members and weighting; a weight change drops the whole
  /// cache first (clearACache), so a hit can never serve a stale-weight mass.
  /// Bypassed (the force-rebuild control) returns freshA without touching the
  /// cache. patch-on-accept and rollback-stability are entirely the memcmp's: an
  /// accepted move permutes/reslices the members and misses; a rejected move that
  /// restores the order hits; there is no per-move cache hook.
  double lookupOrStoreA(const Tree& tree, std::int32_t leafId, std::size_t begin,
                        std::size_t length, double freshA) {
    if (aCacheBypass) return freshA;
    TreeACache& cache = aCacheForTree(tree);
    std::size_t index = static_cast<std::size_t>(leafId);
    if (index >= cache.nodes.size()) cache.nodes.resize(index + 1);
    CachedLeafA& entry = cache.nodes[index];
    if (entry.members.size() == length &&
        std::memcmp(entry.members.data(), members + begin,
                    length * sizeof(std::size_t)) == 0) {
      ++aCacheHits;
      return entry.a;
    }
    ++aCacheMisses;
    storeLeafA(entry, begin, length, freshA);
    return freshA;
  }

  /// Drop the whole cross-sweep A cache. A is weight-dependent, so the sweep
  /// calls this whenever the working weights change - a BCF forest reweights by
  /// m^2 each sweep, a Polya-Gamma family redraws its weights, setWeights swaps
  /// them - the atom analog of LinearGaussianLeaf::invalidateStatistics.
  /// Membership changes need no clear; the per-leaf memcmp catches them.
  void clearACache() {
    aCaches_.clear();
    aCacheByAtom_.clear();
    aCacheUsedBytes_ = 0;
  }

  /// Write the aggregated (A, G, Q) SoA into each leaf's Node suffstat cache
  /// {sumWeights, sumWeightedResponse, sumWeightedResponseSq} - the single seam
  /// every constant-leaf consumer (move scoring, the leaf draw, k/scatter)
  /// reads. Walks fillBottom in the SAME DFS (left-child-first) order
  /// aggregateTree filled the SoA in, so the running atomId over non-empty
  /// leaves selects the SoA entry that owns the current leaf (pin 1). An empty
  /// leaf holds no atom and gets a forced 0.0/0.0/0.0 - byte-identical to what
  /// computeLeafStats writes over a zero-length slice (the misc kernels return
  /// +0.0 for length 0), so the whole tree's caches equal setNodeAverages
  /// BITWISE (block-fusion-stage-a.md 1.2, pin 8). aggregateTree must have run
  /// against the current residual first.
  void writeNodeCaches(Tree& tree) {
    bottomScratch.clear();
    tree.fillBottom(0, bottomScratch);
    std::uint32_t atomId = 0;
    for (std::int32_t leaf : bottomScratch) {
      Node& node(tree.at(leaf));
      if (node.numObservations() == 0) {
        node.sumWeights = 0.0;
        node.sumWeightedResponse = 0.0;
        node.sumWeightedResponseSq = 0.0;
        continue;
      }
      node.sumWeights = A[atomId];
      node.sumWeightedResponse = G[atomId];
      node.sumWeightedResponseSq = Q[atomId];
      ++atomId;
    }
  }

  /// Shipped b=1 per-tree entry: buildForTree + aggregateTree + writeNodeCaches
  /// fused into ONE fillBottom traversal, so the per-tree cost matches legacy
  /// setNodeAverages' single walk instead of paying three tree walks and three
  /// leaf loops. For each leaf in DFS (left-first) order: an empty leaf gets a
  /// forced 0/0/0 cache (byte-identical to computeLeafStats over a zero-length
  /// slice); a non-empty leaf records its atom topology, aggregates (A, G, Q)
  /// through the same misc kernel dispatch, and writes the node suffstat cache.
  /// Every value and its order match the three-call sequence exactly (each
  /// atom's aggregation reads only its own just-set topology and writes only its
  /// own cache), so this is bitwise the separate path - it is a pure per-tree
  /// pass-fusion. The separate methods stay for the component tests + Stage B.
  void buildAggregateWrite(Tree& tree, const ColumnStore& data, const double* g,
                           const double* weights) {
    members = tree.indices;
    std::size_t n = data.numObservations;
    if (trackAtomOf && atomOf.size() != n) atomOf.assign(n, 0u);

    bottomScratch.clear();
    tree.fillBottom(0, bottomScratch);

    numAtoms = 0;
    for (std::int32_t leaf : bottomScratch)
      if (tree.at(leaf).numObservations() > 0) ++numAtoms;
    ensureAtomCapacity(numAtoms);

    std::uint32_t atomId = 0;
    for (std::int32_t leaf : bottomScratch) {
      Node& node(tree.at(leaf));
      if (node.numObservations() == 0) {
        node.sumWeights = 0.0;
        node.sumWeightedResponse = 0.0;
        node.sumWeightedResponseSq = 0.0;
        continue;
      }
      atomBegin[atomId] = node.begin;
      atomEnd[atomId] = node.end;
      leafTuple[atomId] = leaf;
      if (trackAtomOf)
        for (std::size_t k = node.begin; k < node.end; ++k)
          atomOf[members[k]] = atomId;
      aggregateAtom(tree, g, weights, atomId);
      node.sumWeights = A[atomId];
      node.sumWeightedResponse = G[atomId];
      node.sumWeightedResponseSq = Q[atomId];
      ++atomId;
    }
  }

  /// Record each atom's in-block fit S(c) = mu, the constant this sweep's draw
  /// assigned the atom's leaf (design 3.3 at b=1: S is the atom's single-tree
  /// fit). `paramByNode` is the draw's per-leaf parameter keyed by node arena
  /// index, so leafTuple[c] indexes it directly. Touches no RNG and no node
  /// cache: S feeds nothing at b=1 and is inert, carried only so Stage B's S
  /// carry is in place and testable. Must run AFTER the draw.
  void setInBlockFits(const std::vector<double>& paramByNode) {
    for (std::size_t c = 0; c < numAtoms; ++c)
      S[c] = paramByNode[static_cast<std::size_t>(leafTuple[c])];
  }

  /// The atom id currently mapped to a leaf, or invalidAtom if the leaf holds
  /// no atom (it is empty). At b=1 leafTuple is a bijection onto the non-empty
  /// leaves, so a linear scan over the small atom list suffices.
  std::uint32_t atomForLeaf(std::int32_t leafId) const {
    for (std::size_t c = 0; c < numAtoms; ++c)
      if (leafTuple[c] == leafId) return static_cast<std::uint32_t>(c);
    return invalidAtom;
  }

  /// Point atom `atomId` at leaf `node` (topology) and aggregate its (A, G, Q)
  /// over that leaf's member slice through the same kernel dispatch the live
  /// computeLeafStats uses. S is a draw-time quantity and is not touched here.
  void bindAtomToLeaf(const Tree& tree, const double* g, const double* weights,
                      std::uint32_t atomId, std::int32_t node) {
    leafTuple[atomId] = node;
    atomBegin[atomId] = tree.at(node).begin;
    atomEnd[atomId] = tree.at(node).end;
    aggregateAtom(tree, g, weights, atomId);
  }

  /// Copy an atom's aggregated (A, G, Q) into its leaf's node suffstat cache -
  /// the single seam every constant-leaf consumer reads. The values are the
  /// misc-kernel result over the leaf's slice, so they are byte-for-byte the
  /// live computeLeafStats writer's.
  static void writeAtomCacheStatic(Tree& tree, std::int32_t node, double a,
                                   double g, double q) {
    Node& leaf(tree.at(node));
    leaf.sumWeights = a;
    leaf.sumWeightedResponse = g;
    leaf.sumWeightedResponseSq = q;
  }
  void writeAtomCache(Tree& tree, std::uint32_t atomId) {
    writeAtomCacheStatic(tree, leafTuple[atomId], A[atomId], G[atomId],
                         Q[atomId]);
  }

  /// Split the atom of leaf L on an accepted-or-proposed birth. L was just
  /// birthStructure'd (children attached, rule set) but NOT partitioned; this
  /// (block-fusion-stage-a.md 3(iv), pins 4 + 7):
  ///  1. reuses the tree's OWN partitionChildren primitive over the member
  ///     slice - at b=1 members alias tree.indices, so the tree's partition IS
  ///     the atom member partition; no second partitioner is forked;
  ///  2. aggregates the LEFT child then the RIGHT child through the same
  ///     kernel and writes their two node caches, bitwise-identical to birth's
  ///     live computeLeafStats(left) then (right), so the acceptance ratio and
  ///     the draw see identical inputs;
  ///  3. updates topology - the parent slot is recycled for the first non-empty
  ///     child, a fresh atom is appended for the second, atomOf follows the
  ///     moved members, leafTuple/atomBegin/atomEnd track the children.
  /// A pre-split snapshot is recorded so undoSplit can roll a rejected birth
  /// back bitwise. `g` is the current residual and `weights` the working
  /// weights (or null), matching the block-entry aggregateTree call.
  void splitAtom(Tree& tree, const ColumnStore& data, const double* g,
                 const double* weights, std::int32_t parentNode) {
    tree.partitionChildren(data, parentNode);
    std::int32_t leftNode = tree.at(parentNode).leftChild;
    std::int32_t rightNode = leftNode + 1;
    bool leftEmpty = tree.at(leftNode).numObservations() == 0;
    bool rightEmpty = tree.at(rightNode).numObservations() == 0;

    std::uint32_t parent = atomForLeaf(parentNode);
    undo_.active = true;
    undo_.hasParent = parent != invalidAtom;
    undo_.createdAtom = false;
    undo_.numAtomsBefore = numAtoms;
    undo_.movedMembers.clear();

    if (parent == invalidAtom) {
      // The split leaf held no atom, so it was empty and both children are too;
      // the birth cannot be accepted (empty-leaf veto). Force 0.0 caches to
      // match computeLeafStats over the zero-length child slices, create no
      // atoms. This path does not arise for a chain-state birth target.
      writeAtomCacheStatic(tree, leftNode, 0.0, 0.0, 0.0);
      writeAtomCacheStatic(tree, rightNode, 0.0, 0.0, 0.0);
      return;
    }

    undo_.parentAtom = parent;
    undo_.A = A[parent];
    undo_.G = G[parent];
    undo_.Q = Q[parent];
    undo_.S = S[parent];
    undo_.begin = atomBegin[parent];
    undo_.end = atomEnd[parent];
    undo_.leaf = leafTuple[parent];

    if (leftEmpty) {
      // every member went right: recycle the parent slot for the right child;
      // atomOf is unchanged (all members still belong to `parent`).
      bindAtomToLeaf(tree, g, weights, parent, rightNode);
      writeAtomCacheStatic(tree, leftNode, 0.0, 0.0, 0.0);
      writeAtomCache(tree, parent);
    } else if (rightEmpty) {
      // every member went left: recycle the parent slot for the left child.
      bindAtomToLeaf(tree, g, weights, parent, leftNode);
      writeAtomCache(tree, parent);
      writeAtomCacheStatic(tree, rightNode, 0.0, 0.0, 0.0);
    } else {
      // both children full: parent slot -> left child (pin 4: left first), a
      // fresh appended atom -> right child.
      bindAtomToLeaf(tree, g, weights, parent, leftNode);
      writeAtomCache(tree, parent);

      std::uint32_t rightAtom = static_cast<std::uint32_t>(numAtoms);
      ensureAtomCapacity(numAtoms + 1);
      ++numAtoms;
      undo_.createdAtom = true;
      undo_.rightAtom = rightAtom;
      S[rightAtom] = 0.0;
      bindAtomToLeaf(tree, g, weights, rightAtom, rightNode);
      writeAtomCache(tree, rightAtom);

      // atomOf + the moved-member roster are Stage-B bookkeeping; undoSplit's
      // atomOf restore is likewise gated, so the shipped path skips both.
      if (trackAtomOf)
        for (std::size_t k = tree.at(rightNode).begin;
             k < tree.at(rightNode).end; ++k) {
          std::size_t obs = members[k];
          atomOf[obs] = rightAtom;
          undo_.movedMembers.push_back(obs);
        }
    }
  }

  /// Roll a rejected birth's atom split back bitwise: restore the parent atom's
  /// pre-split stats + slice + leaf, drop the appended child atom, and return
  /// the moved members' atomOf to the parent. The node caches are the tree's
  /// concern (the move restores the parent cache from its saved node; the child
  /// caches die with the released pair), so this touches only the SoA. Leaves
  /// the live map byte-for-byte as it was before splitAtom.
  void undoSplit() {
    if (!undo_.active) return;
    undo_.active = false;
    if (!undo_.hasParent) return;
    std::uint32_t parent = undo_.parentAtom;
    A[parent] = undo_.A;
    G[parent] = undo_.G;
    Q[parent] = undo_.Q;
    S[parent] = undo_.S;
    atomBegin[parent] = undo_.begin;
    atomEnd[parent] = undo_.end;
    leafTuple[parent] = undo_.leaf;
    if (trackAtomOf)
      for (std::size_t obs : undo_.movedMembers) atomOf[obs] = parent;
    numAtoms = undo_.numAtomsBefore;
  }

  /// Remove atom `slot` by swapping the last live atom into its place so
  /// [0, numAtoms) stays dense. Atom ids are matched by leaf, never by index
  /// (the fuzzer's mapMatchesRebuild keys on leafTuple), so any dense filling is
  /// legal; swap-and-pop keeps it O(moved-atom members) and repoints the moved
  /// atom's atomOf entries. Used only on the death merge path.
  void removeAtomSlot(std::uint32_t slot) {
    std::uint32_t last = static_cast<std::uint32_t>(numAtoms - 1);
    if (slot != last) {
      A[slot] = A[last];
      G[slot] = G[last];
      Q[slot] = Q[last];
      S[slot] = S[last];
      leafTuple[slot] = leafTuple[last];
      atomBegin[slot] = atomBegin[last];
      atomEnd[slot] = atomEnd[last];
      if (trackAtomOf)
        for (std::size_t k = atomBegin[slot]; k < atomEnd[slot]; ++k)
          atomOf[members[k]] = slot;
    }
    --numAtoms;
  }

  /// Merge the death target's two child atoms back into one parent atom,
  /// mirroring Tree::orphanChildren at the SoA level (block-fusion-stage-a.md
  /// 1.4, pin 5). orphanChildren already forms the parent NODE cache as
  /// left + right of the (atom-written) child caches, so the cache math is left
  /// UNCHANGED and the sweep's death scoring is untouched; this patches ONLY the
  /// atom SoA topology: the two child atoms become one parent atom whose
  /// (A, G, Q) are the children's summed LEFT-then-RIGHT (pin 5, so the SoA
  /// matches the orphaned cache bitwise), leafTuple/begin/end track the parent
  /// leaf, atomOf follows the merged members, and the vacated right slot is
  /// filled by the last atom so [0, numAtoms) stays dense.
  ///
  /// The design carries A, G additively across a death (block-fusion.md 4.1); a
  /// from-scratch full-slice rescan would regroup the kernel's mod-5 sum and
  /// round differently, so the merged atom is intentionally NON-canonical until
  /// the next block-entry rebuild re-scans it - exactly as the node cache is.
  ///
  /// `leftChildNode` is the death target's pre-orphan left child (right =
  /// leftChildNode + 1); the caller passes it because orphanChildren has already
  /// detached the pair (tree.at(parentNode).leftChild == invalidNode).
  void mergeAtoms(Tree& tree, std::int32_t parentNode,
                  std::int32_t leftChildNode) {
    std::uint32_t leftAtom = atomForLeaf(leftChildNode);
    std::uint32_t rightAtom = atomForLeaf(leftChildNode + 1);
    // A chain-state death target has two non-empty leaf children, so both atoms
    // exist; tolerate a degenerate empty child defensively.
    if (leftAtom == invalidAtom && rightAtom == invalidAtom) return;
    if (leftAtom == invalidAtom) {  // only the right child held members
      leftAtom = rightAtom;
      rightAtom = invalidAtom;
    }

    std::uint32_t keep = leftAtom;
    if (rightAtom != invalidAtom) {
      A[keep] = A[keep] + A[rightAtom];  // left + right (pin 5); explicit order
      G[keep] = G[keep] + G[rightAtom];
      Q[keep] = Q[keep] + Q[rightAtom];
    }
    const Node& parent(tree.at(parentNode));
    leafTuple[keep] = parentNode;
    atomBegin[keep] = parent.begin;
    atomEnd[keep] = parent.end;
    if (trackAtomOf)
      for (std::size_t k = parent.begin; k < parent.end; ++k)
        atomOf[members[k]] = keep;
    if (rightAtom != invalidAtom) removeAtomSlot(rightAtom);
  }

  /// Rollback record for one pending change/swap subtree refresh (mirrors
  /// Tree::SubtreeSnapshot for the SoA). snapshotSubtree fills it before the
  /// rule mutates; restoreSubtree memcpys the SoA + the subtree's atomOf back so
  /// a rejected change/swap leaves the map byte-for-byte as it was. The full
  /// per-atom SoA is captured (K is #leaves, tiny) so an appended transient atom
  /// or a swap-and-pop is undone unconditionally; atomOf is captured only over
  /// the subtree's obs, keyed by obs id so the restore is independent of whether
  /// the tree's own index buffer has been restored yet.
  struct SubtreeUndo {
    bool active = false;
    std::size_t numAtoms = 0;
    std::vector<double> A, G, Q, S;
    std::vector<std::size_t> begin, end;
    std::vector<std::int32_t> leafTuple;
    std::vector<std::size_t> obsList;    // subtree obs, pre-move
    std::vector<std::uint32_t> obsAtom;  // atomOf of each subtree obs, pre-move
  } subtreeUndo_;

  /// Snapshot the SoA + the subtree's atomOf before a change/swap mutates the
  /// rule, so a rejected move restores via restoreSubtree. Call BEFORE setting
  /// the new rule / running refreshSubtree.
  void snapshotSubtree(const Tree& tree, std::int32_t subtreeRoot) {
    subtreeUndo_.active = true;
    subtreeUndo_.numAtoms = numAtoms;
    subtreeUndo_.A.assign(A.begin(), A.begin() + numAtoms);
    subtreeUndo_.G.assign(G.begin(), G.begin() + numAtoms);
    subtreeUndo_.Q.assign(Q.begin(), Q.begin() + numAtoms);
    subtreeUndo_.S.assign(S.begin(), S.begin() + numAtoms);
    subtreeUndo_.begin.assign(atomBegin.begin(), atomBegin.begin() + numAtoms);
    subtreeUndo_.end.assign(atomEnd.begin(), atomEnd.begin() + numAtoms);
    subtreeUndo_.leafTuple.assign(leafTuple.begin(),
                                  leafTuple.begin() + numAtoms);
    const Node& root(tree.at(subtreeRoot));
    subtreeUndo_.obsList.clear();
    subtreeUndo_.obsAtom.clear();
    // atomOf is Stage-B bookkeeping; when it is off there is nothing to save and
    // restoreSubtree's atomOf pass is likewise gated.
    if (trackAtomOf)
      for (std::size_t k = root.begin; k < root.end; ++k) {
        std::size_t obs = members[k];
        subtreeUndo_.obsList.push_back(obs);
        subtreeUndo_.obsAtom.push_back(atomOf[obs]);
      }
  }

  /// Restore the SoA + subtree atomOf a rejected change/swap saved. Any atom
  /// appended by a transient empty->non-empty leaf is discarded by resetting
  /// numAtoms; the [0, numAtoms) prefix is copied back bitwise. Independent of
  /// the tree's own index-buffer restore (atomOf is keyed by obs id), so the two
  /// may run in either order.
  void restoreSubtree() {
    if (!subtreeUndo_.active) return;
    subtreeUndo_.active = false;
    numAtoms = subtreeUndo_.numAtoms;
    std::copy(subtreeUndo_.A.begin(), subtreeUndo_.A.end(), A.begin());
    std::copy(subtreeUndo_.G.begin(), subtreeUndo_.G.end(), G.begin());
    std::copy(subtreeUndo_.Q.begin(), subtreeUndo_.Q.end(), Q.begin());
    std::copy(subtreeUndo_.S.begin(), subtreeUndo_.S.end(), S.begin());
    std::copy(subtreeUndo_.begin.begin(), subtreeUndo_.begin.end(),
              atomBegin.begin());
    std::copy(subtreeUndo_.end.begin(), subtreeUndo_.end.end(), atomEnd.begin());
    std::copy(subtreeUndo_.leafTuple.begin(), subtreeUndo_.leafTuple.end(),
              leafTuple.begin());
    if (trackAtomOf)
      for (std::size_t j = 0; j < subtreeUndo_.obsList.size(); ++j)
        atomOf[subtreeUndo_.obsList[j]] = subtreeUndo_.obsAtom[j];
  }

  /// Re-slice + re-aggregate one leaf's atom and write its node suffstat cache.
  /// A non-empty leaf keeps (or, for a previously-empty leaf turning non-empty
  /// in a to-be-rejected proposal, gains) exactly one atom over its current
  /// slice, aggregated through the same kernel computeLeafStats picks. An empty
  /// leaf holds no live atom: it gets a forced 0/0/0 cache (byte-identical to
  /// computeLeafStats over a zero-length slice) and any pre-existing atom is
  /// collapsed to a zero-length slice (never removed, so the change stays
  /// confined to the subtree and a later restore finds the slot intact). In a
  /// valid ACCEPTED change/swap every leaf is non-empty before and after, so
  /// this neither creates nor collapses - each subtree leaf is simply re-sliced.
  void refreshLeafAtom(Tree& tree, const double* g, const double* weights,
                       std::int32_t leafNode) {
    Node& node(tree.at(leafNode));
    std::uint32_t c = atomForLeaf(leafNode);
    if (node.numObservations() == 0) {
      node.sumWeights = 0.0;
      node.sumWeightedResponse = 0.0;
      node.sumWeightedResponseSq = 0.0;
      if (c != invalidAtom) bindAtomToLeaf(tree, g, weights, c, leafNode);
      return;
    }
    if (c == invalidAtom) {
      c = static_cast<std::uint32_t>(numAtoms);
      ensureAtomCapacity(numAtoms + 1);
      ++numAtoms;
      S[c] = 0.0;
    }
    bindAtomToLeaf(tree, g, weights, c, leafNode);
    writeAtomCache(tree, c);
    if (trackAtomOf)
      for (std::size_t k = node.begin; k < node.end; ++k) atomOf[members[k]] = c;
  }

  /// The DFS half of AtomMap::refreshSubtree: partition each internal node with
  /// the tree's OWN partitionChildren, then recurse LEFT then RIGHT, aggregating
  /// each bottom node's atom on the way down (pin 6). The call sequence is a
  /// byte-for-byte mirror of Tree::refreshSubtree, so the index-buffer
  /// permutation - and hence every leaf's member order and the kernel sum over
  /// it - is identical to the live path.
  void refreshSubtreeRec(Tree& tree, const ColumnStore& data, const double* g,
                         const double* weights, std::int32_t node) {
    if (tree.at(node).isBottom()) {
      refreshLeafAtom(tree, g, weights, node);
      return;
    }
    tree.partitionChildren(data, node);
    std::int32_t left = tree.at(node).leftChild;
    refreshSubtreeRec(tree, data, g, weights, left);
    refreshSubtreeRec(tree, data, g, weights, left + 1);
  }

  /// Repartition and re-aggregate the atoms of the subtree rooted at
  /// `subtreeRoot` after its rule changed (change/swap), mirroring
  /// Tree::refreshSubtree's DFS EXACTLY (block-fusion-stage-a.md 3(v), pin 6):
  /// the tree's own partitionChildren re-routes the members (DESIGN A: members
  /// alias tree.indices, so the tree partition IS the atom member partition; no
  /// second partitioner is forked), and each leaf's (A, G, Q) is re-aggregated
  /// through the same kernel and written into its node cache - byte-for-byte
  /// Tree::refreshSubtree's computeLeafStats, so change/swap scoring and the
  /// draw read identical inputs and the equivalence anchor holds. This REPLACES
  /// tree.refreshSubtree under the flag; running both would partition twice and
  /// scramble the member order. Snapshot with snapshotSubtree first.
  void refreshSubtree(Tree& tree, const ColumnStore& data, const double* g,
                      const double* weights, std::int32_t subtreeRoot) {
    refreshSubtreeRec(tree, data, g, weights, subtreeRoot);
  }

  // ===================================================================
  // The joint b>1 atom map (block-fusion.md 2.3, 3, 4.1-4.5).
  //
  // At blockWidth > 1 an atom is a cell of the JOINT partition of the b block
  // trees; leafOf(c, j) is its leaf in block tree j. `members` points at the
  // OWNED buffer atomMembersOwned. None of the code below runs at blockWidth
  // == 1, so the shipped path stays byte-for-byte the b=1 machinery above.
  // ===================================================================

  std::int32_t leafOf(std::size_t atomId, std::size_t j) const {
    return leafTuple[atomId * blockWidth + j];
  }
  void setLeafOf(std::size_t atomId, std::size_t j, std::int32_t leaf) {
    leafTuple[atomId * blockWidth + j] = leaf;
  }

  // Reused scratch so the b>1 build/move path does not allocate per call.
  std::vector<std::int32_t> obsTupleScratch_;  // per-obs leaf b-tuple, row-major
  std::vector<std::size_t> memberSortScratch_; // ascending-member agg buffer

  /// Aggregate atom c over its member SET: reduce g (and w) through the same
  /// misc kernel the b=1 path uses, but over the members sorted ascending into
  /// a scratch, so (A, G) depend only on WHICH observations are in the atom,
  /// not on the owned buffer's order. The birth split leaves children in the
  /// unstable two-pointer order while a from-scratch buildForBlock groups them
  /// ascending; the sort makes both land bitwise-equal (A, G), which is what
  /// the b>1 fuzzer's rebuild oracle checks. Q is carried for the harness; the
  /// interior regroup (commit iv) drops it. A is served through the atom-keyed
  /// cross-sweep cache, validated against the owned member slice.
  void aggregateAtomBlock(const double* g, const double* weights,
                          std::uint32_t atomId) {
    std::size_t begin = atomBegin[atomId];
    std::size_t length = atomEnd[atomId] - begin;
    memberSortScratch_.assign(members + begin, members + begin + length);
    std::sort(memberSortScratch_.begin(), memberSortScratch_.end());
    double a = 0.0, gMass = 0.0, q = 0.0;
    if (length > 0) {
      if (weights == nullptr)
        misc_computeIndexedSufficientStatisticsFast(g, memberSortScratch_.data(),
                                                    length, &a, &gMass, &q);
      else
        misc_computeIndexedWeightedSufficientStatisticsFast(
          g, memberSortScratch_.data(), length, weights, &a, &gMass, &q);
    }
    A[atomId] = lookupOrStoreABlock(atomId, begin, length, a);
    G[atomId] = gMass;
    Q[atomId] = q;
  }

  /// Serve/store atom c's weight mass A. Keyed by atom id, gated by a memcmp of
  /// the cached member list against the atom's OWNED slice members[begin..begin
  /// + length) -- NOT tree.indices, which at b>1 is a different (aliased) buffer
  /// the joint map does not own. On a hit the served value is byte-identical to
  /// the fresh mass (same members, same weighting, same kernel).
  double lookupOrStoreABlock(std::uint32_t atomId, std::size_t begin,
                             std::size_t length, double freshA) {
    if (aCacheBypass) return freshA;
    if (atomId >= aCacheByAtom_.size()) aCacheByAtom_.resize(atomId + 1);
    CachedLeafA& entry = aCacheByAtom_[atomId];
    if (entry.members.size() == length &&
        std::memcmp(entry.members.data(), members + begin,
                    length * sizeof(std::size_t)) == 0) {
      ++aCacheHits;
      return entry.a;
    }
    ++aCacheMisses;
    storeLeafA(entry, begin, length, freshA);
    return freshA;
  }

  /// Lay out the owned buffer + SoA from a per-observation leaf b-tuple. Groups
  /// the n observations by tuple (first-appearance atom ids, members ascending
  /// within an atom), fills atomBegin/atomEnd/leafTuple/atomOf, and aggregates
  /// every atom. This is the back half shared by the full build (buildForBlock)
  /// and the death/change patches, which differ only in how they derive the
  /// tuples: buildForBlock re-routes every tree (the independent oracle), the
  /// patches derive them incrementally from the current map plus the one
  /// tree's mutation, so a wrong patch produces different tuples and the
  /// rebuild oracle catches it. The std::map grouping is the fuzzer-scale
  /// build; the interior regroup (commits iii-iv) replaces it for run().
  void regroupByObsTuples(const std::int32_t* obsTuple, std::size_t n,
                          const double* g, const double* weights) {
    std::size_t b = blockWidth;
    if (atomOf.size() != n) atomOf.assign(n, 0u);
    atomMembersOwned.resize(n);
    members = atomMembersOwned.data();

    std::map<std::vector<std::int32_t>, std::uint32_t> ids;
    std::vector<std::vector<std::int32_t>> keys;
    numAtoms = 0;
    for (std::size_t i = 0; i < n; ++i) {
      std::vector<std::int32_t> key(obsTuple + i * b, obsTuple + i * b + b);
      auto found = ids.find(key);
      std::uint32_t id;
      if (found == ids.end()) {
        id = static_cast<std::uint32_t>(numAtoms++);
        ids.emplace(key, id);
        keys.push_back(std::move(key));
      } else {
        id = found->second;
      }
      atomOf[i] = id;
    }
    ensureAtomCapacity(numAtoms);

    std::vector<std::size_t> cursor(numAtoms, 0);
    for (std::size_t i = 0; i < n; ++i) ++cursor[atomOf[i]];
    std::size_t offset = 0;
    for (std::size_t c = 0; c < numAtoms; ++c) {
      atomBegin[c] = offset;
      offset += cursor[c];
      atomEnd[c] = offset;
      cursor[c] = atomBegin[c];
      for (std::size_t j = 0; j < b; ++j) setLeafOf(c, j, keys[c][j]);
      S[c] = 0.0;
    }
    for (std::size_t i = 0; i < n; ++i)
      members[cursor[atomOf[i]]++] = i;
    for (std::uint32_t c = 0; c < numAtoms; ++c)
      aggregateAtomBlock(g, weights, c);
  }

  /// Build the joint b>1 map from scratch: route every observation through the
  /// b block trees into its leaf b-tuple, then group (block-fusion.md 4.5, the
  /// O(bn) block-init / rebuild path). This is the fuzzer's independent oracle;
  /// it derives the tuples purely from the trees, sharing no state with the
  /// incrementally-patched map it is checked against.
  void buildForBlock(const std::vector<const Tree*>& trees,
                     const ColumnStore& data, const double* g,
                     const double* weights) {
    std::size_t n = data.numObservations;
    blockWidth = trees.size();
    obsTupleScratch_.assign(n * blockWidth, 0);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < blockWidth; ++j)
        obsTupleScratch_[i * blockWidth + j] =
          trees[j]->findBottomNodeForObservation(data, i);
    regroupByObsTuples(obsTupleScratch_.data(), n, g, weights);
  }

  /// Record atom c's in-block fit S(c) = sum over the b block trees of the
  /// tree's constant fit on the atom's leaf (block-fusion.md 2.2, the
  /// Gauss-Seidel coupling carried per atom). paramByTree[j] is block tree j's
  /// per-node parameter, keyed by arena id. Inert while the b>1 draw path is
  /// unwired (nothing reads S yet); carried so the S coupling is testable.
  void setInBlockFitsBlock(const std::vector<std::vector<double>>& paramByTree) {
    for (std::size_t c = 0; c < numAtoms; ++c) {
      double s = 0.0;
      for (std::size_t j = 0; j < blockWidth; ++j)
        s += paramByTree[j][static_cast<std::size_t>(leafOf(c, j))];
      S[c] = s;
    }
  }

  /// Block-entry static field (block-fusion.md 2.1): the frozen outside-block
  /// fit O_i = F_i - sum over the b block trees of treeFits[(t0 + j)*n + i],
  /// read off the running full-forest fit F so no O((m - b) n) rescan is paid,
  /// and g_i = w_i (y_i - O_i) (unweighted: y_i - O_i) into `g` (size n). This
  /// is the block-static per-observation field the joint map aggregates into
  /// the per-atom (A, G) masses; each g_i is one per-observation formula, so it
  /// is EXACT (no reduction regroup). `O` receives O_i when non-null.
  static void blockStaticField(const double* y, const double* weights,
                               const double* F, const double* treeFitsBase,
                               std::size_t n, std::size_t t0, std::size_t b,
                               double* g, double* O = nullptr) {
    for (std::size_t i = 0; i < n; ++i) {
      double inBlock = 0.0;
      for (std::size_t j = 0; j < b; ++j)
        inBlock += treeFitsBase[(t0 + j) * n + i];
      double o = F[i] - inBlock;
      if (O != nullptr) O[i] = o;
      g[i] = (weights == nullptr ? 1.0 : weights[i]) * (y[i] - o);
    }
  }

  /// Seed each atom's in-block fit S(c) = sum over the b block trees of the
  /// tree's CURRENT (old) fit on the atom (block-fusion.md 2.2, 3), read from
  /// the per-observation treeFits at a representative member: every member of c
  /// shares tree t_j's leaf, so its constant fit is one value over the atom.
  /// Starts the Gauss-Seidel coupling from the correct in-block fit at block
  /// entry. The block spans trees [t0, t0 + blockWidth); members alias the owned
  /// buffer after buildForBlock.
  void seedInBlockFitsFromTreeFits(const double* treeFitsBase, std::size_t n,
                                   std::size_t t0) {
    for (std::size_t c = 0; c < numAtoms; ++c) {
      std::size_t rep = members[atomBegin[c]];
      double s = 0.0;
      for (std::size_t j = 0; j < blockWidth; ++j)
        s += treeFitsBase[(t0 + j) * n + rep];
      S[c] = s;
    }
  }

  /// Block-exit fit scatter (block-fusion.md 3.5): materialize the b block
  /// trees' drawn leaf means into the per-observation treeFits over each atom's
  /// owned members - one walk of the owned buffer writes all b trees' fits - and
  /// carry the running full-forest fit F incrementally (F += new block fit - old
  /// block fit) so the next block's O needs no rescan. paramByTree[j] is block
  /// tree j's per-node parameter, arena-indexed. The scatter is an assignment,
  /// so it is EXACT; F is a running accumulation and only needs to stay a tight
  /// tolerance from a fresh full sum (the sweep-end residual is rebuilt exactly
  /// from treeFits, never from F; block-fusion.md 3.7).
  void scatterInBlockFits(double* treeFitsBase, double* F, std::size_t n,
                          std::size_t t0,
                          const std::vector<std::vector<double>>& paramByTree) {
    for (std::size_t c = 0; c < numAtoms; ++c)
      for (std::size_t m = atomBegin[c]; m < atomEnd[c]; ++m) {
        std::size_t i = members[m];
        for (std::size_t j = 0; j < blockWidth; ++j) {
          double newFit = paramByTree[j][static_cast<std::size_t>(leafOf(c, j))];
          double* slot = treeFitsBase + (t0 + j) * n + i;
          F[i] += newFit - *slot;
          *slot = newFit;
        }
      }
  }

  /// Partition the OWNED slice [begin, end) of the block by tree t_j's rule on
  /// `node`, REUSING the tree's own partitionChildren over atomMembersOwned
  /// (block-fusion.md 5.2: no second partitioner). The node is borrowed: its
  /// begin/end are pointed at the atom slice, its parent forced non-root so the
  /// dense branch takes the INDEXED kernel (misc_partitionRange assumes an
  /// identity slice, which the owned buffer is not), tree.indices is swapped to
  /// the owned buffer, and all of it is restored after. Returns the left count.
  std::size_t partitionOwnedSlice(Tree& tree, const ColumnStore& data,
                                  std::int32_t node, std::size_t begin,
                                  std::size_t end) {
    Node saved = tree.at(node);
    std::size_t* savedIndices = tree.indices;
    tree.at(node).begin = begin;
    tree.at(node).end = end;
    if (tree.at(node).parent == invalidNode) tree.at(node).parent = node;
    tree.indices = atomMembersOwned.data();
    tree.partitionChildren(data, node);
    std::size_t numOnLeft = tree.at(tree.at(node).leftChild).numObservations();
    tree.indices = savedIndices;
    tree.at(node) = saved;
    return numOnLeft;
  }

  /// Split, for a birth on block tree t_j at leaf L (already birthStructure'd),
  /// every atom whose t_j slot is L, slicing ONLY that coordinate and holding
  /// the other b-1 fixed (block-fusion.md 4.2). Each such atom's owned slice is
  /// partitioned in place by L's rule, so children stay CONTIGUOUS within the
  /// parent slice: the atom is relabelled to L_left and a fresh atom appended
  /// for L_right (or just relabelled when all members fall one side). (A, G) are
  /// re-aggregated over each child's member set. Rollback is the whole-map
  /// snapshotBlock/restoreBlock (the b>1 analogue of undoSplit).
  void splitAtomBlock(Tree& tree, const ColumnStore& data, std::size_t j,
                      std::int32_t parentLeaf, const double* g,
                      const double* weights) {
    std::int32_t leftNode = tree.at(parentLeaf).leftChild;
    std::int32_t rightNode = leftNode + 1;
    std::size_t originalNumAtoms = numAtoms;
    for (std::uint32_t c = 0; c < originalNumAtoms; ++c) {
      if (leafOf(c, j) != parentLeaf) continue;
      std::size_t begin = atomBegin[c];
      std::size_t end = atomEnd[c];
      std::size_t numOnLeft =
        partitionOwnedSlice(tree, data, parentLeaf, begin, end);
      if (numOnLeft == 0) {                      // all members went right
        setLeafOf(c, j, rightNode);
        aggregateAtomBlock(g, weights, c);
      } else if (numOnLeft == end - begin) {     // all members went left
        setLeafOf(c, j, leftNode);
        aggregateAtomBlock(g, weights, c);
      } else {
        std::size_t mid = begin + numOnLeft;
        atomEnd[c] = mid;
        setLeafOf(c, j, leftNode);
        aggregateAtomBlock(g, weights, c);

        std::uint32_t r = static_cast<std::uint32_t>(numAtoms);
        ensureAtomCapacity(numAtoms + 1);
        ++numAtoms;
        atomBegin[r] = mid;
        atomEnd[r] = end;
        for (std::size_t k = 0; k < blockWidth; ++k) setLeafOf(r, k, leafOf(c, k));
        setLeafOf(r, j, rightNode);
        S[r] = S[c];
        for (std::size_t m = mid; m < end; ++m) atomOf[members[m]] = r;
        aggregateAtomBlock(g, weights, r);
      }
    }
  }

  /// Merge, for a death on block tree t_j (leaf L reabsorbing children L_left,
  /// L_right), every atom whose t_j slot is a child back onto L, then regroup
  /// (block-fusion.md 4.1). The coordinate relabel is derived from the current
  /// map (trusting the death); the rebuild oracle re-routes t_j and catches a
  /// wrong relabel. Cross-atom merges (two atoms agreeing on the other b-1 slots
  /// once the t_j slot collapses to L) are non-adjacent in the owned buffer, so
  /// the regroup re-lays it out contiguously -- the O(n) re-canonicalization the
  /// design admits for the merge axis (4.5); the interior perf path is commit iv.
  void mergeAtomsBlock(std::size_t j, std::int32_t parentLeaf,
                       std::int32_t leftChild, std::int32_t rightChild,
                       const double* g, const double* weights) {
    std::size_t n = atomOf.size();
    std::size_t b = blockWidth;
    obsTupleScratch_.assign(n * b, 0);
    for (std::size_t i = 0; i < n; ++i) {
      std::uint32_t c = atomOf[i];
      for (std::size_t k = 0; k < b; ++k) {
        std::int32_t t = leafOf(c, k);
        if (k == j && (t == leftChild || t == rightChild)) t = parentLeaf;
        obsTupleScratch_[i * b + k] = t;
      }
    }
    regroupByObsTuples(obsTupleScratch_.data(), n, g, weights);
  }

  /// Re-slice, for a change/swap on block tree t_j's subtree rooted at P, every
  /// atom whose t_j slot is a leaf under P: re-route only those observations
  /// through t_j (the other b-1 slots are unchanged, read from the current map),
  /// then regroup (block-fusion.md 4.1, the surviving O(n_subtree) scan). A
  /// re-slice both splits an atom (its members scatter to several new t_j
  /// leaves) and merges across atoms (two old atoms landing on the same new
  /// tuple), so it goes through the contiguous regroup rather than an in-place
  /// partition. The oracle re-routes ALL trees and catches a wrong affected set
  /// or coordinate.
  void refreshSubtreeBlock(const Tree& tree, const ColumnStore& data,
                           std::size_t j, std::int32_t subtreeRoot,
                           const double* g, const double* weights) {
    std::size_t n = atomOf.size();
    std::size_t b = blockWidth;
    std::vector<std::int32_t> underLeaves;
    tree.fillBottom(subtreeRoot, underLeaves);
    std::vector<char> underSubtree(tree.nodes.size(), 0);
    for (std::int32_t leaf : underLeaves)
      underSubtree[static_cast<std::size_t>(leaf)] = 1;
    obsTupleScratch_.assign(n * b, 0);
    for (std::size_t i = 0; i < n; ++i) {
      std::uint32_t c = atomOf[i];
      for (std::size_t k = 0; k < b; ++k) obsTupleScratch_[i * b + k] = leafOf(c, k);
      std::int32_t curJ = leafOf(c, j);
      if (underSubtree[static_cast<std::size_t>(curJ)])
        obsTupleScratch_[i * b + j] =
          tree.findBottomNodeForObservation(data, i);
    }
    regroupByObsTuples(obsTupleScratch_.data(), n, g, weights);
  }

  /// Whole-map snapshot for b>1 move rollback: the SoA prefix, leafTuple (b ids
  /// per atom), atomOf, and the owned buffer. restoreBlock returns the map
  /// byte-for-byte, so a rejected birth/death/change/swap is undone bitwise.
  /// The b=1 undoSplit/restoreSubtree are targeted; at b>1 the atom count is
  /// tiny and the merge/change paths touch scattered atoms, so one snapshot is
  /// simplest and always correct (a targeted undo is a commit-iv perf lever).
  struct BlockUndo {
    bool active = false;
    std::size_t numAtoms = 0;
    std::vector<double> A, G, Q, S;
    std::vector<std::size_t> begin, end;
    std::vector<std::int32_t> leafTuple;
    std::vector<std::uint32_t> atomOf;
    std::vector<std::size_t> membersBuf;
  } blockUndo_;

  void snapshotBlock(std::size_t n) {
    blockUndo_.active = true;
    blockUndo_.numAtoms = numAtoms;
    blockUndo_.A.assign(A.begin(), A.begin() + numAtoms);
    blockUndo_.G.assign(G.begin(), G.begin() + numAtoms);
    blockUndo_.Q.assign(Q.begin(), Q.begin() + numAtoms);
    blockUndo_.S.assign(S.begin(), S.begin() + numAtoms);
    blockUndo_.begin.assign(atomBegin.begin(), atomBegin.begin() + numAtoms);
    blockUndo_.end.assign(atomEnd.begin(), atomEnd.begin() + numAtoms);
    blockUndo_.leafTuple.assign(leafTuple.begin(),
                                leafTuple.begin() + numAtoms * blockWidth);
    blockUndo_.atomOf.assign(atomOf.begin(), atomOf.begin() + n);
    blockUndo_.membersBuf.assign(atomMembersOwned.begin(),
                                 atomMembersOwned.begin() + n);
  }

  void restoreBlock(std::size_t n) {
    if (!blockUndo_.active) return;
    blockUndo_.active = false;
    numAtoms = blockUndo_.numAtoms;
    ensureAtomCapacity(numAtoms);
    std::copy(blockUndo_.A.begin(), blockUndo_.A.end(), A.begin());
    std::copy(blockUndo_.G.begin(), blockUndo_.G.end(), G.begin());
    std::copy(blockUndo_.Q.begin(), blockUndo_.Q.end(), Q.begin());
    std::copy(blockUndo_.S.begin(), blockUndo_.S.end(), S.begin());
    std::copy(blockUndo_.begin.begin(), blockUndo_.begin.end(),
              atomBegin.begin());
    std::copy(blockUndo_.end.begin(), blockUndo_.end.end(), atomEnd.begin());
    std::copy(blockUndo_.leafTuple.begin(), blockUndo_.leafTuple.end(),
              leafTuple.begin());
    std::copy(blockUndo_.atomOf.begin(), blockUndo_.atomOf.end(),
              atomOf.begin());
    atomMembersOwned.resize(n);
    std::copy(blockUndo_.membersBuf.begin(), blockUndo_.membersBuf.end(),
              atomMembersOwned.begin());
    members = atomMembersOwned.data();
  }
};

}  // namespace bartcore

#endif  // BARTCORE_ATOMS_HPP
