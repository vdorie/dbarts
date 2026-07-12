#ifndef BARTCORE_ATOMS_HPP
#define BARTCORE_ATOMS_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include <misc/stats.h>

#include "data.hpp"
#include "tree.hpp"

/// Block-fusion Stage A build knob (block-fusion-stage-a.md 5.1). A compile-time
/// switch, DEFAULT OFF: the shipped constant-leaf sweep sources each sweep's
/// per-leaf sufficient statistics from the live setNodeAverages/computeLeafStats
/// writer, so the released engine is byte-for-byte today's. A dev/test build
/// defines BARTCORE_BLOCK_FUSION to non-zero (e.g. -DBARTCORE_BLOCK_FUSION=1) to
/// route that per-sweep suffstat SOURCE through the b=1 atom path
/// (aggregateTree + writeNodeCaches) for the equivalence gate; because the atom
/// aggregation is bitwise the kernel and lands the same values in the same node
/// cache, the draws are unchanged. Both writers stay compiled regardless of the
/// switch, so tests/cpp can drive them side by side. Commit (vii) flips the
/// default ON for the constant-leaf steady-state sweep; this commit keeps the
/// shipped default OFF (clean abort: define it to 0 or leave it undefined).
#ifndef BARTCORE_BLOCK_FUSION
#define BARTCORE_BLOCK_FUSION 0
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

  // --- member buffer (DESIGN A: aliases the block tree's index buffer) ---
  std::size_t* members = nullptr;

  /// Live atom count == non-empty-leaf count at b=1.
  std::size_t numAtoms = 0;

  /// Reused fillBottom scratch so aggregation does not allocate per call.
  std::vector<std::int32_t> bottomScratch;

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
    numAtoms = 0;
  }

  /// Grow every per-atom vector to hold k atoms.
  void ensureAtomCapacity(std::size_t k) {
    atomBegin.resize(k);
    atomEnd.resize(k);
    leafTuple.resize(k);
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
    A[atomId] = a;
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
};

}  // namespace bartcore

#endif  // BARTCORE_ATOMS_HPP
