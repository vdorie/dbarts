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

  /// Sentinel for "no atom" (an empty leaf holds none).
  static constexpr std::uint32_t invalidAtom = 0xFFFFFFFFu;

  /// Reused fillBottom scratch so aggregation does not allocate per call.
  std::vector<std::int32_t> bottomScratch;

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

      for (std::size_t k = tree.at(rightNode).begin; k < tree.at(rightNode).end;
           ++k) {
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
    for (std::size_t obs : undo_.movedMembers) atomOf[obs] = parent;
    numAtoms = undo_.numAtomsBefore;
  }
};

}  // namespace bartcore

#endif  // BARTCORE_ATOMS_HPP
