#ifndef BARTCORE_TREE_HPP
#define BARTCORE_TREE_HPP

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include <external/io.h>
#include <misc/stats.h>
#include <misc/linearAlgebra.h>
#include <misc/partition.h>

#include "data.hpp"

namespace bartcore {

using std::int32_t;
using std::size_t;

constexpr int32_t invalidVariable = -1;
constexpr int32_t invalidNode = -1;

/// A split rule: ordinal columns split by cut-point threshold (code >
/// splitIndex goes right), categorical columns by category subset (bit c of
/// categoryDirections set sends category c right). Which member is live is
/// determined by the column type of variableIndex.
struct Rule {
  int32_t variableIndex = invalidVariable;
  union {
    int32_t splitIndex;
    std::uint32_t categoryDirections;
  };

  Rule() : splitIndex(0) {}

  bool categoryGoesRight(xint_t code) const {
    return ((categoryDirections >> code) & 1u) != 0;
  }
  bool sendsRight(const ColumnStore& data, xint_t code) const {
    return data.types[static_cast<size_t>(variableIndex)] ==
             ColumnType::categorical
      ? categoryGoesRight(code)
      : static_cast<int32_t>(code) > splitIndex;
  }

  bool equals(const Rule& other) const {
    return variableIndex == other.variableIndex && splitIndex == other.splitIndex;
  }
};

/// One node of a value-encoded flattened tree, in pre-order (parent, left
/// subtree, right subtree). Internal nodes store the split as data values -
/// the cut point for an ordinal rule, the direction mask for a categorical
/// one - so a flattened tree can be replayed against raw predictors without
/// the store that quantized them; leaves store their parameter. The same
/// format serves saved-tree storage, external reporting, and state
/// serialization: a cut value maps back to its index exactly because cut
/// points are unique and stored as the doubles they were computed as.
struct FlatNode {
  int32_t variable = invalidVariable;  // invalidVariable for a leaf
  double value = 0.0;
};

/// Flat-arena node. Children are allocated as adjacent pairs, so
/// rightChild == leftChild + 1 always. Observation indices live in the tree's
/// external buffer; a node's subtree owns exactly [begin, end).
struct Node {
  int32_t parent = invalidNode;
  int32_t leftChild = invalidNode;
  Rule rule;
  size_t begin = 0, end = 0;
  double average = 0.0;
  double numEffectiveObservations = 0.0;

  bool isBottom() const { return leftChild == invalidNode; }
  size_t numObservations() const { return end - begin; }
};

class Tree {
public:
  std::vector<Node> nodes;
  size_t* indices = nullptr;  // external buffer, length numObservations

  void initialize(size_t* indexBuffer, size_t numObservations) {
    indices = indexBuffer;
    nodes.clear();
    freePairs.clear();
    Node root;
    root.begin = 0;
    root.end = numObservations;
    nodes.push_back(root);
    for (size_t i = 0; i < numObservations; ++i) indices[i] = i;
  }

  Node& at(int32_t i) { return nodes[static_cast<size_t>(i)]; }
  const Node& at(int32_t i) const { return nodes[static_cast<size_t>(i)]; }
  int32_t rightChildOf(int32_t i) const { return at(i).leftChild + 1; }
  bool hasSingleNode() const { return at(0).isBottom(); }

  bool childrenAreBottom(int32_t i) const {
    const Node& node(at(i));
    return !node.isBottom() && at(node.leftChild).isBottom() &&
           at(node.leftChild + 1).isBottom();
  }

  size_t depthOf(int32_t i) const {
    size_t result = 0;
    while (at(i).parent != invalidNode) {
      ++result;
      i = at(i).parent;
    }
    return result;
  }

  // Node collectors walk left-first, matching the reference engine's order.
  void fillBottom(int32_t i, std::vector<int32_t>& out) const {
    if (at(i).isBottom()) { out.push_back(i); return; }
    fillBottom(at(i).leftChild, out);
    fillBottom(at(i).leftChild + 1, out);
  }
  void fillNotBottom(int32_t i, std::vector<int32_t>& out) const {
    if (at(i).isBottom()) return;
    out.push_back(i);
    fillNotBottom(at(i).leftChild, out);
    fillNotBottom(at(i).leftChild + 1, out);
  }
  void fillNoGrand(int32_t i, std::vector<int32_t>& out) const {
    if (at(i).isBottom()) return;
    if (childrenAreBottom(i)) { out.push_back(i); return; }
    fillNoGrand(at(i).leftChild, out);
    fillNoGrand(at(i).leftChild + 1, out);
  }
  // Swappable: internal node with at least one internal child.
  void fillSwappable(int32_t i, std::vector<int32_t>& out) const {
    if (at(i).isBottom() || childrenAreBottom(i)) return;
    out.push_back(i);
    fillSwappable(at(i).leftChild, out);
    fillSwappable(at(i).leftChild + 1, out);
  }
  void fillSubtree(int32_t i, std::vector<int32_t>& out) const {
    out.push_back(i);
    if (at(i).isBottom()) return;
    fillSubtree(at(i).leftChild, out);
    fillSubtree(at(i).leftChild + 1, out);
  }

  /// Ancestor-constrained cut interval for an ordinal variable at a node;
  /// the split nearest the node on each side wins, exactly as
  /// setSplitInterval does.
  void splitInterval(const ColumnStore& data, int32_t nodeIndex,
                     int32_t variableIndex, int32_t* left, int32_t* right) const {
    *left = 0;
    *right = static_cast<int32_t>(data.numCuts[static_cast<size_t>(variableIndex)]) - 1;
    bool leftFound = false, rightFound = false;

    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode && !(leftFound && rightFound)) {
      bool isRightChild = current == at(at(current).parent).leftChild + 1;
      current = at(current).parent;
      if (at(current).rule.variableIndex == variableIndex) {
        if (isRightChild && !leftFound) {
          leftFound = true;
          *left = at(current).rule.splitIndex + 1;
        }
        if (!isRightChild && !rightFound) {
          rightFound = true;
          *right = at(current).rule.splitIndex - 1;
        }
      }
    }
  }

  /// The categories of a categorical variable that can reach a node, as a
  /// bitmask: every ancestor rule on the variable filters by the side the
  /// path descends.
  std::uint32_t reachableCategories(const ColumnStore& data, int32_t nodeIndex,
                                    int32_t variableIndex) const {
    std::uint32_t numCategories =
      data.numCuts[static_cast<size_t>(variableIndex)];
    std::uint32_t mask = numCategories >= 32
      ? 0xffffffffu
      : (1u << numCategories) - 1u;

    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode) {
      bool isRightChild = current == at(at(current).parent).leftChild + 1;
      current = at(current).parent;
      if (at(current).rule.variableIndex == variableIndex) {
        mask &= isRightChild ? at(current).rule.categoryDirections
                             : ~at(current).rule.categoryDirections;
      }
    }
    return mask;
  }

  bool variableAvailable(const ColumnStore& data, int32_t nodeIndex,
                         int32_t variableIndex) const {
    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical)
      return std::popcount(reachableCategories(data, nodeIndex,
                                               variableIndex)) >= 2;
    int32_t left, right;
    splitInterval(data, nodeIndex, variableIndex, &left, &right);
    return right >= left;
  }

  size_t numVariablesAvailable(const ColumnStore& data, int32_t nodeIndex) const {
    size_t result = 0;
    for (size_t j = 0; j < data.numPredictors; ++j)
      if (variableAvailable(data, nodeIndex, static_cast<int32_t>(j))) ++result;
    return result;
  }

  int32_t findIthAvailableVariable(const ColumnStore& data, int32_t nodeIndex,
                                   size_t i) const {
    size_t count = 0;
    for (size_t j = 0; j < data.numPredictors; ++j) {
      if (variableAvailable(data, nodeIndex, static_cast<int32_t>(j))) {
        if (count == i) return static_cast<int32_t>(j);
        ++count;
      }
    }
    return invalidVariable;
  }

  /// Leaf sufficient statistics. The root intentionally uses the non-indexed
  /// kernels like the reference engine (identical values, cheaper access).
  /// The htm entry points with a null manager run the fast serial unrolled
  /// accumulators; the plain misc_compute* versions use a slower online
  /// algorithm.
  void computeLeafStats(int32_t nodeIndex, const double* y, const double* weights) {
    Node& node(at(nodeIndex));
    bool isRoot = node.parent == invalidNode;
    if (isRoot) {
      if (weights == nullptr) {
        node.average = misc_htm_computeMean(nullptr, 0, y, node.numObservations());
        node.numEffectiveObservations = static_cast<double>(node.numObservations());
      } else {
        node.average = misc_htm_computeWeightedMean(
          nullptr, 0, y, node.numObservations(), weights,
          &node.numEffectiveObservations);
      }
    } else {
      if (weights == nullptr) {
        node.average = misc_htm_computeIndexedMean(nullptr, 0, y,
                                                   indices + node.begin,
                                                   node.numObservations());
        node.numEffectiveObservations = static_cast<double>(node.numObservations());
      } else {
        node.average = misc_htm_computeIndexedWeightedMean(
          nullptr, 0, y, indices + node.begin, node.numObservations(), weights,
          &node.numEffectiveObservations);
      }
    }
  }

  double computeVariance(int32_t nodeIndex, const double* y,
                         const double* weights) const {
    const Node& node(at(nodeIndex));
    bool isRoot = node.parent == invalidNode;
    if (isRoot) {
      return weights == nullptr
        ? misc_htm_computeVarianceForKnownMean(nullptr, 0, y,
                                               node.numObservations(), node.average)
        : misc_htm_computeWeightedVarianceForKnownMean(
            nullptr, 0, y, node.numObservations(), weights, node.average);
    }
    return weights == nullptr
      ? misc_htm_computeIndexedVarianceForKnownMean(nullptr, 0, y,
                                                    indices + node.begin,
                                                    node.numObservations(),
                                                    node.average)
      : misc_htm_computeIndexedWeightedVarianceForKnownMean(
          nullptr, 0, y, indices + node.begin, node.numObservations(), weights,
          node.average);
  }

  void setNodeAverages(const double* y, const double* weights) {
    bottomScratch.clear();
    fillBottom(0, bottomScratch);
    for (int32_t i : bottomScratch) computeLeafStats(i, y, weights);
  }

  /// Two-pointer in-place partition by category mask: bit-clear codes go
  /// left. The mask analogue of misc_partitionIndices, sans SIMD.
  static size_t partitionIndicesByMask(const xint_t* column,
                                       std::uint32_t directions,
                                       size_t* indices, size_t length) {
    size_t lo = 0, hi = length;
    // invariant: [0, lo) is left-bound, [hi, length) is right-bound
    while (true) {
      while (lo < hi && ((directions >> column[indices[lo]]) & 1u) == 0) ++lo;
      while (lo < hi && ((directions >> column[indices[hi - 1]]) & 1u) != 0)
        --hi;
      if (hi - lo < 2) break;
      size_t temp = indices[lo];
      indices[lo] = indices[hi - 1];
      indices[hi - 1] = temp;
      ++lo;
      --hi;
    }
    return lo;
  }

  /// Partition a node's observations between its children by its rule.
  void partitionChildren(const ColumnStore& data, int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    Node& left(at(node.leftChild));
    Node& right(at(node.leftChild + 1));

    size_t numOnLeft = 0;
    if (node.numObservations() > 0) {
      const xint_t* column = data.column(static_cast<size_t>(node.rule.variableIndex));
      if (data.types[static_cast<size_t>(node.rule.variableIndex)] ==
          ColumnType::categorical) {
        numOnLeft = partitionIndicesByMask(column,
                                           node.rule.categoryDirections,
                                           indices + node.begin,
                                           node.numObservations());
      } else {
        bool isRoot = node.parent == invalidNode;
        numOnLeft = isRoot
          ? misc_partitionRange(column, static_cast<misc_xint_t>(node.rule.splitIndex),
                                indices + node.begin, node.numObservations())
          : misc_partitionIndices(column, static_cast<misc_xint_t>(node.rule.splitIndex),
                                  indices + node.begin, node.numObservations());
      }
    }
    left.begin = node.begin;
    left.end = node.begin + numOnLeft;
    right.begin = left.end;
    right.end = node.end;
  }

  /// For state restore: indices already hold each node's segment as a left
  /// block followed by a right block; set the children's ranges accordingly
  /// without disturbing the stored order, whose floating-point accumulation
  /// history a bitwise-exact continuation depends on. Returns false when a
  /// segment is not actually partitioned by its node's rule.
  bool setPartitionsFromOrderedIndices(const ColumnStore& data,
                                       int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    if (node.isBottom()) return true;

    const xint_t* column =
      data.column(static_cast<size_t>(node.rule.variableIndex));
    size_t numOnLeft = 0;
    while (numOnLeft < node.numObservations() &&
           !node.rule.sendsRight(data, column[indices[node.begin + numOnLeft]]))
      ++numOnLeft;
    for (size_t k = numOnLeft; k < node.numObservations(); ++k)
      if (!node.rule.sendsRight(data, column[indices[node.begin + k]]))
        return false;

    Node& left(at(node.leftChild));
    Node& right(at(node.leftChild + 1));
    left.begin = node.begin;
    left.end = node.begin + numOnLeft;
    right.begin = left.end;
    right.end = node.end;
    return setPartitionsFromOrderedIndices(data, node.leftChild) &&
           setPartitionsFromOrderedIndices(data, node.leftChild + 1);
  }

  /// Structure-only re-route of a subtree's observations, for predictor
  /// mutation; leaf stats are left stale and refreshed by the next run().
  void repartitionSubtree(const ColumnStore& data, int32_t nodeIndex) {
    if (at(nodeIndex).isBottom()) return;
    partitionChildren(data, nodeIndex);
    repartitionSubtree(data, at(nodeIndex).leftChild);
    repartitionSubtree(data, at(nodeIndex).leftChild + 1);
  }

  /// The classic engine's validity criterion after a predictor change: no
  /// bottom node may be left without observations.
  bool bottomNodesAreOccupied() const {
    return bottomNodesAreOccupiedBelow(0);
  }

  /// Repartition a subtree after its rule changed, recomputing leaf stats.
  void refreshSubtree(const ColumnStore& data, int32_t nodeIndex, const double* y,
                      const double* weights) {
    if (at(nodeIndex).isBottom()) {
      computeLeafStats(nodeIndex, y, weights);
      return;
    }
    partitionChildren(data, nodeIndex);
    refreshSubtree(data, at(nodeIndex).leftChild, y, weights);
    refreshSubtree(data, at(nodeIndex).leftChild + 1, y, weights);
  }

  /// Split a leaf: acquire a child pair, partition, and compute child stats.
  void birth(const ColumnStore& data, int32_t nodeIndex, const Rule& rule,
             const double* y, const double* weights) {
    int32_t pair = acquirePair();
    Node& node(at(nodeIndex));  // acquirePair may reallocate; reference after
    node.rule = rule;
    node.leftChild = pair;
    at(pair).parent = nodeIndex;
    at(pair).leftChild = invalidNode;
    at(pair + 1).parent = nodeIndex;
    at(pair + 1).leftChild = invalidNode;
    partitionChildren(data, nodeIndex);
    computeLeafStats(pair, y, weights);
    computeLeafStats(pair + 1, y, weights);
  }

  void undoBirth(int32_t nodeIndex) {
    releasePair(at(nodeIndex).leftChild);
    at(nodeIndex).leftChild = invalidNode;
    at(nodeIndex).rule = Rule();
  }

  /// Merge children into their parent (death proposal). Does not free the
  /// pair so the move can be rejected; call releasePair on acceptance.
  void orphanChildren(int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    const Node& left(at(node.leftChild));
    const Node& right(at(node.leftChild + 1));
    double numEffective =
      left.numEffectiveObservations + right.numEffectiveObservations;
    node.average =
      left.average * (left.numEffectiveObservations / numEffective) +
      right.average * (right.numEffectiveObservations / numEffective);
    node.numEffectiveObservations = numEffective;
    node.leftChild = invalidNode;
  }

  int32_t acquirePair() {
    if (!freePairs.empty()) {
      int32_t result = freePairs.back();
      freePairs.pop_back();
      return result;
    }
    int32_t result = static_cast<int32_t>(nodes.size());
    nodes.emplace_back();
    nodes.emplace_back();
    return result;
  }

  void releasePair(int32_t pair) { freePairs.push_back(pair); }

  /// Snapshot/restore of a subtree for change/swap rollback: node structs
  /// plus the index-segment content the repartition scrambles.
  struct SubtreeSnapshot {
    std::vector<int32_t> nodeIds;
    std::vector<Node> nodeCopies;
    std::vector<size_t> indexSegment;
    size_t begin = 0;
  };

  void snapshotSubtree(int32_t nodeIndex, SubtreeSnapshot& snapshot) const {
    snapshot.nodeIds.clear();
    fillSubtree(nodeIndex, snapshot.nodeIds);
    snapshot.nodeCopies.clear();
    for (int32_t i : snapshot.nodeIds) snapshot.nodeCopies.push_back(at(i));
    const Node& node(at(nodeIndex));
    snapshot.begin = node.begin;
    snapshot.indexSegment.assign(indices + node.begin, indices + node.end);
  }

  void restoreSubtree(const SubtreeSnapshot& snapshot) {
    for (size_t i = 0; i < snapshot.nodeIds.size(); ++i)
      at(snapshot.nodeIds[i]) = snapshot.nodeCopies[i];
    std::memcpy(indices + snapshot.begin, snapshot.indexSegment.data(),
                snapshot.indexSegment.size() * sizeof(size_t));
  }

  /// Descend by row-major codes (test prediction).
  int32_t findBottomNodeForRow(const ColumnStore& data, const xint_t* xt) const {
    int32_t current = 0;
    while (!at(current).isBottom()) {
      const Rule& rule(at(current).rule);
      current = rule.sendsRight(data, xt[rule.variableIndex])
                  ? at(current).leftChild + 1
                  : at(current).leftChild;
    }
    return current;
  }

  /// Descend by a training observation's column-major codes, independently of
  /// the current partitions, optionally overriding one variable's code.
  int32_t findBottomNodeForObservation(const ColumnStore& data, size_t i,
                                       int32_t overrideVariable = invalidVariable,
                                       xint_t overrideCode = 0) const {
    int32_t current = 0;
    while (!at(current).isBottom()) {
      const Rule& rule(at(current).rule);
      xint_t code = rule.variableIndex == overrideVariable
        ? overrideCode
        : data.column(static_cast<size_t>(rule.variableIndex))[i];
      current = rule.sendsRight(data, code) ? at(current).leftChild + 1
                                            : at(current).leftChild;
    }
    return current;
  }

  /// Collapse any node with an unoccupied child into a leaf whose parameter
  /// is the effective-observation-weighted mean of its subtree's leaf
  /// parameters, for forced predictor updates. paramByNode is indexed by
  /// arena id; a subtree with no observations at all gets the plain mean.
  void collapseEmptyNodes(const double* weights, std::vector<double>& paramByNode) {
    collapseEmptyNodesBelow(0, weights, paramByNode);
  }

  /// Point the tree at a new observation buffer (whole-data replacement,
  /// possibly of a different size): identity indices under the root, all
  /// partitions stale until repartitionSubtree.
  void resetObservations(size_t* indexBuffer, size_t numObservations) {
    indices = indexBuffer;
    at(0).begin = 0;
    at(0).end = numObservations;
    for (size_t i = 0; i < numObservations; ++i) indices[i] = i;
  }

  /// After a from-scratch cut rebuild, remap every rule's splitIndex onto the
  /// new cut whose value is nearest the old cut, restricted to the ancestor-
  /// constrained interval; a subtree whose interval empties collapses to a
  /// leaf with the plain mean of its leaf parameters (the reference engine's
  /// mapOldCutPointsOntoNew). paramByNode is indexed by arena id.
  void mapOldCutPointsOntoNew(const ColumnStore& data,
                              const std::vector<std::vector<double>>& oldCutPoints,
                              std::vector<double>& paramByNode) {
    std::vector<int32_t> minIndices(data.numPredictors, 0);
    std::vector<int32_t> maxIndices(data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j)
      maxIndices[j] = static_cast<int32_t>(data.numCuts[j]);
    mapCutPointsBelow(0, data, oldCutPoints, paramByNode, minIndices.data(),
                      maxIndices.data());
  }

  void countVariableUses(std::uint32_t* counts) const {
    countVariableUsesBelow(0, counts);
  }

  /// Info dump in the reference engine's Node::print format: one line per
  /// node in pre-order, indented by depth, with occupancy, top/bottom flags,
  /// per-variable availability, and the rule or the leaf parameter (taken
  /// from paramByNode, indexed by arena id, on the internal scale).
  void print(const ColumnStore& data, const double* paramByNode,
             int indentation, int32_t nodeIndex = 0) const {
    const Node& node(at(nodeIndex));
    ext_printf("%*s", indentation + static_cast<int>(depthOf(nodeIndex)), "");
    ext_printf("n: %lu ", static_cast<unsigned long>(node.numObservations()));
    ext_printf("TBN: %u%u%u ", node.parent == invalidNode ? 1u : 0u,
               node.isBottom() ? 1u : 0u,
               childrenAreBottom(nodeIndex) ? 1u : 0u);
    ext_printf("Avail: ");
    for (size_t j = 0; j < data.numPredictors; ++j)
      ext_printf("%u",
                 variableAvailable(data, nodeIndex, static_cast<int32_t>(j))
                   ? 1u : 0u);

    if (!node.isBottom()) {
      ext_printf(" var: %d ", node.rule.variableIndex);
      size_t variableIndex = static_cast<size_t>(node.rule.variableIndex);
      if (data.types[variableIndex] == ColumnType::categorical) {
        ext_printf("CATRule: ");
        for (std::uint32_t i = 0; i < data.numCuts[variableIndex]; ++i)
          ext_printf(" %u", (node.rule.categoryDirections >> i) & 1);
      } else {
        ext_printf("ORDRule: (%d)=%f", node.rule.splitIndex,
                   data.cutPoints[variableIndex]
                                 [static_cast<size_t>(node.rule.splitIndex)]);
      }
    } else {
      ext_printf(" ave: %f", paramByNode[static_cast<size_t>(nodeIndex)]);
    }
    ext_printf("\n");

    if (!node.isBottom()) {
      print(data, paramByNode, indentation, node.leftChild);
      print(data, paramByNode, indentation, node.leftChild + 1);
    }
  }

  /// Flatten to pre-order value-encoded records, splits resolved against the
  /// store's cuts and leaf parameters taken from paramByNode (indexed by
  /// arena id). counts, when non-null, receives each node's current
  /// observation count in the same order.
  void flatten(const ColumnStore& data, const double* paramByNode,
               std::vector<FlatNode>& nodes,
               std::vector<std::uint32_t>* counts = nullptr) const {
    nodes.clear();
    if (counts != nullptr) counts->clear();
    flattenBelow(0, data, paramByNode, nodes, counts);
  }

  /// Rebuild structure from a flattened form into a freshly initialized
  /// (single-root) tree. Split values map back onto rules exactly: an
  /// ordinal value must equal one of its variable's cuts, a categorical mask
  /// must be a canonical-gauge assignment of the categories reachable at its
  /// node. Partitions are left stale (repartitionSubtree) and paramByNode
  /// receives leaf parameters by arena id. Returns false - possibly
  /// half-built - on a malformed input; validate on a scratch tree before
  /// building into live state.
  bool buildFromFlat(const ColumnStore& data, const FlatNode* flatNodes,
                     size_t numNodes, std::vector<double>& paramByNode) {
    paramByNode.clear();
    size_t pos = 0;
    if (!buildFromFlatBelow(0, data, flatNodes, numNodes, pos, paramByNode))
      return false;
    if (pos != numNodes) return false;
    paramByNode.resize(nodes.size(), 0.0);
    return true;
  }

  std::vector<int32_t> bottomScratch;  // reused across iterations

private:
  void countVariableUsesBelow(int32_t i, std::uint32_t* counts) const {
    if (at(i).isBottom()) return;
    ++counts[at(i).rule.variableIndex];
    countVariableUsesBelow(at(i).leftChild, counts);
    countVariableUsesBelow(at(i).leftChild + 1, counts);
  }

  bool bottomNodesAreOccupiedBelow(int32_t i) const {
    if (at(i).isBottom()) return at(i).numObservations() > 0;
    return bottomNodesAreOccupiedBelow(at(i).leftChild) &&
           bottomNodesAreOccupiedBelow(at(i).leftChild + 1);
  }

  void collapseSubtreeToLeaf(int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    if (node.isBottom()) return;
    collapseSubtreeToLeaf(node.leftChild);
    collapseSubtreeToLeaf(node.leftChild + 1);
    releasePair(node.leftChild);
    node.leftChild = invalidNode;
    node.rule = Rule();
  }

  /// minIndices are inclusive, maxIndices exclusive; both are saved and
  /// restored around the recursion, exactly as the reference walker does.
  /// Categorical rules have nothing to remap (category counts are fixed
  /// across data replacement) and pass through to their children.
  void mapCutPointsBelow(int32_t nodeIndex, const ColumnStore& data,
                         const std::vector<std::vector<double>>& oldCutPoints,
                         std::vector<double>& paramByNode,
                         int32_t* minIndices, int32_t* maxIndices) {
    if (at(nodeIndex).isBottom()) return;

    int32_t varIndex = at(nodeIndex).rule.variableIndex;

    if (data.types[static_cast<size_t>(varIndex)] == ColumnType::categorical) {
      mapCutPointsBelow(at(nodeIndex).leftChild, data, oldCutPoints,
                        paramByNode, minIndices, maxIndices);
      mapCutPointsBelow(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                        paramByNode, minIndices, maxIndices);
      return;
    }

    int32_t minIndex = minIndices[varIndex];
    int32_t maxIndex = maxIndices[varIndex];

    if (minIndex > maxIndex - 1) {
      // no split of this variable remains below the ancestors: the node is
      // fundamentally invalid, so its subtree's parameters carry little
      // information and the merge is a plain mean
      std::vector<int32_t> bottoms;
      fillBottom(nodeIndex, bottoms);
      double param = 0.0;
      for (int32_t i : bottoms) param += paramByNode[static_cast<size_t>(i)];
      paramByNode[static_cast<size_t>(nodeIndex)] =
        param / static_cast<double>(bottoms.size());
      collapseSubtreeToLeaf(nodeIndex);
      return;
    }

    double oldCut =
      oldCutPoints[static_cast<size_t>(varIndex)]
                  [static_cast<size_t>(at(nodeIndex).rule.splitIndex)];
    const double* cuts = data.cutPoints[static_cast<size_t>(varIndex)].data();

    // the first new cut below the old cut's value, then the nearer neighbor
    int32_t firstLessThan = at(nodeIndex).rule.splitIndex < maxIndex
      ? at(nodeIndex).rule.splitIndex
      : maxIndex - 1;
    while (firstLessThan < maxIndex && cuts[firstLessThan] < oldCut)
      ++firstLessThan;
    if (firstLessThan < maxIndex)
      while (firstLessThan >= minIndex && cuts[firstLessThan] >= oldCut)
        --firstLessThan;

    int32_t newIndex;
    if (firstLessThan >= maxIndex - 1) newIndex = maxIndex - 1;
    else if (firstLessThan < minIndex) newIndex = minIndex;
    else if (oldCut - cuts[firstLessThan] < cuts[firstLessThan + 1] - oldCut)
      newIndex = firstLessThan;
    else newIndex = firstLessThan + 1;  // includes an exact value match

    at(nodeIndex).rule.splitIndex = newIndex;

    maxIndices[varIndex] = newIndex;
    mapCutPointsBelow(at(nodeIndex).leftChild, data, oldCutPoints, paramByNode,
                      minIndices, maxIndices);
    maxIndices[varIndex] = maxIndex;

    minIndices[varIndex] = newIndex + 1;
    mapCutPointsBelow(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                      paramByNode, minIndices, maxIndices);
    minIndices[varIndex] = minIndex;
  }

  void collapseEmptyNodesBelow(int32_t nodeIndex, const double* weights,
                               std::vector<double>& paramByNode) {
    if (at(nodeIndex).isBottom()) return;

    if (at(at(nodeIndex).leftChild).numObservations() == 0 ||
        at(at(nodeIndex).leftChild + 1).numObservations() == 0) {
      std::vector<int32_t> bottoms;
      fillBottom(nodeIndex, bottoms);

      double weightTotal = 0.0, paramTotal = 0.0, paramSum = 0.0;
      for (int32_t i : bottoms) {
        const Node& leaf(at(i));
        double weight = weights == nullptr
          ? static_cast<double>(leaf.numObservations())
          : misc_sumIndexedVectorElements(weights, indices + leaf.begin,
                                          leaf.numObservations());
        weightTotal += weight;
        paramTotal += weight * paramByNode[static_cast<size_t>(i)];
        paramSum += paramByNode[static_cast<size_t>(i)];
      }
      paramByNode[static_cast<size_t>(nodeIndex)] = weightTotal > 0.0
        ? paramTotal / weightTotal
        : paramSum / static_cast<double>(bottoms.size());

      collapseSubtreeToLeaf(nodeIndex);
    } else {
      collapseEmptyNodesBelow(at(nodeIndex).leftChild, weights, paramByNode);
      collapseEmptyNodesBelow(at(nodeIndex).leftChild + 1, weights, paramByNode);
    }
  }

  void flattenBelow(int32_t nodeIndex, const ColumnStore& data,
                    const double* paramByNode, std::vector<FlatNode>& out,
                    std::vector<std::uint32_t>* counts) const {
    const Node& node(at(nodeIndex));
    if (counts != nullptr)
      counts->push_back(static_cast<std::uint32_t>(node.numObservations()));

    FlatNode flat;
    if (node.isBottom()) {
      flat.value = paramByNode[static_cast<size_t>(nodeIndex)];
      out.push_back(flat);
      return;
    }

    flat.variable = node.rule.variableIndex;
    flat.value =
      data.types[static_cast<size_t>(flat.variable)] == ColumnType::categorical
        ? static_cast<double>(node.rule.categoryDirections)
        : data.cutPoints[static_cast<size_t>(flat.variable)]
                        [static_cast<size_t>(node.rule.splitIndex)];
    out.push_back(flat);
    flattenBelow(node.leftChild, data, paramByNode, out, counts);
    flattenBelow(node.leftChild + 1, data, paramByNode, out, counts);
  }

  bool buildFromFlatBelow(int32_t nodeIndex, const ColumnStore& data,
                          const FlatNode* flatNodes, size_t numNodes,
                          size_t& pos, std::vector<double>& paramByNode) {
    if (pos >= numNodes) return false;
    const FlatNode& flat(flatNodes[pos++]);

    if (flat.variable == invalidVariable) {
      size_t i = static_cast<size_t>(nodeIndex);
      if (paramByNode.size() <= i) paramByNode.resize(i + 1, 0.0);
      paramByNode[i] = flat.value;
      return true;
    }

    if (flat.variable < 0 ||
        static_cast<size_t>(flat.variable) >= data.numPredictors)
      return false;

    Rule rule;
    rule.variableIndex = flat.variable;
    if (data.types[static_cast<size_t>(flat.variable)] ==
        ColumnType::categorical) {
      if (!(flat.value >= 1.0) || flat.value > 4294967295.0) return false;
      std::uint32_t directions = static_cast<std::uint32_t>(flat.value);
      if (static_cast<double>(directions) != flat.value) return false;
      std::uint32_t reachable =
        reachableCategories(data, nodeIndex, flat.variable);
      // canonical gauge: bits confined to reachable, neither side empty
      if ((directions & ~reachable) != 0 || directions == reachable)
        return false;
      rule.categoryDirections = directions;
    } else {
      const std::vector<double>& cuts(
        data.cutPoints[static_cast<size_t>(flat.variable)]);
      std::uint32_t numCuts = data.numCuts[static_cast<size_t>(flat.variable)];
      std::uint32_t k = 0;
      while (k < numCuts && cuts[k] < flat.value) ++k;
      if (k >= numCuts || cuts[k] != flat.value) return false;
      rule.splitIndex = static_cast<int32_t>(k);
    }

    int32_t pair = acquirePair();
    Node& node(at(nodeIndex));  // acquirePair may reallocate; reference after
    node.rule = rule;
    node.leftChild = pair;
    at(pair).parent = nodeIndex;
    at(pair).leftChild = invalidNode;
    at(pair + 1).parent = nodeIndex;
    at(pair + 1).leftChild = invalidNode;

    return buildFromFlatBelow(pair, data, flatNodes, numNodes, pos,
                              paramByNode) &&
           buildFromFlatBelow(pair + 1, data, flatNodes, numNodes, pos,
                              paramByNode);
  }

  std::vector<int32_t> freePairs;
};

/// Partition indices[lo, hi) of raw column-major predictors around a
/// flattened split so left-bound rows precede right-bound ones, returning
/// the boundary: ordinal rows go left when x <= value, categorical rows when
/// the mask's direction bit for their code is clear. Order within the halves
/// is not preserved (nor needed; replays only count or accumulate).
inline size_t partitionFlatIndices(const FlatNode& flat, const ColumnType* types,
                                   const double* x, size_t numRows,
                                   size_t* indices, size_t lo, size_t hi) {
  const double* column = x + static_cast<size_t>(flat.variable) * numRows;
  size_t mid = lo;
  if (types[flat.variable] == ColumnType::categorical) {
    std::uint32_t directions = static_cast<std::uint32_t>(flat.value);
    for (size_t k = lo; k < hi; ++k) {
      if (((directions >> static_cast<std::uint32_t>(column[indices[k]])) & 1u)
          == 0) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  } else {
    for (size_t k = lo; k < hi; ++k) {
      if (column[indices[k]] <= flat.value) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  }
  return mid;
}

/// Route indices[lo, hi) of a raw column-major matrix through a flattened
/// subtree, writing each node's routed count in pre-order; indices is
/// scrambled. Returns the number of flattened nodes consumed.
inline size_t countFlatObservationsBelow(const FlatNode* flatNodes,
                                         const ColumnType* types,
                                         const double* x, size_t numRows,
                                         size_t* indices, size_t lo, size_t hi,
                                         std::uint32_t* counts) {
  counts[0] = static_cast<std::uint32_t>(hi - lo);
  if (flatNodes[0].variable == invalidVariable) return 1;

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi);
  size_t numNodes = 1;
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, types, x,
                                         numRows, indices, lo, mid,
                                         counts + numNodes);
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, types, x,
                                         numRows, indices, mid, hi,
                                         counts + numNodes);
  return numNodes;
}

/// Add each routed row's leaf parameter into fits (one slot per row).
/// Returns the number of flattened nodes consumed.
inline size_t addFlatPredictionsBelow(const FlatNode* flatNodes,
                                      const ColumnType* types, const double* x,
                                      size_t numRows, size_t* indices,
                                      size_t lo, size_t hi, double* fits) {
  if (flatNodes[0].variable == invalidVariable) {
    for (size_t k = lo; k < hi; ++k) fits[indices[k]] += flatNodes[0].value;
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi);
  size_t numNodes = 1;
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, types, x, numRows,
                                      indices, lo, mid, fits);
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, types, x, numRows,
                                      indices, mid, hi, fits);
  return numNodes;
}

/// Structural well-formedness of a flattened subtree - complete pre-order,
/// variables in range, categorical masks integral and nonzero - without the
/// cut-correspondence and gauge conditions live restoration demands (saved
/// trees replay against raw values, so any split value routes). Returns the
/// number of nodes consumed, 0 when malformed.
inline size_t flatSubtreeIsWellFormed(const ColumnStore& data,
                                      const FlatNode* flatNodes,
                                      size_t numNodes, size_t pos) {
  if (pos >= numNodes) return 0;
  const FlatNode& flat(flatNodes[pos]);
  if (flat.variable == invalidVariable) return 1;
  if (flat.variable < 0 ||
      static_cast<size_t>(flat.variable) >= data.numPredictors)
    return 0;
  if (data.types[static_cast<size_t>(flat.variable)] ==
      ColumnType::categorical) {
    if (!(flat.value >= 1.0) || flat.value > 4294967295.0 ||
        static_cast<double>(static_cast<std::uint32_t>(flat.value)) !=
          flat.value)
      return 0;
  }
  size_t numOnLeft =
    flatSubtreeIsWellFormed(data, flatNodes, numNodes, pos + 1);
  if (numOnLeft == 0) return 0;
  size_t numOnRight =
    flatSubtreeIsWellFormed(data, flatNodes, numNodes, pos + 1 + numOnLeft);
  if (numOnRight == 0) return 0;
  return 1 + numOnLeft + numOnRight;
}

inline bool flatTreeIsWellFormed(const ColumnStore& data,
                                 const FlatNode* flatNodes, size_t numNodes) {
  return numNodes > 0 &&
         flatSubtreeIsWellFormed(data, flatNodes, numNodes, 0) == numNodes;
}

/// Number of records a well-formed flattened subtree occupies.
inline size_t flatSubtreeLength(const FlatNode* flatNodes) {
  if (flatNodes[0].variable == invalidVariable) return 1;
  size_t numOnLeft = flatSubtreeLength(flatNodes + 1);
  return 1 + numOnLeft + flatSubtreeLength(flatNodes + 1 + numOnLeft);
}

/// Info dump of a flattened (saved) tree in the reference engine's
/// SavedNode::print format: no occupancy or availability, ordinal splits by
/// value, leaf predictions on the internal scale.
inline void printFlatSubtree(const ColumnStore& data, const FlatNode* flatNodes,
                             int indentation, size_t depth = 0) {
  const FlatNode& flat(flatNodes[0]);
  bool isBottom = flat.variable == invalidVariable;

  size_t numOnLeft = 0;
  bool childrenAreBottom = false;
  if (!isBottom) {
    numOnLeft = flatSubtreeLength(flatNodes + 1);
    childrenAreBottom =
      flatNodes[1].variable == invalidVariable &&
      flatNodes[1 + numOnLeft].variable == invalidVariable;
  }

  ext_printf("%*s", indentation + static_cast<int>(depth), "");
  ext_printf("TBN: %u%u%u ", depth == 0 ? 1u : 0u, isBottom ? 1u : 0u,
             childrenAreBottom ? 1u : 0u);
  if (!isBottom) {
    ext_printf(" var: %d ", flat.variable);
    if (data.types[static_cast<size_t>(flat.variable)] ==
        ColumnType::categorical) {
      std::uint32_t directions = static_cast<std::uint32_t>(flat.value);
      ext_printf("CATRule: ");
      for (std::uint32_t i = 0;
           i < data.numCuts[static_cast<size_t>(flat.variable)]; ++i)
        ext_printf(" %u", (directions >> i) & 1);
    } else {
      ext_printf("ORDRule: %f", flat.value);
    }
    ext_printf("\n");
    printFlatSubtree(data, flatNodes + 1, indentation, depth + 1);
    printFlatSubtree(data, flatNodes + 1 + numOnLeft, indentation, depth + 1);
  } else {
    ext_printf(" pred: %f", flat.value);
    ext_printf("\n");
  }
}

}  // namespace bartcore

#endif  // BARTCORE_TREE_HPP
