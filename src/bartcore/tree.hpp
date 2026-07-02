#ifndef BARTCORE_TREE_HPP
#define BARTCORE_TREE_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include <misc/stats.h>
#include <misc/linearAlgebra.h>
#include <misc/partition.h>

#include "data.hpp"

namespace bartcore {

using std::int32_t;
using std::size_t;

constexpr int32_t invalidVariable = -1;
constexpr int32_t invalidNode = -1;

/// Ordinal split; the rule vocabulary grows in later phases.
struct Rule {
  int32_t variableIndex = invalidVariable;
  int32_t splitIndex = 0;

  bool equals(const Rule& other) const {
    return variableIndex == other.variableIndex && splitIndex == other.splitIndex;
  }
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

  /// Ancestor-constrained cut interval for a variable at a node; the split
  /// nearest the node on each side wins, exactly as setSplitInterval does.
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

  bool variableAvailable(const ColumnStore& data, int32_t nodeIndex,
                         int32_t variableIndex) const {
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

  /// Partition a node's observations between its children by its rule.
  void partitionChildren(const ColumnStore& data, int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    Node& left(at(node.leftChild));
    Node& right(at(node.leftChild + 1));

    size_t numOnLeft = 0;
    if (node.numObservations() > 0) {
      const xint_t* column = data.column(static_cast<size_t>(node.rule.variableIndex));
      bool isRoot = node.parent == invalidNode;
      numOnLeft = isRoot
        ? misc_partitionRange(column, static_cast<misc_xint_t>(node.rule.splitIndex),
                              indices + node.begin, node.numObservations())
        : misc_partitionIndices(column, static_cast<misc_xint_t>(node.rule.splitIndex),
                                indices + node.begin, node.numObservations());
    }
    left.begin = node.begin;
    left.end = node.begin + numOnLeft;
    right.begin = left.end;
    right.end = node.end;
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
  int32_t findBottomNodeForRow(const xint_t* xt) const {
    int32_t current = 0;
    while (!at(current).isBottom()) {
      const Rule& rule(at(current).rule);
      current = xt[rule.variableIndex] > rule.splitIndex
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
      current = code > rule.splitIndex ? at(current).leftChild + 1
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
  void mapCutPointsBelow(int32_t nodeIndex, const ColumnStore& data,
                         const std::vector<std::vector<double>>& oldCutPoints,
                         std::vector<double>& paramByNode,
                         int32_t* minIndices, int32_t* maxIndices) {
    if (at(nodeIndex).isBottom()) return;

    int32_t varIndex = at(nodeIndex).rule.variableIndex;
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

  std::vector<int32_t> freePairs;
};

}  // namespace bartcore

#endif  // BARTCORE_TREE_HPP
