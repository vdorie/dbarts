#ifndef BARTCORE_TREE_HPP
#define BARTCORE_TREE_HPP

#include <bit>
#include <cmath>
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

// Multi-word mask helpers for pooled categorical rules: position b lives at
// bit b % 64 of word b / 64.
inline bool maskTestBit(const std::uint64_t* words, std::uint32_t bit) {
  return ((words[bit >> 6] >> (bit & 63u)) & 1u) != 0;
}
inline void maskSetBit(std::uint64_t* words, std::uint32_t bit) {
  words[bit >> 6] |= 1ull << (bit & 63u);
}
inline size_t maskPopcount(const std::uint64_t* words, size_t numWords) {
  size_t result = 0;
  for (size_t w = 0; w < numWords; ++w)
    result += static_cast<size_t>(std::popcount(words[w]));
  return result;
}
inline bool maskIsZero(const std::uint64_t* words, size_t numWords) {
  for (size_t w = 0; w < numWords; ++w)
    if (words[w] != 0) return false;
  return true;
}
inline bool maskEquals(const std::uint64_t* a, const std::uint64_t* b,
                       size_t numWords) {
  for (size_t w = 0; w < numWords; ++w)
    if (a[w] != b[w]) return false;
  return true;
}
inline bool maskIsSubsetOf(const std::uint64_t* a, const std::uint64_t* b,
                           size_t numWords) {
  for (size_t w = 0; w < numWords; ++w)
    if ((a[w] & ~b[w]) != 0) return false;
  return true;
}
/// out = a & ~b, the left-branch reachable set below a rule.
inline void maskAndNot(const std::uint64_t* a, const std::uint64_t* b,
                       std::uint64_t* out, size_t numWords) {
  for (size_t w = 0; w < numWords; ++w) out[w] = a[w] & ~b[w];
}

/// A split rule: ordinal columns split by cut-point threshold (code >
/// splitIndex goes right), categorical columns by category subset (bit c of
/// categoryDirections set sends category c right); which is live is
/// determined by the column type of variableIndex. Both kinds share one
/// wide word - a categorical rule's full mask, an ordinal rule's split
/// index in the low half - so bit 63 (naCategory's position) holds the
/// missing direction for either kind: part of the mask arithmetic for
/// categorical rules, read explicitly for ordinal ones. Accessors keep the
/// packing endianness-independent, and equals() comparing the word makes
/// the missing direction part of rule identity.
struct Rule {
  int32_t variableIndex = invalidVariable;
  std::uint64_t bits = 0;

  static constexpr std::uint64_t missingDirectionBit = 1ull << naCategory;

  int32_t splitIndex() const {
    return static_cast<int32_t>(static_cast<std::uint32_t>(bits));
  }
  void setSplitIndex(int32_t index) {
    bits = (bits & missingDirectionBit) | static_cast<std::uint32_t>(index);
  }
  bool missingGoesRight() const { return (bits & missingDirectionBit) != 0; }
  void setMissingGoesRight(bool goesRight) {
    bits = goesRight ? (bits | missingDirectionBit)
                     : (bits & ~missingDirectionBit);
  }
  std::uint64_t categoryDirections() const { return bits; }
  void setCategoryDirections(std::uint64_t directions) { bits = directions; }

  // pooled categorical rules (a column of more than 63 categories) store an
  // offset into their tree's mask pool in place of the mask itself; the
  // word-based accessors and the missing bit are meaningless for them
  size_t maskOffset() const { return static_cast<size_t>(bits); }
  void setMaskOffset(size_t offset) {
    bits = static_cast<std::uint64_t>(offset);
  }

  bool categoryGoesRight(xint_t code) const {
    return ((bits >> code) & 1u) != 0;
  }
  bool sendsRight(const ColumnStore& data, xint_t code) const {
    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical)
      return categoryGoesRight(code);
    if (code == naCode) return missingGoesRight();
    return static_cast<int32_t>(code) > splitIndex();
  }

  bool equals(const Rule& other) const {
    return variableIndex == other.variableIndex && bits == other.bits;
  }
};

/// One node of a value-encoded flattened tree, in pre-order (parent, left
/// subtree, right subtree). Internal nodes store the split as data values -
/// the cut point for an ordinal rule, the direction mask for a categorical
/// one - so a flattened tree can be replayed against raw predictors without
/// the store that quantized them; leaves store their parameter. The same
/// format serves saved-tree storage, external reporting, and state
/// serialization: a cut value maps back to its index exactly because cut
/// points are unique and stored as the doubles they were computed as. A
/// rule on a wide categorical column (more than 53 levels, past double
/// exactness) instead stores in value its word offset into a side channel
/// of mask words holding the category bits, pre-order sequential.
struct FlatNode {
  int32_t variable = invalidVariable;  // invalidVariable for a leaf
  double value = 0.0;
  // bit 0: the rule sends missing values right; kept out of value so a
  // categorical mask stays within the 53 double-exact bits
  std::uint8_t flags = 0;
};

constexpr std::uint8_t flatMissingGoesRight = 0x1u;

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

constexpr size_t minMaskPoolCompactionSize = 256;

class Tree {
public:
  std::vector<Node> nodes;
  size_t* indices = nullptr;  // external buffer, length numObservations
  // Pooled categorical masks: rules on columns of more than 63 categories
  // store offsets into this pool. Entries are immutable once a rule
  // references them, so rule copies (snapshots, node restores) alias
  // freely; moves truncate to a mark on rejection and the chain compacts
  // between moves (docs/design/pooled-masks.md).
  std::vector<std::uint64_t> maskPool;

  void initialize(size_t* indexBuffer, size_t numObservations) {
    indices = indexBuffer;
    nodes.clear();
    freePairs.clear();
    maskPool.clear();
    maskPoolHighWater_ = minMaskPoolCompactionSize;
    Node root;
    root.begin = 0;
    root.end = numObservations;
    nodes.push_back(root);
    for (size_t i = 0; i < numObservations; ++i) indices[i] = i;
  }

  size_t maskPoolMark() const { return maskPool.size(); }
  void truncateMaskPool(size_t mark) { maskPool.resize(mark); }
  size_t allocateMask(size_t numWords) {
    size_t offset = maskPool.size();
    maskPool.resize(offset + numWords, 0);
    return offset;
  }
  const std::uint64_t* maskWordsFor(const Rule& rule) const {
    return maskPool.data() + rule.maskOffset();
  }
  std::uint64_t* mutableMaskWordsFor(size_t offset) {
    return maskPool.data() + offset;
  }

  /// Copy live pooled masks into a fresh pool once garbage from accepted
  /// changes and deaths accumulates past the high-water mark. Safe whenever
  /// no rule copies are held outside the tree (between moves); involves no
  /// generator, so draws are unaffected.
  void compactMaskPoolIfNeeded(const ColumnStore& data) {
    if (maskPool.size() <= maskPoolHighWater_) return;
    compactScratch_.clear();
    compactMaskPoolBelow(0, data, compactScratch_);
    maskPool.swap(compactScratch_);
    maskPoolHighWater_ = 4 * maskPool.size() > minMaskPoolCompactionSize
      ? 4 * maskPool.size() : minMaskPoolCompactionSize;
  }

  /// Rule dispatch for pooled columns; equality compares pool words (two
  /// offsets can hold equal masks).
  bool ruleSendsRight(const ColumnStore& data, const Rule& rule,
                      xint_t code) const {
    if (data.columnIsPooled(static_cast<size_t>(rule.variableIndex)))
      return maskTestBit(maskWordsFor(rule), code);
    return rule.sendsRight(data, code);
  }
  bool ruleMissingGoesRight(const ColumnStore& data, const Rule& rule) const {
    size_t j = static_cast<size_t>(rule.variableIndex);
    if (data.columnIsPooled(j))
      return maskTestBit(maskWordsFor(rule), data.numCuts[j]);
    return rule.missingGoesRight();
  }
  bool rulesAreEqual(const ColumnStore& data, const Rule& a,
                     const Rule& b) const {
    if (a.variableIndex != b.variableIndex) return false;
    if (a.variableIndex != invalidVariable &&
        data.columnIsPooled(static_cast<size_t>(a.variableIndex)))
      return maskEquals(
        maskWordsFor(a), maskWordsFor(b),
        maskWordsForCount(data.numCuts[static_cast<size_t>(a.variableIndex)]));
    return a.bits == b.bits;
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
          *left = at(current).rule.splitIndex() + 1;
        }
        if (!isRightChild && !rightFound) {
          rightFound = true;
          *right = at(current).rule.splitIndex() - 1;
        }
      }
    }
  }

  /// The categories of a categorical variable that can reach a node, as a
  /// bitmask: every ancestor rule on the variable filters by the side the
  /// path descends.
  std::uint64_t reachableCategories(const ColumnStore& data, int32_t nodeIndex,
                                    int32_t variableIndex) const {
    std::uint32_t numCategories =
      data.numCuts[static_cast<size_t>(variableIndex)];
    std::uint64_t mask = numCategories >= 64
      ? ~0ull
      : (1ull << numCategories) - 1ull;
    // a missing value is one more category the ancestor rules filter
    if (data.hasMissing[static_cast<size_t>(variableIndex)])
      mask |= Rule::missingDirectionBit;

    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode) {
      bool isRightChild = current == at(at(current).parent).leftChild + 1;
      current = at(current).parent;
      if (at(current).rule.variableIndex == variableIndex) {
        mask &= isRightChild ? at(current).rule.categoryDirections()
                             : ~at(current).rule.categoryDirections();
      }
    }
    return mask;
  }

  /// The pooled-column analogue of reachableCategories, writing
  /// maskWordsForCount(K) words into out: ones over the K category bits,
  /// the missing position K when the column can route one, then each
  /// ancestor rule's filter by the side the path descends.
  void reachableCategoriesWide(const ColumnStore& data, int32_t nodeIndex,
                               int32_t variableIndex,
                               std::uint64_t* out) const {
    size_t j = static_cast<size_t>(variableIndex);
    std::uint32_t numCategories = data.numCuts[j];
    size_t numWords = maskWordsForCount(numCategories);
    for (size_t w = 0; w < numWords; ++w) out[w] = ~0ull;
    for (std::uint32_t bit = numCategories;
         bit < static_cast<std::uint32_t>(64 * numWords); ++bit)
      out[bit >> 6] &= ~(1ull << (bit & 63u));
    if (data.hasMissing[j]) maskSetBit(out, numCategories);

    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode) {
      bool isRightChild = current == at(at(current).parent).leftChild + 1;
      current = at(current).parent;
      if (at(current).rule.variableIndex == variableIndex) {
        const std::uint64_t* directions = maskWordsFor(at(current).rule);
        if (isRightChild) {
          for (size_t w = 0; w < numWords; ++w) out[w] &= directions[w];
        } else {
          for (size_t w = 0; w < numWords; ++w) out[w] &= ~directions[w];
        }
      }
    }
  }

  bool variableAvailable(const ColumnStore& data, int32_t nodeIndex,
                         int32_t variableIndex) const {
    if (data.types[static_cast<size_t>(variableIndex)] ==
        ColumnType::categorical) {
      size_t j = static_cast<size_t>(variableIndex);
      if (data.columnIsPooled(j)) {
        size_t numWords = maskWordsForCount(data.numCuts[j]);
        reachableScratch_.resize(numWords);
        reachableCategoriesWide(data, nodeIndex, variableIndex,
                                reachableScratch_.data());
        return maskPopcount(reachableScratch_.data(), numWords) >= 2;
      }
      return std::popcount(reachableCategories(data, nodeIndex,
                                               variableIndex)) >= 2;
    }
    int32_t left, right;
    splitInterval(data, nodeIndex, variableIndex, &left, &right);
    return right >= left;
  }

  /// Whether at least one predictor can still split at nodeIndex. Early-exits
  /// at the first available variable, so callers that only need the boolean
  /// skip the full O(p * depth) count.
  bool hasAnyAvailableVariable(const ColumnStore& data, int32_t nodeIndex) const {
    for (size_t j = 0; j < data.numPredictors; ++j)
      if (variableAvailable(data, nodeIndex, static_cast<int32_t>(j)))
        return true;
    return false;
  }

  /// Per-variable availability at nodeIndex in a SINGLE ancestor walk
  /// (O(p + depth)) rather than one walk per predictor (O(p * depth)). Ordinal
  /// variables narrow their [left, right] cut interval and inline categoricals
  /// their reachable mask along the one walk; pooled categoricals (rare, more
  /// than 63 levels) fall back to a direct variableAvailable check. Writes a
  /// 0/1 byte per predictor into `available` and returns the number available.
  /// The bounds are nested along a root path, so the extremum over every
  /// ancestor rule equals the nearest-wins bound of splitInterval /
  /// reachableCategories: availability is bitwise identical to p separate walks.
  size_t collectAvailableVariables(const ColumnStore& data, int32_t nodeIndex,
                                   std::uint8_t* available) const {
    const size_t p = data.numPredictors;
    availLeftScratch_.resize(p);
    availRightScratch_.resize(p);
    availMaskScratch_.resize(p);
    for (size_t j = 0; j < p; ++j) {
      if (data.types[j] == ColumnType::categorical) {
        if (data.columnIsPooled(j)) continue;  // resolved after the walk
        std::uint32_t numCategories = data.numCuts[j];
        availMaskScratch_[j] =
          numCategories >= 64 ? ~0ull : (1ull << numCategories) - 1ull;
        if (data.hasMissing[j])
          availMaskScratch_[j] |= Rule::missingDirectionBit;
      } else {
        availLeftScratch_[j] = 0;
        availRightScratch_[j] = static_cast<int32_t>(data.numCuts[j]) - 1;
      }
    }

    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode) {
      bool isRightChild = current == at(at(current).parent).leftChild + 1;
      current = at(current).parent;
      size_t j = static_cast<size_t>(at(current).rule.variableIndex);
      if (data.types[j] == ColumnType::categorical) {
        if (data.columnIsPooled(j)) continue;
        availMaskScratch_[j] &= isRightChild
          ? at(current).rule.categoryDirections()
          : ~at(current).rule.categoryDirections();
      } else {
        int32_t split = at(current).rule.splitIndex();
        if (isRightChild) {
          if (split + 1 > availLeftScratch_[j]) availLeftScratch_[j] = split + 1;
        } else {
          if (split - 1 < availRightScratch_[j]) availRightScratch_[j] = split - 1;
        }
      }
    }

    size_t count = 0;
    for (size_t j = 0; j < p; ++j) {
      bool avail;
      if (data.types[j] == ColumnType::categorical)
        avail = data.columnIsPooled(j)
          ? variableAvailable(data, nodeIndex, static_cast<int32_t>(j))
          : std::popcount(availMaskScratch_[j]) >= 2;
      else
        avail = availRightScratch_[j] >= availLeftScratch_[j];
      available[j] = avail ? 1 : 0;
      count += avail ? 1 : 0;
    }
    return count;
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
                                       std::uint64_t directions,
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

  /// The pooled-mask sibling of partitionIndicesByMask: the direction bit
  /// for a code lives in the rule's pool words.
  static size_t partitionIndicesByWideMask(const xint_t* column,
                                           const std::uint64_t* directions,
                                           size_t* indices, size_t length) {
    size_t lo = 0, hi = length;
    while (true) {
      while (lo < hi && !maskTestBit(directions, column[indices[lo]])) ++lo;
      while (lo < hi && maskTestBit(directions, column[indices[hi - 1]]))
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

  /// partitionIndicesMIA over rank-bitmap storage: the sparse sibling of
  /// the dense MIA fallback (misc_partitionIndicesSparse handles NA-free
  /// sparse columns).
  static size_t partitionIndicesSparseMIA(const SparseColumnData& column,
                                          const Rule& rule, size_t* indices,
                                          size_t length) {
    int32_t splitIndex = rule.splitIndex();
    bool missingGoesRight = rule.missingGoesRight();
    auto goesRight = [&](size_t i) {
      xint_t code = column.at(i);
      return code == naCode ? missingGoesRight
                            : static_cast<int32_t>(code) > splitIndex;
    };
    size_t lo = 0, hi = length;
    while (true) {
      while (lo < hi && !goesRight(indices[lo])) ++lo;
      while (lo < hi && goesRight(indices[hi - 1])) --hi;
      if (hi - lo < 2) break;
      size_t temp = indices[lo];
      indices[lo] = indices[hi - 1];
      indices[hi - 1] = temp;
      ++lo;
      --hi;
    }
    return lo;
  }

  /// Two-pointer ordinal partition aware of the reserved missing code:
  /// codes at or below the split go left, missing codes go by the rule's
  /// direction. The scalar fallback for columns containing NAs.
  static size_t partitionIndicesMIA(const xint_t* column, const Rule& rule,
                                    size_t* indices, size_t length) {
    int32_t splitIndex = rule.splitIndex();
    bool missingGoesRight = rule.missingGoesRight();
    auto goesRight = [=](xint_t code) {
      return code == naCode ? missingGoesRight
                            : static_cast<int32_t>(code) > splitIndex;
    };
    size_t lo = 0, hi = length;
    while (true) {
      while (lo < hi && !goesRight(column[indices[lo]])) ++lo;
      while (lo < hi && goesRight(column[indices[hi - 1]])) --hi;
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
      size_t variable = static_cast<size_t>(node.rule.variableIndex);
      if (data.types[variable] == ColumnType::categorical) {
        const xint_t* column = data.column(variable);
        if (data.columnIsPooled(variable)) {
          numOnLeft = partitionIndicesByWideMask(column,
                                                 maskWordsFor(node.rule),
                                                 indices + node.begin,
                                                 node.numObservations());
        } else {
          numOnLeft = partitionIndicesByMask(column,
                                             node.rule.categoryDirections(),
                                             indices + node.begin,
                                             node.numObservations());
        }
      } else if (data.columnIsSparse(variable)) {
        // in-place partition at the root too: misc_partitionRange assumes
        // identity index content, which only the dense path maintains
        const SparseColumnData& column = data.sparseColumn(variable);
        if (data.hasMissing[variable]) {
          numOnLeft = partitionIndicesSparseMIA(column, node.rule,
                                                indices + node.begin,
                                                node.numObservations());
        } else {
          numOnLeft = misc_partitionIndicesSparse(
            column.bits.data(), column.wordRanks.data(),
            column.nzCodes.data(), column.zeroCode,
            static_cast<misc_xint_t>(node.rule.splitIndex()),
            indices + node.begin, node.numObservations());
        }
      } else {
        const xint_t* column = data.column(variable);
        bool isRoot = node.parent == invalidNode;
        if (data.hasMissing[variable]) {
          numOnLeft = partitionIndicesMIA(column, node.rule,
                                          indices + node.begin,
                                          node.numObservations());
        } else {
          numOnLeft = isRoot
            ? misc_partitionRange(column, static_cast<misc_xint_t>(node.rule.splitIndex()),
                                  indices + node.begin, node.numObservations())
            : misc_partitionIndices(column, static_cast<misc_xint_t>(node.rule.splitIndex()),
                                    indices + node.begin, node.numObservations());
        }
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

    size_t variable = static_cast<size_t>(node.rule.variableIndex);
    size_t numOnLeft = 0;
    while (numOnLeft < node.numObservations() &&
           !ruleSendsRight(data, node.rule,
                           data.codeAt(variable,
                                       indices[node.begin + numOnLeft])))
      ++numOnLeft;
    for (size_t k = numOnLeft; k < node.numObservations(); ++k)
      if (!ruleSendsRight(data, node.rule,
                          data.codeAt(variable, indices[node.begin + k])))
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
      current = ruleSendsRight(data, rule, xt[rule.variableIndex])
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
        : data.codeAt(static_cast<size_t>(rule.variableIndex), i);
      current = ruleSendsRight(data, rule, code) ? at(current).leftChild + 1
                                                 : at(current).leftChild;
    }
    return current;
  }

  /// Collapse any node with an unoccupied child into a leaf whose parameter
  /// is the effective-observation-weighted mean of its subtree's leaf
  /// parameters, for forced predictor updates. paramByNode is indexed by
  /// arena id, paramStride doubles per node, merged per coordinate; a
  /// subtree with no observations at all gets the plain mean.
  void collapseEmptyNodes(const double* weights, std::vector<double>& paramByNode,
                          size_t paramStride = 1) {
    collapseEmptyNodesBelow(0, weights, paramByNode, paramStride);
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
  /// mapOldCutPointsOntoNew). paramByNode is indexed by arena id,
  /// paramStride doubles per node, merged per coordinate.
  void mapOldCutPointsOntoNew(const ColumnStore& data,
                              const std::vector<std::vector<double>>& oldCutPoints,
                              std::vector<double>& paramByNode,
                              size_t paramStride = 1) {
    std::vector<int32_t> minIndices(data.numPredictors, 0);
    std::vector<int32_t> maxIndices(data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j)
      maxIndices[j] = static_cast<int32_t>(data.numCuts[j]);
    mapCutPointsBelow(0, data, oldCutPoints, paramByNode, minIndices.data(),
                      maxIndices.data(), paramStride);
  }

  void countVariableUses(std::uint32_t* counts) const {
    countVariableUsesBelow(0, counts);
  }

  /// Info dump in the reference engine's Node::print format: one line per
  /// node in pre-order, indented by depth, with occupancy, top/bottom flags,
  /// per-variable availability, and the rule or the leaf parameter (taken
  /// from paramByNode, indexed by arena id with paramStride doubles per
  /// node, on the internal scale); vector-parameter leaves append their
  /// slopes after the intercept.
  void print(const ColumnStore& data, const double* paramByNode,
             int indentation, int32_t nodeIndex = 0,
             size_t paramStride = 1) const {
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
        if (data.columnIsPooled(variableIndex)) {
          const std::uint64_t* directions = maskWordsFor(node.rule);
          for (std::uint32_t i = 0; i < data.numCuts[variableIndex]; ++i)
            ext_printf(" %u", maskTestBit(directions, i) ? 1u : 0u);
        } else {
          for (std::uint32_t i = 0; i < data.numCuts[variableIndex]; ++i)
            ext_printf(" %u", static_cast<unsigned int>(
                                (node.rule.categoryDirections() >> i) & 1));
        }
      } else {
        ext_printf("ORDRule: (%d)=%f", node.rule.splitIndex(),
                   data.cutPoints[variableIndex]
                                 [static_cast<size_t>(node.rule.splitIndex())]);
      }
      if (data.hasMissing[variableIndex])
        ext_printf(" NA: %c", ruleMissingGoesRight(data, node.rule) ? 'R' : 'L');
    } else {
      const double* params =
        paramByNode + static_cast<size_t>(nodeIndex) * paramStride;
      ext_printf(" ave: %f", params[0]);
      if (paramStride > 1) {
        ext_printf(" b:");
        for (size_t j = 1; j < paramStride; ++j) ext_printf(" %f", params[j]);
      }
    }
    ext_printf("\n");

    if (!node.isBottom()) {
      print(data, paramByNode, indentation, node.leftChild, paramStride);
      print(data, paramByNode, indentation, node.leftChild + 1, paramStride);
    }
  }

  /// Flatten to pre-order value-encoded records, splits resolved against the
  /// store's cuts and leaf parameters taken from paramByNode (indexed by
  /// arena id, paramStride doubles per node). counts, when non-null,
  /// receives each node's current observation count in the same order. For
  /// vector-parameter leaves a leaf's record keeps its leading coordinate
  /// (the intercept) in value and appends the remaining paramStride - 1
  /// slopes per leaf, in pre-order, to slopes. masks receives wide
  /// categorical rules' words (see FlatNode) and must be non-null when the
  /// store has any wide categorical column.
  void flatten(const ColumnStore& data, const double* paramByNode,
               std::vector<FlatNode>& nodes,
               std::vector<std::uint32_t>* counts = nullptr,
               size_t paramStride = 1,
               std::vector<double>* slopes = nullptr,
               std::vector<std::uint64_t>* masks = nullptr) const {
    nodes.clear();
    if (counts != nullptr) counts->clear();
    if (slopes != nullptr) slopes->clear();
    if (masks != nullptr) masks->clear();
    flattenBelow(0, data, paramByNode, nodes, counts, paramStride, slopes,
                 masks);
  }

  /// Rebuild structure from a flattened form into a freshly initialized
  /// (single-root) tree. Split values map back onto rules exactly: an
  /// ordinal value must equal one of its variable's cuts, a categorical mask
  /// must be a canonical-gauge assignment of the categories reachable at its
  /// node. Partitions are left stale (repartitionSubtree) and paramByNode
  /// receives leaf parameters by arena id, paramStride doubles per node -
  /// the record's value leading, then that leaf's paramStride - 1 entries of
  /// slopes (pre-order by leaf; the caller validates its length). masks is
  /// the wide categorical side channel, whose offsets must be pre-order
  /// sequential and fully consumed. Returns false - possibly half-built -
  /// on a malformed input; validate on a scratch tree before building into
  /// live state.
  bool buildFromFlat(const ColumnStore& data, const FlatNode* flatNodes,
                     size_t numNodes, std::vector<double>& paramByNode,
                     size_t paramStride = 1, const double* slopes = nullptr,
                     const std::uint64_t* masks = nullptr,
                     size_t numMaskWords = 0) {
    paramByNode.clear();
    size_t pos = 0, leafPos = 0, maskPos = 0;
    if (!buildFromFlatBelow(0, data, flatNodes, numNodes, pos, paramByNode,
                            paramStride, slopes, leafPos, masks, numMaskWords,
                            maskPos))
      return false;
    if (pos != numNodes) return false;
    if (maskPos != numMaskWords) return false;
    paramByNode.resize(nodes.size() * paramStride, 0.0);
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
                         int32_t* minIndices, int32_t* maxIndices,
                         size_t paramStride) {
    if (at(nodeIndex).isBottom()) return;

    int32_t varIndex = at(nodeIndex).rule.variableIndex;

    if (data.types[static_cast<size_t>(varIndex)] == ColumnType::categorical) {
      mapCutPointsBelow(at(nodeIndex).leftChild, data, oldCutPoints,
                        paramByNode, minIndices, maxIndices, paramStride);
      mapCutPointsBelow(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                        paramByNode, minIndices, maxIndices, paramStride);
      return;
    }

    int32_t minIndex = minIndices[varIndex];
    int32_t maxIndex = maxIndices[varIndex];

    if (minIndex > maxIndex - 1) {
      // no split of this variable remains below the ancestors: the node is
      // fundamentally invalid, so its subtree's parameters carry little
      // information and the merge is a plain mean, per coordinate
      std::vector<int32_t> bottoms;
      fillBottom(nodeIndex, bottoms);
      std::vector<double> paramSums(paramStride, 0.0);
      for (int32_t i : bottoms)
        for (size_t j = 0; j < paramStride; ++j)
          paramSums[j] += paramByNode[static_cast<size_t>(i) * paramStride + j];
      for (size_t j = 0; j < paramStride; ++j)
        paramByNode[static_cast<size_t>(nodeIndex) * paramStride + j] =
          paramSums[j] / static_cast<double>(bottoms.size());
      collapseSubtreeToLeaf(nodeIndex);
      return;
    }

    double oldCut =
      oldCutPoints[static_cast<size_t>(varIndex)]
                  [static_cast<size_t>(at(nodeIndex).rule.splitIndex())];
    const double* cuts = data.cutPoints[static_cast<size_t>(varIndex)].data();

    // the first new cut below the old cut's value, then the nearer neighbor
    int32_t firstLessThan = at(nodeIndex).rule.splitIndex() < maxIndex
      ? at(nodeIndex).rule.splitIndex()
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

    at(nodeIndex).rule.setSplitIndex(newIndex);

    maxIndices[varIndex] = newIndex;
    mapCutPointsBelow(at(nodeIndex).leftChild, data, oldCutPoints, paramByNode,
                      minIndices, maxIndices, paramStride);
    maxIndices[varIndex] = maxIndex;

    minIndices[varIndex] = newIndex + 1;
    mapCutPointsBelow(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                      paramByNode, minIndices, maxIndices, paramStride);
    minIndices[varIndex] = minIndex;
  }

  void collapseEmptyNodesBelow(int32_t nodeIndex, const double* weights,
                               std::vector<double>& paramByNode,
                               size_t paramStride) {
    if (at(nodeIndex).isBottom()) return;

    if (at(at(nodeIndex).leftChild).numObservations() == 0 ||
        at(at(nodeIndex).leftChild + 1).numObservations() == 0) {
      std::vector<int32_t> bottoms;
      fillBottom(nodeIndex, bottoms);

      double weightTotal = 0.0;
      std::vector<double> paramTotals(paramStride, 0.0);
      std::vector<double> paramSums(paramStride, 0.0);
      for (int32_t i : bottoms) {
        const Node& leaf(at(i));
        double weight = weights == nullptr
          ? static_cast<double>(leaf.numObservations())
          : misc_sumIndexedVectorElements(weights, indices + leaf.begin,
                                          leaf.numObservations());
        weightTotal += weight;
        const double* params =
          paramByNode.data() + static_cast<size_t>(i) * paramStride;
        for (size_t j = 0; j < paramStride; ++j) {
          paramTotals[j] += weight * params[j];
          paramSums[j] += params[j];
        }
      }
      double* merged =
        paramByNode.data() + static_cast<size_t>(nodeIndex) * paramStride;
      for (size_t j = 0; j < paramStride; ++j)
        merged[j] = weightTotal > 0.0
          ? paramTotals[j] / weightTotal
          : paramSums[j] / static_cast<double>(bottoms.size());

      collapseSubtreeToLeaf(nodeIndex);
    } else {
      collapseEmptyNodesBelow(at(nodeIndex).leftChild, weights, paramByNode,
                              paramStride);
      collapseEmptyNodesBelow(at(nodeIndex).leftChild + 1, weights,
                              paramByNode, paramStride);
    }
  }

  void flattenBelow(int32_t nodeIndex, const ColumnStore& data,
                    const double* paramByNode, std::vector<FlatNode>& out,
                    std::vector<std::uint32_t>* counts, size_t paramStride,
                    std::vector<double>* slopes,
                    std::vector<std::uint64_t>* masks) const {
    const Node& node(at(nodeIndex));
    if (counts != nullptr)
      counts->push_back(static_cast<std::uint32_t>(node.numObservations()));

    FlatNode flat;
    if (node.isBottom()) {
      const double* params =
        paramByNode + static_cast<size_t>(nodeIndex) * paramStride;
      flat.value = params[0];
      if (slopes != nullptr)
        for (size_t j = 1; j < paramStride; ++j) slopes->push_back(params[j]);
      out.push_back(flat);
      return;
    }

    flat.variable = node.rule.variableIndex;
    size_t j = static_cast<size_t>(flat.variable);
    if (data.types[j] != ColumnType::categorical) {
      flat.value = data.cutPoints[j]
                                 [static_cast<size_t>(node.rule.splitIndex())];
    } else if (data.columnHasWideMask(j)) {
      // side-channel encoding: value holds the rule's word offset, the
      // words carry category bits only - the missing direction stays in
      // flags for every rule kind
      std::uint32_t numCategories = data.numCuts[j];
      size_t numWords = maskWordsForCount(numCategories);
      size_t offset = masks->size();
      flat.value = static_cast<double>(offset);
      masks->resize(offset + numWords, 0);
      std::uint64_t* words = masks->data() + offset;
      if (data.columnIsPooled(j)) {
        const std::uint64_t* pooled = maskWordsFor(node.rule);
        for (size_t w = 0; w < numWords; ++w) words[w] = pooled[w];
        words[numCategories >> 6] &= ~(1ull << (numCategories & 63u));
      } else {
        words[0] = node.rule.categoryDirections() & ~Rule::missingDirectionBit;
      }
    } else {
      flat.value = static_cast<double>(node.rule.categoryDirections() &
                                       ~Rule::missingDirectionBit);
    }
    flat.flags = ruleMissingGoesRight(data, node.rule) ? flatMissingGoesRight
                                                       : 0;
    out.push_back(flat);
    flattenBelow(node.leftChild, data, paramByNode, out, counts, paramStride,
                 slopes, masks);
    flattenBelow(node.leftChild + 1, data, paramByNode, out, counts,
                 paramStride, slopes, masks);
  }

  bool buildFromFlatBelow(int32_t nodeIndex, const ColumnStore& data,
                          const FlatNode* flatNodes, size_t numNodes,
                          size_t& pos, std::vector<double>& paramByNode,
                          size_t paramStride, const double* slopes,
                          size_t& leafPos, const std::uint64_t* masks,
                          size_t numMaskWords, size_t& maskPos) {
    if (pos >= numNodes) return false;
    const FlatNode& flat(flatNodes[pos++]);

    if (flat.variable == invalidVariable) {
      if (flat.flags != 0) return false;
      size_t i = static_cast<size_t>(nodeIndex);
      if (paramByNode.size() <= (i + 1) * paramStride - 1)
        paramByNode.resize((i + 1) * paramStride, 0.0);
      paramByNode[i * paramStride] = flat.value;
      if (slopes != nullptr)
        for (size_t j = 1; j < paramStride; ++j)
          paramByNode[i * paramStride + j] =
            slopes[leafPos * (paramStride - 1) + (j - 1)];
      ++leafPos;
      return true;
    }

    if (flat.variable < 0 ||
        static_cast<size_t>(flat.variable) >= data.numPredictors ||
        (flat.flags & ~flatMissingGoesRight) != 0)
      return false;

    Rule rule;
    rule.variableIndex = flat.variable;
    size_t variable = static_cast<size_t>(flat.variable);
    if (data.types[variable] == ColumnType::categorical &&
        data.columnHasWideMask(variable)) {
      // side-channel masks: the value must be the running pre-order word
      // cursor, the words category bits only, and the assembled mask a
      // canonical-gauge assignment like the narrow path's
      std::uint32_t numCategories = data.numCuts[variable];
      size_t numWords = maskWordsForCount(numCategories);
      if (masks == nullptr || flat.value != static_cast<double>(maskPos) ||
          maskPos + numWords > numMaskWords)
        return false;
      const std::uint64_t* words = masks + maskPos;
      maskPos += numWords;
      for (std::uint32_t bit = numCategories;
           bit < static_cast<std::uint32_t>(64 * numWords); ++bit)
        if (maskTestBit(words, bit)) return false;
      if (data.columnIsPooled(variable)) {
        size_t offset = allocateMask(numWords);
        std::uint64_t* directions = mutableMaskWordsFor(offset);
        std::memcpy(directions, words, numWords * sizeof(std::uint64_t));
        if ((flat.flags & flatMissingGoesRight) != 0)
          maskSetBit(directions, numCategories);
        reachableScratch_.resize(numWords);
        reachableCategoriesWide(data, nodeIndex, flat.variable,
                                reachableScratch_.data());
        if (maskIsZero(directions, numWords) ||
            !maskIsSubsetOf(directions, reachableScratch_.data(), numWords) ||
            maskEquals(directions, reachableScratch_.data(), numWords))
          return false;
        rule.setMaskOffset(offset);
      } else {
        // the inline band: 54..63 categories in one word
        std::uint64_t directions = words[0];
        if ((flat.flags & flatMissingGoesRight) != 0)
          directions |= Rule::missingDirectionBit;
        std::uint64_t reachable =
          reachableCategories(data, nodeIndex, flat.variable);
        if (directions == 0 || (directions & ~reachable) != 0 ||
            directions == reachable)
          return false;
        rule.setCategoryDirections(directions);
      }
    } else if (data.types[variable] == ColumnType::categorical) {
      // stored masks cover the observed categories, at most 2^53 - 1 (the
      // widest double-exact value); the missing direction arrives in flags
      if (flat.value < 0.0 || flat.value > 9007199254740991.0)
        return false;
      std::uint64_t directions = static_cast<std::uint64_t>(flat.value);
      if (static_cast<double>(directions) != flat.value) return false;
      if ((flat.flags & flatMissingGoesRight) != 0)
        directions |= Rule::missingDirectionBit;
      std::uint64_t reachable =
        reachableCategories(data, nodeIndex, flat.variable);
      // canonical gauge: bits confined to reachable, neither side empty
      if (directions == 0 || (directions & ~reachable) != 0 ||
          directions == reachable)
        return false;
      rule.setCategoryDirections(directions);
    } else {
      const std::vector<double>& cuts(
        data.cutPoints[static_cast<size_t>(flat.variable)]);
      std::uint32_t numCuts = data.numCuts[static_cast<size_t>(flat.variable)];
      std::uint32_t k = 0;
      while (k < numCuts && cuts[k] < flat.value) ++k;
      if (k >= numCuts || cuts[k] != flat.value) return false;
      rule.setSplitIndex(static_cast<int32_t>(k));
      if ((flat.flags & flatMissingGoesRight) != 0) {
        if (!data.hasMissing[static_cast<size_t>(flat.variable)]) return false;
        rule.setMissingGoesRight(true);
      }
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
                              paramByNode, paramStride, slopes, leafPos,
                              masks, numMaskWords, maskPos) &&
           buildFromFlatBelow(pair + 1, data, flatNodes, numNodes, pos,
                              paramByNode, paramStride, slopes, leafPos,
                              masks, numMaskWords, maskPos);
  }

  void compactMaskPoolBelow(int32_t nodeIndex, const ColumnStore& data,
                            std::vector<std::uint64_t>& fresh) {
    Node& node(at(nodeIndex));
    if (node.isBottom()) return;
    size_t j = static_cast<size_t>(node.rule.variableIndex);
    if (data.columnIsPooled(j)) {
      size_t numWords = maskWordsForCount(data.numCuts[j]);
      const std::uint64_t* words = maskWordsFor(node.rule);
      size_t offset = fresh.size();
      fresh.insert(fresh.end(), words, words + numWords);
      node.rule.setMaskOffset(offset);
    }
    compactMaskPoolBelow(node.leftChild, data, fresh);
    compactMaskPoolBelow(node.leftChild + 1, data, fresh);
  }

  std::vector<int32_t> freePairs;
  size_t maskPoolHighWater_ = minMaskPoolCompactionSize;
  std::vector<std::uint64_t> compactScratch_;
  // wide-reachable scratch for the compute-check-discard call sites
  // (variableAvailable, buildFromFlat's gauge check); mutable because
  // availability queries are logically const
  mutable std::vector<std::uint64_t> reachableScratch_;
  // per-variable narrowing state for collectAvailableVariables' single walk
  mutable std::vector<int32_t> availLeftScratch_;
  mutable std::vector<int32_t> availRightScratch_;
  mutable std::vector<std::uint64_t> availMaskScratch_;
};

/// Partition indices[lo, hi) of raw column-major predictors around a
/// flattened split so left-bound rows precede right-bound ones, returning
/// the boundary: ordinal rows go left when x <= value, categorical rows when
/// the mask's direction bit for their code is clear. Order within the halves
/// is not preserved (nor needed; replays only count or accumulate).
/// numCategories/maskWords resolve wide categorical rules (side-channel
/// masks indexed by the record's value); both may be null for stores
/// without wide columns.
inline size_t partitionFlatIndices(const FlatNode& flat, const ColumnType* types,
                                   const double* x, size_t numRows,
                                   size_t* indices, size_t lo, size_t hi,
                                   const std::uint32_t* numCategories = nullptr,
                                   const std::uint64_t* maskWords = nullptr) {
  const double* column = x + static_cast<size_t>(flat.variable) * numRows;
  size_t mid = lo;
  bool missingGoesLeft = (flat.flags & flatMissingGoesRight) == 0;
  if (types[flat.variable] == ColumnType::categorical &&
      numCategories != nullptr &&
      numCategories[flat.variable] > maxValueEncodableCategories) {
    const std::uint64_t* directions =
      maskWords + static_cast<size_t>(flat.value);
    for (size_t k = lo; k < hi; ++k) {
      double value = column[indices[k]];
      bool goesLeft = isNA(value)
        ? missingGoesLeft
        : !maskTestBit(directions, static_cast<std::uint32_t>(value));
      if (goesLeft) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  } else if (types[flat.variable] == ColumnType::categorical) {
    std::uint64_t directions = static_cast<std::uint64_t>(flat.value);
    for (size_t k = lo; k < hi; ++k) {
      double value = column[indices[k]];
      bool goesLeft = isNA(value)
        ? missingGoesLeft
        : ((directions >> static_cast<std::uint32_t>(value)) & 1u) == 0;
      if (goesLeft) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  } else {
    for (size_t k = lo; k < hi; ++k) {
      double value = column[indices[k]];
      // a NaN comparison is false, which would silently send it right
      bool goesLeft = isNA(value) ? missingGoesLeft : value <= flat.value;
      if (goesLeft) {
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
                                         std::uint32_t* counts,
                                         const std::uint32_t* numCategories = nullptr,
                                         const std::uint64_t* maskWords = nullptr) {
  counts[0] = static_cast<std::uint32_t>(hi - lo);
  if (flatNodes[0].variable == invalidVariable) return 1;

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi,
                         numCategories, maskWords);
  size_t numNodes = 1;
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, types, x,
                                         numRows, indices, lo, mid,
                                         counts + numNodes, numCategories,
                                         maskWords);
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, types, x,
                                         numRows, indices, mid, hi,
                                         counts + numNodes, numCategories,
                                         maskWords);
  return numNodes;
}

/// Add each routed row's leaf parameter into fits (one slot per row).
/// Returns the number of flattened nodes consumed.
inline size_t addFlatPredictionsBelow(const FlatNode* flatNodes,
                                      const ColumnType* types, const double* x,
                                      size_t numRows, size_t* indices,
                                      size_t lo, size_t hi, double* fits,
                                      const std::uint32_t* numCategories = nullptr,
                                      const std::uint64_t* maskWords = nullptr) {
  if (flatNodes[0].variable == invalidVariable) {
    for (size_t k = lo; k < hi; ++k) fits[indices[k]] += flatNodes[0].value;
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi,
                         numCategories, maskWords);
  size_t numNodes = 1;
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, types, x, numRows,
                                      indices, lo, mid, fits, numCategories,
                                      maskWords);
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, types, x, numRows,
                                      indices, mid, hi, fits, numCategories,
                                      maskWords);
  return numNodes;
}

/// The linear-leaf analogue of addFlatPredictionsBelow: a routed row's fit
/// is the leaf's intercept (its record's value) plus its slopes dotted with
/// the row's designated columns, standardized on the fly by the training
/// constants (missing values enter at the standardized mean, zero). slopes
/// holds numSlopes doubles per leaf in pre-order; leafOffset counts the
/// leaves consumed before this subtree. Returns the number of flattened
/// nodes consumed.
inline size_t addFlatLinearPredictionsBelow(
    const FlatNode* flatNodes, const ColumnType* types, const double* x,
    size_t numRows, size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    size_t numSlopes, const double* slopes, size_t leafOffset = 0,
    const std::uint32_t* numCategories = nullptr,
    const std::uint64_t* maskWords = nullptr) {
  if (flatNodes[0].variable == invalidVariable) {
    const double* leafSlopes = slopes + leafOffset * numSlopes;
    for (size_t k = lo; k < hi; ++k) {
      size_t row = indices[k];
      double fit = flatNodes[0].value;
      for (size_t j = 0; j < numSlopes; ++j) {
        double value = x[columns[j] * numRows + row];
        fit += leafSlopes[j] *
               (isNA(value) ? 0.0 : (value - means[j]) / sds[j]);
      }
      fits[row] += fit;
    }
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi,
                         numCategories, maskWords);
  size_t numOnLeft = addFlatLinearPredictionsBelow(
    flatNodes + 1, types, x, numRows, indices, lo, mid, fits, columns, means,
    sds, numSlopes, slopes, leafOffset, numCategories, maskWords);
  size_t numNodes = 1 + numOnLeft;
  numNodes += addFlatLinearPredictionsBelow(
    flatNodes + numNodes, types, x, numRows, indices, mid, hi, fits, columns,
    means, sds, numSlopes, slopes, leafOffset + (numOnLeft + 1) / 2,
    numCategories, maskWords);
  return numNodes;
}

/// Function-valued leaves' saved side channel holds one variable-length
/// block per leaf in pre-order: [count, constant] when count is zero
/// (over-cap and empty leaves replay a constant), otherwise [count,
/// alpha (count), standardized covariate rows (count x numCovariates,
/// row-major)]. Walks the channel computing each leaf's block offset;
/// false when the counts are malformed or the walk does not consume the
/// channel exactly.
inline bool computeFunctionBlockOffsets(const double* blocks,
                                        size_t blocksLength, size_t numLeaves,
                                        size_t numCovariates,
                                        std::vector<size_t>& offsets) {
  offsets.assign(numLeaves, 0);
  size_t cursor = 0;
  for (size_t b = 0; b < numLeaves; ++b) {
    if (cursor >= blocksLength) return false;
    offsets[b] = cursor;
    double count = blocks[cursor];
    if (!(count >= 0.0) || count != std::floor(count) || count > 1.0e8)
      return false;
    size_t width = count == 0.0
      ? 2
      : 1 + static_cast<size_t>(count) * (1 + numCovariates);
    if (cursor + width > blocksLength) return false;
    cursor += width;
  }
  return cursor == blocksLength;
}

/// Stack scratch bound for the function-leaf replay below; the factory
/// validates designations against it.
constexpr size_t maxFunctionLeafCovariates = 8;

/// The function-leaf analogue of addFlatLinearPredictionsBelow: a routed
/// row's fit is its leaf's conditional mean c(x*)' alpha under the
/// squared-exponential kernel, the row standardized on the fly by the
/// training constants (missing values enter at the standardized mean,
/// zero) and distances scaled by the lengthscales - the identical
/// arithmetic order as the live evaluation, so replays bit-match recorded
/// test fits. Zero-count blocks add their stored constant. blockOffsets
/// indexes blocks per pre-order leaf (computeFunctionBlockOffsets);
/// leafOffset counts the leaves consumed before this subtree. Returns the
/// number of flattened nodes consumed.
inline size_t addFlatFunctionPredictionsBelow(
    const FlatNode* flatNodes, const ColumnType* types, const double* x,
    size_t numRows, size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    const double* lengthscales, size_t numCovariates, const double* blocks,
    const size_t* blockOffsets, size_t leafOffset = 0,
    const std::uint32_t* numCategories = nullptr,
    const std::uint64_t* maskWords = nullptr) {
  if (flatNodes[0].variable == invalidVariable) {
    const double* block = blocks + blockOffsets[leafOffset];
    size_t count = static_cast<size_t>(block[0]);
    if (count == 0) {
      for (size_t k = lo; k < hi; ++k) fits[indices[k]] += block[1];
      return 1;
    }
    const double* alpha = block + 1;
    const double* rows = block + 1 + count;
    double uStar[maxFunctionLeafCovariates];
    for (size_t k = lo; k < hi; ++k) {
      size_t row = indices[k];
      for (size_t j = 0; j < numCovariates; ++j) {
        double value = x[columns[j] * numRows + row];
        uStar[j] = isNA(value) ? 0.0 : (value - means[j]) / sds[j];
      }
      double fit = 0.0;
      for (size_t r = 0; r < count; ++r) {
        double distanceSq = 0.0;
        for (size_t j = 0; j < numCovariates; ++j) {
          double difference =
            (uStar[j] - rows[r * numCovariates + j]) / lengthscales[j];
          distanceSq += difference * difference;
        }
        fit += std::exp(-0.5 * distanceSq) * alpha[r];
      }
      fits[row] += fit;
    }
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], types, x, numRows, indices, lo, hi,
                         numCategories, maskWords);
  size_t numOnLeft = addFlatFunctionPredictionsBelow(
    flatNodes + 1, types, x, numRows, indices, lo, mid, fits, columns, means,
    sds, lengthscales, numCovariates, blocks, blockOffsets, leafOffset,
    numCategories, maskWords);
  size_t numNodes = 1 + numOnLeft;
  numNodes += addFlatFunctionPredictionsBelow(
    flatNodes + numNodes, types, x, numRows, indices, mid, hi, fits, columns,
    means, sds, lengthscales, numCovariates, blocks, blockOffsets,
    leafOffset + (numOnLeft + 1) / 2, numCategories, maskWords);
  return numNodes;
}

/// Structural well-formedness of a flattened subtree - complete pre-order,
/// variables in range, categorical masks integral and nonzero, wide masks'
/// offsets pre-order sequential within their channel - without the
/// cut-correspondence and gauge conditions live restoration demands (saved
/// trees replay against raw values, so any split value routes). Returns the
/// number of nodes consumed, 0 when malformed.
inline size_t flatSubtreeIsWellFormed(const ColumnStore& data,
                                      const FlatNode* flatNodes,
                                      size_t numNodes, size_t pos,
                                      const std::uint64_t* masks = nullptr,
                                      size_t numMaskWords = 0,
                                      size_t* maskCursor = nullptr) {
  if (pos >= numNodes) return 0;
  const FlatNode& flat(flatNodes[pos]);
  if (flat.variable == invalidVariable) return flat.flags == 0 ? 1 : 0;
  if (flat.variable < 0 ||
      static_cast<size_t>(flat.variable) >= data.numPredictors ||
      (flat.flags & ~flatMissingGoesRight) != 0)
    return 0;
  size_t variable = static_cast<size_t>(flat.variable);
  if (data.types[variable] == ColumnType::categorical &&
      data.columnHasWideMask(variable)) {
    size_t numWords = maskWordsForCount(data.numCuts[variable]);
    if (masks == nullptr || maskCursor == nullptr ||
        flat.value != static_cast<double>(*maskCursor) ||
        *maskCursor + numWords > numMaskWords)
      return 0;
    // the mask plus the missing direction must send something right
    if (maskIsZero(masks + *maskCursor, numWords) &&
        (flat.flags & flatMissingGoesRight) == 0)
      return 0;
    *maskCursor += numWords;
  } else if (data.types[variable] == ColumnType::categorical) {
    // the mask over observed categories plus the missing direction must
    // send something right
    if (flat.value < 0.0 || flat.value > 9007199254740991.0 ||
        static_cast<double>(static_cast<std::uint64_t>(flat.value)) !=
          flat.value ||
        (flat.value == 0.0 && (flat.flags & flatMissingGoesRight) == 0))
      return 0;
  }
  size_t numOnLeft =
    flatSubtreeIsWellFormed(data, flatNodes, numNodes, pos + 1, masks,
                            numMaskWords, maskCursor);
  if (numOnLeft == 0) return 0;
  size_t numOnRight =
    flatSubtreeIsWellFormed(data, flatNodes, numNodes, pos + 1 + numOnLeft,
                            masks, numMaskWords, maskCursor);
  if (numOnRight == 0) return 0;
  return 1 + numOnLeft + numOnRight;
}

inline bool flatTreeIsWellFormed(const ColumnStore& data,
                                 const FlatNode* flatNodes, size_t numNodes,
                                 const std::uint64_t* masks = nullptr,
                                 size_t numMaskWords = 0) {
  size_t maskCursor = 0;
  return numNodes > 0 &&
         flatSubtreeIsWellFormed(data, flatNodes, numNodes, 0, masks,
                                 numMaskWords, &maskCursor) == numNodes &&
         maskCursor == numMaskWords;
}

/// Number of records a well-formed flattened subtree occupies.
inline size_t flatSubtreeLength(const FlatNode* flatNodes) {
  if (flatNodes[0].variable == invalidVariable) return 1;
  size_t numOnLeft = flatSubtreeLength(flatNodes + 1);
  return 1 + numOnLeft + flatSubtreeLength(flatNodes + 1 + numOnLeft);
}

/// Info dump of a flattened (saved) tree in the reference engine's
/// SavedNode::print format: no occupancy or availability, ordinal splits by
/// value, leaf predictions on the internal scale. slopes, when non-null,
/// holds numSlopes per leaf in pre-order (vector-parameter leaves);
/// leafOffset counts the leaves consumed before this subtree.
inline void printFlatSubtree(const ColumnStore& data, const FlatNode* flatNodes,
                             int indentation, size_t depth = 0,
                             const double* slopes = nullptr,
                             size_t numSlopes = 0, size_t leafOffset = 0,
                             const std::uint64_t* masks = nullptr) {
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
      ext_printf("CATRule: ");
      if (data.columnHasWideMask(static_cast<size_t>(flat.variable))) {
        const std::uint64_t* directions =
          masks + static_cast<size_t>(flat.value);
        for (std::uint32_t i = 0;
             i < data.numCuts[static_cast<size_t>(flat.variable)]; ++i)
          ext_printf(" %u", maskTestBit(directions, i) ? 1u : 0u);
      } else {
        std::uint64_t directions = static_cast<std::uint64_t>(flat.value);
        for (std::uint32_t i = 0;
             i < data.numCuts[static_cast<size_t>(flat.variable)]; ++i)
          ext_printf(" %u", static_cast<unsigned int>((directions >> i) & 1));
      }
    } else {
      ext_printf("ORDRule: %f", flat.value);
    }
    if (data.hasMissing[static_cast<size_t>(flat.variable)])
      ext_printf(" NA: %c",
                 (flat.flags & flatMissingGoesRight) != 0 ? 'R' : 'L');
    ext_printf("\n");
    printFlatSubtree(data, flatNodes + 1, indentation, depth + 1, slopes,
                     numSlopes, leafOffset, masks);
    printFlatSubtree(data, flatNodes + 1 + numOnLeft, indentation, depth + 1,
                     slopes, numSlopes, leafOffset + (numOnLeft + 1) / 2,
                     masks);
  } else {
    ext_printf(" pred: %f", flat.value);
    if (slopes != nullptr && numSlopes > 0) {
      ext_printf(" b:");
      for (size_t j = 0; j < numSlopes; ++j)
        ext_printf(" %f", slopes[leafOffset * numSlopes + j]);
    }
    ext_printf("\n");
  }
}

}  // namespace bartcore

#endif  // BARTCORE_TREE_HPP
