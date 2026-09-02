#ifndef BARTCORE_TREE_HPP
#define BARTCORE_TREE_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
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
inline void maskClearBit(std::uint64_t* words, std::uint32_t bit) {
  words[bit >> 6] &= ~(1ull << (bit & 63u));
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

/// Per-forest interaction constraint: a cap on the number of DISTINCT split
/// variables along any
/// root-to-leaf path (maxOrder) and/or a symmetric forbidden co-occurrence
/// adjacency (variable pairs that may not share a path). Borrowed by the trees
/// of a single forest through Tree::setInteractionConstraint; a null pointer -
/// or an inactive constraint - leaves the availability path byte-for-byte
/// unchanged. The predicate collapses to one admissibility test p(j | A) over
/// the ancestor split-variable set A: max-order is |distinct(A u {j})| <= K,
/// co-occurrence is "no forbidden (j, a) for a in A".
struct InteractionConstraint {
  size_t numPredictors = 0;
  size_t maxOrder = 0;   // 0 leaves the order uncapped
  size_t numWords = 0;   // 64-bit words per p-bitset row: (numPredictors + 63) / 64
  // numPredictors rows of numWords words; bit a of row j is set iff the pair
  // (j, a) may not co-occur. Empty leaves co-occurrence unrestricted.
  std::vector<std::uint64_t> forbidden;

  /// Materialize from a max-order K and a flat list of forbidden pairs (2 *
  /// numPairs column indices, each pair symmetric). Borrowed inputs; the
  /// forbidden adjacency is copied here.
  void build(size_t numPredictors_, size_t maxOrder_,
             const size_t* forbiddenPairs, size_t numPairs) {
    numPredictors = numPredictors_;
    maxOrder = maxOrder_;
    numWords = (numPredictors + 63) / 64;
    forbidden.clear();
    if (numPairs > 0) {
      forbidden.assign(numPredictors * numWords, 0);
      for (size_t k = 0; k < numPairs; ++k) {
        std::uint32_t a = static_cast<std::uint32_t>(forbiddenPairs[2 * k]);
        std::uint32_t b = static_cast<std::uint32_t>(forbiddenPairs[2 * k + 1]);
        maskSetBit(forbidden.data() + a * numWords, b);
        maskSetBit(forbidden.data() + b * numWords, a);
      }
    }
  }

  bool hasOrderCap() const {
    return maxOrder > 0 && maxOrder < numPredictors;
  }
  bool hasForbidden() const { return !forbidden.empty(); }
  bool active() const { return hasOrderCap() || hasForbidden(); }
  const std::uint64_t* forbiddenRow(size_t j) const {
    return forbidden.data() + j * numWords;
  }
};

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
    if (data.splitsBySubset(static_cast<size_t>(variableIndex)))
      return categoryGoesRight(code);
    if (code == naCode) return missingGoesRight();
    return static_cast<int32_t>(code) > splitIndex();
  }

  bool equals(const Rule& other) const {
    return variableIndex == other.variableIndex && bits == other.bits;
  }
};

/// The split kind a FlatNode's payload carries; kept in flags so the replay
/// family routes raw predictors without the store that typed the columns. It
/// is the MECHANIC axis, so a column kind that splits by threshold tags
/// ordinal whatever its semantics; expectedFlatKind below derives it.
enum class FlatKind : std::uint8_t {
  leaf = 0,            // payload: the leaf parameter (an intercept)
  ordinal,             // payload: the cut point, as the double it was
  categoricalInline,   // payload: the direction mask, category bits only
  categoricalPooled    // payload: a word offset into the mask side channel
};

/// One node of a flattened tree, in pre-order (parent, left subtree, right
/// subtree). Internal nodes carry their split in a tagged payload - an
/// ordinal cut point or a categorical direction mask - so a flattened tree
/// replays against raw predictors without the store that quantized them;
/// leaves carry their parameter. The same format serves saved-tree storage,
/// external reporting, and state serialization: a cut value maps back to its
/// index exactly because cut points are unique and stored as the doubles
/// they were computed as. An inline mask (up to 63 categories) rides in the
/// payload word directly, category bits only; a pooled column's wider mask
/// keeps numMaskWords words at maskOffset in a per-tree side channel,
/// pre-order sequential. The missing direction lives in flags for either
/// kind.
struct FlatNode {
  int32_t variable = invalidVariable;  // invalidVariable for a leaf
  std::uint32_t numMaskWords = 0;      // categoricalPooled: words at maskOffset
  union {
    double value = 0.0;        // leaf parameter or ordinal cut point
    std::uint64_t mask;        // categoricalInline direction bits
    std::uint64_t maskOffset;  // categoricalPooled: word offset in the channel
  };
  std::uint8_t flags = 0;      // bit 0 missing-right; bits 1-2 hold FlatKind
};

constexpr std::uint8_t flatMissingGoesRight = 0x1u;
constexpr std::uint8_t flatKindMask = 0x6u;  // bits 1-2
constexpr std::uint8_t flatKindShift = 1u;

inline FlatKind flatKindOf(const FlatNode& node) {
  return static_cast<FlatKind>((node.flags & flatKindMask) >> flatKindShift);
}
inline void setFlatKind(FlatNode& node, FlatKind kind) {
  node.flags = static_cast<std::uint8_t>(
    (node.flags & ~flatKindMask) |
    (static_cast<std::uint8_t>(kind) << flatKindShift));
}

/// The tag an internal node splitting on \p variable must carry. The tag is
/// what lets the replay family route without a store, so it is written here
/// and checked against here rather than re-derived per ladder: one definition
/// for the flatten, the rebuild and the well-formedness check, each of which
/// then holds only its own payload arm.
inline FlatKind expectedFlatKind(const ColumnStore& data, size_t variable) {
  if (!data.splitsBySubset(variable)) return FlatKind::ordinal;
  return data.columnIsPooled(variable) ? FlatKind::categoricalPooled
                                       : FlatKind::categoricalInline;
}

/// Flat-arena node. Children are allocated as adjacent pairs, so
/// rightChild == leftChild + 1 always. Observation indices live in the tree's
/// external buffer; a node's subtree owns exactly [begin, end).
struct Node {
  int32_t parent = invalidNode;
  int32_t leftChild = invalidNode;
  Rule rule;
  size_t begin = 0, end = 0;
  // Constant-leaf sufficient statistic (sum w, sum wz); sumWeights is the
  // effective count (n unweighted). Raw sum wz^2 cancels in every MH ratio.
  double sumWeights = 0.0;
  double sumWeightedResponse = 0.0;

  bool isBottom() const { return leftChild == invalidNode; }
  size_t numObservations() const { return end - begin; }
};

constexpr size_t minMaskPoolCompactionSize = 256;

/// How a collapse merges the leaf parameters of the subtree it replaces.
/// A policy maps a parameter into the space the weighted mean is taken in and
/// back; the merge itself (weights, the no-weight fallback) is common.
///
/// ArithmeticMerge is the additive rule every mean forest uses and the default
/// argument everywhere, so both merge sites compile to the code they carried
/// before the policy existed - a codegen or bitwise difference on a
/// homoscedastic path is a defect, not a tolerance.
struct ArithmeticMerge {
  static double toMergeSpace(double param) { return param; }
  static double fromMergeSpace(double merged) { return merged; }
};

/// The multiplicative analogue, for the variance forest's scale leaves:
/// exp(sum_b w_b log h_b / sum_b w_b). Merging positive factors arithmetically
/// would average over GM-natural quantities and bias the merged s^2 high;
/// exp of a finite mean of logs is positive by construction, which is what
/// keeps stateIsValid's strict positivity and formMeanWeights' division safe
/// structurally. Every h_b a live tree holds is a drawn scale, hence positive.
struct GeometricMerge {
  static double toMergeSpace(double param) { return std::log(param); }
  static double fromMergeSpace(double merged) { return std::exp(merged); }
};

class Tree {
public:
  std::vector<Node> nodes;
  index_t* indices = nullptr;  // external buffer, length numObservations
  // Pooled categorical masks: rules on columns of more than 63 categories
  // store offsets into this pool. Entries are immutable once a rule
  // references them, so rule copies (snapshots, node restores) alias
  // freely; moves truncate to a mark on rejection and the chain compacts
  // between moves.
  std::vector<std::uint64_t> maskPool;

  void initialize(index_t* indexBuffer, size_t numObservations) {
    indices = indexBuffer;
    nodes.clear();
    freePairs.clear();
    maskPool.clear();
    maskPoolHighWater_ = minMaskPoolCompactionSize;
    Node root;
    root.begin = 0;
    root.end = numObservations;
    nodes.push_back(root);
    for (size_t i = 0; i < numObservations; ++i)
      indices[i] = static_cast<index_t>(i);
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
      return maskTestBit(maskWordsFor(rule), data.categoryCounts[j]);
    return rule.missingGoesRight();
  }
  bool rulesAreEqual(const ColumnStore& data, const Rule& a,
                     const Rule& b) const {
    if (a.variableIndex != b.variableIndex) return false;
    if (a.variableIndex != invalidVariable &&
        data.columnIsPooled(static_cast<size_t>(a.variableIndex)))
      return maskEquals(
        maskWordsFor(a), maskWordsFor(b),
        maskWordsForCount(
          data.categoryCounts[static_cast<size_t>(a.variableIndex)]));
    return a.bits == b.bits;
  }

  Node& at(int32_t i) { return nodes[static_cast<size_t>(i)]; }
  const Node& at(int32_t i) const { return nodes[static_cast<size_t>(i)]; }

  /// Emptiness as the branch log-likelihood's veto means it: a leaf is empty
  /// when no member carries positive weight. With no weight vector installed
  /// this IS the member count, bit for bit the test the veto has always run;
  /// with one, a zero-weight row is absent from the likelihood rather than
  /// downweighted, so a leaf of only such rows carries
  /// nothing to estimate a parameter from and its branch must be vetoed. The
  /// scan stops at the first positive weight, so an ordinary leaf costs one
  /// gather; only a leaf that is about to be vetoed walks its members.
  bool leafHasNoWeight(int32_t i, const double* weights) const {
    const Node& node(at(i));
    if (weights == nullptr) return node.numObservations() == 0;
    for (size_t j = node.begin; j < node.end; ++j)
      if (weights[indices[j]] > 0.0) return false;
    return true;
  }
  /// The leaf's rank under the branch veto: 2 when it holds no member at all,
  /// 1 when it holds members but none of positive weight, 0 when a likelihood
  /// term reaches it. The two vetoed levels stay apart because they answer to
  /// different laws. Level 2 is the MEMBERSHIP law every site outside the move
  /// kernels enforces (bottomNodesAreOccupied), so no move may install one
  /// even from a state that is already vetoed; level 1 is the move kernels'
  /// own law and IS reachable in the current state, since weights do not ride
  /// the tree and any weight install can strand a leaf, so a move out of such
  /// a state must be priced rather than compared against a like penalty.
  int leafVetoRank(int32_t i, const double* weights) const {
    if (at(i).numObservations() == 0) return 2;
    return leafHasNoWeight(i, weights) ? 1 : 0;
  }
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

  // Node collectors walk left-first.
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
  /// the split nearest the node on each side wins.
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
      data.categoryCounts[static_cast<size_t>(variableIndex)];
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
    std::uint32_t numCategories = data.categoryCounts[j];
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

  /// Fill `scratch` with variableIndex's reachable-category words at nodeIndex
  /// and return the word count. Pooled columns only; scratch is the caller's,
  /// so a query on a const tree writes no shared state.
  size_t reachableCategoryWords(const ColumnStore& data, int32_t nodeIndex,
                                int32_t variableIndex,
                                std::vector<std::uint64_t>& scratch) const {
    size_t numWords = maskWordsForCount(
      data.categoryCounts[static_cast<size_t>(variableIndex)]);
    scratch.resize(numWords);
    reachableCategoriesWide(data, nodeIndex, variableIndex, scratch.data());
    return numWords;
  }

  /// How many of variableIndex's categories reach nodeIndex, over either
  /// storage tier. `cachedInline` is an inline column's already-walked mask,
  /// for a caller holding one from collectAvailableVariables, or null to walk
  /// it here; a pooled column ignores it and fills `scratch` instead.
  size_t reachableCategoryCount(const ColumnStore& data, int32_t nodeIndex,
                                int32_t variableIndex,
                                std::vector<std::uint64_t>& scratch,
                                const std::uint64_t* cachedInline) const {
    if (data.columnIsPooled(static_cast<size_t>(variableIndex))) {
      size_t numWords =
        reachableCategoryWords(data, nodeIndex, variableIndex, scratch);
      return maskPopcount(scratch.data(), numWords);
    }
    return static_cast<size_t>(std::popcount(
      cachedInline != nullptr
        ? *cachedInline
        : reachableCategories(data, nodeIndex, variableIndex)));
  }

  /// The same set handed to the tier's action rather than counted:
  /// pooled(words, numWords) for a column with more than 63 categories,
  /// narrow(mask) for one the mask fits inline. Both must return the same
  /// type. Nothing here marks or truncates the mask pool - an action that
  /// allocates from it owns the rollback.
  template <typename Pooled, typename Narrow>
  auto withReachableMask(const ColumnStore& data, int32_t nodeIndex,
                         int32_t variableIndex,
                         std::vector<std::uint64_t>& scratch,
                         const std::uint64_t* cachedInline, Pooled pooled,
                         Narrow narrow) const {
    if (data.columnIsPooled(static_cast<size_t>(variableIndex))) {
      size_t numWords =
        reachableCategoryWords(data, nodeIndex, variableIndex, scratch);
      return pooled(scratch.data(), numWords);
    }
    return narrow(cachedInline != nullptr
                    ? *cachedInline
                    : reachableCategories(data, nodeIndex, variableIndex));
  }

  /// Install a per-forest split-variable restriction: a borrowed 0/1 byte per
  /// predictor (1 = splittable), or null to lift it. The availability queries
  /// short-circuit on the null before touching it, so an unrestricted tree
  /// runs the default path unchanged.
  void setColumnMask(const std::uint8_t* columnMask) { columnMask_ = columnMask; }

  /// Install a per-forest interaction constraint (or null to lift it). Only an
  /// ACTIVE constraint is retained, so an inactive one leaves the availability
  /// path byte-for-byte unchanged; the queries short-circuit on the null.
  void setInteractionConstraint(const InteractionConstraint* interaction) {
    interaction_ =
      (interaction != nullptr && interaction->active()) ? interaction : nullptr;
  }
  bool hasInteractionConstraint() const { return interaction_ != nullptr; }

  /// Whole-subtree, all-variables interaction feasibility: walk the subtree
  /// rooted at subtreeRoot carrying the running distinct-ancestor bitset
  /// (seeded from subtreeRoot's own ancestors) and reject if ANY node's split
  /// variable violates the order cap or a forbidden co-occurrence. Unlike the
  /// per-variable categoricalSubtreeIsValid / ordinalRuleIsValid, this couples
  /// DIFFERENT variables: a swap that re-checked only the swapped pair could
  /// miss a sibling-path violation, so every node must be tested. Trivially
  /// true when the constraint is inactive, so the change/swap moves may call
  /// it unconditionally.
  bool interactionSubtreeIsValid(int32_t subtreeRoot) const {
    if (interaction_ == nullptr) return true;
    interactionWalkScratch_.resize(interaction_->numWords);
    collectAncestorVariables(subtreeRoot, interactionWalkScratch_.data());
    return interactionSubtreeWalk(subtreeRoot, interactionWalkScratch_.data());
  }

  /// Whole-subtree column-mask feasibility: reject if ANY decision node in the
  /// subtree rooted at subtreeRoot splits on a variable this tree's column mask
  /// forbids. Unlike interactionSubtreeIsValid the test is per-node independent
  /// - no ancestor set, no order - so a flat scan of the decision nodes decides
  /// it. Trivially true for an unrestricted tree (columnMask_ null short-circuit),
  /// so the default availability path is byte-for-byte unchanged and a warm start
  /// may call it unconditionally.
  bool columnMaskSubtreeIsValid(int32_t subtreeRoot) const {
    if (columnMask_ == nullptr) return true;
    return columnMaskSubtreeWalk(subtreeRoot);
  }

  bool variableAvailable(const ColumnStore& data, int32_t nodeIndex,
                         int32_t variableIndex) const {
    if (!columnAllowed(static_cast<size_t>(variableIndex))) return false;
    bool cutAvailable;
    if (data.splitsBySubset(static_cast<size_t>(variableIndex))) {
      cutAvailable = reachableCategoryCount(data, nodeIndex, variableIndex,
                                            reachableScratch_, nullptr) >= 2;
    } else {
      int32_t left, right;
      splitInterval(data, nodeIndex, variableIndex, &left, &right);
      cutAvailable = right >= left;
    }
    if (!cutAvailable) return false;
    if (interaction_ != nullptr &&
        !interactionVariableAvailable(nodeIndex,
                                      static_cast<size_t>(variableIndex)))
      return false;
    return true;
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
    // interaction constraint: accumulate the distinct ancestor split-variable
    // set along the SAME walk, adding only O(p) rather than a separate pass.
    // Guarded so the unconstrained path never touches the scratch.
    const size_t interactionWords =
      interaction_ != nullptr ? interaction_->numWords : 0;
    if (interaction_ != nullptr) {
      availAncestorScratch_.resize(interactionWords);
      for (size_t w = 0; w < interactionWords; ++w) availAncestorScratch_[w] = 0;
    }
    for (size_t j = 0; j < p; ++j) {
      if (data.splitsBySubset(j)) {
        if (data.columnIsPooled(j)) continue;  // resolved after the walk
        std::uint32_t numCategories = data.categoryCounts[j];
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
      if (interaction_ != nullptr)
        maskSetBit(availAncestorScratch_.data(), static_cast<std::uint32_t>(j));
      if (data.splitsBySubset(j)) {
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

    const size_t order = interaction_ != nullptr
      ? maskPopcount(availAncestorScratch_.data(), interactionWords) : 0;
    size_t count = 0;
    for (size_t j = 0; j < p; ++j) {
      bool avail;
      if (data.splitsBySubset(j))
        avail = data.columnIsPooled(j)
          ? variableAvailable(data, nodeIndex, static_cast<int32_t>(j))
          : std::popcount(availMaskScratch_[j]) >= 2;
      else
        avail = availRightScratch_[j] >= availLeftScratch_[j];
      if (!columnAllowed(j)) avail = false;  // no-op when unrestricted
      // interaction: drop a variable the ancestor set forbids by order or
      // co-occurrence (idempotent for the pooled fallback above)
      if (avail && interaction_ != nullptr)
        avail = interactionAllows(j, availAncestorScratch_.data(), order);
      available[j] = avail ? 1 : 0;
      count += avail ? 1 : 0;
    }
    return count;
  }

  /// The reachable category masks collectAvailableVariables narrowed on its
  /// single walk, indexed by predictor. Valid for INLINE categorical columns
  /// only (the walk skips pooled ones, which need reachableCategoriesWide) and
  /// only until the next availability query on this tree, so a caller that
  /// already collected availability at a node reads R here instead of paying a
  /// second O(p * depth) walk for it.
  const std::uint64_t* inlineReachableMasks() const {
    return availMaskScratch_.data();
  }

  /// Leaf sufficient statistic (sum w, sum wz) in one pass. The root
  /// intentionally uses the non-indexed kernels (identical values, cheaper
  /// access). Templated on the residual element type: the default double path
  /// selects the fp64 kernels (byte-identical to before), the opt-in fp32
  /// residual (ResidT = float) selects the float-input kernels that load float
  /// and accumulate double.
  template <typename ResidT>
  void computeLeafStats(int32_t nodeIndex, const ResidT* y, const double* weights) {
    Node& node(at(nodeIndex));
    bool isRoot = node.parent == invalidNode;
    if (weights == nullptr) {
      if (isRoot) {
        if constexpr (std::is_same_v<ResidT, float>)
          misc_computeFloatSufficientStatisticsFast(
            y, node.numObservations(), &node.sumWeights,
            &node.sumWeightedResponse);
        else
          misc_computeSufficientStatisticsFast(
            y, node.numObservations(), &node.sumWeights,
            &node.sumWeightedResponse);
      } else {
        if constexpr (std::is_same_v<ResidT, float>)
          misc_computeIndexedFloatSufficientStatisticsFast(
            y, indices + node.begin, node.numObservations(), &node.sumWeights,
            &node.sumWeightedResponse);
        else
          misc_computeIndexedSufficientStatisticsFast(
            y, indices + node.begin, node.numObservations(), &node.sumWeights,
            &node.sumWeightedResponse);
      }
    } else {
      if (isRoot) {
        if constexpr (std::is_same_v<ResidT, float>)
          misc_computeWeightedFloatSufficientStatisticsFast(
            y, node.numObservations(), weights, &node.sumWeights,
            &node.sumWeightedResponse);
        else
          misc_computeWeightedSufficientStatisticsFast(
            y, node.numObservations(), weights, &node.sumWeights,
            &node.sumWeightedResponse);
      } else {
        if constexpr (std::is_same_v<ResidT, float>)
          misc_computeIndexedWeightedFloatSufficientStatisticsFast(
            y, indices + node.begin, node.numObservations(), weights,
            &node.sumWeights, &node.sumWeightedResponse);
        else
          misc_computeIndexedWeightedSufficientStatisticsFast(
            y, indices + node.begin, node.numObservations(), weights,
            &node.sumWeights, &node.sumWeightedResponse);
      }
    }
  }

  template <typename ResidT>
  void setNodeAverages(const ResidT* y, const double* weights) {
    bottomScratch.clear();
    fillBottom(0, bottomScratch);
    for (int32_t i : bottomScratch) computeLeafStats(i, y, weights);
  }

  /// Shared two-pointer in-place partition kernel for the scalar fallbacks:
  /// observations for which goesRight is false are gathered into [0, result),
  /// the rest into [result, length). goesRight receives the observation index;
  /// each typed wrapper below supplies the storage- and rule-specific
  /// predicate. The swap sequence is identical to the previously inlined form,
  /// so every wrapper's index permutation - and thus the downstream
  /// sufficient-statistic summation order - is byte-for-byte unchanged.
  template <typename GoesRight>
  static size_t partitionByPredicate(index_t* indices, size_t length,
                                     GoesRight goesRight) {
    size_t lo = 0, hi = length;
    // invariant: [0, lo) is left-bound, [hi, length) is right-bound
    while (true) {
      while (lo < hi && !goesRight(indices[lo])) ++lo;
      while (lo < hi && goesRight(indices[hi - 1])) --hi;
      if (hi - lo < 2) break;
      index_t temp = indices[lo];
      indices[lo] = indices[hi - 1];
      indices[hi - 1] = temp;
      ++lo;
      --hi;
    }
    return lo;
  }

  /// How a scalar partition reads a column's codes: a contiguous dense array
  /// or the rank-bitmap layout.
  enum class CodeStorage { dense, sparse };
  template <CodeStorage storage>
  using PartitionColumn =
    std::conditional_t<storage == CodeStorage::dense, const xint_t*,
                       const SparseColumnData&>;

  /// What a scalar partition tests a code against, and with it the type of the
  /// rule payload the caller supplies: an inline 64-bit category mask by
  /// value, a pooled rule's mask words, or the Rule itself for the ordinal
  /// split with a direction for the reserved missing code.
  enum class PartitionRule { inlineMask, pooledMask, missingAware };

  /// The missing-aware test's operands, decoded from the rule once rather than
  /// per observation.
  struct MissingAwareSplit {
    int32_t splitIndex;
    bool missingGoesRight;
  };

  /// The scalar partitions the SIMD kernels do not cover, over the (storage,
  /// rule) cross product: category membership on either storage, and the
  /// missing-aware ordinal compare on either. The inline-mask arms carry the
  /// shift's code < 64 invariant, which is what confines them to a column of
  /// at most 63 categories plus the reserved missing position.
  ///
  /// The swap sequence is partitionByPredicate's, so every instantiation
  /// produces the same index permutation - and thus the same downstream
  /// sufficient-statistic summation order - as a hand-written loop would.
  template <CodeStorage storage, PartitionRule rule, typename Payload>
  static size_t partitionIndicesScalar(PartitionColumn<storage> column,
                                       const Payload& payload,
                                       index_t* indices, size_t length) {
    auto test = [&] {
      if constexpr (rule == PartitionRule::missingAware)
        return MissingAwareSplit{payload.splitIndex(),
                                 payload.missingGoesRight()};
      else
        return payload;  // a mask arrives already decoded
    }();
    return partitionByPredicate(indices, length, [&](index_t i) {
      xint_t code;
      if constexpr (storage == CodeStorage::dense) code = column[i];
      else code = column.at(i);
      if constexpr (rule == PartitionRule::inlineMask)
        return ((test >> code) & 1u) != 0;
      else if constexpr (rule == PartitionRule::pooledMask)
        return maskTestBit(test, code);
      else
        return code == naCode ? test.missingGoesRight
                              : static_cast<int32_t>(code) > test.splitIndex;
    });
  }

  // pins the misc partition kernels' integer width; a mismatch would
  // silently truncate the cut-index casts below
  static_assert(std::is_same_v<misc_xint_t, std::uint16_t>);

  // pins the C++ gather-index buffer element type to the C kernel index type,
  // so passing index_t* where the retyped misc kernels expect misc_index_t*
  // (and the sizeof memcpy/memcmp over index segments) stays byte-exact; a
  // width mismatch would silently truncate or overread indices.
  static_assert(sizeof(index_t) == sizeof(misc_index_t));

  /// Partition a node's observations between its children by its rule.
  void partitionChildren(const ColumnStore& data, int32_t nodeIndex) {
    Node& node(at(nodeIndex));
    Node& left(at(node.leftChild));
    Node& right(at(node.leftChild + 1));

    size_t numOnLeft = 0;
    if (node.numObservations() > 0) {
      size_t variable = static_cast<size_t>(node.rule.variableIndex);
      index_t* segment = indices + node.begin;
      size_t numMembers = node.numObservations();
      if (data.splitsBySubset(variable)) {
        bool pooled = data.columnIsPooled(variable);
        if (data.columnIsSparse(variable)) {
          const SparseColumnData& column = data.sparseColumn(variable);
          numOnLeft = pooled
            ? partitionIndicesScalar<CodeStorage::sparse,
                                     PartitionRule::pooledMask>(
                column, maskWordsFor(node.rule), segment, numMembers)
            : partitionIndicesScalar<CodeStorage::sparse,
                                     PartitionRule::inlineMask>(
                column, node.rule.categoryDirections(), segment, numMembers);
        } else {
          const xint_t* column = data.column(variable);
          numOnLeft = pooled
            ? partitionIndicesScalar<CodeStorage::dense,
                                     PartitionRule::pooledMask>(
                column, maskWordsFor(node.rule), segment, numMembers)
            : partitionIndicesScalar<CodeStorage::dense,
                                     PartitionRule::inlineMask>(
                column, node.rule.categoryDirections(), segment, numMembers);
        }
      } else if (data.columnIsSparse(variable)) {
        // in-place partition at the root too: misc_partitionRange assumes
        // identity index content, which only the dense path maintains
        const SparseColumnData& column = data.sparseColumn(variable);
        if (data.hasMissing[variable]) {
          numOnLeft = partitionIndicesScalar<CodeStorage::sparse,
                                             PartitionRule::missingAware>(
            column, node.rule, segment, numMembers);
        } else {
          numOnLeft = misc_partitionIndicesSparse(
            column.bits.data(), column.wordRanks.data(),
            column.nzCodes.data(), column.zeroCode,
            static_cast<misc_xint_t>(node.rule.splitIndex()), segment,
            numMembers);
        }
      } else {
        const xint_t* column = data.column(variable);
        bool isRoot = node.parent == invalidNode;
        if (data.hasMissing[variable]) {
          numOnLeft = partitionIndicesScalar<CodeStorage::dense,
                                             PartitionRule::missingAware>(
            column, node.rule, segment, numMembers);
        } else {
          numOnLeft = isRoot
            ? misc_partitionRange(column, static_cast<misc_xint_t>(node.rule.splitIndex()),
                                  segment, numMembers)
            : misc_partitionIndices(column, static_cast<misc_xint_t>(node.rule.splitIndex()),
                                    segment, numMembers);
        }
      }
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

  /// Validity criterion after a predictor change: no bottom node may be left
  /// without observations.
  bool bottomNodesAreOccupied() const {
    return bottomNodesAreOccupiedBelow(0);
  }

  /// Whether the tree is in the set the branch log-likelihood's veto admits:
  /// no bottom node fails leafHasNoWeight. This is the emptiness law of the
  /// move kernels, not the membership law
  /// bottomNodesAreOccupied answers for state restore; the two agree exactly
  /// when no weight vector is installed. leafVetoRank splits the failure into
  /// the two levels the moves order lexicographically: a tree can leave this
  /// set through a weight install (rank 1, transient - the moves price their
  /// way out) but never through a move (rank 2 is refused absolutely).
  bool bottomNodesHaveWeight(const double* weights) const {
    return bottomNodesHaveWeightBelow(0, weights);
  }

  /// Repartition a subtree after its rule changed, recomputing leaf stats.
  template <typename ResidT>
  void refreshSubtree(const ColumnStore& data, int32_t nodeIndex, const ResidT* y,
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
  template <typename ResidT>
  void birth(const ColumnStore& data, int32_t nodeIndex, const Rule& rule,
             const ResidT* y, const double* weights) {
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
    node.sumWeights = left.sumWeights + right.sumWeights;
    node.sumWeightedResponse =
      left.sumWeightedResponse + right.sumWeightedResponse;
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
    std::vector<index_t> indexSegment;
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
                snapshot.indexSegment.size() * sizeof(index_t));
  }

  /// Descend a test row through the storage-aware test accessor, reading only
  /// the columns on the path (test prediction).
  int32_t findBottomNodeForRow(const ColumnStore& data, size_t i) const {
    int32_t current = 0;
    while (!at(current).isBottom()) {
      const Rule& rule(at(current).rule);
      xint_t code =
        data.testCodeAt(static_cast<size_t>(rule.variableIndex), i);
      current = ruleSendsRight(data, rule, code) ? at(current).leftChild + 1
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

  /// Collapse any node with an unoccupied child, or an ordinal split the
  /// current grid can no longer address (setCutPoints shrank the column below
  /// its index; missing values riding the right child keep both children
  /// occupied so the empty-child test alone would spare it), into a leaf whose
  /// parameter is the effective-observation-weighted mean of its subtree's
  /// leaf parameters, for forced predictor updates. paramByNode is indexed by
  /// arena id, paramStride doubles per node, merged per coordinate; a subtree
  /// with no observations at all gets the plain mean. Merge selects the space
  /// the mean is taken in (arithmetic for an additive leaf parameter,
  /// geometric for a multiplicative scale).
  template <typename Merge = ArithmeticMerge>
  void collapseEmptyNodes(const ColumnStore& data, const double* weights,
                          std::vector<double>& paramByNode,
                          size_t paramStride = 1) {
    collapseEmptyNodesBelow<Merge>(0, data, weights, paramByNode, paramStride);
  }

  /// An ordinal rule whose index no longer addresses a cut after the grid
  /// shrank; flatten would read past cutPoints[j]. Categorical masks stay
  /// representable (setCutPoints leaves category counts alone).
  bool ruleIsUnrepresentable(const ColumnStore& data, const Rule& rule) const {
    size_t j = static_cast<size_t>(rule.variableIndex);
    return data.splitsByThreshold(j) &&
           rule.splitIndex() >= static_cast<int32_t>(data.numCuts[j]);
  }

  /// Point the tree at a new observation buffer (whole-data replacement,
  /// possibly of a different size): identity indices under the root, all
  /// partitions stale until repartitionSubtree.
  void resetObservations(index_t* indexBuffer, size_t numObservations) {
    indices = indexBuffer;
    at(0).begin = 0;
    at(0).end = numObservations;
    for (size_t i = 0; i < numObservations; ++i)
      indices[i] = static_cast<index_t>(i);
  }

  /// After a from-scratch cut rebuild, remap every rule's splitIndex onto the
  /// new cut whose value is nearest the old cut, restricted to the ancestor-
  /// constrained interval; a subtree whose interval empties collapses to a
  /// leaf with the plain mean of its leaf parameters. paramByNode is indexed
  /// by arena id, paramStride doubles per node, merged per coordinate. Merge
  /// selects the space the mean is taken in, as for collapseEmptyNodes.
  template <typename Merge = ArithmeticMerge>
  void mapOldCutPointsOntoNew(const ColumnStore& data,
                              const std::vector<std::vector<double>>& oldCutPoints,
                              std::vector<double>& paramByNode,
                              size_t paramStride = 1) {
    std::vector<int32_t> minIndices(data.numPredictors, 0);
    std::vector<int32_t> maxIndices(data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j)
      maxIndices[j] = static_cast<int32_t>(data.numCuts[j]);
    mapCutPointsBelow<Merge>(0, data, oldCutPoints, paramByNode,
                             minIndices.data(), maxIndices.data(), paramStride);
  }

  /// Drop a rule's missing direction after a data mutation stops its column
  /// routing a missing value: hasMissing false puts the bit outside
  /// reachableCategories, so buildFromFlat's gauge check would reject an
  /// otherwise valid rule. The bit routes nothing without missing
  /// observations, so clearing it moves nothing. Pooled masks keep the bit in
  /// their words under their own scheme, so pass through.
  void dropStaleMissingDirections(const ColumnStore& data) {
    dropStaleMissingDirectionsBelow(0, data);
  }

  void countVariableUses(std::uint32_t* counts) const {
    countVariableUsesBelow(0, counts);
  }

  /// Info dump; the per-node output format is R-visible and pinned by tests.
  /// One line per node in pre-order, indented by depth, with occupancy,
  /// top/bottom flags, per-variable availability, and the rule or the leaf
  /// parameter (taken from paramByNode, indexed by arena id with paramStride
  /// doubles per node, on the internal scale); vector-parameter leaves append
  /// their slopes after the intercept.
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
      if (data.splitsBySubset(variableIndex)) {
        ext_printf("CATRule: ");
        if (data.columnIsPooled(variableIndex)) {
          const std::uint64_t* directions = maskWordsFor(node.rule);
          for (std::uint32_t i = 0; i < data.categoryCounts[variableIndex];
               ++i)
            ext_printf(" %u", maskTestBit(directions, i) ? 1u : 0u);
        } else {
          for (std::uint32_t i = 0; i < data.categoryCounts[variableIndex];
               ++i)
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
  /// slopes per leaf, in pre-order, to slopes. masks receives pooled
  /// categorical rules' words (see FlatNode) and must be non-null when the
  /// store has any pooled categorical column.
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
  /// the pooled categorical side channel, whose offsets must be pre-order
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

  bool bottomNodesHaveWeightBelow(int32_t i, const double* weights) const {
    if (at(i).isBottom()) return !leafHasNoWeight(i, weights);
    return bottomNodesHaveWeightBelow(at(i).leftChild, weights) &&
           bottomNodesHaveWeightBelow(at(i).leftChild + 1, weights);
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

  void dropStaleMissingDirectionsBelow(int32_t nodeIndex,
                                       const ColumnStore& data) {
    Node& node(at(nodeIndex));
    if (node.isBottom()) return;
    size_t j = static_cast<size_t>(node.rule.variableIndex);
    if (!data.hasMissing[j] && !data.columnIsPooled(j) &&
        node.rule.missingGoesRight())
      node.rule.setMissingGoesRight(false);
    dropStaleMissingDirectionsBelow(node.leftChild, data);
    dropStaleMissingDirectionsBelow(node.leftChild + 1, data);
  }

  /// minIndices are inclusive, maxIndices exclusive; both are saved and
  /// restored around the recursion. Categorical rules have nothing to remap
  /// (category counts are fixed across data replacement) and pass through to
  /// their children.
  template <typename Merge>
  void mapCutPointsBelow(int32_t nodeIndex, const ColumnStore& data,
                         const std::vector<std::vector<double>>& oldCutPoints,
                         std::vector<double>& paramByNode,
                         int32_t* minIndices, int32_t* maxIndices,
                         size_t paramStride) {
    if (at(nodeIndex).isBottom()) return;

    int32_t varIndex = at(nodeIndex).rule.variableIndex;

    if (data.splitsBySubset(static_cast<size_t>(varIndex))) {
      mapCutPointsBelow<Merge>(at(nodeIndex).leftChild, data, oldCutPoints,
                               paramByNode, minIndices, maxIndices,
                               paramStride);
      mapCutPointsBelow<Merge>(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                               paramByNode, minIndices, maxIndices,
                               paramStride);
      return;
    }

    int32_t minIndex = minIndices[varIndex];
    int32_t maxIndex = maxIndices[varIndex];

    if (minIndex > maxIndex - 1) {
      // no split of this variable remains below the ancestors: merge the
      // subtree's leaf parameters by effective observation count, matching
      // collapseEmptyNodesBelow; plain mean when the subtree holds no weight
      std::vector<int32_t> bottoms;
      fillBottom(nodeIndex, bottoms);
      double weightTotal = 0.0;
      std::vector<double> paramTotals(paramStride, 0.0);
      std::vector<double> paramSums(paramStride, 0.0);
      for (int32_t i : bottoms) {
        double weight = at(i).sumWeights;
        weightTotal += weight;
        const double* params =
          paramByNode.data() + static_cast<size_t>(i) * paramStride;
        for (size_t j = 0; j < paramStride; ++j) {
          double param = Merge::toMergeSpace(params[j]);
          paramTotals[j] += weight * param;
          paramSums[j] += param;
        }
      }
      double* merged =
        paramByNode.data() + static_cast<size_t>(nodeIndex) * paramStride;
      for (size_t j = 0; j < paramStride; ++j)
        merged[j] = Merge::fromMergeSpace(
          weightTotal > 0.0
            ? paramTotals[j] / weightTotal
            : paramSums[j] / static_cast<double>(bottoms.size()));
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
    mapCutPointsBelow<Merge>(at(nodeIndex).leftChild, data, oldCutPoints,
                             paramByNode, minIndices, maxIndices, paramStride);
    maxIndices[varIndex] = maxIndex;

    minIndices[varIndex] = newIndex + 1;
    mapCutPointsBelow<Merge>(at(nodeIndex).leftChild + 1, data, oldCutPoints,
                             paramByNode, minIndices, maxIndices, paramStride);
    minIndices[varIndex] = minIndex;
  }

  template <typename Merge>
  void collapseEmptyNodesBelow(int32_t nodeIndex, const ColumnStore& data,
                               const double* weights,
                               std::vector<double>& paramByNode,
                               size_t paramStride) {
    if (at(nodeIndex).isBottom()) return;

    if (at(at(nodeIndex).leftChild).numObservations() == 0 ||
        at(at(nodeIndex).leftChild + 1).numObservations() == 0 ||
        ruleIsUnrepresentable(data, at(nodeIndex).rule)) {
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
          double param = Merge::toMergeSpace(params[j]);
          paramTotals[j] += weight * param;
          paramSums[j] += param;
        }
      }
      double* merged =
        paramByNode.data() + static_cast<size_t>(nodeIndex) * paramStride;
      for (size_t j = 0; j < paramStride; ++j)
        merged[j] = Merge::fromMergeSpace(
          weightTotal > 0.0
            ? paramTotals[j] / weightTotal
            : paramSums[j] / static_cast<double>(bottoms.size()));

      collapseSubtreeToLeaf(nodeIndex);
    } else {
      collapseEmptyNodesBelow<Merge>(at(nodeIndex).leftChild, data, weights,
                                     paramByNode, paramStride);
      collapseEmptyNodesBelow<Merge>(at(nodeIndex).leftChild + 1, data, weights,
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
    FlatKind kind = expectedFlatKind(data, j);
    setFlatKind(flat, kind);
    if (kind == FlatKind::ordinal) {
      flat.value = data.cutPoints[j]
                                 [static_cast<size_t>(node.rule.splitIndex())];
    } else if (kind == FlatKind::categoricalPooled) {
      // side channel: maskOffset points at numMaskWords words of category
      // bits only; the missing direction stays in flags for either kind
      std::uint32_t numCategories = data.categoryCounts[j];
      size_t numWords = maskWordsForCount(numCategories);
      size_t offset = masks->size();
      flat.maskOffset = static_cast<std::uint64_t>(offset);
      flat.numMaskWords = static_cast<std::uint32_t>(numWords);
      masks->resize(offset + numWords, 0);
      std::uint64_t* words = masks->data() + offset;
      const std::uint64_t* pooled = maskWordsFor(node.rule);
      for (size_t w = 0; w < numWords; ++w) words[w] = pooled[w];
      words[numCategories >> 6] &= ~(1ull << (numCategories & 63u));
    } else {
      flat.mask = node.rule.categoryDirections() & ~Rule::missingDirectionBit;
    }
    if (ruleMissingGoesRight(data, node.rule))
      flat.flags |= flatMissingGoesRight;
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
        (flat.flags & ~(flatMissingGoesRight | flatKindMask)) != 0)
      return false;

    Rule rule;
    rule.variableIndex = flat.variable;
    size_t variable = static_cast<size_t>(flat.variable);
    FlatKind kind = flatKindOf(flat);
    if (kind != expectedFlatKind(data, variable)) return false;
    if (kind == FlatKind::categoricalPooled) {
      // side channel: the offset must be the running pre-order word cursor,
      // the words category bits only, the assembled mask a canonical gauge
      std::uint32_t numCategories = data.categoryCounts[variable];
      size_t numWords = maskWordsForCount(numCategories);
      if (masks == nullptr || flat.maskOffset != maskPos ||
          maskPos + numWords > numMaskWords)
        return false;
      const std::uint64_t* words = masks + maskPos;
      maskPos += numWords;
      for (std::uint32_t bit = numCategories;
           bit < static_cast<std::uint32_t>(64 * numWords); ++bit)
        if (maskTestBit(words, bit)) return false;
      size_t offset = allocateMask(numWords);
      std::uint64_t* directions = mutableMaskWordsFor(offset);
      std::memcpy(directions, words, numWords * sizeof(std::uint64_t));
      if ((flat.flags & flatMissingGoesRight) != 0)
        maskSetBit(directions, numCategories);
      reachableCategoryWords(data, nodeIndex, flat.variable, reachableScratch_);
      if (maskIsZero(directions, numWords) ||
          !maskIsSubsetOf(directions, reachableScratch_.data(), numWords) ||
          maskEquals(directions, reachableScratch_.data(), numWords))
        return false;
      rule.setMaskOffset(offset);
    } else if (kind == FlatKind::categoricalInline) {
      // an inline mask over the observed categories with no bit past them;
      // the missing direction arrives in flags
      std::uint32_t numCategories = data.categoryCounts[variable];
      if ((flat.mask >> numCategories) != 0) return false;
      std::uint64_t directions = flat.mask;
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
      size_t numWords = maskWordsForCount(data.categoryCounts[j]);
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

  // per-forest split-variable restriction (chain.hpp installs it); null (the
  // default) restricts nothing. The availability queries short-circuit on the
  // null, so the unrestricted path never touches the mask.
  const std::uint8_t* columnMask_ = nullptr;
  bool columnAllowed(size_t j) const {
    return columnMask_ == nullptr || columnMask_[j] != 0;
  }

  // per-forest interaction constraint (chain.hpp installs it, only when
  // active); null on the default path. availAncestorScratch_ carries the
  // ancestor variable-set along collectAvailableVariables' single walk;
  // interactionWalkScratch_ is the running bitset of the subtree-validity walk.
  const InteractionConstraint* interaction_ = nullptr;
  mutable std::vector<std::uint64_t> availAncestorScratch_;
  mutable std::vector<std::uint64_t> interactionWalkScratch_;

  /// The admissibility test p(j | A): a NEW variable j is barred once the
  /// distinct-ancestor order reaches K (max-order), and any variable is barred
  /// when a forbidden partner sits in the ancestor set (co-occurrence).
  bool interactionAllows(size_t j, const std::uint64_t* ancestors,
                         size_t order) const {
    if (interaction_->hasOrderCap()) {
      bool alreadyUsed = maskTestBit(ancestors, static_cast<std::uint32_t>(j));
      if (!alreadyUsed && order >= interaction_->maxOrder) return false;
    }
    if (interaction_->hasForbidden()) {
      const std::uint64_t* row = interaction_->forbiddenRow(j);
      for (size_t w = 0; w < interaction_->numWords; ++w)
        if ((row[w] & ancestors[w]) != 0) return false;
    }
    return true;
  }

  /// Distinct split variables STRICTLY above nodeIndex, written into ancestors
  /// (interaction_->numWords words); returns the count (the node's order).
  size_t collectAncestorVariables(int32_t nodeIndex,
                                  std::uint64_t* ancestors) const {
    size_t numWords = interaction_->numWords;
    for (size_t w = 0; w < numWords; ++w) ancestors[w] = 0;
    int32_t current = nodeIndex;
    while (at(current).parent != invalidNode) {
      current = at(current).parent;
      maskSetBit(ancestors,
                 static_cast<std::uint32_t>(at(current).rule.variableIndex));
    }
    return maskPopcount(ancestors, numWords);
  }

  /// Single-variable interaction feasibility at nodeIndex, computing the
  /// ancestor set on the fly (the per-call path variableAvailable /
  /// hasAnyAvailableVariable take).
  bool interactionVariableAvailable(int32_t nodeIndex, size_t j) const {
    interactionWalkScratch_.resize(interaction_->numWords);
    size_t order =
      collectAncestorVariables(nodeIndex, interactionWalkScratch_.data());
    return interactionAllows(j, interactionWalkScratch_.data(), order);
  }

  /// DFS of interactionSubtreeIsValid: test each internal node's variable
  /// against the running ancestor bitset, then recurse with the variable
  /// added, backtracking on the way out.
  bool interactionSubtreeWalk(int32_t nodeIndex,
                              std::uint64_t* ancestors) const {
    const Node& node(at(nodeIndex));
    if (node.isBottom()) return true;
    std::uint32_t j = static_cast<std::uint32_t>(node.rule.variableIndex);
    size_t order = maskPopcount(ancestors, interaction_->numWords);
    if (!interactionAllows(j, ancestors, order)) return false;
    bool wasSet = maskTestBit(ancestors, j);
    if (!wasSet) maskSetBit(ancestors, j);
    bool ok = interactionSubtreeWalk(node.leftChild, ancestors) &&
              interactionSubtreeWalk(node.leftChild + 1, ancestors);
    if (!wasSet) maskClearBit(ancestors, j);
    return ok;
  }

  /// DFS of columnMaskSubtreeIsValid: a decision node fails when its split
  /// variable lies outside the tree's column mask; leaves always pass. No
  /// ancestor bookkeeping - each node is judged on its own variable alone.
  bool columnMaskSubtreeWalk(int32_t nodeIndex) const {
    const Node& node(at(nodeIndex));
    if (node.isBottom()) return true;
    if (!columnAllowed(static_cast<size_t>(node.rule.variableIndex)))
      return false;
    return columnMaskSubtreeWalk(node.leftChild) &&
           columnMaskSubtreeWalk(node.leftChild + 1);
  }
};

/// The reader a dense column-major block hands the flat replay: a bare
/// pointer at a column's first row, so at(row) is one indexed load.
struct DenseColumnReader {
  const double* values;
  double at(size_t row) const { return values[row]; }
};

/// A borrowed dense column-major block in the Columns shape the flat replay
/// reads predictors through: column(j) yields a reader over store column j.
/// Every raw read in a replay - the routing hoist and the leaf covariate
/// loads, which index by the store column - goes through column(), so a
/// source that maps store columns onto other storage substitutes here
/// without a second code path.
struct DenseColumns {
  const double* values;
  size_t numRows;
  DenseColumnReader column(size_t j) const {
    return DenseColumnReader{values + j * numRows};
  }
};

/// Partition indices[lo, hi) of a Columns predictor source around a
/// flattened split so left-bound rows precede right-bound ones, returning
/// the boundary: ordinal rows go left when x <= value, categorical rows when
/// the mask's direction bit for their code is clear. Order within the halves
/// is not preserved (nor needed; replays only count or accumulate). The tag
/// selects the payload; maskWords holds a pooled rule's numMaskWords words at
/// maskOffset and may be null when no rule is pooled.
template <typename Columns>
inline size_t partitionFlatIndices(const FlatNode& flat, const Columns& x,
                                   size_t* indices, size_t lo, size_t hi,
                                   const std::uint64_t* maskWords = nullptr) {
  auto column = x.column(static_cast<size_t>(flat.variable));
  size_t mid = lo;
  bool missingGoesLeft = (flat.flags & flatMissingGoesRight) == 0;
  FlatKind kind = flatKindOf(flat);
  if (kind == FlatKind::categoricalPooled) {
    const std::uint64_t* directions = maskWords + flat.maskOffset;
    // callers validate codes against the training categories; the clamps
    // here and below keep an out-of-range code a defined lookup instead of
    // undefined behavior, so the tree layer is safe standalone
    std::uint32_t maxCode = 64u * flat.numMaskWords - 1u;
    for (size_t k = lo; k < hi; ++k) {
      double value = column.at(indices[k]);
      bool goesLeft = isNA(value)
        ? missingGoesLeft
        : !maskTestBit(directions,
                       std::min(static_cast<std::uint32_t>(value), maxCode));
      if (goesLeft) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  } else if (kind == FlatKind::categoricalInline) {
    std::uint64_t directions = flat.mask;
    for (size_t k = lo; k < hi; ++k) {
      double value = column.at(indices[k]);
      bool goesLeft = isNA(value)
        ? missingGoesLeft
        : ((directions >>
            (static_cast<std::uint32_t>(value) & 63u)) & 1u) == 0;
      if (goesLeft) {
        size_t temp = indices[mid];
        indices[mid] = indices[k];
        indices[k] = temp;
        ++mid;
      }
    }
  } else {
    for (size_t k = lo; k < hi; ++k) {
      double value = column.at(indices[k]);
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

/// The raw column-major entry: numRows is the block's row stride.
inline size_t partitionFlatIndices(const FlatNode& flat, const double* x,
                                   size_t numRows, size_t* indices, size_t lo,
                                   size_t hi,
                                   const std::uint64_t* maskWords = nullptr) {
  return partitionFlatIndices(flat, DenseColumns{x, numRows}, indices, lo, hi,
                              maskWords);
}

/// Route indices[lo, hi) of a Columns predictor source through a flattened
/// subtree, writing each node's routed count in pre-order; indices is
/// scrambled. Returns the number of flattened nodes consumed.
template <typename Columns>
inline size_t countFlatObservationsBelow(const FlatNode* flatNodes,
                                         const Columns& x, size_t* indices,
                                         size_t lo, size_t hi,
                                         std::uint32_t* counts,
                                         const std::uint64_t* maskWords = nullptr) {
  counts[0] = static_cast<std::uint32_t>(hi - lo);
  if (flatNodes[0].variable == invalidVariable) return 1;

  size_t mid =
    partitionFlatIndices(flatNodes[0], x, indices, lo, hi, maskWords);
  size_t numNodes = 1;
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, x, indices, lo,
                                         mid, counts + numNodes, maskWords);
  numNodes += countFlatObservationsBelow(flatNodes + numNodes, x, indices, mid,
                                         hi, counts + numNodes, maskWords);
  return numNodes;
}

/// The raw column-major entry: numRows is the block's row stride.
inline size_t countFlatObservationsBelow(const FlatNode* flatNodes,
                                         const double* x, size_t numRows,
                                         size_t* indices, size_t lo, size_t hi,
                                         std::uint32_t* counts,
                                         const std::uint64_t* maskWords = nullptr) {
  return countFlatObservationsBelow(flatNodes, DenseColumns{x, numRows},
                                    indices, lo, hi, counts, maskWords);
}

/// Add each routed row's leaf parameter into fits (one slot per row).
/// Returns the number of flattened nodes consumed.
template <typename Columns>
inline size_t addFlatPredictionsBelow(const FlatNode* flatNodes,
                                      const Columns& x, size_t* indices,
                                      size_t lo, size_t hi, double* fits,
                                      const std::uint64_t* maskWords = nullptr) {
  if (flatNodes[0].variable == invalidVariable) {
    for (size_t k = lo; k < hi; ++k) fits[indices[k]] += flatNodes[0].value;
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], x, indices, lo, hi, maskWords);
  size_t numNodes = 1;
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, x, indices, lo,
                                      mid, fits, maskWords);
  numNodes += addFlatPredictionsBelow(flatNodes + numNodes, x, indices, mid,
                                      hi, fits, maskWords);
  return numNodes;
}

/// The raw column-major entry: numRows is the block's row stride.
inline size_t addFlatPredictionsBelow(const FlatNode* flatNodes,
                                      const double* x, size_t numRows,
                                      size_t* indices, size_t lo, size_t hi,
                                      double* fits,
                                      const std::uint64_t* maskWords = nullptr) {
  return addFlatPredictionsBelow(flatNodes, DenseColumns{x, numRows}, indices,
                                 lo, hi, fits, maskWords);
}

/// The linear-leaf analogue of addFlatPredictionsBelow: a routed row's fit
/// is the leaf's intercept (its record's value) plus its slopes dotted with
/// the row's designated columns, standardized on the fly by the training
/// constants (missing values enter at the standardized mean, zero). slopes
/// holds numSlopes doubles per leaf in pre-order; leafOffset counts the
/// leaves consumed before this subtree. Returns the number of flattened
/// nodes consumed.
template <typename Columns>
inline size_t addFlatLinearPredictionsBelow(
    const FlatNode* flatNodes, const Columns& x,
    size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    size_t numSlopes, const double* slopes, size_t leafOffset = 0,
    const std::uint64_t* maskWords = nullptr) {
  if (flatNodes[0].variable == invalidVariable) {
    const double* leafSlopes = slopes + leafOffset * numSlopes;
    for (size_t k = lo; k < hi; ++k) {
      size_t row = indices[k];
      double fit = flatNodes[0].value;
      for (size_t j = 0; j < numSlopes; ++j) {
        double value = x.column(columns[j]).at(row);
        fit += leafSlopes[j] *
               (isNA(value) ? 0.0 : (value - means[j]) / sds[j]);
      }
      fits[row] += fit;
    }
    return 1;
  }

  size_t mid =
    partitionFlatIndices(flatNodes[0], x, indices, lo, hi, maskWords);
  size_t numOnLeft = addFlatLinearPredictionsBelow(
    flatNodes + 1, x, indices, lo, mid, fits, columns, means,
    sds, numSlopes, slopes, leafOffset, maskWords);
  size_t numNodes = 1 + numOnLeft;
  numNodes += addFlatLinearPredictionsBelow(
    flatNodes + numNodes, x, indices, mid, hi, fits, columns,
    means, sds, numSlopes, slopes, leafOffset + (numOnLeft + 1) / 2,
    maskWords);
  return numNodes;
}

/// The raw column-major entry: numRows is the block's row stride.
inline size_t addFlatLinearPredictionsBelow(
    const FlatNode* flatNodes, const double* x,
    size_t numRows, size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    size_t numSlopes, const double* slopes, size_t leafOffset = 0,
    const std::uint64_t* maskWords = nullptr) {
  return addFlatLinearPredictionsBelow(
    flatNodes, DenseColumns{x, numRows}, indices, lo, hi, fits, columns, means,
    sds, numSlopes, slopes, leafOffset, maskWords);
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
template <typename Columns>
inline size_t addFlatFunctionPredictionsBelow(
    const FlatNode* flatNodes, const Columns& x,
    size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    const double* lengthscales, size_t numCovariates, const double* blocks,
    const size_t* blockOffsets, size_t leafOffset = 0,
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
        double value = x.column(columns[j]).at(row);
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
    partitionFlatIndices(flatNodes[0], x, indices, lo, hi, maskWords);
  size_t numOnLeft = addFlatFunctionPredictionsBelow(
    flatNodes + 1, x, indices, lo, mid, fits, columns, means,
    sds, lengthscales, numCovariates, blocks, blockOffsets, leafOffset,
    maskWords);
  size_t numNodes = 1 + numOnLeft;
  numNodes += addFlatFunctionPredictionsBelow(
    flatNodes + numNodes, x, indices, mid, hi, fits, columns,
    means, sds, lengthscales, numCovariates, blocks, blockOffsets,
    leafOffset + (numOnLeft + 1) / 2, maskWords);
  return numNodes;
}

/// The raw column-major entry: numRows is the block's row stride.
inline size_t addFlatFunctionPredictionsBelow(
    const FlatNode* flatNodes, const double* x,
    size_t numRows, size_t* indices, size_t lo, size_t hi, double* fits,
    const size_t* columns, const double* means, const double* sds,
    const double* lengthscales, size_t numCovariates, const double* blocks,
    const size_t* blockOffsets, size_t leafOffset = 0,
    const std::uint64_t* maskWords = nullptr) {
  return addFlatFunctionPredictionsBelow(
    flatNodes, DenseColumns{x, numRows}, indices, lo, hi, fits, columns, means,
    sds, lengthscales, numCovariates, blocks, blockOffsets, leafOffset,
    maskWords);
}

/// Structural well-formedness of a flattened subtree - complete pre-order,
/// variables in range, each node's tag matching its column, categorical
/// masks nonzero and confined to the categories, pooled masks' offsets
/// pre-order sequential within their channel - without the cut-correspondence
/// and gauge conditions live restoration demands (saved trees replay against
/// raw values, so any split value routes). Returns the number of nodes
/// consumed, 0 when malformed.
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
      (flat.flags & ~(flatMissingGoesRight | flatKindMask)) != 0)
    return 0;
  size_t variable = static_cast<size_t>(flat.variable);
  FlatKind kind = flatKindOf(flat);
  if (kind != expectedFlatKind(data, variable)) return 0;
  if (kind == FlatKind::categoricalPooled) {
    size_t numWords = maskWordsForCount(data.categoryCounts[variable]);
    if (masks == nullptr || maskCursor == nullptr ||
        flat.maskOffset != *maskCursor ||
        *maskCursor + numWords > numMaskWords)
      return 0;
    // the mask plus the missing direction must send something right
    if (maskIsZero(masks + *maskCursor, numWords) &&
        (flat.flags & flatMissingGoesRight) == 0)
      return 0;
    *maskCursor += numWords;
  } else if (kind == FlatKind::categoricalInline) {
    // an inline mask with no bit past the categories; the mask plus the
    // missing direction must send something right
    if ((flat.mask >> data.categoryCounts[variable]) != 0 ||
        (flat.mask == 0 && (flat.flags & flatMissingGoesRight) == 0))
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

/// Info dump of a flattened (saved) tree; the output format is R-visible and
/// pinned by tests: no occupancy or availability, ordinal splits by value,
/// leaf predictions on the internal scale. slopes, when non-null, holds
/// numSlopes per leaf in pre-order (vector-parameter leaves); leafOffset
/// counts the leaves consumed before this subtree.
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
    if (data.splitsBySubset(static_cast<size_t>(flat.variable))) {
      ext_printf("CATRule: ");
      std::uint32_t numCats =
        data.categoryCounts[static_cast<size_t>(flat.variable)];
      if (flatKindOf(flat) == FlatKind::categoricalPooled) {
        const std::uint64_t* directions = masks + flat.maskOffset;
        for (std::uint32_t i = 0; i < numCats; ++i)
          ext_printf(" %u", maskTestBit(directions, i) ? 1u : 0u);
      } else {
        for (std::uint32_t i = 0; i < numCats; ++i)
          ext_printf(" %u", static_cast<unsigned int>((flat.mask >> i) & 1));
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
