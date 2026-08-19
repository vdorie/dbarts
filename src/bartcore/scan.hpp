#ifndef BARTCORE_SCAN_HPP
#define BARTCORE_SCAN_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "data.hpp"
#include "model.hpp"

// Leaf-model-templated cut-scan primitive: for one ordinal variable over a
// node's index segment, histogram the members' sufficient statistics per code,
// prefix-scan the histogram into every candidate cut's collapsed left/right
// statistics, and score each cut by the leaf's integrated log-likelihood. The
// shared building block of grow-from-root (grow.hpp) and the informed-proposal
// birth/death of docs/design/parallel-bart-frontier.md; landed once here
// so both consumers include it unchanged (docs/design/grow-from-root.md).
//
// Occupancy-aware: a cut with a zero-count side gets the never-selected
// sentinel, so a scan-based builder never creates an empty leaf and the MH
// empty-leaf veto (logLikelihoodForBranch's branch rank) never fires on this
// path.
// The scan omits sum wz^2: it is dead weight for the constant leaf (additive
// over any partition of a node's fixed member set, so it cancels in every
// within-node comparison, cut-vs-cut and split-vs-no-split), so the bin and
// the marginal both carry only (count, sum w, sum wz).

namespace bartcore {

/// The never-selected weight of an occupancy-empty cut: -inf, so exp of it is
/// exactly zero after the caller's max-shift, and the cut is never drawn.
inline constexpr double cutScanEmptySentinel =
  -std::numeric_limits<double>::infinity();

/// Constant-leaf histogram bin: the (count, sum w, sum wz) reduction the
/// ConstantGaussianLeaf marginal consumes. count is the member tally, the
/// histogram's own census; the occupancy contract reads sumWeights, since the
/// emptiness law the move kernels hold counts positive WEIGHT and a bin of
/// zero-weight members carries positive count with none. A matrix-valued leaf
/// (linear, GP) replaces this triple with its (U'WU, U'Wz) block without
/// touching the scan's control flow.
struct ConstantLeafScanBin {
  double count = 0.0;
  double sumWeights = 0.0;
  double sumWeightedResponse = 0.0;

  void addObservation(double weight, double response) {
    count += 1.0;
    sumWeights += weight;
    sumWeightedResponse += weight * response;
  }
  void addBin(const ConstantLeafScanBin& other) {
    count += other.count;
    sumWeights += other.sumWeights;
    sumWeightedResponse += other.sumWeightedResponse;
  }
};

/// Scan every candidate cut of ordinal variable `variable` over the member ids
/// in indices[0, numMembers). Writes numCuts[variable] entries into
/// logLikelihood: entry k scores the split "codes <= k go left" of the node's
/// NON-MISSING members as leaf.logIntegratedLikelihood(left) +
/// logIntegratedLikelihood(right), or cutScanEmptySentinel when either side
/// carries no positive weight. Missing values (naCode) are excluded from the
/// split scan; grow-from-root routes them by the birth-time missing-direction
/// coin, and occupancy on the non-missing weights alone keeps both children
/// weighted, since a routed member only adds weight to the side it joins.
///
/// binScratch is caller-owned reused storage (numCuts[variable] + 1 bins).
/// The marginal carries no sum wz^2: that per-node total is identical under
/// every cut and under no-split, so it cancels against the no-split term the
/// caller assembles and is never accumulated.
///
/// RNG-free and residual-only: the scan reads y (the current tree residual)
/// and weights exactly as computeLeafStats does, so a builder's cut scores are
/// consistent with the leaf stats tree.birth later caches.
template <ScalarLeafModel L, typename ResidT = double>
void scanOrdinalCuts(const ColumnStore& data, std::size_t variable,
                     const index_t* indices, std::size_t numMembers,
                     const ResidT* y, const double* weights, const L& leaf,
                     double k, double residualVariance,
                     std::vector<ConstantLeafScanBin>& binScratch,
                     double* logLikelihood) {
  std::size_t numCuts = static_cast<std::size_t>(data.numCuts[variable]);
  std::size_t numBins = numCuts + 1;
  binScratch.assign(numBins, ConstantLeafScanBin{});

  // histogram the members' statistics per code, branching on storage once so
  // the dense inner loop stays a plain gather (the common case)
  bool dense = !data.columnIsSparse(variable);
  const xint_t* denseColumn = dense ? data.column(variable) : nullptr;
  const SparseColumnData* sparseColumn =
    dense ? nullptr : &data.sparseColumn(variable);
  for (std::size_t i = 0; i < numMembers; ++i) {
    std::size_t obs = indices[i];
    xint_t code = dense ? denseColumn[obs] : sparseColumn->at(obs);
    if (code == naCode) continue;  // routed by the birth-time coin, not scanned
    double weight = weights == nullptr ? 1.0 : weights[obs];
    binScratch[code].addObservation(weight, y[obs]);
  }

  // node total over the non-missing bins (the left-fold order the prefix scan
  // reproduces, so a cut's right side is bitwise total - left)
  ConstantLeafScanBin total;
  for (std::size_t b = 0; b < numBins; ++b) total.addBin(binScratch[b]);

  ConstantLeafScanBin left;
  for (std::size_t cut = 0; cut + 1 < numBins; ++cut) {
    left.addBin(binScratch[cut]);  // codes 0..cut go left
    double rightWeights = total.sumWeights - left.sumWeights;
    // Occupancy is the MOVES' emptiness law - positive weight on each side,
    // not positive member count (docs/design/empty-leaf-veto.md) - so a cut
    // isolating none but zero-weight members is undrawable rather than a leaf
    // a later birth would have to compare -inf against -inf out of. It
    // subsumes the member count, weights being nonnegative: a side with no
    // member sums to zero. Off an installed weight vector every member counts
    // 1.0 and the two laws are the same numbers, so the candidate set is
    // unchanged there; the subtraction is exact where it decides, a suffix of
    // exactly-zero bins leaving the total untouched.
    if (left.sumWeights <= 0.0 || rightWeights <= 0.0) {
      logLikelihood[cut] = cutScanEmptySentinel;  // never selected
      continue;
    }
    double rightWeightedResponse =
      total.sumWeightedResponse - left.sumWeightedResponse;
    logLikelihood[cut] =
      leaf.logIntegratedLikelihood(k, residualVariance, left.sumWeights,
                                   left.sumWeightedResponse) +
      leaf.logIntegratedLikelihood(k, residualVariance, rightWeights,
                                   rightWeightedResponse);
  }
}

// ---------------------------------------------------------------------------
// The categorical sibling of scanOrdinalCuts.
//
// INIT-ONLY. The candidate set below is not a valid Metropolis-Hastings
// neighborhood and must not be reused as one: above the cap it is truncated to
// a sorted-prefix family and REWEIGHTED so the family carries the whole
// variable's prior mass, and even below it the enumeration is over the
// partitions of the categories PRESENT at the node, with the absent reachable
// positions filled in by post-draw coins - which changes the reverse move's
// reachable sets. The grown forest is handed to the exact MH sweeps, which own
// stationarity; a proposal kernel needs the full 2^R - 2 mask space
// (CGMTreePrior::drawRuleForVariable) instead.

/// One category PRESENT at a node: its code and the (count, sum w, sum wz)
/// reduction over the members carrying it.
struct CategoricalScanEntry {
  std::uint32_t code = 0;
  ConstantLeafScanBin bin;
};

/// Caller-owned scratch for the categorical scan. `bins` and `seen` are
/// indexed by CATEGORY CODE and are never cleared wholesale: a column may
/// declare 65536 levels, so an O(K) clear at every node of every tree would
/// dominate the O(numMembers) histogram it precedes. First sight of a code
/// initializes its bin and pushes the code onto `touched`; the clear at the end
/// of the pass walks `touched` alone. Both arrays only ever grow, so the
/// zero-fill invariant on `seen` survives a variable with fewer categories.
struct CategoricalScanScratch {
  std::vector<ConstantLeafScanBin> bins;
  std::vector<std::uint8_t> seen;
  std::vector<std::uint32_t> touched;
  std::vector<CategoricalScanEntry> present;  // compact, sorted
};

/// Present-category count above which the exact enumeration gives way to the
/// sorted-prefix family. 2^(P-1) - 1 = 511 candidates already cover every
/// factor an R user constructs, and raising the cap relocates the boundary
/// rather than removing it.
inline constexpr std::size_t categoricalExhaustiveCap = 10;

/// How many candidates the enumeration emits for P present categories: the
/// full partition set below the cap, the sorted prefixes above it, and nothing
/// at all when the node holds fewer than two distinct categories (the variable
/// still counts in availableSplitProbability, exactly as an ordinal column all
/// of whose cuts are occupancy-empty already does).
inline std::size_t categoricalNumEmitted(std::size_t numPresent) {
  if (numPresent < 2) return 0;
  return numPresent <= categoricalExhaustiveCap
    ? (static_cast<std::size_t>(1) << (numPresent - 1)) - 1
    : numPresent - 1;
}

/// Exact branch: whether sorted present position `position` lies on the side
/// candidate `s` accumulates. The enumeration is a PLAIN BINARY COUNTER
/// s = 0 .. 2^(P-1) - 2 over the P - 1 non-anchor positions with position 0
/// (the anchor) always in, so every emitted candidate leaves both sides
/// non-empty, no partition repeats, and the full non-anchor mask - the one
/// that would empty a child - is outside the range by construction. Scoring
/// and the post-draw mask reconstruction call this same function, so a grown
/// mask cannot disagree with the candidate that was scored.
inline bool exactPartitionHoldsPosition(std::int32_t candidate,
                                        std::size_t position) {
  return position == 0 ||
         ((candidate >> (position - 1)) & 1) != 0;
}

/// Fisher branch: candidate c names the sorted prefix of length c + 1.
inline bool prefixPartitionHoldsPosition(std::int32_t candidate,
                                         std::size_t position) {
  return position <= static_cast<std::size_t>(candidate);
}

/// Pass 1 of the categorical scan, and the post-draw reconstruction: histogram
/// the members in indices[0, numMembers) by category code, compact the occupied
/// bins into scratch.present, and sort them ascending by the singleton-category
/// leaf posterior mean sum wz / (a s2 + sum w) - the shrunken analog of
/// LightGBM's cat_smooth, with a s2 > 0 always, so a zero-weight category needs
/// no special case - breaking ties on the code so the order is a deterministic
/// function of the inputs. Returns P, the number of distinct categories present.
///
/// Unlike scanOrdinalCuts there is NO naCode skip: a missing value takes the
/// column's reserved pseudo-category (63 inline, K pooled) and a categorical
/// rule routes it by another bit of the same mask, so it is a real category
/// here and a real histogram bin.
///
/// RNG-free and a pure function of the member set, which is what lets the
/// caller rebuild the winning variable's order AFTER the discrete draw rather
/// than carry it through. `leaf.scale` supplies the leaf prior's scale; every
/// constant leaf the grow path instantiates carries one.
template <ScalarLeafModel L, typename ResidT = double>
std::size_t scanCategoryHistogram(const ColumnStore& data, std::size_t variable,
                                  const index_t* indices,
                                  std::size_t numMembers, const ResidT* y,
                                  const double* weights, const L& leaf,
                                  double k, double residualVariance,
                                  CategoricalScanScratch& scratch) {
  std::size_t numBins =
    data.hasMissing[variable]
      ? static_cast<std::size_t>(missingCategoryCode(data.numCuts[variable])) + 1
      : static_cast<std::size_t>(data.numCuts[variable]);
  if (scratch.bins.size() < numBins) scratch.bins.resize(numBins);
  if (scratch.seen.size() < numBins) scratch.seen.resize(numBins, 0);
  scratch.touched.clear();

  bool dense = !data.columnIsSparse(variable);
  const xint_t* denseColumn = dense ? data.column(variable) : nullptr;
  const SparseColumnData* sparseColumn =
    dense ? nullptr : &data.sparseColumn(variable);
  for (std::size_t i = 0; i < numMembers; ++i) {
    std::size_t obs = indices[i];
    std::size_t code = dense ? denseColumn[obs] : sparseColumn->at(obs);
    if (scratch.seen[code] == 0) {
      scratch.seen[code] = 1;
      scratch.touched.push_back(static_cast<std::uint32_t>(code));
      scratch.bins[code] = ConstantLeafScanBin{};
    }
    scratch.bins[code].addObservation(weights == nullptr ? 1.0 : weights[obs],
                                      y[obs]);
  }

  scratch.present.clear();
  scratch.present.reserve(scratch.touched.size());
  for (std::uint32_t code : scratch.touched) {
    scratch.seen[code] = 0;  // the only clear the pass performs
    scratch.present.push_back(CategoricalScanEntry{code, scratch.bins[code]});
  }

  double priorVariance = (k / leaf.scale) * (k / leaf.scale) * residualVariance;
  std::sort(scratch.present.begin(), scratch.present.end(),
            [priorVariance](const CategoricalScanEntry& a,
                            const CategoricalScanEntry& b) {
              double keyA =
                a.bin.sumWeightedResponse / (priorVariance + a.bin.sumWeights);
              double keyB =
                b.bin.sumWeightedResponse / (priorVariance + b.bin.sumWeights);
              return keyA != keyB ? keyA < keyB : a.code < b.code;
            });
  return scratch.present.size();
}

/// Scan the partition candidates of categorical variable `variable` over the
/// member ids in indices[0, numMembers), resizing logLikelihood to the emitted
/// count (data-dependent, unlike an ordinal column's fixed numCuts) and
/// returning it. Entry c scores the partition candidate c decodes to as
/// leaf.logIntegratedLikelihood(left) + logIntegratedLikelihood(right), on
/// exactly the scale scanOrdinalCuts and the caller's no-split term use: the
/// omitted sum wz^2 is additive over ANY partition of the node's member set, so
/// it cancels against both.
///
/// Below the cap the enumeration is exact over the present partitions; above it
/// the sorted prefixes stand in for the family. The sentinel is one compare per
/// candidate on each side's WEIGHT, the moves' emptiness law: the enumeration
/// domain already rules out a member-empty side, so what the compare still
/// catches is a side of members that carry no weight, which scores identically
/// to an empty one (logIntegratedLikelihood returns 0.0 at sumWeights == 0)
/// and which a later birth out of would compare -inf against -inf.
template <ScalarLeafModel L, typename ResidT = double>
std::size_t scanCategoricalPartitions(const ColumnStore& data,
                                      std::size_t variable,
                                      const index_t* indices,
                                      std::size_t numMembers, const ResidT* y,
                                      const double* weights, const L& leaf,
                                      double k, double residualVariance,
                                      CategoricalScanScratch& scratch,
                                      std::vector<double>& logLikelihood) {
  std::size_t numPresent =
    scanCategoryHistogram(data, variable, indices, numMembers, y, weights, leaf,
                          k, residualVariance, scratch);
  std::size_t numEmitted = categoricalNumEmitted(numPresent);
  logLikelihood.resize(numEmitted);
  if (numEmitted == 0) return 0;

  const std::vector<CategoricalScanEntry>& present = scratch.present;
  ConstantLeafScanBin total;
  for (const CategoricalScanEntry& entry : present) total.addBin(entry.bin);

  auto score = [&](const ConstantLeafScanBin& left,
                   const ConstantLeafScanBin& right) {
    return left.sumWeights <= 0.0 || right.sumWeights <= 0.0
      ? cutScanEmptySentinel
      : leaf.logIntegratedLikelihood(k, residualVariance, left.sumWeights,
                                     left.sumWeightedResponse) +
        leaf.logIntegratedLikelihood(k, residualVariance, right.sumWeights,
                                     right.sumWeightedResponse);
  };

  if (numPresent <= categoricalExhaustiveCap) {
    // both sides accumulated from the compact array, so a candidate's suffstat
    // is a pure function of its index: no path-dependent add/subtract sequence,
    // no last-ulp drift into a negative sumWeights, and a trivial decode
    for (std::size_t c = 0; c < numEmitted; ++c) {
      ConstantLeafScanBin left, right;
      for (std::size_t position = 0; position < numPresent; ++position)
        (exactPartitionHoldsPosition(static_cast<std::int32_t>(c), position)
           ? left : right).addBin(present[position].bin);
      logLikelihood[c] = score(left, right);
    }
    return numEmitted;
  }

  // above the cap the prefixes nest, so they accumulate like the ordinal prefix
  // scan; recomputing each in O(P) would be quadratic in a 65536-level column
  ConstantLeafScanBin left;
  for (std::size_t c = 0; c + 1 < numPresent; ++c) {
    left.addBin(present[c].bin);
    ConstantLeafScanBin right;
    right.count = total.count - left.count;
    right.sumWeights = total.sumWeights - left.sumWeights;
    right.sumWeightedResponse =
      total.sumWeightedResponse - left.sumWeightedResponse;
    logLikelihood[c] = score(left, right);
  }
  return numEmitted;
}

}  // namespace bartcore

#endif  // BARTCORE_SCAN_HPP
