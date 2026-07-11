#ifndef BARTCORE_SCAN_HPP
#define BARTCORE_SCAN_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "data.hpp"
#include "model.hpp"

// Leaf-model-templated cut-scan primitive: for one ordinal variable over a
// node's index segment, histogram the members' sufficient statistics per code,
// prefix-scan the histogram into every candidate cut's collapsed left/right
// statistics, and score each cut by the leaf's integrated log-likelihood. The
// shared building block of grow-from-root (grow.hpp) and the informed-proposal
// birth/death of docs/design/parallel-bart-frontier.md 3.1; landed once here
// so both consumers include it unchanged (docs/design/grow-from-root.md).
//
// Occupancy-aware: a cut with a zero-count side gets the never-selected
// sentinel, so a scan-based builder never creates an empty leaf and the MH
// empty-leaf veto (logLikelihoodForBranch's -1e7) never runs on this path.
// The scan omits sum wz^2: it is dead weight for the constant leaf (additive
// over any partition of a node's fixed member set, so it cancels in every
// within-node comparison, cut-vs-cut and split-vs-no-split), which is why the
// bin carries only (count, sum w, sum wz) and the marginal is evaluated with a
// zero sumWeightedResponseSq.

namespace bartcore {

/// The never-selected weight of an occupancy-empty cut: -inf, so exp of it is
/// exactly zero after the caller's max-shift, and the cut is never drawn.
inline constexpr double cutScanEmptySentinel =
  -std::numeric_limits<double>::infinity();

/// Constant-leaf histogram bin: the (count, sum w, sum wz) reduction the
/// ConstantGaussianLeaf marginal consumes. count is the member tally driving
/// the occupancy contract (never the weight, which a zero-weight member does
/// not zero). A matrix-valued leaf (linear, GP) replaces this triple with its
/// (U'WU, U'Wz) block without touching the scan's control flow.
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

/// Histogram one dense column's members per code at its storage width. Missing
/// codes (the width's naCode) are skipped, the same left-fold order the prefix
/// scan reproduces. One typed gather per instantiation, no per-element branch.
template <typename CodeT>
inline void histogramDenseCutScan(const CodeT* column,
                                  const std::size_t* indices,
                                  std::size_t numMembers, const double* y,
                                  const double* weights,
                                  std::vector<ConstantLeafScanBin>& bins) {
  for (std::size_t i = 0; i < numMembers; ++i) {
    std::size_t obs = indices[i];
    CodeT code = column[obs];
    if (code == CodeWidthTraits<CodeT>::na) continue;
    double weight = weights == nullptr ? 1.0 : weights[obs];
    bins[code].addObservation(weight, y[obs]);
  }
}

/// Scan every candidate cut of ordinal variable `variable` over the member ids
/// in indices[0, numMembers). Writes numCuts[variable] entries into
/// logLikelihood: entry k scores the split "codes <= k go left" of the node's
/// NON-MISSING members as leaf.logIntegratedLikelihood(left) +
/// logIntegratedLikelihood(right), or cutScanEmptySentinel when either side
/// has zero members. Missing values (naCode) are excluded from the split scan;
/// grow-from-root routes them by the birth-time missing-direction coin, and
/// occupancy on the non-missing counts alone keeps both children non-empty.
///
/// binScratch is caller-owned reused storage (numCuts[variable] + 1 bins).
/// sumWeightedResponseSq is passed as zero to the marginal: its per-node total
/// is identical under every cut and under no-split, so it cancels against the
/// no-split term the caller assembles and need never be accumulated.
///
/// RNG-free and residual-only: the scan reads y (the current tree residual)
/// and weights exactly as computeLeafStats does, so a builder's cut scores are
/// consistent with the leaf stats tree.birth later caches.
template <ScalarLeafModel L>
void scanOrdinalCuts(const ColumnStore& data, std::size_t variable,
                     const std::size_t* indices, std::size_t numMembers,
                     const double* y, const double* weights, const L& leaf,
                     double k, double residualVariance,
                     std::vector<ConstantLeafScanBin>& binScratch,
                     double* logLikelihood) {
  std::size_t numCuts = static_cast<std::size_t>(data.numCuts[variable]);
  std::size_t numBins = numCuts + 1;
  binScratch.assign(numBins, ConstantLeafScanBin{});

  // histogram the members' statistics per code, branching on storage and width
  // once so each dense inner loop stays a plain typed gather (the common case)
  if (data.columnIsSparse(variable)) {
    const SparseColumnData& sparseColumn = data.sparseColumn(variable);
    for (std::size_t i = 0; i < numMembers; ++i) {
      std::size_t obs = indices[i];
      xint_t code = sparseColumn.at(obs);
      if (code == naCode) continue;  // routed by the birth-time coin, not scanned
      double weight = weights == nullptr ? 1.0 : weights[obs];
      binScratch[code].addObservation(weight, y[obs]);
    }
  } else if (data.codeWidth[variable] == 1) {
    histogramDenseCutScan(data.columnU8(variable), indices, numMembers, y,
                          weights, binScratch);
  } else {
    histogramDenseCutScan(data.columnU16(variable), indices, numMembers, y,
                          weights, binScratch);
  }

  // node total over the non-missing bins (the left-fold order the prefix scan
  // reproduces, so a cut's right side is bitwise total - left)
  ConstantLeafScanBin total;
  for (std::size_t b = 0; b < numBins; ++b) total.addBin(binScratch[b]);

  ConstantLeafScanBin left;
  for (std::size_t cut = 0; cut + 1 < numBins; ++cut) {
    left.addBin(binScratch[cut]);  // codes 0..cut go left
    if (left.count == 0.0 || left.count == total.count) {
      logLikelihood[cut] = cutScanEmptySentinel;  // occupancy: never selected
      continue;
    }
    double rightWeights = total.sumWeights - left.sumWeights;
    double rightWeightedResponse =
      total.sumWeightedResponse - left.sumWeightedResponse;
    logLikelihood[cut] =
      leaf.logIntegratedLikelihood(k, residualVariance, left.sumWeights,
                                   left.sumWeightedResponse, 0.0) +
      leaf.logIntegratedLikelihood(k, residualVariance, rightWeights,
                                   rightWeightedResponse, 0.0);
  }
}

}  // namespace bartcore

#endif  // BARTCORE_SCAN_HPP
