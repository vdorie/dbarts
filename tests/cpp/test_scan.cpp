// Component tests for the cut-scan primitive (scan.hpp): per-cut integrated
// log-likelihoods against an independent brute-force recompute (tolerance),
// against a from-scratch histogram collapse (bitwise), the promoted histogram-
// totals check, and the occupancy sentinel. Needs only the data/model layers
// plus scan.hpp, so a touch to moves/chain/sampler does not recompile this TU.
#include "assert.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include <misc/stats.h>

#include <bartcore/data.hpp>
#include <bartcore/model.hpp>
#include <bartcore/scan.hpp>

using namespace bartcore;

namespace {

// Independent per-cut oracle: for each cut, partition the (non-missing) members
// in member order and score with the same marginal the scan calls. Member-order
// accumulation differs from the scan's per-code order, so this agrees only to a
// tolerance; the bitwise agreement is checked separately below.
void bruteForcePerCut(const ColumnStore& store, size_t variable,
                      const index_t* indices, size_t numMembers, const double* y,
                      const double* weights, const ConstantGaussianLeaf& leaf,
                      double k, double residualVariance,
                      std::vector<double>& out) {
  size_t numCuts = store.numCuts[variable];
  out.assign(numCuts, 0.0);
  for (size_t cut = 0; cut < numCuts; ++cut) {
    double lw = 0.0, lwz = 0.0, lwz2 = 0.0, rw = 0.0, rwz = 0.0, rwz2 = 0.0;
    size_t lc = 0, rc = 0;
    for (size_t i = 0; i < numMembers; ++i) {
      size_t obs = indices[i];
      xint_t code = store.codeAt(variable, obs);
      if (code == naCode) continue;
      double w = weights == nullptr ? 1.0 : weights[obs];
      double wz = w * y[obs];
      if (static_cast<size_t>(code) <= cut) {
        ++lc; lw += w; lwz += wz; lwz2 += wz * y[obs];
      } else {
        ++rc; rw += w; rwz += wz; rwz2 += wz * y[obs];
      }
    }
    if (lc == 0 || rc == 0) {
      out[cut] = cutScanEmptySentinel;
    } else {
      out[cut] =
        leaf.logIntegratedLikelihood(k, residualVariance, lw, lwz) +
        leaf.logIntegratedLikelihood(k, residualVariance, rw, rwz);
    }
  }
}

// From-scratch histogram (member order) then fresh per-cut collapse, the same
// left-fold order the scan uses, so this agrees BITWISE with the scan.
void bitwiseCollapse(const ColumnStore& store, size_t variable,
                     const index_t* indices, size_t numMembers, const double* y,
                     const double* weights, const ConstantGaussianLeaf& leaf,
                     double k, double residualVariance,
                     std::vector<double>& out) {
  size_t numCuts = store.numCuts[variable];
  size_t numBins = numCuts + 1;
  std::vector<ConstantLeafScanBin> bins(numBins);
  for (size_t i = 0; i < numMembers; ++i) {
    size_t obs = indices[i];
    xint_t code = store.codeAt(variable, obs);
    if (code == naCode) continue;
    double w = weights == nullptr ? 1.0 : weights[obs];
    bins[code].addObservation(w, y[obs]);
  }
  ConstantLeafScanBin total;
  for (size_t b = 0; b < numBins; ++b) total.addBin(bins[b]);
  out.assign(numCuts, 0.0);
  for (size_t cut = 0; cut + 1 < numBins; ++cut) {
    ConstantLeafScanBin left;
    for (size_t b = 0; b <= cut; ++b) left.addBin(bins[b]);
    if (left.count == 0.0 || left.count == total.count) {
      out[cut] = cutScanEmptySentinel;
      continue;
    }
    double rw = total.sumWeights - left.sumWeights;
    double rwz = total.sumWeightedResponse - left.sumWeightedResponse;
    out[cut] =
      leaf.logIntegratedLikelihood(k, residualVariance, left.sumWeights,
                                   left.sumWeightedResponse) +
      leaf.logIntegratedLikelihood(k, residualVariance, rw, rwz);
  }
}

void testScanAgreement() {
  const size_t n = 200;
  std::vector<double> x(n), y(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    y[i] = runif01() - 0.5;
    w[i] = 0.5 + runif01();
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 20);

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  ConstantGaussianLeaf leaf{0.5};
  double k = 2.0, residualVariance = 0.7;
  size_t numCuts = store.numCuts[0];

  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> scan(numCuts), brute, bitwise;

  // whole-node scan, unweighted and weighted
  for (int weighted = 0; weighted <= 1; ++weighted) {
    const double* weights = weighted ? w.data() : nullptr;
    scanOrdinalCuts(store, 0, members.data(), n, y.data(), weights, leaf, k,
                    residualVariance, binScratch, scan.data());
    bruteForcePerCut(store, 0, members.data(), n, y.data(), weights, leaf, k,
                     residualVariance, brute);
    bitwiseCollapse(store, 0, members.data(), n, y.data(), weights, leaf, k,
                    residualVariance, bitwise);
    bool okTol = true, okBit = true;
    for (size_t cut = 0; cut < numCuts; ++cut) {
      okBit &= scan[cut] == bitwise[cut];
      okTol &= std::fabs(scan[cut] - brute[cut]) <=
               1e-11 * (1.0 + std::fabs(brute[cut]));
    }
    check(okBit, weighted ? "scan == fresh collapse, weighted (bitwise)"
                          : "scan == fresh collapse, unweighted (bitwise)");
    check(okTol, weighted ? "scan == brute per-cut, weighted"
                          : "scan == brute per-cut, unweighted");
  }

  // scattered subset of members (a non-root leaf's gather access pattern)
  std::vector<index_t> subset;
  for (size_t i = 0; i < n; i += 3) subset.push_back((7 * i + 11) % n);
  scanOrdinalCuts(store, 0, subset.data(), subset.size(), y.data(), w.data(),
                  leaf, k, residualVariance, binScratch, scan.data());
  bruteForcePerCut(store, 0, subset.data(), subset.size(), y.data(), w.data(),
                   leaf, k, residualVariance, brute);
  bitwiseCollapse(store, 0, subset.data(), subset.size(), y.data(), w.data(),
                  leaf, k, residualVariance, bitwise);
  bool okTol = true, okBit = true;
  for (size_t cut = 0; cut < numCuts; ++cut) {
    okBit &= scan[cut] == bitwise[cut];
    okTol &= std::fabs(scan[cut] - brute[cut]) <=
             1e-11 * (1.0 + std::fabs(brute[cut]));
  }
  check(okBit, "scattered-subset scan == fresh collapse (bitwise)");
  check(okTol, "scattered-subset scan == brute per-cut");

  printf("ok: scan agreement\n");
}

// Promote the prototype's correctnessCheck: the scan's per-code histogram
// totals equal a direct sufficient-statistic over the same members.
void testHistogramTotals() {
  const size_t n = 512;
  std::vector<double> x(n), y(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = runif01();
    y[i] = runif01() - 0.5;
    w[i] = 0.5 + runif01();
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 100);

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  double swK, swzK;
  misc_computeIndexedWeightedSufficientStatisticsFast(
    y.data(), members.data(), n, w.data(), &swK, &swzK);

  size_t numBins = store.numCuts[0] + 1;
  std::vector<ConstantLeafScanBin> bins(numBins);
  for (size_t i = 0; i < n; ++i)
    bins[store.codeAt(0, i)].addObservation(w[i], y[i]);
  ConstantLeafScanBin total;
  for (size_t b = 0; b < numBins; ++b) total.addBin(bins[b]);

  checkNear(total.sumWeights, swK, 1e-9 * (1.0 + std::fabs(swK)),
            "histogram total sum w == direct suffstat");
  checkNear(total.sumWeightedResponse, swzK, 1e-9 * (1.0 + std::fabs(swzK)),
            "histogram total sum wz == direct suffstat");
  check(total.count == static_cast<double>(n), "histogram counts every member");

  printf("ok: scan histogram totals\n");
}

// Occupancy gate: a cut whose either side has zero non-missing members gets the
// never-selected sentinel, and its softmax weight is exactly zero. The cut grid
// adapts to the column's range, so one-sidedness is forced by scanning only the
// members that land in a middle code band of a full-range column.
void testOccupancySentinel() {
  const size_t n = 400;
  std::vector<double> x(n), y(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i) / static_cast<double>(n - 1);  // full [0, 1]
    y[i] = runif01() - 0.5;
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 10);
  size_t numCuts = store.numCuts[0];

  // gather only the members whose code is in the interior band [3, 6], so every
  // cut below 3 leaves the left empty and every cut at or above 6 the right
  std::vector<index_t> members;
  for (size_t i = 0; i < n; ++i) {
    xint_t code = store.codeAt(0, i);
    if (code >= 3 && code <= 6) members.push_back(i);
  }
  check(!members.empty(), "interior band is populated");

  ConstantGaussianLeaf leaf{0.5};
  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> scan(numCuts);
  scanOrdinalCuts(store, 0, members.data(), members.size(), y.data(), nullptr,
                  leaf, 2.0, 0.7, binScratch, scan.data());

  // an independent per-code census over the scanned members: a cut is one-sided
  // iff every scanned member falls on one side of it
  bool sawSentinel = false, sawFinite = false, occupancyExact = true;
  for (size_t cut = 0; cut < numCuts; ++cut) {
    size_t leftCount = 0, rightCount = 0;
    for (size_t m : members) {
      xint_t code = store.codeAt(0, m);
      if (static_cast<size_t>(code) <= cut) ++leftCount; else ++rightCount;
    }
    bool oneSided = leftCount == 0 || rightCount == 0;
    bool isSentinel = scan[cut] == cutScanEmptySentinel;
    occupancyExact &= (oneSided == isSentinel);
    sawSentinel |= isSentinel;
    sawFinite |= !isSentinel && std::isfinite(scan[cut]);
  }
  check(occupancyExact,
        "sentinel iff a side is empty (occupancy contract)");
  check(sawSentinel, "the banded members produce one-sided cuts");
  check(sawFinite, "the banded members produce splittable cuts");

  // a sentinel cut is never drawn: after the caller's max-shift its weight is
  // exactly zero
  double maxLog = -std::numeric_limits<double>::infinity();
  for (size_t cut = 0; cut < numCuts; ++cut)
    if (scan[cut] > maxLog) maxLog = scan[cut];
  bool sentinelWeightsAreZero = true;
  for (size_t cut = 0; cut < numCuts; ++cut)
    if (scan[cut] == cutScanEmptySentinel)
      sentinelWeightsAreZero &= std::exp(scan[cut] - maxLog) == 0.0;
  check(sentinelWeightsAreZero, "sentinel cut has exactly zero draw weight");

  printf("ok: scan occupancy sentinel\n");
}

// Missing values are excluded from the split scan and left for the birth-time
// coin; the scan over a column with NAs matches a brute force that also skips
// them, and the histogram counts only the non-missing members.
void testMissingExcluded() {
  const size_t n = 150;
  std::vector<double> x(n), y(n);
  size_t numMissing = 0;
  for (size_t i = 0; i < n; ++i) {
    if (i % 7 == 0) { x[i] = std::nan(""); ++numMissing; }
    else x[i] = runif01();
    y[i] = runif01() - 0.5;
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 15);
  check(store.hasMissing[0] == 1, "column reports missing values");

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  ConstantGaussianLeaf leaf{0.5};
  size_t numCuts = store.numCuts[0];
  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> scan(numCuts), bitwise;
  scanOrdinalCuts(store, 0, members.data(), n, y.data(), nullptr, leaf, 2.0, 0.7,
                  binScratch, scan.data());
  bitwiseCollapse(store, 0, members.data(), n, y.data(), nullptr, leaf, 2.0, 0.7,
                  bitwise);

  size_t numBins = numCuts + 1;
  std::vector<ConstantLeafScanBin> bins(numBins);
  for (size_t i = 0; i < n; ++i) {
    xint_t code = store.codeAt(0, i);
    if (code == naCode) continue;
    bins[code].addObservation(1.0, y[i]);
  }
  ConstantLeafScanBin total;
  for (size_t b = 0; b < numBins; ++b) total.addBin(bins[b]);

  bool okBit = true;
  for (size_t cut = 0; cut < numCuts; ++cut) okBit &= scan[cut] == bitwise[cut];
  check(okBit, "scan with missing == fresh collapse (bitwise)");
  check(total.count == static_cast<double>(n - numMissing),
        "histogram excludes missing members");

  printf("ok: scan missing excluded\n");
}

}  // namespace

void runScanTests() {
  testScanAgreement();
  testHistogramTotals();
  testOccupancySentinel();
  testMissingExcluded();
}
