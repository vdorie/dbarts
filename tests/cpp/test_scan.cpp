// Component tests for the cut-scan primitive (scan.hpp): per-cut integrated
// log-likelihoods against an independent brute-force recompute (tolerance),
// against a from-scratch histogram collapse (bitwise), the promoted histogram-
// totals check, the occupancy sentinel, the routing of a node's missing rows
// into the child a candidate's direction names, the categorical scan's
// enumeration (exhaustive against a brute-forced partition set) and
// per-candidate scores, and the two scans' agreement on a shared partition.
// Needs only the data/model layers plus scan.hpp, so a touch to
// moves/chain/sampler does not recompile this TU.
#include "assert.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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

// The bin's split, pinned at a fixture where it bites: count tallies MEMBERS,
// sumWeights sums their weights, so a bin made entirely of zero-weight members
// carries positive count and zero sumWeights at once (ConstantLeafScanBin,
// scan.hpp). The occupancy gate reads sumWeights, the emptiness law the move
// kernels hold (docs/design/empty-leaf-veto.md), so a cut isolating such a bin
// is VETOED with the sentinel rather than scored - it would otherwise land a
// leaf whose every member is invisible to the likelihood, out of which a birth
// compares -inf against -inf.
void testCountWeightSplit() {
  // a saved snapshot of the shared runif01 stream keeps the seed-pinned suites
  // downstream of this TU bitwise intact
  uint64_t savedRngState = rngState;
  const size_t n = 400;
  std::vector<double> x(n), y(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i) / static_cast<double>(n - 1);  // full [0, 1]
    y[i] = runif01() - 0.5;
    w[i] = (i < n / 4) ? 0.0 : 1.0;  // the low quarter carries no weight
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 10);
  size_t numCuts = store.numCuts[0];

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  // x, hence code, is nondecreasing in i, so every code strictly below the
  // first positive-weight member's code is reached only by zero-weight
  // members - a bin made entirely of them
  xint_t boundaryCode = store.codeAt(0, n / 4);
  check(boundaryCode >= 1,
        "the fixture leaves at least one code entirely zero-weight");

  // direct histogram: every such bin carries positive count and zero
  // sumWeights - the split, isolated
  size_t numBins = numCuts + 1;
  std::vector<ConstantLeafScanBin> bins(numBins);
  for (size_t i = 0; i < n; ++i)
    bins[store.codeAt(0, i)].addObservation(w[i], y[i]);
  bool sawSplit = true;
  for (xint_t code = 0; code < boundaryCode; ++code)
    sawSplit &= bins[code].count > 0.0 && bins[code].sumWeights == 0.0;
  check(sawSplit, "every code below the boundary carries positive count and "
                  "zero sumWeights");

  // the scan over the same fixture: the cut isolating those codes carries no
  // weight on one side and is vetoed, while a cut with weight on both sides
  // still scores - the veto is the weight law, not a blanket refusal
  size_t isolatingCut = static_cast<size_t>(boundaryCode) - 1;
  ConstantGaussianLeaf leaf{0.5};
  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> scan(numCuts);
  scanOrdinalCuts(store, 0, members.data(), n, y.data(), w.data(), leaf, 2.0,
                  0.7, binScratch, scan.data());
  check(scan[isolatingCut] == cutScanEmptySentinel,
        "a zero-weight-only side is vetoed, member count notwithstanding");
  size_t weightedCut = static_cast<size_t>(store.codeAt(0, n / 2));
  check(weightedCut > static_cast<size_t>(boundaryCode) &&
          weightedCut + 1 < numCuts,
        "the fixture leaves a cut with weight on both sides");
  check(std::isfinite(scan[weightedCut]),
        "a cut carrying weight on both sides still scores");

  rngState = savedRngState;
  printf("ok: scan count/sumWeights split\n");
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

// Missing rows are ROUTED, not skipped: under MIA the rule's missing direction
// sends them to one child, so where a node holds them the scan emits both
// directions of every cut and each scores with those rows in the child it
// names. Checked against a brute force that routes them the same way, against a
// bitwise fresh collapse, and against the score the scan would have produced by
// dropping them - the defect this layout exists to close. A member set holding
// no missing row keeps the plain one-entry-per-cut layout.
void testMissingRouted() {
  const size_t n = 150;
  std::vector<double> x(n), y(n);
  size_t numMissing = 0;
  for (size_t i = 0; i < n; ++i) {
    bool isMissing = i % 7 == 0;
    if (isMissing) { x[i] = std::nan(""); ++numMissing; }
    else x[i] = runif01();
    // the missing rows carry a level of their own, so which child receives
    // them changes the score
    y[i] = runif01() - 0.5 + (isMissing ? 1.5 : 0.0);
  }
  ColumnStore store;
  store.build(x.data(), n, 1, 15);
  check(store.hasMissing[0] == 1, "column reports missing values");

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  ConstantGaussianLeaf leaf{0.5};
  double k = 2.0, residualVariance = 0.7;
  size_t numCuts = store.numCuts[0];
  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> scan(2 * numCuts), dropped;
  size_t numEmitted =
    scanOrdinalCuts(store, 0, members.data(), n, y.data(), nullptr, leaf, k,
                    residualVariance, binScratch, scan.data());
  check(numEmitted == 2 * numCuts,
        "a node holding missing rows emits both directions of every cut");

  // the missing rows' own reduction, and the per-cut non-missing sides
  ConstantLeafScanBin missing;
  std::vector<ConstantLeafScanBin> bins(numCuts + 1);
  for (size_t i = 0; i < n; ++i) {
    xint_t code = store.codeAt(0, i);
    if (code == naCode) missing.addObservation(1.0, y[i]);
    else bins[code].addObservation(1.0, y[i]);
  }
  check(missing.count == static_cast<double>(numMissing) &&
          missing.sumWeights > 0.0,
        "every missing member reduces into the routed bin");

  // an independent brute force over the raw members, routing naCode to the side
  // the entry names; and the same, dropping them, which the scan must NOT match
  bitwiseCollapse(store, 0, members.data(), n, y.data(), nullptr, leaf, k,
                  residualVariance, dropped);
  bool routedTol = true, differsFromDropped = true, sentinelPaired = true;
  for (size_t cut = 0; cut < numCuts; ++cut) {
    for (size_t side = 0; side < 2; ++side) {
      double lw = 0.0, lwz = 0.0, rw = 0.0, rwz = 0.0;
      size_t lc = 0, rc = 0;
      for (size_t i = 0; i < n; ++i) {
        xint_t code = store.codeAt(0, i);
        bool left = code == naCode ? side == 0
                                   : static_cast<size_t>(code) <= cut;
        if (left) { ++lc; lw += 1.0; lwz += y[i]; }
        else { ++rc; rw += 1.0; rwz += y[i]; }
      }
      double entry = scan[2 * cut + side];
      if (dropped[cut] == cutScanEmptySentinel) {
        sentinelPaired &= entry == cutScanEmptySentinel;
        continue;
      }
      sentinelPaired &= entry != cutScanEmptySentinel;
      check(lc > 0 && rc > 0, "a scored entry leaves both children occupied");
      double expected =
        leaf.logIntegratedLikelihood(k, residualVariance, lw, lwz) +
        leaf.logIntegratedLikelihood(k, residualVariance, rw, rwz);
      routedTol &= std::fabs(entry - expected) <=
                   1e-11 * (1.0 + std::fabs(expected));
      differsFromDropped &= std::fabs(entry - dropped[cut]) > 1e-6;
    }
  }
  check(routedTol, "every entry scores its cut with the missing rows in the "
                   "child its direction names");
  check(sentinelPaired,
        "the occupancy sentinel is on the non-missing sides, so it applies to "
        "both directions of a cut or to neither");
  check(differsFromDropped,
        "routing the missing rows moves every score off the dropped-rows one");

  // a member set with no missing row keeps the plain layout, bitwise
  std::vector<index_t> observed;
  for (size_t i = 0; i < n; ++i)
    if (store.codeAt(0, i) != naCode) observed.push_back(i);
  std::vector<double> plain(2 * numCuts), plainBitwise;
  size_t plainEmitted =
    scanOrdinalCuts(store, 0, observed.data(), observed.size(), y.data(),
                    nullptr, leaf, k, residualVariance, binScratch,
                    plain.data());
  bitwiseCollapse(store, 0, observed.data(), observed.size(), y.data(), nullptr,
                  leaf, k, residualVariance, plainBitwise);
  bool plainBit = plainEmitted == numCuts;
  for (size_t cut = 0; cut < numCuts; ++cut)
    plainBit &= plain[cut] == plainBitwise[cut];
  check(plainBit,
        "a node holding no missing row keeps one entry per cut (bitwise)");

  // A one-sided cut on a member set that DOES hold missing rows: occupancy is
  // read off the non-missing sides, so BOTH direction entries take the
  // sentinel. The pairing is what this asserts - an unwritten second entry
  // would leave a drawable candidate scoring a rule outside the ancestor
  // interval, which is how the scan subsumes that constraint. The buffer is
  // poisoned first so the assertion is about what the scan WRITES, not about
  // what a previous call happened to leave behind.
  std::vector<index_t> banded;
  for (size_t i = 0; i < n; ++i) {
    xint_t code = store.codeAt(0, i);
    if (code == naCode || (code >= 3 && code <= 6)) banded.push_back(i);
  }
  const double poison = -1.0;  // finite, and no marginal sum reaches it
  std::vector<double> bandedScan(2 * numCuts, poison);
  size_t bandedEmitted =
    scanOrdinalCuts(store, 0, banded.data(), banded.size(), y.data(), nullptr,
                    leaf, k, residualVariance, binScratch, bandedScan.data());
  check(bandedEmitted == 2 * numCuts,
        "the banded member set still holds missing rows");
  bool pairedSentinel = true, sawOneSided = false, sawScored = false;
  for (size_t cut = 0; cut < numCuts; ++cut) {
    size_t leftCount = 0, rightCount = 0;
    for (index_t m : banded) {
      xint_t code = store.codeAt(0, m);
      if (code == naCode) continue;  // routed, not scanned for occupancy
      if (static_cast<size_t>(code) <= cut) ++leftCount; else ++rightCount;
    }
    bool oneSided = leftCount == 0 || rightCount == 0;
    sawOneSided |= oneSided;
    sawScored |= !oneSided;
    for (size_t side = 0; side < 2; ++side)
      pairedSentinel &= (bandedScan[2 * cut + side] == cutScanEmptySentinel) ==
                        oneSided;
  }
  check(sawOneSided && sawScored,
        "the banded members leave both one-sided and splittable cuts");
  check(pairedSentinel,
        "a one-sided cut sentinels BOTH missing directions, and a splittable "
        "one neither");

  printf("ok: scan missing routed\n");
}

// The two scans agree where they can be made to answer the same question: one
// ordinal column and one categorical column carrying the SAME grouping and the
// SAME missing rows induce the same partitions of a node's members, so a scored
// ordinal entry and the categorical candidate over the same partition must
// carry the same score. They route missing rows by one rule - into the child
// the candidate names - and a divergence here is a defect in one of them.
void testOrdinalCategoricalScanAgreement() {
  uint64_t savedRngState = rngState;
  const size_t n = 96, numLevels = 4;
  std::vector<double> x(2 * n), y(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    double level = static_cast<double>(i % numLevels);
    x[i] = level;      // ordinal
    x[n + i] = level;  // categorical, same grouping
    y[i] = runif01() - 0.5 + 0.3 * level;
    w[i] = 0.5 + runif01();
  }
  for (size_t m = 0; m < n; m += 13) {  // the same rows missing in both
    x[m] = std::nan("");
    x[n + m] = std::nan("");
    y[m] += 1.2;  // and carrying a level of their own, so routing them matters
  }
  ColumnType types[] = {ColumnType::ordinal, ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 2, numLevels, false, types);
  check(store.hasMissing[0] == 1 && store.hasMissing[1] == 1,
        "both columns of the agreement fixture carry missing values");
  size_t numCuts = store.numCuts[0];
  check(numCuts >= numLevels - 1,
        "the ordinal column cuts between its four levels");

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  ConstantGaussianLeaf leaf{0.5};
  double k = 2.0, residualVariance = 0.7;
  std::vector<ConstantLeafScanBin> binScratch;
  std::vector<double> ordinalScan(2 * numCuts);
  size_t numEmitted =
    scanOrdinalCuts(store, 0, members.data(), n, y.data(), w.data(), leaf, k,
                    residualVariance, binScratch, ordinalScan.data());
  check(numEmitted == 2 * numCuts, "the ordinal scan routes its missing rows");

  CategoricalScanScratch scratch;
  std::vector<double> categoricalScan;
  size_t numCandidates =
    scanCategoricalPartitions(store, 1, members.data(), n, y.data(), w.data(),
                              leaf, k, residualVariance, scratch,
                              categoricalScan);
  size_t numPresent = scratch.present.size();
  check(numPresent == numLevels + 1 && numPresent <= categoricalExhaustiveCap,
        "four levels plus the missing pseudo-category, under the cap");

  // the left member set of every categorical candidate, as a row mask
  std::vector<std::vector<char>> categoricalLeft(numCandidates);
  for (size_t c = 0; c < numCandidates; ++c) {
    uint64_t leftCodes = 0;
    for (size_t position = 0; position < numPresent; ++position)
      if (exactPartitionHoldsPosition(static_cast<int32_t>(c), position))
        leftCodes |= 1ull << scratch.present[position].code;
    categoricalLeft[c].assign(n, 0);
    for (size_t i = 0; i < n; ++i)
      categoricalLeft[c][i] =
        ((leftCodes >> store.codeAt(1, i)) & 1ull) != 0 ? 1 : 0;
  }

  bool everyEntryMatched = true, scoresAgree = true;
  size_t matched = 0, scored = 0;
  for (size_t entry = 0; entry < numEmitted; ++entry) {
    if (ordinalScan[entry] == cutScanEmptySentinel) continue;
    ++scored;
    size_t cut = entry >> 1, side = entry & 1u;
    std::vector<char> left(n, 0);
    for (size_t i = 0; i < n; ++i) {
      xint_t code = store.codeAt(0, i);
      left[i] = (code == naCode ? side == 0
                                : static_cast<size_t>(code) <= cut) ? 1 : 0;
    }
    // the partition is unordered, so a candidate matches by either orientation
    bool found = false;
    for (size_t c = 0; c < numCandidates && !found; ++c) {
      bool same = true, complement = true;
      for (size_t i = 0; i < n; ++i) {
        same &= categoricalLeft[c][i] == left[i];
        complement &= categoricalLeft[c][i] != left[i];
      }
      if (!same && !complement) continue;
      found = true;
      ++matched;
      scoresAgree &=
        std::fabs(ordinalScan[entry] - categoricalScan[c]) <=
        1e-11 * (1.0 + std::fabs(categoricalScan[c]));
    }
    everyEntryMatched &= found;
  }
  check(scored > 0 && matched == scored && everyEntryMatched,
        "every scored ordinal entry names a partition the categorical scan "
        "also enumerates");
  check(scoresAgree,
        "the two scans score the same partition of the same members alike");

  rngState = savedRngState;
  printf("ok: scan ordinal/categorical agreement (%zu of %zu entries, %zu cuts, %zu candidates)\n", scored, numEmitted, numCuts, numCandidates);
}

// ---------------------------------------------------------------------------
// The categorical enumeration.
//
// The recipe is the plain binary counter s = 0 .. 2^(P-1) - 2 over the P - 1
// non-anchor positions. This checks it exhaustively against an independently
// brute-forced enumeration of the unordered both-sides-non-empty partitions,
// for P = 2..12 - a range that deliberately runs past the exhaustive cap, since
// the recipe is a combinatorial fact about the counter and does not stop being
// one where the scan stops using it. It goes red against the reflected Gray map
// G(g) = g ^ (g >> 1) over the same range, which emits the FULL non-anchor mask
// (an empty right child) for every P >= 3 and drops one legitimate partition
// while keeping the count right.
void testCategoricalEnumeration() {
  bool countExact = true, uniqueEverywhere = true, bothSidesNonEmpty = true;
  bool bijective = true;
  for (size_t numPresent = 2; numPresent <= 12; ++numPresent) {
    size_t numEmitted = (static_cast<size_t>(1) << (numPresent - 1)) - 1;

    // the emitted set, canonicalized as the side holding position 0
    std::vector<uint32_t> emitted;
    for (size_t c = 0; c < numEmitted; ++c) {
      uint32_t side = 0, other = 0;
      for (size_t position = 0; position < numPresent; ++position) {
        if (exactPartitionHoldsPosition(static_cast<int32_t>(c), position))
          side |= 1u << position;
        else
          other |= 1u << position;
      }
      bothSidesNonEmpty &= side != 0 && other != 0;
      bothSidesNonEmpty &= (side & 1u) != 0;  // the anchor is always held
      emitted.push_back(side);
    }
    std::vector<uint32_t> sorted(emitted);
    std::sort(sorted.begin(), sorted.end());
    uniqueEverywhere &=
      std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();

    // brute force: every subset of the P positions that holds position 0 and
    // leaves the complement non-empty is exactly one unordered partition
    std::vector<uint32_t> brute;
    for (uint32_t mask = 0; mask < (1u << numPresent); ++mask) {
      if ((mask & 1u) == 0) continue;
      if (mask == (1u << numPresent) - 1u) continue;
      brute.push_back(mask);
    }
    std::sort(brute.begin(), brute.end());
    countExact &= emitted.size() == brute.size();
    bijective &= sorted == brute;

    // the shipped emission count agrees below the cap and switches to the
    // sorted-prefix family above it
    size_t shipped = categoricalNumEmitted(numPresent);
    countExact &= shipped == (numPresent <= categoricalExhaustiveCap
                                ? numEmitted : numPresent - 1);
  }
  check(countExact, "the counter emits 2^(P-1) - 1 candidates, P = 2..12");
  check(bothSidesNonEmpty, "every emitted candidate leaves both sides occupied");
  check(uniqueEverywhere, "no induced partition is emitted twice");
  check(bijective,
        "the emitted set is the brute-forced partition set, P = 2..12");
  check(categoricalNumEmitted(0) == 0 && categoricalNumEmitted(1) == 0,
        "a node with fewer than two present categories emits nothing");

  // the Fisher branch's decode is the sorted prefix of length c + 1
  bool prefixExact = true;
  for (size_t c = 0; c + 1 < 12; ++c)
    for (size_t position = 0; position < 12; ++position)
      prefixExact &= prefixPartitionHoldsPosition(static_cast<int32_t>(c),
                                                  position) ==
                     (position <= c);
  check(prefixExact, "the prefix decode names the sorted prefix of length c + 1");

  printf("ok: scan categorical enumeration\n");
}

// The categorical scan itself, against an independent recompute over the
// partition each emitted candidate decodes to: both branches, the missing
// pseudo-category as a real bin, the shrunken sort order, and the P < 2 node
// that emits nothing while staying available.
void testCategoricalScan() {
  // a saved snapshot of the shared runif01 stream keeps the seed-pinned suites
  // downstream of this TU bitwise intact
  uint64_t savedRngState = rngState;
  const size_t n = 240;
  const uint32_t numLevels = 6;
  std::vector<double> x(n * 2), y(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % numLevels);
    // a second categorical column with 12 levels, past the exhaustive cap
    x[i + n] = static_cast<double>(i % 12);
    y[i] = runif01() - 0.5 + 0.4 * static_cast<double>(i % numLevels);
    w[i] = 0.5 + runif01();
  }
  x[3] = std::nan("");  // one missing value: a real category, not a skip
  ColumnType types[] = {ColumnType::categorical, ColumnType::categorical};
  ColumnStore store;
  store.build(x.data(), n, 2, 20, false, types);
  check(store.hasMissing[0] == 1 && store.numCuts[0] == numLevels,
        "the categorical fixture declares six levels and carries a missing one");

  std::vector<index_t> members(n);
  for (size_t i = 0; i < n; ++i) members[i] = i;

  ConstantGaussianLeaf leaf{0.5};
  double k = 2.0, residualVariance = 0.7;
  CategoricalScanScratch scratch;
  std::vector<double> scan;

  for (size_t variable = 0; variable < 2; ++variable) {
    size_t numEmitted =
      scanCategoricalPartitions(store, variable, members.data(), n, y.data(),
                                w.data(), leaf, k, residualVariance, scratch,
                                scan);
    size_t numPresent = scratch.present.size();
    bool exact = numPresent <= categoricalExhaustiveCap;
    check(numEmitted == categoricalNumEmitted(numPresent) && numEmitted > 0,
          variable == 0 ? "the exhaustive branch emits 2^(P-1) - 1 candidates"
                        : "the prefix branch emits P - 1 candidates");
    check(exact == (variable == 0),
          variable == 0 ? "six levels plus a missing one stay under the cap"
                        : "twelve levels run past the cap");

    // the histogram is a partition of the members, missing row included
    double totalCount = 0.0;
    for (const CategoricalScanEntry& entry : scratch.present)
      totalCount += entry.bin.count;
    check(totalCount == static_cast<double>(n),
          "every member lands in exactly one category bin");

    // the sort key is the singleton-category leaf posterior mean, ascending
    double priorVariance = (k / leaf.scale) * (k / leaf.scale) * residualVariance;
    bool ordered = true;
    for (size_t i = 0; i + 1 < numPresent; ++i) {
      double a = scratch.present[i].bin.sumWeightedResponse /
                 (priorVariance + scratch.present[i].bin.sumWeights);
      double b = scratch.present[i + 1].bin.sumWeightedResponse /
                 (priorVariance + scratch.present[i + 1].bin.sumWeights);
      ordered &= a < b || (a == b && scratch.present[i].code <
                                     scratch.present[i + 1].code);
    }
    check(ordered, "present categories sort by the shrunken category mean");

    // an independent recompute of every emitted candidate, over the raw
    // members rather than the compacted bins
    bool scored = true, anySentinel = false;
    std::vector<CategoricalScanEntry> present(scratch.present);
    for (size_t c = 0; c < numEmitted; ++c) {
      uint64_t rightCodes = 0;
      for (size_t position = 0; position < numPresent; ++position) {
        bool holds = exact
          ? exactPartitionHoldsPosition(static_cast<int32_t>(c), position)
          : prefixPartitionHoldsPosition(static_cast<int32_t>(c), position);
        if (!holds) rightCodes |= 1ull << present[position].code;
      }
      double lw = 0.0, lwz = 0.0, rw = 0.0, rwz = 0.0;
      size_t lc = 0, rc = 0;
      for (size_t i = 0; i < n; ++i) {
        xint_t code = store.codeAt(variable, i);
        bool right = ((rightCodes >> code) & 1ull) != 0;
        if (right) { ++rc; rw += w[i]; rwz += w[i] * y[i]; }
        else { ++lc; lw += w[i]; lwz += w[i] * y[i]; }
      }
      anySentinel |= scan[c] == cutScanEmptySentinel;
      double expected =
        leaf.logIntegratedLikelihood(k, residualVariance, lw, lwz) +
        leaf.logIntegratedLikelihood(k, residualVariance, rw, rwz);
      scored &= lc > 0 && rc > 0 &&
                std::fabs(scan[c] - expected) <=
                  1e-11 * (1.0 + std::fabs(expected));
    }
    check(scored, "every emitted candidate scores its own partition");
    check(!anySentinel,
          "the enumeration domain never reaches the occupancy sentinel");
  }

  // a node holding one category emits nothing and spends no candidate
  std::vector<index_t> single;
  for (size_t i = 0; i < n; ++i)
    if (store.codeAt(0, i) == 2) single.push_back(i);
  check(!single.empty(), "the single-category subset is populated");
  size_t emptyEmission =
    scanCategoricalPartitions(store, 0, single.data(), single.size(), y.data(),
                              nullptr, leaf, k, residualVariance, scratch, scan);
  check(emptyEmission == 0 && scratch.present.size() == 1,
        "a node with one present category emits no candidate");

  rngState = savedRngState;
  printf("ok: scan categorical partitions\n");
}

}  // namespace

void runScanTests() {
  testScanAgreement();
  testHistogramTotals();
  testCountWeightSplit();
  testOccupancySentinel();
  testMissingRouted();
  testCategoricalEnumeration();
  testCategoricalScan();
  testOrdinalCategoricalScanAgreement();
}
