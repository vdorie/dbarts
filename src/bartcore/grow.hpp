#ifndef BARTCORE_GROW_HPP
#define BARTCORE_GROW_HPP

#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <external/random.h>

#include "data.hpp"
#include "model.hpp"
#include "scan.hpp"
#include "tree.hpp"

// XBART-style root-down stochastic tree construction (He, Yalov and Hahn 2019),
// the warm-start forest builder of docs/design/grow-from-root.md. Parallel to
// moves.hpp's metropolisJumpForTree, but a builder rather than a reversible MH
// kernel: at each node it scans every ordinal cut and every categorical
// partition candidate of every available variable (scan.hpp) and draws one
// outcome from the weighted set {no-split} U {(var, candidate)}, the CGM tree
// prior's own factors times the integrated leaf likelihood.
// It is NOT a posterior kernel; the exact MH sweeps that follow own
// stationarity, so the grown forest need only be a LEGAL chain state (non-empty
// leaves, rules in gauge, availability/depth vetoes respected).

namespace bartcore {

/// Per-node reusable scratch for growTreeFromRoot; one instance per chain,
/// reused down the whole recursion (a node's assembly is over before its
/// children recurse, so nothing aliases).
struct GrowScratch {
  std::vector<ConstantLeafScanBin> binScratch;  // scan.hpp ordinal histogram
  CategoricalScanScratch categoryScan;          // scan.hpp category histogram
  std::vector<std::uint8_t> available;          // per-predictor availability
  std::vector<double> cutLogLikelihood;         // per-variable candidate scores
  std::vector<double> candidateLogWeights;      // {no-split} U {valid (var,cut)}
  std::vector<double> candidateWeights;         // exp'd, normalized for the draw
  std::vector<std::int32_t> candidateVariable;  // parallel to the weights
  std::vector<std::int32_t> candidateCut;
  // parallel too: the missing direction an ordinal candidate already names
  // (0 left, 1 right), or -1 where nothing routed and a coin decides
  std::vector<std::int8_t> candidateMissing;
  // pooled categorical masks; grow owns these rather than borrowing
  // Tree::reachableScratch_, which variableAvailable and buildFromFlatBelow
  // also use
  std::vector<std::uint64_t> reachableMask;
  std::vector<std::uint64_t> presentMask;
};

/// Turn the winning categorical candidate into a rule, in the fixed order the
/// draw discipline documents: rebuild the variable's present categories and
/// their sort order (deterministic and RNG-free, so nothing has to survive the
/// discrete draw), decode the candidate into a side of the present set, flip a
/// single orientation coin, then walk the REACHABLE positions ascending and
/// flip one coin per ABSENT one. `rule.variableIndex` is already set; this
/// fills in the mask, inline or through the tree's pool.
template <ScalarLeafModel L, typename ResidT = double>
void growCategoricalRule(const ColumnStore& data, const L& leaf, ext_rng* rng,
                         Tree& tree, std::int32_t nodeIndex, const ResidT* y,
                         const double* weights, double k,
                         double residualVariance, std::size_t begin,
                         std::size_t numMembers, std::int32_t candidate,
                         GrowScratch& scratch, Rule& rule) {
  std::size_t j = static_cast<std::size_t>(rule.variableIndex);
  std::size_t numPresent =
    scanCategoryHistogram(data, j, tree.indices + begin, numMembers, y, weights,
                          leaf, k, residualVariance, scratch.categoryScan);
  const std::vector<CategoricalScanEntry>& present =
    scratch.categoryScan.present;
  bool exact = numPresent <= categoricalExhaustiveCap;
  bool orientation = ext_rng_simulateBernoulli(rng, 0.5) == 1;
  auto goesRight = [&](std::size_t position) {
    bool holds = exact ? exactPartitionHoldsPosition(candidate, position)
                       : prefixPartitionHoldsPosition(candidate, position);
    return holds != orientation;  // the complement when the coin says so
  };

  if (data.columnIsPooled(j)) {
    std::size_t numWords = maskWordsForCount(data.numCuts[j]);
    scratch.reachableMask.resize(numWords);
    tree.reachableCategoriesWide(data, nodeIndex,
                                 static_cast<std::int32_t>(j),
                                 scratch.reachableMask.data());
    scratch.presentMask.assign(numWords, 0);
    // grow never rejects, so the pool needs no mark and no truncate;
    // growForestFromRoot compacts it per grown tree
    std::size_t offset = tree.allocateMask(numWords);
    std::uint64_t* directions = tree.mutableMaskWordsFor(offset);
    for (std::size_t position = 0; position < numPresent; ++position) {
      maskSetBit(scratch.presentMask.data(), present[position].code);
      if (goesRight(position)) maskSetBit(directions, present[position].code);
    }
    for (std::size_t w = 0; w < numWords; ++w) {
      std::uint64_t absent =
        scratch.reachableMask[w] & ~scratch.presentMask[w];
      while (absent != 0) {
        std::uint32_t bit = static_cast<std::uint32_t>(std::countr_zero(absent));
        absent &= absent - 1;
        if (ext_rng_simulateBernoulli(rng, 0.5) == 1)
          maskSetBit(directions, bit + 64u * static_cast<std::uint32_t>(w));
      }
    }
    rule.setMaskOffset(offset);
    return;
  }

  std::uint64_t presentBits = 0, directions = 0;
  for (std::size_t position = 0; position < numPresent; ++position) {
    presentBits |= 1ull << present[position].code;
    if (goesRight(position)) directions |= 1ull << present[position].code;
  }
  std::uint64_t absent = tree.inlineReachableMasks()[j] & ~presentBits;
  while (absent != 0) {
    std::uint32_t category = static_cast<std::uint32_t>(std::countr_zero(absent));
    absent &= absent - 1;
    if (ext_rng_simulateBernoulli(rng, 0.5) == 1) directions |= 1ull << category;
  }
  rule.setCategoryDirections(directions);
}

/// Recursively split `tree` from nodeIndex against the residual y (per-tree
/// working response) and case weights. The leaf's integrated likelihood scores
/// each candidate cut through the scan; the discrete outcome is drawn on rng.
///
/// Draw discipline (documented so the equivalence picture stays predictable):
///  - A node whose prior growth probability is zero - the depth veto or no
///    available variable - draws NOTHING and stays a leaf.
///  - Every other node draws exactly ONE discrete outcome over {no-split} plus
///    the occupancy-nonempty ordinal cuts and the categorical partition
///    candidates of every available variable. A no-split outcome (always
///    representable) ends the node.
///  - An ORDINAL candidate on a column that routes missing values covers those
///    rows too, since under MIA the rule's missing direction sends them to one
///    child. Where the node HOLDS missing members the direction is therefore
///    part of the candidate: the scan emits both directions of every cut, each
///    scored with the missing rows in the child it names, each carrying CGM's
///    mass for ONE rule (ruleForVariableLogProbability's, the interval halved
///    by the two directions), and NO coin is spent. Where the node holds none,
///    the direction moves no statistic: the two rules score identically, the
///    candidate carries their combined mass whole, and ONE symmetric coin picks
///    within the pair after the draw, matching
///    CGMTreePrior::drawRuleForVariable. That is the same rule the categorical
///    branch follows, where a present missing pseudo-category is a real
///    histogram bin the partition routes and an absent one gets a coin.
///    test_grow.cpp chi-squares the realized root-rule frequencies against the
///    exact law over the full rule set and against the law the scan realized
///    while it dropped the missing rows from a split likelihood;
///    docs/design/grow-from-root.md section 7 records what it measured.
///  - A CATEGORICAL split draws exactly 1 + A more symmetric coins and NEVER
///    rejects: one orientation coin, then one per category REACHABLE at the
///    node but ABSENT from its members (A = R - P), taken over the reachable
///    positions ascending, the bit-by-bit convention and generator-granularity
///    reason of CGMTreePrior::drawCategoryPattern. There is no rejection loop
///    because the present side already leaves both children occupied, and the
///    absent positions are drawn rather than pinned left because pinning would
///    bias where a category unseen at the node routes at PREDICT time, in a
///    direction the MH sweeps do not correct. The candidate carries its whole
///    enumerable family's mass spread over what the enumeration emits
///    (CGMTreePrior::categoricalGroupLogProbability), which is the same
///    "candidate stands for a rule group" reading the ordinal missing
///    convention above is measured against.
///
/// Occupancy (scan.hpp) subsumes the ancestor split interval - members
/// violating an ancestor cut have codes confined to the interval, so
/// out-of-interval cuts have an empty side and are already zeroed - and keeps
/// both children carrying positive WEIGHT against the vector passed here, so
/// the built tree satisfies the same structural invariant the MH moves veto on
/// (docs/design/empty-leaf-veto.md). That vector is the caller's composed one,
/// per-forest under a coupling, which is why the law can be read off the scan's
/// own sums rather than off member counts. The categorical branch's enumeration
/// domain rules out an empty side by construction, and its mask comes out
/// nonzero, a strict subset of the node's reachable set and unequal to it,
/// which is buildFromFlat's gauge.
template <ScalarLeafModel L, typename ResidT = double>
void growTreeFromRoot(const ColumnStore& data, const CGMTreePrior& treePrior,
                      const L& leaf, ext_rng* rng, Tree& tree,
                      std::int32_t nodeIndex, const ResidT* y,
                      const double* weights, double k, double sigma,
                      GrowScratch& scratch) {
  double growth = treePrior.growthProbability(tree, data, nodeIndex);
  if (growth <= 0.0) return;  // depth / availability veto: a leaf, no draw

  double residualVariance = sigma * sigma;
  std::size_t begin = tree.at(nodeIndex).begin;
  std::size_t numMembers = tree.at(nodeIndex).numObservations();
  double nodeLogLikelihood =
    leaf.logIntegratedLikelihood(k, residualVariance,
                                 tree.at(nodeIndex).sumWeights,
                                 tree.at(nodeIndex).sumWeightedResponse);

  scratch.candidateLogWeights.clear();
  scratch.candidateVariable.clear();
  scratch.candidateCut.clear();
  scratch.candidateMissing.clear();
  // candidate 0 is always no-split: (1 - growth) * L(node), always finite
  scratch.candidateLogWeights.push_back(std::log(1.0 - growth) +
                                        nodeLogLikelihood);
  scratch.candidateVariable.push_back(invalidVariable);
  scratch.candidateCut.push_back(-1);
  scratch.candidateMissing.push_back(-1);

  scratch.available.resize(data.numPredictors);
  std::size_t numAvailable =
    tree.collectAvailableVariables(data, nodeIndex, scratch.available.data());

  double logGrowth = std::log(growth);
  const double* splitProbabilities = treePrior.splitProbabilities;
  double availableSplitProbability = 0.0;
  if (splitProbabilities != nullptr)
    for (std::size_t j = 0; j < data.numPredictors; ++j)
      if (scratch.available[j]) availableSplitProbability += splitProbabilities[j];

  for (std::size_t j = 0; j < data.numPredictors; ++j) {
    if (!scratch.available[j]) continue;
    double logSplitVariable =
      splitProbabilities == nullptr
        ? -std::log(static_cast<double>(numAvailable))
        : std::log(splitProbabilities[j] / availableSplitProbability);

    if (data.types[j] == ColumnType::categorical) {
      // R, the categories that reach this node: read off the mask
      // collectAvailableVariables already narrowed for an inline column,
      // walked here for a pooled one (that pass skips those)
      std::size_t numReachable;
      if (data.columnIsPooled(j)) {
        std::size_t numWords = maskWordsForCount(data.numCuts[j]);
        scratch.reachableMask.resize(numWords);
        tree.reachableCategoriesWide(data, nodeIndex,
                                     static_cast<std::int32_t>(j),
                                     scratch.reachableMask.data());
        numReachable = maskPopcount(scratch.reachableMask.data(), numWords);
      } else {
        numReachable = static_cast<std::size_t>(
          std::popcount(tree.inlineReachableMasks()[j]));
      }

      std::size_t numEmitted = scanCategoricalPartitions(
        data, j, tree.indices + begin, numMembers, y, weights, leaf, k,
        residualVariance, scratch.categoryScan, scratch.cutLogLikelihood);
      if (numEmitted == 0) continue;  // fewer than two categories present
      std::size_t numPresent = scratch.categoryScan.present.size();
      // a member arrives only by satisfying every ancestor rule, so its
      // category survived every AND: present implies reachable, A = R - P >= 0
      assert(numPresent <= numReachable);

      double splitBase =
        logGrowth + logSplitVariable +
        CGMTreePrior::categoricalGroupLogProbability(
          numReachable, numPresent, static_cast<double>(numEmitted));
      for (std::size_t candidate = 0; candidate < numEmitted; ++candidate) {
        if (scratch.cutLogLikelihood[candidate] == cutScanEmptySentinel)
          continue;
        scratch.candidateLogWeights.push_back(
          splitBase + scratch.cutLogLikelihood[candidate]);
        scratch.candidateVariable.push_back(static_cast<std::int32_t>(j));
        scratch.candidateCut.push_back(static_cast<std::int32_t>(candidate));
        scratch.candidateMissing.push_back(-1);  // the mask routes its own
      }
      continue;
    }

    // P(rule): uniform over the ancestor-constrained interval, and over the two
    // missing directions where the column carries them. The scan emits both
    // directions of a cut exactly when the node's missing rows make them score
    // differently, so a candidate is ONE rule when it does and the whole pair
    // when it does not, and the mass follows: the halving
    // CGMTreePrior::ruleForVariableLogProbability applies to one rule of a
    // missing-bearing column is spent where the candidate is one rule, and
    // cancels against the post-draw coin where it is the pair.
    std::int32_t left, right;
    tree.splitInterval(data, nodeIndex, static_cast<std::int32_t>(j), &left,
                       &right);
    double logCut = -std::log(static_cast<double>(right - left + 1));

    std::size_t numCuts = data.numCuts[j];
    scratch.cutLogLikelihood.resize(data.hasMissing[j] ? 2 * numCuts : numCuts);
    std::size_t numEmitted =
      scanOrdinalCuts(data, j, tree.indices + begin, numMembers, y, weights,
                      leaf, k, residualVariance, scratch.binScratch,
                      scratch.cutLogLikelihood.data());
    bool routesMissing = numEmitted > numCuts;

    double splitBase = logGrowth + logSplitVariable + logCut;
    if (routesMissing) splitBase -= std::log(2.0);
    for (std::size_t entry = 0; entry < numEmitted; ++entry) {
      if (scratch.cutLogLikelihood[entry] == cutScanEmptySentinel) continue;
      scratch.candidateLogWeights.push_back(splitBase +
                                            scratch.cutLogLikelihood[entry]);
      scratch.candidateVariable.push_back(static_cast<std::int32_t>(j));
      // the doubled layout is (cut, direction) interleaved
      scratch.candidateCut.push_back(
        static_cast<std::int32_t>(routesMissing ? entry >> 1 : entry));
      scratch.candidateMissing.push_back(
        routesMissing ? static_cast<std::int8_t>(entry & std::size_t{1})
                      : std::int8_t{-1});
    }
  }

  // draw one outcome; assemble in log space with a max-shift so the leaf
  // marginals (which can be large and negative) do not underflow before the
  // exponentials, and normalize for ext_rng_drawFromDiscreteDistribution
  std::size_t numCandidates = scratch.candidateLogWeights.size();
  double maxLogWeight = -std::numeric_limits<double>::infinity();
  for (double logWeight : scratch.candidateLogWeights)
    if (logWeight > maxLogWeight) maxLogWeight = logWeight;
  scratch.candidateWeights.resize(numCandidates);
  double sum = 0.0;
  for (std::size_t i = 0; i < numCandidates; ++i) {
    double weight = std::exp(scratch.candidateLogWeights[i] - maxLogWeight);
    scratch.candidateWeights[i] = weight;
    sum += weight;
  }
  for (std::size_t i = 0; i < numCandidates; ++i)
    scratch.candidateWeights[i] /= sum;

  std::size_t choice = ext_rng_drawFromDiscreteDistribution(
    rng, scratch.candidateWeights.data(), numCandidates);
  // choice 0 is no-split; a failure (only under a degenerate all-zero weight
  // vector, which the finite no-split term precludes) also ends the node
  if (choice == 0 || choice == EXT_DISCRETE_DRAW_FAILURE) return;

  Rule rule;
  rule.variableIndex = scratch.candidateVariable[choice];
  std::size_t winner = static_cast<std::size_t>(rule.variableIndex);
  if (data.types[winner] == ColumnType::categorical)
    growCategoricalRule(data, leaf, rng, tree, nodeIndex, y, weights, k,
                        residualVariance, begin, numMembers,
                        scratch.candidateCut[choice], scratch, rule);
  else {
    rule.setSplitIndex(scratch.candidateCut[choice]);
    std::int8_t routed = scratch.candidateMissing[choice];
    if (routed >= 0)
      rule.setMissingGoesRight(routed != 0);  // the candidate already named it
    else if (data.hasMissing[winner])
      rule.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);
  }

  tree.birth(data, nodeIndex, rule, y, weights);
  std::int32_t leftChild = tree.at(nodeIndex).leftChild;
  growTreeFromRoot(data, treePrior, leaf, rng, tree, leftChild, y, weights, k,
                   sigma, scratch);
  growTreeFromRoot(data, treePrior, leaf, rng, tree, leftChild + 1, y, weights,
                   k, sigma, scratch);
}

}  // namespace bartcore

#endif  // BARTCORE_GROW_HPP
