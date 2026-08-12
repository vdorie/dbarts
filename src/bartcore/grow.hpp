#ifndef BARTCORE_GROW_HPP
#define BARTCORE_GROW_HPP

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
// kernel: at each node it scans every ordinal cut of every available variable
// (scan.hpp) and draws one outcome from the weighted set {no-split} U {(var,
// cut)}, the CGM tree prior's own factors times the integrated leaf likelihood.
// It is NOT a posterior kernel; the exact MH sweeps that follow own
// stationarity, so the grown forest need only be a LEGAL chain state (non-empty
// leaves, rules in gauge, availability/depth vetoes respected).

namespace bartcore {

/// Per-node reusable scratch for growTreeFromRoot; one instance per chain,
/// reused down the whole recursion (a node's assembly is over before its
/// children recurse, so nothing aliases).
struct GrowScratch {
  std::vector<ConstantLeafScanBin> binScratch;  // scan.hpp histogram
  std::vector<std::uint8_t> available;          // per-predictor availability
  std::vector<double> cutLogLikelihood;         // per-variable cut scores
  std::vector<double> candidateLogWeights;      // {no-split} U {valid (var,cut)}
  std::vector<double> candidateWeights;         // exp'd, normalized for the draw
  std::vector<std::int32_t> candidateVariable;  // parallel to the weights
  std::vector<std::int32_t> candidateCut;
};

/// Recursively split `tree` from nodeIndex against the residual y (per-tree
/// working response) and case weights. The leaf's integrated likelihood scores
/// each candidate cut through the scan; the discrete outcome is drawn on rng.
///
/// Draw discipline (documented so the equivalence picture stays predictable):
///  - A node whose prior growth probability is zero - the depth veto or no
///    available variable - draws NOTHING and stays a leaf.
///  - Every other node draws exactly ONE discrete outcome over {no-split} plus
///    the occupancy-nonempty ordinal cuts of every available variable. A
///    no-split outcome (always representable) ends the node.
///  - A split on a column with missing values draws ONE additional symmetric
///    missing-direction coin, matching CGMTreePrior::drawRuleForVariable.
///    Convention on such a column, measured rather than assumed: an enumerated
///    cut candidate STANDS FOR the two rules {cut, missing left} and {cut,
///    missing right} the coin picks between, but CARRIES the prior mass of one
///    of them (the `- log 2` below), so the pair enters the discrete draw at
///    half its group's prior mass and the remainder accrues to no-split.
///    test_grow.cpp chi-squares the realized root-rule frequencies against
///    that law, against the group's own mass, and against the exact law over
///    the full rule set; docs/design/grow-from-root.md section 7 records what
///    it measured.
///
/// Categorical predictors are ordinal-only in v1: they are never scanned and
/// never split here (that keeps the draw count above exact), so their structure
/// is left for the MH sweeps that refine the warm start to discover. Occupancy
/// (scan.hpp) subsumes the ancestor split interval - members violating an
/// ancestor cut have codes confined to the interval, so out-of-interval cuts
/// have an empty side and are already zeroed - and keeps both children
/// non-empty, so the built tree satisfies MH's structural invariants.
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
  // candidate 0 is always no-split: (1 - growth) * L(node), always finite
  scratch.candidateLogWeights.push_back(std::log(1.0 - growth) +
                                        nodeLogLikelihood);
  scratch.candidateVariable.push_back(invalidVariable);
  scratch.candidateCut.push_back(-1);

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
    if (!scratch.available[j] || data.types[j] == ColumnType::categorical)
      continue;
    double logSplitVariable =
      splitProbabilities == nullptr
        ? -std::log(static_cast<double>(numAvailable))
        : std::log(splitProbabilities[j] / availableSplitProbability);

    // P(cut): uniform over the ancestor-constrained interval, widened by the
    // missing-direction coin when the column can route one (CGMTreePrior's
    // ruleForVariableLogProbability)
    std::int32_t left, right;
    tree.splitInterval(data, nodeIndex, static_cast<std::int32_t>(j), &left,
                       &right);
    double logCut = -std::log(static_cast<double>(right - left + 1));
    if (data.hasMissing[j]) logCut -= std::log(2.0);

    std::size_t numCuts = data.numCuts[j];
    scratch.cutLogLikelihood.resize(numCuts);
    scanOrdinalCuts(data, j, tree.indices + begin, numMembers, y, weights, leaf,
                    k, residualVariance, scratch.binScratch,
                    scratch.cutLogLikelihood.data());

    double splitBase = logGrowth + logSplitVariable + logCut;
    for (std::size_t cut = 0; cut < numCuts; ++cut) {
      if (scratch.cutLogLikelihood[cut] == cutScanEmptySentinel) continue;
      scratch.candidateLogWeights.push_back(splitBase +
                                            scratch.cutLogLikelihood[cut]);
      scratch.candidateVariable.push_back(static_cast<std::int32_t>(j));
      scratch.candidateCut.push_back(static_cast<std::int32_t>(cut));
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
  rule.setSplitIndex(scratch.candidateCut[choice]);
  if (data.hasMissing[static_cast<std::size_t>(rule.variableIndex)])
    rule.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);

  tree.birth(data, nodeIndex, rule, y, weights);
  std::int32_t leftChild = tree.at(nodeIndex).leftChild;
  growTreeFromRoot(data, treePrior, leaf, rng, tree, leftChild, y, weights, k,
                   sigma, scratch);
  growTreeFromRoot(data, treePrior, leaf, rng, tree, leftChild + 1, y, weights,
                   k, sigma, scratch);
}

}  // namespace bartcore

#endif  // BARTCORE_GROW_HPP
