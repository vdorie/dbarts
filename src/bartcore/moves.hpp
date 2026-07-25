#ifndef BARTCORE_MOVES_HPP
#define BARTCORE_MOVES_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <external/random.h>

#include "data.hpp"
#include "model.hpp"
#include "tree.hpp"

// The conjugate Metropolis-Hastings tree moves: birth/death, change, and swap
// proposals with their acceptance ratios.

namespace bartcore {

/// Per-chain scratch for the moves; reused across calls to avoid allocation.
struct MoveScratch {
  std::vector<int32_t> nodeScratch;
  std::vector<double> probabilityScratch;
  Tree::SubtreeSnapshot snapshot;
  // pooled categorical masks: the changed node's reachable set, the pattern
  // draw, and a depth-indexed arena for the validity walk's left-side masks
  std::vector<std::uint64_t> reachableWords;
  std::vector<std::uint64_t> patternWords;
  std::vector<std::uint64_t> maskArena;
};

struct MoveContext {
  const ColumnStore& data;
  const CGMTreePrior& treePrior;
  double birthOrDeathProbability;
  double swapProbability;
  double birthProbability;
  const double* weights;
  double k;
  MoveScratch& scratch;
  // the tree's current leaf-parameter vector M, handed to the scoring of a
  // ParamScoringLeafModel (monotone reads frozen neighbor mu here); null and
  // unread for the conjugate leaves, which integrate every leaf out.
  const double* leafParams = nullptr;
};

template <MoveScorableLeafModel L, typename ResidT = double>
double logLikelihoodForBranch(const MoveContext& ctx, const L& leaf, Tree& tree,
                              int32_t branchIndex, const ResidT* y, double sigma) {
  // a leaf whose score reads M owns the branch marginal outright (the
  // constrained joint over the touched leaves given frozen neighbors); the
  // conjugate leaves fall through to the per-leaf marginal sum unchanged.
  if constexpr (ParamScoringLeafModel<L>)
    return leaf.logLikelihoodForBranchWithParams(tree, branchIndex, y,
                                                 ctx.weights, ctx.k,
                                                 sigma * sigma, ctx.leafParams);

  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(branchIndex, bottoms);

  double result = 0.0;
  for (int32_t i : bottoms) {
    // vetoing any branch with an empty leaf keeps empty leaves out of the
    // chain state (docs/design/empty-leaf-veto.md). The penalty must beat any
    // finite branch score: a valid branch's log-likelihood is unbounded below
    // (it scales with centeredSumOfSquares / residualVariance, large for a big
    // node or a small sigma), so a fixed constant would be out-penalized at
    // scale and the empty leaf accepted. -HUGE_VAL vetoes unconditionally;
    // vetoed-vs-vetoed never arises, so exp of the difference stays 0, not NaN.
    if (tree.at(i).numObservations() == 0) return -HUGE_VAL;
    result += leaf.logIntegratedLikelihoodForNode(tree, y, ctx.weights, ctx.k,
                                                  sigma * sigma, i);
  }
  return result;
}

inline double probabilityOfBirthStep(const MoveContext& ctx, const Tree& tree,
                                     bool birthableNodeExists) {
  if (!birthableNodeExists) return 0.0;
  if (tree.hasSingleNode()) return 1.0;
  return ctx.birthProbability;
}

inline bool birthableNodeExists(const MoveContext& ctx, Tree& tree) {
  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);
  for (int32_t i : bottoms)
    if (tree.hasAnyAvailableVariable(ctx.data, i)) return true;
  return false;
}

inline double probabilityOfSelectingNodeForDeath(Tree& tree,
                                                 std::vector<int32_t>& scratch) {
  scratch.clear();
  tree.fillNoGrand(0, scratch);
  if (scratch.empty()) return 0.0;
  return 1.0 / static_cast<double>(scratch.size());
}

inline double probabilityOfSelectingNodeForBirth(const MoveContext& ctx,
                                                 Tree& tree) {
  if (tree.hasSingleNode()) return 1.0;

  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);

  double totalProbability = 0.0;
  for (int32_t i : bottoms)
    totalProbability += tree.hasAnyAvailableVariable(ctx.data, i) ? 1.0 : 0.0;

  if (totalProbability <= 0.0) return 0.0;
  return 1.0 / totalProbability;
}

inline int32_t drawBirthableNode(const MoveContext& ctx, ext_rng* rng, Tree& tree,
                                 double* nodeSelectionProbability) {
  if (tree.hasSingleNode()) {
    *nodeSelectionProbability = 1.0;
    return 0;
  }

  std::vector<int32_t>& bottoms(ctx.scratch.nodeScratch);
  bottoms.clear();
  tree.fillBottom(0, bottoms);

  std::vector<double>& probabilities(ctx.scratch.probabilityScratch);
  probabilities.resize(bottoms.size());
  double totalProbability = 0.0;
  for (size_t i = 0; i < bottoms.size(); ++i) {
    probabilities[i] =
      tree.hasAnyAvailableVariable(ctx.data, bottoms[i]) ? 1.0 : 0.0;
    totalProbability += probabilities[i];
  }

  if (totalProbability <= 0.0) {
    *nodeSelectionProbability = 0.0;
    return invalidNode;
  }

  misc_scalarMultiplyVectorInPlace(probabilities.data(), probabilities.size(),
                                   1.0 / totalProbability);
  size_t index = ext_rng_drawFromDiscreteDistribution(rng, probabilities.data(),
                                                      probabilities.size());
  *nodeSelectionProbability = probabilities[index];
  return bottoms[index];
}

template <MoveScorableLeafModel L, typename ResidT = double>
double birthOrDeathMove(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                        Tree& tree, const ResidT* y, double sigma,
                        bool* stepTaken, bool* stepWasBirth,
                        int32_t* changedNode = nullptr) {
  // A root-only tree whose lone leaf admits no split variable can neither birth
  // (no rule to draw) nor die (no children); its move is a no-op this sweep.
  // The single-node branch below would otherwise force a birth and draw a rule
  // for invalidVariable. A movable tree never reaches here, so RNG is untouched.
  if (tree.hasSingleNode() && !birthableNodeExists(ctx, tree)) {
    *stepTaken = false;
    *stepWasBirth = false;
    return 0.0;
  }

  double ratio;

  double transitionProbabilityOfSelectingNodeForBirth;
  int32_t nodeToChange =
    drawBirthableNode(ctx, rng, tree, &transitionProbabilityOfSelectingNodeForBirth);

  double transitionProbabilityOfBirthStep =
    probabilityOfBirthStep(ctx, tree, nodeToChange != invalidNode);

  if (ext_rng_simulateBernoulli(rng, transitionProbabilityOfBirthStep) == 1) {
    *stepWasBirth = true;

    double parentPriorGrowthProbability =
      ctx.treePrior.growthProbability(tree, ctx.data, nodeToChange);
    double oldPriorProbability = 1.0 - parentPriorGrowthProbability;
    double oldLogLikelihood =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    // The proposal rule is drawn from the prior, so its density cancels with
    // the prior's in the acceptance ratio.
    Node oldNode = tree.at(nodeToChange);
    size_t maskPoolMark = tree.maskPoolMark();
    Rule newRule = ctx.treePrior.drawRuleAndVariable(tree, ctx.data, rng, nodeToChange);
    tree.birth(ctx.data, nodeToChange, newRule, y, ctx.weights);

    double leftPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild);
    double rightPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild + 1);
    double newPriorProbability = parentPriorGrowthProbability *
                                 (1.0 - leftPriorGrowthProbability) *
                                 (1.0 - rightPriorGrowthProbability);

    double newLogLikelihood =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    double transitionProbabilityOfDeathStep =
      1.0 - probabilityOfBirthStep(ctx, tree, birthableNodeExists(ctx, tree));
    double transitionProbabilityOfSelectingNodeForDeath =
      probabilityOfSelectingNodeForDeath(tree, ctx.scratch.nodeScratch);

    double priorRatio = newPriorProbability / oldPriorProbability;
    double transitionRatio =
      (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath) /
      (transitionProbabilityOfBirthStep * transitionProbabilityOfSelectingNodeForBirth);
    double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);

    ratio = priorRatio * likelihoodRatio * transitionRatio;

    if (ext_rng_simulateContinuousUniform(rng) < ratio) {
      *stepTaken = true;
      if (changedNode != nullptr) *changedNode = nodeToChange;
    } else {
      // Reference behavior: the index segment stays permuted; only structure
      // and cached leaf stats are restored. A rejected pooled draw is the
      // last pool allocation, so the mark reclaims it.
      tree.undoBirth(nodeToChange);
      tree.truncateMaskPool(maskPoolMark);
      tree.at(nodeToChange).sumWeights = oldNode.sumWeights;
      tree.at(nodeToChange).sumWeightedResponse = oldNode.sumWeightedResponse;
      *stepTaken = false;
    }
  } else {
    *stepWasBirth = false;

    double transitionProbabilityOfDeathStep = 1.0 - transitionProbabilityOfBirthStep;

    double transitionProbabilityOfSelectingNodeForDeath;
    std::vector<int32_t>& noGrand(ctx.scratch.nodeScratch);
    noGrand.clear();
    tree.fillNoGrand(0, noGrand);
    size_t index =
      ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, noGrand.size());
    transitionProbabilityOfSelectingNodeForDeath =
      1.0 / static_cast<double>(noGrand.size());
    nodeToChange = noGrand[index];

    double parentPriorGrowthProbability =
      ctx.treePrior.growthProbability(tree, ctx.data, nodeToChange);
    double leftPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild);
    double rightPriorGrowthProbability = ctx.treePrior.growthProbability(
      tree, ctx.data, tree.at(nodeToChange).leftChild + 1);
    double oldLogLikelihood =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

    Node oldNode = tree.at(nodeToChange);
    tree.orphanChildren(nodeToChange);

    double newLogLikelihood =
      logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);
    double transitionProbabilityOfBirthStepReverse =
      probabilityOfBirthStep(ctx, tree, true);
    double reverseTransitionProbabilityOfSelectingNodeForBirth =
      probabilityOfSelectingNodeForBirth(ctx, tree);

    double oldPriorProbability = parentPriorGrowthProbability *
                                 (1.0 - leftPriorGrowthProbability) *
                                 (1.0 - rightPriorGrowthProbability);
    double newPriorProbability = 1.0 - parentPriorGrowthProbability;

    double priorRatio = newPriorProbability / oldPriorProbability;
    double transitionRatio =
      (transitionProbabilityOfBirthStepReverse * reverseTransitionProbabilityOfSelectingNodeForBirth) /
      (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath);
    double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);

    ratio = priorRatio * likelihoodRatio * transitionRatio;

    if (ext_rng_simulateContinuousUniform(rng) < ratio) {
      tree.releasePair(oldNode.leftChild);
      *stepTaken = true;
      if (changedNode != nullptr) *changedNode = nodeToChange;
    } else {
      tree.at(nodeToChange) = oldNode;  // reattaches children unchanged
      *stepTaken = false;
    }
  }

  return ratio < 1.0 ? ratio : 1.0;
}

/// Ordinal cut range at nodeIndex that keeps every descendant split on the
/// same variable logically satisfiable (findGoodOrdinalRules).
inline void findGoodOrdinalRules(const MoveContext& ctx, const Tree& tree,
                                 int32_t nodeIndex, int32_t variableIndex,
                                 int32_t* lower, int32_t* upper) {
  int32_t leftIndex, rightIndex;
  tree.splitInterval(ctx.data, nodeIndex, variableIndex, &leftIndex, &rightIndex);

  // min/max split index used for the variable in each child's subtree
  struct Walker {
    static void minMax(const Tree& tree, int32_t i, int32_t variableIndex,
                       int32_t* min, int32_t* max) {
      const Node& node(tree.at(i));
      if (node.isBottom()) return;
      if (node.rule.variableIndex == variableIndex) {
        if (node.rule.splitIndex() < *min) *min = node.rule.splitIndex();
        if (node.rule.splitIndex() > *max) *max = node.rule.splitIndex();
      }
      minMax(tree, node.leftChild, variableIndex, min, max);
      minMax(tree, node.leftChild + 1, variableIndex, min, max);
    }
  };

  int32_t leftMin = rightIndex + 1, leftMaxOut = leftIndex - 1;
  int32_t rightMinOut = rightIndex + 1, rightMax = leftIndex - 1;
  Walker::minMax(tree, tree.at(nodeIndex).leftChild, variableIndex, &leftMin,
                 &leftMaxOut);
  Walker::minMax(tree, tree.at(nodeIndex).leftChild + 1, variableIndex,
                 &rightMinOut, &rightMax);
  int32_t leftMax = leftMaxOut;
  int32_t rightMin = rightMinOut;

  *lower = std::max(leftIndex, leftMax + 1);
  *upper = std::min(rightIndex, rightMin - 1);
}

inline bool categoricalSubtreeIsValid(const Tree& tree, int32_t nodeIndex,
                                      int32_t variableIndex,
                                      std::uint64_t reachable);
inline bool categoricalSubtreeIsValidWide(const Tree& tree, int32_t nodeIndex,
                                          int32_t variableIndex,
                                          const std::uint64_t* reachable,
                                          size_t numWords,
                                          std::uint64_t* arena, size_t depth);

/// Draw a categorical assignment straight from the node prior (the
/// propose-from-prior mechanism; see docs/design/change-move-balance.md): a
/// single unrestricted gauge-pattern draw over the
/// reachable set, with no descendant-validity rejection loop. Returns true and
/// fills newRule when the draw leaves every same-variable descendant
/// satisfiable; returns false (pi(T') = 0, an automatic no-op) otherwise. The
/// prior density 1/(2^R - 2) cancels the node's rule prior exactly, so the
/// acceptance keeps only the subtree-below and likelihood ratios.
inline bool drawCategoricalRuleFromPrior(const MoveContext& ctx, ext_rng* rng,
                                         Tree& tree, int32_t nodeToChange,
                                         int32_t newVariableIndex, Rule& newRule,
                                         size_t maskPoolMark) {
  int32_t leftChild = tree.at(nodeToChange).leftChild;
  if (ctx.data.columnIsPooled(static_cast<size_t>(newVariableIndex))) {
    size_t numWords = maskWordsForCount(
      ctx.data.numCuts[static_cast<size_t>(newVariableIndex)]);
    std::vector<std::uint64_t>& reachable(ctx.scratch.reachableWords);
    reachable.resize(numWords);
    tree.reachableCategoriesWide(ctx.data, nodeToChange, newVariableIndex,
                                 reachable.data());
    size_t numReachable = maskPopcount(reachable.data(), numWords);
    ctx.scratch.patternWords.resize(numWords);
    ctx.scratch.maskArena.resize((tree.nodes.size() + 1) * numWords);
    size_t offset = tree.allocateMask(numWords);
    CGMTreePrior::drawCategoryPatternWide(
      rng, numReachable, ctx.scratch.patternWords.data(), numWords);
    std::uint64_t* directions = tree.mutableMaskWordsFor(offset);
    CGMTreePrior::categoryDirectionsForPatternWide(
      reachable.data(), ctx.scratch.patternWords.data(), directions, numWords);
    std::uint64_t* leftReachable = ctx.scratch.maskArena.data();
    maskAndNot(reachable.data(), directions, leftReachable, numWords);
    if (!categoricalSubtreeIsValidWide(tree, leftChild, newVariableIndex,
                                       leftReachable, numWords,
                                       ctx.scratch.maskArena.data(), 1) ||
        !categoricalSubtreeIsValidWide(tree, leftChild + 1, newVariableIndex,
                                       directions, numWords,
                                       ctx.scratch.maskArena.data(), 1)) {
      tree.truncateMaskPool(maskPoolMark);
      return false;
    }
    newRule.setMaskOffset(offset);
    return true;
  }

  std::uint64_t reachable =
    tree.reachableCategories(ctx.data, nodeToChange, newVariableIndex);
  int numReachable = std::popcount(reachable);
  std::uint64_t pattern = CGMTreePrior::drawCategoryPattern(rng, numReachable);
  std::uint64_t directions =
    CGMTreePrior::categoryDirectionsForPattern(reachable, pattern);
  if (!categoricalSubtreeIsValid(tree, leftChild, newVariableIndex,
                                 reachable & ~directions) ||
      !categoricalSubtreeIsValid(tree, leftChild + 1, newVariableIndex,
                                 directions))
    return false;  // narrow columns allocate no pool words to reclaim
  newRule.setCategoryDirections(directions);
  return true;
}

/// Change-move proposal kernel: redraw the split variable and rule at an
/// internal node, keeping the subtree below in place. The acceptance satisfies
/// detailed balance,
///   alpha = exp( B(T') - B(T) + yLogL - xLogL + logProposalCorrection ),
/// where B is the tree-prior log-probability of the subtree STRICTLY BELOW the
/// changed node; every prior factor at or above the node cancels between T and
/// T', as in birth/death. The correction is the surviving proposal-density
/// ratio q(T|T')/q(T'|T) and composes PER SIDE, because the forward density
/// uses the new variable v''s mechanism and the reverse the old variable v's:
///   logProposalCorrection =
///       (v' ordinal ? log|Valid_T(v')| - log|SI(v')| : 0)
///     + (v  ordinal ? log|SI(v)| - log|Valid_T'(v)| : 0).
/// An ordinal side draws uniformly over the descendant-valid good set while
/// the node's rule prior normalizes over the ancestor-only interval, leaving
/// the counted ratio (|SI| the interval size, |Valid| the good-set count; the
/// variable prior and missing coin cancel within the side). A categorical side
/// draws straight from the node prior (drawCategoricalRuleFromPrior), whose
/// density cancels its side's prior factor exactly and contributes nothing; a
/// forward draw that strands a descendant split is an automatic no-op
/// (pi(T') = 0). findGoodOrdinalRules and splitInterval both ignore the node's
/// OWN rule, so the reverse ordinal count, re-enumerated on the current tree,
/// equals its value on the changed tree and always contains the old rule
/// (>= 1, never zero); a same-variable redraw gives correction 1. Omitting the
/// correction is the CGM-lineage defect this repairs
/// (docs/design/change-move-balance.md).
template <MoveScorableLeafModel L, typename ResidT = double>
double changeMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                  const ResidT* y, double sigma, bool* stepTaken,
                  int32_t* changedNode = nullptr) {
  *stepTaken = false;

  std::vector<int32_t>& notBottom(ctx.scratch.nodeScratch);
  notBottom.clear();
  tree.fillNotBottom(0, notBottom);
  if (notBottom.empty()) return -1.0;

  size_t nodeNumber =
    ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, notBottom.size());
  int32_t nodeToChange = notBottom[nodeNumber];

  int32_t newVariableIndex =
    ctx.treePrior.drawSplitVariable(tree, ctx.data, rng, nodeToChange);

  Rule newRule;
  newRule.variableIndex = newVariableIndex;
  // covers every exit: a pooled proposal's words are the last allocation,
  // so aborts and MH rejections truncate back; acceptance keeps them (the
  // replaced rule's words become garbage the chain compacts later)
  size_t maskPoolMark = tree.maskPoolMark();

  const bool newIsCategorical =
    ctx.data.types[static_cast<size_t>(newVariableIndex)] ==
    ColumnType::categorical;
  int32_t oldVariableIndex = tree.at(nodeToChange).rule.variableIndex;
  const bool oldIsCategorical =
    ctx.data.types[static_cast<size_t>(oldVariableIndex)] ==
    ColumnType::categorical;
  int32_t forwardValid = 0, forwardInterval = 0;  // new-side ordinal counts
  int32_t reverseValid = 0, reverseInterval = 0;  // old-side ordinal counts

  if (newIsCategorical) {
    // counting descendant-valid gauge patterns is exponential for wide masks,
    // so a categorical proposal draws from the prior; the density cancels the
    // node's rule prior and its side contributes no correction
    if (!drawCategoricalRuleFromPrior(ctx, rng, tree, nodeToChange,
                                      newVariableIndex, newRule, maskPoolMark))
      return -1.0;  // pi(T') = 0: an unsatisfiable prior draw is a no-op
  } else {
    int32_t left, right;
    tree.splitInterval(ctx.data, nodeToChange, newVariableIndex, &left, &right);
    int32_t lower, upper;
    findGoodOrdinalRules(ctx, tree, nodeToChange, newVariableIndex, &lower, &upper);
    if (upper - lower + 1 <= 0) return -1.0;

    newRule.setSplitIndex(static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, lower, upper + 1)));
    // like the birth draw: the missing direction is a fresh symmetric coin
    // whenever the column can route a missing value
    if (ctx.data.hasMissing[static_cast<size_t>(newVariableIndex)])
      newRule.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);

    forwardValid = upper - lower + 1;
    forwardInterval = right - left + 1;
  }

  if (!oldIsCategorical) {
    // the reverse count re-enumerates the OLD variable on the current tree;
    // splitInterval and findGoodOrdinalRules both ignore the node's own rule,
    // so it always contains the old rule and never vanishes. Never run these
    // ordinal counters on a categorical column - a categorical reverse side
    // is a prior draw whose density cancels, contributing nothing.
    int32_t oldLeft, oldRight, oldLower, oldUpper;
    tree.splitInterval(ctx.data, nodeToChange, oldVariableIndex, &oldLeft,
                       &oldRight);
    findGoodOrdinalRules(ctx, tree, nodeToChange, oldVariableIndex, &oldLower,
                         &oldUpper);
    reverseValid = oldUpper - oldLower + 1;
    reverseInterval = oldRight - oldLeft + 1;
  }

  // interaction constraint: the variable was drawn feasibly against
  // nodeToChange's ANCESTORS, but the unchanged subtree below may now strand a
  // descendant (a co-occurrence or order break the drawn variable introduces).
  // treeLogProbability delegates to the availability primitives, which cannot
  // self-detect a node whose OWN variable is barred, so score the -1.0 no-op
  // (pi(T') = 0) directly, exactly as the unsatisfiable categorical draw does.
  if (tree.hasInteractionConstraint()) {
    Rule savedRule = tree.at(nodeToChange).rule;
    tree.at(nodeToChange).rule = newRule;
    bool valid = tree.interactionSubtreeIsValid(nodeToChange);
    tree.at(nodeToChange).rule = savedRule;
    if (!valid) {
      tree.truncateMaskPool(maskPoolMark);
      return -1.0;
    }
  }

  double logProposalCorrection = 0.0;
  if (!newIsCategorical && !oldIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(reverseInterval)) -
      std::log(static_cast<double>(forwardInterval)) +
      std::log(static_cast<double>(forwardValid)) -
      std::log(static_cast<double>(reverseValid));
  } else if (!newIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(forwardValid)) -
      std::log(static_cast<double>(forwardInterval));
  } else if (!oldIsCategorical) {
    logProposalCorrection =
      std::log(static_cast<double>(reverseInterval)) -
      std::log(static_cast<double>(reverseValid));
  }

  // the node's own split-variable and rule-prior factors cancel against the
  // proposal (or are carried by logProposalCorrection), so the pi ratio
  // reduces to the subtree strictly below the changed node
  int32_t leftChild = tree.at(nodeToChange).leftChild;
  double belowX =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  double xLogL = logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  tree.snapshotSubtree(nodeToChange, ctx.scratch.snapshot);

  tree.at(nodeToChange).rule = newRule;
  tree.refreshSubtree(ctx.data, nodeToChange, y, ctx.weights);

  double belowY =
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild) +
    ctx.treePrior.treeLogProbability(tree, ctx.data, leftChild + 1);
  double yLogL = logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  double alpha =
    std::exp((belowY - belowX) + (yLogL - xLogL) + logProposalCorrection);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
    if (changedNode != nullptr) *changedNode = nodeToChange;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
    tree.truncateMaskPool(maskPoolMark);
  }
  return alpha;
}

/// ordinalRuleIsValid: every descendant split on variableIndex must fall
/// inside the interval implied by its ancestors.
inline bool ordinalRuleIsValid(const Tree& tree, int32_t nodeIndex,
                               int32_t variableIndex, int32_t leftIndex,
                               int32_t rightIndex) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    int32_t splitIndex = node.rule.splitIndex();
    if (splitIndex < leftIndex || splitIndex > rightIndex) return false;
    return ordinalRuleIsValid(tree, node.leftChild, variableIndex, leftIndex,
                              splitIndex - 1) &&
           ordinalRuleIsValid(tree, node.leftChild + 1, variableIndex,
                              splitIndex + 1, rightIndex);
  }

  return ordinalRuleIsValid(tree, node.leftChild, variableIndex, leftIndex,
                            rightIndex) &&
         ordinalRuleIsValid(tree, node.leftChild + 1, variableIndex, leftIndex,
                            rightIndex);
}

/// categoricalSubtreeIsValid: every split on variableIndex in the subtree
/// must stay in the canonical gauge (its directions confined to the
/// categories reaching it) and keep at least one reachable category on each
/// side; reachable is the mask entering the subtree.
inline bool categoricalSubtreeIsValid(const Tree& tree, int32_t nodeIndex,
                                      int32_t variableIndex,
                                      std::uint64_t reachable) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    std::uint64_t directions = node.rule.categoryDirections();
    if ((directions & ~reachable) != 0 || directions == 0 ||
        directions == reachable)
      return false;
    return categoricalSubtreeIsValid(tree, node.leftChild, variableIndex,
                                     reachable & ~directions) &&
           categoricalSubtreeIsValid(tree, node.leftChild + 1, variableIndex,
                                     directions);
  }

  return categoricalSubtreeIsValid(tree, node.leftChild, variableIndex,
                                   reachable) &&
         categoricalSubtreeIsValid(tree, node.leftChild + 1, variableIndex,
                                   reachable);
}

/// The pooled-column analogue: reachable spans numWords words, the left
/// branch's filtered set lives in the arena's slot at depth (deeper matching
/// rules use higher slots, so a frame's mask survives its subtree), and the
/// right branch reads the rule's own immutable pool words.
inline bool categoricalSubtreeIsValidWide(const Tree& tree, int32_t nodeIndex,
                                          int32_t variableIndex,
                                          const std::uint64_t* reachable,
                                          size_t numWords,
                                          std::uint64_t* arena, size_t depth) {
  const Node& node(tree.at(nodeIndex));
  if (node.isBottom()) return true;

  if (node.rule.variableIndex == variableIndex) {
    const std::uint64_t* directions = tree.maskWordsFor(node.rule);
    if (!maskIsSubsetOf(directions, reachable, numWords) ||
        maskIsZero(directions, numWords) ||
        maskEquals(directions, reachable, numWords))
      return false;
    std::uint64_t* leftReachable = arena + depth * numWords;
    maskAndNot(reachable, directions, leftReachable, numWords);
    return categoricalSubtreeIsValidWide(tree, node.leftChild, variableIndex,
                                         leftReachable, numWords, arena,
                                         depth + 1) &&
           categoricalSubtreeIsValidWide(tree, node.leftChild + 1,
                                         variableIndex, directions, numWords,
                                         arena, depth + 1);
  }

  return categoricalSubtreeIsValidWide(tree, node.leftChild, variableIndex,
                                       reachable, numWords, arena, depth) &&
         categoricalSubtreeIsValidWide(tree, node.leftChild + 1, variableIndex,
                                       reachable, numWords, arena, depth);
}

inline bool ruleIsValid(const MoveContext& ctx, const Tree& tree, int32_t nodeIndex,
                        int32_t variableIndex) {
  if (ctx.data.types[static_cast<size_t>(variableIndex)] ==
      ColumnType::categorical) {
    size_t j = static_cast<size_t>(variableIndex);
    if (ctx.data.columnIsPooled(j)) {
      size_t numWords = maskWordsForCount(ctx.data.numCuts[j]);
      std::vector<std::uint64_t>& reachable(ctx.scratch.reachableWords);
      reachable.resize(numWords);
      tree.reachableCategoriesWide(ctx.data, nodeIndex, variableIndex,
                                   reachable.data());
      ctx.scratch.maskArena.resize((tree.nodes.size() + 1) * numWords);
      return categoricalSubtreeIsValidWide(tree, nodeIndex, variableIndex,
                                           reachable.data(), numWords,
                                           ctx.scratch.maskArena.data(), 0);
    }
    return categoricalSubtreeIsValid(
      tree, nodeIndex, variableIndex,
      tree.reachableCategories(ctx.data, nodeIndex, variableIndex));
  }

  int32_t leftIndex, rightIndex;
  tree.splitInterval(ctx.data, nodeIndex, variableIndex, &leftIndex, &rightIndex);
  return ordinalRuleIsValid(tree, nodeIndex, variableIndex, leftIndex, rightIndex);
}

template <MoveScorableLeafModel L, typename ResidT = double>
double swapMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                const ResidT* y, double sigma, bool* stepTaken,
                int32_t* changedNode = nullptr) {
  *stepTaken = false;

  std::vector<int32_t>& swappable(ctx.scratch.nodeScratch);
  swappable.clear();
  tree.fillSwappable(0, swappable);
  if (swappable.empty()) return -1.0;

  size_t nodeNumber =
    ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, swappable.size());
  int32_t parent = swappable[nodeNumber];
  int32_t leftChild = tree.at(parent).leftChild;
  int32_t rightChild = leftChild + 1;

  bool leftHasRule = !tree.at(leftChild).isBottom();
  bool rightHasRule = !tree.at(rightChild).isBottom();
  bool childrenHaveSameRule = leftHasRule && rightHasRule &&
    tree.rulesAreEqual(ctx.data, tree.at(leftChild).rule,
                       tree.at(rightChild).rule);

  double alpha;

  // The swap gives the parent a child's rule and each swapped child the
  // parent's. When the two children share a rule both are swapped; otherwise
  // one is picked (a fair coin when both carry a rule, else the only one that
  // does). Either way every swapped child's original rule is childRule, so the
  // two cases differ only in which children move - captured in swapChildren.
  Rule parentRule = tree.at(parent).rule;
  Rule childRule;
  int32_t swapChildren[2];
  int numSwapChildren;
  if (!childrenHaveSameRule) {
    int32_t child;
    if (leftHasRule && rightHasRule) {
      child = ext_rng_simulateBernoulli(rng, 0.5) == 1 ? leftChild : rightChild;
    } else {
      child = leftHasRule ? leftChild : rightChild;
    }
    childRule = tree.at(child).rule;
    swapChildren[0] = child;
    numSwapChildren = 1;
  } else {
    childRule = tree.at(leftChild).rule;
    swapChildren[0] = leftChild;
    swapChildren[1] = rightChild;
    numSwapChildren = 2;
  }

  auto applySwap = [&]() {
    tree.at(parent).rule = childRule;
    for (int i = 0; i < numSwapChildren; ++i)
      tree.at(swapChildren[i]).rule = parentRule;
  };
  auto undoSwap = [&]() {
    tree.at(parent).rule = parentRule;
    for (int i = 0; i < numSwapChildren; ++i)
      tree.at(swapChildren[i]).rule = childRule;
  };

  // test the swap for logical consistency before scoring it
  applySwap();
  bool swapIsSensible = ruleIsValid(ctx, tree, parent, childRule.variableIndex);
  if (childRule.variableIndex != parentRule.variableIndex && swapIsSensible)
    swapIsSensible = ruleIsValid(ctx, tree, parent, parentRule.variableIndex);
  // interaction is a WHOLE-subtree, all-variables property the per-variable
  // ruleIsValid checks above cannot see (the swap sibling-strand break): a
  // swap that lifts x2 above x3 co-occurs a forbidden pair with neither
  // swapped variable equal to x3. Score it the -1.0 no-op (pi(T') = 0).
  if (swapIsSensible) swapIsSensible = tree.interactionSubtreeIsValid(parent);
  undoSwap();

  if (!swapIsSensible) return -1.0;

  // as in changeMove, prior terms outside the swapped subtree cancel
  double xLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data, parent);
  double xLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

  tree.snapshotSubtree(parent, ctx.scratch.snapshot);
  applySwap();
  tree.refreshSubtree(ctx.data, parent, y, ctx.weights);

  double yLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data, parent);
  double yLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

  alpha = std::exp(yLogPi + yLogL - xLogPi - xLogL);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
    if (changedNode != nullptr) *changedNode = parent;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
  }

  return alpha;
}

enum class StepType { birth, death, swap, change };

/// changedNode, when non-null, receives the index of the node whose subtree an
/// ACCEPTED move repartitioned (the birthed/died node, or the changed/swapped
/// subtree root); untouched on rejection or no-op, so gate reads on stepTaken.
template <MoveScorableLeafModel L, typename ResidT = double>
double metropolisJumpForTree(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                             Tree& tree, const ResidT* y, double sigma,
                             bool* stepTaken, StepType* stepType,
                             int32_t* changedNode = nullptr) {
  double alpha;
  double u = ext_rng_simulateContinuousUniform(rng);

  if (u < ctx.birthOrDeathProbability) {
    bool birthed;
    alpha = birthOrDeathMove(ctx, leaf, rng, tree, y, sigma, stepTaken, &birthed,
                             changedNode);
    *stepType = birthed ? StepType::birth : StepType::death;
  } else if (u < ctx.birthOrDeathProbability + ctx.swapProbability) {
    alpha = swapMove(ctx, leaf, rng, tree, y, sigma, stepTaken, changedNode);
    *stepType = StepType::swap;
  } else {
    alpha = changeMove(ctx, leaf, rng, tree, y, sigma, stepTaken, changedNode);
    *stepType = StepType::change;
  }

  return alpha;
}

}  // namespace bartcore

#endif  // BARTCORE_MOVES_HPP
