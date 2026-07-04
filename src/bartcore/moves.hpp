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

// Faithful port of the classic conjugate Metropolis-Hastings tree moves
// (birth/death, change, swap) from src/dbarts/{birthDeathRule,changeRule,
// swapRule}.cpp. Deviations from the reference are noted inline; everything
// else preserves the reference's math and proposal structure.

namespace bartcore {

/// Per-chain scratch for the moves; reused across calls to avoid allocation.
struct MoveScratch {
  std::vector<int32_t> nodeScratch;
  std::vector<double> probabilityScratch;
  Tree::SubtreeSnapshot snapshot;
};

struct MoveContext {
  const ColumnStore& data;
  const CGMTreePrior& treePrior;
  double birthOrDeathProbability;
  double swapProbability;
  double changeProbability;
  double birthProbability;
  const double* weights;
  double k;
  MoveScratch& scratch;
};

template <IntegrableLeafModel L>
double logLikelihoodForBranch(const MoveContext& ctx, const L& leaf, Tree& tree,
                              int32_t branchIndex, const double* y, double sigma) {
  std::vector<int32_t>& bottoms(tree.bottomScratch);
  bottoms.clear();
  tree.fillBottom(branchIndex, bottoms);

  double result = 0.0;
  for (int32_t i : bottoms) {
    // the reference engine vetoes any branch containing an empty leaf, which
    // keeps empty leaves out of the chain state entirely
    if (tree.at(i).numObservations() == 0) return -10000000.0;
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
    if (tree.numVariablesAvailable(ctx.data, i) > 0) return true;
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
    totalProbability += tree.numVariablesAvailable(ctx.data, i) > 0 ? 1.0 : 0.0;

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
      tree.numVariablesAvailable(ctx.data, bottoms[i]) > 0 ? 1.0 : 0.0;
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

template <IntegrableLeafModel L>
double birthOrDeathMove(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                        Tree& tree, const double* y, double sigma,
                        bool* stepTaken, bool* stepWasBirth) {
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
    } else {
      // Reference behavior: the index segment stays permuted; only structure
      // and cached leaf stats are restored.
      tree.undoBirth(nodeToChange);
      tree.at(nodeToChange).average = oldNode.average;
      tree.at(nodeToChange).numEffectiveObservations =
        oldNode.numEffectiveObservations;
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
    double transitionProbabilityOfSelectingNodeForBirth =
      probabilityOfSelectingNodeForBirth(ctx, tree);

    double oldPriorProbability = parentPriorGrowthProbability *
                                 (1.0 - leftPriorGrowthProbability) *
                                 (1.0 - rightPriorGrowthProbability);
    double newPriorProbability = 1.0 - parentPriorGrowthProbability;

    double priorRatio = newPriorProbability / oldPriorProbability;
    double transitionRatio =
      (transitionProbabilityOfBirthStepReverse * transitionProbabilityOfSelectingNodeForBirth) /
      (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath);
    double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);

    ratio = priorRatio * likelihoodRatio * transitionRatio;

    if (ext_rng_simulateContinuousUniform(rng) < ratio) {
      tree.releasePair(oldNode.leftChild);
      *stepTaken = true;
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

  int32_t leftMax = leftIndex - 1;
  int32_t rightMin = rightIndex + 1;

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
  leftMax = leftMaxOut;
  rightMin = rightMinOut;

  *lower = std::max(leftIndex, leftMax + 1);
  *upper = std::min(rightIndex, rightMin - 1);
}

inline bool categoricalSubtreeIsValid(const Tree& tree, int32_t nodeIndex,
                                      int32_t variableIndex,
                                      std::uint64_t reachable);

template <IntegrableLeafModel L>
double changeMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                  const double* y, double sigma, bool* stepTaken) {
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

  if (ctx.data.types[static_cast<size_t>(newVariableIndex)] ==
      ColumnType::categorical) {
    std::uint64_t reachable =
      tree.reachableCategories(ctx.data, nodeToChange, newVariableIndex);
    int numReachable = std::popcount(reachable);
    // Rejection-sample an assignment that keeps the descendant splits on
    // this variable satisfiable: uniform over that good set, which depends
    // only on the tree above and below the node, so the forward and reverse
    // proposals cancel; the abort probability is likewise invariant, so an
    // exhausted draw backs the move out symmetrically.
    bool found = false;
    for (int attempt = 0; attempt < 64 && !found; ++attempt) {
      std::uint64_t pattern =
        CGMTreePrior::drawCategoryPattern(rng, numReachable);
      std::uint64_t directions =
        CGMTreePrior::categoryDirectionsForPattern(reachable, pattern);
      if (categoricalSubtreeIsValid(tree, tree.at(nodeToChange).leftChild,
                                    newVariableIndex, reachable & ~directions) &&
          categoricalSubtreeIsValid(tree, tree.at(nodeToChange).leftChild + 1,
                                    newVariableIndex, directions)) {
        newRule.setCategoryDirections(directions);
        found = true;
      }
    }
    if (!found) return -1.0;
  } else {
    int32_t lower, upper;
    findGoodOrdinalRules(ctx, tree, nodeToChange, newVariableIndex, &lower, &upper);
    if (upper - lower + 1 <= 0) return -1.0;

    newRule.setSplitIndex(static_cast<int32_t>(
      ext_rng_simulateIntegerUniformInRange(rng, lower, upper + 1)));
    // like the birth draw: the missing direction is a fresh symmetric coin
    // whenever the column can route a missing value
    if (ctx.data.hasMissing[static_cast<size_t>(newVariableIndex)])
      newRule.setMissingGoesRight(ext_rng_simulateBernoulli(rng, 0.5) == 1);
  }

  double xLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
  double xLogL = logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  tree.snapshotSubtree(nodeToChange, ctx.scratch.snapshot);

  tree.at(nodeToChange).rule = newRule;
  tree.refreshSubtree(ctx.data, nodeToChange, y, ctx.weights);

  double yLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
  double yLogL = logLikelihoodForBranch(ctx, leaf, tree, nodeToChange, y, sigma);

  double alpha = std::exp(yLogPi + yLogL - xLogPi - xLogL);
  alpha = alpha > 1.0 ? 1.0 : alpha;

  if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
    *stepTaken = true;
  } else {
    tree.restoreSubtree(ctx.scratch.snapshot);
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

inline bool ruleIsValid(const MoveContext& ctx, const Tree& tree, int32_t nodeIndex,
                        int32_t variableIndex) {
  if (ctx.data.types[static_cast<size_t>(variableIndex)] ==
      ColumnType::categorical)
    return categoricalSubtreeIsValid(
      tree, nodeIndex, variableIndex,
      tree.reachableCategories(ctx.data, nodeIndex, variableIndex));

  int32_t leftIndex, rightIndex;
  tree.splitInterval(ctx.data, nodeIndex, variableIndex, &leftIndex, &rightIndex);
  return ordinalRuleIsValid(tree, nodeIndex, variableIndex, leftIndex, rightIndex);
}

template <IntegrableLeafModel L>
double swapMove(const MoveContext& ctx, const L& leaf, ext_rng* rng, Tree& tree,
                const double* y, double sigma, bool* stepTaken) {
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
    tree.at(leftChild).rule.equals(tree.at(rightChild).rule);

  double alpha;

  if (!childrenHaveSameRule) {
    int32_t child;
    if (leftHasRule && rightHasRule) {
      child = ext_rng_simulateBernoulli(rng, 0.5) == 1 ? leftChild : rightChild;
    } else {
      child = leftHasRule ? leftChild : rightChild;
    }

    Rule parentRule = tree.at(parent).rule;
    Rule childRule = tree.at(child).rule;

    // test the swap for logical consistency before scoring it
    tree.at(parent).rule = childRule;
    tree.at(child).rule = parentRule;
    bool swapIsSensible = ruleIsValid(ctx, tree, parent, childRule.variableIndex);
    if (childRule.variableIndex != parentRule.variableIndex && swapIsSensible)
      swapIsSensible = ruleIsValid(ctx, tree, parent, parentRule.variableIndex);
    tree.at(parent).rule = parentRule;
    tree.at(child).rule = childRule;

    if (!swapIsSensible) return -1.0;

    double xLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
    double xLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

    tree.snapshotSubtree(parent, ctx.scratch.snapshot);

    tree.at(parent).rule = childRule;
    tree.at(child).rule = parentRule;
    tree.refreshSubtree(ctx.data, parent, y, ctx.weights);

    double yLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
    double yLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

    alpha = std::exp(yLogPi + yLogL - xLogPi - xLogL);
    alpha = alpha > 1.0 ? 1.0 : alpha;

    if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
      *stepTaken = true;
    } else {
      tree.restoreSubtree(ctx.scratch.snapshot);
    }
  } else {
    // both children share a rule: swap it with the parent's in both
    Rule parentRule = tree.at(parent).rule;
    Rule childRule = tree.at(leftChild).rule;

    tree.at(parent).rule = childRule;
    tree.at(leftChild).rule = parentRule;
    tree.at(rightChild).rule = parentRule;
    bool swapIsSensible = ruleIsValid(ctx, tree, parent, childRule.variableIndex);
    if (childRule.variableIndex != parentRule.variableIndex && swapIsSensible)
      swapIsSensible = ruleIsValid(ctx, tree, parent, parentRule.variableIndex);
    tree.at(parent).rule = parentRule;
    tree.at(leftChild).rule = childRule;
    tree.at(rightChild).rule = childRule;

    if (!swapIsSensible) return -1.0;

    double xLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
    double xLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

    tree.snapshotSubtree(parent, ctx.scratch.snapshot);

    tree.at(parent).rule = childRule;
    tree.at(leftChild).rule = parentRule;
    tree.at(rightChild).rule = parentRule;
    tree.refreshSubtree(ctx.data, parent, y, ctx.weights);

    double yLogPi = ctx.treePrior.treeLogProbability(tree, ctx.data);
    double yLogL = logLikelihoodForBranch(ctx, leaf, tree, parent, y, sigma);

    alpha = std::exp(yLogPi + yLogL - xLogPi - xLogL);
    alpha = alpha > 1.0 ? 1.0 : alpha;

    if (ext_rng_simulateBernoulli(rng, alpha) == 1) {
      *stepTaken = true;
    } else {
      tree.restoreSubtree(ctx.scratch.snapshot);
    }
  }

  return alpha;
}

enum class StepType { birth, death, swap, change };

template <IntegrableLeafModel L>
double metropolisJumpForTree(const MoveContext& ctx, const L& leaf, ext_rng* rng,
                             Tree& tree, const double* y, double sigma,
                             bool* stepTaken, StepType* stepType) {
  double alpha;
  double u = ext_rng_simulateContinuousUniform(rng);

  if (u < ctx.birthOrDeathProbability) {
    bool birthed;
    alpha = birthOrDeathMove(ctx, leaf, rng, tree, y, sigma, stepTaken, &birthed);
    *stepType = birthed ? StepType::birth : StepType::death;
  } else if (u < ctx.birthOrDeathProbability + ctx.swapProbability) {
    alpha = swapMove(ctx, leaf, rng, tree, y, sigma, stepTaken);
    *stepType = StepType::swap;
  } else {
    alpha = changeMove(ctx, leaf, rng, tree, y, sigma, stepTaken);
    *stepType = StepType::change;
  }

  return alpha;
}

}  // namespace bartcore

#endif  // BARTCORE_MOVES_HPP
