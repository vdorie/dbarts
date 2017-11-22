#include "config.hpp"
#include "birthDeathRule.hpp"

#include <cstddef> // size_t
#include <cmath>   // exp
#include <cstring> // memcpy

#include <external/alloca.h>
#include <external/linearAlgebra.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>
#include "likelihood.hpp"
#include "node.hpp"
#include "tree.hpp"

using std::size_t;


namespace {
  using namespace dbarts;
  struct State {
    Node node;
    
    void store(const Node& other);
    void destroy();
    void restore(Node& other) const;
  };
}

namespace dbarts {
  
  Node* drawBirthableNode(const BARTFit& fit, ext_rng* rng, const Tree& tree, double* nodeSelectionProbability);
  Node* drawChildrenKillableNode(const BARTFit& fit, ext_rng* rng, const Tree& tree, double* nodeSelectionProbability);
  
  double computeUnnormalizedNodeBirthProbability(const BARTFit& fit, const Node& node);
  double computeProbabilityOfBirthStep(const BARTFit& fit, const Tree& tree); // same as below but that has a step cached
  double computeProbabilityOfBirthStep(const BARTFit& fit, const Tree& tree, bool birthableNodeExists);
  double computeProbabilityOfSelectingNodeForDeath(const Tree& tree);
  double computeProbabilityOfSelectingNodeForBirth(const BARTFit& fit, const Tree& tree);
    
  // returns probability of jump
  double birthOrDeathNode(const BARTFit& fit, size_t chainNum, Tree& tree, const double* y, bool* stepWasTaken, bool* stepWasBirth)
  {
    dbarts::State& state(fit.state[chainNum]);
    
    double ratio;
    
    ::State* oldStatePtr = ext_stackAllocate(1, ::State);
    ::State& oldState(*oldStatePtr);
    
    
    // Rather than flipping a coin to see if birth or death, we have to first check that either is possible.
    // Since that involves pretty much finding a node to give birth, we just do that and then possibly ignore
    // it.

    double transitionProbabilityOfSelectingNodeForBirth;
    Node* nodeToChangePtr = drawBirthableNode(fit, state.rng, tree, &transitionProbabilityOfSelectingNodeForBirth);
    
    double transitionProbabilityOfBirthStep = computeProbabilityOfBirthStep(fit, tree, nodeToChangePtr != NULL);
    
    if (ext_rng_simulateBernoulli(state.rng, transitionProbabilityOfBirthStep) == 1) {
      *stepWasBirth = true;
      
      Node& nodeToChange(*nodeToChangePtr);
      
      double parentPriorGrowthProbability = fit.model.treePrior->computeGrowthProbability(fit, nodeToChange);
      double oldPriorProbability = 1.0 - parentPriorGrowthProbability;
      double oldLogLikelihood = computeLogLikelihoodForBranch(fit, chainNum, nodeToChange, y, state.sigma);
      
      // now perform birth;
      oldState.store(nodeToChange);

      bool exhaustedLeftSplits, exhaustedRightSplits;
      Rule newRule = fit.model.treePrior->drawRuleAndVariable(fit, state.rng, nodeToChange, &exhaustedLeftSplits, &exhaustedRightSplits);
      nodeToChange.split(fit, chainNum, newRule, y, exhaustedLeftSplits, exhaustedRightSplits);
      
      // determine how to go backwards
      double leftPriorGrowthProbability  = fit.model.treePrior->computeGrowthProbability(fit, *nodeToChange.getLeftChild());
      double rightPriorGrowthProbability = fit.model.treePrior->computeGrowthProbability(fit, *nodeToChange.getRightChild());
      double newPriorProbability = parentPriorGrowthProbability * (1.0 - leftPriorGrowthProbability) * (1.0 - rightPriorGrowthProbability);

      double newLogLikelihood = computeLogLikelihoodForBranch(fit, chainNum, nodeToChange, y, state.sigma);

      double transitionProbabilityOfDeathStep = 1.0 - computeProbabilityOfBirthStep(fit, tree);
      double transitionProbabilityOfSelectingNodeForDeath = computeProbabilityOfSelectingNodeForDeath(tree);
      
      // compute ratios
      double priorRatio = newPriorProbability / oldPriorProbability;
      double transitionRatio = (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath) /
                               (transitionProbabilityOfBirthStep * transitionProbabilityOfSelectingNodeForBirth);
      
      double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);
      
      ratio = priorRatio * likelihoodRatio * transitionRatio;
      
      if (ext_rng_simulateContinuousUniform(state.rng) < ratio) {
        oldState.destroy();
        
        *stepWasTaken = true;
      } else {
        oldState.restore(nodeToChange);
        
        *stepWasTaken = false;
      }
    } else {
      *stepWasBirth = false;
      
      double transitionProbabilityOfDeathStep = 1.0 - transitionProbabilityOfBirthStep;
      
      double transitionProbabilityOfSelectingNodeForDeath;
      nodeToChangePtr = drawChildrenKillableNode(fit, state.rng, tree, &transitionProbabilityOfSelectingNodeForDeath);
      
      Node& nodeToChange(*nodeToChangePtr);
      
      double parentPriorGrowthProbability = fit.model.treePrior->computeGrowthProbability(fit, nodeToChange);
      double leftPriorGrowthProbability   = fit.model.treePrior->computeGrowthProbability(fit, *nodeToChange.getLeftChild());
      double rightPriorGrowthProbability  = fit.model.treePrior->computeGrowthProbability(fit, *nodeToChange.getRightChild());
      double oldLogLikelihood = computeLogLikelihoodForBranch(fit, chainNum, nodeToChange, y, state.sigma);
      
      oldState.store(nodeToChange);
      
      // now figure out how the node could have given birth
      nodeToChange.orphanChildren();
      
      double newLogLikelihood = computeLogLikelihoodForBranch(fit, chainNum, nodeToChange, y, state.sigma);
      transitionProbabilityOfBirthStep = computeProbabilityOfBirthStep(fit, tree, true);
#ifdef MATCH_BAYES_TREE
      ext_simulateContinuousUniform();
#endif
      double transitionProbabilityOfSelectingNodeForBirth = computeProbabilityOfSelectingNodeForBirth(fit, tree);
      
      double oldPriorProbability = parentPriorGrowthProbability * (1.0 - leftPriorGrowthProbability) * (1.0 - rightPriorGrowthProbability);
      double newPriorProbability = 1.0 - parentPriorGrowthProbability;
      
      double priorRatio = newPriorProbability / oldPriorProbability;
      double transitionRatio = (transitionProbabilityOfBirthStep * transitionProbabilityOfSelectingNodeForBirth) /
                               (transitionProbabilityOfDeathStep * transitionProbabilityOfSelectingNodeForDeath);
      
      double likelihoodRatio = std::exp(newLogLikelihood - oldLogLikelihood);
      
      ratio = priorRatio * likelihoodRatio * transitionRatio;
      
      if (ext_rng_simulateContinuousUniform(state.rng) < ratio) {
        oldState.destroy();
        
        *stepWasTaken = true;
      } else {
        oldState.restore(nodeToChange);
        
        *stepWasTaken = false;
      }
    }
    
    ext_stackFree(oldStatePtr);
    
    return ratio < 1.0 ? ratio : 1.0;
  }
  
  // transition mechanism
  double computeProbabilityOfBirthStep(const BARTFit& fit, const Tree& tree)
  {
    NodeVector bottomNodes(tree.getBottomNodes());
    size_t numBottomNodes = bottomNodes.size();
    
    bool birthableNodeExists = false;
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      if (computeUnnormalizedNodeBirthProbability(fit, *bottomNodes[i]) > 0.0) {
        birthableNodeExists = true;
        break;
      }
    }
    
#ifdef MATCH_BAYES_TREE
    if (birthableNodeExists) ext_simulateContinuousUniform();
#endif
    
    return computeProbabilityOfBirthStep(fit, tree, birthableNodeExists);
  }
  
  double computeProbabilityOfBirthStep(const BARTFit& fit, const Tree& tree, bool birthableNodeExists)
  {
    if (!birthableNodeExists) return 0.0;
    if (tree.hasSingleNode()) return 1.0;
    
    return fit.model.birthProbability;
  }
  
  double computeProbabilityOfSelectingNodeForDeath(const Tree& tree)
  {
    size_t numNodesWhoseChildrenAreBottom = tree.getNumNodesWhoseChildrenAreBottom();
    if (numNodesWhoseChildrenAreBottom == 0) return 0.0;
    
    return 1.0 / static_cast<double>(numNodesWhoseChildrenAreBottom);
  }
                                                                                                      
  double computeProbabilityOfSelectingNodeForBirth(const BARTFit& fit, const Tree& tree)
  {
    if (tree.hasSingleNode()) return 1.0;
    
    NodeVector bottomNodes(tree.getBottomNodes());
    size_t numBottomNodes = bottomNodes.size();
    
    double totalProbability = 0.0;
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      totalProbability += computeUnnormalizedNodeBirthProbability(fit, *bottomNodes[i]);
    }
    
    if (totalProbability <= 0.0) return 0.0;
    
    return 1.0 / totalProbability;
  }
  
  Node* drawBirthableNode(const BARTFit& fit, ext_rng* rng, const Tree& tree, double* nodeSelectionProbability)
  {
    Node* result = NULL;
    
#ifndef MATCH_BAYES_TREE
    if (tree.hasSingleNode()) {
      *nodeSelectionProbability = 1.0;
      return tree.getTop();
    }
#endif
    
    NodeVector bottomNodes(tree.getBottomNodes());
    size_t numBottomNodes = bottomNodes.size();
    
    double* nodeBirthProbabilities = ext_stackAllocate(numBottomNodes, double);
    double totalProbability = 0.0;
        
    for (size_t i = 0; i < numBottomNodes; ++i) {
      nodeBirthProbabilities[i] = computeUnnormalizedNodeBirthProbability(fit, *bottomNodes[i]);
      totalProbability += nodeBirthProbabilities[i];
    }
    
    if (totalProbability > 0.0) {
      ext_scalarMultiplyVectorInPlace(nodeBirthProbabilities, numBottomNodes, 1.0 / totalProbability);

      size_t index = ext_rng_drawFromDiscreteDistribution(rng, nodeBirthProbabilities, numBottomNodes);

      result = bottomNodes[index];
      *nodeSelectionProbability = nodeBirthProbabilities[index];
    } else {
      *nodeSelectionProbability = 0.0;
    }
    
    ext_stackFree(nodeBirthProbabilities);
    
    return result;
  }
  
  Node* drawChildrenKillableNode(const BARTFit& fit, ext_rng* rng, const Tree& tree, double* nodeSelectionProbability)
  {
    NodeVector nodesWhoseChildrenAreBottom(tree.getNodesWhoseChildrenAreAtBottom());
    size_t numNodesWhoseChildrenAreBottom = nodesWhoseChildrenAreBottom.size();
    
    if (numNodesWhoseChildrenAreBottom == 0) {
      *nodeSelectionProbability = 0.0;
      return NULL;
    }
    
    size_t index = ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, numNodesWhoseChildrenAreBottom);
    *nodeSelectionProbability = 1.0 / static_cast<double>(numNodesWhoseChildrenAreBottom);
    
    return nodesWhoseChildrenAreBottom[index];
  }
  
  double computeUnnormalizedNodeBirthProbability(const BARTFit& fit, const Node& node)
  {
    bool hasVariablesAvailable = node.getNumVariablesAvailableForSplit(fit.data.numPredictors) > 0;
    
    return hasVariablesAvailable ? 1.0 : 0.0;
  }
}

namespace {
  using namespace dbarts;
  
  void ::State::store(const Node& other) {
    std::memcpy(&node, &other, sizeof(Node));
  }
  
  void ::State::destroy() {
    if (node.getLeftChild() != NULL) {
      // successful death step
      delete node.getLeftChild();
      delete node.getRightChild();
    }
  }
  
  void ::State::restore(Node& other) const {
    if (node.getLeftChild() == NULL) {
      // failed birth step
      if (other.getLeftChild() != NULL) {
        // TODO: clean this up
        delete other.leftChild; other.leftChild = NULL;
        delete other.p.rightChild; other.p.rightChild = NULL;
      }
    }
    std::memcpy(&other, &node, sizeof(Node));
  }
}
