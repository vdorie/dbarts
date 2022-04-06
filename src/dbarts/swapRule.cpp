#include "config.hpp"
#include "swapRule.hpp"

#include <cstddef>
#include <cmath>
#include <dbarts/cstdint.hpp>
#include <cstring>

#include <misc/alloca.h>
#include <external/io.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>
#include <dbarts/types.hpp>
#include "functions.hpp"
#include "likelihood.hpp"
#include "node.hpp"
#include "tree.hpp"

using std::size_t;
using std::uint32_t;

// note: I got real tired of fixing the unreadable code that went before and haven't managed to
// re-write this yet

namespace {
  using namespace dbarts;
  
  // State doesn't handle the sloppiness of swapping with both children, only
  // cares about one. When both are swapped, deals w/left child's rule only.
  struct State {
    Rule parentRule;
    
    double* averages;
    double* numEffectiveObservations;
    
    size_t numNodesInSubtree;
    bool* variablesAvailable;
    
    size_t** observationIndicesPtrs; // duplicates where original was pointing
    size_t* numObservations;  // duplicates length of original
    size_t** observationIndices;     // duplicates content of original
    
    void store(const BARTFit& fit, const Node& node);
    void destroy();
    void restore(const BARTFit& fit, Node& node); // destroys state's internals
  };
}

namespace dbarts {
  bool categoricalRuleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* catGoesRight);
  bool ordinalRuleIsValid(const Node& node, int32_t variableIndex, int32_t leftIndex, int32_t rightIndex);
  bool ruleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex);
  
  double swapRule(const BARTFit& fit, size_t chainNum, Tree& tree, const double* y, double sigma, bool* stepTaken)
  // step which tries swapping rules
  {
    State& state(fit.state[chainNum]);
    
    double alpha = 0.0; // note backout = (alpha = -1)
    *stepTaken = false;
    
    NodeVector& swappableNodes(fit.chainScratch[chainNum].nodeVector);
    swappableNodes.clear();
    tree.fillSwappableNodesVector(swappableNodes);
    size_t numSwappableNodes = swappableNodes.size();
    
    //if there are no swappable rule back out
    if (numSwappableNodes == 0) return -1.0;
    
    // randomly choose a node with a swappable rule = parent
    uint32_t nodeIndex = static_cast<uint32_t>(ext_rng_simulateUnsignedIntegerUniformInRange(state.rng, 0, numSwappableNodes));
    
    Node& parent(*swappableNodes[nodeIndex]);
    Node& leftChild(*parent.getLeftChild());
    Node& rightChild(*parent.getRightChild());
    
    bool leftHasRule = false, rightHasRule = false;
    if ( !leftChild.isBottom() &&  leftChild.p.rule.variableIndex != DBARTS_INVALID_RULE_VARIABLE) leftHasRule  = true;
    if (!rightChild.isBottom() && rightChild.p.rule.variableIndex != DBARTS_INVALID_RULE_VARIABLE) rightHasRule = true;
    
    if (!leftHasRule && !rightHasRule) ext_throwError("error in SwapRule: neither child of parent has a rule\n");
    
    bool childrenHaveSameRule = leftHasRule && rightHasRule && leftChild.p.rule.equals(rightChild.p.rule);
    
    if (!childrenHaveSameRule) {
      //find out which children have rules and pick one
      
      Node* childPtr;
      
      if (leftHasRule && rightHasRule) {
        if (ext_rng_simulateBernoulli(state.rng, 0.5) == 1) {
          childPtr = &leftChild;
        } else {
          childPtr = &rightChild;
        }
      } else if (leftHasRule) {
        childPtr = &leftChild;
      } else {
        childPtr = &rightChild;
      }
      
      Node& child(*childPtr);
      
      // swap rules between parent and child and test that no conflicts arise from doing so
      parent.p.rule.swapWith(child.p.rule);
      
      int32_t parentVariableIndex = parent.p.rule.variableIndex;
      int32_t childVariableIndex  = child.p.rule.variableIndex;
      
      bool swapIsSensible = ruleIsValid(fit, parent, parentVariableIndex);
      if (parentVariableIndex != childVariableIndex && swapIsSensible) swapIsSensible = ruleIsValid(fit, parent, childVariableIndex);
      
      // swap back to calculate probabilities
      parent.p.rule.swapWith(child.p.rule);
      
      //if the swap was ok (rules made sense)
      if (swapIsSensible) {
        ::State oldState;
        oldState.store(fit, parent);
        
        double XLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double XLogL = computeLogLikelihoodForBranch(fit, chainNum, parent, y, sigma);
        
        parent.p.rule.swapWith(child.p.rule);
        
        parent.addObservationsToChildren(fit, chainNum, y);
        
        //  fix VarAvail
        parentVariableIndex = parent.p.rule.variableIndex;
        childVariableIndex  =  child.p.rule.variableIndex;
        updateVariablesAvailable(fit, parent, parentVariableIndex);
        if (parentVariableIndex != childVariableIndex) updateVariablesAvailable(fit, parent, childVariableIndex);
        
        //get logpri and logL from current tree (X)
        double YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double YLogL = computeLogLikelihoodForBranch(fit, chainNum, parent, y, sigma);
                
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_rng_simulateBernoulli(state.rng, alpha) == 1) {
          oldState.destroy();
          
          *stepTaken = true;
        } else {
          oldState.restore(fit, parent);
        }
      } else {
        alpha = -1.0; //not a legal swap	
      }
    } else {
      Rule oldRightChildRule = rightChild.p.rule;
      // std::memcpy(&oldRightChildRule, &parent.rightChild->rule, sizeof(Rule));
      
      parent.p.rule.swapWith(leftChild.p.rule);
      // temporarily just copy in left rule; give ownership over memory if step not rejected
      rightChild.p.rule = leftChild.p.rule;
      // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(Rule));
      
      //check if rule is ok
      int32_t parentVariableIndex = parent.p.rule.variableIndex;
      int32_t childVariableIndex  = leftChild.p.rule.variableIndex;
      
      bool swapIsSensible = ruleIsValid(fit, parent, parentVariableIndex);
      if (parentVariableIndex != childVariableIndex && swapIsSensible) swapIsSensible = ruleIsValid(fit, parent, childVariableIndex);
      
      if (swapIsSensible) {
        // swap back to calculate probabilities
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = leftChild.p.rule;
        // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(Rule));
        
        ::State oldState;
        oldState.store(fit, parent);
        
        double XLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double XLogL = computeLogLikelihoodForBranch(fit, chainNum, parent, y, sigma);
        
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = leftChild.p.rule;
        // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(Rule));
        
        parent.addObservationsToChildren(fit, chainNum, y);
        
        //  fix VarAvail
        childVariableIndex = leftChild.p.rule.variableIndex;
        parentVariableIndex = parent.p.rule.variableIndex;
        updateVariablesAvailable(fit, parent, parentVariableIndex);
        if (parentVariableIndex != childVariableIndex) updateVariablesAvailable(fit, parent, childVariableIndex);
        
        double YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double YLogL = computeLogLikelihoodForBranch(fit, chainNum, parent, y, sigma);
        
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_rng_simulateBernoulli(state.rng, alpha) == 1) {
          oldState.destroy();
          // accept, so make right rule copy deep and trash old
          rightChild.p.rule.copyFrom(leftChild.p.rule);
          
          *stepTaken = true;
        } else {
          oldState.restore(fit, parent);
          // reject, so copy back in old right rule
          rightChild.p.rule = oldRightChildRule;
          // std::memcpy(&parent.rightChild->rule, &oldRightChildRule, sizeof(Rule));
          
          *stepTaken = false;
        }
      } else {
        // checkrule failed, swap back
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = oldRightChildRule;
        // std::memcpy(&parent.rightChild->rule, &oldRightChildRule, sizeof(Rule));
        
        alpha = -1.0;
        *stepTaken = false;
      }
    }
    
    return alpha;
  }
  
  bool categoricalRuleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* catGoesRight)
  {
    if (node.isBottom()) return true;
    
    uint32_t numCategories = fit.numCutsPerVariable[variableIndex];
    
    bool* leftChildCategories  = misc_stackAllocate(numCategories, bool);
    bool* rightChildCategories = misc_stackAllocate(numCategories, bool);
    
    for(uint32_t i = 0; i < numCategories; ++i) {
      leftChildCategories[i]  = catGoesRight[i];
      rightChildCategories[i] = catGoesRight[i];
    }
    
    if (node.p.rule.variableIndex == variableIndex) {
      for (uint32_t i = 0; i < numCategories; ++i) {
        if (catGoesRight[i] == true) {
          if (node.p.rule.categoryGoesRight(i)) {
            leftChildCategories[i] = false;
          } else {
            rightChildCategories[i] = false;
          }
        }
      }
    }
    
    if (countTrueValues( leftChildCategories, numCategories) == 0 ||
        countTrueValues(rightChildCategories, numCategories) == 0)
    {
      misc_stackFree(leftChildCategories);
      misc_stackFree(rightChildCategories);
      
      return false;
    }
    
    if (!categoricalRuleIsValid(fit, *node.getLeftChild(), variableIndex,  leftChildCategories) ||
        !categoricalRuleIsValid(fit, *node.getRightChild(), variableIndex, rightChildCategories))
    {
      misc_stackFree(leftChildCategories);
      misc_stackFree(rightChildCategories);
      
      return false;
    }
    
    return true;
  }
  
  
  bool ordinalRuleIsValid(const Node& node, int32_t variableIndex, int32_t leftIndex, int32_t rightIndex)
  {
    if (node.isBottom()) return true;

    int32_t ruleVariableIndex = node.p.rule.variableIndex;

    if (ruleVariableIndex == variableIndex) {
      int32_t splitIndex = node.p.rule.splitIndex;
      
      if (splitIndex < leftIndex || splitIndex > rightIndex) return false;
      
      
       if (!ordinalRuleIsValid( *node.getLeftChild(), variableIndex, leftIndex, splitIndex - 1) ||
           !ordinalRuleIsValid(*node.getRightChild(), variableIndex, splitIndex + 1, rightIndex))
       {
         return false;
       }
      
      return true;
    }
    
    if (!ordinalRuleIsValid( *node.getLeftChild(), variableIndex, leftIndex, rightIndex) ||
        !ordinalRuleIsValid(*node.getRightChild(), variableIndex, leftIndex, rightIndex))
    {
      return false;
    }
    return true;
  }
  
  bool ruleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex)
  //starting at node n, check rules using VarI to see if they make sense
  {
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      bool* catGoesRight = misc_stackAllocate(fit.numCutsPerVariable[variableIndex], bool);
      setCategoryReachability(fit, node, variableIndex, catGoesRight);
      
      bool result = categoricalRuleIsValid(fit, node, variableIndex, catGoesRight);
      misc_stackFree(catGoesRight);
      return result;
    }
    
    int32_t leftIndex, rightIndex;
    setSplitInterval(fit, node, variableIndex, &leftIndex, &rightIndex);
    return ordinalRuleIsValid(node, variableIndex, leftIndex, rightIndex);
  }
}

// see comments in changeRule.cpp for what the heck I'm doing here
namespace {
  using namespace dbarts;
  
  void storeTree(::State& state, const BARTFit& fit, const Node& node, size_t& nodeIndex, size_t& bottomNodeIndex) {
    // copy variables available w/brute force
    std::memcpy(state.variablesAvailable + nodeIndex * fit.data.numPredictors, node.variablesAvailableForSplit, fit.data.numPredictors * sizeof(bool));

    state.observationIndicesPtrs[nodeIndex] = node.observationIndices;
    state.numObservations[nodeIndex] = node.numObservations;
    state.observationIndices[nodeIndex] = new size_t[node.numObservations];
    std::memcpy(state.observationIndices[nodeIndex], const_cast<const size_t*>(node.observationIndices), node.numObservations * sizeof(size_t));

    ++nodeIndex;
    
    if (node.isBottom()) {
      state.averages[bottomNodeIndex] = node.getAverage();
      state.numEffectiveObservations[bottomNodeIndex++] = node.getNumEffectiveObservations();
      return;
    }
    
    storeTree(state, fit, *node.getLeftChild(), nodeIndex, bottomNodeIndex);
    storeTree(state, fit, *node.getRightChild(), nodeIndex, bottomNodeIndex);
  }
  
  void restoreTree(::State& state, const BARTFit& fit, Node& node, size_t& nodeIndex, size_t& bottomNodeIndex) {
    std::memcpy(node.variablesAvailableForSplit, state.variablesAvailable + nodeIndex * fit.data.numPredictors, fit.data.numPredictors * sizeof(bool));

    node.observationIndices = state.observationIndicesPtrs[nodeIndex];
    node.numObservations = state.numObservations[nodeIndex];
    std::memcpy(node.observationIndices, const_cast<const size_t*>(state.observationIndices[nodeIndex]), state.numObservations[nodeIndex] * sizeof(size_t));
    
    ++nodeIndex;
    
    if (node.isBottom()) {
      node.setAverage(state.averages[bottomNodeIndex]);
      node.setNumEffectiveObservations(state.numEffectiveObservations[bottomNodeIndex++]);
      return;
    }
    
    restoreTree(state, fit, *node.getLeftChild(), nodeIndex, bottomNodeIndex);
    restoreTree(state, fit, *node.getRightChild(), nodeIndex, bottomNodeIndex);
  }
  
  void ::State::store(const BARTFit& fit, const Node& node) {
    parentRule = node.p.rule;
                
    size_t numBottomNodes = node.getNumBottomNodes();
    
    averages = new double[numBottomNodes];
    numEffectiveObservations = new double[numBottomNodes];
    
    numNodesInSubtree = 1 + node.getNumNodesBelow();
    variablesAvailable = new bool[numNodesInSubtree * fit.data.numPredictors];
    
    observationIndicesPtrs = new size_t*[numNodesInSubtree];
    numObservations = new size_t[numNodesInSubtree];
    observationIndices = new size_t*[numNodesInSubtree];
    
    size_t nodeIndex = 0, bottomNodeIndex = 0;
    storeTree(*this, fit, node, nodeIndex, bottomNodeIndex);
  }
  
  void State::destroy() {
    delete [] averages;
    delete [] numEffectiveObservations;
    
    delete [] variablesAvailable;
    
    delete [] observationIndicesPtrs;
    delete [] numObservations;
    for (size_t i = 0; i < numNodesInSubtree; ++i) delete [] observationIndices[i];
    delete [] observationIndices;
  }
  
  void State::restore(const BARTFit& fit, Node& node) {
    bool leftWasSwapped = parentRule.equals(node.getLeftChild()->p.rule);
    node.p.rule.swapWith(leftWasSwapped ? node.getLeftChild()->p.rule : node.getRightChild()->p.rule);
    
    size_t nodeIndex = 0, bottomNodeIndex = 0;
    restoreTree(*this, fit, node, nodeIndex, bottomNodeIndex);
    
    delete [] averages;
    delete [] numEffectiveObservations;
    
    delete [] variablesAvailable;
    
    delete [] observationIndicesPtrs;
    delete [] numObservations;
    for (size_t i = 0; i < numNodesInSubtree; ++i) delete [] observationIndices[i];
    delete [] observationIndices;
  }
}
