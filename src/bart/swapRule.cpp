#include "config.hpp"
#include "swapRule.hpp"

#include <cstddef>
#include <cmath>
#include <bart/cstdint>
#include <cstring>

#include <external/stats.h>
#include <external/io.h>
#include <external/alloca.h>

#include <bart/bartFit.hpp>
#include <bart/model.hpp>
#include <bart/scratch.hpp>
#include <bart/types.hpp>
#include "functions.hpp"
#include "likelihood.hpp"
#include "node.hpp"
#include "tree.hpp"

using std::size_t;
using std::uint32_t;

// note: I got real tired of fixing the unreadable code that went before and haven't managed to
// re-write this yet

namespace {
  // State doesn't handle the sloppiness of swapping with both children, only
  // cares about one. When both are swapped, deals w/left child's rule only.
  struct State {
    bart::Rule parentRule;
    
    double* averages;
    double* numEffectiveObservations;
    
    size_t numNodesInSubtree;
    bool* variablesAvailable;
    
    size_t** observationIndicesPtrs; // duplicates where original was pointing
    size_t* numObservations;  // duplicates length of original
    size_t** observationIndices;     // duplicates content of original
    
    void store(const bart::BARTFit& fit, const bart::Node& node);
    void destroy();
    void restore(const bart::BARTFit& fit, bart::Node& node); // destroys state's internals
  };
}

namespace bart {
  bool categoricalRuleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* catGoesRight);
  bool ordinalRuleIsValid(const Node& node, int32_t variableIndex, int32_t leftIndex, int32_t rightIndex);
  bool ruleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex);
  
  double swapRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken)
  // step which tries swapping rules
  {
    double alpha = 0.0; // note backout = (alpha = -1)
    *stepTaken = false;
    
    NodeVector swappableNodes(tree.getSwappableNodes());
    size_t numSwappableNodes = swappableNodes.size();
    
    //if there are no swappable rule back out
    if (numSwappableNodes == 0) return -1.0;
    
    // randomly choose a node with a swappable rule = parent
    uint32_t nodeIndex = (uint32_t) ext_simulateUnsignedIntegerUniformInRange(0, numSwappableNodes);
    
    Node& parent(*swappableNodes[nodeIndex]);
    Node& leftChild(*parent.getLeftChild());
    Node& rightChild(*parent.getRightChild());
    
    bool leftHasRule = false, rightHasRule = false;
    if ( !leftChild.isBottom() &&  leftChild.p.rule.variableIndex != BART_INVALID_RULE_VARIABLE) leftHasRule  = true;
    if (!rightChild.isBottom() && rightChild.p.rule.variableIndex != BART_INVALID_RULE_VARIABLE) rightHasRule = true;
    
    if (!leftHasRule && !rightHasRule) ext_throwError("error in SwapRule: neither child of parent has a rule\n");
    
    bool childrenHaveSameRule = leftHasRule && rightHasRule && leftChild.p.rule.equals(rightChild.p.rule);
    
    if (!childrenHaveSameRule) {
      //find out which children have rules and pick one
      
      Node* childPtr;
      
      if (leftHasRule && rightHasRule) {
        if (ext_simulateBernoulli(0.5) == 1) {
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
        double XLogL = computeLogLikelihoodForBranch(fit, parent, y);
        
        parent.p.rule.swapWith(child.p.rule);
        
        parent.addObservationsToChildren(fit, y);
        
        //  fix VarAvail
        parentVariableIndex = parent.p.rule.variableIndex;
        childVariableIndex  =  child.p.rule.variableIndex;
        updateVariablesAvailable(fit, parent, parentVariableIndex);
        if (parentVariableIndex != childVariableIndex) updateVariablesAvailable(fit, parent, childVariableIndex);
        
        //get logpri and logL from current tree (X)
        double YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double YLogL = computeLogLikelihoodForBranch(fit, parent, y);
                
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_simulateBernoulli(alpha) == 1) {
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
      // std::memcpy(&oldRightChildRule, &parent.rightChild->rule, sizeof(bart::Rule));
      
      parent.p.rule.swapWith(leftChild.p.rule);
      // temporarily just copy in left rule; give ownership over memory if step not rejected
      rightChild.p.rule = leftChild.p.rule;
      // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(bart::Rule));
      
      //check if rule is ok
      int32_t parentVariableIndex = parent.p.rule.variableIndex;
      int32_t childVariableIndex  = leftChild.p.rule.variableIndex;
      
      bool swapIsSensible = ruleIsValid(fit, parent, parentVariableIndex);
      if (parentVariableIndex != childVariableIndex && swapIsSensible) swapIsSensible = ruleIsValid(fit, parent, childVariableIndex);
      
      if (swapIsSensible) {
        // swap back to calculate probabilities
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = leftChild.p.rule;
        // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(bart::Rule));
        
        ::State oldState;
        oldState.store(fit, parent);
        
        double XLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double XLogL = computeLogLikelihoodForBranch(fit, parent, y);
        
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = leftChild.p.rule;
        // std::memcpy(&parent.rightChild->rule, &parent.leftChild->rule, sizeof(bart::Rule));
        
        parent.addObservationsToChildren(fit, y);
        
        //  fix VarAvail
        childVariableIndex = leftChild.p.rule.variableIndex;
        parentVariableIndex = parent.p.rule.variableIndex;
        updateVariablesAvailable(fit, parent, parentVariableIndex);
        if (parentVariableIndex != childVariableIndex) updateVariablesAvailable(fit, parent, childVariableIndex);
        
        double YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        double YLogL = computeLogLikelihoodForBranch(fit, parent, y);
        
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_simulateBernoulli(alpha) == 1) {
          oldState.destroy();
          // accept, so make right rule copy deep and trash old
          rightChild.p.rule.copyFrom(leftChild.p.rule);
          
          *stepTaken = true;
        } else {
          oldState.restore(fit, parent);
          // reject, so copy back in old right rule
          rightChild.p.rule = oldRightChildRule;
          // std::memcpy(&parent.rightChild->rule, &oldRightChildRule, sizeof(bart::Rule));
          
          *stepTaken = false;
        }
      } else {
        // checkrule failed, swap back
        parent.p.rule.swapWith(leftChild.p.rule);
        rightChild.p.rule = oldRightChildRule;
        // std::memcpy(&parent.rightChild->rule, &oldRightChildRule, sizeof(bart::Rule));
        
        alpha = -1.0;
        *stepTaken = false;
      }
    }
    
    return alpha;
  }
  
  bool categoricalRuleIsValid(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* catGoesRight)
  {
    if (node.isBottom()) return true;
    
    uint32_t numCategories = fit.scratch.numCutsPerVariable[variableIndex];
    
    bool* leftChildCategories  = ext_stackAllocate(numCategories, bool);
    bool* rightChildCategories = ext_stackAllocate(numCategories, bool);
    
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
      ext_stackFree(leftChildCategories);
      ext_stackFree(rightChildCategories);
      
      return false;
    }
    
    if (!categoricalRuleIsValid(fit, *node.getLeftChild(), variableIndex,  leftChildCategories) ||
        !categoricalRuleIsValid(fit, *node.getRightChild(), variableIndex, rightChildCategories))
    {
      ext_stackFree(leftChildCategories);
      ext_stackFree(rightChildCategories);
      
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
      bool* catGoesRight = ext_stackAllocate(fit.scratch.numCutsPerVariable[variableIndex], bool);
      setCategoryReachability(fit, node, variableIndex, catGoesRight);
      
      bool result = categoricalRuleIsValid(fit, node, variableIndex, catGoesRight);
      ext_stackFree(catGoesRight);
      return result;
    }
    
    int32_t leftIndex, rightIndex;
    setSplitInterval(fit, node, variableIndex, &leftIndex, &rightIndex);
    return ordinalRuleIsValid(node, variableIndex, leftIndex, rightIndex);
  }
}

// see comments in changeRule.cpp for what the heck I'm doing here
namespace {
  void storeTree(State& state, const bart::BARTFit& fit, const bart::Node& node, size_t& nodeIndex, size_t& bottomNodeIndex) {
    // copy variables available w/brute force
    std::memcpy(state.variablesAvailable + nodeIndex * fit.data.numPredictors, node.variablesAvailableForSplit, fit.data.numPredictors * sizeof(bool));

    state.observationIndicesPtrs[nodeIndex] = node.observationIndices;
    state.numObservations[nodeIndex] = node.numObservations;
    state.observationIndices[nodeIndex] = new size_t[node.numObservations];
    std::memcpy(state.observationIndices[nodeIndex], (const size_t*) node.observationIndices, node.numObservations * sizeof(size_t));

    ++nodeIndex;
    
    if (node.isBottom()) {
      state.averages[bottomNodeIndex] = node.getAverage();
      state.numEffectiveObservations[bottomNodeIndex++] = node.getNumEffectiveObservations();
      return;
    }
    
    storeTree(state, fit, *node.getLeftChild(), nodeIndex, bottomNodeIndex);
    storeTree(state, fit, *node.getRightChild(), nodeIndex, bottomNodeIndex);
  }
  
  void restoreTree(State& state, const bart::BARTFit& fit, bart::Node& node, size_t& nodeIndex, size_t& bottomNodeIndex) {
    std::memcpy(node.variablesAvailableForSplit, state.variablesAvailable + nodeIndex * fit.data.numPredictors, fit.data.numPredictors * sizeof(bool));

    node.observationIndices = state.observationIndicesPtrs[nodeIndex];
    node.numObservations = state.numObservations[nodeIndex];
    std::memcpy(node.observationIndices, (const size_t*) state.observationIndices[nodeIndex], state.numObservations[nodeIndex] * sizeof(size_t));
    
    ++nodeIndex;
    
    if (node.isBottom()) {
      node.setAverage(state.averages[bottomNodeIndex]);
      node.setNumEffectiveObservations(state.numEffectiveObservations[bottomNodeIndex++]);
      return;
    }
    
    restoreTree(state, fit, *node.getLeftChild(), nodeIndex, bottomNodeIndex);
    restoreTree(state, fit, *node.getRightChild(), nodeIndex, bottomNodeIndex);
  }
  
  void State::store(const bart::BARTFit& fit, const bart::Node& node) {
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
  
  void State::restore(const bart::BARTFit& fit, bart::Node& node) {
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
