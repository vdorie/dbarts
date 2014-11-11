#include "config.hpp"
#include "functions.hpp"

#include <vector>

#include <external/io.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/types.hpp>
#include "birthDeathRule.hpp"
#include "changeRule.hpp"
#include "node.hpp"
#include "swapRule.hpp"

using std::int32_t;
using std::uint32_t;

namespace dbarts {
  void updateCategoricalVariablesAvailable(const BARTFit& fit, Node* node, int32_t variableIndex, bool* catGoesRight);
  void updateOrdinalVariablesAvailable(const BARTFit& fit, Node* node, int32_t variableIndex, int32_t leftIndex, int32_t rightIndex);
  
  void updateVariablesAvailable(const BARTFit& fit, Node& node, int32_t variableIndex)
  {
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      bool* catGoesRight = new bool[fit.scratch.numCutsPerVariable[variableIndex]];
      setCategoryReachability(fit, node, variableIndex, catGoesRight);
      
      // note catGoesRight is deleted in here
      updateCategoricalVariablesAvailable(fit, &node, variableIndex, catGoesRight);
    } else {
      int32_t leftIndex, rightIndex;
      setSplitInterval(fit, node, variableIndex, &leftIndex, &rightIndex);
      updateOrdinalVariablesAvailable(fit, &node, variableIndex, leftIndex, rightIndex);
    }
  }
  
  void updateOrdinalVariablesAvailable(const BARTFit& fit, Node* node, int32_t variableIndex, int32_t leftIndex, int32_t rightIndex)
  {
    int32_t numSplits = rightIndex - leftIndex + 1;
    node->variablesAvailableForSplit[variableIndex] = (numSplits >= 1);
    
    if (!node->isBottom()) {
      int32_t leftChildLeftIndex, leftChildRightIndex, rightChildLeftIndex, rightChildRightIndex;
      leftChildLeftIndex = leftIndex;
      rightChildLeftIndex = leftIndex;
      leftChildRightIndex = rightIndex;
      rightChildRightIndex = rightIndex;
      
      if (node->p.rule.variableIndex == variableIndex) {
        leftChildRightIndex = node->p.rule.splitIndex - 1;
        rightChildLeftIndex = node->p.rule.splitIndex + 1;
      }
      
      updateOrdinalVariablesAvailable(fit, node->getLeftChild(), variableIndex, leftChildLeftIndex, leftChildRightIndex);
      updateOrdinalVariablesAvailable(fit, node->getRightChild(), variableIndex, rightChildLeftIndex, rightChildRightIndex);
    }
  }
  
  void updateCategoricalVariablesAvailable(const BARTFit& fit, Node* node, int32_t variableIndex, bool* catGoesRight)
  {
    uint32_t numCategories = fit.scratch.numCutsPerVariable[variableIndex];


    node->variablesAvailableForSplit[variableIndex] = countTrueValues(catGoesRight, numCategories) >= 2;
        
    if (!node->isBottom()) {
      bool* leftChildCategories  = new bool [numCategories];
      bool* rightChildCategories = new bool [numCategories];
      
      for (size_t i = 0; i < numCategories; ++i) {
        leftChildCategories[i] = catGoesRight[i];
        rightChildCategories[i] = catGoesRight[i];
      }
      
      if (node->p.rule.variableIndex == variableIndex) {
        uint32_t categoryMask = 1u;
        for (size_t i = 0; i < numCategories; ++i) {
          if (catGoesRight[i] == true) {
            if ((node->p.rule.categoryDirections & categoryMask) != 0) {
              leftChildCategories[i] = false;
            } else {
              rightChildCategories[i] = false;
            }
            categoryMask <<= 1;
          }
        }
      }
      
      updateCategoricalVariablesAvailable(fit, node->getLeftChild(), variableIndex, leftChildCategories);
      updateCategoricalVariablesAvailable(fit, node->getRightChild(), variableIndex, rightChildCategories);
    }
    
    delete [] catGoesRight;
  }
  
  double metropolisJumpForTree(const BARTFit& fit, Tree& tree, const double* y,
                               bool* stepTaken, StepType* stepType)
  { 
    double alpha;
    bool birthedTree;
    
    
    
    double u = ext_rng_simulateContinuousUniform(fit.control.rng);
    // ext_printf("type: %s; ", u < fit.model.birthOrDeathProbability ? "birth/death" : (u < fit.model.birthOrDeathProbability + fit.model.swapProbability ? "swap" : "change"));
    if (u < fit.model.birthOrDeathProbability) {
      alpha = birthOrDeathNode(fit, tree, y, stepTaken, &birthedTree);
      if (birthedTree == true) {
        *stepType = BIRTH;
      } else {
        *stepType = DEATH;
      }
    } else if (u < fit.model.birthOrDeathProbability + fit.model.swapProbability) { 
      alpha = swapRule(fit, tree, y, stepTaken);
      *stepType = SWAP;
    } else {
      alpha = changeRule(fit, tree, y, stepTaken);
      *stepType = CHANGE;
    }
    // const char * const jumpNames[] = { "birth", "death", "swap", "change" };
    // ext_printf("jump: %s, succ: %s, u: %f\n", jumpNames[*stepType], *stepTaken ? "true" : "false", u);
    

    return alpha;
  }
  
  void setBinaryRepresentation(uint32_t length, uint32_t ind, bool* d)
  {
    if (length > 64) ext_throwError("attempt to get binary representation for more than 32 categories not supported.");
    for (uint32_t i = 0; i < length; ++i) {
      d[i] = ((ind & 1) == true);
      ind >>= 1;
    }
  }
  
  size_t countTrueValues(bool* v, size_t length)
  {
    size_t result = 0;
    for (size_t i = 0; i < length; ++i) if (v[i] == true) ++result;
    return result;
  }

  void setCategoryReachability(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* categoriesCanReachNode)
  {
    if (fit.data.variableTypes[variableIndex] != CATEGORICAL) ext_throwError("error in setCategoryBranching: not a categorical variable\n");
    
    uint32_t numCategories = fit.scratch.numCutsPerVariable[variableIndex];
    for (uint32_t i = 0; i < numCategories; ++i) categoriesCanReachNode[i] = true;
    
    const Node* curr = &node;
    const Node* parent;
    while (!curr->isTop()) {
      parent = curr->parent;
      if (parent->p.rule.variableIndex == variableIndex) {
        if (curr == parent->getLeftChild()) {
          for (uint32_t i = 0; i < numCategories; ++i) { // parent woulda taken this right, but we're left child
            if (parent->p.rule.categoryGoesRight(i)) categoriesCanReachNode[i] = false;
          }
        } else {
          for (uint32_t i = 0; i < numCategories; ++i) {
            if (parent->p.rule.categoryGoesRight(i) == false) categoriesCanReachNode[i] = false;

          }
        }
      }
      curr = parent;
    }
  }
  
  int32_t findIndexOfIthPositiveValue(bool* values, size_t numValues, size_t i)
  {
    size_t positiveValueCount = 0;
    
    for (uint32_t j = 0; j < numValues; ++j) {
      if (values[j] == true) {
        if (positiveValueCount == i) return static_cast<int32_t>(j);
        ++positiveValueCount;
      }
    }
    
    return DBARTS_INVALID_RULE_VARIABLE;
  }
  
  // get interval of available splits for ordered variable
  //
  // we go up from bottom of tree, and when we find our variable:
  //   if node is a right child, can only split to right of it
  //   if left child, split to left of it
  void setSplitInterval(const BARTFit& fit, const Node& startNode, int32_t variableIndex, int32_t* leftIndex, int32_t* rightIndex)
  {
    if (variableIndex == DBARTS_INVALID_RULE_VARIABLE) {
      ext_throwError("error in getSplitInterval: variable index invalid\n");
    }
    if (fit.data.variableTypes[variableIndex] != ORDINAL) {
      ext_throwError("error in getSplitInterval: variable not ordered\n");
    }
    
    bool leftFound = false;
    bool rightFound = false;
    
    *leftIndex = 0; // left value if you top out
    *rightIndex = static_cast<int32_t>(fit.scratch.numCutsPerVariable[variableIndex]) - 1; // right value if you top out
    
    bool isRightChild;
    
    const Node* curr = &startNode;
    
    // move up tree until you have topped out or found both left and right
    while (!curr->isTop() && !(leftFound && rightFound)) {
      
      if (curr == curr->parent->getRightChild()) {
        isRightChild = true;
      } else {
        isRightChild = false;
      }
      
      curr = curr->parent;
      
      // if you find the variable set the left or right
      if (curr->p.rule.variableIndex == variableIndex) {
        if (isRightChild && !leftFound) {
          leftFound = true;
          *leftIndex = curr->p.rule.splitIndex + 1;
        }
        if (!isRightChild && !rightFound) {
          rightFound = true;
          *rightIndex = curr->p.rule.splitIndex - 1;
        }
      }
    }
  }
}
