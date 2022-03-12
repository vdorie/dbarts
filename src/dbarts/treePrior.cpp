#include "config.hpp"
#include <dbarts/model.hpp>

#include <cmath>
#include <dbarts/cstdint.hpp>

#include <misc/alloca.h>

#include <external/io.h>
#include <external/random.h>
#include <external/stats.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/data.hpp>
#include <dbarts/types.hpp>
#include "functions.hpp"
#include "node.hpp"
#include "tree.hpp"

using std::uint32_t;
using std::uint64_t;
using std::int32_t;

namespace {
  using namespace dbarts;
  
  double computeTreeLogProbability(const CGMPrior& prior, const BARTFit& fit, const Node& node);
}

namespace dbarts {
  double CGMPrior::computeGrowthProbability(const BARTFit& fit, const Node& node) const
  {
    if (node.getNumVariablesAvailableForSplit(fit.data.numPredictors) == 0) return 0.0;

#ifdef MATCH_BAYES_TREE
    if (node.getNumEffectiveObservations() < 5.0) {
      return 0.001 * base / std::pow(1.0 + node.getDepth(), power);
    }
#endif
    
    return base / std::pow(1.0 + static_cast<double>(node.getDepth()), power);
  }
  
  double CGMPrior::computeTreeLogProbability(const BARTFit& fit, const Tree& tree) const
  {
    return ::computeTreeLogProbability(*this, fit, tree.top);
  }
  
  double CGMPrior::computeSplitVariableLogProbability(const BARTFit& fit, const Node& node) const
  {
    if (splitProbabilities == NULL)
      return -std::log(static_cast<double>(node.getNumVariablesAvailableForSplit(fit.data.numPredictors)));

    double totalProbability = 0.0;
    for (size_t i = 0; i < fit.data.numPredictors; ++i)
      if (node.variablesAvailableForSplit[i]) totalProbability += splitProbabilities[i];

    return std::log(splitProbabilities[node.p.rule.variableIndex] / totalProbability);
  }
  
  double CGMPrior::computeRuleForVariableLogProbability(const BARTFit& fit, const Node& node) const
  {
    int32_t variableIndex = node.p.rule.variableIndex;
    
    double result;
    
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      uint32_t numCategories = fit.numCutsPerVariable[variableIndex];
      
      bool* categoriesCanReachNode = misc_stackAllocate(numCategories, bool);
      
      setCategoryReachability(fit, node, variableIndex, categoriesCanReachNode);
      
      uint32_t numCategoriesCanReachNode = 0;
      for (size_t i = 0; i < numCategories; ++i) if (categoriesCanReachNode[i]) ++numCategoriesCanReachNode;
      
      result  = std::log(std::pow(2.0, static_cast<double>(numCategoriesCanReachNode) - 1.0) - 1.0);
      result -= std::log(std::pow(2.0, static_cast<double>(numCategories - numCategoriesCanReachNode)));
      
      misc_stackFree(categoriesCanReachNode);
    } else {
      int32_t leftCutIndex, rightCutIndex;
      setSplitInterval(fit, node, variableIndex, &leftCutIndex, &rightCutIndex);
      result = -std::log(static_cast<double>(rightCutIndex - leftCutIndex + 1));
    }
    
    return result;
  }
  
  Rule CGMPrior::drawRuleAndVariable(const BARTFit& fit, ext_rng* rng, const Node& node, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const
  {
    int32_t variableIndex = drawSplitVariable(fit, rng, node);
    return drawRuleForVariable(fit, rng, node, variableIndex, exhaustedLeftSplits, exhaustedRightSplits);
  }
  
  int32_t CGMPrior::drawSplitVariable(const BARTFit& fit, ext_rng* rng, const Node& node) const
  {
    if (splitProbabilities == NULL) {
      size_t numGoodVariables = node.getNumVariablesAvailableForSplit(fit.data.numPredictors);
      
      size_t variableNumber = ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, numGoodVariables);
      
      return findIndexOfIthPositiveValue(node.variablesAvailableForSplit, fit.data.numPredictors, variableNumber);
    } else {
      double totalProbability = 0.0;
      for (int32_t i = 0; i < static_cast<int32_t>(fit.data.numPredictors); ++i) {
        if (node.variablesAvailableForSplit[i]) {
          totalProbability += splitProbabilities[i];
        }
      }

      double cutoff = ext_rng_simulateContinuousUniform(rng) * totalProbability;

      double runningProbability = 0.0;
      for (int32_t i = 0; i < static_cast<int32_t>(fit.data.numPredictors); ++i) {
        if (node.variablesAvailableForSplit[i]) {
          runningProbability += splitProbabilities[i];
          if (runningProbability >= cutoff)
            return i;
        }
      }

      ext_throwError("drawSplitVariable went beyond array extent without selecting a variable");
    }
  }
  
  Rule CGMPrior::drawRuleForVariable(const BARTFit& fit, ext_rng* rng, const Node& node, int32_t variableIndex, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const
  {
    Rule result = { DBARTS_INVALID_RULE_VARIABLE, { DBARTS_INVALID_RULE_VARIABLE } };
    
    result.variableIndex = variableIndex;
    
    *exhaustedLeftSplits = false;
    *exhaustedRightSplits = false;
    
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      uint32_t numCategories = fit.numCutsPerVariable[variableIndex];
      
      bool* categoriesCanReachNode = misc_stackAllocate(numCategories, bool);
      // result.categoryDirections = new CategoryBranchingType[numCategories];
      result.categoryDirections = 0l;
      
      setCategoryReachability(fit, node, variableIndex, categoriesCanReachNode);
      
      uint32_t numCategoriesCanReachNode = 0;
      for (uint32_t i = 0; i < numCategories; ++i) if (categoriesCanReachNode[i]) ++numCategoriesCanReachNode;
      if (numCategoriesCanReachNode < 2) {
        misc_stackFree(categoriesCanReachNode);
        ext_throwError("error in TreePrior::drawRule: less than 2 values left for cat var\n");
      }
      
      bool* sendCategoriesRight = misc_stackAllocate(numCategoriesCanReachNode, bool);
      sendCategoriesRight[0] = true; // the first value always goes right so that at least one does
      
      uint64_t categoryIndex = ext_rng_simulateUnsignedIntegerUniformInRange(rng, 0, static_cast<uint64_t>(std::pow(2.0, static_cast<double>(numCategoriesCanReachNode) - 1.0) - 1.0));
      setBinaryRepresentation(numCategoriesCanReachNode - 1, static_cast<uint32_t>(categoryIndex), sendCategoriesRight + 1);
      
      uint32_t sendIndex = 0;
      for (uint32_t i = 0; i < numCategories; ++i) {
        if (categoriesCanReachNode[i]) {
          // result.categoryDirections[i] = sendCategoriesRight[sendIndex++] == true ? BART_CAT_RIGHT : BART_CAT_LEFT;
          if (sendCategoriesRight[sendIndex++]) result.setCategoryGoesRight(i);
          else result.setCategoryGoesLeft(i);
        } else {
          // result.categoryDirections[i] = ext_simulateBernoulli(0.5) == 1 ? BART_CAT_RIGHT : BART_CAT_LEFT;
          if (ext_rng_simulateBernoulli(rng, 0.5) == 1) result.setCategoryGoesRight(i);
          else result.setCategoryGoesLeft(i);
        }
      }
      
      uint32_t numSentCategories = 0;
      for (uint32_t i = 0; i < numCategoriesCanReachNode; ++i) {
        if (sendCategoriesRight[i] == true) ++numSentCategories;
      }
      
      if (numCategoriesCanReachNode - numSentCategories == 1) *exhaustedLeftSplits = true;
      if (numSentCategories == 1) *exhaustedRightSplits = true;
      
      misc_stackFree(sendCategoriesRight);
      misc_stackFree(categoriesCanReachNode);
    } else {
      int32_t leftIndex, rightIndex;
      setSplitInterval(fit, node, variableIndex, &leftIndex, &rightIndex);
      
      int32_t numSplits = rightIndex - leftIndex + 1;
      if (numSplits == 0) {
        ext_printf("error in drawRuleFromPrior: no splits left for ordered var\n");
      }
      
      result.splitIndex = static_cast<int32_t>(ext_rng_simulateIntegerUniformInRange(rng, leftIndex, rightIndex + 1));
      
      if (result.splitIndex == leftIndex) *exhaustedLeftSplits = true;
      if (result.splitIndex == rightIndex) *exhaustedRightSplits = true;
    }
    
    return result;
  } 
}

namespace {
  using namespace dbarts;
  
  double computeTreeLogProbability(const CGMPrior& prior, const BARTFit& fit, const Node& node)
  {
    double probabilityNodeIsNotTerminal = prior.computeGrowthProbability(fit, node);
    
    if (node.isBottom()) return std::log(1.0 - probabilityNodeIsNotTerminal);
    
    double result;
    
    result  = std::log(probabilityNodeIsNotTerminal);
    result += prior.computeSplitVariableLogProbability(fit, node);
    result += prior.computeRuleForVariableLogProbability(fit, node);
    
    result = result + ::computeTreeLogProbability(prior, fit, *node.getLeftChild()) + ::computeTreeLogProbability(prior, fit, *node.getRightChild());
    
    return result;
  }
}
