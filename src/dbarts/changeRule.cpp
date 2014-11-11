#include "config.hpp"
#include "changeRule.hpp"

#include <dbarts/cstdint.hpp>
#include <cmath>
#include <cstddef>
#include <cstring>

#include <algorithm>

#include <external/alloca.h>
#include <external/io.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/types.hpp>
#include "functions.hpp"
#include "likelihood.hpp"
#include "node.hpp"
#include "tree.hpp"

using std::size_t;
using std::int32_t;
using std::uint32_t;
using std::uint64_t;

// note: I got real tired of fixing the unreadable code that went before and haven't managed to
// re-write this yet

namespace {
  using namespace dbarts;
  // Since a change rule, well, changes the rule at a node it does not influence the 
  // tree structure. Thus, we store the innards of nodes of the tree flattened by
  // imposing a particular tree-walking order.
  struct State {
    Rule rule;
    
    double* averages;
    double* numEffectiveObservations;
    
    size_t numNodesInSubtree;
    bool* variablesAvailable;

    size_t** observationIndicesPtrs; // duplicates where original was pointing
    size_t* numObservations;         // duplicates length of original
    size_t** observationIndices;     // duplicates content of original
    
    void store(const BARTFit& fit, const Node& node);
    void destroy();
    void restore(const BARTFit& fit, Node& node); // invalidates afterwards
  };
}

namespace dbarts {
  void findReachableBottomNodesForCategory(const Node* curr, int32_t variableIndex, size_t categoryIndex, NodeVector& bottomVector, bool* nodesAreReachable);
  void findGoodOrdinalRules(const BARTFit& fit, const Node& node, int32_t variableIndex, int32_t* lowerIndex, int32_t* upperIndex);
  void findOrdinalMinMaxSplitIndices(const BARTFit& fit, const Node& node, int32_t variableIndex, int32_t* min, int32_t* max);
  void findGoodCategoricalRules(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* categoryCombinationsAreGood, uint32_t* firstGoodCategory);
  bool allTrue(bool* v, size_t length);
  size_t getIndexOfFirstTrueValue(bool* v, size_t length);
  
  
  double changeRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken)
  // step which tries changing the rule 
  {
    double XLogPi, XLogL, YLogPi, YLogL;
    
    double alpha;
    
    *stepTaken = false;
    
    NodeVector notBottomNodes(tree.getNotBottomNodes());
    size_t numNotBottomNodes = notBottomNodes.size();
    
    // get list of nodes with rule = nodes which are not bottom
    if (numNotBottomNodes == 0) return -1.0;
    
    // randomly choose a notBottom node = nodeToChange
    //u=ran1(&idum);
    size_t nodeIndex = ext_rng_simulateUnsignedIntegerUniformInRange(fit.control.rng, 0, numNotBottomNodes);
    Node& nodeToChange(*notBottomNodes[nodeIndex]);
    
    //given the node, choose a new variable for the new rule
    int32_t newVariableIndex = fit.model.treePrior->drawSplitVariable(fit, nodeToChange);
    
    if (fit.data.variableTypes[newVariableIndex] == CATEGORICAL) {
      // get the list of good cat rules given var choice
      uint32_t firstGoodCategory;
      
      uint32_t numCategories = fit.scratch.numCutsPerVariable[newVariableIndex];
      size_t numCategoryCombinations = (1 << (numCategories - 1)) - 1;
      // bool* categoryCombinationsAreGood = new bool[numCategoryCombinations];
      bool* categoryCombinationsAreGood = ext_stackAllocate(numCategoryCombinations, bool);
      
      findGoodCategoricalRules(fit, nodeToChange, newVariableIndex, categoryCombinationsAreGood, &firstGoodCategory);
      uint64_t numGoodRules = countTrueValues(categoryCombinationsAreGood, numCategoryCombinations);
      
      //if there are any good cat rules
      if (numGoodRules > 0) {
        // draw the rule from list of good ones
        //u=ran1(&idum);
        uint32_t goodCategoryNumber = static_cast<uint32_t>(ext_rng_simulateUnsignedIntegerUniformInRange(fit.control.rng, 0, numGoodRules));
        uint32_t categoryCombinationNumber = static_cast<uint32_t>(findIndexOfIthPositiveValue(categoryCombinationsAreGood, numCategoryCombinations, goodCategoryNumber));
        
        
        //get logpri and logL from current tree (X)
        XLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        XLogL = computeLogLikelihoodForBranch(fit, nodeToChange, y);
        
        // copy old rule
        ::State oldState;
        oldState.store(fit, nodeToChange);
        
        // change rule at nodeToChange to the new one
        bool* sel = ext_stackAllocate(numCategories - 1, bool);
        setBinaryRepresentation(numCategories - 1, categoryCombinationNumber, sel);
        
        nodeToChange.p.rule.variableIndex = newVariableIndex;
        nodeToChange.p.rule.categoryDirections = 0u;
        for (size_t j = 0; j < firstGoodCategory; ++j) {
          if (sel[j] == true) nodeToChange.p.rule.setCategoryGoesRight(static_cast<uint32_t>(j));
          else nodeToChange.p.rule.setCategoryGoesLeft(static_cast<uint32_t>(j));
        }
        nodeToChange.p.rule.setCategoryGoesRight(firstGoodCategory);
        for (size_t j = firstGoodCategory + 1; j < numCategories; ++j) {
          if (sel[j - 1] == true) nodeToChange.p.rule.setCategoryGoesRight(static_cast<uint32_t>(j));
          else nodeToChange.p.rule.setCategoryGoesLeft(static_cast<uint32_t>(j));
        }
        
        // fix data at nodes below nodeToChange given new rule
        nodeToChange.addObservationsToChildren(fit, y);
        
        // fix VarAvail
        updateVariablesAvailable(fit, nodeToChange, newVariableIndex);
        if (newVariableIndex != oldState.rule.variableIndex) updateVariablesAvailable(fit, nodeToChange, oldState.rule.variableIndex);
        
        //get logpri and logL from candidate tree (Y)
        YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        YLogL = computeLogLikelihoodForBranch(fit, nodeToChange, y);
        
        //draw go nogo
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_rng_simulateBernoulli(fit.control.rng, alpha) == 1) {
          oldState.destroy();
          
          *stepTaken = true;
        } else {
          oldState.restore(fit, nodeToChange);
          
          *stepTaken = false;
        }
        
        ext_stackFree(sel);
      } else {
        // if no rules for that var abort step
        alpha = -1.0;
      }
      
      ext_stackFree(categoryCombinationsAreGood);
    } else {
      
      //ORD variable
      
      // get the set of good rules = [l,r]
      int32_t left, right;
      findGoodOrdinalRules(fit, nodeToChange, newVariableIndex, &left, &right);
      int32_t numSplitVariables = right - left + 1;
      
      // if there are any rules
      if (numSplitVariables > 0) {
        int32_t newRuleIndex = static_cast<int32_t>(ext_rng_simulateIntegerUniformInRange(fit.control.rng, left, right + 1));
        
        //get logpri and logL from current tree (X)
        XLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        XLogL  = computeLogLikelihoodForBranch(fit, nodeToChange, y);
        
        // copy old rule
        ::State oldState;
        oldState.store(fit, nodeToChange);
        
        // change rule at nodeToChange to the new one
        nodeToChange.p.rule.variableIndex = newVariableIndex;
        nodeToChange.p.rule.splitIndex    = newRuleIndex;
        
        nodeToChange.addObservationsToChildren(fit, y);
        
        updateVariablesAvailable(fit, nodeToChange, newVariableIndex);
        if (newVariableIndex != oldState.rule.variableIndex) updateVariablesAvailable(fit, nodeToChange, oldState.rule.variableIndex);
        
        //get logpri and logL from candidate tree (Y)
        YLogPi = fit.model.treePrior->computeTreeLogProbability(fit, tree);
        YLogL  = computeLogLikelihoodForBranch(fit, nodeToChange, y);
        
        alpha = std::exp(YLogPi + YLogL - XLogPi - XLogL);
        alpha = (alpha > 1.0 ? 1.0 : alpha);
        
        if (ext_rng_simulateBernoulli(fit.control.rng, alpha) == 1) {	
          oldState.destroy();
          *stepTaken = true;
        } else {
          oldState.restore(fit, nodeToChange);
          
          *stepTaken = false;
        }
      } else {
        // if no rules for that var abort step
        alpha = -1.0;
      }
    }
    return alpha; // note -1 means backed out
  }
  
  
  void findReachableBottomNodesForCategory(const Node* curr, int32_t variableIndex, size_t categoryIndex, NodeVector& bottomVector, bool* nodesAreReachable)
  //adds 1 to fcount i if category cat associated with variable VarI
  // can get from node curr to the ith bottom node
  // 
  {
    if (curr->isBottom()) {
      size_t index = 0;
      while (curr != bottomVector[index]) ++index;
      nodesAreReachable[index] = true;
    } else {
      if (curr->p.rule.variableIndex == variableIndex) {
//        if (curr->rule.categoryDirections[categoryIndex] == BART_CAT_RIGHT) {
        if (curr->p.rule.categoryGoesRight(static_cast<uint32_t>(categoryIndex))) {
          findReachableBottomNodesForCategory(curr->getRightChild(), variableIndex, categoryIndex, bottomVector, nodesAreReachable);
        } else {
          findReachableBottomNodesForCategory(curr->getLeftChild(), variableIndex, categoryIndex, bottomVector, nodesAreReachable);
        }
      } else {
        findReachableBottomNodesForCategory(curr->getRightChild(), variableIndex, categoryIndex, bottomVector, nodesAreReachable);
        findReachableBottomNodesForCategory(curr->getLeftChild(), variableIndex, categoryIndex, bottomVector, nodesAreReachable);
      }
    }
  }
  
  void findGoodOrdinalRules(const BARTFit& fit, const Node& node, int32_t variableIndex, int32_t* lowerIndex, int32_t* upperIndex)
  //good rule have splits in [l,u]
  {
    int32_t leftIndex, rightIndex; 
    leftIndex = 0; // left value if you top out
    rightIndex = static_cast<int32_t>(fit.scratch.numCutsPerVariable[variableIndex]) - 1; // right value if you top out
    
    int32_t leftMin, leftMax, rightMin, rightMax;
    
    setSplitInterval(fit, node, variableIndex, &leftIndex, &rightIndex);
    
    leftMin  = rightIndex + 1;
    rightMin = rightIndex + 1;
    leftMax  = leftIndex - 1;
    rightMax = leftIndex - 1;
    
    findOrdinalMinMaxSplitIndices(fit, *node.getLeftChild(), variableIndex, &leftMin, &leftMax);
    findOrdinalMinMaxSplitIndices(fit, *node.getRightChild(), variableIndex, &rightMin, &rightMax);
    
    *lowerIndex = std::max(leftIndex, leftMax + 1);
    *upperIndex = std::min(rightIndex, rightMin - 1);
  }
  
  void findOrdinalMinMaxSplitIndices(const BARTFit& fit, const Node& node, int32_t variableIndex, int32_t* min, int32_t* max)
  //used to find good ord rule, go down tree adjusting min and max whenever VarI is used
  {
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) ext_throwError("error in findOrdinalMinMaxSplitIndices, called on CATEGORICAL var");
    
    if (!node.isBottom()) {
      if (variableIndex == node.p.rule.variableIndex) {
        if (node.p.rule.splitIndex < *min) *min = node.p.rule.splitIndex;
        if (node.p.rule.splitIndex > *max) *max = node.p.rule.splitIndex;
      }
      findOrdinalMinMaxSplitIndices(fit, *node.getLeftChild(), variableIndex, min, max);
      findOrdinalMinMaxSplitIndices(fit, *node.getRightChild(), variableIndex, min, max);
    }
  }
  
  void findGoodCategoricalRules(const BARTFit& fit, const Node& node, int32_t variableIndex, bool* categoryCombinationsAreGood, uint32_t* firstGoodCategory)
  // finds out which categorical rule using VarI are good.
  // a good rule is one that does not result in logically empty bottom nodes
  //n: the node at which the rule is to be set
  //VarI: the variable ~ the rule
  //RuleInd: integer vector whose length is the number of possible rules 2^(NR-1) - 1
  //	on exit 1 if rule ok 0 otherwise, already allocated
  //firstone: first category still "alive" at node n, depends on tree above n
  {
    uint32_t numCategories = fit.scratch.numCutsPerVariable[variableIndex];
    bool* sel = ext_stackAllocate(numCategories, bool);
    
    bool* categoriesGoRight = ext_stackAllocate(numCategories, bool);
    setCategoryReachability(fit, node, variableIndex, categoriesGoRight);
    
    *firstGoodCategory = static_cast<uint32_t>(getIndexOfFirstTrueValue(categoriesGoRight, numCategories));
    
    if (*firstGoodCategory == numCategories) ext_printf("error in findGoodCategoricalRule: no available categories\n");
    
    sel[*firstGoodCategory] = true;
    
    bool* sel1 = ext_stackAllocate(numCategories - 1, bool);
    
    NodeVector leftBottomVector(node.getLeftChild()->getBottomVector());
    size_t numLeftBottomNodes = leftBottomVector.size();
    bool* leftNodesAreReachable = ext_stackAllocate(numLeftBottomNodes, bool);
    
    NodeVector rightBottomVector(node.getRightChild()->getBottomVector());
    size_t numRightBottomNodes = rightBottomVector.size();
    bool* rightNodesAreReachable = ext_stackAllocate(numRightBottomNodes, bool);
    
    
    // 2^(numCategories - 1) - 1
    size_t numCategoryCombinations = (1 << (numCategories - 1)) - 1;
    for (size_t i = 0; i < numCategoryCombinations; ++i) categoryCombinationsAreGood[i] = false;
    
    // enumerate through all possible category branches and test them for validity
    for (uint64_t i = 0; i < numCategoryCombinations; ++i) {
      setBinaryRepresentation(numCategories - 1, static_cast<uint32_t>(i), sel1);
      
      for (size_t j = 0; j < *firstGoodCategory; ++j) sel[j] = sel1[j];
      for (size_t j = *firstGoodCategory + 1; j < numCategories; ++j) sel[j] = sel1[j - 1];
      
      for (size_t j = 0; j < numLeftBottomNodes; ++j) leftNodesAreReachable[j] = false;
      for (size_t j = 0; j < numRightBottomNodes; ++j) rightNodesAreReachable[j] = true;
      
      for (size_t j = 0; j < numCategories; ++j) {
        if (categoriesGoRight[j] == true) {
          if (sel[j] == true) {
            findReachableBottomNodesForCategory(node.getRightChild(), variableIndex, j, rightBottomVector, rightNodesAreReachable);
          } else {
            findReachableBottomNodesForCategory(node.getLeftChild(), variableIndex, j, leftBottomVector, leftNodesAreReachable);
          }
        }
        if (allTrue(leftNodesAreReachable, numLeftBottomNodes) &&
            allTrue(rightNodesAreReachable, numRightBottomNodes))
        {
          categoryCombinationsAreGood[i] = true;
          break;
        }
      }
    }
    
    ext_stackFree(sel);
    ext_stackFree(sel1);
    ext_stackFree(categoriesGoRight);
    ext_stackFree(leftNodesAreReachable);
    ext_stackFree(rightNodesAreReachable);
  }
  
  bool allTrue(bool* v, size_t length) {
    for (size_t i = 0; i < length; ++i) {
      if (v[i] == false) return false;
    }
    return true;
  }
  
  // result will equal length if can't find any
  size_t getIndexOfFirstTrueValue(bool* v, size_t length)
  {
    size_t i;
    for (i = 0; i < length; ++i) {
      if (v[i] == true) break;
    }
    
    return i;
  }
}

// This is a bit of a mess since I want to use copy constructors as often as
// possible but I have an array(ish) of IndexVectors to dupe. Rather than
// allocate an array of IndexVectors locally which then call default constructors
// and then use assignment, instead I placement new each one copying from the
// appropriate place.
//
// Also, it is a bit of a mess since I flatten the tree below.

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
    rule = node.p.rule;
    
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
  
  void ::State::destroy() {
    delete [] averages;
    delete [] numEffectiveObservations;
    
    delete [] variablesAvailable;
    
    delete [] observationIndicesPtrs;
    delete [] numObservations;
    for (size_t i = 0; i < numNodesInSubtree; ++i) delete [] observationIndices[i];
    delete [] observationIndices;
  }
  
  void ::State::restore(const BARTFit& fit, Node& node) {
    node.p.rule = rule;
        
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
