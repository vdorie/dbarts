#include "config.hpp"
#include "node.hpp"

#include <cstring>    // memcpy
#include <algorithm>  // int max

#include <external/alloca.h>
#include <external/io.h>
#include <external/linearAlgebra.h>
#include <external/stats.h>
#include <external/stats_mt.h>

#include <bart/bartFit.hpp>
#include <bart/data.hpp>
#include <bart/model.hpp>
#include <bart/scratch.hpp>
#include "functions.hpp"

using std::uint64_t;

namespace bart {
  
  double Rule::getSplitValue(const BARTFit& fit) const
  {
    if (variableIndex < 0) return -1000.0;
    if (fit.data.variableTypes[variableIndex] != ORDINAL) return -2000.0;
    
    return fit.scratch.cutPoints[variableIndex][splitIndex];
  }
  
  Rule::Rule() : variableIndex(BART_INVALID_RULE_VARIABLE), splitIndex(BART_INVALID_RULE_VARIABLE)
  {
  }
  
  Rule::Rule(const Rule& other) : variableIndex(other.variableIndex),
    splitIndex(other.splitIndex)
  {

  }
  
  void Rule::invalidate() {
    variableIndex = BART_INVALID_RULE_VARIABLE;
    splitIndex = BART_INVALID_RULE_VARIABLE;
  }
  
  bool Rule::goesRight(const BARTFit& fit, const double* x) const
  {
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      // x is a double, but that is 64 bits wide, and as such we can treat it as
      // a 64 bit integer
      uint32_t categoryId = (uint32_t) *((const uint64_t*) (x + variableIndex));
      
      return categoryGoesRight(categoryId);
    } else {
      const double* splitValues = fit.scratch.cutPoints[variableIndex];
      
      return x[variableIndex] > splitValues[splitIndex];
    }
  }
  
  /*bool Rule::goesRight(const BARTFit& fit, const double* x) const
  {
    if (variableIndex == BART_INVALID_RULE_VARIABLE) ext_throwError("invalid rule for in Rule::goesRight");
    
    if (fit.variableTypes[variableIndex] == CATEGORICAL) {
      
      size_t numCategories = fit.numCutsPerVariable[splitIndex];

      
      uint32_t categoryId;
      for (categoryId = 0; categoryId < numCategories; ++categoryId) if (x[variableIndex] == categoryValues[categoryId]) break;
      
      if (categoryId == numCategories) ext_throwError("error in Rule::goesRight, no match found for cat variable");
      
      return categoryGoesRight(categoryId);
    } else {
      const double*& splitValues(fit.cutPoints[variableIndex]);
    
      return x[variableIndex] > splitValues[splitIndex];
    }
  } */
  
  void Rule::copyFrom(const BARTFit& fit, const Rule& other)
  {
    if (other.variableIndex == BART_INVALID_RULE_VARIABLE) {
      variableIndex = BART_INVALID_RULE_VARIABLE;
      splitIndex    = BART_INVALID_RULE_VARIABLE;
      return;
    }
    
    variableIndex = other.variableIndex;
    
    if (fit.data.variableTypes[variableIndex] == ORDINAL) {
      splitIndex = other.splitIndex;
    } else {
      categoryDirections = other.categoryDirections;
    }
  }
  
  void Rule::swapWith(Rule& other)
  {
    Rule temp;
    std::memcpy(&temp, (const void*) &other, sizeof(Rule));
    std::memcpy(&other, this,  sizeof(Rule));
    std::memcpy(this, &temp, sizeof(Rule));
  }

  bool Rule::equals(const BARTFit& fit, const Rule& other) const {
    if (variableIndex != other.variableIndex) return false;
    
    if (fit.data.variableTypes[variableIndex] == CATEGORICAL) {
      return categoryDirections == other.categoryDirections;
    }
    
    return splitIndex == other.splitIndex;
  }
  
  
  void Node::clearObservations()
  {
    if (!isTop()) {
      observationIndices = NULL;
      numObservationsInNode = 0;
    }
    if (!isBottom()) {
      leftChild->clearObservations();
      rightChild->clearObservations();
    } else {
      average = 0.0;
    }
  }
  
  void Node::clear()
  {
    if (!isBottom()) {
      delete leftChild;
      delete rightChild;
      
      leftChild = NULL;
    }
    clearObservations();
    rule.invalidate();
  }
  
  Node::Node(size_t* observationIndices, size_t numObservationsInNode, size_t numPredictors) :
    parent(NULL), leftChild(NULL), average(0.0), enumerationIndex(BART_INVALID_NODE_ENUM), variablesAvailableForSplit(NULL),
    observationIndices(observationIndices), numObservationsInNode(numObservationsInNode)
  {
    variablesAvailableForSplit = new bool[numPredictors];
    for (size_t i = 0; i < numPredictors; ++i) variablesAvailableForSplit[i] = true;
  }
  
  Node::Node(const Node& parent, size_t numPredictors) :
    parent(const_cast<Node*>(&parent)), leftChild(NULL), average(0.0), enumerationIndex(BART_INVALID_NODE_ENUM),
    variablesAvailableForSplit(NULL), observationIndices(NULL), numObservationsInNode(0)
  {
    variablesAvailableForSplit = new bool[numPredictors];
    std::memcpy(variablesAvailableForSplit, parent.variablesAvailableForSplit, sizeof(bool) * numPredictors);
  }
  
  Node::~Node()
  {
    if (leftChild != NULL) {
      delete leftChild; leftChild = NULL;
      delete rightChild; rightChild = NULL;
    }
    delete [] variablesAvailableForSplit; variablesAvailableForSplit = NULL;
  }
  
  void Node::copyFrom(const BARTFit& fit, const Node& other)
  {
    parent = other.parent;
    leftChild = other.leftChild;
    if (leftChild != NULL) rightChild = other.rightChild;
    else average = other.average;
    
    rule.copyFrom(fit, other.rule);
    
    enumerationIndex = other.enumerationIndex;
    std::memcpy(variablesAvailableForSplit, other.variablesAvailableForSplit, sizeof(bool) * fit.data.numPredictors);
    
    observationIndices = other.observationIndices;
    numObservationsInNode = other.numObservationsInNode;
  }
  
  void Node::print(const BARTFit& fit) const
  {
    size_t depth = getDepth();
        
    for (size_t i = 0; i < depth; ++i) ext_printf("  ");
    
    ext_printf("node:");
    ext_printf(" n: %lu", getNumObservationsInNode());
    ext_printf(" TBN: %u%u%u", isTop(), isBottom(), childrenAreBottom());
    ext_printf(" Avail: ");
    
    for (size_t i = 0; i < fit.data.numPredictors; ++i) ext_printf("%u", variablesAvailableForSplit[i]);
    
    if (!isBottom()) {
      ext_printf(" var: %d ", rule.variableIndex);
      
      if (fit.data.variableTypes[rule.variableIndex] == CATEGORICAL) {
        ext_printf("CATRule: ");
        for (size_t i = 0; 0 < fit.scratch.numCutsPerVariable[rule.variableIndex]; ++i) ext_printf(" %u", (rule.categoryDirections >> i) & 1);
      } else {
        ext_printf("ORDRule: (%d)=%f", rule.splitIndex, rule.getSplitValue(fit));
      }
    } else {
      ext_printf(" ave: %f", average);
    }
    ext_printf("\n");
    
    if (!isBottom()) {
      leftChild->print(fit);
      rightChild->print(fit);
    }
  }
  
  size_t Node::getNumBottomNodes() const
  {
    if (isBottom()) {
      return 1;
    } else {
      return leftChild->getNumBottomNodes() + rightChild->getNumBottomNodes();
    }
  }
  
  size_t Node::getNumNotBottomNodes() const
  {
    if (isBottom()) return 0;
    
    return leftChild->getNumNotBottomNodes() + rightChild->getNumNotBottomNodes() + 1;
  }
  
  size_t Node::getNumNoGrandNodes() const
  {
    if (isBottom()) return 0;
    if (childrenAreBottom()) return 1;
    return (leftChild->getNumNoGrandNodes() + rightChild->getNumNoGrandNodes());
  }
  
  size_t Node::getNumSwappableNodes() const
  {
    if (isBottom() || childrenAreBottom()) return 0;
    if ((leftChild->isBottom()  || leftChild->childrenAreBottom()) &&
        (rightChild->isBottom() || rightChild->childrenAreBottom())) return 1;
    
    return (leftChild->getNumSwappableNodes() + rightChild->getNumSwappableNodes() + 1);
  }
}

// NOTE: below assumes that walk tree on left
namespace {
  void fillBottomVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom()) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }
    
    fillBottomVector(*node.leftChild, result);
    fillBottomVector(*node.rightChild, result);
  }
  
  void fillAndEnumerateBottomVector(bart::Node& node, bart::NodeVector& result, size_t& index)
  {
    if (node.isBottom()) {
      result.push_back(&node);
      node.enumerationIndex = index++;
      return;
    }
    
    fillAndEnumerateBottomVector(*node.leftChild, result, index);
    fillAndEnumerateBottomVector(*node.rightChild, result, index);
  }
  
  void fillNoGrandVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom()) return;
    if (node.childrenAreBottom()) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }

    fillNoGrandVector(*node.leftChild, result);
    fillNoGrandVector(*node.rightChild, result);
  }
  
  void fillNotBottomVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom()) return;
    if (node.childrenAreBottom()) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }
  
    fillNotBottomVector(*node.leftChild, result);
    fillNotBottomVector(*node.rightChild, result);
    
    result.push_back(const_cast<bart::Node*>(&node));
  }
  
  void fillSwappableVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom() || node.childrenAreBottom()) return;
    if ((node.leftChild->isBottom()  || node.leftChild->childrenAreBottom()) && 
        (node.rightChild->isBottom() || node.rightChild->childrenAreBottom())) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }
    
    fillSwappableVector(*node.leftChild, result);
    fillSwappableVector(*node.rightChild, result);
    
    result.push_back(const_cast<bart::Node*>(&node));
  }
}
namespace bart {  
  NodeVector Node::getBottomVector() const
  {
    NodeVector result;
    fillBottomVector(*this, result);
    return result;
  }
  
  NodeVector Node::getAndEnumerateBottomVector()
  {
    size_t index = 0;
    NodeVector result;
    fillAndEnumerateBottomVector(*this, result, index);
    return result;
  }
  
  NodeVector Node::getNoGrandVector() const
  {
    NodeVector result;
    fillNoGrandVector(*this, result);
    return result;
  }
  
  NodeVector Node::getNotBottomVector() const
  {
    NodeVector result;
    fillNotBottomVector(*this, result);
    return result;
  }
  
  NodeVector Node::getSwappableVector() const
  {
    NodeVector result;
    fillSwappableVector(*this, result);
    return result;
  }
  
  Node* Node::findBottomNode(const BARTFit& fit, const double *x) const
  {
    if (isBottom()) return const_cast<Node*>(this);
    
    if (rule.goesRight(fit, x)) return rightChild->findBottomNode(fit, x);
    
    return leftChild->findBottomNode(fit, x);
  }
}


namespace {
  using namespace bart;
  
  struct IndexOrdering {
    const BARTFit& fit;
    const Rule &rule;
    
    IndexOrdering(const BARTFit& fit, const Rule &rule) : fit(fit), rule(rule) { }
    
    bool operator()(size_t i) const { return rule.goesRight(fit, fit.scratch.Xt + i * fit.data.numPredictors); }
  };
  
  // returns how many observations are on the "left"
  size_t partitionRange(size_t* restrict indices, size_t startIndex, size_t length, IndexOrdering& restrict indexGoesRight) {
    size_t lengthOfLeft;
    
    size_t lh = 0, rh = length - 1;
    size_t i = startIndex;
    while (lh <= rh && rh > 0) {
      if (indexGoesRight(i)) {
        indices[rh] = i;
        i = startIndex + rh--;
      } else {
        indices[lh] = i;
        i = startIndex + ++lh;
      }
    }
    if (lh == 0 && rh == 0) { // ugliness w/wrapping around at 0 makes an off-by-one when all obs go right
      indices[startIndex] = i;
      if (indexGoesRight(i)) {
        lengthOfLeft = 0;
      } else {
        lengthOfLeft = 1;
      }
    } else {
      lengthOfLeft = lh;
    }
    return lengthOfLeft;
  }
  
  size_t partitionIndices(size_t* restrict indices, size_t length, IndexOrdering& restrict indexGoesRight) {
    if (length == 0) return 0;
    
    size_t lengthOfLeft;
    
    size_t lh = 0, rh = length - 1;
    while (lh <= rh && rh > 0) {
      if (indexGoesRight(indices[lh])) {
        size_t temp = indices[rh];
        indices[rh] = indices[lh];
        indices[lh] = temp;
        --rh;
      } else {
        ++lh;
      }
    }
    if (lh == 0 && rh == 0) {
      if (indexGoesRight(indices[0])) {
        lengthOfLeft = 0;
      } else {
        lengthOfLeft = 1;
      }
    } else {
      lengthOfLeft = lh;
    }
    
    return lengthOfLeft;
  }
  
  /*
   // http://en.wikipedia.org/wiki/XOR_swap_algorithm
   void ext_swapVectors(size_t* restrict x, size_t* restrict y, size_t length)
   {
   if (length == 0) return;
   
   size_t lengthMod5 = length % 5;
   
   if (lengthMod5 != 0) {
   for (size_t i = 0; i < lengthMod5; ++i) {
   x[i] ^= y[i];
   y[i] ^= x[i];
   x[i] ^= y[i];
   }
   if (length < 5) return;
   }
   
   for (size_t i = lengthMod5; i < length; i += 5) {
   x[i    ] ^= y[i    ]; y[i    ] ^= x[i    ]; x[i    ] ^= y[i    ];
   x[i + 1] ^= y[i + 1]; y[i + 1] ^= x[i + 1]; x[i + 1] ^= y[i + 1];
   x[i + 2] ^= y[i + 2]; y[i + 2] ^= x[i + 2]; x[i + 2] ^= y[i + 2];
   x[i + 3] ^= y[i + 3]; y[i + 3] ^= x[i + 3]; x[i + 3] ^= y[i + 3];
   x[i + 4] ^= y[i + 4]; y[i + 4] ^= x[i + 4]; x[i + 4] ^= y[i + 4];
   }
   }
   
   // merges adjacent partitions of the form:
   // [ l1 r1 l2 r2 ]
   size_t mergeAdjacentPartitions(size_t* array, size_t firstTotalLength, size_t firstLeftLength,
   size_t secondLeftLength)
   {
   // size_t* l1 = array;
   size_t* r1 = array + firstLeftLength;
   size_t* l2 = array + firstTotalLength;
   // size_t* r2 = array + firstTotalLength + secondLeftLength;
   
   size_t firstRightLength = firstTotalLength - firstLeftLength;
   
   if (secondLeftLength <= firstRightLength) {
   ext_swapVectors(r1, l2, secondLeftLength);
   // end up w/[ l1 l2 r1_2 r1_1 r2 ]
   } else {
   ext_swapVectors(r1, l2 + (secondLeftLength - firstRightLength), firstRightLength);
   // end up w/[ l1 l2_2 l2_1 r1 r2 ]
   }
   
   return firstLeftLength + secondLeftLength;
   }
  
  struct PartitionThreadData {
    size_t* indices;
    size_t startIndex;
    size_t length;
    IndexOrdering* ordering;
    size_t numOnLeft;
  };
  
  size_t mergePartitions(PartitionThreadData* data, size_t numThreads)
  {
    while (numThreads > 1) {
      if (numThreads % 2 == 1) {
        // if odd number, merge last two
        PartitionThreadData* left = &data[numThreads - 2];
        PartitionThreadData* right = &data[numThreads - 1];
        
        left->numOnLeft = mergeAdjacentPartitions(left->indices, left->length, left->numOnLeft, right->numOnLeft);
        left->length += right->length;
        
        --numThreads;
      }
        
      for (size_t i = 0; i < numThreads / 2; ++i) {
        PartitionThreadData* left = &data[2 * i];
        PartitionThreadData* right = &data[2 * i + 1];
        
        left->numOnLeft = mergeAdjacentPartitions(left->indices, left->length, left->numOnLeft, right->numOnLeft);
        left->length += right->length;
        
        // now shift down in array so that valid stuffs always occupy the beginning
        if (i > 0) {
          right = &data[i];
          std::memcpy(right, (const PartitionThreadData*) left, sizeof(PartitionThreadData));
        }
      }
      numThreads /= 2;
    }
    return data[0].numOnLeft;
  }
  
  void partitionTask(void* v_data) {
    PartitionThreadData& data(*static_cast<PartitionThreadData*>(v_data));

    data.numOnLeft = (data.startIndex != ((size_t) -1) ? 
                      partitionRange(data.indices, data.startIndex, data.length, *data.ordering) :
                      partitionIndices(data.indices, data.length, *data.ordering));
  } */
} // anon namespace
// MT not worth it for this, apparently
// #define MIN_NUM_OBSERVATIONS_IN_NODE_PER_THREAD 5000

namespace bart {
  void Node::addObservationsToChildren(const BARTFit& fit, const double* y) {
    if (isBottom()) {
      if (isTop())
        average = ext_mt_computeMean(fit.threadManager, y, numObservationsInNode);
      else
        average = ext_mt_computeIndexedMean(fit.threadManager, y, observationIndices, numObservationsInNode);
      return;
    }
    
    leftChild->clearObservations();
    rightChild->clearObservations();
    
    /*size_t numThreads, numElementsPerThread;
    ext_mt_getNumThreadsForJob(fit.threadManager, numObservationsInNode, MIN_NUM_OBSERVATIONS_IN_NODE_PER_THREAD,
                               &numThreads, &numElementsPerThread); */
    
    size_t numOnLeft;
    IndexOrdering ordering(fit, rule);
    
    //if (numThreads <= 1) {
      numOnLeft = (isTop() ?
                   partitionRange(observationIndices, 0, numObservationsInNode, ordering) :
                   partitionIndices(observationIndices, numObservationsInNode, ordering));
    /*} else {
      PartitionThreadData* threadData = ext_stackAllocate(numThreads, PartitionThreadData);
      void** threadDataPtrs = ext_stackAllocate(numThreads, void*);
      
      size_t i;
      for (i = 0; i < numThreads - 1; ++i) {
        threadData[i].indices = observationIndices + i * numElementsPerThread;
        threadData[i].startIndex = isTop() ? i * numElementsPerThread : ((size_t) -1);
        threadData[i].length = numElementsPerThread;
        threadData[i].ordering = &ordering;
        threadDataPtrs[i] = &threadData[i];
      }
      threadData[i].indices = observationIndices + i * numElementsPerThread;
      threadData[i].startIndex = isTop() ? i * numElementsPerThread : ((size_t) -1);
      threadData[i].length = numObservationsInNode - i * numElementsPerThread;
      threadData[i].ordering = &ordering;
      threadDataPtrs[i] = &threadData[i];
     
      
      
      ext_mt_runTasks(fit.threadManager, &partitionTask, threadDataPtrs, numThreads);
      
      
      
      numOnLeft = mergePartitions(threadData, numThreads);
      
      ext_stackFree(threadDataPtrs);
      ext_stackFree(threadData);
    } */
    
    
    leftChild->observationIndices = observationIndices;
    leftChild->numObservationsInNode = numOnLeft;
    rightChild->observationIndices = observationIndices + numOnLeft;
    rightChild->numObservationsInNode = numObservationsInNode - numOnLeft;
    
    
    leftChild->addObservationsToChildren(fit, y);
    rightChild->addObservationsToChildren(fit, y);
  }
  
  void Node::addObservationsToChildren(const BARTFit& fit) {
    if (isBottom()) {
      average = 0.0;
      return;
    }
    
    leftChild->clearObservations();
    rightChild->clearObservations();
    
    size_t numOnLeft;
    IndexOrdering ordering(fit, rule);
    
    numOnLeft = (isTop() ?
                 partitionRange(observationIndices, 0, numObservationsInNode, ordering) :
                 partitionIndices(observationIndices, numObservationsInNode, ordering));
    
    
    leftChild->observationIndices = observationIndices;
    leftChild->numObservationsInNode = numOnLeft;
    rightChild->observationIndices = observationIndices + numOnLeft;
    rightChild->numObservationsInNode = numObservationsInNode - numOnLeft;
    
    
    leftChild->addObservationsToChildren(fit);
    rightChild->addObservationsToChildren(fit);
  }
	
  void Node::setAverage(const BARTFit& fit, const double* y)
  {
    leftChild = NULL;
    
    if (isTop()) {
      average = ext_mt_computeMean(fit.threadManager, y, numObservationsInNode);
    } else {
      average = ext_mt_computeIndexedMean(fit.threadManager, y, observationIndices, numObservationsInNode);
    }
  }
  
  void Node::setAverages(const BARTFit& fit, const double* y)
  {
    if (isBottom()) {
      setAverage(fit, y);
      return;
    }
    
    leftChild->setAverages(fit, y);
    rightChild->setAverages(fit, y);
  }
  
  double Node::computeVariance(const BARTFit& fit, const double* y) const
  {
    if (isTop()) {
      return ext_mt_computeVarianceForKnownMean(fit.threadManager, y, numObservationsInNode, average);
    } else {
      return ext_mt_computeIndexedVarianceForKnownMean(fit.threadManager, y, observationIndices, numObservationsInNode, average);
    }
  }
  
  double Node::drawFromPosterior(const EndNodePrior& endNodePrior, double residualVariance) const
  {
    size_t numObservationsInNode = getNumObservationsInNode();
    
    if (numObservationsInNode == 0) return 0.0;
      
    return endNodePrior.drawFromPosterior(getAverage(), numObservationsInNode, residualVariance);
  }
  
  void Node::setPredictions(double* y_hat, double prediction) const
  {
    if (isTop()) {
      ext_setVectorToConstant(y_hat, numObservationsInNode, prediction);
      return;
    }
    
    for (size_t i = 0; i < numObservationsInNode; ++i) y_hat[observationIndices[i]] = prediction;
  }
  
  size_t Node::getDepth() const
  {
    size_t result = 0;
    const Node* node = this;

    while (!node->isTop()) {
      ++result;
      node = node->parent;
    }
    
    return result;
  }
  
  size_t Node::getDepthBelow() const
  {
    if (childrenAreBottom()) return 1;
    if (isBottom()) return 0;
    return (1 + (size_t) std::max(leftChild->getDepthBelow(), rightChild->getDepthBelow()));
  }
  
  size_t Node::getNumNodesBelow() const
  {
    if (isBottom()) return 0;
    return 2 + leftChild->getNumNodesBelow() + rightChild->getNumNodesBelow();
  }
  
  size_t Node::getNumVariablesAvailableForSplit(size_t numVariables) const {
    return countTrueValues(variablesAvailableForSplit, numVariables);
  }

  void Node::split(const BARTFit& fit, const Rule& newRule, const double* y, bool exhaustedLeftSplits, bool exhaustedRightSplits) {
    if (newRule.variableIndex < 0) ext_throwError("error in split: rule not set\n");
    
    std::memcpy(&rule, &newRule, sizeof(Rule));
    
    leftChild  = new Node(*this, fit.data.numPredictors);
    rightChild = new Node(*this, fit.data.numPredictors);
    
    if (exhaustedLeftSplits)   leftChild->variablesAvailableForSplit[rule.variableIndex] = false;
    if (exhaustedRightSplits) rightChild->variablesAvailableForSplit[rule.variableIndex] = false;
    
    addObservationsToChildren(fit, y);
  }

  void Node::orphanChildren() {
    double numObservationsInNode = (double) getNumObservationsInNode();
    double mu = leftChild->getAverage() * ((double) leftChild->getNumObservationsInNode() / numObservationsInNode) +
                rightChild->getAverage() * ((double) rightChild->getNumObservationsInNode() / numObservationsInNode);
    leftChild = NULL;
    setAverage(mu);
    
    rule.invalidate();
  }
  
  void Node::countVariableUses(uint32_t* variableCounts) const
  {
    if (isBottom()) return;
    
    ++variableCounts[rule.variableIndex];
    
    leftChild->countVariableUses(variableCounts);
    rightChild->countVariableUses(variableCounts);
  }
}
