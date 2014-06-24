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
  
  void Rule::copyFrom(const Rule& other)
  {
    if (other.variableIndex == BART_INVALID_RULE_VARIABLE) {
      variableIndex = BART_INVALID_RULE_VARIABLE;
      splitIndex    = BART_INVALID_RULE_VARIABLE;
      return;
    }
    
    variableIndex = other.variableIndex;
    splitIndex    = other.splitIndex;
  }
  
  void Rule::swapWith(Rule& other)
  {
    Rule temp(other);
    other = *this;
    *this = temp;
  }

  bool Rule::equals(const Rule& other) const {
    if (variableIndex != other.variableIndex) return false;
    
    // since is a union of variables of the same width, bit-wise equality is sufficient
    return splitIndex == other.splitIndex;
  }
  
  
  void Node::clearObservations()
  {
    if (!isTop()) {
      observationIndices = NULL;
      numObservations = 0;
    }
    if (!isBottom()) {
      leftChild->clearObservations();
      p.rightChild->clearObservations();
    } else {
      m.average = 0.0;
    }
  }
  
  void Node::clear()
  {
    if (!isBottom()) {
      delete leftChild;
      delete p.rightChild;
      
      leftChild = NULL;
      p.rule.invalidate();
    }
    clearObservations();
  }
  
  Node::Node(size_t* observationIndices, size_t numObservations, size_t numPredictors) :
    parent(NULL), leftChild(NULL), enumerationIndex(BART_INVALID_NODE_ENUM), variablesAvailableForSplit(NULL),
    observationIndices(observationIndices), numObservations(numObservations)
  {
    variablesAvailableForSplit = new bool[numPredictors];
    for (size_t i = 0; i < numPredictors; ++i) variablesAvailableForSplit[i] = true;
  }
  
  Node::Node(const Node& parent, size_t numPredictors) :
    parent(const_cast<Node*>(&parent)), leftChild(NULL), enumerationIndex(BART_INVALID_NODE_ENUM),
    variablesAvailableForSplit(NULL), observationIndices(NULL), numObservations(0)
  {
    variablesAvailableForSplit = new bool[numPredictors];
    std::memcpy(variablesAvailableForSplit, parent.variablesAvailableForSplit, sizeof(bool) * numPredictors);
  }
  
  Node::~Node()
  {
    if (leftChild != NULL) {
      delete leftChild; leftChild = NULL;
      delete p.rightChild; p.rightChild = NULL;
    }
    delete [] variablesAvailableForSplit; variablesAvailableForSplit = NULL;
  }
  
  void Node::copyFrom(const BARTFit& fit, const Node& other)
  {
    parent = other.parent;
    leftChild = other.leftChild;
    if (leftChild != NULL) {
      p.rightChild = other.p.rightChild;
      p.rule.copyFrom(other.p.rule);
    }
    else {
      m.average = other.m.average;
      m.numEffectiveObservations = other.m.numEffectiveObservations;
    }
    
    enumerationIndex = other.enumerationIndex;
    std::memcpy(variablesAvailableForSplit, other.variablesAvailableForSplit, sizeof(bool) * fit.data.numPredictors);
    
    observationIndices = other.observationIndices;
    numObservations = other.numObservations;
  }
  
  void Node::print(const BARTFit& fit) const
  {
    size_t depth = getDepth();
        
    for (size_t i = 0; i < depth; ++i) ext_printf("  ");
    
    ext_printf("node:");
    ext_printf(" n: %lu", getNumObservations());
    ext_printf(" TBN: %u%u%u", isTop(), isBottom(), childrenAreBottom());
    ext_printf(" Avail: ");
    
    for (size_t i = 0; i < fit.data.numPredictors; ++i) ext_printf("%u", variablesAvailableForSplit[i]);
    
    if (!isBottom()) {
      ext_printf(" var: %d ", p.rule.variableIndex);
      
      if (fit.data.variableTypes[p.rule.variableIndex] == CATEGORICAL) {
        ext_printf("CATRule: ");
        for (size_t i = 0; 0 < fit.scratch.numCutsPerVariable[p.rule.variableIndex]; ++i) ext_printf(" %u", (p.rule.categoryDirections >> i) & 1);
      } else {
        ext_printf("ORDRule: (%d)=%f", p.rule.splitIndex, p.rule.getSplitValue(fit));
      }
    } else {
      ext_printf(" ave: %f", m.average);
    }
    ext_printf("\n");
    
    if (!isBottom()) {
      leftChild->print(fit);
      p.rightChild->print(fit);
    }
  }
  
  size_t Node::getNumBottomNodes() const
  {
    if (isBottom()) {
      return 1;
    } else {
      return leftChild->getNumBottomNodes() + p.rightChild->getNumBottomNodes();
    }
  }
  
  size_t Node::getNumNotBottomNodes() const
  {
    if (isBottom()) return 0;
    
    return leftChild->getNumNotBottomNodes() + p.rightChild->getNumNotBottomNodes() + 1;
  }
  
  size_t Node::getNumNoGrandNodes() const
  {
    if (isBottom()) return 0;
    if (childrenAreBottom()) return 1;
    return (leftChild->getNumNoGrandNodes() + p.rightChild->getNumNoGrandNodes());
  }
  
  size_t Node::getNumSwappableNodes() const
  {
    if (isBottom() || childrenAreBottom()) return 0;
    if ((leftChild->isBottom()  || leftChild->childrenAreBottom()) &&
        (p.rightChild->isBottom() || p.rightChild->childrenAreBottom())) return 1;
    
    return (leftChild->getNumSwappableNodes() + p.rightChild->getNumSwappableNodes() + 1);
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
    fillBottomVector(*node.p.rightChild, result);
  }
  
  void fillAndEnumerateBottomVector(bart::Node& node, bart::NodeVector& result, size_t& index)
  {
    if (node.isBottom()) {
      result.push_back(&node);
      node.enumerationIndex = index++;
      return;
    }
    
    fillAndEnumerateBottomVector(*node.getLeftChild(), result, index);
    fillAndEnumerateBottomVector(*node.getRightChild(), result, index);
  }
  
  void fillNoGrandVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom()) return;
    if (node.childrenAreBottom()) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }

    fillNoGrandVector(*node.getLeftChild(), result);
    fillNoGrandVector(*node.getRightChild(), result);
  }
  
  void fillNotBottomVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom()) return;
    if (node.childrenAreBottom()) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }
  
    fillNotBottomVector(*node.leftChild, result);
    fillNotBottomVector(*node.p.rightChild, result);
    
    result.push_back(const_cast<bart::Node*>(&node));
  }
  
  void fillSwappableVector(const bart::Node& node, bart::NodeVector& result)
  {
    if (node.isBottom() || node.childrenAreBottom()) return;
    if ((node.leftChild->isBottom()  || node.leftChild->childrenAreBottom()) && 
        (node.p.rightChild->isBottom() || node.p.rightChild->childrenAreBottom())) {
      result.push_back(const_cast<bart::Node*>(&node));
      return;
    }
    
    fillSwappableVector(*node.leftChild, result);
    fillSwappableVector(*node.p.rightChild, result);
    
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
    
    if (p.rule.goesRight(fit, x)) return p.rightChild->findBottomNode(fit, x);
    
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
      if (isTop()) {
        if (fit.data.weights == NULL) {
          m.average = ext_mt_computeMean(fit.threadManager, y, numObservations);
          m.numEffectiveObservations = (double) numObservations;
        } else {
          m.average = ext_mt_computeWeightedMean(fit.threadManager, y, numObservations, fit.data.weights, &m.numEffectiveObservations);
        }
      } else {
        if (fit.data.weights == NULL) {
          m.average = ext_mt_computeIndexedMean(fit.threadManager, y, observationIndices, numObservations);
          m.numEffectiveObservations = (double) numObservations;
        } else {
          m.average = ext_mt_computeIndexedWeightedMean(fit.threadManager, y, observationIndices, numObservations, fit.data.weights, &m.numEffectiveObservations);
        }
      }
      
      return;
    }
    
    leftChild->clearObservations();
    p.rightChild->clearObservations();
    
    /*size_t numThreads, numElementsPerThread;
    ext_mt_getNumThreadsForJob(fit.threadManager, numObservations, MIN_NUM_OBSERVATIONS_IN_NODE_PER_THREAD,
                               &numThreads, &numElementsPerThread); */
    
    size_t numOnLeft;
    IndexOrdering ordering(fit, p.rule);
    
    //if (numThreads <= 1) {
      numOnLeft = (isTop() ?
                   partitionRange(observationIndices, 0, numObservations, ordering) :
                   partitionIndices(observationIndices, numObservations, ordering));
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
      threadData[i].length = numObservations - i * numElementsPerThread;
      threadData[i].ordering = &ordering;
      threadDataPtrs[i] = &threadData[i];
     
      
      
      ext_mt_runTasks(fit.threadManager, &partitionTask, threadDataPtrs, numThreads);
      
      
      
      numOnLeft = mergePartitions(threadData, numThreads);
      
      ext_stackFree(threadDataPtrs);
      ext_stackFree(threadData);
    } */
    
    
    leftChild->observationIndices = observationIndices;
    leftChild->numObservations = numOnLeft;
    p.rightChild->observationIndices = observationIndices + numOnLeft;
    p.rightChild->numObservations = numObservations - numOnLeft;
    
    
    leftChild->addObservationsToChildren(fit, y);
    p.rightChild->addObservationsToChildren(fit, y);
  }
  
  void Node::addObservationsToChildren(const BARTFit& fit) {
    if (isBottom()) {
      m.average = 0.0;
      return;
    }
    
    leftChild->clearObservations();
    p.rightChild->clearObservations();
    
    size_t numOnLeft;
    IndexOrdering ordering(fit, p.rule);
    
    numOnLeft = (isTop() ?
                 partitionRange(observationIndices, 0, numObservations, ordering) :
                 partitionIndices(observationIndices, numObservations, ordering));
    
    
    leftChild->observationIndices = observationIndices;
    leftChild->numObservations = numOnLeft;
    p.rightChild->observationIndices = observationIndices + numOnLeft;
    p.rightChild->numObservations = numObservations - numOnLeft;
    
    
    leftChild->addObservationsToChildren(fit);
    p.rightChild->addObservationsToChildren(fit);
  }
	
  void Node::setAverage(const BARTFit& fit, const double* y)
  {
    leftChild = NULL;
    
    size_t numObservations = getNumObservations();
    
    if (isTop()) {
      if (fit.data.weights == NULL) {
        m.average = ext_mt_computeMean(fit.threadManager, y, numObservations);
        m.numEffectiveObservations = (double) numObservations;
      }
      else m.average = ext_mt_computeWeightedMean(fit.threadManager, y, numObservations, fit.data.weights, &m.numEffectiveObservations);
    } else {
      if (fit.data.weights == NULL) {
        m.average = ext_mt_computeIndexedMean(fit.threadManager, y, observationIndices, numObservations);
        m.numEffectiveObservations = (double) numObservations;
      }
      else m.average = ext_mt_computeIndexedWeightedMean(fit.threadManager, y, observationIndices, numObservations, fit.data.weights, &m.numEffectiveObservations);
    }
  }
  
  void Node::setAverages(const BARTFit& fit, const double* y)
  {
    if (isBottom()) {
      setAverage(fit, y);
      return;
    }
    
    leftChild->setAverages(fit, y);
    p.rightChild->setAverages(fit, y);
  }
  
  double Node::computeVariance(const BARTFit& fit, const double* y) const
  {
    if (isTop()) {
      if (fit.data.weights == NULL) {
        return ext_mt_computeVarianceForKnownMean(fit.threadManager, y, numObservations, getAverage());
      } else {
        return ext_mt_computeWeightedVarianceForKnownMean(fit.threadManager, y, numObservations, fit.data.weights, getAverage());
      }
    } else {
      if (fit.data.weights == NULL) {
        return ext_mt_computeIndexedVarianceForKnownMean(fit.threadManager, y, observationIndices, numObservations, getAverage());
      } else {
        return ext_mt_computeIndexedWeightedVarianceForKnownMean(fit.threadManager, y, observationIndices, numObservations, fit.data.weights, getAverage());
      }
    }
  }
  
  double Node::drawFromPosterior(const EndNodePrior& endNodePrior, double residualVariance) const
  {
    if (getNumObservations() == 0) return 0.0;
      
    return endNodePrior.drawFromPosterior(getAverage(), getNumEffectiveObservations(), residualVariance);
  }
  
  void Node::setPredictions(double* restrict y_hat, double prediction) const restrict
  {
    size_t numObservations = getNumObservations();
    
    if (isTop()) {
      ext_setVectorToConstant(y_hat, numObservations, prediction);
      return;
    }
    
    size_t i = 0;
    size_t lengthMod5 = numObservations % 5;
    
    if (lengthMod5 != 0) {
      for ( ; i < lengthMod5; ++i) y_hat[observationIndices[i]] = prediction;
      if (numObservations < 5) return;
    }
    
    for ( ; i < numObservations; i += 5) {
      y_hat[observationIndices[i]]     = prediction;
      y_hat[observationIndices[i + 1]] = prediction;
      y_hat[observationIndices[i + 2]] = prediction;
      y_hat[observationIndices[i + 3]] = prediction;
      y_hat[observationIndices[i + 4]] = prediction;
    }
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
    return (1 + (size_t) std::max(leftChild->getDepthBelow(), p.rightChild->getDepthBelow()));
  }
  
  size_t Node::getNumNodesBelow() const
  {
    if (isBottom()) return 0;
    return 2 + leftChild->getNumNodesBelow() + p.rightChild->getNumNodesBelow();
  }
  
  size_t Node::getNumVariablesAvailableForSplit(size_t numVariables) const {
    return countTrueValues(variablesAvailableForSplit, numVariables);
  }

  void Node::split(const BARTFit& fit, const Rule& newRule, const double* y, bool exhaustedLeftSplits, bool exhaustedRightSplits) {
    if (newRule.variableIndex < 0) ext_throwError("error in split: rule not set\n");
    
    p.rule = newRule;
    
    leftChild    = new Node(*this, fit.data.numPredictors);
    p.rightChild = new Node(*this, fit.data.numPredictors);
    
    if (exhaustedLeftSplits)     leftChild->variablesAvailableForSplit[p.rule.variableIndex] = false;
    if (exhaustedRightSplits) p.rightChild->variablesAvailableForSplit[p.rule.variableIndex] = false;
    
    addObservationsToChildren(fit, y);
  }

  void Node::orphanChildren() {
    // do this w/o clobbering children pointers until details are nailed down
    double numEffectiveObservations = leftChild->m.numEffectiveObservations + p.rightChild->m.numEffectiveObservations;
    
    double average = leftChild->m.average * (leftChild->m.numEffectiveObservations / numEffectiveObservations) +
                     p.rightChild->m.average * (p.rightChild->m.numEffectiveObservations / numEffectiveObservations);
    
    leftChild = NULL;
    m.average = average;
    m.numEffectiveObservations = numEffectiveObservations;
  }
  
  void Node::countVariableUses(uint32_t* variableCounts) const
  {
    if (isBottom()) return;
    
    ++variableCounts[p.rule.variableIndex];
    
    leftChild->countVariableUses(variableCounts);
    p.rightChild->countVariableUses(variableCounts);
  }
}
