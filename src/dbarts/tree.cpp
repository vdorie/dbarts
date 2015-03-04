#include "config.hpp"
#include "tree.hpp"

#include <cstring>
#include <cstdio>

#include <external/alloca.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/state.hpp>

namespace {
  using namespace dbarts;
  
  // multithread me!
  size_t* createObservationToNodeIndexMap(const BARTFit& fit, const Node& top,
                                          const double* Xt, size_t numObservations)
  {
    if (numObservations == 0) return NULL;
    
    size_t* map = new size_t[numObservations];
        
    for (size_t i = 0; i < numObservations; ++i) {
      const Node* bottomNode = top.findBottomNode(fit, Xt + i * fit.data.numPredictors);
      
      map[i] = bottomNode->enumerationIndex;
    }
    
    return map;
  }
}

namespace dbarts {
  void Tree::setNodeAverages(const BARTFit& fit, const double* y) {
    NodeVector bottomNodes(getBottomNodes());
    
    size_t numBottomNodes = bottomNodes.size();
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      bottomNodes[i]->setAverage(fit, y);
    }
  }
  
  void Tree::sampleAveragesAndSetFits(const BARTFit& fit, double* trainingFits, double* testFits)
  {
    NodeVector bottomNodes(top.getAndEnumerateBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    double* nodePosteriorPredictions = NULL;
    
    if (testFits != NULL) nodePosteriorPredictions = ext_stackAllocate(numBottomNodes, double);
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      const Node& bottomNode(*bottomNodes[i]);
      
      double posteriorPrediction = bottomNode.drawFromPosterior(fit.control.rng, *fit.model.muPrior, fit.state.sigma * fit.state.sigma);
      bottomNode.setPredictions(trainingFits, posteriorPrediction);
      
      if (testFits != NULL) nodePosteriorPredictions[i] = posteriorPrediction;
    }
    
    if (testFits != NULL) {
      size_t* observationNodeMap = createObservationToNodeIndexMap(fit, top, fit.scratch.Xt_test, fit.data.numTestObservations);
      for (size_t i = 0; i < fit.data.numTestObservations; ++i) testFits[i] = nodePosteriorPredictions[observationNodeMap[i]];
      delete [] observationNodeMap;
      
      ext_stackFree(nodePosteriorPredictions);
    }
  }
  
  double* Tree::recoverAveragesFromFits(const BARTFit&, const double* treeFits)
  {
    NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    double* result = new double[numBottomNodes];
    for (size_t i = 0; i < numBottomNodes; ++i) {
      if (bottomNodes[i]->isTop()) {
        result[i] = treeFits[0];
      } else if (bottomNodes[i]->getNumObservations() > 0) {
        result[i] = treeFits[bottomNodes[i]->observationIndices[0]];
      } else {
        result[i] = 0.0;
      }
    }
    
    return(result);
  }
  
  void Tree::setCurrentFitsFromAverages(const BARTFit& fit, const double* posteriorPredictions, double* trainingFits, double* testFits)
  {
    NodeVector bottomNodes(top.getAndEnumerateBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    if (trainingFits != NULL) {
      for (size_t i = 0; i < numBottomNodes; ++i) {
        const Node& bottomNode(*bottomNodes[i]);
        
        bottomNode.setPredictions(trainingFits, posteriorPredictions[i]);
      }
    }
    
    if (testFits != NULL) {
      size_t* observationNodeMap = createObservationToNodeIndexMap(fit, top, fit.scratch.Xt_test, fit.data.numTestObservations);
      for (size_t i = 0; i < fit.data.numTestObservations; ++i) testFits[i] = posteriorPredictions[observationNodeMap[i]];
      delete [] observationNodeMap;
    }
  }
}

namespace {
  using namespace dbarts;
  void mapCutPoints(Node& n, const BARTFit& fit, const double* const* oldCutPoints);
  void collapseEmptyNodes(Node& n, const BARTFit& fit, double* posteriorPredictions);
}

namespace dbarts {
  void Tree::mapOldCutPointsOntoNew(const BARTFit& fit, const double* const* oldCutPoints)
  {
    mapCutPoints(top, fit, oldCutPoints);
  }
  
  void Tree::collapseEmptyNodes(const BARTFit&fit, double* posteriorPredictions)
  {
    top.enumerateBottomNodes();
    ::collapseEmptyNodes(top, fit, posteriorPredictions);
    
    NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    for (size_t i = 0; i < numBottomNodes; ++i)
      posteriorPredictions[i] = posteriorPredictions[bottomNodes[i]->enumerationIndex];
  }
  
  
  void Tree::countVariableUses(uint32_t* variableCounts) {
    top.countVariableUses(variableCounts);
  }
  
  bool Tree::isValid() const {
    const NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    for (size_t j = 0; j < numBottomNodes; ++j) {
      if (bottomNodes[j]->getNumObservations() == 0) return false;
    }
    
    return true;
  }
}

namespace {
  using namespace dbarts;
  
  void mapCutPoints(Node& n, const BARTFit& fit, const double* const* oldCutPoints)
  {
    if (n.isBottom() || n.p.rule.variableIndex == DBARTS_INVALID_RULE_VARIABLE) return;
    
    int32_t varIndex = n.p.rule.variableIndex;
    
    if (fit.data.variableTypes[varIndex] == ORDINAL) {
      int32_t numCuts = static_cast<int32_t>(fit.scratch.numCutsPerVariable[varIndex]);
      double oldCut = oldCutPoints[varIndex][n.p.rule.splitIndex];
      const double* cutPoints_i = fit.scratch.cutPoints[varIndex];
      
      if (numCuts == 0) {
        n.p.rule.invalidate();
        return;
      } else {
        int32_t firstLessThan = n.p.rule.splitIndex;
        // if it starts out below, move it above
        while (firstLessThan < numCuts && cutPoints_i[firstLessThan] < oldCut) ++firstLessThan;
        // now nudge it back down
        if (firstLessThan < numCuts) while (firstLessThan >= 0 && cutPoints_i[firstLessThan] >= oldCut) --firstLessThan;
  
        int32_t newIndex;
        if (firstLessThan >= numCuts - 1) newIndex = numCuts - 1;
        else if (firstLessThan < 0) newIndex = 0;
        else if (cutPoints_i[firstLessThan + 1] == oldCut) newIndex = firstLessThan + 1;
        else if (oldCut - cutPoints_i[firstLessThan] < cutPoints_i[firstLessThan + 1] - oldCut) newIndex = firstLessThan;
        else newIndex = firstLessThan + 1;
        
        n.p.rule.splitIndex = newIndex;
      }        
    }
    
    mapCutPoints(*n.leftChild, fit, oldCutPoints);
    mapCutPoints(*n.p.rightChild, fit, oldCutPoints);
  }
  
  void collapseEmptyNodes(Node& n, const BARTFit& fit, double* posteriorPredictions)
  {
    if (n.isBottom()) return; // only happens if is top and bottom
    
    if (n.getLeftChild()->getNumObservations() == 0 || n.getRightChild()->getNumObservations() == 0) {
      const NodeVector bottomNodes(n.getBottomVector());
      size_t numBottomNodes = bottomNodes.size();
      double* weights = ext_stackAllocate(numBottomNodes, double);
      double* params  = ext_stackAllocate(numBottomNodes, double);
      
      for (size_t i = 0; i < numBottomNodes; ++i) {
        Node& bottomNode(*bottomNodes[i]);
        weights[i] = fit.data.weights == NULL ? (double) bottomNode.getNumObservations() : ext_sumIndexedVectorElements(fit.data.weights, bottomNode.observationIndices, bottomNode.getNumObservations());
        params[i] = posteriorPredictions[bottomNodes[i]->enumerationIndex];
      }
      size_t leftMostEnumerationIndex = bottomNodes[0]->enumerationIndex;
      delete n.getLeftChild();
      delete n.getRightChild();
      
      posteriorPredictions[leftMostEnumerationIndex] = ext_computeWeightedMean(params, numBottomNodes, weights, NULL);
      n.enumerationIndex = leftMostEnumerationIndex;
      
      ext_stackFree(params);
      ext_stackFree(weights);
    } else {
      if (!n.getLeftChild()->isBottom()) collapseEmptyNodes(*n.getLeftChild(), fit, posteriorPredictions);
      if (!n.getRightChild()->isBottom()) collapseEmptyNodes(*n.getRightChild(), fit, posteriorPredictions);
    }
  }
}
