#ifndef DBARTS_TREE_HPP
#define DBARTS_TREE_HPP

#include <cstddef>
#include <dbarts/cstdint.hpp>

#include "node.hpp"

namespace dbarts {
  using std::size_t;
  using std::uint32_t;
  
  struct BARTFit;
  
  struct Tree {
    Node top;
    
    Tree(size_t* indices, size_t numObservations, size_t numPredictors) : top(indices, numObservations, numPredictors) { }
    
    void sampleAveragesAndSetFits(const BARTFit& fit, double* trainingFits, double* testFits);
    double* recoverAveragesFromFits(const BARTFit& fit, const double* treeFits); // allocates response; are ordered as bottom nodes are
    void setCurrentFitsFromAverages(const BARTFit& fit, const double* posteriorPredictions, double* trainingFits, double* testFits);
    
    void mapOldCutPointsOntoNew(const BARTFit& fit, const double* const* oldCutPoints, double* posteriorPredictions);
    void collapseEmptyNodes(const BARTFit& fit, double* posteriorPredictions);
    
    void sampleFromPrior(const BARTFit& fit);
    
    Node* getTop() const;
    bool hasSingleNode() const;
    
    size_t getNumBottomNodes() const;
    size_t getNumNotBottomNodes() const;
    size_t getNumNodesWhoseChildrenAreBottom() const;
    size_t getNumSwappableNodes() const;
    
    NodeVector getBottomNodes() const;
    NodeVector getNotBottomNodes() const;
    NodeVector getNodesWhoseChildrenAreAtBottom() const;
    NodeVector getSwappableNodes() const;
    
    void setNodeAverages(const BARTFit& fit, const double* y);
    
    void countVariableUses(uint32_t* variableCounts);
    
    const char* createString() const;
    
    bool isValid() const;
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  inline Node* Tree::getTop() const { return const_cast<Node*>(&top); }
  inline bool Tree::hasSingleNode() const { return top.isBottom(); }
  
  inline size_t Tree::getNumBottomNodes() const { return top.getNumBottomNodes(); }
  inline size_t Tree::getNumNotBottomNodes() const { return top.getNumNotBottomNodes(); }
  inline size_t Tree::getNumNodesWhoseChildrenAreBottom() const { return top.getNumNoGrandNodes(); }
  inline size_t Tree::getNumSwappableNodes() const { return top.getNumSwappableNodes(); }
  
  inline NodeVector Tree::getBottomNodes() const { return top.getBottomVector(); }
  inline NodeVector Tree::getNotBottomNodes() const { return top.getNotBottomVector(); }
  inline NodeVector Tree::getNodesWhoseChildrenAreAtBottom() const { return top.getNoGrandVector(); }
  inline NodeVector Tree::getSwappableNodes() const { return top.getSwappableVector(); }
}

#endif
