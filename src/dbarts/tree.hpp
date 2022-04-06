#ifndef DBARTS_TREE_HPP
#define DBARTS_TREE_HPP

#include <cstddef>
#include <dbarts/cstdint.hpp>
#include <dbarts/types.hpp>

#include "node.hpp"

extern "C" struct ext_rng;

namespace dbarts {
  struct BARTFit;
  struct SavedTree;
  
  struct Tree {
    Node top;
    
    Tree(std::size_t* indices, std::size_t numObservations, std::size_t numPredictors) : top(indices, numObservations, numPredictors) { }
    
    void sampleParametersAndSetFits(const BARTFit& fit, std::size_t chainNum, double* trainingFits, double* testFits);
    double* recoverParametersFromFits(const BARTFit& fit, const double* treeFits); // allocates result; are ordered as bottom nodes are
    double* recoverParametersFromFits(const BARTFit& fit, const double* treeFits, std::size_t* numBottomNodes); // allocates result; are ordered as bottom nodes are
    void setCurrentFitsFromParameters(const BARTFit& fit, const double* nodeParams, double* trainingFits, double* testFits);
    void setCurrentFitsFromParameters(const BARTFit& fit, const double* nodeParams, const xint_t* xt, std::size_t numObservations, double* fits);
    
    // deals largely if there are a different number of cut points, since a tree could then conceivably have
    // splits out of range
    void mapOldCutPointsOntoNew(const BARTFit& fit, const double* const* oldCutPoints, double* nodeParams);
    
    void collapseEmptyNodes(); // this ignores parameters and should be used on an invalid tree
    void collapseEmptyNodes(const BARTFit& fit, double* nodeParams); // this combines the parameters in the nodes
    
    void sampleStructureFromPrior(const BARTFit& fit, ext_rng* rng);
    void sampleParametersFromPrior(const BARTFit& fit, std::size_t chainNum, double* trainingFits, double* testFits);
    
    Node* getTop() const;
    bool hasSingleNode() const;
    
    std::size_t getNumBottomNodes() const;
    std::size_t getNumNotBottomNodes() const;
    std::size_t getNumNodesWhoseChildrenAreBottom() const;
    std::size_t getNumSwappableNodes() const;
    
    NodeVector getBottomNodes() const;
    void fillBottomNodesVector(NodeVector& nodeVector) const;
    NodeVector getNotBottomNodes() const;
    void fillNotBottomNodesVector(NodeVector& nodeVector) const;
    NodeVector getNodesWhoseChildrenAreAtBottom() const;
    NodeVector getSwappableNodes() const;
    void fillSwappableNodesVector(NodeVector& nodeVector) const;
    
    void setNodeAverages(const BARTFit& fit, std::size_t chainNum, const double* y);
    
    void countVariableUses(std::uint32_t* variableCounts) const;
    
    const char* createString() const;
    
    std::size_t getSerializedLength(const BARTFit& fit) const;
    std::size_t serialize(const BARTFit& fit, void* state) const;
    std::size_t deserialize(const BARTFit& fit, const void* state);
    
    bool isValid() const;
  };
  
  struct SavedTree {
    SavedNode top;
    
    SavedTree() : top() { }
    void copyStructureFrom(const BARTFit& fit, const Tree& other, const double* treeFits);
    
    void getPredictions(const BARTFit& fit, const double* xt, std::size_t numTestObservations, double* result);
    
    std::size_t getSerializedLength() const;
    std::size_t serialize(void* state) const;
    std::size_t deserialize(const void* state);
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  inline Node* Tree::getTop() const { return const_cast<Node*>(&top); }
  inline bool Tree::hasSingleNode() const { return top.isBottom(); }
  
  inline std::size_t Tree::getNumBottomNodes() const { return top.getNumBottomNodes(); }
  inline std::size_t Tree::getNumNotBottomNodes() const { return top.getNumNotBottomNodes(); }
  inline std::size_t Tree::getNumNodesWhoseChildrenAreBottom() const { return top.getNumNoGrandNodes(); }
  inline std::size_t Tree::getNumSwappableNodes() const { return top.getNumSwappableNodes(); }
  
  inline NodeVector Tree::getBottomNodes() const { return top.getBottomVector(); }
  inline void Tree::fillBottomNodesVector(NodeVector& nodeVector) const { return top.fillBottomVector(nodeVector); }
  inline NodeVector Tree::getNotBottomNodes() const { return top.getNotBottomVector(); }
  inline void Tree::fillNotBottomNodesVector(NodeVector& nodeVector) const { return top.fillNotBottomVector(nodeVector); }
  inline NodeVector Tree::getNodesWhoseChildrenAreAtBottom() const { return top.getNoGrandVector(); }
  inline NodeVector Tree::getSwappableNodes() const { return top.getSwappableVector(); }
  inline void Tree::fillSwappableNodesVector(NodeVector& nodeVector) const { return top.fillSwappableVector(nodeVector); }
}

#endif
