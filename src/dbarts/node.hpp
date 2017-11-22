#ifndef DBARTS_NODE_HPP
#define DBARTS_NODE_HPP

#include <dbarts/cstdint.hpp>
#include <cstddef>
#include <vector>

#include <dbarts/types.hpp>

struct ext_rng;

namespace dbarts {
  struct BARTFit;
  struct EndNodePrior;
  
#define DBARTS_INVALID_RULE_VARIABLE -1
  struct Rule {
    std::int32_t variableIndex;
    
    union {
      std::int32_t splitIndex;
      std::uint32_t categoryDirections;
    };
    
    void invalidate();
    
    bool goesRight(const BARTFit& fit, const double* x) const;
    bool categoryGoesRight(std::uint32_t categoryId) const;
    void setCategoryGoesRight(std::uint32_t categoryId);
    void setCategoryGoesLeft(std::uint32_t categoryId);
    double getSplitValue(const BARTFit& fit) const;
    
    bool equals(const Rule& other) const;
    void copyFrom(const Rule& other);
    void swapWith(Rule& other);
  };
  
  struct VarUsage {
    std::uint32_t depth;
    std::size_t nodeIndex;
    std::uint32_t variableIndex;
  };
  
  struct Node;
  
  struct ParentMembers {
    Node* rightChild;
    Rule rule;
  };
  
  struct EndNodeMembers {
    double average;
    double numEffectiveObservations;
  };
  
  typedef std::vector<Node*> NodeVector;
  
  struct Node {
    Node* parent;
    Node* leftChild;
    
    union {
      ParentMembers p;
      EndNodeMembers m;
    };
    
#define BART_INVALID_NODE_ENUM static_cast<std::size_t>(-1)
    std::size_t enumerationIndex;
    bool* variablesAvailableForSplit;
    
    std::size_t* observationIndices;
    std::size_t numObservations;
    
    Node(std::size_t* observationIndices, std::size_t numObservations, std::size_t numPredictors); // node is assumed at top
    Node(const Node& parent, std::size_t numPredictors); // node attaches to parent; parent should add observations
    ~Node();
    
    void copyFrom(const BARTFit& fit, const Node& other);
    
    
    bool isTop() const;
    bool isBottom() const;
    bool childrenAreBottom() const;
    
    Node* getParent() const;
    Node* getLeftChild() const;
    Node* getRightChild() const;
    
    std::size_t getNumBottomNodes() const;
    std::size_t getNumNotBottomNodes() const;
    std::size_t getNumNoGrandNodes() const;
    std::size_t getNumSwappableNodes() const;
    
    NodeVector getBottomVector() const;
    NodeVector getNoGrandVector() const;
    NodeVector getNotBottomVector() const;
    NodeVector getSwappableVector() const;
    
    void enumerateBottomNodes();
    NodeVector getAndEnumerateBottomVector(); // the nodes will have their enumeration indices set to their array index
    
    Node* findBottomNode(const BARTFit& fit, const double* x) const;
        
    void print(const BARTFit& fit) const;
    
    void setAverage(double average);                       // call these only on bottom nodes
    void setAverage(const BARTFit& fit, std::size_t chainNum, const double* y);  //
    void setAverages(const BARTFit& fit, std::size_t chainNum, const double* y); // call anywhere and it'll recurse
    void setNumEffectiveObservations(double n);
    
    double getAverage() const;
    double getNumEffectiveObservations() const;
    double computeVariance(const BARTFit& fit, std::size_t chainNum, const double* y) const;
    
    std::size_t getNumObservations() const;
    void addObservationsToChildren(const BARTFit& fit);
    void addObservationsToChildren(const BARTFit& fit, std::size_t chainNum, const double* y); // computes averages in bottom nodes as it goes
    void setObservationIndices(std::size_t* indices);
    void clearObservations();
    void clear();
    
    double drawFromPosterior(ext_rng* rng, const EndNodePrior& endNodePrior, double residualVariance) const;
    void setPredictions(double* y_hat, double prediction) const;
        
    std::size_t getDepth() const;
    std::size_t getDepthBelow() const;
    
    std::size_t getNumNodesBelow() const;
    std::size_t getNumVariablesAvailableForSplit(std::size_t numVariables) const;
    
    void split(const BARTFit& fit, const Rule& rule, bool exhaustedLeftSplits, bool exhaustedRightSplits);
    void split(const BARTFit& fit, std::size_t chainNum, const Rule& rule, const double* y, bool exhaustedLeftSplits, bool exhaustedRightSplits);
    void orphanChildren();
    
    void countVariableUses(std::uint32_t* variableCounts) const;
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // here for inlining purposes but contain implementation details and shouldn't be
  // relied on
  inline bool Node::isTop() const { return parent == NULL; }
  inline bool Node::isBottom() const { return leftChild == NULL; }
  inline bool Node::childrenAreBottom() const { return leftChild != NULL && leftChild->leftChild == NULL && p.rightChild->leftChild == NULL; }
  
  inline Node* Node::getParent() const { return const_cast<Node*>(parent); }
  inline Node* Node::getLeftChild() const { return const_cast<Node*>(leftChild); }
  inline Node* Node::getRightChild() const { return const_cast<Node*>(p.rightChild); }

  inline std::size_t Node::getNumObservations() const { return numObservations; }
  inline double Node::getAverage() const { return m.average; }

#ifdef MATCH_BAYES_TREE
  // This only means something if weights are supplied, which BayesTree didn't have.
  // It is also only meaningful on non-end nodes when using MATCH_BAYES_TREE.
  inline double Node::getNumEffectiveObservations() const { if (leftChild == NULL) return m.numEffectiveObservations; else return leftChild->getNumEffectiveObservations() + p.rightChild->getNumEffectiveObservations(); }
#else
  inline double Node::getNumEffectiveObservations() const { return m.numEffectiveObservations; }
#endif
  inline void Node::setAverage(double newAverage) { leftChild = NULL; m.average = newAverage; }
  inline void Node::setNumEffectiveObservations(double n) { leftChild = NULL; m.numEffectiveObservations = n; }
  inline void Node::setObservationIndices(std::size_t* indices) { observationIndices = indices; }
  
  inline bool Rule::categoryGoesRight(std::uint32_t categoryId) const { return ((1u << categoryId) & categoryDirections) != 0; }
  inline void Rule::setCategoryGoesRight(std::uint32_t categoryId) { categoryDirections |= (1u << categoryId); }
  inline void Rule::setCategoryGoesLeft(std::uint32_t categoryId) { categoryDirections &= ~(1u << categoryId); }
}

#endif
