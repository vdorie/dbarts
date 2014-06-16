#ifndef BART_NODE_HPP
#define BART_NODE_HPP

#include <bart/cstdint>
#include <cstddef>
using std::size_t;
using std::uint32_t;

#include <vector>

#include <bart/types.hpp>

namespace bart {
  struct BARTFit;
  struct EndNodePrior;
  
#define BART_INVALID_RULE_VARIABLE -1
  struct Rule {
    int32_t variableIndex;
    
    union {
      int32_t splitIndex;
      uint32_t categoryDirections;
    };
    
    Rule();
    Rule(const Rule& other); // copy constructor; duplicates
    
    void invalidate();
    
    bool goesRight(const BARTFit& fit, const double* x) const;
    bool categoryGoesRight(uint32_t categoryId) const;
    void setCategoryGoesRight(uint32_t categoryId);
    void setCategoryGoesLeft(uint32_t categoryId);
    double getSplitValue(const BARTFit& fit) const;
    
    bool equals(const BARTFit& fit, const Rule& other) const;
    void copyFrom(const BARTFit& fit, const Rule& other);
    void swapWith(Rule& other);
  };
  
  struct VarUsage {
    uint32_t depth;
    size_t nodeIndex;
    uint32_t variableIndex;
  };
  
  struct Node;
  typedef std::vector<Node*> NodeVector;
  
  struct Node {
    Node* parent;
    Node* leftChild;
    union {
      Node* rightChild; // average only applies to nodes w/o children
      double average;
    };
    
    Rule rule;
#define BART_INVALID_NODE_ENUM ((size_t) -1)
    size_t enumerationIndex;
    bool* variablesAvailableForSplit;
    
    size_t* observationIndices;
    size_t numObservationsInNode;
    
    Node(size_t* observationIndices, size_t numObservationsInNode, size_t numPredictors); // node is assumed at top
    Node(const Node& parent, size_t numPredictors); // node attaches to parent; parent should add observations
    ~Node();
    
    void copyFrom(const BARTFit& fit, const Node& other);
    
    
    bool isTop() const;
    bool isBottom() const;
    bool childrenAreBottom() const;
    
    
    size_t getNumBottomNodes() const;
    size_t getNumNotBottomNodes() const;
    size_t getNumNoGrandNodes() const;
    size_t getNumSwappableNodes() const;
    
    NodeVector getBottomVector() const;
    NodeVector getNoGrandVector() const;
    NodeVector getNotBottomVector() const;
    NodeVector getSwappableVector() const;
    
    NodeVector getAndEnumerateBottomVector(); // the nodes will have their enumeration indices set to their array index
    
    Node* findBottomNode(const BARTFit& fit, const double* x) const;
        
    void print(const BARTFit& fit) const;
    
    void setAverage(double average);                       // call these only on bottom nodes
    void setAverage(const BARTFit& fit, const double* y);  //
    void setAverages(const BARTFit& fit, const double* y); // call anywhere and it'll recurse
    double getAverage() const;
    double computeVariance(const BARTFit& fit, const double* y) const;
    
    size_t getNumObservationsInNode() const;
    void addObservationsToChildren(const BARTFit& fit);
    void addObservationsToChildren(const BARTFit& fit, const double* y); // computes averages in bottom nodes as it goes
    void clearObservations();
    void clear();
    
    double drawFromPosterior(const EndNodePrior& endNodePrior, double residualVariance) const;
    void setPredictions(double* y_hat, double prediction) const;
        
    size_t getDepth() const;
    size_t getDepthBelow() const;
    
    size_t getNumNodesBelow() const;
    size_t getNumVariablesAvailableForSplit(size_t numVariables) const;
    
    void split(const BARTFit& fit, const Rule& rule, const double* y, bool exhaustedLeftSplits, bool exhaustedRightSplits);
    void orphanChildren();
    
    void countVariableUses(uint32_t* variableCounts) const;
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // here for inlining purposes but contain implementation details and shouldn't be
  // relied on
  inline bool Node::isTop() const { return parent == NULL; }
  inline bool Node::isBottom() const { return leftChild == NULL; }
  inline bool Node::childrenAreBottom() const { return leftChild != NULL && leftChild->leftChild == NULL && rightChild->leftChild == NULL; }

  inline size_t Node::getNumObservationsInNode() const { return numObservationsInNode; }
  inline double Node::getAverage() const { return average; }
  inline void Node::setAverage(double newAverage) { leftChild = NULL; average = newAverage; }
  
  inline bool Rule::categoryGoesRight(uint32_t categoryId) const { return ((1u << categoryId) & categoryDirections) != 0; }
  inline void Rule::setCategoryGoesRight(uint32_t categoryId) { categoryDirections |= (1u << categoryId); }
  inline void Rule::setCategoryGoesLeft(uint32_t categoryId) { categoryDirections &= ~(1u << categoryId); }
}

#endif
