#ifndef BART_STATE_HPP
#define BART_STATE_HPP

#include <cstddef>

namespace bart {
  struct Tree;
  struct BARTFit;
  
  struct State {
    Tree* trees;
    size_t* treeIndices;
    
    double* treeFits;      // numObs x numTrees;     vals for tree <=> x + i * numObs
    double* totalFits;
    double* totalTestFits; // numTestObs x 1

    double sigma;
    
    const char* const* createTreeStrings(const BARTFit& fit) const;
    void recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings);
  };
} // namespace bart

#endif // BART_STATE_HPP
