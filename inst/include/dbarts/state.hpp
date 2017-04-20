#ifndef DBARTS_STATE_HPP
#define DBARTS_STATE_HPP

#include <cstddef>

#include <external/random.h>

namespace dbarts {
  struct Tree;
  struct BARTFit;
  
  struct State {
    Tree* trees;
    std::size_t* treeIndices; // numObs * numTrees
    
    double* treeFits;      // numObs x numTrees;     vals for tree <=> x + i * numObs
    double* totalFits;     // numObs
    double* totalTestFits; // numTestObs x 1

    double sigma;
    
    double runningTime;
    
    ext_rng* rng;
    
    const char* const* createTreeStrings(const BARTFit& fit) const;
    void recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings);
  };
} // namespace dbarts

#endif // DBARTS_STATE_HPP
