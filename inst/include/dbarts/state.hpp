#ifndef DBARTS_STATE_HPP
#define DBARTS_STATE_HPP

#include <cstddef>

#include <external/random.h>

namespace dbarts {
  struct Control;
  struct Data;
  struct Tree;
  struct BARTFit;
  
  struct State {
    std::size_t* treeIndices; // numObs x numTrees
    Tree* trees;              // numTrees
    double* treeFits;         // numObs x numTrees; vals for tree <=> obsNum + treeNum * numObs
    
    std::size_t* savedTreeIndices; // numObs x numTrees x numSamples
    Tree* savedTrees;              // numTrees x numSamples
    double* savedTreeFits;         // numObs x numTrees x numSamples; vals for tree <=> obsNum + treeNum * numObs + sampleNum * numTrees * numSamples

    double sigma;
    
    ext_rng* rng;
    
    State(const Control& control, const Data& data);
    void invalidate(std::size_t numTrees, std::size_t numSamples);
    
    // returns true if resize was necessary
    bool resize(const BARTFit& fit, const Control& newControl);
    bool resize(const BARTFit& fit, std::size_t numSamples);
    
    const char* const* createTreeStrings(const BARTFit& fit, bool useSavedTrees) const;
    void recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings, bool useSavedTrees);
  };
} // namespace dbarts

#endif // DBARTS_STATE_HPP
