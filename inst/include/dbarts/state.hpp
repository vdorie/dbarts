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
    Tree* trees;              // numTrees x numSamples
    std::size_t* treeIndices; // numObs x numTrees x numSamples
    
    double* treeFits;      // numObs x numTrees x numSamples; vals for tree <=> obsNum + treeNum * numObs + sampleNum * numTrees * numSamples

    double* sigma; // 1 x numSamples
    
    ext_rng* rng;
    
    State(const Control& control, const Data& data);
    void invalidate(const Control& contol, std::size_t numSamples);
    ~State();
    
    bool resize(const BARTFit& fit, const Control& newControl);
    bool resize(const BARTFit& fit, std::size_t numSamples);
    
    const char* const* createTreeStrings(const BARTFit& fit) const;
    void recreateTreesFromStrings(const BARTFit& fit, const char* const* treeStrings);
  };
} // namespace dbarts

#endif // DBARTS_STATE_HPP
