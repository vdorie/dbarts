#ifndef DBARTS_STATE_HPP
#define DBARTS_STATE_HPP

#include <cstddef>

#include <dbarts/random.hpp>

namespace dbarts {
  struct Control;
  struct Data;
  struct Tree;
  struct SavedTree;
  struct BARTFit;
  
  struct State {
    std::size_t* treeIndices; // numObs x numTrees
    Tree* trees;              // numTrees
    double* treeFits;         // numObs x numTrees; vals for tree <=> obsNum + treeNum * numObs
    
    SavedTree* savedTrees;              // numTrees x numSamples
    
    double sigma;
    double k;
    
    ext_rng* rng;
    
    State(const Control& control, const Data& data, size_t treeFitsStride, bool allocateTreeFitsAligned);
    void invalidate(std::size_t numTrees, std::size_t numSamples, bool treeFitsWereAllocatedAligned);
    
    // returns true if resize was necessary
    bool resize(const BARTFit& fit, const Control& newControl);
    bool resize(const BARTFit& fit, std::size_t numSamples);
    
    std::size_t getSerializedTreesLength(const BARTFit& fit) const;
    void serializeTrees(const BARTFit& fit, void* state) const;
    void deserializeTrees(const BARTFit& fit, const void* state);
    
    std::size_t getSerializedSavedTreesLength(const BARTFit& fit) const;
    void serializeSavedTrees(const BARTFit& fit, void* state) const;
    void deserializeSavedTrees(const BARTFit& fit, const void* state);
  };
} // namespace dbarts

#endif // DBARTS_STATE_HPP
